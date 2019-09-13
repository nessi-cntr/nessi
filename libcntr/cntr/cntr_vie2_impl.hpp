#ifndef CNTR_VIE2_IMPL_H
#define CNTR_VIE2_IMPL_H

#include "cntr_vie2_decl.hpp"
#include "fourier.hpp"
#include "cntr_elements.hpp"
#include "cntr_herm_matrix_timestep_decl.hpp"
#include "cntr_herm_matrix_decl.hpp"
#include "cntr_matsubara_impl.hpp"
#include "cntr_convolution_impl.hpp"

namespace cntr {

#define CPLX std::complex<T>
/* #######################################################################################
#
#  [1+F]G=Q, G,Q hermitian
#
###########################################################################################*/

/*###########################################################################################
#
#   RETARDED FUNCTION:
#
###########################################################################################*/

/** \brief <b> One step VIE solver \f$(1+F)*G=Q\f$ for the retarded component of the Green's function \f$G\f$ at a given timestep</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for the retarded component \f$G^R(t, t^\prime)\f$ at a given timestep, for given:
* > input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$, the source term \f$Q(t, t^\prime)\f$, and
* > the integrator class 'I'.
*
* \note: F,G and Q are instances of the template class `GG`, representing `herm_matrix`.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] time step
* @param &G
* > [GG] Solution
* @param &F
* > [GG] Green's function on left-hand side
* @param &Fcc
* > [GG] Hermitian conjugate of F
* @param &Q
* > [GG] Green's function on right-hand side
* @param I
* > [Integrator] integrator class
* @param h
* > [double] time interval
*/
template <typename T, class GG, int SIZE1>
void vie2_timestep_ret(int n, GG &G, GG &F, GG &Fcc, GG &Q, integration::Integrator<T> &I,
                       T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1;
    int ss, sg, n1, l, j, i, p, q, m, j1, j2, size1 = G.size1();
    cplx *gret, *gtemp, *mm, *qq, *qqj, *stemp, *one, cplx_i, cweight, *diffw;
    cplx *sret, w0;
    T weight;
    cplx_i = cplx(0, 1);

    ss = F.element_size();
    sg = G.element_size();
    gtemp = new cplx[k * sg];
    diffw = new cplx[k1 + 1];
    qq = new cplx[(n + 1) * sg];
    one = new cplx[sg];
    mm = new cplx[k * k * sg];
    stemp = new cplx[sg]; // sic
    element_set<T, SIZE1>(size1, one, 1);
    // check consistency:
    assert(n > k);
    assert(F.nt() >= n);
    assert(G.nt() >= n);
    assert(G.sig() == F.sig());
    // SET ENTRIES IN TIMESTEP TO 0
    gret = G.retptr(n, 0);
    n1 = (n + 1) * sg;
    for (i = 0; i < n1; i++)
        gret[i] = 0;
    // INITIAL VALUE t' = n
    element_set<T, SIZE1>(size1, G.retptr(n, n), Q.retptr(n, n));
    // START VALUES  t' = n-j, j = 1...k: solve a kxk problem
    for (i = 0; i < k * k * sg; i++)
        mm[i] = 0;
    for (i = 0; i < k * sg; i++)
        qq[i] = 0;
    for (j = 1; j <= k; j++) {
        p = j - 1;
        element_incr<T, SIZE1>(size1, qq + p * sg, Q.retptr(n, n - j));
        element_incr<T, SIZE1>(size1, mm + sg * (p + k * p), one);
        // integral
        for (l = 0; l <= k; l++) {
            weight = -h * I.gregory_weights(j, l);
            if (n - l >= n - j) {
                element_set<T, SIZE1>(
                    size1, stemp, F.retptr(n - l, n - j)); // stemp is element of type G!!
            } else {
                element_set<T, SIZE1>(size1, stemp, Fcc.retptr(n - j, n - l));
                element_conj<T, SIZE1>(size1, stemp);
                weight *= -1;
            }
            if (l == 0) {
                element_incr<T, SIZE1>(size1, qq + p * sg, weight, G.retptr(n, n), stemp);
            } else {
                q = l - 1;
                element_incr<T, SIZE1>(size1, mm + sg * (p * k + q), -weight, stemp);
            }
        }
    }
    element_linsolve_left<T, SIZE1>(size1, k, gtemp, mm, qq); // gtemp * mm = qq
    for (j = 1; j <= k; j++)
        element_set<T, SIZE1>(size1, G.retptr(n, n - j), gtemp + (j - 1) * sg);
    // Compute the contribution to the convolution int dm G(n,m)Sigma(m,j)
    // from G(n,m=n-k..n) to m=0...n-k, store into qq(j), without factor h!!
    for (i = 0; i < sg * (n + 1); i++)
        qq[i] = 0;
    for (m = n - k; m <= n; m++) {
        weight = -I.gregory_omega(n - m);
        gret = G.retptr(n, m);
        sret = F.retptr(m, 0);
        for (j = 0; j <= n - k2; j++) {
            element_incr<T, SIZE1>(size1, qq + j * sg, -I.gregory_weights(n - j, n - m),
                                   gret, F.retptr(m, j));
            sret += ss;
        }
        j1 = (n - k2 + 1 < 0 ? 0 : n - k2 + 1);
        for (j = j1; j < n - k; j++) { // start weights
            weight = I.gregory_weights(n - j, n - m);
            element_incr<T, SIZE1>(size1, qq + j * sg, -I.gregory_weights(n - j, n - m),
                                   gret, F.retptr(m, j));
            sret += ss;
        }
    }
    // DER REST VOM SCHUETZENFEST: t' = n-l, l = k+1 ... n:
    w0 = -h * I.gregory_omega(0);
    for (l = k + 1; l <= n; l++) {
        j = n - l;
        // set up mm and qqj for 1x1 problem:
        qqj = qq + j * sg;
        for (i = 0; i < sg; i++) {
            qqj[i] *= h;
        }
        element_set<T, SIZE1>(size1, stemp, F.retptr(j, j));
        for (i = 0; i < sg; i++)
            mm[i] = one[i] - w0 * stemp[i];
        element_incr<T, SIZE1>(size1, qqj, Q.retptr(n, j));
        element_linsolve_left<T, SIZE1>(size1, G.retptr(n, j), mm, qqj);
        // compute the contribution of Gret(n,j) to
        // int dm G(n,j)Sigma(j,j2)  for j2<j, store into qq(j2), without h
        gret = G.retptr(n, j);
        sret = F.retptr(j, 0);
        for (j2 = 0; j2 < j - k; j2++) {
            element_incr<T, SIZE1>(size1, qq + j2 * sg, -I.gregory_weights(n - j2, n - j),
                                   gret, sret);
            sret += ss;
        }
        j1 = (j - k < 0 ? 0 : j - k);
        for (j2 = j1; j2 < j; j2++) { // start weights
            weight = -I.gregory_omega(j - j2);
            element_incr<T, SIZE1>(size1, qq + j2 * sg, -I.gregory_weights(n - j2, n - j),
                                   gret, sret);
            sret += ss;
        }
    }
    delete[] diffw;
    delete[] stemp;
    delete[] qq;
    delete[] mm;
    delete[] gtemp;
    delete[] one;
    return;
}
/** \brief <b> VIE solver \f$(1+F)*G=Q\f$ for a left-mixing component of the Green's function \f$G\f$ for the first k timesteps</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for a left-mixing component of the Green's function \f$G^\rceil(t, t^\prime)\f$
* > for the first k timesteps (given by the integrator class 'I').
* > Here, are given: input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$,
* > and a source term \f$Q(t, t^\prime)\f$.
*
* \note: F,G and Q are instances of the template class `GG`, representing `herm_matrix`.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [GG] solution
* @param &F
* > [GG] green's function  on left-hand side
* @param &Fcc
* > [GG] Complex conjugate of F
* @param &Q
* > [GG] green's function  on right-hand side
* @param I
* > [Integrator] integrator class
* @param h
* > [double] time interval
*/
template <typename T, class GG, int SIZE1>
void vie2_start_ret(GG &G, GG &F, GG &Fcc, GG &Q, integration::Integrator<T> &I, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1;
    int sg, l, n, j, p, q, i, dim, size1 = G.size1();
    cplx ih, *gtemp, *mm, *qq, *stemp, minusi, *one;
    T weight;

    sg = G.element_size();
    assert(F.nt() >= k);
    assert(G.nt() >= k);
    assert(G.sig() == F.sig());

    // temporary storage
    qq = new cplx[k1 * sg];
    mm = new cplx[k1 * k1 * sg];
    one = new cplx[sg];
    gtemp = new cplx[k1 * sg];
    stemp = new cplx[sg]; // sic

    element_set<T, SIZE1>(size1, one, 1.0);
    minusi = cplx(0, -1);

    // set initial values:
    for (n = 0; n <= k; n++)
        element_set<T, SIZE1>(size1, G.retptr(n, n), Q.retptr(n, n));

    for (j = 0; j < k; j++) { // determine G(n,j), n=j+1..k
        dim = k - j;
        for (i = 0; i < k * sg; i++)
            qq[i] = 0;
        for (i = 0; i < k * k * sg; i++)
            mm[i] = 0;
        // fill the matrix mm, indices p,q
        for (n = j + 1; n <= k; n++) {
            p = n - (j + 1);
            element_incr<T, SIZE1>(size1, qq + p * sg, Q.retptr(n, j));
            for (l = 0; l <= k; l++) {
                q = l - (j + 1);
                if (l <= j) { // this involves G which is known and goes into qq(p)
                    element_conj<T, SIZE1>(size1, gtemp, G.retptr(j, l));
                    element_smul<T, SIZE1>(size1, gtemp, -1);
                    element_set_zero<T, SIZE1>(size1, stemp);
                    element_incr<T, SIZE1>(size1, stemp, F.retptr(n, l), gtemp);
                    for (i = 0; i < sg; i++)
                        qq[p * sg + i] += -h * I.poly_integration(j, n, l) * stemp[i];
                } else { // this goes into mm(p,q)
                    for (i = 0; i < sg; i++)
                        mm[sg * (p * dim + q) + i] = 0.0;
                    if (l == n) {
                        for (i = 0; i < sg; i++) {
                            mm[sg * (p * dim + q) + i] += one[i];
                        }
                    }
                    weight = -h * I.poly_integration(j, n, l);
                    if (n >= l) {
                        element_set<T, SIZE1>(size1, stemp, F.retptr(n, l));
                    } else {
                        element_set<T, SIZE1>(size1, stemp, Fcc.retptr(l, n));
                        element_conj<T, SIZE1>(size1, stemp);
                        weight *= -1;
                    }
                    for (i = 0; i < sg; i++)
                        mm[sg * (p * dim + q) + i] -= weight * stemp[i];
                }
            }
        }
        element_linsolve_right<T, SIZE1>(size1, dim, gtemp, mm, qq);
        for (n = j + 1; n <= k; n++) {
            p = n - (j + 1);
            element_set<T, SIZE1>(size1, G.retptr(n, j), gtemp + p * sg);
        }
    }
    delete[] stemp;
    delete[] qq;
    delete[] mm;
    delete[] gtemp;
    delete[] one;
    return;
}
/*####################################################################################
#
#   MATSUBARA FUNCTION g(tau) = 1/beta sum_n e^{-iomn tau} (iomn-H-Sigma(iomn))^{-1}
#
####################################################################################*/
/// @private
// TODO:
// For Sigma=0 the boundary values seem to be a very badly approximated, although
// the jump=1 is accounted far correctly. Sigma!=0 somehow makes the performance
// better
template <typename T, class GG, int SIZE1>
void vie2_mat_fourier_dispatch(GG &G, GG &F, GG &Fcc, GG &Q, T beta, int pcf = 20, int order = 3) {
    typedef std::complex<T> cplx;
    cplx *fmdft, *fiomn, *qiomn, *qiomn1, *qmdft, *qmasy, *z1, *z2, *z3, *zinv, *one;
    cplx *expfac, *xmat;
    int ntau, m, r, p, m2, sig, size1 = G.size1(), l, sg, ss, matsub_one;
    T dtau, omn;
    // T arg;

    assert(G.ntau() == F.ntau());
    sig = G.sig();
    assert(sig == -1 || sig == 1); // Tested for both Bosons and Fermions (StrandH)
    assert(G.sig() == F.sig());
    sg = G.element_size();
    ss = F.element_size();
    ntau = G.ntau();
    dtau = beta / ntau;

    matsub_one = (sig == -1 ? 1.0 : 0.0); // 1 for Fermions, 0 for Bosons

    // pcf=10;
    if (ntau % 2 == 1) {
        std::cerr << "matsubara_inverse: ntau odd" << std::endl;
        abort();
    }
    m2 = ntau / 2;
    xmat = new cplx[(ntau + 1) * sg];
    fmdft = new cplx[(ntau + 1) * sg];
    qmdft = new cplx[(ntau + 1) * sg];
    expfac = new cplx[ntau + 1];
    qmasy = new cplx[sg];
    z1 = new cplx[sg];
    z2 = new cplx[sg];
    z3 = new cplx[sg];
    zinv = new cplx[sg];
    fiomn = new cplx[sg];
    qiomn = new cplx[sg];
    qiomn1 = new cplx[sg];
    one = new cplx[sg];
    element_set<T, SIZE1>(size1, one, 1.0);


    // First order tail correction factor
    // qmasy = -1 * ( Q(0) + Q(\beta) ) NB! Only valid for Fermions
    // The general expression reads: qmasy = sig * Q(\beta) - Q(0)

    element_set<T, SIZE1>(size1, qmasy, Q.matptr(0));
    element_smul<T, SIZE1>(size1, qmasy, -1.0);
    element_set<T, SIZE1>(size1, z1, Q.matptr(ntau));
    element_smul<T, SIZE1>(size1, z1, sig);
    element_incr<T, SIZE1>(size1, qmasy, z1);

    // Fermion specific code
    // element_set<T,SIZE1>(size1,qmasy,Q.matptr(0));
    // element_incr<T,SIZE1>(size1,qmasy,Q.matptr(ntau));
    // element_smul<T,SIZE1>(size1,qmasy,-1.0);

    set_first_order_tail<T, SIZE1>(xmat, qmasy, beta, sg, ntau, sig, size1);

    // Raw Discrete Fourier Transform of F(\tau) & G(\tau)
    // Sig determines whether the Matsubara frequencies are Fermionic or Bosonic
    // F(i\omega_n) computed for \omega_n \ge 0 (i.e. including zero for Bosons)
    matsubara_dft<T, GG, SIZE1>(fmdft, F, sig);
    matsubara_dft<T, GG, SIZE1>(qmdft, Q, sig);

    // Oversampling loop over Matsubara freqencies
    // nomega = (2*pcf + 1) * ntau
    for (p = -pcf; p <= pcf; p++) {

        // For Bosons we sample the middle interval including the zero frequency.
        // (For Fermions there is no need to fiddle with the frequency interval)
        int m_low = 0, m_high = 0;
        if (sig == 1) {
            m_high = (p < 0 ? 0 : 1);
            m_low = (p > 0 ? -1 : 0);
        }

        // Loop over matsubara frequencies
        for (m = -m2 - m_low; m <= m2 - 1 + m_high; m++) {

            // omn=(2*(m+p*ntau)+1)*PI/(ntau*dtau); // Fermionic Matsubara freq.
            omn = (2 * (m + p * ntau) + matsub_one) * PI /
                  (ntau * dtau); // General Matsubara freq.

            matsubara_ft<T, GG, SIZE1>(fiomn, m + p * ntau, F, fmdft, sig, beta, order);
            matsubara_ft<T, GG, SIZE1>(qiomn, m + p * ntau, Q, qmdft, sig, beta, order);

            // This is our tail correction in frequency space
            // z3 = -i/omn * qmasy = 1/(i*omn) * qmasy
            element_set<T, SIZE1>(size1, z3, qmasy);

            // For Bosons we need to skip the zeroth frequency
            if (sig == 1 && m + p * ntau == 0)
                element_smul<T, SIZE1>(size1, z3, cplx(0.0, 0.0));
            else
                element_smul<T, SIZE1>(size1, z3, cplx(0.0, -1.0 / omn));

            // qiomn1 = qiomn - 1/(i*omn) * qmasy, remove tail correction from q
            for (l = 0; l < sg; l++)
                qiomn1[l] = qiomn[l] - z3[l];

            // z1 = 1 + f
            for (l = 0; l < sg; l++)
                z1[l] = fiomn[l] + one[l];
            // z2 = f * z3 = f * qmasy/(i*omn) ?? multiply f with tail correction
            element_mult<T, SIZE1>(size1, z2, fiomn, z3);

            // add correction to z2 = qi - z2 = q - qmasy/(i*omn) - f * qmasy/(i*omn)
            // so z2 = q - (1 + f) * qmasy/(i*omn)
            for (l = 0; l < sg; l++)
                z2[l] = qiomn1[l] - z2[l];

            // zinv = 1./z1 = (1 + f)^{-1}
            element_inverse<T, SIZE1>(size1, zinv, z1);
            element_mult<T, SIZE1>(size1, z1, zinv, z2); // REIHENFOLGE ?? Oroa dig inte ;)
            element_smul<T, SIZE1>(size1, z1, 1.0 / beta);

            // z1 = 1/beta * zinv * z2 = 1/beta (1+f)^{-1} * z2
            // => z1 = 1/beta (1+f)^{-1} * [q - (1+f)*qmasy/(i*omn)]

            // So actually we could simplify this to
            // z1 = 1/beta * [ (1+f)^{-1} * q - qmasy/(i*omn) ] ?

            for (r = 0; r <= ntau; r++) { // Loop over tau

                double arg = omn * r * dtau;
                std::complex<double> expfactor = cplx(cos(arg), -sin(arg));
                element_set<T, SIZE1>(size1, z2, z1);
                element_smul<T, SIZE1>(size1, z2, expfactor);
                element_incr<T, SIZE1>(size1, xmat + r * sg, z2);

            } // End tau loop
        }     // End matsubara freq loop
    }         // End oversampling loop

    for (r = 0; r <= ntau; r++)
        element_set<T, SIZE1>(size1, G.matptr(r), xmat + r * sg);

    delete[] xmat;
    delete[] fmdft;
    delete[] qmdft;
    delete[] expfac;
    delete[] qmasy;
    delete[] z1;
    delete[] z2;
    delete[] z3;
    delete[] zinv;
    delete[] qiomn;
    delete[] qiomn1;
    delete[] fiomn;
    delete[] one;
    return;
}

/// @private
template <typename T, class GG, int SIZE1>
void vie2_mat_fixpoint_dispatch(GG &G, GG &F, GG &Fcc, GG &Q, T beta,
                             integration::Integrator<T> &I, int nfixpoint, int pcf = 5,
                             int order = 3) {
    int fixpoint, ntau = G.ntau(), size1 = G.size1(), r;

    vie2_mat_fourier_dispatch<T, GG, SIZE1>(G, F, Fcc, Q, beta, pcf, order);

    if (nfixpoint > 0) {
        GG Qn(-1, ntau, size1, Q.sig()), dG(-1, ntau, size1, G.sig());
        for (fixpoint = 0; fixpoint < nfixpoint; fixpoint++) {
            convolution_matsubara(Qn, F, G, I, beta);
            for (r = 0; r <= ntau; r++) {
                element_incr<T, SIZE1>(size1, Qn.matptr(r), G.matptr(r));
                element_incr<T, SIZE1>(size1, Qn.matptr(r), -1.0, Q.matptr(r));
            }
            vie2_mat_fourier_dispatch<T, GG, SIZE1>(dG, F, Fcc, Qn, beta, pcf, order);
            for (r = 0; r <= ntau; r++) {
                element_incr<T, SIZE1>(size1, G.matptr(r), -1.0, dG.matptr(r));
            }
        }
    }
    return;
}

/// @private
template <typename T, class GG, int SIZE1>
void vie2_mat_steep_dispatch(GG &G, GG &F, GG &Fcc, GG &Q, T beta,
			    integration::Integrator<T> &I, int maxiter, T maxerr, int pcf = 3,
			    int order = 3) {
    int iter, ntau = G.ntau(), size1 = G.size1(), r;
    std::complex<T> *temp = new std::complex<T>[size1 * size1];
    std::complex<T> *Rcc = new std::complex<T>[size1 * size1];
    std::complex<T> *Pcc = new std::complex<T>[size1 * size1];

    vie2_mat_fourier_dispatch<T, GG, SIZE1>(G, F, Fcc, Q, beta, pcf, order);

    if (maxiter > 0) {
      double alpha,c1,c2,err,rsold,rsnew,weight;
      GG R(-1, ntau, size1, Q.sig()), AR(-1, ntau, size1, Q.sig());
      GG P(-1, ntau, size1, Q.sig()), AP(-1, ntau, size1, Q.sig());
      GG C1(-1, ntau, size1, Q.sig());
      GG C2(-1, ntau, size1, Q.sig());

      //initial guees for residue
      convolution_matsubara(C1, F, G, I, beta);
      convolution_matsubara(C2, G, Fcc, I, beta);
      for (r = 0; r <= ntau; r++) {
          element_set<T, SIZE1>(size1, R.matptr(r), C1.matptr(r));
          element_incr<T, SIZE1>(size1, R.matptr(r), 1.0, C2.matptr(r));
          element_smul<T, SIZE1>(size1, R.matptr(r), 0.5);
	        element_incr<T, SIZE1>(size1, R.matptr(r), 1.0, G.matptr(r));
	        element_incr<T, SIZE1>(size1, R.matptr(r), -1.0, Q.matptr(r));
      }
      R.smul(-1,-1.0);
      for (r=0; r <= ntau; r++)
        element_set<T, SIZE1>(size1, P.matptr(r), R.matptr(r));


      // norm of residue
      rsold = 0.0;
      for (r = 0; r <= ntau; r++) {
         if(r <= I.get_k()) {
           weight = I.gregory_omega(r);
         } else if(r >= ntau - I.get_k()) {
           weight = I.gregory_omega(ntau-r);
         } else {
           weight = 1.0;
          }
          element_conj<T, SIZE1>(size1, Rcc, R.matptr(r));
          element_mult<T, SIZE1>(size1, temp, Rcc, R.matptr(r));
          for(int s=0; s<size1*size1; s++)
            rsold += weight * (temp[s]).real();
      }

      // if norm of residue < threshold, no iteration
      if(rsold > maxerr) {
        for (iter = 0; iter < maxiter; iter++) {

          // compute A.P
          convolution_matsubara(C1, F, P, I, beta);
          convolution_matsubara(C2, P, Fcc, I, beta);
          for (r = 0; r <= ntau; r++) {
            element_set<T, SIZE1>(size1, AP.matptr(r), C1.matptr(r));
            element_incr<T, SIZE1>(size1, AP.matptr(r), 1.0, C2.matptr(r));
            element_smul<T, SIZE1>(size1, AP.matptr(r), 0.5);
            element_incr<T, SIZE1>(size1, AP.matptr(r), 1.0, P.matptr(r));
          }

          // compute P.A.P
          c2 = 0.0;
          for (r = 0; r <= ntau; r++) {
            if(r <= I.get_k()) {
              weight = I.gregory_omega(r);
            } else if(r >= ntau - I.get_k()) {
              weight = I.gregory_omega(ntau-r);
            } else {
              weight = 1.0;
            }
            element_conj<T, SIZE1>(size1, Pcc, P.matptr(r));
            element_mult<T, SIZE1>(size1, temp, Pcc, AP.matptr(r));
            for(int s=0; s<size1*size1; s++)
              c2 += weight * (temp[s]).real();
          }

          alpha = rsold / c2;

          for (r = 0; r <= ntau; r++) {
            element_incr<T, SIZE1>(size1, G.matptr(r), alpha, P.matptr(r));
            element_incr<T, SIZE1>(size1, R.matptr(r), -alpha, AP.matptr(r));
          }

          rsnew = 0.0;
          for (r = 0; r <= ntau; r++) {
            if(r <= I.get_k()) {
              weight = I.gregory_omega(r);
            } else if(r >= ntau - I.get_k()) {
              weight = I.gregory_omega(ntau-r);
            } else {
              weight = 1.0;
            }
            element_conj<T, SIZE1>(size1, Rcc, R.matptr(r));
            element_mult<T, SIZE1>(size1, temp, Rcc, R.matptr(r));
            for(int s=0; s<size1*size1; s++)
              rsnew += weight * (temp[s]).real();
          }

          for (r = 0; r <= ntau; r++) {
            element_set<T, SIZE1>(size1, temp, R.matptr(r));
            element_incr<T, SIZE1>(size1, temp, rsnew/rsold, P.matptr(r));
            element_set<T, SIZE1>(size1, P.matptr(r), temp);
          }

          rsold = rsnew;

          if (sqrt(rsold) < maxerr) break;
        }
      }
      delete[] temp;
      delete[] Rcc;
      delete[] Pcc;
      return;
    }
}

/*###########################################################################################
#
#   tv-FUNCTION:
#
###########################################################################################*/
/** \brief <b> One step VIE solver \f$(1+F)*G=Q\f$ for a left-mixing component of the Green's function \f$G\f$ at a given timestep</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for a left-mixing component of \f$G^\rceil(t, t^\prime)\f$ at a given timestep, for given:
* > input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$, the source term \f$Q(t, t^\prime)\f$, and
* > the integrator class 'I'.
*
* \note: F,G and Q are instances of the template class `GG`, representing `herm_matrix`.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] time step
* @param &G
* > [GG] solution
* @param &F
* > [GG] green's function  on left-hand side
* @param &Fcc
* > [GG] Complex conjugate of F
* @param &Q
* > [GG] green's function  on right-hand side
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/
template <typename T, class GG, int SIZE1>
void vie2_timestep_tv(int n, GG &G, GG &F, GG &Fcc, GG &Q, integration::Integrator<T> &I,
                      T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k();
    int sg, n1, l, j, ntau, size1 = G.size1();
    cplx *gtv, cweight, ih, *mm, *qq, *stemp, minusi, *one;
    T weight;

    sg = G.element_size();
    ntau = G.ntau();
    // check consistency:  (more assertations follow in convolution)
    assert(n > k);
    assert(F.nt() >= n);
    assert(G.nt() >= n);
    assert(G.sig() == F.sig());

    one = new cplx[sg];
    qq = new cplx[sg];
    mm = new cplx[sg];
    stemp = new cplx[sg]; // sic
    element_set<T, SIZE1>(size1, one, 1.0);
    // SET ENTRIES IN TIMESTEP(TV) TO 0
    gtv = G.tvptr(n, 0);
    n1 = (ntau + 1) * sg;
    for (l = 0; l < n1; l++)
        gtv[l] = 0;
    // CONVOLUTION SIGMA*G:  --->  Gtv(n,m)
    convolution_timestep_tv<T, herm_matrix<T>, SIZE1>(n, G, F, Fcc, G, G, I, beta,
                                                      h); // note: this sets only tv
    // Now solve
    // [ 1 - h w(n,0) Sigma(n,n) ] G(n,m)  = Q(m),
    // where Q is initially stored in G(n,m)
    element_set<T, SIZE1>(size1, stemp, F.retptr(n, n));
    weight = h * I.gregory_weights(n, 0);
    for (l = 0; l < sg; l++)
        mm[l] = one[l] + weight * stemp[l];
    for (j = 0; j <= ntau; j++) {
        for (l = 0; l < sg; l++)
            qq[l] = -G.tvptr(n, j)[l] + Q.tvptr(n, j)[l];
        element_linsolve_right<T, SIZE1>(size1, G.tvptr(n, j), mm, qq);
    }
    delete[] stemp;
    delete[] qq;
    delete[] mm;
    delete[] one;
    return;
}

/** \brief <b> VIE solver \f$(1+F)*G=Q\f$ for a left-mixing component of the Green's function \f$G\f$ for the first k timesteps</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for a left-mixing component of \f$G^\rceil(t, t^\prime)\f$
* > for the first k timesteps (given by the integrator class 'I').
* > Here, are given: input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$,
* > and the source term \f$Q(t, t^\prime)\f$.
*
* \note: F,G and Q are instances of the template class `GG`, representing `herm_matrix`.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [GG] solution
* @param &F
* > [GG] green's function  on left-hand side
* @param &Fcc
* > [GG] Complex conjugate of F
* @param &Q
* > [GG] green's function  on right-hand side
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/
template <typename T, class GG, int SIZE1>
void vie2_start_tv(GG &G, GG &F, GG &Fcc, GG &Q, integration::Integrator<T> &I, T beta,
                   T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1;
    int sg, l, m, j, ntau, p, q, n, sig, size1 = G.size1();
    cplx cweight, *mm, *qq, *stemp, *one, *gtemp;
    cplx cplx_i = cplx(0.0, 1.0);
    T dtau;
    sg = G.element_size();
    ntau = G.ntau();
    dtau = beta / ntau;
    sig = G.sig();
    // check consistency:  (more assertations follow in convolution)
    assert(F.nt() >= k);
    assert(G.nt() >= k);
    assert(F.ntau() == ntau);
    assert(G.sig() == F.sig());

    one = new cplx[sg];
    qq = new cplx[k1 * sg];
    mm = new cplx[k1 * k1 * sg];
    stemp = new cplx[sg]; // sic
    gtemp = new cplx[k1 * sg];
    element_set<T, SIZE1>(size1, one, 1.0);

    // CONVOLUTION  -i int dtau Sigma^tv(n,tau)G^mat(tau,m) ---> Gtv(n,m)
    for (n = 0; n <= k; n++) {
        for (m = 0; m <= ntau; m++) {
            matsubara_integral_2<T, SIZE1>(size1, m, ntau, gtemp, F.tvptr(n, 0), G.matptr(0),
                                           I, G.sig());
            for (l = 0; l < sg; l++)
                G.tvptr(n, m)[l] = -dtau * gtemp[l];
        }
    }
    // loop over m --> determine G^tv(n,m) for n=0...k
    for (m = 0; m <= ntau; m++) {
        for (l = 0; l < k1 * k1 * sg; l++)
            mm[l] = 0;
        for (l = 0; l < k1 * sg; l++)
            qq[l] = 0;
        // derive linear equations
        // mm(p,q)*G(q)=Q(p) for p=n-1=0...k-1, q=n-1=0...k
        // G(p)=G(p+1,m)
        for (n = 0; n <= k; n++) {
            // 1 -- goes into m(p,p)
            element_incr<T, SIZE1>(size1, mm + sg * (n * k1 + n), one);
            // integral 0..n
            for (j = 0; j <= k; j++) {
                cweight = h * I.gregory_weights(n, j);
                if (n >= j) {
                    element_set<T, SIZE1>(size1, stemp, F.retptr(n, j));
                } else {
                    element_set<T, SIZE1>(size1, stemp, Fcc.retptr(j, n));
                    element_conj<T, SIZE1>(size1, stemp);
                    element_smul<T, SIZE1>(size1, stemp, -1);
                }
                for (l = 0; l < sg; l++)
                    mm[sg * (n * k1 + j) + l] += cweight * stemp[l];
            }
            // integral Sigmatv*Gmat --> take from Gtv(n,m), write into qq
            element_incr<T, SIZE1>(size1, qq + n * sg, G.tvptr(n, m));
            element_incr<T, SIZE1>(size1, qq + n * sg, Q.tvptr(n, m));
        }
        element_linsolve_right<T, SIZE1>(size1, k1, gtemp, mm,
                                         qq); // solve kXk problem mm*gtemp=qq
        // write elements into Gtv
        for (n = 0; n <= k; n++)
            element_set<T, SIZE1>(size1, G.tvptr(n, m), gtemp + n * sg);
    }
    delete[] stemp;
    delete[] gtemp;
    delete[] qq;
    delete[] mm;
    delete[] one;
    return;
}
/*###########################################################################################
#
#   les-FUNCTION: use id_t G(t,nh)= xxx  , t=0...n
#   note: G(t,n) is not continuous in memory!!!
#   for n<k: start routine
#
###########################################################################################*/
/** \brief <b> One step VIE solver \f$(1+F)*G=Q\f$ for a lesser component of the Green's function \f$G\f$ at a given timestep</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for a lesser component of \f$G^<(t, t^\prime)\f$ at a given timestep, for given:
* > input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$, the source term \f$Q(t, t^\prime)\f$, and
* > an integrator class 'I'. One uses use \f$ i \partial t G(t,nh)= xxx \f$  , t=0...n
*
* \note: F,G and Q are instances of the template class `GG`, representing `herm_matrix`.
* \note: G(t,n) is not continuous in memory! for n<k: start routine
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] time step
* @param &G
* > [GG] solution
* @param &F
* > [GG] green's function  on left-hand side
* @param &Fcc
* > [GG] Complex conjugate of F
* @param &Q
* > [GG] green's function  on right-hand side
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/
template <typename T, class GG, int SIZE1>
void vie2_timestep_les(int n, GG &G, GG &F, GG &Fcc, GG &Q, integration::Integrator<T> &I,
                       T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1;
    int sg, n1, l, m, j, ntau, p, q, sig, size1 = G.size1();
    cplx *gles, cweight, ih, *mm, *qq, *stemp, cplx_i = cplx(0, 1), *one, *gtemp, *qtemp;

    n1 = (n > k ? n : k);
    sg = G.element_size();
    ntau = G.ntau();
    sig = G.sig();
    // check consistency:  (more assertations follow in convolution)
    assert(F.nt() >= n1);
    assert(G.nt() >= n1);
    assert(G.sig() == F.sig());

    one = new cplx[sg];
    qq = new cplx[k1 * sg];
    mm = new cplx[k1 * k1 * sg];
    gtemp = new cplx[sg];
    qtemp = new cplx[sg];
    stemp = new cplx[sg]; // sic
    gles = new cplx[(n1 + 1) * sg];
    for (j = 0; j <= n1; j++)
        element_set_zero<T, SIZE1>(size1, gles + j * sg);
    element_set<T, SIZE1>(size1, one, 1.0);

    // CONVOLUTION SIGMA*G:  --->  G^les(j,n)
    // Note: this is only the tv*vt + les*adv part, Gles is not adressed
    // Note: gles is determioned for j=0...max(k,n)!!!
    convolution_timestep_les_tvvt<T, GG, SIZE1>(n, gles, G, F, Fcc, G, G, I, beta,
                                                h); // sig needed!!!
    convolution_timestep_les_lesadv<T, GG, SIZE1>(n, gles, G, F, Fcc, G, G, I, beta, h);
    for (j = 0; j <= n1; j++)
        element_smul<T, SIZE1>(size1, gles + j * sg, -1.0);

    // Start for integrodifferential equation: j=1...k
    // .... the usual mess:
    // derive linear equations
    // mm(p,q)*G(q)=Q(p) for p=j-1=0...k-1, q=j-1=0...k
    // G(p)=G(p,m)
    for (l = 0; l < k1 * k1 * sg; l++)
        mm[l] = 0;
    for (l = 0; l < k1 * sg; l++)
        qq[l] = 0;
    for (j = 0; j <= k; j++) {
        p = j;
        if (j <= n) {
            element_incr<T, SIZE1>(size1, qq + p * sg, Q.lesptr(j, n));
        } else {
            element_conj<T, SIZE1>(size1, qtemp, Q.lesptr(n, j));
            element_smul<T, SIZE1>(size1, qtemp, -1.0);
            element_incr<T, SIZE1>(size1, qq + p * sg, qtemp);
        }
        // integral Sigmatv*Gvt + Sigma^les*G^adv --> take from gles(j)
        // after that, gles can be overwritten
        element_incr<T, SIZE1>(size1, qq + p * sg, gles + j * sg);
        //  -- goes into m(p,p)
        element_incr<T, SIZE1>(size1, mm + sg * (p * k1 + p), one);
        // integral Sigma^ret(j,m)G^les(m,n)
        for (m = 0; m <= k; m++) {
            cweight = -h * I.gregory_weights(j, m);
            if (0 && m == 0) { // goes into qq(p)
                element_incr<T, SIZE1>(size1, qq + p * sg, cweight, F.retptr(j, 0),
                                       gles + 0 * sg);
            } else { // goes into mm(p,q)
                // q=m-1;
                q = m;
                if (j >= m) {
                    element_set<T, SIZE1>(size1, stemp, F.retptr(j, m));
                } else {
                    element_set<T, SIZE1>(size1, stemp, Fcc.retptr(m, j));
                    element_conj<T, SIZE1>(size1, stemp);
                    element_smul<T, SIZE1>(size1, stemp, -1);
                }
                for (l = 0; l < sg; l++)
                    mm[sg * (p * k1 + q) + l] += -cweight * stemp[l];
            }
        }
    }
    element_linsolve_right<T, SIZE1>(size1, k1, gles, mm,
                                     qq); // solve k1Xk1 problem mm*gtemp=qq
    // integrodifferential equation k+1...n
    for (j = k + 1; j <= n; j++) {
        element_set_zero<T, SIZE1>(size1, mm);
        // CONTRIBUTION FROM INTEGRAL tv*vt+les*adv
        element_set<T, SIZE1>(size1, qq, gles + j * sg); // now gles(j) may be overwritten
        element_incr<T, SIZE1>(size1, qq, Q.lesptr(j, n));
        element_incr<T, SIZE1>(size1, mm, one);
        // CONTRIBUTION FROM INTEGRAL Sigma^ret(j,m)*G^les(m,n)
        element_set<T, SIZE1>(size1, stemp, F.retptr(j, j));
        for (l = 0; l < sg; l++)
            mm[l] += h * stemp[l] * I.gregory_weights(j, j);
        for (m = 0; m < j; m++) {
            cweight = -h * I.gregory_weights(j, m);
            element_incr<T, SIZE1>(size1, qq, cweight, F.retptr(j, m), gles + m * sg);
        }
        element_linsolve_right<T, SIZE1>(size1, gles + j * sg, mm, qq);
    }
    // write elements into Gles
    for (j = 0; j <= n; j++)
        element_set<T, SIZE1>(size1, G.lesptr(j, n), gles + j * sg);
    delete[] stemp;
    delete[] gtemp;
    delete[] qtemp;
    delete[] gles;
    delete[] qq;
    delete[] mm;
    delete[] one;
    return;
}
/** \brief <b> VIE solver \f$(1+F)*G=Q\f$ for a lesser component of the Green's function \f$G\f$ for the first k timesteps</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for a lesser component of \f$G^(t, t^\prime)<\f$
* > for the first k timesteps (given by the integrator class 'I').
* > Here, are given: input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$,
* > and the source term \f$Q(t, t^\prime)\f$. Basically, the function calls 'vie2_timestep_les'
*
* \note: F,G and Q are instances of the template class `GG`, representing `herm_matrix`.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [GG] solution
* @param &F
* > [GG] green's function  on left-hand side
* @param &Fcc
* > [GG] Complex conjugate of F
* @param &Q
* > [GG] green's function  on right-hand side
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/
template <typename T, class GG, int SIZE1>
void vie2_start_les(GG &G, GG &F, GG &Fcc, GG &Q, integration::Integrator<T> &I, T beta,
                    T h) {
    int k = I.get_k(), n;
    for (n = 0; n <= k; n++)
        vie2_timestep_les<T, GG, SIZE1>(n, G, F, Fcc, Q, I, beta, h);
    return;
}
/** \brief <b> VIE solver \f$(1+F)*G=Q\f$ using Fourier method on matsubara axis</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$
* > on a Matsubara axis using the Fourier method.
* > Here, are given: input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$,
* > and the source term \f$Q(t, t^\prime)\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix] solution
* @param &F
* > [herm_matrix] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix] Complex conjugate of F
* @param &Q
* > [herm_matrix] green's function  on right-hand side
* @param beta
* > [double] inverse temperature
* @param order
* > [double] order
*/
template <typename T>
void vie2_mat_fourier(herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc, herm_matrix<T> &Q,
              T beta, int order) {
    int size1 = G.size1(), pcf = 20;
    assert(G.size1() == F.size1());
    assert(G.ntau() == F.ntau());
    switch (size1) {
    case 1:
        vie2_mat_fourier_dispatch<T, herm_matrix<T>, 1>(G, F, Fcc, Q, beta, pcf, order);
        break;
    case 2:
        vie2_mat_fourier_dispatch<T, herm_matrix<T>, 2>(G, F, Fcc, Q, beta, pcf, order);
        break;
    case 3:
        vie2_mat_fourier_dispatch<T, herm_matrix<T>, 3>(G, F, Fcc, Q, beta, pcf, order);
        break;
    case 4:
        vie2_mat_fourier_dispatch<T, herm_matrix<T>, 4>(G, F, Fcc, Q, beta, pcf, order);
        break;
    case 5:
        vie2_mat_fourier_dispatch<T, herm_matrix<T>, 5>(G, F, Fcc, Q, beta, pcf, order);
        break;
    case 8:
        vie2_mat_fourier_dispatch<T, herm_matrix<T>, 8>(G, F, Fcc, Q, beta, pcf, order);
        break;
    default:
        vie2_mat_fourier_dispatch<T, herm_matrix<T>, LARGESIZE>(G, F, Fcc, Q, beta, pcf, order);
        break;
    }
}
/** \brief <b> VIE solver \f$(1+F)*G=Q\f$ using fixpoint iteration method on matsubara axis</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$
* > on a Matsubara axis using the fixpoint iteration method.
* > Here, are given: input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$,
* > and the source term \f$Q(t, t^\prime)\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix] solution
* @param &F
* > [herm_matrix] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix] Complex conjugate of F
* @param &Q
* > [herm_matrix] green's function  on right-hand side
* @param beta
* > [double] inverse temperature
* @param I
* > [Integrator] integrator class
* @param nfixpoint
* > [double]
*/
template <typename T>
void vie2_mat_fixpoint(herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc,
                    herm_matrix<T> &Q, T beta, integration::Integrator<T> &I, int nfixpoint) {
    int size1 = G.size1(), pcf = 5, order = 3;
    assert(G.size1() == F.size1());
    assert(G.ntau() == F.ntau());
    switch (size1) {
    case 1:
        vie2_mat_fixpoint_dispatch<T, herm_matrix<T>, 1>(G, F, Fcc, Q, beta, I, nfixpoint, pcf,
                                                      order);
        break;
    case 2:
        vie2_mat_fixpoint_dispatch<T, herm_matrix<T>, 2>(G, F, Fcc, Q, beta, I, nfixpoint, pcf,
                                                      order);
        break;
    case 3:
        vie2_mat_fixpoint_dispatch<T, herm_matrix<T>, 3>(G, F, Fcc, Q, beta, I, nfixpoint, pcf,
                                                      order);
        break;
    case 4:
        vie2_mat_fixpoint_dispatch<T, herm_matrix<T>, 4>(G, F, Fcc, Q, beta, I, nfixpoint, pcf,
                                                      order);
        break;
    case 5:
        vie2_mat_fixpoint_dispatch<T, herm_matrix<T>, 5>(G, F, Fcc, Q, beta, I, nfixpoint, pcf,
                                                      order);
        break;
    case 8:
        vie2_mat_fixpoint_dispatch<T, herm_matrix<T>, 8>(G, F, Fcc, Q, beta, I, nfixpoint, pcf,
                                                      order);
        break;
    default:
        vie2_mat_fixpoint_dispatch<T, herm_matrix<T>, LARGESIZE>(G, F, Fcc, Q, beta, I, nfixpoint,
                                                              pcf, order);
        break;
    }
}
/** \brief <b> VIE solver \f$(1+F)*G=Q\f$ using steepest descent minimization on matsubara axis</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$
* > on a Matsubara axis using the steepest descent minimization method.
* > Here, are given: input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$,
* > and the source term \f$Q(t, t^\prime)\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix] solution
* @param &F
* > [herm_matrix] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix] Complex conjugate of F
* @param &Q
* > [herm_matrix] green's function  on right-hand side
* @param beta
* > [double] inverse temperature
* @param I
* > [Integrator] integrator class
* @param maxiter
* > [int]
* @param tol
* >
*/
template <typename T>
void vie2_mat_steep(herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc,
		 herm_matrix<T> &Q, T beta, integration::Integrator<T> &I, int maxiter, T tol) {
    int size1 = G.size1(), pcf = 3, order = 3;
    assert(G.size1() == F.size1());
    assert(G.ntau() == F.ntau());
    switch (size1) {
    case 1:
      vie2_mat_steep_dispatch<T, herm_matrix<T>, 1>(G, F, Fcc, Q, beta, I, maxiter, tol, pcf,
						 order);
      break;
    case 2:
        vie2_mat_steep_dispatch<T, herm_matrix<T>, 2>(G, F, Fcc, Q, beta, I, maxiter, tol, pcf,
                                                      order);
        break;
    case 3:
        vie2_mat_steep_dispatch<T, herm_matrix<T>, 3>(G, F, Fcc, Q, beta, I, maxiter, tol, pcf,
                                                      order);
        break;
    case 4:
        vie2_mat_steep_dispatch<T, herm_matrix<T>, 4>(G, F, Fcc, Q, beta, I, maxiter, tol, pcf,
                                                      order);
        break;
    case 5:
        vie2_mat_steep_dispatch<T, herm_matrix<T>, 5>(G, F, Fcc, Q, beta, I, maxiter, tol, pcf,
                                                      order);
        break;
    case 8:
        vie2_mat_steep_dispatch<T, herm_matrix<T>, 8>(G, F, Fcc, Q, beta, I, maxiter, tol, pcf,
                                                      order);
        break;
    default:
        vie2_mat_steep_dispatch<T, herm_matrix<T>, LARGESIZE>(G, F, Fcc, Q, beta, I, maxiter, tol,
                                                              pcf, order);
        break;
    }
}

/** \brief <b> VIE solver \f$(1+F)*G=Q\f$ for a Green's function \f$G\f$ on the Matsubara axis</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for a hermitian matrix \f$G(t, t^\prime)\f$ on a Matsubara axis, for given:
* > input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$, and the source term \f$Q(t, t^\prime)\f$.
* > There are 3 possible methods for solution: Fourier, steep, and fixpoint.
* > Fixpoint method is choosen by default.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix] solution
* @param &F
* > [herm_matrix] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix] Complex conjugate of F
* @param &Q
* > [herm_matrix] green's function  on right-hand side
* @param beta
* > [double] inverse temperature
* @param I
* > [Integrator] integrator class
* @param method
* > [const] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
*/
template <typename T>
void vie2_mat(herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc,
	 herm_matrix<T> &Q, T beta, integration::Integrator<T> &I, const int method){

  const int fourier_order=3;
  int maxiter;
  T tol = 1.0e-12;

  switch(method) {
  case CNTR_MAT_FOURIER:
    vie2_mat_fourier(G, F, Fcc, Q, beta, fourier_order);
    break;
  case CNTR_MAT_CG:
    maxiter = 40;
    vie2_mat_steep(G, F, Fcc, Q, beta, I, maxiter, tol);
    break;
  default:
    maxiter = 6;
    vie2_mat_fixpoint(G, F, Fcc, Q, beta, I, maxiter);
    break;
  }

}

/** \brief <b> VIE solver \f$(1+F)*G=Q\f$ for a Green's function \f$G\f$ on the Matsubara axis</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for a hermitian matrix \f$G(t, t^\prime)\f$ on a Matsubara axis, for given:
* > input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$, and the source term \f$Q(t, t^\prime)\f$.
* > There are 3 possible methods for solution: Fourier, steep, and fixpoint.
* > Fixpoint method is choosen by default.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix] solution
* @param &F
* > [herm_matrix] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix] Complex conjugate of F
* @param &Q
* > [herm_matrix] green's function  on right-hand side
* @param beta
* > [double] inverse temperature
* @param I
* > [Integrator] integrator class
* @param method
* > [const] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
*/
template <typename T>
void vie2_mat(T beta, herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc,
     herm_matrix<T> &Q, const int kt, const int method){

  const int fourier_order=3;
  int maxiter;
  T tol = 1.0e-12;

  switch(method) {
  case CNTR_MAT_FOURIER:
    vie2_mat_fourier(G, F, Fcc, Q, beta, fourier_order);
    break;
  case CNTR_MAT_CG:
    maxiter = 40;
    vie2_mat_steep(G, F, Fcc, Q, beta, integration::I<T>(kt), maxiter, tol);
    break;
  default:
    maxiter = 6;
    vie2_mat_fixpoint(G, F, Fcc, Q, beta, integration::I<T>(kt), maxiter);
    break;
  }

}

/** \brief <b> VIE solver \f$(1+F)*G=Q\f$ for a Green's function \f$G\f$ for the first k timesteps</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for a hermitian matrix \f$G(t, t^\prime)\f$
* > for the first k timesteps (given by the integrator class 'I').
* > Here, are given: input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$,
* > and the source term \f$Q(t, t^\prime)\f$. Basically, the function calls 'vie2_start_ret', 'vie2_start_tv', and 'vie2_start_les'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix] solution
* @param &F
* > [herm_matrix] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix] Complex conjugate of F
* @param &Q
* > [herm_matrix] green's function  on right-hand side
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/
template <typename T>
void vie2_start(herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc, herm_matrix<T> &Q,
                integration::Integrator<T> &I, T beta, T h) {
    int size1 = G.size1(), k = I.k();
    assert(G.size1() == F.size1());
    assert(G.ntau() == F.ntau());
    assert(G.nt() >= k);
    assert(F.nt() >= k);
    assert(G.size1() == Fcc.size1());
    assert(G.ntau() == Fcc.ntau());
    assert(Fcc.nt() >= k);

    switch (size1) {
    case 1:
        vie2_start_ret<T, herm_matrix<T>, 1>(G, F, Fcc, Q, I, h);
        vie2_start_tv<T, herm_matrix<T>, 1>(G, F, Fcc, Q, I, beta, h);
        vie2_start_les<T, herm_matrix<T>, 1>(G, F, Fcc, Q, I, beta, h);
        break;
    case 2:
        vie2_start_ret<T, herm_matrix<T>, 2>(G, F, Fcc, Q, I, h);
        vie2_start_tv<T, herm_matrix<T>, 2>(G, F, Fcc, Q, I, beta, h);
        vie2_start_les<T, herm_matrix<T>, 2>(G, F, Fcc, Q, I, beta, h);
        break;
    case 3:
        vie2_start_ret<T, herm_matrix<T>, 3>(G, F, Fcc, Q, I, h);
        vie2_start_tv<T, herm_matrix<T>, 3>(G, F, Fcc, Q, I, beta, h);
        vie2_start_les<T, herm_matrix<T>, 3>(G, F, Fcc, Q, I, beta, h);
        break;
    case 4:
        vie2_start_ret<T, herm_matrix<T>, 4>(G, F, Fcc, Q, I, h);
        vie2_start_tv<T, herm_matrix<T>, 4>(G, F, Fcc, Q, I, beta, h);
        vie2_start_les<T, herm_matrix<T>, 4>(G, F, Fcc, Q, I, beta, h);
        break;
    case 5:
        vie2_start_ret<T, herm_matrix<T>, 5>(G, F, Fcc, Q, I, h);
        vie2_start_tv<T, herm_matrix<T>, 5>(G, F, Fcc, Q, I, beta, h);
        vie2_start_les<T, herm_matrix<T>, 5>(G, F, Fcc, Q, I, beta, h);
        break;
    case 8:
        vie2_start_ret<T, herm_matrix<T>, 8>(G, F, Fcc, Q, I, h);
        vie2_start_tv<T, herm_matrix<T>, 8>(G, F, Fcc, Q, I, beta, h);
        vie2_start_les<T, herm_matrix<T>, 8>(G, F, Fcc, Q, I, beta, h);
        break;
    default:
        vie2_start_ret<T, herm_matrix<T>, LARGESIZE>(G, F, Fcc, Q, I, h);
        vie2_start_tv<T, herm_matrix<T>, LARGESIZE>(G, F, Fcc, Q, I, beta, h);
        vie2_start_les<T, herm_matrix<T>, LARGESIZE>(G, F, Fcc, Q, I, beta, h);
        break;
    }
}

/** \brief <b> VIE solver \f$(1+F)*G=Q\f$ for a Green's function \f$G\f$ for the first k timesteps</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for a hermitian matrix \f$G(t, t^\prime)\f$
* > for the first k timesteps (given by the integrator class 'I').
* > Here, are given: input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$,
* > and the source term \f$Q(t, t^\prime)\f$. Basically, the function calls 'vie2_start_ret', 'vie2_start_tv', and 'vie2_start_les'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix] solution
* @param &F
* > [herm_matrix] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix] Complex conjugate of F
* @param &Q
* > [herm_matrix] green's function  on right-hand side
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/
template <typename T>
void vie2_start(T beta, T h, herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc, herm_matrix<T> &Q,
                const int kt) {
    int size1 = G.size1();
    assert(G.size1() == F.size1());
    assert(G.ntau() == F.ntau());
    assert(G.nt() >= kt);
    assert(F.nt() >= kt);
    assert(G.size1() == Fcc.size1());
    assert(G.ntau() == Fcc.ntau());
    assert(Fcc.nt() >= kt);

    switch (size1) {
    case 1:
        vie2_start_ret<T, herm_matrix<T>, 1>(G, F, Fcc, Q, integration::I<T>(kt), h);
        vie2_start_tv<T, herm_matrix<T>, 1>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        vie2_start_les<T, herm_matrix<T>, 1>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        break;
    case 2:
        vie2_start_ret<T, herm_matrix<T>, 2>(G, F, Fcc, Q, integration::I<T>(kt), h);
        vie2_start_tv<T, herm_matrix<T>, 2>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        vie2_start_les<T, herm_matrix<T>, 2>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        break;
    case 3:
        vie2_start_ret<T, herm_matrix<T>, 3>(G, F, Fcc, Q, integration::I<T>(kt), h);
        vie2_start_tv<T, herm_matrix<T>, 3>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        vie2_start_les<T, herm_matrix<T>, 3>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        break;
    case 4:
        vie2_start_ret<T, herm_matrix<T>, 4>(G, F, Fcc, Q, integration::I<T>(kt), h);
        vie2_start_tv<T, herm_matrix<T>, 4>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        vie2_start_les<T, herm_matrix<T>, 4>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        break;
    case 5:
        vie2_start_ret<T, herm_matrix<T>, 5>(G, F, Fcc, Q, integration::I<T>(kt), h);
        vie2_start_tv<T, herm_matrix<T>, 5>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        vie2_start_les<T, herm_matrix<T>, 5>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        break;
    case 8:
        vie2_start_ret<T, herm_matrix<T>, 8>(G, F, Fcc, Q, integration::I<T>(kt), h);
        vie2_start_tv<T, herm_matrix<T>, 8>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        vie2_start_les<T, herm_matrix<T>, 8>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        break;
    default:
        vie2_start_ret<T, herm_matrix<T>, LARGESIZE>(G, F, Fcc, Q, integration::I<T>(kt), h);
        vie2_start_tv<T, herm_matrix<T>, LARGESIZE>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        vie2_start_les<T, herm_matrix<T>, LARGESIZE>(G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        break;
    }
}

/** \brief <b> One step VIE solver \f$(1+F)*G=Q\f$ for a Green's function \f$G\f$ at a given timestep</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for \f$G(t, t^\prime)\f$ at a given timestep, for given:
* > input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$, the source term \f$Q(t, t^\prime)\f$, and
* > the integrator class 'I'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] time step
* @param &G
* > [herm_matrix<T>] solution
* @param &F
* > [herm_matrix<T>] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix<T>] Complex conjugate of F
* @param &Q
* > [herm_matrix<T>] green's function  on right-hand side
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
* @param matsubara_method
* > [const] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
*/
template <typename T>
void vie2_timestep(int n, herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc,
                   herm_matrix<T> &Q, integration::Integrator<T> &I, T beta, T h, const int matsubara_method) {
    int size1 = G.size1(), k = I.k();
    assert(G.size1() == F.size1());
    assert(G.size1() == Fcc.size1());
    assert(G.ntau() == F.ntau());
    assert(G.ntau() == Fcc.ntau());
    assert(G.nt() >= n);
    assert(F.nt() >= n);
    assert(Fcc.nt() >= n);

    if (n==-1){
        vie2_mat(G, F, Fcc, Q, beta, I, matsubara_method);
    }else if(n<=k){
        vie2_start(G,F,Fcc,Q,integration::I<T>(n),beta,h);
    }else{
        switch (size1) {
        case 1:
            vie2_timestep_ret<T, herm_matrix<T>, 1>(n, G, Fcc, F, Q, I, h);
            vie2_timestep_tv<T, herm_matrix<T>, 1>(n, G, F, Fcc, Q, I, beta, h);
            vie2_timestep_les<T, herm_matrix<T>, 1>(n, G, F, Fcc, Q, I, beta, h);
            break;
        case 2:
            vie2_timestep_ret<T, herm_matrix<T>, 2>(n, G, Fcc, F, Q, I, h);
            vie2_timestep_tv<T, herm_matrix<T>, 2>(n, G, F, Fcc, Q, I, beta, h);
            vie2_timestep_les<T, herm_matrix<T>, 2>(n, G, F, Fcc, Q, I, beta, h);
            break;
        case 3:
            vie2_timestep_ret<T, herm_matrix<T>, 3>(n, G, Fcc, F, Q, I, h);
            vie2_timestep_tv<T, herm_matrix<T>, 3>(n, G, F, Fcc, Q, I, beta, h);
            vie2_timestep_les<T, herm_matrix<T>, 3>(n, G, F, Fcc, Q, I, beta, h);
            break;
        case 4:
            vie2_timestep_ret<T, herm_matrix<T>, 4>(n, G, Fcc, F, Q, I, h);
            vie2_timestep_tv<T, herm_matrix<T>, 4>(n, G, F, Fcc, Q, I, beta, h);
            vie2_timestep_les<T, herm_matrix<T>, 4>(n, G, F, Fcc, Q, I, beta, h);
            break;
        case 5:
            vie2_timestep_ret<T, herm_matrix<T>, 5>(n, G, Fcc, F, Q, I, h);
            vie2_timestep_tv<T, herm_matrix<T>, 5>(n, G, F, Fcc, Q, I, beta, h);
            vie2_timestep_les<T, herm_matrix<T>, 5>(n, G, F, Fcc, Q, I, beta, h);
            break;
        case 8:
            vie2_timestep_ret<T, herm_matrix<T>, 8>(n, G, Fcc, F, Q, I, h);
            vie2_timestep_tv<T, herm_matrix<T>, 8>(n, G, F, Fcc, Q, I, beta, h);
            vie2_timestep_les<T, herm_matrix<T>, 8>(n, G, F, Fcc, Q, I, beta, h);
            break;
        default:
            vie2_timestep_ret<T, herm_matrix<T>, LARGESIZE>(n, G, Fcc, F, Q, I, h);
            vie2_timestep_tv<T, herm_matrix<T>, LARGESIZE>(n, G, F, Fcc, Q, I, beta, h);
            vie2_timestep_les<T, herm_matrix<T>, LARGESIZE>(n, G, F, Fcc, Q, I, beta, h);
        break;
        }
    }
}


/** \brief <b> One step VIE solver \f$(1+F)*G=Q\f$ for a Green's function \f$G\f$ at a given timestep</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for \f$G(t, t^\prime)\f$ at a given timestep, for given:
* > input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$, the source term \f$Q(t, t^\prime)\f$, and
* > the integrator class 'I'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] time step
* @param &G
* > [herm_matrix<T>] solution
* @param &F
* > [herm_matrix<T>] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix<T>] Complex conjugate of F
* @param &Q
* > [herm_matrix<T>] green's function  on right-hand side
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
* @param matsubara_method
* > [const] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
*/
template <typename T>
void vie2_timestep(int n, T beta, T h, herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc,
                   herm_matrix<T> &Q, const int kt,  const int matsubara_method) {
    int size1 = G.size1(), k = I.k();
    assert(G.size1() == F.size1());
    assert(G.size1() == Fcc.size1());
    assert(G.ntau() == F.ntau());
    assert(G.ntau() == Fcc.ntau());
    assert(G.nt() >= n);
    assert(F.nt() >= n);
    assert(Fcc.nt() >= n);

    if (n==-1){
        vie2_mat(G, F, Fcc, Q, beta, integration::I<T>(kt), matsubara_method);
    }else if(n<=k){
        vie2_start(G,F,Fcc,Q,integration::I<T>(n),beta,h);
    }else{
        switch (size1) {
        case 1:
            vie2_timestep_ret<T, herm_matrix<T>, 1>(n, G, Fcc, F, Q, integration::I<T>(kt), h);
            vie2_timestep_tv<T, herm_matrix<T>, 1>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            vie2_timestep_les<T, herm_matrix<T>, 1>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            break;
        case 2:
            vie2_timestep_ret<T, herm_matrix<T>, 2>(n, G, Fcc, F, Q, integration::I<T>(kt), h);
            vie2_timestep_tv<T, herm_matrix<T>, 2>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            vie2_timestep_les<T, herm_matrix<T>, 2>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            break;
        case 3:
            vie2_timestep_ret<T, herm_matrix<T>, 3>(n, G, Fcc, F, Q, integration::I<T>(kt), h);
            vie2_timestep_tv<T, herm_matrix<T>, 3>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            vie2_timestep_les<T, herm_matrix<T>, 3>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            break;
        case 4:
            vie2_timestep_ret<T, herm_matrix<T>, 4>(n, G, Fcc, F, Q, integration::I<T>(kt), h);
            vie2_timestep_tv<T, herm_matrix<T>, 4>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            vie2_timestep_les<T, herm_matrix<T>, 4>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            break;
        case 5:
            vie2_timestep_ret<T, herm_matrix<T>, 5>(n, G, Fcc, F, Q, integration::I<T>(kt), h);
            vie2_timestep_tv<T, herm_matrix<T>, 5>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            vie2_timestep_les<T, herm_matrix<T>, 5>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            break;
        case 8:
            vie2_timestep_ret<T, herm_matrix<T>, 8>(n, G, Fcc, F, Q, integration::I<T>(kt), h);
            vie2_timestep_tv<T, herm_matrix<T>, 8>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            vie2_timestep_les<T, herm_matrix<T>, 8>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            break;
        default:
            vie2_timestep_ret<T, herm_matrix<T>, LARGESIZE>(n, G, Fcc, F, Q, integration::I<T>(kt), h);
            vie2_timestep_tv<T, herm_matrix<T>, LARGESIZE>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
            vie2_timestep_les<T, herm_matrix<T>, LARGESIZE>(n, G, F, Fcc, Q, integration::I<T>(kt), beta, h);
        break;
        }
    }
}

/** \brief <b> One step VIE solver \f$(1+F)*G=Q\f$ for a Green's function \f$G\f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for a hermitian matrix \f$G(t, t^\prime)\f$
* > with given \f$F(t, t^\prime)\f$ and \f$Q(t, t^\prime)\f$.
* > Here, one calls the routines 'vie2_mat()', 'vie2_start()', 'vie2_timestep'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix<T>] solution
* @param &F
* > [herm_matrix<T>] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix<T>] Complex conjugate of F
* @param &Q
* > [herm_matrix<T>] green's function  on right-hand side
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
* @param matsubara_method
* > [const] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
*/
template <typename T>
void vie2(herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc, herm_matrix<T> &Q,
          integration::Integrator<T> &I, T beta, T h, const int matsubara_method) {
    int tstp, k = I.k();
    vie2_mat(G, F, Fcc, Q, beta, I, matsubara_method);
    if (G.nt() >= 0)
        vie2_start(G, F, Fcc, Q, I, beta, h);
    for (tstp = k + 1; tstp <= G.nt(); tstp++)
        vie2_timestep(tstp, G, F, Fcc, Q, I, beta, h);
}

/** \brief <b> One step VIE solver \f$(1+F)*G=Q\f$ for a Green's function \f$G\f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the linear equation \f$(1+F)*G=Q\f$ for a hermitian matrix \f$G(t, t^\prime)\f$
* > with given \f$F(t, t^\prime)\f$ and \f$Q(t, t^\prime)\f$.
* > Here, one calls the routines 'vie2_mat()', 'vie2_start()', 'vie2_timestep'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix<T>] solution
* @param &F
* > [herm_matrix<T>] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix<T>] Complex conjugate of F
* @param &Q
* > [herm_matrix<T>] green's function  on right-hand side
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
* @param matsubara_method
* > [const] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
*/
template <typename T>
void vie2(T beta, T h, herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc, herm_matrix<T> &Q,
          const int kt, const int matsubara_method) {
    int tstp, k = I.k();
    vie2_mat(beta, G, F, Fcc, Q, kt, matsubara_method);
    if (G.nt() >= 0)
        vie2_start(beta, h, G, F, Fcc, Q, kt);
    for (tstp = k + 1; tstp <= G.nt(); tstp++)
        vie2_timestep(tstp, beta, h, G, F, Fcc, Q, kt);
}

/** \brief <b> One step VIE solver \f$(1+F)*G=Q\f$ for Green's function with instantaneous contributions for given integration order. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solvs the linear equation \f$(1+F)*G=Q\f$ for \f$G(t, t^\prime)\f$
* > with the given input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$, and the source term \f$Q(t, t^\prime)\f$
* > for the given integration order 'kt' at a given timestep.
* > The Green's functions \f$G(t, t^\prime)\f$ and \f$Q(t, t^\prime)\f$ have instantaneous contributions.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] time step
* @param &G
* > [herm_matrix<T>] solution
* @param &Gsin
* > [function<T>] singular component fo G
* @param &F
* > [herm_matrix<T>] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix<T>] Complex conjugate of F
* @param &Fsin
* > [function<T>] singular component fo F
* @param &Q
* > [herm_matrix<T>] green's function  on right-hand side
* @param &Qsin
* > [function<T>] singular component of Q
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
* @param kt
* > [int] integration order, 'I'
*/
template < typename T>
void vie2_timestep_sin(int n,herm_matrix<T> &G,function<T> &Gsin,herm_matrix<T> &F,herm_matrix<T> &Fcc, function<T> &Fsin ,herm_matrix<T> &Q,function<T> &Qsin,T beta,T h,int kt){

    int n1=(n<=kt && n>=0 ? 0 : n);
    int n2=(n<=kt && n>=0 ? kt : n);
    int nt=F.nt(),ntau=F.ntau(),size1=F.size1(),sig=F.sig();

    assert(n >= -1);
    assert(ntau > 0);
    assert(kt > 0 && kt <= 5);
    assert(kt <= 2 * ntau + 2);
    assert(ntau == Fcc.ntau());
    assert(ntau == F.ntau());
    assert(ntau == Q.ntau());
    assert(F.sig()== G.sig());
    assert(Fcc.sig()== G.sig());
    assert(Q.sig()== G.sig());
    assert(F.size1()== size1);
    assert(F.size2()== size1);
    assert(Fcc.size1()== size1);
    assert(Fcc.size2()== size1);
    assert(Q.size1()== size1);
    assert(Q.size2()== size1);
    assert(n1 <= F.nt());
    assert(n1 <= Fcc.nt());
    assert(n1 <= G.nt());
    assert(n1 <= Q.nt());

    function<T> funFinv(n2,size1);

    cntr::herm_matrix<T> tmpF(nt,ntau,size1,sig),tmpFcc(nt,ntau,size1,sig),tmpF1(nt,ntau,size1,sig),tmpQ(nt,ntau,size1,sig);
    cdmatrix cdF,cdFinv,cdQ,cdG;
    //Check consistency
    assert(G.sig()==F.sig());
    assert(G.nt()==F.nt());
    assert(G.nt()==Q.nt());

    for(int n=-1;n<=n2;n++){
      Fsin.get_value(n,cdF);
      Qsin.get_value(n,cdQ);
      cdF=cdF+cdmatrix::Identity(size1,size1);
      cdFinv=cdF.inverse();  //Expected  that these matrices are small
      cdG=cdFinv*cdQ;
      funFinv.set_value(n,cdFinv);
      Gsin.set_value(n,cdG);
    }

    tmpF=F;
    tmpFcc=Fcc;
    tmpQ=Q;

    //Set new F and Q
    for(int n=-1;n<=n2;n++){
      tmpF.left_multiply(n,funFinv,1.0);
      tmpFcc.right_multiply(n,funFinv,1.0);
      tmpQ.left_multiply(n,funFinv,1.0);
    }
    tmpF1=tmpF;

    for(int n=-1;n<=n2;n++){
      tmpF1.right_multiply(n,Gsin,1.0);
    }
    for(int i=-1;i<=n2;i++){
      tmpQ.incr_timestep(i,tmpF1,-1.0);
    }
    vie2_timestep(n,G,tmpF,tmpFcc,tmpQ,integration::I<T>(kt),beta,h);
  }


/// @private
// function calls with alpha and possibly function
#if CNTR_USE_OMP == 1
template <typename T, class GG, int SIZE1>
void vie2_timestep_omp_dispatch(int omp_num_threads, int tstp, GG &B, CPLX alpha, GG &A,
                                GG &Acc, CPLX *f0, CPLX *ft, GG &Q,
                                integration::Integrator<T> &I, T beta, T h) {
    int kt = I.get_k();
    int ntau = A.ntau();
    int size1 = A.size1();
    int sc = size1 * size1;
    bool func = (ft == NULL ? false : true);
    int n1 = (tstp == -1 || tstp > kt ? tstp : kt);
    CPLX *ftcc, *f0cc;
    herm_matrix_timestep<T> Qtstp(tstp, ntau, size1);
    // save Q
    Q.get_timestep(tstp, Qtstp);
    B.set_timestep_zero(tstp);
    if (func) {
        ftcc = new CPLX[sc * n1];
        f0cc = new CPLX[sc];
        element_conj<T, SIZE1>(size1, f0cc, f0);
        for (int n = 0; n <= n1; n++)
            element_conj<T, SIZE1>(size1, ftcc + n * sc, ft + n * sc);
    } else {
        ftcc = NULL;
        f0cc = NULL;
    }
    {
        CPLX *mtemp = new CPLX[sc];
        CPLX *one = new CPLX[sc];
        // mtemp= A.retptr(tstp,tstp)*ft(tstp)*alpha*h  [NB1 x NB1]
        if (func) {
            element_mult<T, SIZE1>(size1, mtemp, A.retptr(tstp, tstp), ft + sc * tstp);
            element_smul<T, SIZE1>(size1, mtemp, h * alpha);
        } else {
            element_set<T, SIZE1>(size1, mtemp, A.retptr(tstp, tstp));
            element_smul<T, SIZE1>(size1, mtemp, h * alpha);
        }
        // one = diag(1)  [NB1 x NB1]
        element_set<T, SIZE1>(size1, one, CPLX(1, 0));
#pragma omp parallel num_threads(omp_num_threads)
        {
            int nomp = omp_get_num_threads();
            int tid = omp_get_thread_num();
            int i, n, m;
            CPLX *mm = new CPLX[sc];
            T wt;
            std::vector<bool> mask_tv, mask_les, mask_ret;
            mask_tv = std::vector<bool>(ntau + 1, false);
            mask_ret = std::vector<bool>(tstp + 1, false);
            mask_les = std::vector<bool>(tstp + 1, false);
            for (i = 0; i <= ntau; i++)
                if (i % nomp == tid)
                    mask_tv[i] = true;
            for (i = 0; i <= tstp; i++)
                if (i % nomp == tid)
                    mask_ret[i] = true;
            for (int i = 0; i <= tstp; i++)
                if (mask_ret[i])
                    mask_les[tstp - i] = true;
            mask_les[tstp] = false; // computed separately at the end
            // retarded part:
            // Q1 = Q - alpha*A*ft*B, with Bret(tstp,n)=0
            incr_convolution_ret<T, GG, SIZE1>(tstp, mask_ret, -alpha, Q, A, Acc, ft, B, B,
                                               I, h);
            for (n = 0; n <= tstp; n++) {
                if (mask_ret[n]) {
                    wt = (tstp < kt ? I.poly_integration(n, tstp, tstp)
                                    : I.gregory_weights(tstp - n, 0));
                    element_set<T, SIZE1>(size1, mm, one);
                    element_incr<T, SIZE1>(size1, mm, wt, mtemp);
                    element_linsolve_right<T, SIZE1>(size1, 1, B.retptr(tstp, n), mm,
                                                     Q.retptr(tstp, n));
                }
            }
            // tv part:
            // Q -> Q - alpha*A*ft*B, with Btv(tstp,n)=0
            incr_convolution_tv<T, GG, SIZE1>(tstp, mask_tv, -alpha, Q, A, Acc, f0, ft, B, B,
                                              I, beta, h);
            for (m = 0; m <= ntau; m++) {
                if (mask_tv[m]) {
                    wt = I.gregory_weights(tstp, tstp);
                    element_set<T, SIZE1>(size1, mm, one);
                    element_incr<T, SIZE1>(size1, mm, wt, mtemp);
                    element_linsolve_right<T, SIZE1>(size1, 1, B.tvptr(tstp, m), mm,
                                                     Q.tvptr(tstp, m));
                }
            }
            // les part:
            // Q -> Q - conj(alpha)*B*ftcc*Acc, with Bles(n,tstp)=0, n=0...tstp-1 (Bret
            // already known!)
            mask_les[tstp] = false;
            incr_convolution_les<T, GG, SIZE1>(tstp, mask_les, -conj(alpha), Q, B, B, f0cc,
                                               ftcc, Acc, A, I, beta, h);
            for (n = 0; n < tstp; n++) {
                if (mask_les[n]) {
                    wt = I.gregory_weights(tstp, tstp);
                    element_set<T, SIZE1>(size1, mm, one);
                    element_incr<T, SIZE1>(size1, mm, wt, mtemp);
                    element_conj<T, SIZE1>(size1, mm);
                    element_linsolve_left<T, SIZE1>(size1, 1, B.lesptr(n, tstp), mm,
                                                    Q.lesptr(n, tstp));
                }
            }
        }
        // Bles(tstp,tstp) (Bles(n<tstp,tstp) etc. enters the convolution!)
        {
            int n;
            CPLX *mm = new CPLX[sc];
            std::vector<bool> mask_les;
            T wt;
            n = tstp;
            mask_les = std::vector<bool>(tstp + 1, false);
            mask_les[n] = true;
            incr_convolution_les<T, GG, SIZE1>(tstp, mask_les, -conj(alpha), Q, B, B, f0cc,
                                               ftcc, Acc, A, I, beta, h);
            wt = I.gregory_weights(tstp, tstp);
            element_set<T, SIZE1>(size1, mm, one);
            element_incr<T, SIZE1>(size1, mm, wt, mtemp);
            element_conj<T, SIZE1>(size1, mm);
            element_linsolve_left<T, SIZE1>(size1, 1, B.lesptr(n, tstp), mm,
                                            Q.lesptr(n, tstp));
            delete[] mm;
        }
        delete[] one;
        delete[] mtemp;
        if (func) {
            delete[] f0cc;
            delete[] ftcc;
        }
        // restore Q
        Q.set_timestep(tstp, Qtstp);
    }
}
// same function call as above: anyway, only for  timestep
/** \brief <b> One step VIE solver \f$(1+F)*G=Q\f$ for Green's function for given timestep. OpenMP parallelized </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > OpenMP version of 'vie2_timestep'.
* > One solvs the linear equation \f$(1+F)*G=Q\f$ for \f$G(t, t^\prime)\f$
* > with the given input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$, the source term \f$Q(t, t^\prime)\f$, and
* > the integrator class 'I'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of threads for omp parallelization
* @param tstp
* > [int] time step
* @param &G
* > [herm_matrix<T>] solution
* @param &F
* > [herm_matrix<T>] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix<T>] Complex conjugate of F
* @param &Q
* > [herm_matrix<T>] green's function  on right-hand side
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
* @param matsubara_method
* > [const] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
*/
template <typename T>
void vie2_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &G, herm_matrix<T> &F,
                       herm_matrix<T> &Fcc, herm_matrix<T> &Q, integration::Integrator<T> &I,
                       T beta, T h, const int matsubara_method) {
    int ntau = G.ntau();
    int size1 = G.size1(), kt = I.k();
    int n1 = (tstp >= kt ? tstp : kt);

    assert(tstp >= 0);
    assert(ntau > 0);
    assert(kt > 0 && kt <= 5);
    assert(kt <= 2 * ntau + 2);
    assert(ntau == Fcc.ntau());
    assert(ntau == F.ntau());
    assert(ntau == Q.ntau());
    assert(F.sig()== G.sig());
    assert(Fcc.sig()== G.sig());
    assert(Q.sig()== G.sig());
    assert(F.size1()== size1);
    assert(F.size2()== size1);
    assert(Fcc.size1()== size1);
    assert(Fcc.size2()== size1);
    assert(Q.size1()== size1);
    assert(Q.size2()== size1);
    assert(n1 <= F.nt());
    assert(n1 <= Fcc.nt());
    assert(n1 <= G.nt());
    assert(n1 <= Q.nt());

    if (tstp==-1){
        vie2_mat(G, F, Fcc, Q, beta, I, matsubara_method);
    }else if(tstp<=kt){
        cntr::vie2_start(G,F,Fcc,Q,integration::I<double>(tstp),beta,h);
    }else{
        switch (size1){
        case 1:
            vie2_timestep_omp_dispatch<T, herm_matrix<T>, 1>(
            omp_num_threads, tstp, G, CPLX(1, 0), F, Fcc, NULL, NULL, Q, I, beta, h);
	break;
        case 2:
            vie2_timestep_omp_dispatch<T, herm_matrix<T>, 2>(
            omp_num_threads, tstp, G, CPLX(1, 0), F, Fcc, NULL, NULL, Q, I, beta, h);
        break;
        case 3:
            vie2_timestep_omp_dispatch<T, herm_matrix<T>, 3>(
            omp_num_threads, tstp, G, CPLX(1, 0), F, Fcc, NULL, NULL, Q, I, beta, h);
        break;
        case 4:
            vie2_timestep_omp_dispatch<T, herm_matrix<T>, 4>(
            omp_num_threads, tstp, G, CPLX(1, 0), F, Fcc, NULL, NULL, Q, I, beta, h);
        break;
        case 5:
            vie2_timestep_omp_dispatch<T, herm_matrix<T>, 5>(
            omp_num_threads, tstp, G, CPLX(1, 0), F, Fcc, NULL, NULL, Q, I, beta, h);
        break;
        case 8:
            vie2_timestep_omp_dispatch<T, herm_matrix<T>, 8>(
            omp_num_threads, tstp, G, CPLX(1, 0), F, Fcc, NULL, NULL, Q, I, beta, h);
        break;
        }

    }
}


/** \brief <b> One step VIE solver \f$(1+F)*G=Q\f$ for Green's function with instantaneous contributions for given integration order. OpenMP parallelized </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > OpenMP version of 'vie2_timestep_sin'.
* > One solves the linear equation \f$(1+F)*G=Q\f$ for a hermitian matrix \f$G(t, t^\prime)\f$.
* > Here, are given: input kernel \f$F(t, t^\prime)\f$, its hermitian conjugate \f$F^\ddagger(t, t^\prime)\f$, the source term \f$Q(t, t^\prime)\f$,
* > and the integrator class 'I'.
* > The Green's functions \f$G(t, t^\prime)\f$ and \f$Q(t, t^\prime)\f$ have instantaneous contributions.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of threads for omp parallelization
* @param tstp
* > [int] time step
* @param &G
* > [herm_matrix<T>] solution
* @param &Gsin
* > [function<T>] singular component fo G
* @param &F
* > [herm_matrix<T>] green's function  on left-hand side
* @param &Fcc
* > [herm_matrix<T>] Complex conjugate of F
* @param &Fsin
* > [function<T>] singular component fo F
* @param &Q
* > [herm_matrix<T>] green's function  on right-hand side
* @param &Qsin
* > [function<T>] singular component of Q
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/

template < typename T>
void vie2_timestep_sin_omp(int omp_num_threads, int tstp,herm_matrix<T> &G,function<T> &Gsin,
                            herm_matrix<T> &F,herm_matrix<T> &Fcc, function<T> &Fsin ,
                            herm_matrix<T> &Q,function<T> &Qsin,integration::Integrator<T> &I,T beta,T h){

    int kt=I.k();
    int n1=(tstp<=kt && tstp>=0 ? 0 : tstp);
    int n2=(tstp<=kt && tstp>=0 ? kt : tstp);
    int nt=F.nt(),ntau=F.ntau(),size1=F.size1(),sig=F.sig();

    assert(tstp >= 0);
    assert(ntau > 0);
    assert(kt > 0 && kt <= 5);
    assert(kt <= 2 * ntau + 2);
    assert(ntau == Fcc.ntau());
    assert(ntau == F.ntau());
    assert(ntau == Q.ntau());
    assert(F.sig()== G.sig());
    assert(Fcc.sig()== G.sig());
    assert(Q.sig()== G.sig());
    assert(F.size1()== size1);
    assert(F.size2()== size1);
    assert(Fcc.size1()== size1);
    assert(Fcc.size2()== size1);
    assert(Q.size1()== size1);
    assert(Q.size2()== size1);
    assert(n1 <= F.nt());
    assert(n1 <= Fcc.nt());
    assert(n1 <= G.nt());
    assert(n1 <= Q.nt());

    function<T> funFinv(n2,size1);

    cntr::herm_matrix<T> tmpF(nt,ntau,size1,sig),tmpFcc(nt,ntau,size1,sig),tmpF1(nt,ntau,size1,sig),tmpQ(nt,ntau,size1,sig);
    cdmatrix cdF,cdFinv,cdQ,cdG;
    //Check consistency
    assert(G.sig()==F.sig());
    assert(G.nt()==F.nt());
    assert(G.nt()==Q.nt());

    for(int n=-1;n<=n2;n++){
      Fsin.get_value(n,cdF);
      Qsin.get_value(n,cdQ);
      cdF=cdF+cdmatrix::Identity(size1,size1);
      cdFinv=cdF.inverse();  //Expected  that these matrices are small
      cdG=cdFinv*cdQ;
      funFinv.set_value(n,cdFinv);
      Gsin.set_value(n,cdG);
    }

    tmpF=F;
    tmpFcc=Fcc;
    tmpQ=Q;

    //Set new F and Q
    for(int n=-1;n<=n2;n++){
      tmpF.left_multiply(n,funFinv,1.0);
      tmpFcc.right_multiply(n,funFinv,1.0);
      tmpQ.left_multiply(n,funFinv,1.0);
    }
    tmpF1=tmpF;

    for(int n=-1;n<=n2;n++){
      tmpF1.right_multiply(n,Gsin,1.0);
    }
    for(int i=-1;i<=n2;i++){
      tmpQ.incr_timestep(i,tmpF1,-1.0);
    }

    vie2_timestep_omp(omp_num_threads,tstp,G,tmpF,tmpFcc,tmpQ,I,beta,h);
  }

  /** \brief <b> One step VIE solver \f$(1+F)*G=Q\f$ for a Green's function \f$G\f$. OpenMP parallelized</b>
  *
  * <!-- ====== DOCUMENTATION ====== -->
  *
  *   \par Purpose
  * <!-- ========= -->
  *
  * > OpenMP version of 'vie2'.
  * > One solves the linear equation \f$(1+F)*G=Q\f$ for a hermitian matrix \f$G(t, t^\prime)\f$ with given \f$F(t, t^\prime)\f$ and \f$Q(t, t^\prime)\f$.
  * > Here, one calls the routines 'vie2_mat()', 'vie2_start()', 'vie2_timestep'.
  *
  *
  * <!-- ARGUMENTS
  *      ========= -->
  *
  * @param omp_num_threads
  * > [int] number of threads for omp parallelization
  * @param &G
  * > [herm_matrix<T>] solution
  * @param &F
  * > [herm_matrix<T>] green's function  on left-hand side
  * @param &Fcc
  * > [herm_matrix<T>] Complex conjugate of F
  * @param &Q
  * > [herm_matrix<T>] green's function  on right-hand side
  * @param I
  * > [Integrator] integrator class
  * @param beta
  * > [double] inverse temperature
  * @param h
  * > [double] time interval
  * @param matsubara_method
  * > [const] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
  */
template <typename T>
void vie2_omp(int omp_num_threads, herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc,
              herm_matrix<T> &Q, integration::Integrator<T> &I, T beta, T h, const int matsubara_method) {
    int tstp, k = I.k();
    vie2_mat(G, F, Fcc, Q, beta, I, matsubara_method);
    // vie2_mat(G,F,Fcc,Q,beta,3);
    if (G.nt() >= 0)
        vie2_start(G, F, Fcc, Q, I, beta, h);
    for (tstp = k + 1; tstp <= G.nt(); tstp++)
        vie2_timestep_omp(omp_num_threads, tstp, G, F, Fcc, Q, I, beta, h);
}

#endif // CNTR_USE_OMP

#undef CPLX

} // namespace cntr

#endif  // CNTR_VIE2_IMPL_H
