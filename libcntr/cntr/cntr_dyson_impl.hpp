#ifndef CNTR_DYSON_IMPL_H
#define CNTR_DYSON_IMPL_H

#include "cntr_dyson_decl.hpp"
//#include "cntr_exception.hpp"
#include "cntr_elements.hpp"
#include "cntr_function_decl.hpp"
#include "cntr_herm_matrix_decl.hpp"
#include "cntr_matsubara_impl.hpp"
#include "cntr_convolution_impl.hpp"
#include "cntr_equilibrium_decl.hpp"
#include "cntr_vie2_decl.hpp"
#include "cntr_utilities_decl.hpp"

namespace cntr {

/* #######################################################################################
#      DYSON EQUATION:
#
#          [ id/dt + mu - H(t) ] G(t,t') - [Sigma*G](t,t') = 1(t,t')
#          G(t,t')[-id/dt' + mu - H(t')] - [G*Sigma](t,t') = 1(t,t')
#
#  (*) "_timestep" computes G(t,t') at timestep n, i.e.,  G^ret(nh,t'<=nh),
#      G^les(t<=nh,nh), G^tv(nt,tau=0..beta). Timestep must be >k.
#
#  (*) For the computation of timestep n, G(t,t') and Sigma(t,t') are adressed at
#      times t,t' <= max(n,k), where k is the Integration order (see Integrator).
#      I.e., the timesteps n=0..k must be computed seperately, using the routine
#      "_start", which assumes that the Matsubara component of G and Sigma(t,t')
#      for t,t'<=k are given.
#
#  (*)  The following types are allowed for G and Sigma:
#         - G=scalar,  Sigma=scalar
#         - G=matrix,  Sigma=scalar,matrix, (sparse matrix)
#       Htype is some type that can be converted into G by G.element_set
#       The Dyson equation assumes H,G, and Sigma to be hermitian.
#
###########################################################################################*/

/*###########################################################################################
#
#   RETARDED FUNCTION:
#
###########################################################################################*/

/// @private
/** \brief <b> One step Dyson solver (integral-differential form) for the retarded component of the Green's function \f$G\f$ </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ i\frac{d}{dt} + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for the retarded component of the Green's function \f$G^R(t, t^\prime)\f$ at a given timestep.
* > Here, are given: \f$\Sigma(t, t^\prime)\f$, \f$\mu\f$, and \f$H(t)\f$.
*
* \note: \f$G\f$ and \f$\Sigma\f$ are instances of the template class `GG`, representing `herm_matrix`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] time step
* @param G
* > [GG] solution
* @param mu
* > chemical potential
* @param H
* > [complex] complex function
* @param Sigma
* > [GG] self-energy, \f$\Sigma\f$
* @param I
* > [Integrator] integrator class
* @param h
* > time interval
*/
template <typename T, class GG, int SIZE1>
void dyson_timestep_ret(int n, GG &G, T mu, std::complex<T> *H, GG &Sigma,
                        integration::Integrator<T> &I, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1;
    int ss, sg, n1, l, j, i, p, q, m, j1, j2, size1 = G.size1();
    cplx *gret, *gtemp, *mm, *qq, *qqj, *stemp, *one, cplx_i, cweight, *diffw;
    cplx *sret, w0, *hj;
    T weight;
    cplx_i = cplx(0, 1);

    ss = Sigma.element_size();
    sg = G.element_size();
    gtemp = new cplx[k * sg];
    diffw = new cplx[k1 + 1];
    qq = new cplx[(n + 1) * sg];
    one = new cplx[sg];
    mm = new cplx[k * k * sg];
    hj = new cplx[sg];
    stemp = new cplx[sg]; // sic
    element_set<T, SIZE1>(size1, one, 1);
    // check consistency:
    assert(n > k);
    assert(Sigma.nt() >= n);
    assert(G.nt() >= n);
    assert(G.sig() == Sigma.sig());
    // SET ENTRIES IN TIMESTEP TO 0
    gret = G.retptr(n, 0);
    n1 = (n + 1) * sg;
    for (i = 0; i < n1; i++)
        gret[i] = 0;
    // INITIAL VALUE t' = n
    element_set<T, SIZE1>(size1, G.retptr(n, n), -cplx_i);
    // START VALUES  t' = n-j, j = 1...k: solve a kxk problem
    for (i = 0; i < k * k * sg; i++)
        mm[i] = 0;
    for (i = 0; i < k * sg; i++)
        qq[i] = 0;
    for (j = 1; j <= k; j++) {
        p = j - 1;
        // derivatives:
        for (l = 0; l <= k; l++) {
            cweight = cplx_i / h * I.poly_differentiation(j, l);
            if (l == 0) {
                element_incr<T, SIZE1>(size1, qq + p * sg, -cweight, G.retptr(n, n));
            } else {
                q = l - 1;
                element_incr<T, SIZE1>(size1, mm + sg * (p * k + q), cweight);
            }
        }
        // H
        element_set<T, SIZE1>(size1, gtemp, H + (n - j) * sg);
        element_smul<T, SIZE1>(size1, gtemp, -1);
        for (i = 0; i < sg; i++)
            gtemp[i] += mu * one[i];
        element_incr<T, SIZE1>(size1, mm + sg * (p + k * p), gtemp);
        // integral
        for (l = 0; l <= k; l++) {
            weight = h * I.gregory_weights(j, l);
            if (n - l >= n - j) {
                element_set<T, SIZE1>(
                    size1, stemp,
                    Sigma.retptr(n - l, n - j)); // stemp is element of type G!!
            } else {
                element_set<T, SIZE1>(size1, stemp, Sigma.retptr(n - j, n - l));
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
        weight = I.gregory_omega(n - m);
        gret = G.retptr(n, m);
        sret = Sigma.retptr(m, 0);
        for (j = 0; j <= n - k2; j++) {
            element_incr<T, SIZE1>(size1, qq + j * sg, I.gregory_weights(n - j, n - m), gret,
                                   Sigma.retptr(m, j));
            sret += ss;
        }
        j1 = (n - k2 + 1 < 0 ? 0 : n - k2 + 1);
        for (j = j1; j < n - k; j++) { // start weights
            weight = I.gregory_weights(n - j, n - m);
            element_incr<T, SIZE1>(size1, qq + j * sg, I.gregory_weights(n - j, n - m), gret,
                                   Sigma.retptr(m, j));
            sret += ss;
        }
    }
    // DER REST VOM SCHUETZENFEST: t' = n-l, l = k+1 ... n:
    for (p = 0; p <= k1; p++)
        diffw[p] = I.bd_weights(p) * cplx_i / h; // use BD(k+1!!)
    w0 = h * I.gregory_omega(0);
    for (l = k + 1; l <= n; l++) {
        j = n - l;
        element_set<T, SIZE1>(size1, hj, H + j * sg);
        element_conj<T, SIZE1>(size1, hj);
        // set up mm and qqj for 1x1 problem:
        qqj = qq + j * sg;
        for (i = 0; i < sg; i++) {
            qqj[i] *= h;
            for (p = 1; p <= k1; p++)
                qqj[i] += -diffw[p] * G.retptr(n, n - l + p)[i];
        }
        element_set<T, SIZE1>(size1, stemp, Sigma.retptr(j, j));
        for (i = 0; i < sg; i++)
            mm[i] = diffw[0] * one[i] + mu * one[i] - hj[i] - w0 * stemp[i];
        element_linsolve_left<T, SIZE1>(size1, G.retptr(n, j), mm, qqj);
        // compute the contribution of Gret(n,j) to
        // int dm G(n,j)Sigma(j,j2)  for j2<j, store into qq(j2), without h
        gret = G.retptr(n, j);
        sret = Sigma.retptr(j, 0);
        for (j2 = 0; j2 < j - k; j2++) {
            element_incr<T, SIZE1>(size1, qq + j2 * sg, I.gregory_weights(n - j2, n - j),
                                   gret, sret);
            sret += ss;
        }
        j1 = (j - k < 0 ? 0 : j - k);
        for (j2 = j1; j2 < j; j2++) { // start weights
            weight = I.gregory_omega(j - j2);
            element_incr<T, SIZE1>(size1, qq + j2 * sg, I.gregory_weights(n - j2, n - j),
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
    delete[] hj;
    return;
}

/// @private
/** \brief <b> Start-up procedure for calculation of the retarded component of the Green's function \f$G\f$ from the Dyson equition (integral-differential form) </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ i\frac{d}{dt} + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for the retarded component of the Green's function \f$G^R(t, t^\prime)\f$
* > for the first k timesteps (given by the integrator class 'I').
* > Here, are given: \f$\Sigma(t, t^\prime)\f$, \f$\mu\f$, and \f$H(t)\f$.
*
* \note: \f$G\f$ and \f$\Sigma\f$ are instances of the template class `GG`, representing `herm_matrix`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param G
* > [GG] solution
* @param mu
* > chemical potential
* @param H
* > [complex] time-dependent complex function
* @param Sigma
* > [GG] self-energy, \f$\Sigma\f$
* @param I
* > [Integrator] integrator class
* @param h
* > time interval
*/
template <typename T, class GG, int SIZE1>
void dyson_start_ret(GG &G, T mu, std::complex<T> *H, GG &Sigma,
                     integration::Integrator<T> &I, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k();
    int sg, l, n, j, p, q, i, dim, size1 = G.size1();
    cplx ih, *gtemp, *mm, *qq, *stemp, minusi, *one, *hj;
    T weight;

    sg = G.element_size();
    assert(Sigma.nt() >= k);
    assert(G.nt() >= k);
    assert(G.sig() == Sigma.sig());

    // temporary storage
    qq = new cplx[k * sg];
    mm = new cplx[k * k * sg];
    one = new cplx[sg];
    hj = new cplx[sg];
    gtemp = new cplx[k * sg];
    stemp = new cplx[sg]; // sic

    element_set<T, SIZE1>(size1, one, 1.0);
    minusi = cplx(0, -1);

    // set initial values:
    for (n = 0; n <= k; n++)
        element_set<T, SIZE1>(size1, G.retptr(n, n), cplx(0, -1.0));

    for (j = 0; j < k; j++) { // determine G(n,j), n=j+1..k
        dim = k - j;
        for (i = 0; i < k * sg; i++)
            qq[i] = 0;
        for (i = 0; i < k * k * sg; i++)
            mm[i] = 0;
        // fill the matrix mm, indices p,q
        for (n = j + 1; n <= k; n++) {
            p = n - (j + 1);
            for (l = 0; l <= k; l++) {
                q = l - (j + 1);
                if (l <= j) { // this involves G which is known and goes into qq(p)
                    element_conj<T, SIZE1>(size1, gtemp, G.retptr(j, l));
                    element_smul<T, SIZE1>(size1, gtemp, -1);
                    element_set_zero<T, SIZE1>(size1, stemp);
                    element_incr<T, SIZE1>(size1, stemp, Sigma.retptr(n, l), gtemp);
                    for (i = 0; i < sg; i++)
                        qq[p * sg + i] +=
                            minusi / h * I.poly_differentiation(n, l) * gtemp[i] +
                            h * I.poly_integration(j, n, l) * stemp[i];
                } else { // this goes into mm(p,q)
                    for (i = 0; i < sg; i++)
                        mm[sg * (p * dim + q) + i] =
                            -minusi / h * I.poly_differentiation(n, l) * one[i];
                    if (l == n) {
                        element_set<T, SIZE1>(size1, hj, H + l * sg);
                        for (i = 0; i < sg; i++) {
                            mm[sg * (p * dim + q) + i] += mu * one[i] - hj[i];
                        }
                    }
                    weight = h * I.poly_integration(j, n, l);
                    if (n >= l) {
                        element_set<T, SIZE1>(size1, stemp, Sigma.retptr(n, l));
                    } else {
                        element_set<T, SIZE1>(size1, stemp, Sigma.retptr(l, n));
                        element_conj<T, SIZE1>(size1, stemp);
                        weight *= -1;
                    }
                    for (i = 0; i < sg; i++)
                        mm[sg * (p * dim + q) + i] -= weight * stemp[i];
                }
            }
        }
        // solve mm*gtemp=qq
        element_linsolve_right<T, SIZE1>(size1, dim, gtemp, mm, qq);
        for (n = j + 1; n <= k; n++) {
            p = n - (j + 1);
            //      std::cerr << n << " " << j << " " << gtemp[p] << std:: endl;
            element_set<T, SIZE1>(size1, G.retptr(n, j), gtemp + p * sg);
        }
    }
    delete[] stemp;
    delete[] qq;
    delete[] mm;
    delete[] gtemp;
    delete[] one;
    delete[] hj;
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
void dyson_mat_fourier_dispatch(GG &G, GG &Sigma, T mu, std::complex<T> *H0, T beta, int order = 3) {
    typedef std::complex<double> cplx;
    cplx *sigmadft, *sigmaiomn, *z1, *z2, *one;
    cplx *expfac, *gmat, *hj, iomn, *zinv;
    int ntau, m, r, pcf, p, m2, sg, ss, l, sig, size1 = G.size1();
    double dtau;

    assert(G.ntau() == Sigma.ntau());
    sig = G.sig();
    assert(G.sig() == Sigma.sig());
    sg = G.element_size();
    ss = Sigma.element_size();
    ntau = G.ntau();
    dtau = beta / ntau;
    if (ntau % 2 == 1) {
        std::cerr << "matsubara_inverse: ntau odd" << std::endl;
        abort();
    }
    sigmadft = new cplx[(ntau + 1) * ss];
    sigmaiomn = new cplx[ss];
    expfac = new cplx[ntau + 1];
    gmat = new cplx[(ntau + 1) * sg];
    z1 = new cplx[sg];
    z2 = new cplx[sg];
    hj = new cplx[sg];
    one = new cplx[sg];
    zinv = new cplx[sg];
    element_set<T, SIZE1>(size1, one, 1.0);
    element_set<T, SIZE1>(size1, hj, H0);
    for (l = 0; l < sg; l++)
        hj[l] -= mu * one[l];
    pcf = 10;
    m2 = ntau / 2;
    matsubara_dft<T, GG, SIZE1>(sigmadft, Sigma, sig);
    set_first_order_tail<T, SIZE1>(gmat, one, beta, sg, ntau, sig, size1);

    for (m = -m2; m <= m2 - 1; m++) {
        for (p = -pcf; p <= pcf; p++) {

            iomn = cplx(0, get_omega(m + p * ntau, beta, sig));
            matsubara_ft<T, GG, SIZE1>(sigmaiomn, m + p * ntau, Sigma, sigmadft, sig, beta,
                                       order);

            element_set<T, SIZE1>(size1, z1, sigmaiomn); // convert
            element_incr<T, SIZE1>(size1, z1, hj);       // z1=H+Sigma(iomn)

            if (sig == 1 && m + p * ntau == 0) {
                // For Bosons we need to special treat the zeroth frequency
                // as 1/iwn diverges.

                // z2 = 1./(-H - S)
                element_inverse<T, SIZE1>(size1, zinv, z1);
                element_set<T, SIZE1>(size1, z2, zinv);
                element_smul<T, SIZE1>(size1, z2, -1.0);

            } else {

                // z1 = H + S
                // zinv = 1/( iwn * ( iwn - H - S ))
                // z2 = 1/(iwn - H - S) - 1/iwn = (H + S) / (iwn * (iwn - H - S))

                for (l = 0; l < sg; l++)
                    z2[l] = iomn * one[l] - z1[l];
                element_inverse<T, SIZE1>(size1, zinv, z2);

                element_smul<T, SIZE1>(size1, zinv, 1. / iomn);
                element_mult<T, SIZE1>(size1, z2, zinv, z1);
            }

            element_smul<T, SIZE1>(size1, z2, 1 / beta);
            for (r = 0; r <= ntau; r++) {
                cplx expfac(std::exp(-get_tau(r, beta, ntau) * iomn));
                for (l = 0; l < sg; l++)
                    gmat[r * sg + l] += z2[l] * expfac;
            }
        }
    }
    for (r = 0; r <= ntau; r++) {
        element_set<T, SIZE1>(size1, G.matptr(r), gmat + r * sg);
    }

    delete[] sigmadft;
    delete[] sigmaiomn;
    delete[] expfac;
    delete[] gmat;
    delete[] z1;
    delete[] z2;
    delete[] hj;
    delete[] one;
    delete[] zinv;
    return;
}

/// @private
template <typename T, class GG, int SIZE1>
void dyson_mat_fixpoint_dispatch(GG &G, GG &Sigma, T mu, cdmatrix &H0,
                 integration::Integrator<T> &I, T beta, int fixpiter){

  int k = I.get_k();
  int ntau = G.ntau(), size1=G.size1();
  int iter;
  T hdummy=1.0;
  herm_matrix<T> G0(-1,ntau,size1,G.sig());
  herm_matrix<T> G0xSGM(-1,ntau,size1,G.sig());

  green_from_H(G0, mu, H0, beta, hdummy);

  convolution_matsubara(G0xSGM, G0, Sigma, I, beta);
  G0xSGM.smul(-1,-1);

  vie2_mat_fixpoint(G, G0xSGM, G0xSGM, G0, beta, I, fixpiter);

}

/// @private
template <typename T, class GG, int SIZE1>
void dyson_mat_fixpoint_dispatch(GG &G, GG &Sigma, T mu, cdmatrix &H0, function<T> &SigmaMF,
                 integration::Integrator<T> &I, T beta, int fixpiter){

  int k = I.get_k();
  int ntau = G.ntau(), size1=G.size1();
  T hdummy=1.0;
  herm_matrix<T> G0(-1,ntau,size1,G.sig());
  herm_matrix<T> G0xSGM(-1,ntau,size1,G.sig());
  herm_matrix<T> G0xMF(-1,ntau,size1,G.sig());

  green_from_H(G0, mu, H0, beta, hdummy);

  convolution_matsubara(G0xSGM, G0, Sigma, I, beta);
  G0xMF.set_timestep(-1,G0);
  G0xMF.right_multiply(-1,SigmaMF);
  G0xSGM.incr_timestep(-1,G0xMF);
  G0xSGM.smul(-1,-1);

  vie2_mat_fixpoint(G, G0xSGM, G0xSGM, G0, beta, I, fixpiter);

}

/// @private
template <typename T, class GG, int SIZE1>
void dyson_mat_fixpoint_dispatch(GG &G, GG &Sigma, T mu, cdmatrix &H0, cdmatrix &SigmaMF,
                 integration::Integrator<T> &I, T beta, int fixpiter){

  int k = I.get_k();
  int ntau = G.ntau(), size1=G.size1();
  T hdummy=1.0;
  herm_matrix<T> G0(-1,ntau,size1,G.sig());
  herm_matrix<T> G0xSGM(-1,ntau,size1,G.sig());
  herm_matrix<T> G0xMF(-1,ntau,size1,G.sig());
  function<T> sgmf(-1,size1);

  green_from_H(G0, mu, H0, beta, hdummy);

  convolution_matsubara(G0xSGM, G0, Sigma, I, beta);
  sgmf.set_value(-1,SigmaMF);
  G0xMF.set_timestep(-1,G0);
  G0xMF.right_multiply(-1,sgmf);
  G0xSGM.incr_timestep(-1,G0xMF);
  G0xSGM.smul(-1,-1);

  vie2_mat_fixpoint(G, G0xSGM, G0xSGM, G0, beta, I, fixpiter);

}

/// @private
template <typename T, class GG, int SIZE1>
void dyson_mat_steep_dispatch(GG &G, GG &Sigma, T mu, cdmatrix &H0,
                  integration::Integrator<T> &I, T beta, int maxiter, T tol){

  int k = I.get_k();
  int ntau = G.ntau(), size1=G.size1();
  T hdummy=1.0;
  herm_matrix<T> G0(-1,ntau,size1,G.sig());
  herm_matrix<T> G0xSGM(-1,ntau,size1,G.sig());
  herm_matrix<T> SGMxG0(-1,ntau,size1,G.sig());

  green_from_H(G0, mu, H0, beta, hdummy);

  convolution_matsubara(G0xSGM, G0, Sigma, I, beta);
  convolution_matsubara(SGMxG0, Sigma, G0, I, beta);
  G0xSGM.smul(-1,-1);
  SGMxG0.smul(-1,-1);

  vie2_mat_steep(G, G0xSGM, SGMxG0, G0, beta, I, maxiter, tol);

}

/// @private
template <typename T, class GG, int SIZE1>
void dyson_mat_steep_dispatch(GG &G, GG &Sigma, T mu, cdmatrix &H0, function<T> &SigmaMF,
                  integration::Integrator<T> &I, T beta, int maxiter, T tol){

  int k = I.get_k();
  int ntau = G.ntau(), size1=G.size1();
  T hdummy=1.0;
  herm_matrix<T> G0(-1,ntau,size1,G.sig());
  herm_matrix<T> G0xSGM(-1,ntau,size1,G.sig());
  herm_matrix<T> SGMxG0(-1,ntau,size1,G.sig());
  herm_matrix<T> G0xMF(-1,ntau,size1,G.sig());
  herm_matrix<T> MFxG0(-1,ntau,size1,G.sig());

  green_from_H(G0, mu, H0, beta, hdummy);

  G0xMF.set_timestep(-1,G0);
  G0xMF.right_multiply(-1,SigmaMF);
  MFxG0.set_timestep(-1,G0);
  MFxG0.left_multiply(-1,SigmaMF);

  convolution_matsubara(G0xSGM, G0, Sigma, I, beta);
  convolution_matsubara(SGMxG0, Sigma, G0, I, beta);
  G0xSGM.incr_timestep(-1,G0xMF);
  G0xSGM.smul(-1,-1);
  SGMxG0.incr_timestep(-1,MFxG0);
  SGMxG0.smul(-1,-1);

  vie2_mat_steep(G, G0xSGM, SGMxG0, G0, beta, I, maxiter, tol);

}

/// @private
template <typename T, class GG, int SIZE1>
void dyson_mat_steep_dispatch(GG &G, GG &Sigma, T mu, cdmatrix &H0, cdmatrix &SigmaMF,
                 integration::Integrator<T> &I, T beta, int maxiter, T tol){

  int k = I.get_k();
  int ntau = G.ntau(), size1=G.size1();
  T hdummy=1.0;
  herm_matrix<T> G0(-1,ntau,size1,G.sig());
  herm_matrix<T> G0xSGM(-1,ntau,size1,G.sig());
  herm_matrix<T> SGMxG0(-1,ntau,size1,G.sig());
  herm_matrix<T> G0xMF(-1,ntau,size1,G.sig());
  herm_matrix<T> MFxG0(-1,ntau,size1,G.sig());
  function<T> sgmf(-1,size1);

  green_from_H(G0, mu, H0, beta, hdummy);

  sgmf.set_value(-1,SigmaMF);
  G0xMF.set_timestep(-1,G0);
  G0xMF.right_multiply(-1,sgmf);
  MFxG0.set_timestep(-1,G0);
  MFxG0.left_multiply(-1,sgmf);

  convolution_matsubara(G0xSGM, G0, Sigma, I, beta);
  convolution_matsubara(SGMxG0, G0, Sigma, I, beta);
  G0xSGM.incr_timestep(-1,G0xMF);
  G0xSGM.smul(-1,-1);
  SGMxG0.incr_timestep(-1,MFxG0);
  SGMxG0.smul(-1,-1);

  vie2_mat_steep(G, G0xSGM, SGMxG0, G0, beta, I, maxiter, tol);

}


/*###########################################################################################
#
#   tv-FUNCTION:
#
###########################################################################################*/
/// @private
/** \brief <b> One step Dyson solver (integral-differential form) for the left-mixing component of the Green's function \f$G\f$ </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a left-mixing component of the Green's function \f$G^\rceil(t, t^\prime)\f$ at a given timestep.
* > Here, are given: \f$\Sigma(t, t^\prime)\f$, \f$\mu\f$, and \f$H(t)\f$.
*
* \note: $G$ and \f$\Sigma\f$ are instances of the template class `GG`, representing `herm_matrix`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] time step
* @param &G
* > [GG] solution
* @param mu
* > [T] chemical potential
* @param &Hn
* > [complex<T>] complex function at a time step 'n'
* @param &Sigma
* > [GG] self-energy
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/
template <typename T, class GG, int SIZE1>
void dyson_timestep_tv(int n, GG &G, T mu, std::complex<T> *Hn, GG &Sigma,
                       integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k();
    int sg, n1, l, m, j, ntau, size1 = G.size1();
    cplx *gtv, *gtv1, cweight, ih, *mm, *qq, *stemp, minusi, *one, *htemp;
    T weight;

    sg = G.element_size();
    ntau = G.ntau();
    // check consistency:  (more assertations follow in convolution)
    assert(n > k);
    assert(Sigma.nt() >= n);
    assert(G.nt() >= n);
    assert(G.sig() == Sigma.sig());

    one = new cplx[sg];
    qq = new cplx[sg];
    mm = new cplx[sg];
    stemp = new cplx[sg]; // sic
    htemp = new cplx[sg];
    element_set<T, SIZE1>(size1, one, 1.0);
    // SET ENTRIES IN TIMESTEP(TV) TO 0
    gtv = G.tvptr(n, 0);
    n1 = (ntau + 1) * sg;
    for (l = 0; l < n1; l++)
        gtv[l] = 0;
    // CONVOLUTION SIGMA*G:  --->  Gtv(n,m)
    convolution_timestep_tv<T, herm_matrix<T>, SIZE1>(n, G, Sigma, Sigma, G, G, I, beta,
                                                      h); // note: this sets only tv
    // ACCUMULATE CONTRIBUTION TO id/dt G(t,t') FROM t=mh, m=n-k..n-1
    ih = cplx(0, 1 / h);
    for (m = n - k - 1; m < n; m++) {
        cweight = ih * I.bd_weights(n - m); // use BD(k+1!!)
        // G(n,j) -= cweight*G(m,j), for j=0...m
        n1 = (ntau + 1) * sg;
        gtv = G.tvptr(m, 0);
        gtv1 = G.tvptr(n, 0);
        for (l = 0; l < n1; l++)
            gtv1[l] -= cweight * gtv[l];
    }
    // Now solve
    // [ i/h bd(0) - H - h w(n,0) Sigma(n,n) ] G(n,m)  = Q(m),
    // where Q is initially stored in G(n,m)
    element_set<T, SIZE1>(size1, stemp, Sigma.retptr(n, n));
    weight = -h * I.gregory_weights(n, 0);
    element_set<T, SIZE1>(size1, htemp, Hn);
    for (l = 0; l < sg; l++)
        mm[l] = ih * I.bd_weights(0) * one[l] + weight * stemp[l] + mu * one[l] - htemp[l];
    for (j = 0; j <= ntau; j++) {
        for (l = 0; l < sg; l++)
            qq[l] = G.tvptr(n, j)[l];
        element_linsolve_right<T, SIZE1>(size1, G.tvptr(n, j), mm, qq);
    }
    delete[] stemp;
    delete[] htemp;
    delete[] qq;
    delete[] mm;
    delete[] one;
    return;
}
/// @private
/** \brief <b> Start-up procedure for calculation of the left-mixing component of the Green's function \f$G\f$ from the Dyson equition (integral-differential form) </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for the left-mixing component of the Green's function \f$G^\rceil(t, t^\prime)\f$
* > for the first k timesteps (given by the integrator class 'I').
* > Here, are given: \f$\Sigma(t, t^\prime)\f$, \f$\mu\f$, and \f$H(t)\f$.
*
* \note: $G$ and \f$\Sigma\f$ are instances of the template class `GG`, representing `herm_matrix`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [GG] solution
* @param mu
* > [T] chemical potential
* @param &H
* > [complex<T>] time-dependent complex function
* @param &Sigma
* > [GG] self-energy
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/
template <typename T, class GG, int SIZE1>
void dyson_start_tv(GG &G, T mu, std::complex<T> *H, GG &Sigma,
                    integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k();
    int sg, l, m, j, ntau, p, q, n, sig, size1 = G.size1();
    cplx cweight, *mm, *qq, *stemp, *one, *gtemp;
    cplx cplx_i = cplx(0.0, 1.0);
    T dtau;
    sg = G.element_size();
    ntau = G.ntau();
    dtau = beta / ntau;
    sig = G.sig();
    // check consistency:  (more assertations follow in convolution)
    assert(Sigma.nt() >= k);
    assert(G.nt() >= k);
    assert(Sigma.ntau() == ntau);
    assert(G.sig() == Sigma.sig());

    one = new cplx[sg];
    qq = new cplx[k * sg];
    mm = new cplx[k * k * sg];
    stemp = new cplx[sg]; // sic
    gtemp = new cplx[k * sg];
    element_set<T, SIZE1>(size1, one, 1.0);
    // INITIAL VALUE: G^tv(0,tau) = i sgn G^mat(beta-tau)  (sgn=Bose/Fermi)
    for (m = 0; m <= ntau; m++)
        for (l = 0; l < sg; l++)
            G.tvptr(0, m)[l] = ((T)sig) * cplx_i * G.matptr(ntau - m)[l];
    // CONVOLUTION  -i int dtau Sigma^tv(n,tau)G^mat(tau,m) ---> Gtv(n,m)
    for (n = 1; n <= k; n++) {
        for (m = 0; m <= ntau; m++) {
            matsubara_integral_2<T, SIZE1>(size1, m, ntau, gtemp, Sigma.tvptr(n, 0),
                                           G.matptr(0), I, G.sig());
            for (l = 0; l < sg; l++)
                G.tvptr(n, m)[l] = dtau * gtemp[l];
        }
    }
    // loop over m --> determine G^tv(n,m) for n=1...k
    for (m = 0; m <= ntau; m++) {
        for (l = 0; l < k * k * sg; l++)
            mm[l] = 0;
        for (l = 0; l < k * sg; l++)
            qq[l] = 0;
        // derive linear equations
        // mm(p,q)*G(q)=Q(p) for p=n-1=0...k-1, q=n-1=0...k
        // G(p)=G(p+1,m)
        for (n = 1; n <= k; n++) {
            p = n - 1;
            // derivative id/dt Gtv(n,m)
            for (j = 0; j <= k; j++) {
                cweight = cplx_i / h * I.poly_differentiation(n, j);
                if (j == 0) {
                    for (l = 0; l < sg; l++)
                        qq[p * sg + l] -= cweight * G.tvptr(0, m)[l];
                } else { // goes into mm(p,q)
                    q = j - 1;
                    for (l = 0; l < sg; l++)
                        mm[sg * (p * k + q) + l] += cweight * one[l];
                }
            }
            // H -- goes into m(p,p)
            element_set<T, SIZE1>(size1, gtemp, H + sg * n);
            element_smul<T, SIZE1>(size1, gtemp, -1.0);
            for (l = 0; l < sg; l++)
                gtemp[l] += mu * one[l];
            element_incr<T, SIZE1>(size1, mm + sg * (p * k + p), gtemp);
            // integral 0..n
            for (j = 0; j <= k; j++) {
                cweight = h * I.gregory_weights(n, j);
                if (j == 0) { // goes into qq(p)
                    element_incr<T, SIZE1>(size1, qq + p * sg, cweight, Sigma.retptr(n, 0),
                                           G.tvptr(0, m));
                } else { // goes into mm(p,q)
                    q = j - 1;
                    if (n >= j) {
                        element_set<T, SIZE1>(size1, stemp, Sigma.retptr(n, j));
                    } else {
                        element_set<T, SIZE1>(size1, stemp, Sigma.retptr(j, n));
                        element_conj<T, SIZE1>(size1, stemp);
                        element_smul<T, SIZE1>(size1, stemp, -1);
                    }
                    for (l = 0; l < sg; l++)
                        mm[sg * (p * k + q) + l] += -cweight * stemp[l];
                }
            }
            // integral Sigmatv*Gmat --> take from Gtv(n,m), write into qq
            element_incr<T, SIZE1>(size1, qq + p * sg, G.tvptr(n, m));
        }
        element_linsolve_right<T, SIZE1>(size1, k, gtemp, mm,
                                         qq); // solve kXk problem mm*gtemp=qq
        // write elements into Gtv
        for (n = 1; n <= k; n++)
            element_set<T, SIZE1>(size1, G.tvptr(n, m), gtemp + (n - 1) * sg);
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
/// @private
/** \brief <b> One step Dyson solver (integral-differential form) for the lesser component of the Green's function \f$G\f$ </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a lesser component of \f$G(t, t^\prime)\f$ at a given timestep.
* > Here, are given: \f$\Sigma(t, t^\prime)\f$, \f$\mu\f$, and \f$H(t)\f$.
*
* \note: $G$ and \f$\Sigma\f$ are instances of the template class `GG`, representing `herm_matrix`.
*
* \note: \f$G\f$(t,n) is not continuous in memory!
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] time step
* @param &G
* > [GG] solution
* @param mu
* > [T] chemical potential
* @param &H
* > [complex<T>] time-dependent complex function
* @param &Sigma
* > [GG] self-energy
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/
template <typename T, class GG, int SIZE1>
void dyson_timestep_les(int n, GG &G, T mu, std::complex<T> *H, GG &Sigma,
                        integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1;
    int sg, n1, l, m, j, ntau, p, q, sig, size1 = G.size1();
    cplx *gles, cweight, ih, *mm, *qq, *stemp, cplx_i = cplx(0, 1), *one, *gtemp;

    n1 = (n > k ? n : k);
    sg = G.element_size();
    ntau = G.ntau();
    sig = G.sig();
    // check consistency:  (more assertations follow in convolution)
    assert(Sigma.nt() >= n1);
    assert(G.nt() >= n1);
    assert(G.sig() == Sigma.sig());

    one = new cplx[sg];
    qq = new cplx[k * sg];
    mm = new cplx[k * k * sg];
    gtemp = new cplx[sg];
    stemp = new cplx[sg]; // sic
    gles = new cplx[(n1 + 1) * sg];
    for (j = 0; j <= n1; j++)
        element_set_zero<T, SIZE1>(size1, gles + j * sg);
    element_set<T, SIZE1>(size1, one, 1.0);

    // CONVOLUTION SIGMA*G:  --->  G^les(j,n)
    // Note: this is only the tv*vt + les*adv part, Gles is not adressed
    // Note: gles is determioned for j=0...max(k,n)!!!
    convolution_timestep_les_tvvt<T, GG, SIZE1>(n, gles, G, Sigma, Sigma, G, G, I, beta,
                                                h); // sig needed!!!
    convolution_timestep_les_lesadv<T, GG, SIZE1>(n, gles, G, Sigma, Sigma, G, G, I, beta,
                                                  h);
    // INITIAL VALUE  G^les(0,n) = -G^tv(n,0)^*
    element_conj<T, SIZE1>(size1, gles, G.tvptr(n, 0));
    element_smul<T, SIZE1>(size1, gles, -1);
    // Start for integrodifferential equation: j=1...k
    // .... the usual mess:
    // derive linear equations
    // mm(p,q)*G(q)=Q(p) for p=j-1=0...k-1, q=j-1=0...k
    // G(p)=G(p+1,m)
    for (l = 0; l < k * k * sg; l++)
        mm[l] = 0;
    for (l = 0; l < k * sg; l++)
        qq[l] = 0;
    for (j = 1; j <= k; j++) {
        p = j - 1;
        // integral Sigmatv*Gvt + Sigma^les*G^adv --> take from gles(j)
        // after that, gles can be overwritten
        element_incr<T, SIZE1>(size1, qq + p * sg, gles + j * sg);
        // derivative id/dt Gtv(n,m)
        for (m = 0; m <= k; m++) {
            cweight = cplx_i / h * I.poly_differentiation(j, m);
            if (m == 0) {
                for (l = 0; l < sg; l++)
                    qq[p * sg + l] -= cweight * gles[0 * sg + l];
            } else {
                q = m - 1;
                for (l = 0; l < sg; l++)
                    mm[sg * (p * k + q) + l] += cweight * one[l];
            }
        }
        // H -- goes into m(p,p)
        element_set<T, SIZE1>(size1, gtemp, H + sg * j);
        element_smul<T, SIZE1>(size1, gtemp, -1);
        for (l = 0; l < sg; l++)
            gtemp[l] += mu * one[l];
        element_incr<T, SIZE1>(size1, mm + sg * (p * k + p), gtemp);
        // integral Sigma^ret(j,m)G^les(m,n)
        for (m = 0; m <= k; m++) {
            cweight = h * I.gregory_weights(j, m);
            if (m == 0) { // goes into qq(p)
                element_incr<T, SIZE1>(size1, qq + p * sg, cweight, Sigma.retptr(j, 0),
                                       gles + 0 * sg);
            } else { // goes into mm(p,q)
                q = m - 1;
                if (j >= m) {
                    element_set<T, SIZE1>(size1, stemp, Sigma.retptr(j, m));
                } else {
                    element_set<T, SIZE1>(size1, stemp, Sigma.retptr(m, j));
                    element_conj<T, SIZE1>(size1, stemp);
                    element_smul<T, SIZE1>(size1, stemp, -1);
                }
                for (l = 0; l < sg; l++)
                    mm[sg * (p * k + q) + l] += -cweight * stemp[l];
            }
        }
    }
    element_linsolve_right<T, SIZE1>(size1, k, gles + 1 * sg, mm,
                                     qq); // solve kXk problem mm*gtemp=qq
    // integrodifferential equation k+1...n
    for (j = k + 1; j <= n; j++) {
        element_set_zero<T, SIZE1>(size1, mm);
        // CONTRIBUTION FROM INTEGRAL tv*vt+les*adv
        element_set<T, SIZE1>(size1, qq, gles + j * sg); // now gles(j) may be overwritten
        // ACCUMULATE CONTRIBUTION TO id/dt G(j-p,n) p=1...k1 into qq
        for (p = 1; p <= k1; p++) { // use BD(k+1) !!!
            cweight = -cplx_i / h * I.bd_weights(p);
            for (l = 0; l < sg; l++)
                qq[l] += cweight * gles[sg * (j - p) + l];
        }
        for (l = 0; l < sg; l++)
            mm[l] += cplx_i / h * I.bd_weights(0) * one[l];
        // CONTRIBUTION FROM H
        element_set<T, SIZE1>(size1, gtemp, H + sg * j);
        element_smul<T, SIZE1>(size1, gtemp, -1);
        for (l = 0; l < sg; l++)
            gtemp[l] += mu * one[l];
        element_incr<T, SIZE1>(size1, mm, gtemp);
        // CONTRIBUTION FROM INTEGRAL Sigma^ret(j,m)*G^les(m,j)
        element_set<T, SIZE1>(size1, stemp, Sigma.retptr(j, j));
        for (l = 0; l < sg; l++)
            mm[l] += -h * stemp[l] * I.gregory_weights(j, j);
        for (m = 0; m < j; m++) {
            cweight = h * I.gregory_weights(j, m);
            element_incr<T, SIZE1>(size1, qq, cweight, Sigma.retptr(j, m), gles + m * sg);
        }
        element_linsolve_right<T, SIZE1>(size1, gles + j * sg, mm, qq);
    }
    // write elements into Gles
    for (j = 0; j <= n; j++)
        element_set<T, SIZE1>(size1, G.lesptr(j, n), gles + j * sg);
    delete[] stemp;
    delete[] gtemp;
    delete[] gles;
    delete[] qq;
    delete[] mm;
    delete[] one;
    return;
}

/// @private
/** \brief <b> Start-up procedure for calculation of the lesser component of the Green's function \f$G\f$ from the Dyson equition (integral-differential form) </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a lesser component of \f$G(t, t^\prime)\f$
* > for the first k timesteps (given by the integrator class 'I').
* > Here, are given: \f$\Sigma(t, t^\prime)\f$, \f$\mu\f$, and \f$H(t)\f$.
*
* \note: $G$ and \f$\Sigma\f$ are instances of the template class `GG`, representing `herm_matrix`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [GG] solution
* @param mu
* > [T] chemical potential
* @param &H
* > [complex<T>] time-dependent complex function
* @param &Sigma
* > [GG] self-energy
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/
template <typename T, class GG, int SIZE1>
void dyson_start_les(GG &G, T mu, std::complex<T> *H, GG &Sigma,
                     integration::Integrator<T> &I, T beta, T h) {
    int k = I.get_k(), n;
    for (n = 0; n <= k; n++)
        dyson_timestep_les<T, GG, SIZE1>(n, G, mu, H, Sigma, I, beta, h);
    return;
}
/// @private
// preferred: H passed as a cntr::function object:
template <typename T>
void dyson_mat_fourier(herm_matrix<T> &G, herm_matrix<T> &Sigma, T mu, function<T> &H, T beta,
               int order) {
    int size1 = G.size1();
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    if (size1 == 1)
      dyson_mat_fourier_dispatch<T, herm_matrix<T>, 1>(G, Sigma, mu, H.ptr(-1), beta, order);
    else
      dyson_mat_fourier_dispatch<T, herm_matrix<T>, LARGESIZE>(G, Sigma, mu, H.ptr(-1), beta,
                                                         order);
}
/// @private
template <typename T>
void dyson_mat_fourier(herm_matrix<T> &G, herm_matrix<T> &Sigma, T mu, function<T> &H,
               function<T> &SigmaMF, T beta, int order) {
    int size1 = G.size1();
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    std::complex<T> *hmf;
    hmf = new std::complex<T>[size1*size1];

    if (size1 == 1){
      element_set<T, 1>(size1, hmf, H.ptr(-1));
      element_incr<T, 1>(size1, hmf, 1.0, SigmaMF.ptr(-1));
      dyson_mat_fourier_dispatch<T, herm_matrix<T>, 1>(G, Sigma, mu, hmf, beta, order);
    } else {
      element_set<T, LARGESIZE>(size1, hmf, H.ptr(-1));
      element_incr<T, LARGESIZE>(size1, hmf, 1.0, SigmaMF.ptr(-1));
      dyson_mat_fourier_dispatch<T, herm_matrix<T>, LARGESIZE>(G, Sigma, mu, hmf, beta,
                                                         order);
    }
    delete hmf;
}
/// @private
template <typename T>
void dyson_mat_fixpoint(herm_matrix<T> &G, herm_matrix<T> &Sigma, T mu, function<T> &H,
            integration::Integrator<T> &I,T beta, int fixpiter) {
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    int size1 = G.size1();
    cdmatrix h0(size1,size1);

    H.get_value(-1,h0);
    if (size1 == 1)
      dyson_mat_fixpoint_dispatch<T, herm_matrix<T>, 1>(G, Sigma, mu, h0, I, beta, fixpiter);
    else
      dyson_mat_fixpoint_dispatch<T, herm_matrix<T>, LARGESIZE>(G, Sigma, mu, h0, I, beta,
                                fixpiter);
}
/// @private
template <typename T>
void dyson_mat_fixpoint(herm_matrix<T> &G, herm_matrix<T> &Sigma, T mu, function<T> &H,
            function<T> &SigmaMF, integration::Integrator<T> &I,T beta, int fixpiter) {
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    int size1 = G.size1();
    cdmatrix h0(size1,size1);

    H.get_value(-1,h0);
    if (size1 == 1)
      dyson_mat_fixpoint_dispatch<T, herm_matrix<T>, 1>(G, Sigma, mu, h0, SigmaMF, I, beta, fixpiter);
    else
      dyson_mat_fixpoint_dispatch<T, herm_matrix<T>, LARGESIZE>(G, Sigma, mu, h0, SigmaMF, I, beta,
                                fixpiter);
}


/// @private
template <typename T>
void dyson_mat_steep(herm_matrix<T> &G, herm_matrix<T> &Sigma, T mu, function<T> &H,
             integration::Integrator<T> &I,T beta, int maxiter, T tol) {
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert(H.size1() == G.size1() );
    int size1 = G.size1();
    cdmatrix h0(size1,size1);

    H.get_value(-1,h0);
    if (size1 == 1)
      dyson_mat_steep_dispatch<T, herm_matrix<T>, 1>(G, Sigma, mu, h0, I, beta, maxiter, tol);
    else
      dyson_mat_steep_dispatch<T, herm_matrix<T>, LARGESIZE>(G, Sigma, mu, h0, I, beta,
                                 maxiter, tol);
}

/// @private
template <typename T>
void dyson_mat_steep(herm_matrix<T> &G, herm_matrix<T> &Sigma, T mu, function<T> &H,
             function<T> &SigmaMF, integration::Integrator<T> &I,T beta, int maxiter, T tol) {
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert(H.size1() == G.size1() );
    int size1 = G.size1();
    cdmatrix h0(size1,size1);

    H.get_value(-1,h0);
    if (size1 == 1)
      dyson_mat_steep_dispatch<T, herm_matrix<T>, 1>(G, Sigma, mu, h0, SigmaMF, I, beta, maxiter, tol);
    else
      dyson_mat_steep_dispatch<T, herm_matrix<T>, LARGESIZE>(G, Sigma, mu, h0, SigmaMF, I, beta,
                                 maxiter, tol);
}

/// @private
/** \brief <b> Dyson solver (integral-differential form) for a Green's function \f$G\f$. Global interface</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a hermitian matrix \f$G(t, t^\prime)\f$ on a Matsubara axis.
* > There are 3 possible methods for solution: Fourier, steep, and fixpoint.
* > Fixpoint method is choosen by default.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix<T>] solution
* @param &Sigma
* > [herm_matrix<T>] self-energy
* @param mu
* > [T] chemical potential
* @param &H
* > [function<T>] time-dependent function
* @param &SigmaMF
* > [herm_matrix<T>] mean-field self-energy
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param method
* > [const] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
* @param force_hermitian
* > [const bool] force hermitian solution, if 'true'
*/
template <typename T>
void dyson_mat(herm_matrix<T> &G, herm_matrix<T> &Sigma, T mu, function<T> &H,
           function<T> &SigmaMF, integration::Integrator<T> &I,T beta, const int method,
         const bool force_hermitian){
  assert(method <= 2 && "UNKNOWN CNTR_MAT_METHOD");


  const int fourier_order = 3;
  const double tol=1.0e-12;
  int maxiter;
  switch(method){
  case CNTR_MAT_FOURIER:
    dyson_mat_fourier(G, Sigma, mu, H, SigmaMF, beta, fourier_order);
    break;
  case CNTR_MAT_CG:
    maxiter = 40;
    dyson_mat_steep(G, Sigma, mu, H, SigmaMF, I, beta, maxiter, tol);
    break;
  default:
    maxiter = 6;
    dyson_mat_fixpoint(G, Sigma, mu, H, SigmaMF, I, beta, maxiter);
    break;
  }
  if(force_hermitian){
    force_matsubara_hermitian(G);
  }
}


// global interface
/** \brief <b> Dyson solver (integral-differential form) for a Green's function \f$G\f$. Global interface</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a hermitian matrix \f$G(t, t^\prime)\f$ on a Matsubara axis.
* > There are 3 possible methods for solution: Fourier, steep, and fixpoint.
* > Fixpoint method is choosen by default.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix<T>] solution
* @param mu
* > [T] chemical potential
* @param &H
* > [function<T>] time-dependent function
* @param &SigmaMF
* > [herm_matrix<T>] mean-field self-energy
* @param &Sigma
* > [herm_matrix<T>] self-energy
* @param beta
* > [double] inverse temperature
* @param SolveOrder
* > [int] integrator order
* @param method
* > [const] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
* @param force_hermitian
* > [const bool] force hermitian solution, if 'true'
*/
template <typename T>
void dyson_mat(herm_matrix<T> &G, T mu, function<T> &H,
           function<T> &SigmaMF, herm_matrix<T> &Sigma, T beta, const int SolveOrder,const int method,
         const bool force_hermitian){
  assert(method <= 2 && "UNKNOWN CNTR_MAT_METHOD");
  assert(SolveOrder <= MAX_SOLVE_ORDER);

  const int fourier_order = 3;
  const double tol=1.0e-12;
  int maxiter;
  switch(method){
  case CNTR_MAT_FOURIER:
    dyson_mat_fourier(G, Sigma, mu, H, SigmaMF, beta, fourier_order);
    break;
  case CNTR_MAT_CG:
    maxiter = 40;
    dyson_mat_steep(G, Sigma, mu, H, SigmaMF, integration::I<T>(SolveOrder), beta, maxiter, tol);
    break;
  default:
    maxiter = 6;
    dyson_mat_fixpoint(G, Sigma, mu, H, SigmaMF, integration::I<T>(SolveOrder), beta, maxiter);
    break;
  }
  if(force_hermitian){
    force_matsubara_hermitian(G);
  }
}
/// @private
/** \brief <b> Dyson solver (integral-differential form) for a Green's function \f$G\f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a hermitian matrix \f$G(t, t^\prime)\f$ on a Matsubara axis.
* > There are 3 possible methods for solution: Fourier, steep, and fixpoint.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix<T>] solution
* @param &Sigma
* > [herm_matrix<T>] self-energy
* @param mu
* > [T] chemical potential
* @param &H
* > [function<T>] time-dependent function
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param method
* > [const] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
* @param force_hermitian
* > [const bool] force hermitian solution, if 'true'
*/
template <typename T>
void dyson_mat(herm_matrix<T> &G, herm_matrix<T> &Sigma, T mu, function<T> &H,
           integration::Integrator<T> &I,T beta, const int method,const bool force_hermitian){
  assert(method <= 2 && "UNKNOWN CNTR_MAT_METHOD");

  const int fourier_order = 3;
  const double tol=1.0e-12;
  int maxiter;

  switch(method){
  case 0:
    dyson_mat_fourier(G, Sigma, mu, H, beta, fourier_order);
    break;
  case 1:
    maxiter = 40;
    dyson_mat_steep(G, Sigma, mu, H, I, beta, maxiter, tol);
    break;
  case 2:
    maxiter=6;
    dyson_mat_fixpoint(G, Sigma, mu, H, I, beta, maxiter);
    break;
  }
  if(force_hermitian){
    force_matsubara_hermitian(G);
  }
}

/** \brief <b> Dyson solver (integral-differential form) for a Green's function \f$G\f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a hermitian matrix \f$G(t, t^\prime)\f$ on a Matsubara axis.
* > There are 3 possible methods for solution: Fourier, steep, and fixpoint.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix<T>] solution
* @param mu
* > [T] chemical potential
* @param &Sigma
* > [herm_matrix<T>] self-energy
* @param &H
* > [function<T>] time-dependent function
* @param SolveOrder
* > [int] integrator order
* @param beta
* > [double] inverse temperature
* @param method
* > [int] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
* @param force_hermitian
* > [bool] force hermitian solution, if 'true'
*/
template <typename T>
void dyson_mat(herm_matrix<T> &G, T mu, function<T> &H, herm_matrix<T> &Sigma,
           T beta, const int SolveOrder, const int method,const bool force_hermitian){
  assert(method <= 2 && "UNKNOWN CNTR_MAT_METHOD");
  assert(SolveOrder <= MAX_SOLVE_ORDER);

  const int fourier_order = 3;
  const double tol=1.0e-12;
  int maxiter;

  switch(method){
  case 0:
    dyson_mat_fourier(G, Sigma, mu, H, beta, fourier_order);
    break;
  case 1:
    maxiter = 40;
    dyson_mat_steep(G, Sigma, mu, H, integration::I<T>(SolveOrder), beta, maxiter, tol);
    break;
  case 2:
    maxiter=6;
    dyson_mat_fixpoint(G, Sigma, mu, H, integration::I<T>(SolveOrder), beta, maxiter);
    break;
  }
  if(force_hermitian){
    force_matsubara_hermitian(G);
  }
}
/// @private
/** \brief <b> Start-up procedure for solving the Dyson equation of the integral-differential form for a Green's function \f$G\f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a hermitian matrix \f$G(t, t^\prime)\f$ for the first k timesteps (given by the integrator class 'I').
* > One assumes that the Matsubara component of \f$G\f$ and \f$\Sigma(t,t^\prime)\f$
* > for \f$t,t^\prime\f$<=k are given.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix<T>] solution
* @param mu
* > [T] chemical potential
* @param &H
* > [function<T>] time-dependent function
* @param &Sigma
* > [herm_matrix<T>] self-energy
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/
template <typename T>
void dyson_start(herm_matrix<T> &G, T mu, function<T> &H, herm_matrix<T> &Sigma,
                 integration::Integrator<T> &I, T beta, T h) {
    int size1 = G.size1(), k = I.k();
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert(G.nt() >= k);
    assert(Sigma.nt() >= k);
    if (size1 == 1) {
        dyson_start_ret<T, herm_matrix<T>, 1>(G, mu, H.ptr(0), Sigma, I, h);
        dyson_start_tv<T, herm_matrix<T>, 1>(G, mu, H.ptr(0), Sigma, I, beta, h);
        dyson_start_les<T, herm_matrix<T>, 1>(G, mu, H.ptr(0), Sigma, I, beta, h);
    } else {
        dyson_start_ret<T, herm_matrix<T>, LARGESIZE>(G, mu, H.ptr(0), Sigma, I, h);
        dyson_start_tv<T, herm_matrix<T>, LARGESIZE>(G, mu, H.ptr(0), Sigma, I, beta, h);
        dyson_start_les<T, herm_matrix<T>, LARGESIZE>(G, mu, H.ptr(0), Sigma, I, beta, h);
    }
}

/** \brief <b> Start-up procedure for solving the Dyson equation of the integral-differential form for a Green's function \f$G\f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a hermitian matrix \f$G(t, t^\prime)\f$ for the first k timesteps (given by the integrator class 'I').
* > One assumes that the Matsubara component of \f$G\f$ and \f$\Sigma(t,t^\prime)\f$
* > for \f$t,t^\prime\f$<=k are given.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix<T>] solution
* @param mu
* > [T] chemical potential
* @param &H
* > [function<T>] time-dependent function
* @param &Sigma
* > [herm_matrix<T>] self-energy
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T>
void dyson_start(herm_matrix<T> &G, T mu, function<T> &H, herm_matrix<T> &Sigma,
                 T beta, T h, const int SolveOrder) {
    int size1 = G.size1();
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert(G.nt() >= SolveOrder);
    assert(Sigma.nt() >= SolveOrder);
    if (size1 == 1) {
        dyson_start_ret<T, herm_matrix<T>, 1>(G, mu, H.ptr(0), Sigma, integration::I<T>(SolveOrder), h);
        dyson_start_tv<T, herm_matrix<T>, 1>(G, mu, H.ptr(0), Sigma, integration::I<T>(SolveOrder), beta, h);
        dyson_start_les<T, herm_matrix<T>, 1>(G, mu, H.ptr(0), Sigma, integration::I<T>(SolveOrder), beta, h);
    } else {
        dyson_start_ret<T, herm_matrix<T>, LARGESIZE>(G, mu, H.ptr(0), Sigma, integration::I<T>(SolveOrder), h);
        dyson_start_tv<T, herm_matrix<T>, LARGESIZE>(G, mu, H.ptr(0), Sigma, integration::I<T>(SolveOrder), beta, h);
        dyson_start_les<T, herm_matrix<T>, LARGESIZE>(G, mu, H.ptr(0), Sigma, integration::I<T>(SolveOrder), beta, h);
    }
}
/// @private
/** \brief <b> One step Dyson solver (integral-differential form) for a Green's function \f$G\f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a hermitian matrix \f$G(t, t^\prime)\f$ at a given timestep 'n',
* > i.e.,  G^ret(nh,t'<=nh), G^les(t<=nh,nh), G^tv(nt,tau=0..beta). Timestep must be >k,
* > where k is the Integration order 'I'.
* > The timesteps n=0..k must be computed seperately, using the routine
* > "_start", which assumes that the Matsubara component of \f$G\f$ and \f$\Sigma(t,t^\prime)\f$
* > for \f$t,t^\prime\f$<=k are given.
* > Here, are given: \f$\Sigma(t, t^\prime)\f$, \f$\mu\f$, and \f$H(t)\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] time step
* @param &G
* > [herm_matrix<T>] solution
* @param mu
* > [T] chemical potential
* @param &H
* > [function<T>] time-dependent function
* @param &Sigma
* > [herm_matrix<T>] self-energy
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
*/
template <typename T>
void dyson_timestep(int n, herm_matrix<T> &G, T mu, function<T> &H, herm_matrix<T> &Sigma,
                    integration::Integrator<T> &I, T beta, T h) {
    int size1 = G.size1(), k = I.k();
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert(G.nt() >= n);
    assert(Sigma.nt() >= n);
    assert(n > k);
    if (size1 == 1) {
        dyson_timestep_ret<T, herm_matrix<T>, 1>(n, G, mu, H.ptr(0), Sigma, I, h);
        dyson_timestep_tv<T, herm_matrix<T>, 1>(n, G, mu, H.ptr(n), Sigma, I, beta, h);
        dyson_timestep_les<T, herm_matrix<T>, 1>(n, G, mu, H.ptr(0), Sigma, I, beta, h);
    } else {
        dyson_timestep_ret<T, herm_matrix<T>, LARGESIZE>(n, G, mu, H.ptr(0), Sigma, I, h);
        dyson_timestep_tv<T, herm_matrix<T>, LARGESIZE>(n, G, mu, H.ptr(n), Sigma, I, beta,
                                                        h);
        dyson_timestep_les<T, herm_matrix<T>, LARGESIZE>(n, G, mu, H.ptr(0), Sigma, I, beta,
                                                         h);
    }
}


/** \brief <b> One step Dyson solver (integral-differential form) for a Green's function \f$G\f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a hermitian matrix \f$G(t, t^\prime)\f$ at a given timestep 'n',
* > i.e.,  G^ret(nh,t'<=nh), G^les(t<=nh,nh), G^tv(nt,tau=0..beta). Timestep must be >k,
* > where k is the Integration order 'I'.
* > The timesteps n=0..k must be computed seperately, using the routine
* > "_start", which assumes that the Matsubara component of \f$G\f$ and \f$\Sigma(t,t^\prime)\f$
* > for \f$t,t^\prime\f$<=k are given.
* > Here, are given: \f$\Sigma(t, t^\prime)\f$, \f$\mu\f$, and \f$H(t)\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] time step
* @param &G
* > [herm_matrix<T>] solution
* @param mu
* > [T] chemical potential
* @param &H
* > [function<T>] time-dependent function
* @param &Sigma
* > [herm_matrix<T>] self-energy
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T>
void dyson_timestep(int n, herm_matrix<T> &G, T mu, function<T> &H, herm_matrix<T> &Sigma,
                    T beta, T h, const int SolveOrder) {
    int size1 = G.size1();
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert(G.nt() >= n);
    assert(Sigma.nt() >= n);
    assert(n > SolveOrder);
    if (size1 == 1) {
        dyson_timestep_ret<T, herm_matrix<T>, 1>(n, G, mu, H.ptr(0), Sigma, integration::I<T>(SolveOrder), h);
        dyson_timestep_tv<T, herm_matrix<T>, 1>(n, G, mu, H.ptr(n), Sigma, integration::I<T>(SolveOrder), beta, h);
        dyson_timestep_les<T, herm_matrix<T>, 1>(n, G, mu, H.ptr(0), Sigma, integration::I<T>(SolveOrder), beta, h);
    } else {
        dyson_timestep_ret<T, herm_matrix<T>, LARGESIZE>(n, G, mu, H.ptr(0), Sigma, integration::I<T>(SolveOrder), h);
        dyson_timestep_tv<T, herm_matrix<T>, LARGESIZE>(n, G, mu, H.ptr(n), Sigma, integration::I<T>(SolveOrder), beta,
                                                        h);
        dyson_timestep_les<T, herm_matrix<T>, LARGESIZE>(n, G, mu, H.ptr(0), Sigma, integration::I<T>(SolveOrder), beta,
                                                         h);
    }
}
/// @private
/** \brief <b> Solver of the Dyson equation in the integral-differential form for a Green's function \f$G\f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a hermitian matrix \f$G(t, t^\prime)\f$
* > with given \f$\Sigma(t, t^\prime)\f$, \f$\mu\f$, and \f$H(t)\f$.
* > Here, one calls the routines 'dyson_mat()', 'dyson_start()', 'dyson_timestep'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix<T>] solution
* @param mu
* > [T] chemical potential
* @param &H
* > [function<T>] time-dependent function
* @param &Sigma
* > [herm_matrix<T>] self-energy
* @param I
* > [Integrator] integrator class
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
* @param matsubara_method
* > [int] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
* @param force_hermitian
* > [bool] force hermitian solution
*/
template <typename T>
void dyson(herm_matrix<T> &G, T mu, function<T> &H, herm_matrix<T> &Sigma,
           integration::Integrator<T> &I, T beta, T h, const int matsubara_method,
           const bool force_hermitian) {
    int n, k = I.k(), nt = G.nt();
    dyson_mat(G, Sigma, mu, H, I, beta,  matsubara_method, force_hermitian);
    if (nt >= 0)
        dyson_start(G, mu, H, Sigma, I, beta, h);
    for (n = k + 1; n <= nt; n++)
        dyson_timestep(n, G, mu, H, Sigma, I, beta, h);
}


/** \brief <b> Solver of the Dyson equation in the integral-differential form for a Green's function \f$G\f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > One solves the Dyson equation of the following form:
* > \f$ [ id/dt + \mu - H(t) ] G(t,t^\prime) - [\Sigma*G](t,t^\prime) = \delta(t,t^\prime)\f$
* > for a hermitian matrix \f$G(t, t^\prime)\f$
* > with given \f$\Sigma(t, t^\prime)\f$, \f$\mu\f$, and \f$H(t)\f$.
* > Here, one calls the routines 'dyson_mat()', 'dyson_start()', 'dyson_timestep'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param &G
* > [herm_matrix<T>] solution
* @param mu
* > [T] chemical potential
* @param &H
* > [function<T>] time-dependent function
* @param &Sigma
* > [herm_matrix<T>] self-energy
* @param beta
* > [double] inverse temperature
* @param h
* > [double] time interval
* @param SolveOrder
* > [int] integrator order
* @param matsubara_method
* > [int] Solution method on the Matsubara axis with 0: Fourier, 1: steep, 2: fixpoint
* @param force_hermitian
* > [bool] force hermitian solution
*/
template <typename T>
void dyson(herm_matrix<T> &G, T mu, function<T> &H, herm_matrix<T> &Sigma,
           T beta, T h, const int SolveOrder, const int matsubara_method,
           const bool force_hermitian) {
    int n, nt = G.nt();
    dyson_mat(G, mu, H, Sigma, beta, SolveOrder, matsubara_method, force_hermitian);
    if (nt >= 0)
        dyson_start(G, mu, H, Sigma, beta, h, SolveOrder);
    for (n = SolveOrder + 1; n <= nt; n++)
        dyson_timestep(n, G, mu, H, Sigma, beta, h, SolveOrder);
}

}
#endif  // CNTR_DYSON_IMPL_H