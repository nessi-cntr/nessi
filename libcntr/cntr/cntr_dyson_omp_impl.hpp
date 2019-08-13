#ifndef CNTR_DYSON_OMP_IMPL_H
#define CNTR_DYSON_OMP_IMPL_H

#include "cntr_dyson_omp_decl.hpp"
//#include "cntr_exception.hpp"
#include "cntr_elements.hpp"
#include "cntr_dyson_decl.hpp"
#include "cntr_dyson_impl.hpp"
#include "cntr_function_decl.hpp"
#include "cntr_herm_matrix_decl.hpp"
#include "cntr_herm_pseudo_decl.hpp"
#include "cntr_pseudo_convolution_decl.hpp"
#include "cntr_pseudodyson_decl.hpp"
#include "cntr_convolution_impl.hpp"

namespace cntr {

#if CNTR_USE_OMP == 1

/* #######################################################################################
#      DYSON EQUATION:
#
#          [ id/dt + mu - H(t) ] G(t,t') - [Sigma*G](t,t') = 1(t,t')
#          G(t,t')[-id/dt' + mu - H(t')] - [G*Sigma](t,t') = 1(t,t')
#
#      explanation see cntr_dyson
#      for omp paralellizaion I use conjugate equation, but this does not work for
#
#
###########################################################################################*/

/*###########################################################################################
#   RETARDED FUNCTION: (GG = herm_matrix or herm_pseudo)
###########################################################################################*/
template <typename T, class GG, int SIZE1>
void dyson_timestep_ret_omp(int omp_num_threads, int n, GG &G, T mu, std::complex<T> *H,
                            GG &Sigma, integration::Integrator<T> &I, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1;
    cplx cplx_i = cplx(0, 1);
    int size1 = G.size1();
    int sg = G.element_size();
    // check consistency:
    assert(k + 1<= n);
    assert(n<= Sigma.nt());
    assert(n<= G.nt());
    assert(G.sig()== Sigma.sig());

    ///////////////////////////////////////////////////////////////////////////////////////
    // SET ENTRIES IN TIMESTEP TO 0
    {
        cplx *gret = G.retptr(n, 0);
        int n1 = (n + 1) * sg, i;
        for (i = 0; i < n1; i++)
            gret[i] = 0;
    }
    ///////////////////////////////////////////////////////////////////////////////////////
    // INITIAL VALUE t' = n
    element_set<T, SIZE1>(size1, G.retptr(n, n), -cplx_i);
    ///////////////////////////////////////////////////////////////////////////////////////
    // START VALUES  t' = n-j, j = 1...k: solve a kxk problem
    // this is the same as before, i.e., use conjugate equation ii*d/dt1 G(t,t1) = ...
    // using the other equation ii*d/dt G(t,t1) = ... seems to be unstable
    // (may be interesting to investigate this instability in general)
    {
        int i, j, p, l, q;
        cplx w0, cweight;
        T weight;
        cplx *gtemp = new cplx[k * sg];
        cplx *diffw = new cplx[k1 + 1];
        cplx *qq = new cplx[(n + 1) * sg];
        cplx *one = new cplx[sg];
        cplx *mm = new cplx[k * k * sg];
        cplx *hj = new cplx[sg];
        cplx *stemp = new cplx[sg]; // sic
        element_set<T, SIZE1>(size1, one, 1);
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
                    element_incr<T, SIZE1>(size1, qq + p * sg, weight, G.retptr(n, n),
                                           stemp);
                } else {
                    q = l - 1;
                    element_incr<T, SIZE1>(size1, mm + sg * (p * k + q), -weight, stemp);
                }
            }
        }
        element_linsolve_left<T, SIZE1>(size1, k, gtemp, mm, qq); // gtemp * mm = qq
        for (j = 1; j <= k; j++)
            element_set<T, SIZE1>(size1, G.retptr(n, n - j), gtemp + (j - 1) * sg);
        delete[] gtemp;
        delete[] diffw;
        delete[] qq;
        delete[] one;
        delete[] mm;
        delete[] hj;
        delete[] stemp;
    }
///////////////////////////////////////////////////////////////////////////////////////
// now use equation ii*d/dt G(t,t1) = ... to compute G(n*h,j*h),j=0 ... n-k-1
// OMP parallelization over j

#pragma omp parallel num_threads(omp_num_threads)
    {
        // convolution Sigma*G ->> written to G, on
        int j, p, i;
        int nomp = omp_get_num_threads();
        int tid = omp_get_thread_num();
        std::vector<bool> mask_ret(n + 1, false);
        cplx w0 = h * I.gregory_omega(0);
        cplx *diffw = new cplx[k1 + 1];
        cplx *qq = new cplx[sg];
        cplx *mm = new cplx[sg];
        for (i = 0; i < n - k; i++)
            if (i % nomp == tid)
                mask_ret[i] = true;
        incr_convolution_ret<T, GG, SIZE1>(n, mask_ret, cplx(1.0, 0.0), G, Sigma, Sigma,
                                           NULL, G, G, I, h);
        for (p = 0; p <= k1; p++)
            diffw[p] = I.bd_weights(p) * cplx_i / h; // use BD(k+1!!)
        for (j = 0; j < n - k; j++) {
            if (mask_ret[j]) {
                element_set<T, SIZE1>(size1, qq, G.retptr(n, j)); // << Sigma*G(n,j)
                // set up mm and qqj for 1x1 problem:
                for (p = 1; p <= k1; p++)
                    element_incr<T, SIZE1>(size1, qq, -diffw[p], G.retptr(n - p, j));
                element_set<T, SIZE1>(size1, mm, diffw[0] + mu);
                element_incr<T, SIZE1>(size1, mm, -w0, Sigma.retptr(j, j));
                element_incr<T, SIZE1>(size1, mm, cplx(-1.0, 0.0), H + n * sg);
                element_linsolve_right<T, SIZE1>(size1, G.retptr(n, j), mm, qq); // mm*G=qq
            }
        }
        delete[] diffw;
        delete[] qq;
        delete[] mm;
    }
    return;
}
// retarded start is not omp-paralellized
/*###########################################################################################
#   TV FUNCTION:
###########################################################################################*/
// GG = herm_matrix
template <typename T, class GG, int SIZE1>
void dyson_timestep_tv_omp(int omp_num_threads, int n, GG &G, T mu, std::complex<T> *Hn,
                           GG &Sigma, integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int size1 = G.size1();
    int k = I.get_k(), k1 = k + 1;
    int sg = G.element_size();
    int ntau = G.ntau();
    cplx ih = cplx(0, 1.0 / h);
    assert(k + 1<= n);
    assert(n<= Sigma.nt());
    assert(n<= G.nt());
    assert(G.sig()== Sigma.sig());

    for (int j = 0; j <= ntau; j++)
        element_set_zero<T, SIZE1>(size1, G.tvptr(n, j));
#pragma omp parallel num_threads(omp_num_threads)
    {
        // convolution Sigma*G ->> written to G.tv
        int j, p, i;
        int nomp = omp_get_num_threads();
        int tid = omp_get_thread_num();
        std::vector<bool> mask(ntau + 1, false);
        cplx cweight;
        cplx *diffw = new cplx[k1 + 1];
        cplx *qq = new cplx[sg];
        cplx *mm = new cplx[sg];
        for (i = 0; i <= ntau; i++)
            if (i % nomp == tid)
                mask[i] = true;
        incr_convolution_tv<T, GG, SIZE1>(n, mask, cplx(1.0, 0.0), G, Sigma, Sigma, NULL,
                                          NULL, G, G, I, beta, h);
        // Now solve
        // [ i/h bd(0) - H - h w(n,0) Sigma(n,n) ] G(n,m)  = Q(m),
        // where Q is initially stored in G(n,m)
        element_set<T, SIZE1>(size1, mm, ih * I.bd_weights(0) + mu);
        element_incr<T, SIZE1>(size1, mm, cplx(-1.0, 0.0), Hn);
        cweight = -h * I.gregory_weights(n, 0);
        element_incr<T, SIZE1>(size1, mm, cweight, Sigma.retptr(n, n));
        // ACCUMULATE CONTRIBUTION TO id/dt G(t,t') FROM t=mh, m=n-k..n-1
        for (p = 0; p <= k1; p++)
            diffw[p] = ih * I.bd_weights(p); // use BD(k+1!!)
        for (j = 0; j <= ntau; j++) {
            if (mask[j]) {
                element_set<T, SIZE1>(size1, qq, G.tvptr(n, j));
                for (p = 1; p <= k1; p++)
                    element_incr<T, SIZE1>(size1, qq, -diffw[p], G.tvptr(n - p, j));
                element_linsolve_right<T, SIZE1>(size1, G.tvptr(n, j), mm, qq);
            }
        }
        delete[] qq;
        delete[] mm;
        delete[] diffw;
    }
    return;
}
// GG = pseudo_matrix the only differebnce is the convolution!
/// @private
template <typename T, class GG, int SIZE1>
void pseudodyson_timestep_tv_omp(int omp_num_threads, int n, GG &G, T mu,
                                 std::complex<T> *Hn, GG &Sigma,
                                 integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int size1 = G.size1();
    int k = I.get_k(), k1 = k + 1;
    int sg = G.element_size();
    int ntau = G.ntau();
    cplx ih = cplx(0, 1.0 / h);
    assert(k + 1<= n);
    assert(n<= Sigma.nt());
    assert(n<= G.nt());
    assert(G.sig()== Sigma.sig());

    for (int j = 0; j <= ntau; j++)
        element_set_zero<T, SIZE1>(size1, G.tvptr(n, j));
#pragma omp parallel num_threads(omp_num_threads)
    {
        // convolution Sigma*G ->> written to G.tv
        int j, p, i;
        int nomp = omp_get_num_threads();
        int tid = omp_get_thread_num();
        std::vector<bool> mask(ntau + 1, false);
        cplx cweight;
        cplx *diffw = new cplx[k1 + 1];
        cplx *qq = new cplx[sg];
        cplx *mm = new cplx[sg];
        for (i = 0; i <= ntau; i++)
            if (i % nomp == tid)
                mask[i] = true;
        incr_pseudo_convolution_tv<T, GG, SIZE1>(n, mask, cplx(1.0, 0.0), G, Sigma, Sigma,
                                                 NULL, NULL, G, G, I, beta, h);
        // Now solve
        // [ i/h bd(0) - H - h w(n,0) Sigma(n,n) ] G(n,m)  = Q(m),
        // where Q is initially stored in G(n,m)
        element_set<T, SIZE1>(size1, mm, ih * I.bd_weights(0) + mu);
        element_incr<T, SIZE1>(size1, mm, cplx(-1.0, 0.0), Hn);
        cweight = -h * I.gregory_weights(n, 0);
        element_incr<T, SIZE1>(size1, mm, cweight, Sigma.retptr(n, n));
        // ACCUMULATE CONTRIBUTION TO id/dt G(t,t') FROM t=mh, m=n-k..n-1
        for (p = 0; p <= k1; p++)
            diffw[p] = ih * I.bd_weights(p); // use BD(k+1!!)
        for (j = 0; j <= ntau; j++) {
            if (mask[j]) {
                element_set<T, SIZE1>(size1, qq, G.tvptr(n, j));
                for (p = 1; p <= k1; p++)
                    element_incr<T, SIZE1>(size1, qq, -diffw[p], G.tvptr(n - p, j));
                element_linsolve_right<T, SIZE1>(size1, G.tvptr(n, j), mm, qq);
            }
        }
        delete[] qq;
        delete[] mm;
        delete[] diffw;
    }
    return;
}
/*###########################################################################################
#   LESSER FUNCTION: (GG = herm_matrix or herm_pseudo)
###########################################################################################*/
template <typename T, class GG, int SIZE1>
void dyson_timestep_les_omp(int omp_num_threads, int n, GG &G, T mu, std::complex<T> *H,
                            GG &Sigma, integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    cplx cplx_i = cplx(0, 1);
    int k = I.get_k(), k1 = k + 1;
    int size1 = G.size1();
    int sg = G.element_size();
    int n1 = (n > k ? n : k);
    //////////////////////////////////////////////////////////////////////////////////////////
    // check consistency:  (more assertations follow in convolution)
    assert(k + 1<= n);
    assert(G.ntau()<= Sigma.ntau());
    assert(n1<= Sigma.nt());
    assert(n1<= G.nt());
    assert(G.sig()== Sigma.sig());
    // OMP PARALELLIZARION STARTS ONLY FOR n>=2*k+1
    if (n < 2 * k + 1) {
        return dyson_timestep_les<T, GG, SIZE1>(n, G, mu, H, Sigma, I, beta, h);
    }
    for (int j = 0; j <= n; j++)
        element_set_zero<T, SIZE1>(size1, G.lesptr(j, n));
///////////////////////////////////////////////////////////////////////////////////////////
// get G(j,n), j=0...n-k-1 from d/dt' G(t,t') equation
#pragma omp parallel num_threads(omp_num_threads)
    {
        int j, p, i;
        int nomp = omp_get_num_threads();
        int tid = omp_get_thread_num();
        std::vector<bool> mask_les(n + 1, false);
        cplx w0 = h * I.gregory_omega(0);
        cplx *diffw = new cplx[k1 + 1];
        cplx *qq = new cplx[sg];
        cplx *mm = new cplx[sg];
        cplx *stemp = new cplx[sg];
        // convolution Sigma*G ->> written to G
        for (i = 0; i < n - k; i++)
            if (i % nomp == tid)
                mask_les[i] = true;
        incr_convolution_les<T, GG, SIZE1>(n, mask_les, cplx(1.0, 0.0), G, G, G, NULL, NULL,
                                           Sigma, Sigma, I, beta, h);
        for (p = 0; p <= k1; p++)
            diffw[p] = I.bd_weights(p) * cplx_i / h; // use BD(k+1!!)
        for (j = 0; j < n - k; j++) {
            if (mask_les[j]) {
                element_set<T, SIZE1>(size1, qq, G.lesptr(j, n)); // << G*Sigma(j,n)
                // set up mm and qqj for 1x1 problem:
                for (p = 1; p <= k1; p++)
                    element_incr<T, SIZE1>(size1, qq, diffw[p], G.lesptr(j, n - p));
                element_set<T, SIZE1>(size1, mm, -diffw[0] + mu);
                element_conj<T, SIZE1>(size1, stemp, Sigma.retptr(j, j));
                element_incr<T, SIZE1>(size1, mm, -w0, stemp);
                element_incr<T, SIZE1>(size1, mm, cplx(-1.0, 0.0), H + n * sg);
                element_linsolve_left<T, SIZE1>(size1, G.lesptr(j, n), mm, qq);
            }
        }
        delete[] diffw;
        delete[] qq;
        delete[] mm;
        delete[] stemp;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    // get G(j,n), j=n-k...n from d/dt G(t,t') equation (old implementation)
    // currently not paralellized
    {
        int j, p, m;
        cplx *gles = new cplx[(n + 1) * sg];
        cplx *qq = new cplx[k * sg];
        cplx *mm = new cplx[k * k * sg];
        cplx cweight;
// CONVOLUTION SIGMA*G:  --->  G^les(j,n) j=n-k...n
// Note: this is only the tv*vt + les*adv part, Gles is not adressed
// coupld parallelize this:
#pragma omp parallel num_threads(omp_num_threads)
        {
            int j1;
            int nomp = omp_get_num_threads();
            int tid = omp_get_thread_num();
            for (j1 = n - k; j1 <= n; j1++) {
                if ((n - j1) % nomp == tid) {
                    element_set_zero<T, SIZE1>(size1, gles + j1 * sg);
                    convolution_timestep_les_tvvt<T, GG, SIZE1>(n, j1, j1, gles, G, Sigma,
                                                                Sigma, G, G, I, beta, h);
                    convolution_timestep_les_lesadv<T, GG, SIZE1>(n, j1, j1, gles, G, Sigma,
                                                                  Sigma, G, G, I, beta, h);
                }
            }
        }
        for (j = n - k; j <= n; j++) {
            // CONTRIBUTION FROM INTEGRAL tv*vt+les*adv
            element_set<T, SIZE1>(size1, qq, gles + j * sg);
            // ACCUMULATE CONTRIBUTION TO id/dt G(j-p,n) p=1...k1 into qq
            for (p = 1; p <= k1; p++) { // use BD(k+1) !!!
                cweight = -cplx_i / h * I.bd_weights(p);
                element_incr<T, SIZE1>(size1, qq, cweight, G.lesptr(j - p, n));
            }
            element_set<T, SIZE1>(size1, mm, cplx_i / h * I.bd_weights(0) + mu);
            element_incr<T, SIZE1>(size1, mm, cplx(-1.0, 0.0), H + sg * j);
            cweight = -h * I.gregory_weights(j, j);
            element_incr<T, SIZE1>(size1, mm, cweight, Sigma.retptr(j, j));
            for (m = 0; m < j; m++) {
                cweight = h * I.gregory_weights(j, m);
                element_incr<T, SIZE1>(size1, qq, cweight, Sigma.retptr(j, m),
                                       G.lesptr(m, n));
            }
            element_linsolve_right<T, SIZE1>(size1, G.lesptr(j, n), mm, qq);
        }
        delete[] gles;
        delete[] qq;
        delete[] mm;
    }
    return;
}
//// start les uses d/dt Gles(t,n*h)=... equation ... not OMP paralellized

////////////////////////////////////////////////////////////////////////////////////////////////
// main implementation:
// with function object:
/// @private
template <typename T>
void pseudodyson_timestep_omp(int omp_num_threads, int n, herm_pseudo<T> &G, T lam0,
                              function<T> &H, herm_pseudo<T> &Sigma,
                              integration::Integrator<T> &I, T beta, T h) {
    int size1 = G.size1(), k = I.k();
    int omp_num_threads1 = (omp_num_threads == -1 ? omp_get_max_threads() : omp_num_threads);
    assert(k + 1<= n);
    assert(n<= Sigma.nt());
    assert(n<= G.nt());
    assert(n<= H.nt());
    assert(G.sig()== Sigma.sig());
    assert(G.size1()== Sigma.size1());
    assert(G.size1()== H.size1());
    assert(G.ntau()== Sigma.ntau());
    if (size1 == 1) {
        dyson_timestep_ret_omp<T, herm_pseudo<T>, 1>(omp_num_threads1, n, G, lam0, H.ptr(0),
                                                     Sigma, I, h);
        pseudodyson_timestep_tv_omp<T, herm_pseudo<T>, 1>(omp_num_threads1, n, G, lam0,
                                                          H.ptr(n), Sigma, I, beta, h);
        dyson_timestep_les_omp<T, herm_pseudo<T>, 1>(omp_num_threads1, n, G, lam0, H.ptr(0),
                                                     Sigma, I, beta, h);
    } else {
        dyson_timestep_ret_omp<T, herm_pseudo<T>, LARGESIZE>(omp_num_threads1, n, G, lam0,
                                                             H.ptr(0), Sigma, I, h);
        pseudodyson_timestep_tv_omp<T, herm_pseudo<T>, LARGESIZE>(
            omp_num_threads1, n, G, lam0, H.ptr(n), Sigma, I, beta, h);
        dyson_timestep_les_omp<T, herm_pseudo<T>, LARGESIZE>(omp_num_threads1, n, G, lam0,
                                                             H.ptr(0), Sigma, I, beta, h);
    }
}
// with raw pointer:
/// @private
template <typename T>
void pseudodyson_timestep_omp(int omp_num_threads, int n, herm_pseudo<T> &G, T lam0,
                              std::complex<T> *Ht, herm_pseudo<T> &Sigma,
                              integration::Integrator<T> &I, T beta, T h) {
    int size1 = G.size1(), k = I.k();
    int omp_num_threads1 = (omp_num_threads == -1 ? omp_get_max_threads() : omp_num_threads);
    assert(k + 1<= n);
    assert(n<= Sigma.nt());
    assert(n<= G.nt());
    assert(G.sig()== Sigma.sig());
    assert(G.size1()== Sigma.size1());
    assert(G.ntau()== Sigma.ntau());
    if (size1 == 1) {
        dyson_timestep_ret_omp<T, herm_pseudo<T>, 1>(omp_num_threads1, n, G, lam0, Ht, Sigma,
                                                     I, h);
        pseudodyson_timestep_tv_omp<T, herm_pseudo<T>, 1>(
            omp_num_threads1, n, G, lam0, Ht + n * size1 * size1, Sigma, I, beta, h);
        dyson_timestep_les_omp<T, herm_pseudo<T>, 1>(omp_num_threads1, n, G, lam0, Ht, Sigma,
                                                     I, beta, h);
    } else {
        dyson_timestep_ret_omp<T, herm_pseudo<T>, LARGESIZE>(omp_num_threads1, n, G, lam0,
                                                             Ht, Sigma, I, h);
        pseudodyson_timestep_tv_omp<T, herm_pseudo<T>, LARGESIZE>(
            omp_num_threads1, n, G, lam0, Ht + n * size1 * size1, Sigma, I, beta, h);
        dyson_timestep_les_omp<T, herm_pseudo<T>, LARGESIZE>(omp_num_threads1, n, G, lam0,
                                                             Ht, Sigma, I, beta, h);
    }
}
// only for compatibility, works for size1 only
/// @private
template <typename T>
void pseudodyson_timestep_omp(int omp_num_threads, int n, herm_pseudo<T> &G, T lam0,
                              std::vector<std::complex<T>> &Ht, herm_pseudo<T> &Sigma,
                              integration::Integrator<T> &I, T beta, T h) {
    int size1 = G.size1(), k = I.k();
    std::complex<T> *hh = new std::complex<T>[n + 1];
    int omp_num_threads1 = (omp_num_threads == -1 ? omp_get_max_threads() : omp_num_threads);
    assert(k + 1<= n);
    assert(size1== 1);
    assert(n<= Sigma.nt());
    assert(n<= G.nt());
    assert(G.sig()== Sigma.sig());
    assert(G.size1()== Sigma.size1());
    assert(G.ntau()== Sigma.ntau());
    assert(n + 1<= (int)Ht.size());
    for (int j = 0; j <= n; j++)
        hh[j] = Ht[j];
    dyson_timestep_ret_omp<T, herm_pseudo<T>, 1>(omp_num_threads1, n, G, lam0, hh, Sigma, I,
                                                 h);
    pseudodyson_timestep_tv_omp<T, herm_pseudo<T>, 1>(
        omp_num_threads1, n, G, lam0, hh + n * size1 * size1, Sigma, I, beta, h);
    dyson_timestep_les_omp<T, herm_pseudo<T>, 1>(omp_num_threads1, n, G, lam0, hh, Sigma, I,
                                                 beta, h);
    delete[] hh;
}

template <typename T>
void dyson_timestep_omp(int omp_num_threads, int n, herm_matrix<T> &G, T lam0,
                        function<T> &H, herm_matrix<T> &Sigma, integration::Integrator<T> &I,
                        T beta, T h) {
    int size1 = G.size1(), k = I.k();
    int omp_num_threads1 = (omp_num_threads == -1 ? omp_get_max_threads() : omp_num_threads);
    assert(k + 1<= n);
    assert(n<= Sigma.nt());
    assert(n<= G.nt());
    assert(n<= H.nt());
    assert(G.sig()== Sigma.sig());
    assert(G.size1()== Sigma.size1());
    assert(G.size1()== H.size1());
    assert(G.ntau()== Sigma.ntau());
    if (size1 == 1) {
        dyson_timestep_ret_omp<T, herm_matrix<T>, 1>(omp_num_threads1, n, G, lam0, H.ptr(0),
                                                     Sigma, I, h);
        dyson_timestep_tv_omp<T, herm_matrix<T>, 1>(omp_num_threads1, n, G, lam0, H.ptr(n),
                                                    Sigma, I, beta, h);
        dyson_timestep_les_omp<T, herm_matrix<T>, 1>(omp_num_threads1, n, G, lam0, H.ptr(0),
                                                     Sigma, I, beta, h);
    } else {
        dyson_timestep_ret_omp<T, herm_matrix<T>, LARGESIZE>(omp_num_threads1, n, G, lam0,
                                                             H.ptr(0), Sigma, I, h);
        dyson_timestep_tv_omp<T, herm_matrix<T>, LARGESIZE>(omp_num_threads1, n, G, lam0,
                                                            H.ptr(n), Sigma, I, beta, h);
        dyson_timestep_les_omp<T, herm_matrix<T>, LARGESIZE>(omp_num_threads1, n, G, lam0,
                                                             H.ptr(0), Sigma, I, beta, h);
    }
}

#endif // CNTR_USE_OMP

}  // namespace cntr

#endif  // CNTR_DYSON_OMP_IMPL_H
