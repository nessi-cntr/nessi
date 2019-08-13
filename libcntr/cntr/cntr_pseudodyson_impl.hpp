#ifndef CNTR_PSEUDODYSON_IMPL_H
#define CNTR_PSEUDODYSON_IMPL_H

#include "cntr_pseudodyson_decl.hpp"

namespace cntr {

/*#########################################################################################
#
#  PSEUDOPARTICLE DYSON EQUATION: [ id/dt + lam0 - H -Sigma ] G = 1
#
#  lam0 must be included in H
#
#  ANALOGOUS TO USUAL DYSON, BUT WITH DIFFERENT CONVOLUTION
#
##########################################################################################*/
/// @private
template <typename T, int SIZE1>
void pseudodyson_mat_dispatch(herm_pseudo<T> &G, T lam0, std::complex<T> *H0,
                              herm_pseudo<T> &Sigma, integration::Integrator<T> &I, T beta) {
    typedef std::complex<T> cplx;
    int ntau, k = I.get_k(), k1 = k + 1, sg, l, m, p, q, n, j, size1 = G.size1();
    T dtau;
    cplx *mm, *qq, *gtemp, *one, *stemp, cweight;

    ntau = G.ntau();
    assert(ntau > k);
    assert(Sigma.ntau() == ntau);
    sg = G.element_size();
    dtau = beta / ntau;

    mm = new cplx[k * k * sg];
    qq = new cplx[k * sg];
    gtemp = new cplx[k * sg];
    stemp = new cplx[k * sg];
    one = new cplx[sg];
    element_set<T, SIZE1>(size1, one, 1.0);
    // INITIAL VALUE : GMAT=-1
    for (l = 0; l < sg; l++)
        G.matptr(0)[l] = -one[l];
    // THE FIRST K STEPS Gmat(n) for n=1...k
    for (l = 0; l < k * k * sg; l++)
        mm[l] = 0;
    for (l = 0; l < k * sg; l++)
        qq[l] = 0;
    // derive linear equations for gtemp(p)=Gmat(p+1)
    // mm(p,q)*G(q)=Q(p) for p,q=0...k-1
    for (n = 1; n <= k; n++) {
        p = n - 1;
        // derivative -d/dtau Gmat(n)
        for (m = 0; m <= k; m++) {
            cweight = -1 / dtau * I.poly_differentiation(n, m);
            if (m == 0) {
                for (l = 0; l < sg; l++)
                    qq[p * sg + l] -= cweight * G.matptr(0)[l];
            } else { // goes into mm(p,q)
                q = m - 1;
                for (l = 0; l < sg; l++)
                    mm[sg * (p * k + q) + l] += cweight * one[l];
            }
        }
        // H -- goes into m(p,p)
        element_set<T, SIZE1>(size1, gtemp, H0);
        for (l = 0; l < sg; l++)
            mm[sg * (p * k + p) + l] += lam0 * one[l] - gtemp[l];
        // integral int_0^n dx Sigma(n*dtau-x)G(x) :
        if (n < k) {
            // need the rcorr for n<k:
            // Zitat: R=int_0^n dx a(n-x) b(x) = sum_{m,j=0}^k I.rcorr(n,m,j) a(m)*b(j)
            for (m = 0; m <= k; m++) {
                element_set<T, SIZE1>(size1, stemp, Sigma.matptr(m));
                for (j = 0; j <= k; j++) {
                    cweight = dtau * I.rcorr(n, m, j);
                    if (j == 0) { // goes into qq(p)
                        element_incr<T, SIZE1>(size1, qq + p * sg, cweight, stemp,
                                               G.matptr(0));
                    } else { // goes into mm(p,q)
                        q = j - 1;
                        for (l = 0; l < sg; l++)
                            mm[sg * (p * k + q) + l] -= stemp[l] * cweight;
                    }
                }
            }
        } else {
            // usual integral for n=k
            // R=int_0^k dx a(n-x) b(x) = sum_{m=0}^k I.gregory_weights(k,m) a(k-m)*b(m)
            for (m = 0; m <= k; m++) {
                cweight = dtau * I.gregory_weights(k, m);
                element_set<T, SIZE1>(size1, stemp, Sigma.matptr(k - m));
                if (m == 0) { // goes into qq(p)
                    element_incr<T, SIZE1>(size1, qq + p * sg, cweight, stemp, G.matptr(0));
                } else { // goes into mm(p,q)
                    q = m - 1;
                    for (l = 0; l < sg; l++)
                        mm[sg * (p * k + q) + l] -= stemp[l] * cweight;
                }
            }
        }
    }

    /*   std::cout.precision(15);
        std::cout << "mat: " ;
        for(p=0;p<k;p++){ std::cout << "\t" << qq[p];}
        std::cout << std::endl;
        */

    element_linsolve_right<T, SIZE1>(size1, k, gtemp, mm,
                                     qq); // solve kXk problem mm*gtemp=qq
    for (p = 0; p < k; p++)
        element_set<T, SIZE1>(size1, G.matptr(p + 1), gtemp + sg * p);
    // the remaining steps
    for (n = k1; n <= ntau; n++) {
        for (l = 0; l < sg; l++)
            mm[l] = 0;
        for (l = 0; l < sg; l++)
            qq[l] = 0;
        // integral
        if (n < 2 * k1) {
            for (m = 0; m < n; m++) {
                cweight = I.gregory_weights(n, m);
                element_incr<T, SIZE1>(size1, qq, cweight, Sigma.matptr(n - m), G.matptr(m));
            }
        } else { // some weights are 1
            for (m = 0; m <= k; m++) {
                cweight = I.gregory_weights(n, m);
                element_incr<T, SIZE1>(size1, qq, cweight, Sigma.matptr(n - m), G.matptr(m));
            }
            for (m = k1; m < n - k; m++) {
                element_incr<T, SIZE1>(size1, qq, Sigma.matptr(n - m), G.matptr(m));
            }
            for (m = n - k; m < n; m++) {
                cweight = I.gregory_omega(n - m);
                element_incr<T, SIZE1>(size1, qq, cweight, Sigma.matptr(n - m), G.matptr(m));
            }
        }
        element_smul<T, SIZE1>(size1, qq, dtau);
        // derivative -d/dtau Gmat(n) :use BD(k+1)
        for (m = 1; m <= k1; m++) {
            cweight = 1 / dtau * I.bd_weights(m);
            for (l = 0; l < sg; l++)
                qq[l] += G.matptr(n - m)[l] * cweight;
        }
        // determine mm
        element_set<T, SIZE1>(size1, gtemp, H0);
        element_set<T, SIZE1>(size1, stemp, Sigma.matptr(0));
        for (l = 0; l < sg; l++)
            mm[l] = -I.bd_weights(0) / dtau * one[l] + lam0 * one[l] - gtemp[l] -
                    dtau * I.gregory_weights(n, n) * stemp[l];
        element_linsolve_right<T, SIZE1>(size1, G.matptr(n), mm,
                                         qq); // solve 1X1 problem mm*gtemp=qq
    }
    delete[] stemp;
    delete[] gtemp;
    delete[] mm;
    delete[] qq;
    delete[] one;
    return;
}
/// @private
template <typename T, int SIZE1>
void pseudodyson_timestep_tv(int n, herm_pseudo<T> &G, T lam0, std::complex<T> *H0,
                             herm_pseudo<T> &Sigma, integration::Integrator<T> &I, T beta,
                             T h) {
    // tv: solve conjugate equation: (quite analogous to matsubara)
    // d/dm Gtv(n,m) - Gtv(n,m) H0* = int_{ntau-m}^ntau  dx Gtv(n,x) Smat(x-m)  + int_0^n dj
    // Gret(n,j) Stv(j,m)
    // start this equation from beta and integrate to 0: m -> l=beta-m, X(l)=Gtv(n.beta-l)
    // -d/dl X(l) - X(l) H0* = int_0^l dj X(j) Smat(l-j)  + int_0^n dj Gret(n,j)
    // Stv(j,beta-l)
    // NOTE: restriction of dtau integral to x between m and n on CC
    typedef std::complex<T> cplx;
    int ntau, k = I.get_k(), k1 = k + 1, sg, l, m, p, q, i, j, n1, sig, size1 = G.size1();
    T dtau;
    cplx *mm, *qq, *gtemp, *htemp, *one, *stemp, *xx, cweight;

    n1 = (n > k ? n : k);
    ntau = G.ntau();
    assert(ntau > k);
    assert(Sigma.ntau() == ntau);
    assert(G.nt() >= n1);
    assert(Sigma.nt() >= n1);
    sig = G.sig();
    assert(sig == Sigma.sig());

    sg = G.element_size();
    dtau = beta / ntau;

    mm = new cplx[k * k * sg];
    qq = new cplx[k * sg];
    xx = new cplx[(ntau + 1) * sg];
    stemp = new cplx[k * sg];
    htemp = new cplx[sg];
    gtemp = new cplx[sg];
    one = new cplx[sg];
    element_set<T, SIZE1>(size1, one, 1.0);
    element_set<T, SIZE1>(size1, htemp, H0);
    element_conj<T, SIZE1>(size1, htemp);
    for (i = 0; i < sg * (ntau + 1); i++)
        xx[i] = 0;
    // int_0^n dj Gret(n,j) Stv(j,beta-i)  --->  X(i) for i=1..ntau
    for (j = 0; j <= n1; j++) {
        if (j <= n) {
            element_set<T, SIZE1>(size1, gtemp, G.retptr(n, j));
        } else {
            element_conj<T, SIZE1>(size1, gtemp, G.retptr(j, n));
            element_smul<T, SIZE1>(size1, gtemp, -1);
        }
        element_smul<T, SIZE1>(size1, gtemp, I.gregory_weights(n, j) * h);
        for (i = 1; i <= ntau; i++) {
            element_incr<T, SIZE1>(size1, xx + i * sg, gtemp, Sigma.tvptr(j, ntau - i));
        }
    }
    // INITIAL VALUE : X(0) = Gtv(t,beta) = -sig Gvt(beta-beta,t)* =  sig Gret(t,0)
    element_set<T, SIZE1>(size1, xx, G.retptr(n, 0));
    element_smul<T, SIZE1>(size1, xx, sig);
    // THE FIRST K STEPS X(i) for i=1...k
    for (l = 0; l < k * k * sg; l++)
        mm[l] = 0;
    for (l = 0; l < k * sg; l++)
        qq[l] = 0;
    // derive linear equations for
    // X(q+1) mm(p,q)=Q(p) for p,q=0...k-1
    for (i = 1; i <= k; i++) {
        p = i - 1;
        element_set<T, SIZE1>(size1, qq + p * sg, xx + i * sg); /// CHECK: xx=0 for n=0???
        // derivative -d/di X(i)
        for (m = 0; m <= k; m++) {
            cweight = -1 / dtau * I.poly_differentiation(i, m);
            if (m == 0) {
                for (l = 0; l < sg; l++)
                    qq[p * sg + l] -= cweight * xx[sg * 0 + l];
            } else { // goes into mm(p,q)
                q = m - 1;
                for (l = 0; l < sg; l++)
                    mm[sg * (p * k + q) + l] += cweight * one[l];
            }
        }
        // H -- goes into m(p,p)
        for (l = 0; l < sg; l++)
            mm[sg * (p * k + p) + l] += lam0 * one[l] - htemp[l];
        // integral int_0^i dx X(x) Smat(i-x):
        if (i < k) {
            // need the rcorr:
            // Zitat: R=int_0^i dx a(i-x) b(x) = sum_{m,j=0}^k I.rcorr(i,m,j) a(m)*b(j)
            for (m = 0; m <= k; m++) {
                element_set<T, SIZE1>(size1, stemp, Sigma.matptr(m));
                for (j = 0; j <= k; j++) {
                    cweight = dtau * I.rcorr(i, m, j);
                    if (j == 0) { // goes into qq(p)
                        element_incr<T, SIZE1>(size1, qq + p * sg, cweight, xx + 0 * sg,
                                               stemp);
                    } else { // goes into mm(p,q)
                        q = j - 1;
                        for (l = 0; l < sg; l++)
                            mm[sg * (p * k + q) + l] -= stemp[l] * cweight;
                    }
                }
            }
        } else { // i=k
                 // usual integral for n=k
            // R=int_0^k dx a(n-x) b(x) = sum_{m=0}^k I.gregory_weights(k,m) a(k-m)*b(m)
            for (m = 0; m <= k; m++) {
                element_set<T, SIZE1>(size1, stemp, Sigma.matptr(k - m));
                cweight = dtau * I.gregory_weights(i, m);
                if (m == 0) { // goes into qq(p)
                    element_incr<T, SIZE1>(size1, qq + p * sg, cweight, xx + 0 * sg, stemp);
                } else { // goes into mm(p,q)
                    q = m - 1;
                    for (l = 0; l < sg; l++)
                        mm[sg * (p * k + q) + l] -= stemp[l] * cweight;
                }
            }
        }
    }
    /*   std::cout.precision(15);
     std::cout << "tv: " ;
        for(p=0;p<k;p++){ std::cout << "\t" << qq[p];}
        std::cout << std::endl;
        */

    element_linsolve_left<T, SIZE1>(size1, k, xx + 1 * sg, mm,
                                    qq); // solve kXk problem xx*mm=qq
    // the remaining steps
    for (i = k1; i <= ntau; i++) {
        for (l = 0; l < sg; l++)
            mm[l] = 0;
        for (l = 0; l < sg; l++)
            qq[l] = 0;
        // integral
        if (i < 2 * k1) {
            for (m = 0; m < i; m++) {
                cweight = I.gregory_weights(i, m);
                element_incr<T, SIZE1>(size1, qq, cweight, xx + m * sg, Sigma.matptr(i - m));
            }
        } else { // some weights are 1
            for (m = 0; m <= k; m++) {
                cweight = I.gregory_weights(i, m);
                element_incr<T, SIZE1>(size1, qq, cweight, xx + m * sg, Sigma.matptr(i - m));
            }
            for (m = k1; m < i - k; m++) {
                element_incr<T, SIZE1>(size1, qq, xx + m * sg, Sigma.matptr(i - m));
            }
            for (m = i - k; m < i; m++) {
                cweight = I.gregory_omega(i - m);
                element_incr<T, SIZE1>(size1, qq, cweight, xx + m * sg, Sigma.matptr(i - m));
            }
        }
        element_smul<T, SIZE1>(size1, qq, dtau);
        element_incr<T, SIZE1>(size1, qq, xx + i * sg);
        // derivative -d/dtau Gmat(i) :use BD(k+1)
        for (m = 1; m <= k1; m++) {
            cweight = 1 / dtau * I.bd_weights(m);
            for (l = 0; l < sg; l++)
                qq[l] += xx[sg * (i - m) + l] * cweight;
        }
        // determine mm
        element_set<T, SIZE1>(size1, stemp, Sigma.matptr(0));
        for (l = 0; l < sg; l++)
            mm[l] = -I.bd_weights(0) / dtau * one[l] + (lam0 * one[l] - htemp[l]) -
                    dtau * I.gregory_weights(i, i) * stemp[l];
        element_linsolve_left<T, SIZE1>(size1, xx + i * sg, mm,
                                        qq); // solve 1X1 problem gtemp*mm=qq
    }
    for (l = 0; l <= ntau; l++)
        element_set<T, SIZE1>(size1, G.tvptr(n, ntau - l), xx + sg * l);
    delete[] stemp;
    delete[] htemp;
    delete[] gtemp;
    delete[] xx;
    delete[] mm;
    delete[] qq;
    delete[] one;
    return;
}
/// @private
template <typename T, int SIZE1>
void pseudodyson_start_tv(herm_pseudo<T> &G, T lam0, std::complex<T> *H0,
                          herm_pseudo<T> &Sigma, integration::Integrator<T> &I, T beta,
                          T h) {
    int n, k = I.k();
    for (n = 0; n <= k; n++)
        pseudodyson_timestep_tv<T, SIZE1>(n, G, lam0, H0, Sigma, I, beta, h);
}
// for matrix geenfunctions: H passed as Matrix
/// @private
template <typename T, class Matrix>
void pseudodyson_mat(herm_pseudo<T> &G, T lam0, Matrix &H0, herm_pseudo<T> &Sigma,
                     integration::Integrator<T> &I, T beta) {
    int size1 = G.size1(), l, m;
    std::complex<T> *h0 = new std::complex<T>[size1 * size1];
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert((int)H0.size1() == size1 && (int)H0.size2() == size1);
    for (l = 0; l < size1; l++)
        for (m = 0; m < size1; m++)
            h0[l * size1 + m] = H0(l, m);
    if (size1 == 1)
        pseudodyson_mat_dispatch<T, 1>(G, lam0, h0, Sigma, I, beta);
    else
        pseudodyson_mat_dispatch<T, LARGESIZE>(G, lam0, h0, Sigma, I, beta);
    delete[] h0;
}
/// @private
template <typename T, class Matrix>
void pseudodyson_start(herm_pseudo<T> &G, T lam0, Matrix &H0, std::vector<Matrix> &H,
                       herm_pseudo<T> &Sigma, integration::Integrator<T> &I, T beta, T h) {
    int size1 = G.size1(), l, m, n, k = I.k(), k1 = k + 1, sg = size1 * size1;
    std::complex<T> *hh = new std::complex<T>[k1 * size1 * size1];
    std::complex<T> *h0 = new std::complex<T>[size1 * size1];
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert(G.nt() >= k);
    assert(Sigma.nt() >= k);
    assert((int)H.size() >= k1);
    for (n = 0; n <= k; n++) {
        assert((int)H[n].size1() == size1 && (int)H[n].size2() == size1);
        for (l = 0; l < size1; l++)
            for (m = 0; m < size1; m++)
                hh[n * sg + l * size1 + m] = H[n](l, m);
    }
    assert((int)H0.size1() == size1 && (int)H0.size2() == size1);
    for (l = 0; l < size1; l++)
        for (m = 0; m < size1; m++)
            h0[l * size1 + m] = H0(l, m);
    if (size1 == 1) {
        dyson_start_ret<T, herm_pseudo<T>, 1>(G, lam0, hh, Sigma, I, h);
        pseudodyson_start_tv<T, 1>(G, lam0, h0, Sigma, I, beta, h);
        dyson_start_les<T, herm_pseudo<T>, 1>(G, lam0, hh, Sigma, I, beta, h);
    } else {
        dyson_start_ret<T, herm_pseudo<T>, LARGESIZE>(G, lam0, hh, Sigma, I, h);
        pseudodyson_start_tv<T, LARGESIZE>(G, lam0, h0, Sigma, I, beta, h);
        dyson_start_les<T, herm_pseudo<T>, LARGESIZE>(G, lam0, hh, Sigma, I, beta, h);
    }
    delete[] h0;
    delete[] hh;
}
/// @private
template <typename T, class Matrix>
void pseudodyson_timestep(int n, herm_pseudo<T> &G, T lam0, Matrix &H0,
                          std::vector<Matrix> &H, herm_pseudo<T> &Sigma,
                          integration::Integrator<T> &I, T beta, T h) {
    int size1 = G.size1(), l, m, j, k = I.k(), sg = size1 * size1;
    std::complex<T> *hh = new std::complex<T>[(n + 1) * size1 * size1];
    std::complex<T> *h0 = new std::complex<T>[size1 * size1];
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert(G.nt() >= n);
    assert(Sigma.nt() >= n);
    assert(n > k);
    assert((int)H.size() >= n + 1);
    for (j = 0; j <= n; j++) {
        assert((int)H[j].size1() == size1 && (int)H[j].size2() == size1);
        for (l = 0; l < size1; l++)
            for (m = 0; m < size1; m++)
                hh[j * sg + l * size1 + m] = H[j](l, m);
    }
    assert((int)H0.size1() == size1 && (int)H0.size2() == size1);
    for (l = 0; l < size1; l++)
        for (m = 0; m < size1; m++)
            h0[l * size1 + m] = H0(l, m);
    if (size1 == 1) {
        dyson_timestep_ret<T, herm_pseudo<T>, 1>(n, G, lam0, hh, Sigma, I, h);
        pseudodyson_timestep_tv<T, 1>(n, G, lam0, h0, Sigma, I, beta, h);
        dyson_timestep_les<T, herm_pseudo<T>, 1>(n, G, lam0, hh, Sigma, I, beta, h);
    } else {
        dyson_timestep_ret<T, herm_pseudo<T>, LARGESIZE>(n, G, lam0, hh, Sigma, I, h);
        pseudodyson_timestep_tv<T, LARGESIZE>(n, G, lam0, h0, Sigma, I, beta, h);
        dyson_timestep_les<T, herm_pseudo<T>, LARGESIZE>(n, G, lam0, hh, Sigma, I, beta, h);
    }
    delete[] hh;
    delete[] h0;
}
// with function object:
/// @private
template <typename T>
void pseudodyson_mat(herm_pseudo<T> &G, T lam0, function<T> &H0, herm_pseudo<T> &Sigma,
                     integration::Integrator<T> &I, T beta) {
    int size1 = G.size1();
    std::complex<T> *h0;
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert((int)H0.size1() == size1 && (int)H0.size2() == size1);
    h0 = H0.ptr(-1);
    if (size1 == 1)
        pseudodyson_mat_dispatch<T, 1>(G, lam0, h0, Sigma, I, beta);
    else
        pseudodyson_mat_dispatch<T, LARGESIZE>(G, lam0, h0, Sigma, I, beta);
}
/// @private
template <typename T>
void pseudodyson_start(herm_pseudo<T> &G, T lam0, function<T> &H, herm_pseudo<T> &Sigma,
                       integration::Integrator<T> &I, T beta, T h) {
    int size1 = G.size1(), k = I.k();
    std::complex<T> *hh, *h0;
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert(G.nt() >= k);
    assert(Sigma.nt() >= k);
    assert(H.nt() >= k);
    assert((int)H.size1() == size1 && (int)H.size2() == size1);
    hh = H.ptr(0);
    h0 = H.ptr(-1);
    if (size1 == 1) {
        dyson_start_ret<T, herm_pseudo<T>, 1>(G, lam0, hh, Sigma, I, h);
        pseudodyson_start_tv<T, 1>(G, lam0, h0, Sigma, I, beta, h);
        dyson_start_les<T, herm_pseudo<T>, 1>(G, lam0, hh, Sigma, I, beta, h);
    } else {
        dyson_start_ret<T, herm_pseudo<T>, LARGESIZE>(G, lam0, hh, Sigma, I, h);
        pseudodyson_start_tv<T, LARGESIZE>(G, lam0, h0, Sigma, I, beta, h);
        dyson_start_les<T, herm_pseudo<T>, LARGESIZE>(G, lam0, hh, Sigma, I, beta, h);
    }
}
/// @private
template <typename T>
void pseudodyson_timestep(int n, herm_pseudo<T> &G, T lam0, function<T> &H,
                          herm_pseudo<T> &Sigma, integration::Integrator<T> &I, T beta,
                          T h) {
    int size1 = G.size1(), k = I.k();
    std::complex<T> *hh, *h0;
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert(G.nt() >= n);
    assert(Sigma.nt() >= n);
    assert(n > k);
    assert(H.nt() >= n);
    assert((int)H.size1() == size1 && (int)H.size2() == size1);
    hh = H.ptr(0);
    h0 = H.ptr(-1);
    if (size1 == 1) {
        dyson_timestep_ret<T, herm_pseudo<T>, 1>(n, G, lam0, hh, Sigma, I, h);
        pseudodyson_timestep_tv<T, 1>(n, G, lam0, h0, Sigma, I, beta, h);
        dyson_timestep_les<T, herm_pseudo<T>, 1>(n, G, lam0, hh, Sigma, I, beta, h);
    } else {
        dyson_timestep_ret<T, herm_pseudo<T>, LARGESIZE>(n, G, lam0, hh, Sigma, I, h);
        pseudodyson_timestep_tv<T, LARGESIZE>(n, G, lam0, h0, Sigma, I, beta, h);
        dyson_timestep_les<T, herm_pseudo<T>, LARGESIZE>(n, G, lam0, hh, Sigma, I, beta, h);
    }
}

// the following work for size1=1 only
/// @private
template <typename T>
void pseudodyson_mat(herm_pseudo<T> &G, T lam0, std::complex<T> &H0, herm_pseudo<T> &Sigma,
                     integration::Integrator<T> &I, T beta) {
    int size1 = G.size1();
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    if (size1 == 1)
        pseudodyson_mat_dispatch<T, 1>(G, lam0, &H0, Sigma, I, beta);
    else
        pseudodyson_mat_dispatch<T, LARGESIZE>(G, lam0, &H0, Sigma, I, beta);
}
/// @private
template <typename T>
void pseudodyson_start(herm_pseudo<T> &G, T lam0, std::complex<T> H0,
                       std::vector<std::complex<T>> &H, herm_pseudo<T> &Sigma,
                       integration::Integrator<T> &I, T beta, T h) {
    int size1 = G.size1(), n, k = I.k(), k1 = k + 1;
    std::complex<T> *hh = new std::complex<T>[k1];
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert(G.nt() >= k);
    assert(Sigma.nt() >= k);
    assert(size1 == 1);
    assert((int)H.size() >= k1);
    for (n = 0; n <= k; n++)
        hh[n] = H[n];
    if (size1 == 1) {
        dyson_start_ret<T, herm_pseudo<T>, 1>(G, lam0, hh, Sigma, I, h);
        pseudodyson_start_tv<T, 1>(G, lam0, &H0, Sigma, I, beta, h);
        dyson_start_les<T, herm_pseudo<T>, 1>(G, lam0, hh, Sigma, I, beta, h);
    } else {
        dyson_start_ret<T, herm_pseudo<T>, LARGESIZE>(G, lam0, hh, Sigma, I, h);
        pseudodyson_start_tv<T, LARGESIZE>(G, lam0, &H0, Sigma, I, beta, h);
        dyson_start_les<T, herm_pseudo<T>, LARGESIZE>(G, lam0, hh, Sigma, I, beta, h);
    }
    delete[] hh;
}
/// @private
template <typename T>
void pseudodyson_timestep(int tstp, herm_pseudo<T> &G, T lam0, std::complex<T> H0,
                          std::vector<std::complex<T>> &H, herm_pseudo<T> &Sigma,
                          integration::Integrator<T> &I, T beta, T h) {
    int size1 = G.size1(), n, k = I.k();
    std::complex<T> *hh = new std::complex<T>[tstp + 1];
    assert(G.size1() == Sigma.size1());
    assert(G.ntau() == Sigma.ntau());
    assert(G.nt() >= tstp);
    assert(Sigma.nt() >= tstp);
    assert(size1 == 1);
    assert((int)H.size() >= tstp + 1);
    assert(tstp > k);
    for (n = 0; n <= tstp; n++)
        hh[n] = H[n];
    if (size1 == 1) {
        dyson_timestep_ret<T, herm_pseudo<T>, 1>(tstp, G, lam0, hh, Sigma, I, h);
        pseudodyson_timestep_tv<T, 1>(tstp, G, lam0, &H0, Sigma, I, beta, h);
        dyson_timestep_les<T, herm_pseudo<T>, 1>(tstp, G, lam0, hh, Sigma, I, beta, h);
    } else {
        dyson_timestep_ret<T, herm_pseudo<T>, LARGESIZE>(tstp, G, lam0, hh, Sigma, I, h);
        pseudodyson_timestep_tv<T, LARGESIZE>(tstp, G, lam0, &H0, Sigma, I, beta, h);
        dyson_timestep_les<T, herm_pseudo<T>, LARGESIZE>(tstp, G, lam0, hh, Sigma, I, beta,
                                                         h);
    }
    delete[] hh;
}

} // namespace cntr

#endif  // CNTR_PSEUDODYSON_IMPL_H
