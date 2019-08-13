#ifndef CNTR_PSEUDO_VIE2_IMPL_H
#define CNTR_PSEUDO_VIE2_IMPL_H

#include "cntr_pseudo_vie2_decl.hpp"

namespace cntr {

#define CPLX std::complex<T>
/////////////////////////////////////////////////////////////////////////////////////////
// VIE2 FOR PSEUDO-PARTICLE GREEN'S FUNCTIONS:
// ONLY CONCVOLUTION IS DIFFERENT FROM USUAL VIE2
// IMPLEMENTATION OF PSEUDO_VIE2 BASED ON THE NEW (OMP PARALELLIZED VIE2 SOLVER)
// Matsubara: implementation copied from pseudodyson_mat_dispatch, without derivative term
/// @private
template <typename T, int SIZE1>
void pseudo_vie2_mat_dispatch(herm_pseudo<T> &G, herm_pseudo<T> &F, herm_pseudo<T> &Fcc,
                              herm_pseudo<T> &Q, integration::Integrator<T> &I, T beta) {
    typedef std::complex<T> cplx;
    int ntau, k = I.get_k(), k1 = k + 1, sg, l, m, p, q, n, j, size1 = G.size1();
    T dtau;
    cplx *mm, *qq, *gtemp, *one, *stemp, cweight;

    ntau = G.ntau();
    assert(ntau > k);
    assert(F.ntau() == ntau);
    assert(Fcc.ntau() == ntau);
    sg = G.element_size();
    dtau = beta / ntau;

    mm = new cplx[k * k * sg];
    qq = new cplx[k * sg];
    gtemp = new cplx[k * sg];
    stemp = new cplx[k * sg];
    one = new cplx[sg];
    element_set<T, SIZE1>(size1, one, 1.0);
    // INITIAL VALUE : GMAT=Q
    for (l = 0; l < sg; l++)
        G.matptr(0)[l] = Q.matptr(0)[l]; //-one[l];
    // THE FIRST K STEPS Gmat(n) for n=1...k
    for (l = 0; l < k * k * sg; l++)
        mm[l] = 0;
    for (l = 0; l < k * sg; l++)
        qq[l] = 0;
    for (n = 1; n <= k; n++) {
        p = n - 1;
        element_set<T, SIZE1>(size1, qq + p * sg, Q.matptr(n));
    }
    // derive linear equations for gtemp(p)=Gmat(p+1)
    // mm(p,q)*G(q)=Q(p) for p,q=0...k-1
    for (n = 1; n <= k; n++) {
        p = n - 1;
        // "1" -- goes into m(p,p)
        for (l = 0; l < sg; l++)
            mm[sg * (p * k + p) + l] += one[l];
        // integral int_0^n dx F(n*dtau-x)G(x) :
        if (n < k) {
            // need the rcorr for n<k:
            // Zitat: R=int_0^n dx a(n-x) b(x) = sum_{m,j=0}^k I.rcorr(n,m,j) a(m)*b(j)
            for (m = 0; m <= k; m++) {
                element_set<T, SIZE1>(size1, stemp, F.matptr(m));
                for (j = 0; j <= k; j++) {
                    cweight = dtau * I.rcorr(n, m, j);
                    if (j == 0) { // goes into qq(p)
                        element_incr<T, SIZE1>(size1, qq + p * sg, -cweight, stemp,
                                               G.matptr(0));
                    } else { // goes into mm(p,q)
                        q = j - 1;
                        for (l = 0; l < sg; l++)
                            mm[sg * (p * k + q) + l] += stemp[l] * cweight;
                    }
                }
            }
        } else {
            // usual integral for n=k
            // R=int_0^k dx a(n-x) b(x) = sum_{m=0}^k I.gregory_weights(k,m) a(k-m)*b(m)
            for (m = 0; m <= k; m++) {
                cweight = dtau * I.gregory_weights(k, m);
                element_set<T, SIZE1>(size1, stemp, F.matptr(k - m));
                if (m == 0) { // goes into qq(p)
                    element_incr<T, SIZE1>(size1, qq + p * sg, -cweight, stemp, G.matptr(0));
                } else { // goes into mm(p,q)
                    q = m - 1;
                    for (l = 0; l < sg; l++)
                        mm[sg * (p * k + q) + l] += stemp[l] * cweight;
                }
            }
        }
    }
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
        // element_set<T,SIZE1>(size1,qq,Q.matptr(n));
        element_set_zero<T, SIZE1>(size1, qq);
        // integral
        if (n < 2 * k1) {
            for (m = 0; m < n; m++) {
                cweight = I.gregory_weights(n, m);
                element_incr<T, SIZE1>(size1, qq, -cweight, F.matptr(n - m), G.matptr(m));
            }
        } else { // some weights are 1
            for (m = 0; m <= k; m++) {
                cweight = I.gregory_weights(n, m);
                element_incr<T, SIZE1>(size1, qq, -cweight, F.matptr(n - m), G.matptr(m));
            }
            for (m = k1; m < n - k; m++) {
                cweight = 1.0;
                element_incr<T, SIZE1>(size1, qq, -cweight, F.matptr(n - m), G.matptr(m));
                // element_incr<T,SIZE1>(size1,qq,F.matptr(n-m),G.matptr(m));
            }
            for (m = n - k; m < n; m++) {
                cweight = I.gregory_omega(n - m);
                element_incr<T, SIZE1>(size1, qq, -cweight, F.matptr(n - m), G.matptr(m));
            }
        }
        element_smul<T, SIZE1>(size1, qq, dtau);
        element_incr<T, SIZE1>(size1, qq, Q.matptr(n));
        // determine mm
        element_set<T, SIZE1>(size1, stemp, F.matptr(0));
        for (l = 0; l < sg; l++)
            mm[l] = one[l] + dtau * I.gregory_weights(n, n) * stemp[l];
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
template <typename T>
void pseudo_vie2_mat(herm_pseudo<T> &G, herm_pseudo<T> &F, herm_pseudo<T> &Fcc,
                     herm_pseudo<T> &Q, integration::Integrator<T> &I, T beta) {
    int size1 = G.size1();
    assert(G.size1() == F.size1());
    assert(G.ntau() == F.ntau());
    assert(G.size1() == Fcc.size1());
    assert(G.ntau() == Fcc.ntau());
    if (size1 == 1) {
        pseudo_vie2_mat_dispatch<T, 1>(G, F, Fcc, Q, I, beta);
    } else {
        pseudo_vie2_mat_dispatch<T, LARGESIZE>(G, F, Fcc, Q, I, beta);
    }
}
// TV: for starting, I just copy the pseudodyson_tv routine (without derivative term,
// lam0->1, integral + => -)
/// @private
template <typename T, int SIZE1>
void pseudo_vie2_timestep_tv(int n, herm_pseudo<T> &G, herm_pseudo<T> &F,
                             herm_pseudo<T> &Fcc, herm_pseudo<T> &Q,
                             integration::Integrator<T> &I, T beta, T h) {
    // tv: solve conjugate equation: (quite analogous to matsubara)
    // Gtv(n,m) = Qtv(n,m) - int_{ntau-m}^ntau dx Gtv(n,x) Fccmat(x-m)  - int_0^n dj
    // Gret(n,j) Fcctv(j,m)
    // start this equation from beta and integrate to 0: m -> l=beta-m, X(l)=Gtv(n,beta-l)
    // X(l) = Qtv(n,beta-l) - int_0^l dj X(j) Fccmat(l-j)  - int_0^n dj Gret(n,j)
    // Fcctv(j,beta-l)
    typedef std::complex<T> cplx;
    int ntau, k = I.get_k(), k1 = k + 1, sg, l, m, p, q, i, j, n1, sig, size1 = G.size1();
    T dtau;
    cplx *mm, *qq, *gtemp, *one, *stemp, *xx, cweight;

    n1 = (n > k ? n : k);
    ntau = G.ntau();
    assert(ntau > k);
    assert(F.ntau() == ntau);
    assert(G.nt() >= n1);
    assert(F.nt() >= n1);
    sig = G.sig();
    assert(sig == F.sig());
    sg = G.element_size();
    dtau = beta / ntau;

    mm = new cplx[k * k * sg];
    qq = new cplx[k * sg];
    xx = new cplx[(ntau + 1) * sg];
    stemp = new cplx[k * sg];
    gtemp = new cplx[sg];
    one = new cplx[sg];
    element_set<T, SIZE1>(size1, one, 1.0);
    // for(i=0;i<sg*(ntau+1);i++) xx[i]=0;
    // Qtv(n,beta-i) -int_0^n dj Gret(n,j) Fcctv(j,beta-i)  --->  X(i) for i=1..ntau
    for (i = 1; i <= ntau; i++) {
        // element_incr<T,SIZE1>(size1,xx+i*sg,gtemp,Q.tvptr(n,ntau-i));
        element_set<T, SIZE1>(size1, xx + i * sg, Q.tvptr(n, ntau - i));
    }
    for (j = 0; j <= n1; j++) {
        if (j <= n) {
            element_set<T, SIZE1>(size1, gtemp, G.retptr(n, j));
        } else {
            element_conj<T, SIZE1>(size1, gtemp, G.retptr(j, n));
            element_smul<T, SIZE1>(size1, gtemp, -1);
        }
        element_smul<T, SIZE1>(size1, gtemp, -I.gregory_weights(n, j) * h);
        for (i = 1; i <= ntau; i++) {
            element_incr<T, SIZE1>(size1, xx + i * sg, gtemp, Fcc.tvptr(j, ntau - i));
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
        for (l = 0; l < sg; l++)
            mm[sg * (p * k + p) + l] += one[l];
        // integral int_0^i dx X(x) Fccmat(i-x):
        if (i < k) {
            // need the rcorr:
            // Zitat: R=int_0^i dx a(i-x) b(x) = sum_{m,j=0}^k I.rcorr(i,m,j) a(m)*b(j)
            for (m = 0; m <= k; m++) {
                element_set<T, SIZE1>(size1, stemp, Fcc.matptr(m));
                for (j = 0; j <= k; j++) {
                    cweight = dtau * I.rcorr(i, m, j);
                    if (j == 0) { // goes into qq(p)
                        element_incr<T, SIZE1>(size1, qq + p * sg, -cweight, xx + 0 * sg,
                                               stemp);
                    } else { // goes into mm(p,q)
                        q = j - 1;
                        for (l = 0; l < sg; l++)
                            mm[sg * (p * k + q) + l] += stemp[l] * cweight;
                    }
                }
            }
        } else { // i=k
                 // usual integral for n=k
            // R=int_0^k dx a(n-x) b(x) = sum_{m=0}^k I.gregory_weights(k,m) a(k-m)*b(m)
            for (m = 0; m <= k; m++) {
                element_set<T, SIZE1>(size1, stemp, Fcc.matptr(k - m));
                cweight = dtau * I.gregory_weights(i, m);
                if (m == 0) { // goes into qq(p)
                    element_incr<T, SIZE1>(size1, qq + p * sg, -cweight, xx + 0 * sg, stemp);
                } else { // goes into mm(p,q)
                    q = m - 1;
                    for (l = 0; l < sg; l++)
                        mm[sg * (p * k + q) + l] += stemp[l] * cweight;
                }
            }
        }
    }
    // if(n==5) std::cout << " qq= " << qq[0] << " " << qq[1] << " " << "..." << std::endl;
    // if(n==5) std::cout << " mm= " << mm[0] << " " << mm[1] << " " << "..." << std::endl;
    element_linsolve_left<T, SIZE1>(size1, k, xx + 1 * sg, mm,
                                    qq); // solve kXk problem xx*mm=qq
    // if(n==5) std::cout << " xx= " << xx[0] << " " << xx[1] << " " << "..." << std::endl;
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
                element_incr<T, SIZE1>(size1, qq, -cweight, xx + m * sg, Fcc.matptr(i - m));
            }
        } else { // some weights are 1
            for (m = 0; m <= k; m++) {
                cweight = I.gregory_weights(i, m);
                element_incr<T, SIZE1>(size1, qq, -cweight, xx + m * sg, Fcc.matptr(i - m));
            }
            for (m = k1; m < i - k; m++) {
                cweight = 1.0;
                element_incr<T, SIZE1>(size1, qq, -cweight, xx + m * sg, Fcc.matptr(i - m));
            }
            for (m = i - k; m < i; m++) {
                cweight = I.gregory_omega(i - m);
                element_incr<T, SIZE1>(size1, qq, -cweight, xx + m * sg, Fcc.matptr(i - m));
            }
        }
        element_smul<T, SIZE1>(size1, qq, dtau);
        element_incr<T, SIZE1>(size1, qq, xx + i * sg);
        // determine mm
        element_set<T, SIZE1>(size1, stemp, Fcc.matptr(0));
        for (l = 0; l < sg; l++)
            mm[l] = one[l] + dtau * I.gregory_weights(i, i) * stemp[l];
        element_linsolve_left<T, SIZE1>(size1, xx + i * sg, mm,
                                        qq); // solve 1X1 problem gtemp*mm=qq
    }
    for (l = 0; l <= ntau; l++)
        element_set<T, SIZE1>(size1, G.tvptr(n, ntau - l), xx + sg * l);
    delete[] stemp;
    delete[] gtemp;
    delete[] xx;
    delete[] mm;
    delete[] qq;
    delete[] one;
    return;
}
/// @private
template <typename T, int SIZE1>
void pseudo_vie2_start_tv(herm_pseudo<T> &G, herm_pseudo<T> &F, herm_pseudo<T> &Fcc,
                          herm_pseudo<T> &Q, integration::Integrator<T> &I, T beta, T h) {
    int n, k = I.k();
    for (n = 0; n <= k; n++)
        pseudo_vie2_timestep_tv<T, SIZE1>(n, G, F, Fcc, Q, I, beta, h);
}
/// @private
template <typename T>
void pseudo_vie2_start(herm_pseudo<T> &G, herm_pseudo<T> &F, herm_pseudo<T> &Fcc,
                       herm_pseudo<T> &Q, integration::Integrator<T> &I, T beta, T h) {
    int size1 = G.size1(), k = I.k();
    assert(G.size1() == F.size1());
    assert(G.ntau() == F.ntau());
    assert(G.nt() >= k);
    assert(F.nt() >= k);
    assert(G.size1() == Fcc.size1());
    assert(G.ntau() == Fcc.ntau());
    assert(Fcc.nt() >= k);
    if (size1 == 1) {
        vie2_start_ret<T, herm_pseudo<T>, 1>(
            G, F, Fcc, Q, I, h); // same as usual vie2, because convolution is the same
        pseudo_vie2_start_tv<T, 1>(G, F, Fcc, Q, I, beta, h);
        vie2_start_les<T, herm_pseudo<T>, 1>(G, F, Fcc, Q, I, beta, h);
    } else {
        vie2_start_ret<T, herm_pseudo<T>, LARGESIZE>(G, F, Fcc, Q, I, h);
        pseudo_vie2_start_tv<T, LARGESIZE>(G, F, Fcc, Q, I, beta, h);
        vie2_start_les<T, herm_pseudo<T>, LARGESIZE>(G, F, Fcc, Q, I, beta, h);
    }
}
/// @private
template <typename T, class GG, int SIZE1>
void pseudo_vie2_timestep_new_dispatch(int tstp, GG &B, CPLX alpha, GG &A, GG &Acc, CPLX *f0,
                                       CPLX *ft, GG &Q, integration::Integrator<T> &I,
                                       T beta, T h) {
    int kt = I.get_k();
    int ntau = A.ntau();
    int size1 = A.size1();
    int sc = size1 * size1;
    bool func = (ft == NULL ? false : true);
    int n, m;
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
        for (n = 0; n <= n1; n++)
            element_conj<T, SIZE1>(size1, ftcc + n * sc, ft + n * sc);
    } else {
        ftcc = NULL;
        f0cc = NULL;
    }
    {
        std::vector<bool> mask_ret(tstp + 1, true);
        std::vector<bool> mask_les(tstp + 1, true);
        std::vector<bool> mask_tv(ntau + 1, true);
        CPLX *mtemp = new CPLX[sc];
        CPLX *one = new CPLX[sc];
        CPLX *mm = new CPLX[sc];
        T wt;
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
        // retarded part:
        // Q1 = Q - alpha*A*ft*B, with Bret(tstp,n)=0
        incr_convolution_ret<T, GG, SIZE1>(tstp, mask_ret, -alpha, Q, A, Acc, ft, B, B, I,
                                           h);
        for (n = 0; n <= tstp; n++) {
            wt = (tstp < kt ? I.poly_integration(n, tstp, tstp)
                            : I.gregory_weights(tstp - n, 0));
            element_set<T, SIZE1>(size1, mm, one);
            element_incr<T, SIZE1>(size1, mm, wt, mtemp);
            element_linsolve_right<T, SIZE1>(size1, 1, B.retptr(tstp, n), mm,
                                             Q.retptr(tstp, n));
        }
        // tv part:
        // Q -> Q - alpha*A*ft*B, with Btv(tstp,n)=0
        // incr_convolution_tv<T,GG,SIZE1>(tstp,mask_tv,-alpha,Q,A,Acc,f0,ft,B,B,I,beta,h);
        incr_pseudo_convolution_tv<T, GG, SIZE1>(tstp, mask_tv, -alpha, Q, A, Acc, f0, ft, B,
                                                 B, I, beta, h);
        for (m = 0; m <= ntau; m++) {
            wt = I.gregory_weights(tstp, tstp);
            element_set<T, SIZE1>(size1, mm, one);
            element_incr<T, SIZE1>(size1, mm, wt, mtemp);
            element_linsolve_right<T, SIZE1>(size1, 1, B.tvptr(tstp, m), mm,
                                             Q.tvptr(tstp, m));
        }
        // les part:
        // Q -> Q - conj(alpha)*B*ftcc*Acc, with Bles(n,tstp)=0, n=0...tstp-1 (Bret already
        // known!)
        mask_les[tstp] = false;
        incr_convolution_les<T, GG, SIZE1>(tstp, mask_les, -conj(alpha), Q, B, B, f0cc, ftcc,
                                           Acc, A, I, beta, h);
        for (n = 0; n < tstp; n++) {
            wt = I.gregory_weights(tstp, tstp);
            element_set<T, SIZE1>(size1, mm, one);
            element_incr<T, SIZE1>(size1, mm, wt, mtemp);
            element_conj<T, SIZE1>(size1, mm);
            element_linsolve_left<T, SIZE1>(size1, 1, B.lesptr(n, tstp), mm,
                                            Q.lesptr(n, tstp));
        }
        // Bles(tstp,tstp) (Bles(n<tstp,tstp) etc. enters the convolution!)
        n = tstp;
        mask_les = std::vector<bool>(tstp + 1, false);
        mask_les[n] = true;
        incr_convolution_les<T, GG, SIZE1>(tstp, mask_les, -conj(alpha), Q, B, B, f0cc, ftcc,
                                           Acc, A, I, beta, h);
        wt = I.gregory_weights(tstp, tstp);
        element_set<T, SIZE1>(size1, mm, one);
        element_incr<T, SIZE1>(size1, mm, wt, mtemp);
        element_conj<T, SIZE1>(size1, mm);
        element_linsolve_left<T, SIZE1>(size1, 1, B.lesptr(n, tstp), mm, Q.lesptr(n, tstp));
        delete[] mm;
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
/// @private
template <typename T>
void pseudo_vie2_timestep(int tstp, herm_pseudo<T> &G, herm_pseudo<T> &F,
                          herm_pseudo<T> &Fcc, herm_pseudo<T> &Q,
                          integration::Integrator<T> &I, T beta, T h) {
    int ntau = G.ntau();
    int size1 = G.size1(), kt = I.k();
    int n1 = (tstp >= kt ? tstp : kt);
    assert(tstp >= 0);
    assert(ntau > 0);
    assert(kt>=0 && kt<=5);
    assert(kt<= 2 * ntau + 2);
    assert(ntau==Fcc.ntau());
    assert(ntau==F.ntau());
    assert(ntau==Q.ntau());
    assert(F.sig()==G.sig());
    assert(Fcc.sig()==G.sig());
    assert(Q.sig()==G.sig());
    assert(F.size1()==size1);
    assert(F.size2()==size1);
    assert(Fcc.size1()==size1);
    assert(Fcc.size2()==size1);
    assert(Q.size1()==size1);
    assert(Q.size2()==size1);
    assert(F.nt()>=n1);
    assert(Fcc.nt()>=n1);
    assert(G.nt()>=n1);
    assert(Q.nt()>=n1);

    if (size1 == 1) {
        pseudo_vie2_timestep_new_dispatch<T, herm_pseudo<T>, 1>(tstp, G, CPLX(1, 0), F, Fcc,
                                                                NULL, NULL, Q, I, beta, h);
    } else {
        pseudo_vie2_timestep_new_dispatch<T, herm_pseudo<T>, LARGESIZE>(
            tstp, G, CPLX(1, 0), F, Fcc, NULL, NULL, Q, I, beta, h);
    }
}
/// @private
template <typename T>
void pseudo_vie2(herm_pseudo<T> &G, herm_pseudo<T> &F, herm_pseudo<T> &Fcc,
                 herm_pseudo<T> &Q, integration::Integrator<T> &I, T beta, T h) {
    int tstp, k = I.k();
    pseudo_vie2_mat(G, F, Fcc, Q, I, beta);
    if (G.nt() >= 0)
        pseudo_vie2_start(G, F, Fcc, Q, I, beta, h);
    for (tstp = k + 1; tstp <= G.nt(); tstp++)
        pseudo_vie2_timestep(tstp, G, F, Fcc, Q, I, beta, h);
}
/// @private
#if CNTR_USE_OMP == 1
/// @private
template <typename T, class GG, int SIZE1>
void pseudo_vie2_timestep_omp_dispatch(int omp_num_threads, int tstp, GG &B, CPLX alpha,
                                       GG &A, GG &Acc, CPLX *f0, CPLX *ft, GG &Q,
                                       integration::Integrator<T> &I, T beta, T h) {
    /// identical to vie2_timestep_omp_dispatch, but with pseudo-convolution:
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
            // !!!!PSEUDO-CONVOLUTION --- only different line to vie2_timestep_omp_dispatch
            // !!!!!
            // incr_convolution_tv<T,GG,SIZE1>(tstp,mask_tv,-alpha,Q,A,Acc,f0,ft,B,B,I,beta,h);
            incr_pseudo_convolution_tv<T, GG, SIZE1>(tstp, mask_tv, -alpha, Q, A, Acc, f0,
                                                     ft, B, B, I, beta, h);
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
/// @private
template <typename T>
void pseudo_vie2_timestep_omp(int omp_num_threads, int tstp, herm_pseudo<T> &G,
                              herm_pseudo<T> &F, herm_pseudo<T> &Fcc, herm_pseudo<T> &Q,
                              integration::Integrator<T> &I, T beta, T h) {
    int ntau = G.ntau();
    int size1 = G.size1(), kt = I.k();
    int n1 = (tstp >= kt ? tstp : kt);
    assert(tstp >= 0);
    assert(ntau > 0);
    assert(kt>=0 && kt<=5);
    assert(kt<= 2 * ntau + 2);
    assert(ntau==Fcc.ntau());
    assert(ntau== F.ntau());
    assert(ntau== Q.ntau());
    assert(F.sig()== G.sig());
    assert(Fcc.sig()== G.sig());
    assert(Q.sig()== G.sig());
    assert(F.size1()== size1);
    assert(F.size2()== size1);
    assert(Fcc.size1()== size1);
    assert(Fcc.size2()== size1);
    assert(Q.size1()== size1);
    assert(Q.size2()== size1);
    assert(F.nt()>=n1);
    assert(Fcc.nt()>=n1);
    assert(G.nt()>=n1);
    assert(Q.nt()>=n1);

    if (size1 == 1) {
        pseudo_vie2_timestep_omp_dispatch<T, herm_pseudo<T>, 1>(
            omp_num_threads, tstp, G, CPLX(1, 0), F, Fcc, NULL, NULL, Q, I, beta, h);
    } else {
        pseudo_vie2_timestep_omp_dispatch<T, herm_pseudo<T>, LARGESIZE>(
            omp_num_threads, tstp, G, CPLX(1, 0), F, Fcc, NULL, NULL, Q, I, beta, h);
    }
}
/// @private
template <typename T>
void pseudo_vie2_omp(int omp_num_threads, herm_pseudo<T> &G, herm_pseudo<T> &F,
                     herm_pseudo<T> &Fcc, herm_pseudo<T> &Q, integration::Integrator<T> &I,
                     T beta, T h) {
    int tstp, k = I.k();
    pseudo_vie2_mat(G, F, Fcc, Q, I, beta);
    // vie2_mat(G,F,Fcc,Q,beta,3);
    if (G.nt() >= 0)
        pseudo_vie2_start(G, F, Fcc, Q, I, beta, h);
    for (tstp = k + 1; tstp <= G.nt(); tstp++)
        pseudo_vie2_timestep_omp(omp_num_threads, tstp, G, F, Fcc, Q, I, beta, h);
}
#endif // CNTR_USE_OMP

#undef CPLX

}  // namespace cntr

#endif  // CNTR_PSEUDO_VIE2_IMPL_H
