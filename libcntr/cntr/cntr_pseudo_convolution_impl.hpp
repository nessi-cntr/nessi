#ifndef CNTR_PSEUDO_CONVOLUTION_IMPL_H
#define CNTR_PSEUDO_CONVOLUTION_IMPL_H

#include "cntr_pseudo_convolution_decl.hpp"
//#include "cntr_exception.hpp"
#include "cntr_elements.hpp"
#include "cntr_matsubara_decl.hpp"
#include "cntr_convolution_decl.hpp"
#include "cntr_herm_matrix_decl.hpp"

namespace cntr {

#define CPLX std::complex<T>

/* /////////////////////////////////////////////////////////////////////////////////////////
OVERVIEW:
ROUTINES FOR COMPUTING THE CONVOLUTION  C = gamma * C + alpha * A * ft * B  at one timestep

** ALSO PSEUDO_CONVOLUTION ROUTINES:
///////////////////////////////////////////////////////////////////////////////////////// */

////////////////////////////////////////////////////////////////////////////////////////////////////
/// possibly needed for pseudodyson_tv_timestep: restricted tau-integral:
//  Ctv(t,tau) += alpha*Aret*Btv(t,tau) - ii * alpha * int_t1^beta dt1 Atv(t,t1)*Bmat(t1,tau)
/// @private
template <typename T, class GG, int SIZE1>
void incr_pseudo_convolution_tv(int tstp, std::vector<bool> &mask, CPLX alpha, GG &C, GG &A,
                                GG &Acc, CPLX *f0, CPLX *ft, GG &B, GG &Bcc,
                                integration::Integrator<T> &I, T beta, T h) {
    int ntau = A.ntau();
    int size1 = A.size1();
    int sc = size1 * size1;
    bool func = (ft == NULL ? false : true);
    CPLX adt = alpha * h;
    CPLX adtau = alpha * beta * (1.0 / ntau);
    int kt = I.get_k();
    int n1 = (tstp >= kt ? tstp : kt);
    {
        // CONVOLUTION OF TV SECTION
        int j, m, n, saf = size1 * size1, sfb = size1 * size1;
        CPLX *ctemp1 = new CPLX[sc];
        CPLX *ctemp2 = new CPLX[sc];
        CPLX *bmat;
        CPLX *aret;
        CPLX *bmat1 = 0;
        CPLX *aret1 = 0;
        {
            // aret[j]=Aret(tstp,j)*f(j)   j=0 ... n1
            if (func) {
                CPLX *atemp = new CPLX[sc];
                aret1 = new CPLX[(n1 + 1) * saf];
                aret = aret1;
                for (j = 0; j <= tstp; j++) {
                    element_mult<T, SIZE1>(size1, aret1 + saf * j, A.retptr(tstp, j),
                                           ft + j * sc);
                }
                for (j = tstp + 1; j <= n1; j++) {
                    element_minusconj<T, SIZE1>(size1, atemp, Acc.retptr(j, tstp));
                    element_mult<T, SIZE1>(size1, aret1 + saf * j, atemp, ft + j * sc);
                }
                delete[] atemp;
            } else {
                if (n1 == tstp) {
                    aret = A.retptr(tstp, 0);
                } else {
                    aret1 = new CPLX[(n1 + 1) * saf];
                    aret = aret1;
                    for (j = 0; j <= tstp; j++)
                        element_set<T, SIZE1>(size1, aret1 + sc * j, A.retptr(tstp, j));
                    for (j = tstp + 1; j <= n1; j++)
                        element_minusconj<T, SIZE1>(size1, aret1 + sc * j,
                                                    Acc.retptr(j, tstp));
                }
            }
            // bmat[m]=f0*Bmat(m)
            if (func) {
                bmat1 = new CPLX[(ntau + 1) * sfb];
                bmat = bmat1;
                for (m = 0; m <= ntau; m++)
                    element_mult<T, SIZE1>(size1, bmat1 + m * sfb, f0, B.matptr(m));
            } else {
                bmat = B.matptr(0);
            }
        }
        for (m = 0; m <= ntau; m++) {
            if (mask[m]) {
                // CONTRIBUTION FROM Atv * Bmat : very similar to computing the matsubara
                // convolution
                // !!!!! Here is the only difference between pseudo_convolution and
                // convolution:
                // matsubara_integral_2_2 instead of matsubara_integral_2
                // NOTE: B.sig() really not needed for matsubara_integral_2_2,
                // other than for matsubara_integral_2
                matsubara_integral_2_2<T, SIZE1>(size1, m, ntau, ctemp1, A.tvptr(tstp, 0),
                                                 bmat, I);
                element_smul<T, SIZE1>(size1, ctemp1, adtau);
                // CONTRIBUTION FROM Aret * Btv
                element_set_zero<T, SIZE1>(size1, ctemp2);
                if (tstp <= 2 * kt + 2) {
                    for (n = 0; n <= n1; n++) {
                        element_incr<T, SIZE1>(size1, ctemp2, I.gregory_weights(tstp, n),
                                               aret + n * saf, B.tvptr(n, m));
                    }
                } else {
                    for (n = 0; n <= kt; n++) {
                        element_incr<T, SIZE1>(size1, ctemp2, I.gregory_omega(n),
                                               aret + n * saf, B.tvptr(n, m));
                    }
                    for (n = kt + 1; n < tstp - kt; n++) {
                        element_incr<T, SIZE1>(size1, ctemp2, aret + n * saf, B.tvptr(n, m));
                    }
                    for (n = tstp - kt; n <= tstp; n++) {
                        element_incr<T, SIZE1>(size1, ctemp2, I.gregory_omega(tstp - n),
                                               aret + n * saf, B.tvptr(n, m));
                    }
                }
                element_smul<T, SIZE1>(size1, ctemp2, adt);
                element_incr<T, SIZE1>(size1, C.tvptr(tstp, m), ctemp1);
                element_incr<T, SIZE1>(size1, C.tvptr(tstp, m), ctemp2);
            }
        }
        if (bmat1 != 0)
            delete[] bmat1;
        if (aret1 != 0)
            delete[] aret1;
        delete[] ctemp1;
        delete[] ctemp2;
    }
    return;
}
/// @private
template <typename T, class GG, int SIZE1>
void incr_pseudo_convolution_mat(std::vector<bool> &mask, CPLX alpha, GG &C, GG &A, CPLX *f0,
                                 GG &B, integration::Integrator<T> &I, T beta) {
    int ntau = A.ntau();
    int size1 = A.size1();
    int sc = size1 * size1;
    bool func = (f0 == NULL ? false : true);
    CPLX adtau = alpha * beta * (1.0 / ntau);
    {
        // CONVOLUTION OF MATSUBARA GREENFUNCTIONS
        int m, sfb = sc;
        CPLX *ctemp = new CPLX[sc], *bmat1 = 0;
        CPLX *btemp;
        if (func) {
            bmat1 = new CPLX[(ntau + 1) * sfb];
            for (m = 0; m <= ntau; m++) {
                element_mult<T, SIZE1>(size1, bmat1 + m * sfb, f0, B.matptr(m));
            }
            btemp = bmat1;
        } else {
            btemp = B.matptr(0);
        }
        // t1=omp_get_wtime();
        for (m = 0; m <= ntau; m++) {
            if (mask[m]) {
                // compute cmat(m*dtau) = int_0^tau dx amat(tau-x) f0 b(x)
                // NOTE: HERE IS THE ONLY DIFFERENCE TO THE USUAL MATSUBARA CONVOLUTION
                // matsubara_integral_1<T,SIZE1>(size1,m,ntau,ctemp,A.matptr(0),btemp,I,A.sig());
                matsubara_integral_1_1<T, SIZE1>(size1, m, ntau, ctemp, A.matptr(0), btemp,
                                                 I);
                element_incr<T, SIZE1>(size1, C.matptr(m), adtau, ctemp);
            }
        }
        // t2=omp_get_wtime();
        if (func)
            delete[] bmat1;
        delete[] ctemp;
    }
    // std::cout << "tid " << omp_get_thread_num() << " parallel time " << t2-t1  <<
    // std::endl;
    return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// pseudo-particle convolution: only mat and tv are different
//
///////////////////////////////////////////////////////////////////////////////////////////////////////
/// @private
template <typename T, class GG, int SIZE1>
void incr_pseudo_convolution(int tstp, CPLX alpha, GG &C, GG &A, GG &Acc, CPLX *f0, CPLX *ft,
                             GG &B, GG &Bcc, integration::Integrator<T> &I, T beta, T h) {
    int ntau = A.ntau();
    // this function is still not on top level, so no asserts!
    if (tstp == -1) {
        std::vector<bool> mask(ntau + 1, true);
        incr_pseudo_convolution_mat<T, GG, SIZE1>(mask, alpha, C, A, f0, B, I, beta);
    } else if (tstp >= 0) {
        std::vector<bool> mask_ret(tstp + 1, true);
        std::vector<bool> mask_tv(ntau + 1, true);
        std::vector<bool> mask_les(tstp + 1, true);
        incr_convolution_ret<T, GG, SIZE1>(tstp, mask_ret, alpha, C, A, Acc, ft, B, Bcc, I,
                                           h);
        incr_pseudo_convolution_tv<T, GG, SIZE1>(tstp, mask_tv, alpha, C, A, Acc, f0, ft, B,
                                                 Bcc, I, beta, h);
        incr_convolution_les<T, GG, SIZE1>(tstp, mask_les, alpha, C, A, Acc, f0, ft, B, Bcc,
                                           I, beta, h);
    }
    return;
}
/// @private
template <typename T>
void pseudo_convolution_timestep(int tstp, herm_pseudo<T> &C, herm_pseudo<T> &A,
                                 herm_pseudo<T> &Acc, function<T> &ft, herm_pseudo<T> &B,
                                 herm_pseudo<T> &Bcc, integration::Integrator<T> &I, T beta,
                                 T h) {
    int kt = I.k();
    int ntmin = (tstp == -1 || tstp > kt ? tstp : kt);
    if (tstp < -1)
        return;
    int size1 = A.size1();
    std::complex<T> *fttemp;
    assert(C.size1()==size1);
    assert(ft.size1()==size1);
    assert(B.size1()==size1);
    assert(Bcc.size1()==size1);
    assert(Acc.size1()==size1);
    assert(kt>=0 && kt<=5);
    assert(C.ntau()>=kt);
    assert(C.nt()>=ntmin);
    assert(A.nt()>=ntmin);
    assert(Acc.nt()>=ntmin);
    assert(B.nt()>=ntmin);
    assert(Bcc.nt()>=ntmin);
    assert(ft.nt()>=ntmin);
    assert(C.ntau()==A.ntau());
    assert(C.ntau()==Acc.ntau());
    assert(C.ntau()==B.ntau());
    assert(C.ntau()==Bcc.ntau());
    C.set_timestep_zero(tstp);
    fttemp = (tstp == -1 ? ft.ptr(-1) : ft.ptr(0));
    if (size1 == 1) {
        incr_pseudo_convolution<T, herm_pseudo<T>, 1>(tstp, CPLX(1, 0), C, A, Acc,
                                                      ft.ptr(-1), fttemp, B, Bcc,
                                                      integration::I<T>(kt), beta, h);
    } else {
        incr_pseudo_convolution<T, herm_pseudo<T>, LARGESIZE>(
            tstp, CPLX(1, 0), C, A, Acc, ft.ptr(-1), fttemp, B, Bcc, integration::I<T>(kt),
            beta, h);
    }
}
/// @private
template <typename T>
void pseudo_convolution_timestep(int n, herm_pseudo<T> &C, herm_pseudo<T> &A,
                                 function<T> &ft, herm_pseudo<T> &B,
                                 integration::Integrator<T> &I, T beta, T h) {
    pseudo_convolution_timestep<T>(n, C, A, A, ft, B, B, I, beta, h);
}
/// @private
template <typename T>
void pseudo_convolution_matsubara(herm_pseudo<T> &C, herm_pseudo<T> &A, herm_pseudo<T> &B,
                                  integration::Integrator<T> &I, T beta) {
    pseudo_convolution_timestep<T>(-1, C, A, A, B, B, I, beta, 0.0);
}
/// @private
template <typename T>
void pseudo_convolution_matsubara(herm_pseudo<T> &C, herm_pseudo<T> &A, function<T> &ft,
                                  herm_pseudo<T> &B, integration::Integrator<T> &I, T beta) {
    pseudo_convolution_timestep<T>(-1, C, A, A, ft, B, B, I, beta, 0.0);
}
/// @private
template <typename T>
void pseudo_convolution(herm_pseudo<T> &C, herm_pseudo<T> &A, herm_pseudo<T> &Acc,
                        function<T> &ft, herm_pseudo<T> &B, herm_pseudo<T> &Bcc,
                        integration::Integrator<T> &I, T beta, T h) {
    int tstp;
    for (tstp = -1; tstp <= C.nt(); tstp++)
        pseudo_convolution_timestep<T>(tstp, C, A, Acc, ft, B, Bcc, I, beta, h);
}
/// @private
template <typename T>
void pseudo_convolution_timestep(int tstp, herm_pseudo<T> &C, herm_pseudo<T> &A,
                                 herm_pseudo<T> &Acc, herm_pseudo<T> &B, herm_pseudo<T> &Bcc,
                                 integration::Integrator<T> &I, T beta, T h) {
    int kt = I.k();
    int ntmin = (tstp == -1 || tstp > kt ? tstp : kt);
    if (tstp < -1)
        return;
    int size1 = A.size1();
    assert(C.size1()==size1);
    assert(B.size1()==size1);
    assert(Bcc.size1()==size1);
    assert(Acc.size1()==size1);
    assert(kt>=0 && kt <=5);
    assert(C.ntau()>=kt);
    assert(C.nt()>=ntmin);
    assert(A.nt()>=ntmin);
    assert(Acc.nt()>=ntmin);
    assert(B.nt()>=ntmin);
    assert(Bcc.nt()>=ntmin);
    assert(C.ntau()==A.ntau());
    assert(C.ntau()==Acc.ntau());
    assert(C.ntau()==B.ntau());
    assert(C.ntau()==Bcc.ntau());
    C.set_timestep_zero(tstp);
    if (size1 == 1) {
        incr_pseudo_convolution<T, herm_pseudo<T>, 1>(
            tstp, CPLX(1, 0), C, A, Acc, NULL, NULL, B, Bcc, integration::I<T>(kt), beta, h);
    } else {
        incr_pseudo_convolution<T, herm_pseudo<T>, LARGESIZE>(
            tstp, CPLX(1, 0), C, A, Acc, NULL, NULL, B, Bcc, integration::I<T>(kt), beta, h);
    }
}
/// @private
template <typename T>
void pseudo_convolution_timestep(int n, herm_pseudo<T> &C, herm_pseudo<T> &A,
                                 herm_pseudo<T> &B, integration::Integrator<T> &I, T beta,
                                 T h) {
    pseudo_convolution_timestep<T>(n, C, A, A, B, B, I, beta, h);
}
/// @private
template <typename T>
void pseudo_convolution(herm_pseudo<T> &C, herm_pseudo<T> &A, herm_pseudo<T> &Acc,
                        herm_pseudo<T> &B, herm_pseudo<T> &Bcc,
                        integration::Integrator<T> &I, T beta, T h) {
    int tstp;
    for (tstp = -1; tstp <= C.nt(); tstp++)
        pseudo_convolution_timestep<T>(tstp, C, A, Acc, B, Bcc, I, beta, h);
}

#undef CPLX

}  // namespace cntr

#endif  // CNTR_PSEUDO_CONVOLUTION_IMPL_H
