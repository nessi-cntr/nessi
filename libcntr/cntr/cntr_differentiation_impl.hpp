#ifndef CNTR_DIFFERENTIATION_IMPL_H
#define CNTR_DIFFERENTIATION_IMPL_H

#include "cntr_differentiation_decl.hpp"
//#include "cntr_exception.hpp"
#include "cntr_herm_matrix_decl.hpp"

namespace cntr {

/*###########################################################################################
#
#   DIFFERENTIATION  dA(t,t') = id/dt A(t,t')   or  dA(t,t') = -id/dt' A(t,t')
#
#   ...  currently only implemented for herm_matrix of size1=1 -> solve using ptr
#   ...  Acc is the conjugate function to A (similar as in convolution).
#        Used only in for ret and les. If A is hermitian, simply provide Acc=A.
#   ...  For the computation of timestep tstp, A(t,t') is addressed at t,t' <= max(n,k),
#        where k is the Integartion order (see Integrator). I.e., tstp=0..k can be computed
#        only if A is given for t,t'<=k.
#
#
###########################################################################################*/
// left derivative (1)
/// @private
template <typename T>
void deriv1_timestep(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A, herm_matrix<T> &Acc,
                     integration::Integrator<T> &I, T beta, T h) {
    int tstp2, nt = A.nt();
    if (tstp == -1) {
        deriv1_matsubara(dA, A, I, beta);
    } // matsubara
    else {
        for (tstp2 = 0; tstp2 <= nt; tstp2++)
            deriv1_element(tstp, tstp2, dA, A, Acc, I, h); // ret and less
        deriv1_tv(tstp, dA, A, I, h); // tv
    }
}

// right derivative (2)
/// @private
template <typename T>
void deriv2_timestep(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A, herm_matrix<T> &Acc,
                     integration::Integrator<T> &I, T beta, T h) {
    int tstp2, nt = A.nt();
    if (tstp == -1) {
        deriv1_matsubara(dA, A, I, beta);
    } // matsubara: deriv1=deriv2
    else {
        for (tstp2 = 0; tstp2 <= nt; tstp2++)
            deriv2_element(tstp, tstp2, dA, A, Acc, I, h); // ret and less
        deriv2_tv(tstp, dA, A, I, beta); // tv
    }
}

// matsubara: left and right derivative equal
// deriv1 = -d/dtau
// deriv2 = +d/dtau' A(tau-tau') = -d/dtau A(tau-tau') = deriv1
// computed for all imaginary time directly
/// @private
template <typename T>
void deriv1_matsubara(herm_matrix<T> &dA, herm_matrix<T> &A, integration::Integrator<T> &I,
                      T beta) {
    typedef std::complex<T> cplx;
    int tstp2, l, k = I.get_k(), k1 = k + 1, k2 = int(k1 / 2), ntau = A.ntau();
    cplx a, c;
    double hv = beta / ntau;
    // matsubara: use midpoint to keep possible symmetry, use fwd/bwd near endpoints
    // PROBLEM WITH MIDPOINT, USE FWD AND BWD ONLY
    // for(tstp2=0;tstp2<k2;tstp2++)//imaginary time
    for (tstp2 = 0; tstp2 < k1; tstp2++) // imaginary time
    {
        // fwd diff near tau=0
        c = 0;
        for (l = 0; l <= k1; l++) {
            A.get_mat(tstp2 + l, a);
            c -= a * I.bd_weights(l); // fwd has same weights as bwd with minus sign
        }
        c /= -hv;
        dA.set_mat(tstp2, c);

        // bwd diff near tau=beta
        c = 0;
        for (l = 0; l <= k1; l++) {
            A.get_mat(ntau - tstp2 - l, a);
            c += a * I.bd_weights(l);
        }
        c /= -hv;
        dA.set_mat(ntau - tstp2, c);
    }
    // for the rest use midpoint differentiation
    // PROBLEM WITH MIDPOINT, USE FWD AND BWD ONLY
    // for(tstp2=k2;tstp2<=ntau-k2;tstp2++)//imaginary time
    for (tstp2 = k1; tstp2 <= ntau - k1; tstp2++) // imaginary time
    {
        c = 0;
        /*for(l=-k2;l<=k2;l++)
         {
         A.get_mat(tstp2-l,a);
         //c += a*I.mp_weights(l);
         }*/
        for (l = 0; l <= k1; l++) {
            A.get_mat(tstp2 - l, a);
            c += a * I.bd_weights(l);
        }
        c /= -hv;
        dA.set_mat(tstp2, c);
    }

    /*

     // check symmetry Matsubara function on (0,beta)
     A.get_mat(0,a);
     A.get_mat(ntau,c);
     if( abs(a-c) <= 1e-15 )//symmetric
     {
     dA.set_mat(int(ntau/2),0);//dA(ntau/2)=0 due to symmetry
     for(tstp2=int(ntau/2+1);tstp2<=ntau;tstp2++)//imaginary time
     {
     c=0;
     for(l=0;l<=k1;l++)
     //for(l=-k2,k2,l++)
     {
     A.get_mat(tstp2-l,a);
     c += a*I.bd_weights(l);
     //A.get_mat(tstp2-l,a);
     //c += a*I.mp_weights(l)
     }
     c/=-hv;
     dA.set_mat(tstp2,c);
     dA.set_mat(ntau-tstp2,-c);
     }
     }
     else if( abs(a+c) <= 1e-15 )//antisymmetric
     {
     for(tstp2=int(ntau/2);tstp2<=ntau;tstp2++)//imaginary time
     {
     c=0;
     for(l=0;l<=k1;l++)
     {
     A.get_mat(tstp2-l,a);
     c += a*I.bd_weights(l);
     }
     c/=-hv;
     dA.set_mat(tstp2,c);
     dA.set_mat(ntau-tstp2,c);//note that tstp2=ntau/2 is set twice, they are indeed
     identical
     }

     }
     else{ std::cerr << "matsubara function neither symmetic nor antisymmetric" << std::endl;
     abort();}*/
}

// lesser/greater: deriv1= +i d/dt
// element access only for retarded and lesser
// A not hermitian: keep track of boundaries explicitly and use Acc when needed
// not very efficient (many if/else...) improve when desired
/// @private
template <typename T>
void deriv1_element(int tstp1, int tstp2, herm_matrix<T> &dA, herm_matrix<T> &A,
                    herm_matrix<T> &Acc, integration::Integrator<T> &I, T h) {
    typedef std::complex<T> cplx;
    int l, k = I.get_k(), k1 = k + 1;
    cplx a, c, plusi;
    plusi = cplx(0, 1);
    if ((tstp1 == -1) || (tstp2 == -1)) {
        std::cerr << "deriv1_element only on real times" << std::endl;
        abort();
    } else {
        if (tstp2 <= tstp1) {
            // retarded
            c = 0;
            if (tstp1 < k1) {
                for (l = 0; l <= k; l++) // NB one order less then differentiation!!
                {
                    // Use get_ret but keep track of i<j by hand (non-hermitian)
                    // Later: replace get_ret with pointer
                    // Even later: use function increment
                    if (l < tstp2) {
                        Acc.get_ret(tstp2, l, a);
                        a.real(-a.real()); // make minus conjugate by hand
                    } else {
                        A.get_ret(l, tstp2, a);
                    }
                    c += a * I.poly_differentiation(tstp1, l);
                }
            } else {
                for (l = 0; l <= k1; l++) {
                    if ((tstp1 - l) < tstp2) {
                        Acc.get_ret(tstp2, tstp1 - l, a);
                        a.real(-a.real()); // make minus conjugate by hand
                    } else {
                        A.get_ret(tstp1 - l, tstp2, a);
                    }
                    c += a * I.bd_weights(l);
                }
            }
            c = c * plusi / h;
            dA.set_ret(tstp1, tstp2, c);
        }
        if (tstp2 >= tstp1) {
            // lesser
            c = 0;
            if (tstp1 < k1) {
                for (l = 0; l <= k; l++) // NB one order less then differentiation!!
                {
                    if (l > tstp2) {
                        Acc.get_les(tstp2, l, a);
                        a.real(-a.real()); // make minus conjugate by hand
                    } else {
                        A.get_les(l, tstp2, a);
                    }
                    c += a * I.poly_differentiation(tstp1, l);
                }
            } else {
                for (l = 0; l <= k1; l++) {
                    // get_les automatically returns -conj(Gles(tstp2,tstp-l)) for
                    // tsp1-l>tstp2
                    if (l > tstp2) {
                        Acc.get_les(tstp2, tstp1 - l, a);
                        a.real(-a.real()); // make minus conjugate by hand
                    } else {
                        A.get_les(tstp1 - l, tstp2, a);
                    }
                    c += a * I.bd_weights(l);
                }
            }
            c = c * plusi / h;
            dA.set_les(tstp1, tstp2, c);
        }
    }
}

// lesser/greater: deriv2= -i d/dt
// element access only for retarded and lesser
// not very efficient (many if/else...) improve when needed
/// @private
template <typename T>
void deriv2_element(int tstp1, int tstp2, herm_matrix<T> &dA, herm_matrix<T> &A,
                    herm_matrix<T> &Acc, integration::Integrator<T> &I, T h) {
    typedef std::complex<T> cplx;
    int l, k = I.get_k(), k1 = k + 1;
    cplx a, c, minusi;
    minusi = cplx(0, -1);
    if ((tstp1 == -1) || (tstp2 == -1)) {
        std::cerr << "deriv2_element only on real times" << std::endl;
        abort();
    } else {
        if (tstp2 <= tstp1) {
            // retarded
            c = 0;
            if (tstp2 < k1) {
                for (l = 0; l <= k; l++) // NB one order less then differentiation!!
                {
                    if (tstp1 < l) {
                        Acc.get_ret(l, tstp1, a);
                        a.real(-a.real());
                    } else {
                        A.get_ret(tstp1, l, a);
                    }
                    c += a * I.poly_differentiation(tstp2, l);
                }
            } else {
                for (l = 0; l <= k1; l++) {
                    if (tstp1 < (tstp2 - l)) {
                        Acc.get_ret(tstp2 - l, tstp1, a);
                        a.real(-a.real());
                    } else {
                        A.get_ret(tstp1, tstp2 - l, a);
                    }
                    c += a * I.bd_weights(l);
                }
            }
            c = c * minusi / h;
            dA.set_ret(tstp1, tstp2, c);
        }
        if (tstp2 >= tstp1) {
            // lesser
            c = 0;
            if (tstp2 < k1) {
                for (l = 0; l <= k; l++) // NB one order less then differentiation!!
                {
                    if (tstp1 > l) {
                        Acc.get_les(l, tstp1, a);
                        a.real(-a.real());
                    } else {
                        A.get_les(tstp1, l, a);
                    }
                    c += a * I.poly_differentiation(tstp2, l);
                }
            } else {
                for (l = 0; l <= k1; l++) {
                    // get_les automatically returns -conj(Gles(tstp2,tstp-l)) for
                    // tsp1-l>tstp2
                    if (tstp1 > tstp2 - l) {
                        Acc.get_les(tstp2 - l, tstp1, a);
                        a.real(-a.real());
                    } else {
                        A.get_les(tstp1, tstp2 - l, a);
                    }
                    c += a * I.bd_weights(l);
                }
            }
            c = c * minusi / h;
            dA.set_les(tstp1, tstp2, c);
        }
    }
}

// mixed tv: deriv1 = +i d/dt
// calculated directly whole timestep (sum imaginary time)
/// @private
template <typename T>
void deriv1_tv(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A,
               integration::Integrator<T> &I, T h) {
    typedef std::complex<T> cplx;
    int tstp2, l, k = I.get_k(), k1 = k + 1, ntau = A.ntau();
    cplx a, c, plusi;
    plusi = cplx(0, 1);
    if (tstp == -1) {
        std::cerr << "deriv1_tv requires real time (left)" << std::endl;
        abort();
    } else {
        // tv (mixed: left real, right imag)
        for (tstp2 = 0; tstp2 <= ntau; tstp2++) // imaginary time
        {
            c = 0;
            if (tstp < k1) {
                for (l = 0; l <= k; l++) // NB one order less then differentiation!!
                {
                    // tv well-behaved for all t and tau
                    A.get_tv(l, tstp2, a);
                    c += a * I.poly_differentiation(tstp, l);
                }
            } else {
                for (l = 0; l <= k1; l++) {
                    // tv well-behaved for all t and tau
                    A.get_tv(tstp - l, tstp2, a);
                    c += a * I.bd_weights(l);
                }
            }
            c = c * plusi / h;
            dA.set_tv(tstp, tstp2, c);
        }
    }
}

// mixed tv: deriv2 = +d/dtau'
// calculates directly whole timestep (sum imaginary time)
/// @private
template <typename T>
void deriv2_tv(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A,
               integration::Integrator<T> &I, T beta) {
    typedef std::complex<T> cplx;
    int tstp2, l, k = I.get_k(), k1 = k + 1, ntau = A.ntau();
    cplx a, c;
    double hv = beta / ntau;
    if (tstp == -1) {
        std::cerr << "deriv2_tv requires real time (left)" << std::endl;
        abort();
    } else {
        // tv (mixed: left real, right imag)
        for (tstp2 = 0; tstp2 <= ntau; tstp2++) // imaginary time
        {
            c = 0;
            if (tstp2 < k1) {
                for (l = 0; l <= k; l++) // NB one order less then differentiation!!
                {
                    // tv well-behaved for all t and tau
                    A.get_tv(tstp, l, a);
                    c += a * I.poly_differentiation(tstp2, l);
                }
            } else {
                for (l = 0; l <= k1; l++) {
                    // tv well-behaved for all t and tau
                    A.get_tv(tstp, tstp2 - l, a);
                    c += a * I.bd_weights(l);
                }
            }
            c /= hv;
            dA.set_tv(tstp, tstp2, c);
        }
    }
}


// ============================== new interfaces ==============================

template <typename T>
void deriv1_matsubara(herm_matrix<T> &dA, herm_matrix<T> &A, T beta, int SolveOrder){
    deriv1_matsubara(dA, A, integration::I<T>(SolveOrder), beta);
}
template <typename T>
void deriv1_element(int tstp1, int tstp2, herm_matrix<T> &dA, herm_matrix<T> &A,
                    herm_matrix<T> &Acc, T h, int SolveOrder){
    deriv1_element(tstp1, tstp2, dA, A, Acc, integration::I<T>(SolveOrder), h);
}
template <typename T>
void deriv2_element(int tstp1, int tstp2, herm_matrix<T> &dA, herm_matrix<T> &A,
                    herm_matrix<T> &Acc, T h, int SolveOrder){
    deriv2_element(tstp1, tstp2, dA, A, Acc, integration::I<T>(SolveOrder), h);
}
template <typename T>
void deriv1_tv(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A,
               T h, int SolveOrder){
    deriv1_tv(tstp, dA, A, integration::I<T>(SolveOrder), h);
}
template <typename T>
void deriv2_tv(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A,
               T beta, int SolveOrder){
    deriv2_tv(tstp, dA, A, integration::I<T>(SolveOrder), beta);
}
template <typename T>
void deriv1_timestep(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A, herm_matrix<T> &Acc,
                     T beta, T h, int SolveOrder){
    deriv1_timestep(tstp, dA, A, Acc, integration::I<T>(SolveOrder), beta, h);
}
template <typename T>
void deriv2_timestep(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A, herm_matrix<T> &Acc,
                     T beta, T h, int SolveOrder){
    deriv2_timestep(tstp, dA, A, Acc, integration::I<T>(SolveOrder), beta, h);
}

} // namespace cntr

#endif  // CNTR_DIFFERENTIATION_IMPL_H
