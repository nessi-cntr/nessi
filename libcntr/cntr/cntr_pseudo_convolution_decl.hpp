#ifndef CNTR_PSEUDO_CONVOLUTION_DECL_H
#define CNTR_PSEUDO_CONVOLUTION_DECL_H

#include "cntr_global_settings.hpp"
#include "integration.hpp"

namespace cntr {

template <typename T> class function;
template <typename T> class herm_matrix;
template <typename T> class herm_pseudo;

/* #######################################################################################
#  C = A*f*B or C=A*B
###########################################################################################*/
/// @private
template <typename T>
void pseudo_convolution_timestep(int tstp, herm_pseudo<T> &C, herm_pseudo<T> &A,
                                 herm_pseudo<T> &Acc, function<T> &ft, herm_pseudo<T> &B,
                                 herm_pseudo<T> &Bcc, integration::Integrator<T> &I, T beta,
                                 T h);
/// @private
template <typename T>
void pseudo_convolution_timestep(int n, herm_pseudo<T> &C, herm_pseudo<T> &A,
                                 function<T> &ft, herm_pseudo<T> &B,
                                 integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T>
void pseudo_convolution_matsubara(herm_pseudo<T> &C, herm_pseudo<T> &A, herm_pseudo<T> &B,
                                  integration::Integrator<T> &I, T beta);
/// @private
template <typename T>
void pseudo_convolution_matsubara(herm_pseudo<T> &C, herm_pseudo<T> &A, function<T> &ft,
                                  herm_pseudo<T> &B, integration::Integrator<T> &I, T beta);
/// @private
template <typename T>
void pseudo_convolution(herm_pseudo<T> &C, herm_pseudo<T> &A, herm_pseudo<T> &Acc,
                        function<T> &ft, herm_pseudo<T> &B, herm_pseudo<T> &Bcc,
                        integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T>
void pseudo_convolution_timestep(int tstp, herm_pseudo<T> &C, herm_pseudo<T> &A,
                                 herm_pseudo<T> &Acc, herm_pseudo<T> &B, herm_pseudo<T> &Bcc,
                                 integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T>
void pseudo_convolution_timestep(int n, herm_pseudo<T> &C, herm_pseudo<T> &A,
                                 herm_pseudo<T> &B, integration::Integrator<T> &I, T beta,
                                 T h);
/// @private
template <typename T>
void pseudo_convolution(herm_pseudo<T> &C, herm_pseudo<T> &A, herm_pseudo<T> &Acc,
                        herm_pseudo<T> &B, herm_pseudo<T> &Bcc,
                        integration::Integrator<T> &I, T beta, T h);

/* #######################################################################################
#  Components
###########################################################################################*/

#define CPLX std::complex<T>
/// @private
template <typename T, class GG, int SIZE1>
void incr_pseudo_convolution_tv(int tstp, std::vector<bool> &mask, CPLX alpha, GG &C, GG &A,
                                GG &Acc, CPLX *f0, CPLX *ft, GG &B, GG &Bcc,
                                integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T, class GG, int SIZE1>
void incr_pseudo_convolution_mat(std::vector<bool> &mask, CPLX alpha, GG &C, GG &A, CPLX *f0,
                                 GG &B, integration::Integrator<T> &I, T beta);

#undef CPLX

}  // namespace cntr

#endif  // CNTR_PSEUDO_CONVOLUTION_DECL_H
