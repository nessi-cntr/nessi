#ifndef CNTR_PSEUDODYSON_DECL_H
#define CNTR_PSEUDODYSON_DECL_H

#include "cntr_global_settings.hpp"
#include "integration.hpp"

namespace cntr {

template <typename T> class function;
template <typename T> class herm_pseudo;

/*#########################################################################################
#
#  PSEUDOPARTICLE DYSON EQUATION: [ id/dt + lam0 - H -Sigma ] G = 1
#
#  ANALOGOUS TO USUAL DYSON, BUT WITH DIFFERENT CONVOLUTION
#
##########################################################################################*/
// for matrix geenfunctions: H passed as Matrix
/// @private
template <typename T, class Matrix>
void pseudodyson_mat(herm_pseudo<T> &G, T lam0, Matrix &H0, herm_pseudo<T> &Sigma, T beta1,
                     integration::Integrator<T> &I, T beta2);
/// @private
template <typename T, class Matrix>
void pseudodyson_start(herm_pseudo<T> &G, T lam0, Matrix &H0, std::vector<Matrix> &H,
                       herm_pseudo<T> &Sigma, integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T, class Matrix>
void pseudodyson_timestep(int n, herm_pseudo<T> &G, T lam0, Matrix &H0,
                          std::vector<Matrix> &H, herm_pseudo<T> &Sigma,
                          integration::Integrator<T> &I, T beta, T h);
//
/// @private
template <typename T>
void pseudodyson_mat(herm_pseudo<T> &G, T lam0, function<T> &H0, herm_pseudo<T> &Sigma,
                     T beta1, integration::Integrator<T> &I, T beta2);
/// @private
template <typename T>
void pseudodyson_start(herm_pseudo<T> &G, T lam0, function<T> &H, herm_pseudo<T> &Sigma,
                       integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T>
void pseudodyson_timestep(int n, herm_pseudo<T> &G, T lam0, function<T> &H,
                          herm_pseudo<T> &Sigma, integration::Integrator<T> &I, T beta, T h);
// the following work for size1=1 only
/// @private
template <typename T>
void pseudodyson_mat(herm_pseudo<T> &G, T lam0, std::complex<T> &H0, herm_pseudo<T> &Sigma,
                     T beta1, integration::Integrator<T> &I, T beta2);
/// @private
template <typename T>
void pseudodyson_start(herm_pseudo<T> &G, T lam0, std::complex<T> H0,
                       std::vector<std::complex<T>> &H, herm_pseudo<T> &Sigma,
                       integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T>
void pseudodyson_timestep(int tstp, herm_pseudo<T> &G, T lam0, std::complex<T> H0,
                          std::vector<std::complex<T>> &H, herm_pseudo<T> &Sigma,
                          integration::Integrator<T> &I, T beta, T h);

} // namespace cntr

#endif  // CNTR_PSEUDODYSON_DECL_H
