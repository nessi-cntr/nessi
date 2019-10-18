#ifndef CNTR_MATSUBARA_DECL_H
#define CNTR_MATSUBARA_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

// ----------------------------------------------------------------------

/* /////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for imaginary time and matsubara frequencies
///////////////////////////////////////////////////////////////////////////////////////// */

/// @private
template <typename F, typename I> F get_tau(const I tau_idx, const F beta, const I ntau);
/// @private
template <typename F, typename I> F get_omega(const I m, const F beta, const I sig);

// ----------------------------------------------------------------------
/// @private
template <typename T, int SIZE1>
void set_first_order_tail(std::complex<T> *xmat, std::complex<T> *coeff, T beta, int sg,
                          int ntau, int sig, int size1);

// ----------------------------------------------------------------------

template <typename T, class GG, int SIZE1>
void matsubara_dft(std::complex<T> *mdft, GG &G, int sig);
template <typename T, class GG, int SIZE1>
void matsubara_ft(std::complex<T> *result, int m, GG &G, std::complex<T> *mdft, int sig,
                  T beta, int order);

// ----------------------------------------------------------------------

}  // namespace cntr

#endif  // CNTR_MATSUBARA_DECL_H
