#ifndef CNTR_DYSON_OMP_DECL_H
#define CNTR_DYSON_OMP_DECL_H

#include "cntr_global_settings.hpp"
#include "integration.hpp"

namespace cntr {

template <typename T> class function;
template <typename T> class herm_matrix;
/// @private
template <typename T> class herm_pseudo;

#if CNTR_USE_OMP == 1

template <typename T, class GG, int SIZE1>
void dyson_timestep_ret_omp(int omp_num_threads, int n, GG &G, T mu, std::complex<T> *H,
                            GG &Sigma, integration::Integrator<T> &I, T h);
template <typename T, class GG, int SIZE1>
void dyson_timestep_tv_omp(int omp_num_threads, int n, GG &G, T mu, std::complex<T> *Hn,
                           GG &Sigma, integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T, class GG, int SIZE1>
void pseudodyson_timestep_tv_omp(int omp_num_threads, int n, GG &G, T mu,
                                 std::complex<T> *Hn, GG &Sigma,
                                 integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void dyson_timestep_les_omp(int omp_num_threads, int n, GG &G, T mu, std::complex<T> *H,
                            GG &Sigma, integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T>
void pseudodyson_timestep_omp(int omp_num_threads, int n, herm_pseudo<T> &G, T lam0,
                              function<T> &H, herm_pseudo<T> &Sigma,
                              integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T>
void pseudodyson_timestep_omp(int omp_num_threads, int n, herm_pseudo<T> &G, T lam0,
                              std::complex<T> *Ht, herm_pseudo<T> &Sigma,
                              integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T>
void pseudodyson_timestep_omp(int omp_num_threads, int n, herm_pseudo<T> &G, T lam0,
                              std::vector<std::complex<T>> &Ht, herm_pseudo<T> &Sigma,
                              integration::Integrator<T> &I, T beta, T h);
template <typename T>
void dyson_timestep_omp(int omp_num_threads, int n, herm_matrix<T> &G, T lam0,
                        function<T> &H, herm_matrix<T> &Sigma, integration::Integrator<T> &I,
                        T beta, T h);

#endif // CNTR_USE_OMP

}  // namespace cntr

#endif  // CNTR_DYSON_OMP_DECL_H
