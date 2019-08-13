#ifndef CNTR_CONVOLUTION_DECL_H
#define CNTR_CONVOLUTION_DECL_H

#include "cntr_global_settings.hpp"
#include "integration.hpp"

namespace cntr {

template <typename T> class function;
template <typename T> class herm_matrix;

/*###########################################################################################
#
#   CONVOLUTION  C=A*B   or  C=A*f*B
#
#   ...  currently only implemented for herm_matrix
#   to compute the convolution of non-hermitian greenfunctions, use the second one ...
#   implemented in green_cntr_convolution.hpp
#
###########################################################################################*/
template <typename T, class GG>
void convolution_matsubara(GG &C, GG &A, GG &B, integration::Integrator<T> &I, T beta);
#if CNTR_USE_OMP == 1
template <typename T, class GG>
void convolution_matsubara_omp(int nomp, GG &C, GG &A, GG &B, integration::Integrator<T> &I, T beta);
#endif
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                          herm_matrix<T> &B, herm_matrix<T> &Bcc,
                          integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &B,
                          integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                 herm_matrix<T> &B, herm_matrix<T> &Bcc, integration::Integrator<T> &I,
                 T beta, T h);
// with f
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                          function<T> &ft, herm_matrix<T> &B, herm_matrix<T> &Bcc,
                          integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, function<T> &ft,
                          herm_matrix<T> &B, integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc, function<T> &ft,
                 herm_matrix<T> &B, herm_matrix<T> &Bcc, integration::Integrator<T> &I,
                 T beta, T h);
//
// only density matrix is computed ... defined in cntr_convolution_with_f.hpp
// thus this works with for GG=herm_pseudo and GG=herm_matrix
template <typename T, class GG>
void convolution_density_matrix(int tstp, std::complex<T> *rho, GG &A, GG &Acc,
                                function<T> &ft, GG &B, GG &Bcc,
                                integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG>
void convolution_density_matrix(int tstp, std::complex<T> *rho, GG &A, function<T> &ft,
                                GG &B, integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG>
void convolution_density_matrix(int tstp, std::complex<T> *rho, GG &A, GG &Acc, GG &B,
                                GG &Bcc, integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG>
void convolution_density_matrix(int tstp, std::complex<T> *rho, GG &A, GG &B,
                                integration::Integrator<T> &I, T beta, T h);

template <typename T, class GG>
void convolution_density_matrix(int tstp, cdmatrix &Cles, GG &A, GG &B,
                  integration::Integrator<T> &I, T beta, T h);

template <typename T, class GG>
void convolution_les_timediag(int tstp, cdmatrix &Cles, GG &A, GG &B,
                                integration::Integrator<T> &I, T beta, T h);

template <typename T, class GG>
void convolution_matsubara(GG &C, GG &A, std::complex<T> *f0, GG &B,
                           integration::Integrator<T> &I, T beta);
template <typename T, class GG, int SIZE1>
void convolution_timestep_ret(int n, GG &C, GG &A, GG &Acc, std::complex<T> *ft, GG &B,
                              GG &Bcc, integration::Integrator<T> &I, T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_tv(int n, std::complex<T> *ctv, GG &C, GG &A, GG &Acc, GG &B,
                             GG &Bcc, integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_tv(int n, GG &C, GG &A, GG &Acc, GG &B, GG &Bcc,
                             integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_tv(int n, std::complex<T> *ctv, GG &C, GG &A, GG &Acc,
                             std::complex<T> *f0, std::complex<T> *ft, GG &B, GG &Bcc,
                             integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_tv(int n, GG &C, GG &A, GG &Acc, std::complex<T> *f0,
                             std::complex<T> *ft, GG &B, GG &Bcc,
                             integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_tvvt(int n, int j1, int j2, std::complex<T> *cles, GG &C,
                                   GG &A, GG &Acc, GG &B, GG &Bcc,
                                   integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_tvvt(int n, std::complex<T> *cles, GG &C, GG &A, GG &Acc,
                                   GG &B, GG &Bcc, integration::Integrator<T> &I, T beta,
                                   T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_tvvt(int n, std::complex<T> *cles, GG &C, GG &A, GG &Acc,
                                   std::complex<T> *f0, GG &B, GG &Bcc,
                                   integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_lesadv(int n, int j1, int j2, std::complex<T> *cles, GG &C,
                                     GG &A, GG &Acc, GG &B, GG &Bcc,
                                     integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_lesadv(int n, std::complex<T> *cles, GG &C, GG &A, GG &Acc,
                                     GG &B, GG &Bcc, integration::Integrator<T> &I, T beta,
                                     T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_lesadv(int n, std::complex<T> *cles, GG &C, GG &A, GG &Acc,
                                     std::complex<T> *ft, GG &B, GG &Bcc,
                                     integration::Integrator<T> &I, T beta, T h);
template < typename T, class GG, int SIZE1 >
void convolution_timestep_les_retles(int n,std::complex<T> *cles,GG &C,GG &A, GG &Acc,
std::complex<T> *ft,GG &B,GG &Bcc,integration::Integrator<T> &I,T beta,T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_les(int n, GG &C, GG &A, GG &Acc, std::complex<T> *f0,
                              std::complex<T> *ft, GG &B, GG &Bcc,
                              integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                          std::complex<T> *f0, std::complex<T> *ft, herm_matrix<T> &B,
                          herm_matrix<T> &Bcc, integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, std::complex<T> *f0,
                          std::complex<T> *ft, herm_matrix<T> &B,
                          integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                 std::complex<T> *f0, std::complex<T> *ft, herm_matrix<T> &B,
                 herm_matrix<T> &Bcc, integration::Integrator<T> &I, T beta, T h);

template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                          function<T> &ft, herm_matrix<T> &B, herm_matrix<T> &Bcc,
                          integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, function<T> &ft,
                          herm_matrix<T> &B, integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc, function<T> &ft,
                 herm_matrix<T> &B, herm_matrix<T> &Bcc, integration::Integrator<T> &I,
                 T beta, T h);

template <typename T, class GG, int SIZE1>
void convolution_timestep_les_jn_lesadv(int j, int n, std::complex<T> *cc, int sizec, GG &A,
                                        GG &Acc, std::complex<T> *ft, GG &B, GG &Bcc,
                                        integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_jn_tvvt(int j, int n, std::complex<T> *cc, int sizec, GG &A,
                                      GG &Acc, std::complex<T> *f0, GG &B, GG &Bcc,
                                      integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_jn_retles(int j, int n, std::complex<T> *cc, int sizec, GG &A,
                                        GG &Acc, std::complex<T> *ft, GG &B, GG &Bcc,
                                        integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_jn(int j, int n, std::complex<T> *cc, int sizec, GG &A,
                                 GG &Acc, std::complex<T> *f0, std::complex<T> *ft, GG &B,
                                 GG &Bcc, integration::Integrator<T> &I, T beta, T h);

template <typename T, class GG>
void convolution_density_matrix(int n, std::complex<T> *rho, GG &A, GG &Acc, function<T> &ft,
                                GG &B, GG &Bcc, integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG>
void convolution_density_matrix(int tstp, std::complex<T> *rho, GG &A, function<T> &ft,
                                GG &B, integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG>
void convolution_density_matrix(int n, std::complex<T> *rho, GG &A, GG &Acc, GG &B, GG &Bcc,
                                integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG>
void convolution_density_matrix(int tstp, std::complex<T> *rho, GG &A, GG &B,
                                integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG>
void convolution_density_matrix(int tstp, cdmatrix &Cles, GG &A, GG &B,
                  integration::Integrator<T> &I, T beta, T h);

/*###########################################################################################
#
#   HELPERS
#
###########################################################################################*/
/// @private
template <typename T, int SIZE1>
void matsubara_integral_1(int size1, int m, int ntau, std::complex<T> *C, std::complex<T> *A,
                          std::complex<T> *B, integration::Integrator<T> &I, int sig = -1);
/// @private
template <typename T, int SIZE1>
void matsubara_integral_1_1(int size1, int m, int ntau, std::complex<T> *C,
                            std::complex<T> *A, std::complex<T> *B,
                            integration::Integrator<T> &I);
/// @private
template <typename T, int SIZE1>
void matsubara_integral_2(int size1, int m, int ntau, std::complex<T> *C, std::complex<T> *A,
                          std::complex<T> *B, integration::Integrator<T> &I, int sig = -1);
/// @private
template <typename T, int SIZE1>
void matsubara_integral_2_2(int size1, int m, int ntau, std::complex<T> *C,
                            std::complex<T> *A, std::complex<T> *B,
                            integration::Integrator<T> &I);

/* /////////////////////////////////////////////////////////////////////////////////////////
// INCREMENETAL
///////////////////////////////////////////////////////////////////////////////////////// */
#define CPLX std::complex<T>
template <typename T, class GG, int SIZE1>
void incr_convolution_mat(std::vector<bool> &mask, CPLX alpha, GG &C, GG &A, CPLX *f0, GG &B,
                          integration::Integrator<T> &I, T beta);
template <typename T, class GG, int SIZE1>
void incr_convolution_ret(int tstp, std::vector<bool> &mask, CPLX alpha, GG &C, GG &A,
                          GG &Acc, CPLX *ft, GG &B, GG &Bcc, integration::Integrator<T> &I,
                          T h);
template <typename T, class GG, int SIZE1>
void incr_convolution_tv(int tstp, std::vector<bool> &mask, CPLX alpha, GG &C, GG &A,
                         GG &Acc, CPLX *f0, CPLX *ft, GG &B, GG &Bcc,
                         integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void incr_convolution_les(int tstp, std::vector<bool> &mask, CPLX alpha, GG &C, GG &A,
                          GG &Acc, CPLX *f0, CPLX *ft, GG &B, GG &Bcc,
                          integration::Integrator<T> &I, T beta, T h);
template <typename T, class GG, int SIZE1>
void incr_convolution(int tstp, CPLX alpha, GG &C, GG &A, GG &Acc, CPLX *f0, CPLX *ft, GG &B,
                      GG &Bcc, integration::Integrator<T> &I, T beta, T h);

//////////////////////////////////////////////////////////////////////////////////////////////////////
// NEW VERSIONS
template <typename T>
void convolution_timestep_new(int tstp, herm_matrix<T> &C, herm_matrix<T> &A,
                              herm_matrix<T> &Acc, function<T> &ft, herm_matrix<T> &B,
                              herm_matrix<T> &Bcc, integration::Integrator<T> &I, T beta,
                              T h);
template <typename T>
void convolution_timestep_new(int n, herm_matrix<T> &C, herm_matrix<T> &A, function<T> &ft,
                              herm_matrix<T> &B, integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_matsubara_new(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &B,
                               integration::Integrator<T> &I, T beta);
template <typename T>
void convolution_matsubara_new(herm_matrix<T> &C, herm_matrix<T> &A, function<T> &ft,
                               herm_matrix<T> &B, integration::Integrator<T> &I, T beta);
template <typename T>
void convolution_new(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                     function<T> &ft, herm_matrix<T> &B, herm_matrix<T> &Bcc,
                     integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_timestep_new(int tstp, herm_matrix<T> &C, herm_matrix<T> &A,
                              herm_matrix<T> &Acc, herm_matrix<T> &B, herm_matrix<T> &Bcc,
                              integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_timestep_new(int n, herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &B,
                              integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_new(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                     herm_matrix<T> &B, herm_matrix<T> &Bcc, integration::Integrator<T> &I,
                     T beta, T h);

#undef CPLX

/* #######################################################################################
#  Parallel
###########################################################################################*/
#define CPLX std::complex<T>

#if CNTR_USE_OMP == 1
template <typename T, class GG, int SIZE1>
void incr_convolution_omp(int omp_num_threads, int tstp, CPLX alpha, GG &C, GG &A, GG &Acc,
                          CPLX *f0, CPLX *ft, GG &B, GG &Bcc, integration::Integrator<T> &I,
                          T beta, T h);

template <typename T>
void convolution_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &C,
                              herm_matrix<T> &A, herm_matrix<T> &Acc, function<T> &ft,
                              herm_matrix<T> &B, herm_matrix<T> &Bcc,
                              integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &C,
                              herm_matrix<T> &A, function<T> &ft, herm_matrix<T> &B,
                              integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_matsubara_omp(int omp_num_threads, herm_matrix<T> &C, herm_matrix<T> &A,
                               herm_matrix<T> &B, integration::Integrator<T> &I, T beta);
template <typename T>
void convolution_matsubara_omp(int omp_num_threads, herm_matrix<T> &C, herm_matrix<T> &A,
                               function<T> &ft, herm_matrix<T> &B,
                               integration::Integrator<T> &I, T beta);

template <typename T>
void convolution_omp(int omp_num_threads, herm_matrix<T> &C, herm_matrix<T> &A,
                     herm_matrix<T> &Acc, function<T> &ft, herm_matrix<T> &B,
                     herm_matrix<T> &Bcc, integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &C,
                              herm_matrix<T> &A, herm_matrix<T> &Acc, herm_matrix<T> &B,
                              herm_matrix<T> &Bcc, integration::Integrator<T> &I, T beta,
                              T h);
template <typename T>
void convolution_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &C,
                              herm_matrix<T> &A, herm_matrix<T> &B,
                              integration::Integrator<T> &I, T beta, T h);
template <typename T>
void convolution_omp(int omp_num_threads, herm_matrix<T> &C, herm_matrix<T> &A,
                     herm_matrix<T> &Acc, herm_matrix<T> &B, herm_matrix<T> &Bcc,
                     integration::Integrator<T> &I, T beta, T h);
#undef CPLX
#endif // CNTR_USE_OMP

} // namespace cntr

#endif  // CNTR_CONVOLUTION_DECL_H
