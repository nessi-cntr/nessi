#ifndef CNTR_VIE2_DECL_H
#define CNTR_VIE2_DECL_H

#include "cntr_global_settings.hpp"
#include "integration.hpp"

namespace cntr {

  template <typename T> class herm_matrix;
  template <typename T> class function;

/* #######################################################################################
#
#  [1+F]G=Q, G,Q hermitian
#
###########################################################################################*/
  template <typename T>
  void vie2_mat(herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc, herm_matrix<T> &Q,
		T beta, integration::Integrator<T> &I, const int method=CNTR_MAT_FIXPOINT);

  template <typename T>
  void vie2_mat_fourier(herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc, herm_matrix<T> &Q,
		T beta, int order = 3);
  template <typename T>
  void vie2_mat_fixpoint(herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc,
		      herm_matrix<T> &Q, T beta, integration::Integrator<T> &I,
		      int nfixpoint = 6);
  template <typename T>
  void vie2_mat_steep(herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc,
		   herm_matrix<T> &Q, T beta, integration::Integrator<T> &I,
		   int maxiter = 8, T tol=1.0e-16);
  template <typename T>
  void vie2_start(herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc, herm_matrix<T> &Q,
		  integration::Integrator<T> &I, T beta, T h);
  template <typename T>
  void vie2_timestep(int n, herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc,
		     herm_matrix<T> &Q, integration::Integrator<T> &I, T beta, T h,
		     const int matsubara_method=CNTR_MAT_FIXPOINT);
  template <typename T>
  void vie2(herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc, herm_matrix<T> &Q,
	    integration::Integrator<T> &I, T beta, T h, const int matsubara_method=CNTR_MAT_FIXPOINT);
  template <typename T>
  void vie2_timestep_sin(int n,herm_matrix<T> &G,function<T> &Gsin,herm_matrix<T> &F,herm_matrix<T> &Fcc, function<T> &Fsin ,
      herm_matrix<T> &Q,function<T> &Qsin,T beta,T h,int kt);

/* /////////////////////////////////////////////////////////////////////////////////////////
// OpenMP
///////////////////////////////////////////////////////////////////////////////////////// */

#if CNTR_USE_OMP == 1
  template <typename T>
  void vie2_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &G, herm_matrix<T> &F,
			 herm_matrix<T> &Fcc, herm_matrix<T> &Q, integration::Integrator<T> &I,
			 T beta, T h, const int matsubara_method=CNTR_MAT_FIXPOINT);
  template <typename T>
  void vie2_omp(int omp_num_threads, herm_matrix<T> &G, herm_matrix<T> &F, herm_matrix<T> &Fcc,
		herm_matrix<T> &Q, integration::Integrator<T> &I, T beta, T h,
		const int matsubara_method=CNTR_MAT_FIXPOINT);
#endif // CNTR_USE_OMP

} // namespace cntr

#endif  // CNTR_VIE2_DECL_H
