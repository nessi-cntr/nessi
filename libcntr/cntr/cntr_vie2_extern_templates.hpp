#ifndef CNTR_VIE2_EXTERN_TEMPLATES_H
#define CNTR_VIE2_EXTERN_TEMPLATES_H

#include "cntr_vie2_decl.hpp"

namespace cntr {

  extern template
  void vie2_mat<double>(herm_matrix<double> &G, herm_matrix<double> &F, herm_matrix<double> &Fcc,
	   herm_matrix<double> &Q, double beta, integration::Integrator<double> &I, const int method);
  extern template
  void vie2_mat_fourier<double>(herm_matrix<double> &G,herm_matrix<double> &F,herm_matrix<double> &Fcc,
			herm_matrix<double> &Q,double beta,int order);
  extern template
  void vie2_mat_fixpoint<double>(herm_matrix<double> &G,herm_matrix<double> &F,herm_matrix<double> &Fcc,
			      herm_matrix<double> &Q,double beta,integration::Integrator<double> &I, int nfixpoint);
  extern template
  void vie2_mat_steep<double>(herm_matrix<double> &G, herm_matrix<double> &F, herm_matrix<double> &Fcc,
		   herm_matrix<double> &Q, double beta, integration::Integrator<double> &I,
		   int maxiter, double tol);

  extern template
  void vie2_start<double>(herm_matrix<double> &G,herm_matrix<double> &F,herm_matrix<double> &Fcc,
			  herm_matrix<double> &Q, integration::Integrator<double> &I, double beta,double h);
  extern template
  void vie2_timestep<double>(int n,herm_matrix<double> &G,herm_matrix<double> &F,herm_matrix<double> &Fcc,
			     herm_matrix<double> &Q, integration::Integrator<double> &I, double beta,double h,
			     const int matsubara_method);
  extern template
  void vie2<double>(herm_matrix<double> &G,herm_matrix<double> &F,herm_matrix<double> &Fcc,
		    herm_matrix<double> &Q, integration::Integrator<double> &I, double beta,double h,
		    const int matsubara_method);
  extern template
  void vie2_timestep_sin(int n,herm_matrix<double> &G,function<double> &Gsin,herm_matrix<double> &F,herm_matrix<double> &Fcc,
        function<double> &Fsin,herm_matrix<double> &Q,function<double> &Qsin,double beta,double h,int kt);
#if CNTR_USE_OMP == 1
  extern template
  void vie2_timestep_omp<double>(int omp_num_threads, int n, herm_matrix<double> &G,
				 herm_matrix<double> &F, herm_matrix<double> &Fcc, herm_matrix<double> &Q,
				 integration::Integrator<double> &I, double beta, double h,
				 const int matsubara_method);
  extern template
  void vie2_omp<double>(int omp_num_threads, herm_matrix<double> &G, herm_matrix<double> &F,
			herm_matrix<double> &Fcc, herm_matrix<double> &Q, integration::Integrator<double> &I, double beta, double h,
			const int matsubara_method);
#endif // CNTR_USE_OMP

}  // namespace cntr

#endif  // CNTR_VIE2_EXTERN_TEMPLATES_H
