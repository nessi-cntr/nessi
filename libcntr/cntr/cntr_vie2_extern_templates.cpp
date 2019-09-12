#include "cntr_vie2_extern_templates.hpp"
#include "cntr_vie2_impl.hpp"

namespace cntr {

  template void vie2_mat<double>(herm_matrix<double> &G, herm_matrix<double> &F, herm_matrix<double> &Fcc,
			herm_matrix<double> &Q, double beta, integration::Integrator<double> &I,
			const int method);

  template void vie2_mat<double>(double beta, herm_matrix<double> &G, herm_matrix<double> &F, herm_matrix<double> &Fcc,
			herm_matrix<double> &Q, const int kt, const int method);
  
  template void vie2_mat_fourier<double>(herm_matrix<double> &G,herm_matrix<double> &F,
					 herm_matrix<double> &Fcc,
					 herm_matrix<double> &Q,double beta,int order);
  template void vie2_mat_fixpoint<double>(herm_matrix<double> &G,herm_matrix<double> &F,herm_matrix<double> &Fcc,
				       herm_matrix<double> &Q,double beta,integration::Integrator<double> &I,
				       int nfixpoint);
  template void vie2_mat_steep<double>(herm_matrix<double> &G, herm_matrix<double> &F, herm_matrix<double> &Fcc,
			    herm_matrix<double> &Q, double beta, integration::Integrator<double> &I,
			    int maxiter, double tol);
  
  template void vie2_start<double>(herm_matrix<double> &G,herm_matrix<double> &F,herm_matrix<double>
				   &Fcc,herm_matrix<double> &Q, integration::Integrator<double> &I,
				   double beta,double h);
  template void vie2_start<double>(double beta,double h, herm_matrix<double> &G,herm_matrix<double> &F,herm_matrix<double>
				   &Fcc,herm_matrix<double> &Q, const int kt);

  template void vie2_timestep<double>(int n,herm_matrix<double> &G,herm_matrix<double> &F,herm_matrix<double> &Fcc,
				      herm_matrix<double> &Q, integration::Integrator<double> &I,
				      double beta,double h, const int matsubara_method);
  template void vie2_timestep<double>(int n, double beta, double h, herm_matrix<double> &G,herm_matrix<double> &F,herm_matrix<double> &Fcc,
				      herm_matrix<double> &Q, const int kt,
				      const int matsubara_method);

  template void vie2<double>(herm_matrix<double> &G,herm_matrix<double> &F,herm_matrix<double> &Fcc,
			     herm_matrix<double> &Q, integration::Integrator<double> &I, double beta,double h,
			     const int matsubara_method);
  template void vie2<double>(double beta,double h, herm_matrix<double> &G,herm_matrix<double> &F,herm_matrix<double> &Fcc,
			     herm_matrix<double> &Q, const int kt, const int matsubara_method);


  template void vie2_timestep_sin(int n,herm_matrix<double> &G,function<double> &Gsin,herm_matrix<double> &F,herm_matrix<double> &Fcc, function<double> &Fsin ,
      			 herm_matrix<double> &Q,function<double> &Qsin,double beta,double h,int kt);
#if CNTR_USE_OMP==1
  template void vie2_timestep_omp<double>(int omp_num_threads,int n,herm_matrix<double> &G,
					  herm_matrix<double> &F,herm_matrix<double> &Fcc,herm_matrix<double> &Q,
					  integration::Integrator<double> &I, double beta,double h,
					  const int matsubara_method);
  template void vie2_timestep_omp<double>(int omp_num_threads,int n,double beta,double h,herm_matrix<double> &G,
					  herm_matrix<double> &F,herm_matrix<double> &Fcc,herm_matrix<double> &Q,
					  const int kt, 
					  const int matsubara_method);
  template void vie2_omp<double>(int omp_num_threads, double beta, double h, herm_matrix<double> &G, herm_matrix<double> &F,
				 herm_matrix<double> &Fcc, herm_matrix<double> &Q,
				 const int kt,
				 const int matsubara_method);
#endif // CNTR_USE_OMP

}  // namespace cntr
