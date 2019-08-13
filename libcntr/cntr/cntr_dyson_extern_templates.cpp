#include "cntr_dyson_extern_templates.hpp"
#include "cntr_dyson_impl.hpp"

namespace cntr {
  template
  void dyson_mat_fourier<double>(herm_matrix<double> &G,herm_matrix<double> &Sigma,double mu,function<double> &H,double beta,int order=3);
  template
  void dyson_mat_fourier<double>(herm_matrix<double> &G,herm_matrix<double> &Sigma,double mu,function<double> &H,
  				 function<double> &SigmaMF, double beta,int order=3);

  template
  void dyson_mat_fixpoint<double>(herm_matrix<double> &G,herm_matrix<double> &Sigma,double mu,function<double> &H,
  				  integration::Integrator<double> &I, double beta,int fixpiter=6);
  template
  void dyson_mat_fixpoint<double>(herm_matrix<double> &G,herm_matrix<double> &Sigma,double mu,function<double> &H,
  				  function<double> &SigmaMF,
  				  integration::Integrator<double> &I, double beta,int fixpiter=6);

  template
  void dyson_mat_steep<double>(herm_matrix<double> &G, herm_matrix<double> &Sigma, double mu, function<double> &H,
  			       integration::Integrator<double> &I, double beta, int maxiter = 10, double tol=1.0e-16);
  template
  void dyson_mat_steep<double>(herm_matrix<double> &G, herm_matrix<double> &Sigma, double mu, function<double> &H,
  			       function<double> &SigmaMF, integration::Integrator<double> &I, double beta,
  			       int maxiter = 10, double tol=1.0e-16);
  template
  void dyson_mat<double>(herm_matrix<double> &G, herm_matrix<double> &Sigma, double mu, function<double> &H,
  			 integration::Integrator<double> &I, double beta, const int method=CNTR_MAT_FIXPOINT,
         const bool force_hermitian=true);
  template
  void dyson_mat<double>(herm_matrix<double> &G, herm_matrix<double> &Sigma, double mu, function<double> &H,
  			 function<double> &SigmaMF, integration::Integrator<double> &I, double beta,
			      const int method=CNTR_MAT_FIXPOINT, const bool force_hermitian=true);


  template
  void dyson_start<double>(herm_matrix<double> &G,double mu,function<double> &H, herm_matrix<double> &Sigma,
  			   integration::Integrator<double> &I, double beta,double h);
  template
  void dyson_timestep<double>(int n,herm_matrix<double> &G,double mu,function<double> &H, herm_matrix<double> &Sigma,
  			      integration::Integrator<double> &I, double beta,double h);
  template
  void dyson<double>(herm_matrix<double> &G,double mu,function<double> &H, herm_matrix<double> &Sigma,
  		     integration::Integrator<double> &I, double beta,double h, const int matsubara_method=CNTR_MAT_FIXPOINT,
          const bool force_hermitian=true);

}  // namespace cntr
