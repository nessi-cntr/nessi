#include "cntr_dyson_extern_templates.hpp"
#include "cntr_dyson_impl.hpp"

namespace cntr {

// internal interfaces

  template
  void dyson_mat_fourier<double>(herm_matrix<double> &G,herm_matrix<double> &Sigma,double mu,function<double> &H,double beta,int order=3);
  template
  void dyson_mat_fourier<double>(herm_matrix<double> &G,herm_matrix<double> &Sigma,double mu,function<double> &H,
  				 function<double> &SigmaMF, double beta,int order);

  template
  void dyson_mat_fixpoint<double>(herm_matrix<double> &G,herm_matrix<double> &Sigma,double mu,function<double> &H,
  				  integration::Integrator<double> &I, double beta,int fixpiter);
  template
  void dyson_mat_fixpoint<double>(herm_matrix<double> &G,herm_matrix<double> &Sigma,double mu,function<double> &H,
  				  function<double> &SigmaMF,
  				  integration::Integrator<double> &I, double beta,int fixpiter);

  template
  void dyson_mat_steep<double>(herm_matrix<double> &G, herm_matrix<double> &Sigma, double mu, function<double> &H,
  			       integration::Integrator<double> &I, double beta, int maxiter, double tol);
  template
  void dyson_mat_steep<double>(herm_matrix<double> &G, herm_matrix<double> &Sigma, double mu, function<double> &H,
  			       function<double> &SigmaMF, integration::Integrator<double> &I, double beta,
  			       int maxiter, double tol);
  template
  void dyson_mat<double>(herm_matrix<double> &G, herm_matrix<double> &Sigma, double mu, function<double> &H,
  			 integration::Integrator<double> &I, double beta, const int method,
         const bool force_hermitian);
  template
  void dyson_mat<double>(herm_matrix<double> &G, herm_matrix<double> &Sigma, double mu, function<double> &H,
  			 function<double> &SigmaMF, integration::Integrator<double> &I, double beta,
			      const int method, const bool force_hermitian);


  template
  void dyson_start<double>(herm_matrix<double> &G,double mu,function<double> &H, herm_matrix<double> &Sigma,
  			   integration::Integrator<double> &I, double beta,double h);
  template
  void dyson_timestep<double>(int n,herm_matrix<double> &G,double mu,function<double> &H, herm_matrix<double> &Sigma,
  			      integration::Integrator<double> &I, double beta,double h);
  template
  void dyson<double>(herm_matrix<double> &G,double mu,function<double> &H, herm_matrix<double> &Sigma,
  		     integration::Integrator<double> &I, double beta,double h, const int matsubara_method,
          const bool force_hermitian);


// documented user interfaces
 template
  void dyson_mat<double>(herm_matrix<double> &G, double mu, function<double> &H, herm_matrix<double> &Sigma, double beta, 
     const int SolveOrder, const int method,
     const bool force_hermitian);

 template
  void dyson_mat<double>(herm_matrix<double> &G, double mu, function<double> &H, function<double> &SigmaMF, double beta, 
     herm_matrix<double> &Sigma, const int SolveOrder, const int method,
     const bool force_hermitian);

 template
  void dyson_start<double>(herm_matrix<double> &G,  double mu, function<double> &H, herm_matrix<double> &Sigma, double beta,double h, 
    const int SolveOrder);

 template
  void dyson_timestep<double>(int n, herm_matrix<double> &G, double mu, function<double> &H, herm_matrix<double> &Sigma, double beta, double h, 
    const int SolveOrder);

 template
  void dyson<double>(herm_matrix<double> &G, double mu, function<double> &H, herm_matrix<double> &Sigma, double beta, double h, 
    const int SolveOrder, const int matsubara_method,
    const bool force_hermitian);

}  // namespace cntr
