#ifndef CNTR_DYSON_EXTERN_TEMPLATES_H
#define CNTR_DYSON_EXTERN_TEMPLATES_H

#include "cntr_dyson_decl.hpp"

namespace cntr {

// internal interfaces
  extern template
  void dyson_mat_fourier<double>(herm_matrix<double> &G,herm_matrix<double> &Sigma,double mu,function<double> &H,double beta,int order);
  extern template
  void dyson_mat_fourier<double>(herm_matrix<double> &G,herm_matrix<double> &Sigma,double mu,function<double> &H,
				 function<double> &SigmaMF, double beta,int order);


  extern template
  void dyson_mat_fixpoint<double>(herm_matrix<double> &G,herm_matrix<double> &Sigma,double mu,function<double> &H,
						  integration::Integrator<double> &I, double beta,int fixpiter);
  extern template
  void dyson_mat_fixpoint<double>(herm_matrix<double> &G,herm_matrix<double> &Sigma,double mu,function<double> &H,
						  function<double> &SigmaMF,
						  integration::Integrator<double> &I, double beta,int fixpiter);


  extern template
  void dyson_mat_steep<double>(herm_matrix<double> &G, herm_matrix<double> &Sigma, double mu, function<double> &H,
			       integration::Integrator<double> &I, double beta, int maxiter, double tol);
  extern template
  void dyson_mat_steep<double>(herm_matrix<double> &G, herm_matrix<double> &Sigma, double mu, function<double> &H,
			       function<double> &SigmaMF, integration::Integrator<double> &I, double beta,
			       int maxiter, double tol);

  extern template
  void dyson_mat<double>(herm_matrix<double> &G, herm_matrix<double> &Sigma, double mu, function<double> &H,
			 integration::Integrator<double> &I, double beta, const int method, const bool force_hermitian);
  extern template
  void dyson_mat<double>(herm_matrix<double> &G, herm_matrix<double> &Sigma, double mu, function<double> &H,
			 function<double> &SigmaMF, integration::Integrator<double> &I, double beta, const int method,
      const bool force_hermitian);


  extern template
  void dyson_start<double>(herm_matrix<double> &G,double mu,function<double> &H, herm_matrix<double> &Sigma,
			   integration::Integrator<double> &I, double beta,double h);
  extern template
  void dyson_timestep<double>(int n,herm_matrix<double> &G,double mu,function<double> &H, herm_matrix<double> &Sigma,
			      integration::Integrator<double> &I, double beta,double h);
  extern template
  void dyson<double>(herm_matrix<double> &G,double mu,function<double> &H, herm_matrix<double> &Sigma,
		     integration::Integrator<double> &I, double beta,double h, const int matsubara_method,
        const bool force_hermitian);

// documented user interfaces
  extern template
  void dyson_mat<double>(herm_matrix<double> &G, double mu, function<double> &H, herm_matrix<double> &Sigma, double beta, 
     const int SolveOrder, const int method,
     const bool force_hermitian);

  extern template
  void dyson_mat<double>(herm_matrix<double> &G, double mu, function<double> &H, function<double> &SigmaMF,
     herm_matrix<double> &Sigma, double beta, const int SolveOrder, const int method,
     const bool force_hermitian);

  extern template
  void dyson_start<double>(herm_matrix<double> &G,  double mu, function<double> &H, herm_matrix<double> &Sigma, double beta,double h, 
    const int SolveOrder);

  extern template
  void dyson_timestep<double>(int n, herm_matrix<double> &G, double mu, function<double> &H, herm_matrix<double> &Sigma, double beta, double h, 
    const int SolveOrder);

  extern template
  void dyson<double>(herm_matrix<double> &G, double mu, function<double> &H, herm_matrix<double> &Sigma, double beta, double h, 
    const int SolveOrder, const int matsubara_method,
    const bool force_hermitian);

}  // namespace cntr

#endif  // CNTR_DYSON_EXTERN_TEMPLATES_H
