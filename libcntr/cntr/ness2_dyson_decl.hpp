#ifndef NESS2_DYSON_DECL_H
#define NESS2_DYSON_DECL_H

#include "ness2_global_settings.hpp"

namespace ness2{
    
  enum dyson_bath_type {BATH_NONE=0,BATH_CONST,BATH_GAUSS,BATH_OHMIC};


  /**
   * \brief Solve the Dyson equation for the steady-state Green's function.
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   * \par Purpose
   * > Computes the time-domain solution \f$ G^{R|<}(t) \f$ of the Dyson equation for a steady-state system,
   *   using a frequency-domain method with optional bath regularization.
   * > The method transforms \f$ \Sigma^{R|<}(t) \f$ to frequency space, solves Dyson equations there,
   *   and back-transforms to obtain \f$ G^{R|<}(t) \f$.
   *
   * <!-- ========= ARGUMENTS ========= -->
   *
   * @param[out] G
   *  Output Green's function \f$ G^{R|<}(t) \f$ of type `herm_matrix_ness`, modified in-place.
   *
   * @param[in] mu
   *  Chemical potential \f$ \mu \f$.
   *
   * @param[in] epsilon
   *  Complex-valued single-particle Hamiltonian matrix \f$ \epsilon \f$ (dimension must match `G`).
   *
   * @param[in] Sigma
   *  Input self-energy \f$ \Sigma^{R|<}(t) \f$ of type `herm_matrix_ness`, defined on the same grid as `G`.
   *
   * @param[in] h
   *  Time step size \f$ h \f$ used for time grids. 
   *
   * @param[in] method
   *  Integral transformation method (default = `FFT_TRAPEZ`). Currently only `FFT_TRAPEZ` is implemented.
   *
   * @param[in] bath_type
   *  Type of bath regularization added to the self-energy. 
   *  If `bath_type != REG_NONE`, a bath self-energy \f$ \Sigma_{\rm reg} \f$ is added to \f$ \Sigma \f$:
   *  - `REG_CONST`: Constant density (Lorentzian broadening)
   *  - `REG_GAUSS`: Gaussian cutoff in frequency space.
   *  - `REG_OHMIC`: Ohmic spectral function with cutoff.
   *
   * @param[in] bath_eta
   *  Broadening parameter \f$ \eta \f$ for the bath (used for all regularization types except `BATH_NONE`).
   *
   * @param[in] bath_beta
   *  Inverse temperature \f$ \beta \f$ of the bath (used in `BATH_CONST`, `BATH_GAUSS`, `BATH_OHMIC`).
   *
   * @param[in] bath_mu
   *  Chemical potential \f$ \mu_{\rm bath} \f$ of the bath (used in `BATH_CONST`, `BATH_GAUSS`).
   *
   */ 
  void dyson(herm_matrix_ness &G, double mu, cdmatrix &epsilon, herm_matrix_ness &Sigma, double h,
	     fft_integral_method method=FFT_TRAPEZ,
	     dyson_bath_type bath_type=BATH_NONE,double bath_eta=0,double bath_beta=0,double bath_mu=0);
    
  
  //@private
  void dyson_fft(herm_matrix_ness  & G, herm_matrix_ness  & Sigma, cdmatrix & H,double mu,double h,
		 dyson_bath_type bath_type,double bath_eta,double bath_beta,double bath_mu);

  //@private
  void sigma_bath(double &ret,double &les, double omega,double h,dyson_bath_type bath_type,double bath_eta,double bath_beta,double bath_mu);


    
} //namespace ness
#endif  // NESS_EQUILIBRIUM_IMPL_H
