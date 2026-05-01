#ifndef NESS2_EQUILIBRIUM_IMPL_H
#define NESS2_EQUILIBRIUM_IMPL_H

#include "ness2_global_settings.hpp"
#include "ness2_equilibrium_decl.hpp"

namespace ness2{

  /// @private
  gauss::gauss() {
    mu_ = 0;
    sigma_=1;
    lo_ = mu_-6*sigma_;
    hi_ = mu_+6*sigma_;
  }
  /// @private
  gauss::gauss(double mu, double sigma) {
    mu_ = mu;
    sigma_ = sigma;
    lo_ = mu_-6*sigma;
    hi_ = mu_+6*sigma;
  }
  /// @private
  double gauss::operator()(double x){
    return std::exp(-0.5*((x-mu_)/sigma_)*((x-mu_)/sigma_))/(sigma_*2.506623299999997);
  }

  
  /// @private
  ohmic_sym::ohmic_sym() {
    omegac_ = 1;
    lo_ = -omegac_*20;
    hi_ = omegac_*20;
  }
  /// @private
  ohmic_sym::ohmic_sym(double omegac) {
    omegac_ = omegac;
    lo_=-omegac*20;
    hi_=omegac*20; // value is 10^-6
  }
  /// @private
  double ohmic_sym::operator()(double x) {
    double abs_x = (x >= 0) ? x : -x;
    double result = abs_x  * std::exp(-abs_x / omegac_) / (2.0 * omegac_ * omegac_ );
    return (x >= 0) ? result : -result;
  }
  
  
  bethedos::bethedos() {
    V_ = 1;
    lo_ = -2;
    hi_ = 2;
    E0_ = 0;
  }
  bethedos::bethedos(double a, double b){
    lo_ = a;
    hi_ = b;
    V_ = (b-a)/4.0;
    E0_ = (b+a)/2.0;
  }
  double bethedos::operator()(double x) {
    double arg = 4.0 * V_ * V_ - (x-E0_) * (x-E0_);
    double num = V_ * V_ * 3.14159265358979323846 * 2;
    return (arg < 0 ? 0.0 : sqrt(arg) / num);
  }
  
  
} //namespace ness
#endif  // NESS_EQUILIBRIUM_IMPL_H
