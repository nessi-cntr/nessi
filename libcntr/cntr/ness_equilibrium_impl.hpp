#ifndef NESS_EQUILIBRIUM_IMPL_H
#define CNTR_EQUILIBRIUM_IMPL_H

#include "ness_global_settings.hpp"
#include "ness_equilibrium_decl.hpp"
#include "cntr/cntr.hpp"

using namespace std;
using namespace cntr;

namespace ness
{

/// @private
/** \brief Default constructor of gaussian DOS, mean `mu_=0` and standard deviation `sigma=1`  */
gauss::gauss() {
    mu_ = 0;
    sigma_=1;
    lo_ = -6*sigma_;
    hi_ = 6*sigma_;
}
/// @private
/** \brief  Constructor of gaussian DOS, mean `mu_` and standard deviation `sigma`  */
gauss::gauss(double mu, double sigma) {
    mu_ = mu;
    sigma_ = sigma;
    lo_ = -6*sigma;
    hi_ = 6*sigma;
}
/// @private
double gauss::operator()(double x){
    return std::exp(-0.5*((x-mu_)/sigma_)*((x-mu_)/sigma_))/(sigma_*2.5066);
}


/// @private
/** \brief Default constructor of ohmic DOS, cutoff `omegac_=1`   */
ohmic_sym::ohmic_sym() {
    omegac_ = 1;
    lo_ = 0.01;
    hi_ = omegac_*20;
}
/// @private
/** \brief Constructor of ohmic DOS, cutoff `omegac_`   */
ohmic_sym::ohmic_sym(double omegac) {
    omegac_ = omegac;
    lo_=-omegac*20;
    hi_=omegac*20; // value is 10^-6
}
/// @private
double ohmic_sym::operator()(double x) {
    double abs_x = (x >= 0) ? x : -x;
    double result = abs_x * abs_x * std::exp(-abs_x / omegac_) / (2.0 * omegac_ * omegac_ * omegac_);
    return (x >= 0) ? result : -result;
}

} //namespace ness
#endif  // NESS_EQUILIBRIUM_IMPL_H
