#ifndef CNTR_PSEUDO_VIE2_DECL_H
#define CNTR_PSEUDO_VIE2_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

/* #######################################################################################
#  [1+F]G=Q, G,Q hermitian
###########################################################################################*/
/// @private
template <typename T>
void pseudo_vie2_mat(herm_pseudo<T> &G, herm_pseudo<T> &F, herm_pseudo<T> &Fcc,
                     herm_pseudo<T> &Q, integration::Integrator<T> &I, T beta);
/// @private
template <typename T>
void pseudo_vie2_start(herm_pseudo<T> &G, herm_pseudo<T> &F, herm_pseudo<T> &Fcc,
                       herm_pseudo<T> &Q, integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T>
void pseudo_vie2_timestep(int n, herm_pseudo<T> &G, herm_pseudo<T> &F, herm_pseudo<T> &Fcc,
                          herm_pseudo<T> &Q, integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T>
void pseudo_vie2(herm_pseudo<T> &G, herm_pseudo<T> &F, herm_pseudo<T> &Fcc,
                 herm_pseudo<T> &Q, integration::Integrator<T> &I, T beta, T h);

}  // namespace cntr

#endif  // CNTR_PSEUDO_VIE2_DECL_H
