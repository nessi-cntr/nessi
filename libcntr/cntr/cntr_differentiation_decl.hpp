#ifndef CNTR_DIFFERENTIATION_DECL_H
#define CNTR_DIFFERENTIATION_DECL_H

#include "cntr_global_settings.hpp"
#include "integration.hpp"

namespace cntr {

template <typename T> class herm_matrix;

/*###########################################################################################
#
#   DIFFERENTIATION  dA(t,t') = id/dt A(t,t')   or  dA(t,t') = -id/dt' A(t,t')
#
#   ...  currently only implemented for herm_matrix
#   to compute differentiation of contour functions assuming they are known for all times
#   Matsubara need still working midpoint differentiation (now all bwd except small times
(fwd))
#
###########################################################################################*/
/// @private
template <typename T>
void deriv1_matsubara(herm_matrix<T> &dA, herm_matrix<T> &A, integration::Integrator<T> &I,
                      T beta);
/// @private
template <typename T>
void deriv1_element(int tstp1, int tstp2, herm_matrix<T> &dA, herm_matrix<T> &A,
                    herm_matrix<T> &Acc, integration::Integrator<T> &I, T h);
/// @private
template <typename T>
void deriv2_element(int tstp1, int tstp2, herm_matrix<T> &dA, herm_matrix<T> &A,
                    herm_matrix<T> &Acc, integration::Integrator<T> &I, T h);
/// @private
template <typename T>
void deriv1_tv(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A,
               integration::Integrator<T> &I, T h);
/// @private
template <typename T>
void deriv2_tv(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A,
               integration::Integrator<T> &I, T beta);
/// @private
template <typename T>
void deriv1_timestep(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A, herm_matrix<T> &Acc,
                     integration::Integrator<T> &I, T beta, T h);
/// @private
template <typename T>
void deriv2_timestep(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A, herm_matrix<T> &Acc,
                     integration::Integrator<T> &I, T beta, T h);


template <typename T>
void deriv1_matsubara(herm_matrix<T> &dA, herm_matrix<T> &A, T beta, int SolveOrder=MAX_SOLVE_ORDER);
template <typename T>
void deriv1_element(int tstp1, int tstp2, herm_matrix<T> &dA, herm_matrix<T> &A,
                    herm_matrix<T> &Acc, T h, int SolveOrder=MAX_SOLVE_ORDER);
template <typename T>
void deriv2_element(int tstp1, int tstp2, herm_matrix<T> &dA, herm_matrix<T> &A,
                    herm_matrix<T> &Acc, T h, int SolveOrder=MAX_SOLVE_ORDER);
template <typename T>
void deriv1_tv(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A,
               T h, int SolveOrder=MAX_SOLVE_ORDER);
template <typename T>
void deriv2_tv(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A,
               T beta, int SolveOrder=MAX_SOLVE_ORDER);
template <typename T>
void deriv1_timestep(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A, herm_matrix<T> &Acc,
                     T beta, T h, int SolveOrder=MAX_SOLVE_ORDER);
template <typename T>
void deriv2_timestep(int tstp, herm_matrix<T> &dA, herm_matrix<T> &A, herm_matrix<T> &Acc,
                     T beta, T h, int SolveOrder=MAX_SOLVE_ORDER);

}  // namespace cntr

#endif  // CNTR_DIFFERENTIATION_DECL_H
