#ifndef CNTR_RESPONSE_CONVOLUTION_DECL_H
#define CNTR_RESPONSE_CONVOLUTION_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

/*###########################################################################################
#
#   CONVOLUTION  c=W*f   with a response function
#
#	c(t) =  \int_C  dt' W_{a1,a2}(t,t') f_{b1,b2}(t')
#
#	- Assuming that W is hermitian, i.e.,
#     Wles_{a1,a2}(t',t) = -Wles_{a2,a1}(t,t')^*  etc.
#	- needs timesteps 0...kt of W for t<=kt
#     => the bersion using only a herm_matrix_timestep works only for tstp=-1 and tstp>=kt
#
###########################################################################################*/


/** \brief <b> Evaluate the convolution between a two-time contour function (\f$W\f$) and a real-time function (\f$f\f$) at the time step (\f$t\f$); \f$ c(t) =  \int_C  dt' W_{a1,a2}(t,t') f_{b1,b2}(t') \f$. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * This function evaluates the convolution between a two-time contour function (\f$W\f$) and a real-time function (\f$f\f$) at the time step (\f$t\f$); \f$ c(t) =  \int_C  dt' W_{a1,a2}(t,t') f_{b1,b2}(t') \f$.
 * It can be used, for example, to evaluate a lineaer response by using a response function as \f$W\f$ and an external field as \f$f\f$.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > Time step.
 * @param cc
 * > Value of the conovlution at the time step; \f$ c(t) =  \int_C  dt' W_{a1,a2}(t,t') f_{b1,b2}(t') \f$.
 * @param W
 * > Two-time contour object in the Matrix form with the hermitian symmetry.
 * @param a1
 * > First index of the matrix \f$W\f$.
 * @param a2
 * > Second index of the matrix \f$W\f$.
 * @param f
 * > Function in the Matrix from on the real-time axis.
 * @param b1
 * > First index of the matrix \f$f\f$.
 * @param b2
 * > Second index of the matrix \f$f\f$.
 * @param kt
 * > Integration order
 * @param beta
 * > Inversed temperature
 * @param h
 * > time step interval
 */
template <typename T>
void response_convolution(int tstp, std::complex<T> &cc, herm_matrix<T> &W, int a1, int a2,
                          function<T> &f, int b1, int b2, int kt, T beta, T h);

/** \brief <b> Evaluate the convolution between a two-time contour function (\f$W\f$) and a real-time function (\f$f\f$) at the time step (\f$t\f$); \f$ c(t) =  \int_C  dt' W_{a1,a2}(t,t') f_{b1,b2}(t') \f$. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * This function evaluates the convolution between a two-time contour function (\f$W\f$) and a real-time function (\f$f\f$) at the time step (\f$t\f$); \f$ c(t) =  \int_C  dt' W_{a1,a2}(t,t') f_{b1,b2}(t') \f$.
 * It can be used, for example, to evaluate a lineaer response by using a response function as \f$W\f$ and an external field as \f$f\f$.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > Time step.
 * @param cc
 * > Value of the conovlution at the time step; \f$ c(t) =  \int_C  dt' W_{a1,a2}(t,t') f_{b1,b2}(t') \f$.
 * @param W
 * > Two-time contour object in the Matrix form with the hermitian symmetry.
 * @param a1
 * > First index of the matrix \f$W\f$.
 * @param a2
 * > Second index of the matrix \f$W\f$.
 * @param f
 * > Function in the Matrix from on the real-time axis.
 * @param b1
 * > First index of the matrix \f$f\f$.
 * @param b2
 * > Second index of the matrix \f$f\f$.
 * @param kt
 * > Integration order
 * @param beta
 * > Inversed temperature
 * @param h
 * > time step interval
 */
template <typename T>
void response_convolution(int tstp, std::complex<T> &cc, herm_matrix_timestep<T> &W, int a1,
                          int a2, function<T> &f, int b1, int b2, int kt, T beta, T h);

} // namespace cntr

#endif  // CNTR_RESPONSE_CONVOLUTION_DECL_H
