#ifndef CNTR_BUBBLE_DECL_H
#define CNTR_BUBBLE_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

/*#########################################################################################
 #
 #   ...   useful routines to compute diagrams
 #
 #   BUBBLE1:  C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b2,b1}(t2,t1)
 #
 #   BUBBLE2:  C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b1,b2}(t1,t2)
 #
 #   template GGA,GGB,GGC can be herm_matrix,herm_matrix_timestep,herm_matrix_timestep_view
 #
 #   (implementation in cntr_diagram_utilities.hpp, easy to add further calls with timestep
 #    instead of green functions)
 #########################################################################################*/


/** \brief <b> Evaluate a bubble diagram (\f$C\f$) from two-time contour functions \f$A,B\f$ at the time step; \f$ C_{c_1,c_2}(t_1,t_2) = i  A_{a_1,a_2}(t_t,t_2) * B_{b_2,b_1}(t_2,t_1) \f$. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 * Evaluate the two-time contour function \f$C\f$ represented as a bubble diagram with two-time functions \f$A,B\f$; \f$ C_{c_1,c_2}(t_1,t_2) = i  A_{a_1,a_2}(t_t,t_2) * B_{b_2,b_1}(t_2,t_1) \f$.
 * This evaluation is done at the time step (i.e. \f$ t_1 \f$ or \f$ t_2 \f$ is the time step) for all components (retarded, lesser, left-mixing and Matsubara).
 * The evaluated value of \f$ i  A_{a_1,a_2}(t_t,t_2) * B_{b_2,b_1}(t_2,t_1)\f$ is stored at \f$ C_{c_1,c_2}(t_1,t_2)\f$.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > Time step.
 * @param C
 * > Two-time contour object in the Matrix form defined as \f$ C_{c_1,c_2}(t_1,t_2) = i  A_{a_1,a_2}(t_t,t_2) * B_{b_2,b_1}(t_2,t_1) \f$.
 * @param c1
 * > First index of the matrix \f$C\f$.
 * @param c2
 * > Second index of the matrix \f$C\f$.
 * @param A
 * > Two-time contour object in the Matrix form.
 * @param Acc
 * > Two-time contour object in the Matrix form, which is adjoint of \f$A\f$.
 * @param a1
 * > First index of the matrix \f$A\f$.
 * @param a2
 * > Second index of the matrix \f$A\f$.
 * @param B
 * > Two-time contour object in the Matrix form.
 * @param Bcc
 * > Two-time contour object in the Matrix form, which is adjoint of \f$B\f$.
 * @param b1
 * > 'Second' index of the matrix \f$B\f$.
 * @param b2
 * > 'First' index of the matrix \f$B\f$.
 */

template <class GGC, class GGA, class GGB>
void Bubble1(int tstp, GGC &C, int c1, int c2, GGA &A, GGA &Acc, int a1, int a2, GGB &B,
             GGB &Bcc, int b1, int b2);


/** \brief <b> Evaluate a bubble diagram (\f$C\f$) from two-time contour functions \f$A,B\f$ at the time step; \f$ C_{c_1,c_2}(t_1,t_2) = i  A_{a_1,a_2}(t_t,t_2) * B_{b_2,b_1}(t_2,t_1) \f$. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 * Evaluate the two-time contour function \f$C\f$ represented as a bubble diagram with two-time functions \f$A,B\f$; \f$ C_{c_1,c_2}(t_1,t_2) = i  A_{a_1,a_2}(t_t,t_2) * B_{b_2,b_1}(t_2,t_1) \f$.
 * This evaluation is done at the time step (i.e. \f$ t_1 \f$ or \f$ t_2 \f$ is the time step) for all components (retarded, lesser, left-mixing and Matsubara).
 * The evaluated value of \f$ i  A_{a_1,a_2}(t_t,t_2) * B_{b_2,b_1}(t_2,t_1)\f$ is stored at \f$ C_{c_1,c_2}(t_1,t_2)\f$.
 * Here \f$A\f$ and \f$B\f$ are assumed to have the hermitian symmetry.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > Time step.
 * @param C
 * > Two-time contour object in the Matrix form defined as \f$ C_{c_1,c_2}(t_1,t_2) = i  A_{a_1,a_2}(t_t,t_2) * B_{b_2,b_1}(t_2,t_1) \f$.
 * @param c1
 * > First index of the matrix \f$C\f$.
 * @param c2
 * > Second index of the matrix \f$C\f$.
 * @param A
 * > Two-time contour object in the Matrix form with the hermitian symmetry.
 * @param a1
 * > First index of the matrix \f$A\f$.
 * @param a2
 * > Second index of the matrix \f$A\f$.
 * @param B
 * > Two-time contour object in the Matrix form with the hermitian symmetry.
 * @param b1
 * > 'Second' index of the matrix \f$B\f$.
 * @param b2
 * > 'First' index of the matrix \f$B\f$.
 */

template <class GGC, class GGA, class GGB>
void Bubble1(int tstp, GGC &C, int c1, int c2, GGA &A, int a1, int a2, GGB &B, int b1,
             int b2);

/** \brief <b> Evaluate a bubble diagram (\f$C\f$) from two-time contour functions \f$A,B\f$ at the time step; \f$ C_{0,0}(t_1,t_2) = i  A_{0,0}(t_t,t_2) * B_{0,0}(t_2,t_1) \f$. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 * Evaluate the two-time contour function \f$C\f$ represented as a bubble diagram with two-time functions \f$A,B\f$; \f$ C_{0,0}(t_1,t_2) = i  A_{0,0}(t_t,t_2) * B_{0,0}(t_2,t_1) \f$.
 * Here it is assmued that \f$A\f$, \f$B\f$ and \f$C\f$ are \f$1\times 1\f$ matrices.
 * This evaluation is done at the time step (i.e. \f$ t_1 \f$ or \f$ t_2 \f$ is the time step) for all components (retarded, lesser, left-mixing and Matsubara).
 * The evaluated value of \f$ i  A_{0,0}(t_t,t_2) * B_{0,0}(t_2,t_1)\f$ is stored at \f$ C_{0,0}(t_1,t_2)\f$.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > Time step.
 * @param C
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form defined as \f$ C_{0,0}(t_1,t_2) = i  A_{0,0}(t_t,t_2) * B_{0,0}(t_2,t_1) \f$.
 * @param A
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form.
 * @param Acc
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form, which is adjoint of \f$A\f$.
 * @param B
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form.
 * @param Bcc
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form, which is adjoint of \f$B\f$.
 */

template <class GGC, class GGA, class GGB>
void Bubble1(int tstp, GGC &C, GGA &A, GGA &Acc, GGB &B, GGB &Bcc);

/** \brief <b> Evaluate a bubble diagram (\f$C\f$) from two-time contour functions \f$A,B\f$ at the time step; \f$ C_{0,0}(t_1,t_2) = i  A_{0,0}(t_t,t_2) * B_{0,0}(t_2,t_1) \f$. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 * Evaluate the two-time contour function \f$C\f$ represented as a bubble diagram with two-time functions \f$A,B\f$; \f$ C_{0,0}(t_1,t_2) = i  A_{0,0}(t_t,t_2) * B_{0,0}(t_2,t_1) \f$.
 * Here it is assmued that \f$A\f$, \f$B\f$ and \f$C\f$ are \f$1\times 1\f$ matrices.
 * This evaluation is done at the time step (i.e. \f$ t_1 \f$ or \f$ t_2 \f$ is the time step) for all components (retarded, lesser, left-mixing and Matsubara).
 * The evaluated value of \f$ i  A_{0,0}(t_t,t_2) * B_{0,0}(t_2,t_1)\f$ is stored at \f$ C_{0,0}(t_1,t_2)\f$.
 * Here \f$A\f$ and \f$B\f$ are assumed to have the hermitian symmetry.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > Time step.
 * @param C
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form defined as \f$ C_{0,0}(t_1,t_2) = i  A_{0,0}(t_t,t_2) * B_{0,0}(t_2,t_1) \f$.
 * @param A
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form with the hermitian symmetry.
 * @param B
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form with the hermitian symmetry.
 */
template <class GGC, class GGA, class GGB> void Bubble1(int tstp, GGC &C, GGA &A, GGB &B);





////////////

/** \brief <b> Evaluate a bubble diagram (\f$C\f$) from two-time contour functions \f$A,B\f$ at the time step; \f$ C_{c_1,c_2}(t_1,t_2) = i  A_{a_1,a_2}(t_1,t_2) * B_{b_1,b_2}(t_1,t_2) \f$. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 * Evaluate the two-time contour function \f$C\f$ represented as a bubble diagram with two-time functions \f$A,B\f$; \f$ C_{c_1,c_2}(t_1,t_2) = i  A_{a_1,a_2}(t_1,t_2) * B_{b_1,b_2}(t_1,t_2) \f$.
 * This evaluation is done at the time step (i.e. \f$ t_1 \f$ or \f$ t_2 \f$ is the time step) for all components (retarded, lesser, left-mixing and Matsubara).
 * The evaluated value of \f$ i  A_{a_1,a_2}(t_1,t_2) * B_{b_1,b_2}(t_1,t_2) \f$ is stored at \f$ C_{c_1,c_2}(t_1,t_2)\f$.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > Time step.
 * @param C
 * > Two-time contour object in the Matrix form defined as \f$ C_{c_1,c_2}(t_1,t_2) = i  A_{a_1,a_2}(t_1,t_2) * B_{b_1,b_2}(t_1,t_2) \f$.
 * @param c1
 * > First index of the matrix \f$C\f$.
 * @param c2
 * > Second index of the matrix \f$C\f$.
 * @param A
 * > Two-time contour object in the Matrix form.
 * @param Acc
 * > Two-time contour object in the Matrix form, which is adjoint of \f$A\f$.
 * @param a1
 * > First index of the matrix \f$A\f$.
 * @param a2
 * > Second index of the matrix \f$A\f$.
 * @param B
 * > Two-time contour object in the Matrix form.
 * @param Bcc
 * > Two-time contour object in the Matrix form, which is adjoint of \f$B\f$.
 * @param b1
 * > First index of the matrix \f$B\f$.
 * @param b2
 * > Second index of the matrix \f$B\f$.
 */

template <class GGC, class GGA, class GGB>
void Bubble2(int tstp, GGC &C, int c1, int c2, GGA &A, GGA &Acc, int a1, int a2, GGB &B,
             GGB &Bcc, int b1, int b2);

/** \brief <b> Evaluate a bubble diagram (\f$C\f$) from two-time contour functions \f$A,B\f$ at the time step; \f$ C_{c_1,c_2}(t_1,t_2) = i  A_{a_1,a_2}(t_1,t_2) * B_{b_1,b_2}(t_1,t_2) \f$. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 * Evaluate the two-time contour function \f$C\f$ represented as a bubble diagram with two-time functions \f$A,B\f$; \f$ C_{c_1,c_2}(t_1,t_2) = i  A_{a_1,a_2}(t_1,t_2) * B_{b_1,b_2}(t_1,t_2) \f$.
 * This evaluation is done at the time step (i.e. \f$ t_1 \f$ or \f$ t_2 \f$ is the time step) for all components (retarded, lesser, left-mixing and Matsubara).
 * The evaluated value of \f$ i  A_{a_1,a_2}(t_1,t_2) * B_{b_1,b_2}(t_1,t_2)\f$ is stored at \f$ C_{c_1,c_2}(t_1,t_2)\f$.
 * Here \f$A\f$ and \f$B\f$ are assumed to have the hermitian symmetry.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > Time step.
 * @param C
 * > Two-time contour object in the Matrix form defined as \f$ C_{c_1,c_2}(t_1,t_2) = i  A_{a_1,a_2}(t_1,t_2) * B_{b_1,b_2}(t_1,t_2) \f$.
 * @param c1
 * > First index of the matrix \f$C\f$.
 * @param c2
 * > Second index of the matrix \f$C\f$.
 * @param A
 * > Two-time contour object in the Matrix form with the hermitian symmetry.
 * @param a1
 * > First index of the matrix \f$A\f$.
 * @param a2
 * > Second index of the matrix \f$A\f$.
 * @param B
 * > Two-time contour object in the Matrix form with the hermitian symmetry.
 * @param b1
 * > First index of the matrix \f$B\f$.
 * @param b2
 * > Second index of the matrix \f$B\f$.
 */

template <class GGC, class GGA, class GGB>
void Bubble2(int tstp, GGC &C, int c1, int c2, GGA &A, int a1, int a2, GGB &B, int b1,
             int b2);

/** \brief <b> Evaluate a bubble diagram (\f$C\f$) from two-time contour functions \f$A,B\f$ at the time step; \f$ C_{0,0}(t_1,t_2) = i  A_{0,0}(t_1,t_2) * B_{0,0}(t_1,t_2) \f$. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 * Evaluate the two-time contour function \f$C\f$ represented as a bubble diagram with two-time functions \f$A,B\f$; \f$ C_{0,0}(t_1,t_2) = i  A_{0,0}(t_1,t_2) * B_{0,0}(t_1,t_2) \f$.
 * Here it is assmued that \f$A\f$, \f$B\f$ and \f$C\f$ are \f$1\times 1\f$ matrices.
 * This evaluation is done at the time step (i.e. \f$ t_1 \f$ or \f$ t_2 \f$ is the time step) for all components (retarded, lesser, left-mixing and Matsubara).
 * The evaluated value of \f$ i  A_{0,0}(t_1,t_2) * B_{0,0}(t_1,t_2)\f$ is stored at \f$ C_{0,0}(t_1,t_2)\f$.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > Time step.
 * @param C
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form defined as \f$ C_{0,0}(t_1,t_2) = i  A_{0,0}(t_1,t_2) * B_{0,0}(t_1,t_2) \f$.
 * @param A
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form.
 * @param Acc
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form, which is adjoint of \f$A\f$.
 * @param B
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form.
 * @param Bcc
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form, which is adjoint of \f$B\f$.
 */

template <class GGC, class GGA, class GGB>
void Bubble2(int tstp, GGC &C, GGA &A, GGA &Acc, GGB &B, GGB &Bcc);

/** \brief <b> Evaluate a bubble diagram (\f$C\f$) from two-time contour functions \f$A,B\f$ at the time step; \f$ C_{0,0}(t_1,t_2) = i  A_{0,0}(t_1,t_2) * B_{0,0}(t_1,t_2) \f$. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 * Evaluate the two-time contour function \f$C\f$ represented as a bubble diagram with two-time functions \f$A,B\f$; \f$ C_{0,0}(t_1,t_2) = i  A_{0,0}(t_1,t_2) * B_{0,0}(t_1,t_2) \f$.
 * Here it is assmued that \f$A\f$, \f$B\f$ and \f$C\f$ are \f$1\times 1\f$ matrices.
 * This evaluation is done at the time step (i.e. \f$ t_1 \f$ or \f$ t_2 \f$ is the time step) for all components (retarded, lesser, left-mixing and Matsubara).
 * The evaluated value of \f$ i  A_{0,0}(t_1,t_2) * B_{0,0}(t_1,t_2)\f$ is stored at \f$ C_{0,0}(t_1,t_2)\f$.
 * Here \f$A\f$ and \f$B\f$ are assumed to have the hermitian symmetry.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > Time step.
 * @param C
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form defined as \f$ C_{0,0}(t_1,t_2) = i  A_{0,0}(t_1,t_2) * B_{0,0}(t_1,t_2) \f$.
 * @param A
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form with the hermitian symmetry.
 * @param B
 * > Two-time contour object in the \f$1\times 1\f$-Matrix form with the hermitian symmetry.
 */

template <class GGC, class GGA, class GGB> void Bubble2(int tstp, GGC &C, GGA &A, GGB &B);
///////////

///////////

// #if 0
// template<typename T>
// void Bubble1(int tstp, herm_matrix<T> &C,int c1,int c2, herm_matrix<T> &A,herm_matrix<T> &Acc,int a1,int a2,
// herm_matrix<T> &B,herm_matrix<T> &Bcc,int b1,int b2);
// template<typename T>
// void Bubble1(int tstp, herm_matrix<T> &C,int c1,int c2, herm_matrix<T> &A,int a1,int a2,herm_matrix<T> &B,int b1,int b2);
// template<typename T>
// void Bubble1(int tstp, herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,herm_matrix<T> &B, herm_matrix<T> &Bcc);
// template<typename T>
// void Bubble1(int tstp, herm_matrix<T> &C, herm_matrix<T> &A,herm_matrix<T> &B);
// template<typename T>
// void Bubble1(int tstp, herm_matrix_timestep<T> &C,int c1,int c2, herm_matrix<T> &A,herm_matrix<T> &Acc,int a1,int a2,
// herm_matrix<T> &B,herm_matrix<T> &Bcc,int b1,int b2);
// template<typename T>
// void Bubble1(int tstp, herm_matrix_timestep<T> &C,int c1,int c2, herm_matrix<T> &A,int a1,int a2,herm_matrix<T> &B,int b1,int b2);
// template<typename T>
// void Bubble1(int tstp, herm_matrix_timestep<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,herm_matrix<T> &B, herm_matrix<T> &Bcc);
// template<typename T>
// void Bubble1(int tstp, herm_matrix_timestep<T> &C, herm_matrix<T> &A,herm_matrix<T> &B);
// ///
// template<typename T>
// void Bubble2(int tstp, herm_matrix_timestep_view<T> &C,int c1,int c2, herm_matrix_timestep_view<T> &A, herm_matrix_timestep_view<T> &Acc,int a1,int a2,
// herm_matrix_timestep_view<T> &B, herm_matrix_timestep_view<T> &Bcc,int b1,int b2);
// template<typename T>
// void Bubble2(int tstp, herm_matrix_timestep_view<T> &C,int c1,int c2, herm_matrix_timestep_view<T> &A,int a1,int a2,
// herm_matrix_timestep_view<T> &B,int b1,int b2);
// template<typename T>
// void Bubble2(int tstp, herm_matrix_timestep_view<T> &C, herm_matrix_timestep_view<T> &A, herm_matrix_timestep_view<T> &Acc,
// herm_matrix_timestep_view<T> &B, herm_matrix_timestep_view<T> &Bcc);
// template<typename T>
// void Bubble2(int tstp, herm_matrix_timestep_view<T> &C, herm_matrix_timestep_view<T> &A,herm_matrix_timestep_view<T> &B);
// ////
// template<typename T>
// void Bubble2(int tstp, herm_matrix<T> &C,int c1,int c2, herm_matrix<T> &A, herm_matrix<T> &Acc,int a1,int a2,
// herm_matrix<T> &B, herm_matrix<T> &Bcc,int b1,int b2);
// template<typename T>
// void Bubble2(int tstp, herm_matrix<T> &C,int c1,int c2, herm_matrix<T> &A,int a1,int a2,herm_matrix<T> &B,int b1,int b2);
// template<typename T>
// void Bubble2(int tstp, herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,herm_matrix<T> &B, herm_matrix<T> &Bcc);
// template<typename T>
// void Bubble2(int tstp, herm_matrix<T> &C, herm_matrix<T> &A,herm_matrix<T> &B);
// template<typename T>
// void Bubble2(int tstp, herm_matrix_timestep<T> &C,int c1,int c2, herm_matrix<T> &A, herm_matrix<T> &Acc,int a1,int a2,
// herm_matrix<T> &B, herm_matrix<T> &Bcc,int b1,int b2);
// template<typename T>
// void Bubble2(int tstp, herm_matrix_timestep<T> &C,int c1,int c2, herm_matrix<T> &A,int a1,int a2,herm_matrix<T> &B,int b1,int b2);
// template<typename T>
// void Bubble2(int tstp, herm_matrix_timestep<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,herm_matrix<T> &B, herm_matrix<T> &Bcc);
// template<typename T>
// void Bubble2(int tstp, herm_matrix_timestep<T> &C, herm_matrix<T> &A,herm_matrix<T> &B);
// #endif

} // namespace cntr

#endif  // CNTR_BUBBLE_DECL_H
