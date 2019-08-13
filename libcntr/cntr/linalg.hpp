/*#####################################################################
#
#  LINEAR ALGEBRA INTERFACE:
#
#  FOR A MOMENT, I STILL USE THOSE FROM THE EIGEN LIBRARY,
#  BUT THY MAY BE REPLACED BY DIRECT BLAS/LAPACK CALLS
#  AT A LATER TIME
#
#####################################################################*/

#ifndef _LINALG_H
#define _LINALG_H 1

#include <stdio.h>
#include <complex>
#include <cassert>
#include <cstring>

namespace linalg {

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// interface to simple linear algebra routines:
// implementation (using EIGEN library) defined in linalg_eigen.cpp
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//void cplx_sq_solve(void *a,void  *b,void *x,int n,int d);

/** \brief <b> Evaluate the eigen set of a given Hermitian matrix \f$A\f$. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* Evaluate the eigen system of a given Hermitian matrix \f$A\f$, and put the eigen values into "eval" and
* and the eigen vectors into "evec".
* <!-- ARGUMENTS
*      ========= -->
*
* @param size
* > Size of the matrix \f$A\f$.
* @param A
* > A Hermitian matrix given as a pointer.
* @param eval
* > A complex vector into which the eigen values are set. Given as a pointer.
* @param evec
* > A complex Matrix into which the eigen vectors are set. Given as a pointer.
*/

void eigen_hermv(int size,std::complex<double> *A,double *eval,std::complex<double> *evec);


/** \brief <b> Solve a linear equation ax=b. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* Solve a linear equation ax=b for a given real square matrx (\f$a\f$) and a real vector (\f$b\f$).
* <!-- ARGUMENTS
*      ========= -->
*
* @param ad
* > A real square matrix given as a pointer.
* @param bd
* > A real vector given as a pointer.
* @param xd
* > A real vector into which the solution is set. Given as a pointer.
* @param dim
* > Size of the matrix and the vector
*/
void real_sq_solve(double *ad,double *bd,double *xd,int dim);


/** \brief <b> Solve a linear equation ax=b. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* Solve a linear equation ax=b for a given complex square matrix \f$a\f$ and a complex vector \f$b\f$.
* <!-- ARGUMENTS
*      ========= -->
*
* @param a
* > A complex square matrix given as a pointer
* @param b
* > A complex vector given as a pointer
* @param x
* > A complex vector into which the solution is set.
* @param dim
* > Size of the matrix and the vector
*/
void cplx_sq_solve(void *a,void *b,void *x,int dim);

/** \brief <b> Solve a linear equation \f$ {\rm Diag}[a,a,a,.,a] x=b \f$. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* Solve a linear equation
* \f$ {\rm Diag}[a,a,a,.,a] x=b \f$ for a given complex square matrix \f$a\f$ and a complex vector \f$b\f$.
* The size of \f$a\f$ is \f$n\f$. Here \f$ {\rm Diag}[a,a,a,.,a]\f$ is the block diagonal square matrix, whose size is \f${\rm dim}\times d\f$.
* <!-- ARGUMENTS
*      ========= -->
*
* @param a
* > A complex square matrix, whose size is \f${\rm dim}\f$.
* @param b
* > A complex vector, whose size is \f$n\times d\f$
* @param x
* > A complex vector into which the solution is set.
* @param dim
* > Size of the matrix \f$a\f$
* @param d
* > \f${\rm dim}\times d\f$ is the size of \f$b\f$ and \f$x\f$.
*/

void cplx_sq_solve_many(void *a,void *b,void *x,int dim,int d);

//void cplx_tri_solve(double *a,double *b,double *x,int dim);

/** \brief <b> Evaluate the inverse matrix of a complex matrix \f$a\f$. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* Evaluate the inverse matrix of a given complex matrix \f$a\f$ and put it to \f$x\f$.
* <!-- ARGUMENTS
*      ========= -->
*
* @param a
* > A complex square matrix given as a pointer.
* @param x
* > A complex square matrix into which \f$a^{-1}\f$ is set. Given as a pointer.
* @param n
* > Size of the matrix \f$a\f$ and \f$x\f$.
*/

void cplx_matrix_inverse(void *a,void *x,int n);

/** \brief <b> Evaluate the inverse matrix of a real matrix \f$a\f$. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* Evaluate the inverse matrix of a given real matrix \f$a\f$ and put it to \f$x\f$.
* <!-- ARGUMENTS
*      ========= -->
*
* @param a
* > A real square matrix given as a pointer.
* @param x
* > A real square matrix into which \f$a^{-1}\f$ is set. Given as a pointer.
* @param n
* > Size of the matrix \f$a\f$ and \f$x\f$.
*/
void linalg_matrix_inverse(double *a,double *x,int n);
//void QR_decomposition(double *aa,double *qq,double *rr,int n,int m);
/////////////////////////////////////////////////////////////////


} // namespace



#endif  // linalg
