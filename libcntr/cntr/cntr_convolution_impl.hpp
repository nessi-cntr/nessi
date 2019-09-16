#ifndef CNTR_CONVOLUTION_IMPL_H
#define CNTR_CONVOLUTION_IMPL_H

#include "eigen_map.hpp"
#include "cntr_convolution_decl.hpp"
#include "cntr_elements.hpp"
#include "cntr_function_decl.hpp"
#include "cntr_herm_matrix_decl.hpp"

namespace cntr {

/* #######################################################################################
#
#  CONTOUR CONVOLUTION  C = A*B
#
#  (0) The routine computes C(t,t') = int_CC dt1 A(t,t1)B(t1,t') at timestep n,
#      or matsubara convolution.  i.e.,
#      C^ret(nh,t'<=nh), C^les(t<=nh,nh), C^tv(nt,tau=0..beta). If the timestep n==-1,
#      the Matsubara convolution is performed, and C^ret,C^tv, C^les are untouched
#
#  (1) A and B are contour Greenfunctions of type GA and GB, which may be scalar,
#      matrix etc. The behavior of the elements for the following operations must
#      be specified:
#      -  multiplication:  routine element_incr<T,GC,GA,GB>
#      -  conjugation: routine element_conjugate<T,G>
#
#  (2) Acc and Bcc are the conjugate functions to A and B. If A or B are hermitian,
#      just call convolution with A=Acc or B=Bcc, respectively.
#
#  (3) For the computation of timestep n, A(t,t') and B(t,t') are adressed at times
#      t,t' <= max(n,k), where k is the Integartion order (see Integrator).  I.e.,
#      the timesteps n=0..k can be computed only if A and B are given for t,t'<=k.
#
###########################################################################################*/

//   MATSUBARA INTEGRAL:
//
//   the integral   c = int_0^beta dx a(tau-x)b(x)  for matsubara_integral_1
//   the integral   c = int_0^beta dx a(x)b(x-tau)  for matsubara_integral_2
//
//   for beta=ntau, tau=x,
//   a and b are antiperiodic, i.e., a(-x+beta) =-a(-x) for 0 <= x <= beta
//
//   a,b,c are square-matrix-vauled, dimension size1, element size
//   sa=size1*size1
//   a is a pointer to an array of ntau elements of size sa*cplx:  a + i*sa ->
//   a(i)
//   b is a pointer to an array of ntau elements of size sb*cplx:  b + i*sa ->
//   b(i)
//   c is a pointer to an array of one element of size sa*cplx


/// @private
/** \brief <b> Calculates integral c = \f$\int_0^\beta dx a(\tau-x)b(x)\f$ </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate the integral:
 * > \f$ C = \int_0^\beta dx A(\tau-x)B(x)\f$
 * > The objects A,B and C are of type GG.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param size1
 * > [int] size of the matrix
 * @param m
 * > [int] number of the time step
 * @param ntau
 * > [int] nuber of the tau points
 * @param *C
 * > [std::complex] Pointer to the function, to which the results of calculation is given (this is just a constant matrix)
 * @param *A
 * > [std::complex] Pointer to a function A (this is just a constant matrix)
 * @param B
 * > [std::complex] Pointer to a function B (this is just a constant matrix)
 * @param I
* > [Integrator] integrator class
 * @param sig
 * > [int] Set `sig = -1` for fermions or `sig = +1` for bosons
 */
template <typename T, int SIZE1>
void matsubara_integral_1(int size1, int m, int ntau, std::complex<T> *C,
                          std::complex<T> *A, std::complex<T> *B,
                          integration::Integrator<T> &I, int sig) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1, j, l, sa1, sb1, sc1;
    cplx *amat, *bmat, *ctemp1, *ctemp2;
    T weight;

    sa1 = size1 * size1;
    sb1 = sa1;
    sc1 = sa1;
    // std::cout << "conv " << m << " " << ntau << std::endl;
    ctemp1 = new cplx[sc1];
    ctemp2 = new cplx[sc1];
    // CONTRIBUTION  FROM 0...TAU
    for (l = 0; l < sc1; l++) {
        ctemp1[l] = 0;
    }
    if (m >= k2 - 1) {
        amat = A + m * sa1;
        bmat = B;
        for (j = 0; j <= k; j++) {
            weight = I.gregory_omega(j);
            element_incr<T, SIZE1>(size1, ctemp1, weight, amat, bmat);
            amat -= sa1;
            bmat += sb1;
        }
        for (j = k1; j < m - k; j++) {
            element_incr<T, SIZE1>(size1, ctemp1, amat, bmat);
            amat -= sa1;
            bmat += sb1;
        }
        for (j = m - k; j <= m; j++) {
            weight = I.gregory_omega(m - j);
            element_incr<T, SIZE1>(size1, ctemp1, weight, amat, bmat);
            amat -= sa1;
            bmat += sb1;
        }
    } else if (m >= k) {
        amat = A + m * sa1;
        bmat = B;
        for (j = 0; j <= m; j++) {
            weight = I.gregory_weights(m, j);
            element_incr<T, SIZE1>(size1, ctemp1, weight, amat, bmat);
            amat -= sa1;
            bmat += sb1;
        }
    } else if (m > 0) { // here we need another strage boundary correction
        for (l = 0; l <= k; l++) {
            for (j = 0; j <= k; j++) {
                weight = I.rcorr(m, l, j);
                element_incr<T, SIZE1>(size1, ctemp1, weight, A + l * sa1,
                                       B + j * sb1);
            }
        }
    }
    // CONTRIBUTION  FROM TAU...BETA: mind the minus sign below
    for (l = 0; l < sc1; l++) {
        ctemp2[l] = 0;
    }
    if (ntau - m >= k2 - 1) { // usual gregory integration
        amat = A + sa1 * ntau;
        bmat = B + m * sb1;
        for (j = m; j <= m + k; j++) {
            weight = I.gregory_omega(j - m);
            element_incr<T, SIZE1>(size1, ctemp2, weight, amat, bmat);
            amat -= sa1;
            bmat += sb1;
        }
        for (j = m + k1; j < ntau - k; j++) {
            element_incr<T, SIZE1>(size1, ctemp2, amat, bmat);
            amat -= sa1;
            bmat += sb1;
        }
        for (j = ntau - k; j <= ntau; j++) {
            weight = I.gregory_omega(ntau - j);
            element_incr<T, SIZE1>(size1, ctemp2, weight, amat, bmat);
            amat -= sa1;
            bmat += sb1;
        }
    } else if (ntau - m >= k) {
        amat = A + sa1 * ntau;
        bmat = B + sb1 * m;
        for (j = m; j <= ntau; j++) {
            weight = I.gregory_weights(ntau - m, ntau - j);
            element_incr<T, SIZE1>(size1, ctemp2, weight, amat, bmat);
            amat -= sa1;
            bmat += sb1;
        }
    } else if (ntau - m >
               0) { // here we need another strange boundary correction
        // std::cout << "i am here " << std::endl;
        for (l = 0; l <= k; l++) {
            for (j = 0; j <= k; j++) {
                weight = I.rcorr(ntau - m, l, j);
                element_incr<T, SIZE1>(size1, ctemp2, weight,
                                       A + sa1 * (ntau - l),
                                       B + sb1 * (ntau - j));
            }
        }
    }
    // cmat += ctemp1 + sig * ctemp2:
    for (l = 0; l < sc1; l++)
        C[l] = ctemp1[l] + std::complex<T>(sig, 0.0) * ctemp2[l];
    delete[] ctemp1;
    delete[] ctemp2;
    return;
}
/// @private
/** \brief <b> Calculates contribution from 0 ... \f$\tau\f$ of matsubara integral_1 </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate only contribution from \f$0...\tau\f$ of 'matsubara_integral_1'
 * > ('sig' not needed here), i.e.:
 * > \f$ C = \int_0^\tau A(\tau-x)B(x)\f$
 * > The objects A,B and C are of type GG.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param size1
 * [int] size of the matrix
 * @param m
 * > [int] number of the time step
 * @param ntau
 * > [int] nuber of the tau points
 * @param *C
 * > [std::complex] Pointer to the function, to which the results of calculation is given (this is just a constant matrix)
 * @param *A
 * > [std::complex] Pointer to a function A (this is just a constant matrix)
 * @param B
 * > [std::complex] Pointer to a function B (this is just a constant matrix)
 * @param I
* > [Integrator] integrator class
 */
template <typename T, int SIZE1>
void matsubara_integral_1_1(int size1, int m, int ntau, std::complex<T> *C,
                            std::complex<T> *A, std::complex<T> *B,
                            integration::Integrator<T> &I) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1, j, l, sa1, sb1, sc1;
    cplx *amat, *bmat, *ctemp1;
    T weight;

    sa1 = size1 * size1;
    sb1 = sa1;
    sc1 = sa1;
    // std::cout << "conv " << m << " " << ntau << std::endl;
    ctemp1 = new cplx[sc1];
    // CONTRIBUTION  FROM 0...TAU
    for (l = 0; l < sc1; l++) {
        ctemp1[l] = 0;
    }
    if (m >= k2 - 1) {
        amat = A + m * sa1;
        bmat = B;
        for (j = 0; j <= k; j++) {
            weight = I.gregory_omega(j);
            element_incr<T, SIZE1>(size1, ctemp1, weight, amat, bmat);
            amat -= sa1;
            bmat += sb1;
        }
        for (j = k1; j < m - k; j++) {
            element_incr<T, SIZE1>(size1, ctemp1, amat, bmat);
            amat -= sa1;
            bmat += sb1;
        }
        for (j = m - k; j <= m; j++) {
            weight = I.gregory_omega(m - j);
            element_incr<T, SIZE1>(size1, ctemp1, weight, amat, bmat);
            amat -= sa1;
            bmat += sb1;
        }
    } else if (m >= k) {
        amat = A + m * sa1;
        bmat = B;
        for (j = 0; j <= m; j++) {
            weight = I.gregory_weights(m, j);
            element_incr<T, SIZE1>(size1, ctemp1, weight, amat, bmat);
            amat -= sa1;
            bmat += sb1;
        }
    } else if (m > 0) { // here we need another strage boundary correction
        for (l = 0; l <= k; l++) {
            for (j = 0; j <= k; j++) {
                weight = I.rcorr(m, l, j);
                element_incr<T, SIZE1>(size1, ctemp1, weight, A + l * sa1,
                                       B + j * sb1);
            }
        }
    }
    // cmat += ctemp1:
    for (l = 0; l < sc1; l++)
        C[l] = ctemp1[l];
    delete[] ctemp1;
    return;
}
/// @private
/** \brief <b> Calculates contribution from \f$\tau\f$ ... \f$\beta\f$ of matsubara integral_2 </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate only contribution from \f$\tau...\beta\f$ of 'matsubara_integral_2'
 * > ('sig' not needed here), i.e.:
 * > \f$ C = \int_\tau^\beta dx A(x)B(x-\tau)\f$
 * > The objects A,B and C are of type GG.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param size1
 * [int] size of the matrix
 * @param m
 * > [int] number of the time step
 * @param ntau
 * > [int] nuber of the tau points
 * @param *C
 * > [std::complex] Pointer to the function, to which the results of calculation is given (this is just a constant matrix)
 * @param *A
 * > [std::complex] Pointer to a function A (this is just a constant matrix)
 * @param B
 * > [std::complex] Pointer to a function B (this is just a constant matrix)
 * @param I
* > [Integrator] integrator class
 */
template <typename T, int SIZE1>
void matsubara_integral_2_2(int size1, int m, int ntau, std::complex<T> *C,
                            std::complex<T> *A, std::complex<T> *B,
                            integration::Integrator<T> &I) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1, j, l, sa1, sb1, sc1;
    cplx *amat, *bmat, *ctemp2;
    T weight;
    sa1 = size1 * size1;
    sb1 = sa1;
    sc1 = sa1;
    ctemp2 = new cplx[sc1];
    // CONTRIBUTION FROM TAU ... BETA
    for (l = 0; l < sc1; l++)
        ctemp2[l] = 0;
    if (m == ntau) { // donothing
    } else if (m > ntau - k) {
        for (l = 0; l <= k; l++) {
            for (j = 0; j <= k; j++) {
                weight = I.rcorr(ntau - m, l, j);
                element_incr<T, SIZE1>(size1, ctemp2, weight,
                                       A + sa1 * (ntau - l), B + sb1 * j);
            }
        }
    } else if (m > ntau - k2 + 1) {
        amat = A + sa1 * m;
        bmat = B;
        for (l = 0; l <= ntau - m; l++) {
            weight = I.gregory_weights(ntau - m, l);
            element_incr<T, SIZE1>(size1, ctemp2, weight, amat, bmat);
            amat += sa1;
            bmat += sb1;
        }
    } else {
        amat = A + sa1 * m;
        bmat = B;
        for (l = m; l <= m + k; l++) {
            weight = I.gregory_omega(l - m);
            element_incr<T, SIZE1>(size1, ctemp2, weight, amat, bmat);
            amat += sa1;
            bmat += sb1;
        }
        for (l = m + k1; l < ntau - k; l++) {
            element_incr<T, SIZE1>(size1, ctemp2, amat, bmat);
            amat += sa1;
            bmat += sb1;
        }
        for (l = ntau - k; l <= ntau; l++) {
            weight = I.gregory_omega(ntau - l);
            element_incr<T, SIZE1>(size1, ctemp2, weight, amat, bmat);
            amat += sa1;
            bmat += sb1;
        }
    }
    // ctv += ctemp2 :
    for (l = 0; l < sc1; l++)
        C[l] = ctemp2[l];
    delete[] ctemp2;
}
/// @private
/** \brief <b> Calculates integral \f$\int_0^\beta dx a(x)b(x-\tau)\f$ </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate the integral:
 * > \f$ C = \int_0^\beta dx A(x)B(x-\tau)\f$
 * > The objects A,B and C are of type GG.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param size1
 * [int] size of the matrix
 * @param m
 * > [int] number of the time step
 * @param ntau
 * > [int] nuber of the tau points
 * @param *C
 * > [std::complex] Pointer to the function, to which the results of calculation is given (this is just a constant matrix)
 * @param *A
 * > [std::complex] Pointer to a function A (this is just a constant matrix)
 * @param B
 * > [std::complex] Pointer to a function B (this is just a constant matrix)
 * @param I
* > [Integrator] integrator class
 * @param sig
 * > [int] Set `sig = -1` for fermions or `sig = +1` for bosons
 */
template <typename T, int SIZE1>
void matsubara_integral_2(int size1, int m, int ntau, std::complex<T> *C,
                          std::complex<T> *A, std::complex<T> *B,
                          integration::Integrator<T> &I, int sig) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1, j, l, sa1, sb1, sc1;
    cplx *amat, *bmat, *ctemp1, *ctemp2;
    T weight;

    sa1 = size1 * size1;
    sb1 = sa1;
    sc1 = sa1;

    ctemp1 = new cplx[sc1];
    ctemp2 = new cplx[sc1];
    // CONTRIBUTION FROM 0 ... TAU
    for (l = 0; l < sc1; l++)
        ctemp1[l] = 0;
    if (m == 0) {
        // donothing
    } else if (m < k) {
        for (j = 0; j <= k; j++) {
            for (l = 0; l <= k; l++) {
                weight = I.rcorr(m, l, j);
                element_incr<T, SIZE1>(size1, ctemp1, weight, A + j * sa1,
                                       B + sb1 * (ntau - l));
            }
        }
    } else if (m < k2 - 1) {
        amat = A + sa1 * m;
        bmat = B + sb1 * ntau;
        for (l = 0; l <= m; l++) {
            weight = I.gregory_weights(m, l);
            element_incr<T, SIZE1>(size1, ctemp1, weight, amat, bmat);
            amat -= sa1;
            bmat -= sb1;
        }
    } else {
        amat = A + sa1 * m;
        bmat = B + sb1 * ntau;
        for (l = 0; l <= k; l++) {
            weight = I.gregory_omega(l);
            element_incr<T, SIZE1>(size1, ctemp1, weight, amat, bmat);
            amat -= sa1;
            bmat -= sb1;
        }
        for (l = k1; l < m - k; l++) {
            element_incr<T, SIZE1>(size1, ctemp1, amat, bmat);
            amat -= sa1;
            bmat -= sb1;
        }
        for (l = m - k; l <= m; l++) {
            weight = I.gregory_omega(m - l);
            element_incr<T, SIZE1>(size1, ctemp1, weight, amat, bmat);
            amat -= sa1;
            bmat -= sb1;
        }
    }
    // CONTRIBUTION FROM TAU ... BETA
    for (l = 0; l < sc1; l++)
        ctemp2[l] = 0;
    if (m == ntau) { // donothing
    } else if (m > ntau - k) {
        for (l = 0; l <= k; l++) {
            for (j = 0; j <= k; j++) {
                weight = I.rcorr(ntau - m, l, j);
                element_incr<T, SIZE1>(size1, ctemp2, weight,
                                       A + sa1 * (ntau - l), B + sb1 * j);
            }
        }
    } else if (m > ntau - k2 + 1) {
        amat = A + sa1 * m;
        bmat = B;
        for (l = 0; l <= ntau - m; l++) {
            weight = I.gregory_weights(ntau - m, l);
            element_incr<T, SIZE1>(size1, ctemp2, weight, amat, bmat);
            amat += sa1;
            bmat += sb1;
        }
    } else {
        amat = A + sa1 * m;
        bmat = B;
        for (l = m; l <= m + k; l++) {
            weight = I.gregory_omega(l - m);
            element_incr<T, SIZE1>(size1, ctemp2, weight, amat, bmat);
            amat += sa1;
            bmat += sb1;
        }
        for (l = m + k1; l < ntau - k; l++) {
            element_incr<T, SIZE1>(size1, ctemp2, amat, bmat);
            amat += sa1;
            bmat += sb1;
        }
        for (l = ntau - k; l <= ntau; l++) {
            weight = I.gregory_omega(ntau - l);
            element_incr<T, SIZE1>(size1, ctemp2, weight, amat, bmat);
            amat += sa1;
            bmat += sb1;
        }
    }
    // ctv += ctemp2 + sig * ctemp1:
    for (l = 0; l < sc1; l++)
        C[l] = ctemp2[l] + std::complex<T>(sig, 0.0) * ctemp1[l];

    delete[] ctemp1;
    delete[] ctemp2;
}

/// @private
/** \brief <b> Performs the Matsubara convolution</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Performs Matsubara convolution \f$C^M=A^M*B^M\f$, i.e. we compute
* > \f$ C^M(\tau) = \int_0^\beta dx A^M(\tau-x) B^M(x)\f$
* > The objects A,B and C are of the general class 'GG'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param C
* > [GG] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [GG] contour Green's function
* @param B
* > [GG] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
*/
template <typename T, class GG, int SIZE1>
void convolution_matsubara_dispatch(GG &C, GG &A, GG &B,
                                    integration::Integrator<T> &I, T beta) {
    int ntau, l, m, size1 = C.size1();
    std::complex<T> *cmat;
    T dtau;
    ntau = A.ntau();
    for (m = 0; m <= ntau;
         m++) { // compute cmat(m*dtau) = int_0^beta dx amat(tau-x) b(x)
        matsubara_integral_1<T, SIZE1>(size1, m, ntau, C.matptr(m),
                                       A.matptr(0), B.matptr(0), I, A.sig());
    }
    // multiply by dtau:
    dtau = beta / ntau;
    cmat = C.matptr(0);
    m = (ntau + 1) * C.element_size();
    for (l = 0; l < m; l++)
        cmat[l] *= dtau;
    return;
}

/** \brief <b> Returns the result of the Matsubara convolution of two matrices. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calls Matsubara convolution routines to compute Matsubara convolution \f$C^M=A^M*B^M\f$.
* > \f$C^{R}\f$, \f$C^{\rceil}\f$, and \f$C^{<}\f$ are untouched.
* > The objects A,B and C are of the general class 'GG'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nomp
* > The number of threads.
* @param C
* > [GG] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [GG] contour Green's function
* @param B
* > [GG] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
*/
template <typename T, class GG>
void convolution_matsubara(GG &C, GG &A, GG &B, integration::Integrator<T> &I,
                           T beta) {
    int size1 = A.size1();
    assert(B.ntau() == A.ntau());
    assert(C.ntau() == A.ntau());
    assert(B.size1() == size1);
    assert(C.size1() == size1);
    if (size1 == 1)
      convolution_matsubara_dispatch<T, GG, 1>(C, A, B, I, beta);
    else
      convolution_matsubara_dispatch<T, GG, LARGESIZE>(C, A, B, I, beta);
}

#if CNTR_USE_OMP == 1

/// @private
/** \brief <b> Performs the Matsubara convolution. Uses OpenMP.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Performs Matsubara convolution \f$C^M=A^M*B^M\f$, i.e. we compute
* > \f$ C^M(\tau) = \int_0^\beta dx A^M(\tau-x) B^M(x)\f$
* > The objects A,B and C are of the general class 'GG'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nomp
* > The number of threads.
* @param C
* > [GG] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [GG] contour Green's function
* @param B
* > [GG] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
*/
template <typename T, class GG, int SIZE1>
void convolution_matsubara_omp_dispatch(int nomp, GG &C, GG &A, GG &B,
                                    integration::Integrator<T> &I, T beta) {
    int ntau, l, m, size1 = C.size1();
    std::complex<T> *cmat;
    T dtau;
    ntau = A.ntau();
#pragma omp parallel for num_threads(nomp)
    for (m = 0; m <= ntau;
         m++) { // compute cmat(m*dtau) = int_0^beta dx amat(tau-x) b(x)
        matsubara_integral_1<T, SIZE1>(size1, m, ntau, C.matptr(m),
                                       A.matptr(0), B.matptr(0), I, A.sig());
    }
    // multiply by dtau:
    dtau = beta / ntau;
    cmat = C.matptr(0);
    m = (ntau + 1) * C.element_size();
    for (l = 0; l < m; l++)
        cmat[l] *= dtau;
    return;
}

/** \brief <b> Returns the result of the Matsubara convolution of two matrices. Uses OpenMP.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calls Matsubara convolution routines to compute Matsubara convolution \f$C^M=A^M*B^M\f$.
* > \f$C^{R}\f$, \f$C^{\rceil}\f$, and \f$C^{<}\f$ are untouched.
* > The objects A,B and C are of the general class 'GG'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nomp
* > The number of threads.
* @param C
* > [GG] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [GG] contour Green's function
* @param B
* > [GG] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
*/
template <typename T, class GG>
void convolution_matsubara_nomp(int nomp, GG &C, GG &A, GG &B, integration::Integrator<T> &I,
                           T beta) {
    int size1 = A.size1();
    assert(B.ntau() == A.ntau());
    assert(C.ntau() == A.ntau());
    assert(B.size1() == size1);
    assert(C.size1() == size1);
    if (size1 == 1)
      convolution_matsubara_omp_dispatch<T, GG, 1>(nomp, C, A, B, I, beta);
    else
      convolution_matsubara_omp_dispatch<T, GG, LARGESIZE>(nomp, C, A, B, I, beta);
}

#endif // CNTR_USE_OMP

/** \brief <b> Retarded convolution at a given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate
 * > \f{eqnarray*}{
 *     C &=& A*B\\
 *     C^{R}(t,t') &=& \int_{t'}^t d\bar{t} A^{R}(t,\bar{t}) B^{R}(\bar{t},t')\f}
 * >
 * > At time-step `t = n h` for all \f$t'\f$.
 * > The result is written into `C`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * [int] given time-step
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_ret(int n, GG &C, GG &A, GG &Acc, GG &B, GG &Bcc,
                              integration::Integrator<T> &I, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k();
    int sa, sb, sc, j, m, j1, n1, l, size1 = C.size1();
    cplx *aret, *cret, *bret, *btemp, *atemp, *result;
    T weight;

    // duplicated arguments
    sa = A.element_size();
    sb = B.element_size();
    sc = C.element_size();
    n1 = (n < k ? k : n);
    atemp = new cplx[sa];
    btemp = new cplx[sb];
    // such that G=Sigma*G can be called without creating a mess,
    // data are first written in a temporary variable and then written to C at
    // the end
    result = new cplx[(n + 1) * sc];
    for (l = 0; l < (n + 1) * sc; l++)
        result[l] = 0;
    // check consistency:
    assert(sa * sb * sc != 0);
    assert(Acc.element_size() == sa);
    assert(Bcc.element_size() == sb);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);
    assert(C.nt() >= n);

    if (n >= k) {
        // CONTRIBUTION FROM BRET: loop over lines of Bret
        for (m = 0; m <= n; m++) { // contribution to integral from Bret(m,j)
            aret = A.retptr(n, m);
            for (l = 0; l < sa; l++)
                atemp[l] = aret[l] * h; // here enters h
            // the triangle j <= m
            bret = B.retptr(m, 0);
            cret = result;
            // in the following sector the weights are 1
            if (m < n - k) {
                for (j = 0; j < m - k; j++) {
                    element_incr<T, SIZE1>(size1, cret, atemp,
                                           bret); // cret += aret*bret
                    bret += sb;
                    cret += sc;
                }
            } else { // m>=n-k
                weight = I.gregory_omega(n - m);
                for (j = 0; j < m - k; j++) {
                    element_incr<T, SIZE1>(size1, cret, weight, atemp,
                                           bret); // cret += aret*bret
                    bret += sb;
                    cret += sc;
                }
            }
            // contribution from the stripe m-j <= k, with different weight
            j1 = m - k;
            if (j1 < 0)
                j1 = 0;
            for (j = j1; j <= m; j++) {
                weight = I.gregory_weights(n - j, n - m);
                element_incr<T, SIZE1>(size1, cret, weight, atemp,
                                       bret); // cret += aret*bret
                bret += sb;
                cret += sc;
            }
        }
        // CONTRIBUTION FROM BRET^CONJ:
        for (m = n - k; m < n; m++) {
            aret = A.retptr(n, m);
            for (l = 0; l < sa; l++)
                atemp[l] = aret[l] * h; // here enters h
            for (j = m + 1; j <= n; j++) {
                weight = I.gregory_weights(n - j, n - m);
                element_conj<T, SIZE1>(size1, btemp, Bcc.retptr(j, m));
                // mind the minus sign: B(m,j) continued to -Bcc(j,m)*
                element_incr<T, SIZE1>(size1, result + j * sc, -weight, atemp,
                                       btemp);
            }
        }
    } else { // n < k
        for (j = 0; j <= n; j++) {
            cret = result + j * sc;
            for (m = 0; m <= k; m++) {
                weight = I.poly_integration(j, n, m) * h;
                if (m >= j) {
                    for (l = 0; l < sb; l++)
                        btemp[l] = B.retptr(m, j)[l];
                } else {
                    element_conj<T, SIZE1>(size1, btemp, Bcc.retptr(j, m));
                    weight *= -1;
                }
                if (n >= m) {
                    for (l = 0; l < sa; l++)
                        atemp[l] = A.retptr(n, m)[l];
                } else {
                    element_conj<T, SIZE1>(size1, atemp, Acc.retptr(m, n));
                    weight *= -1;
                }
                // std::cout << n << " " << m << " " << j << std::endl;
                element_incr<T, SIZE1>(size1, cret, weight, atemp, btemp);
            }
        }
    }

    cret = C.retptr(n, 0);
    for (l = 0; l < (n + 1) * sc; l++)
        cret[l] = result[l];
    delete[] result;
    delete[] atemp;
    delete[] btemp;
    return;
}

/** \brief <b> Left-Mixing convolution at given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate
 * > \f{eqnarray*}{
 *     C &=& A*B\\
 *     C^{\rceil}(t,\tau') &=& \int_0^\beta ds  A^{\rceil}(t,s) B^M(s - \tau')
 *                       + \int_0^t ds A^R(t,s) B^{\rceil}(s,\tau') \f}
 * >
 * > At time-step `t = n h` for all \f$\tau'\f$.
 * > The result is written into `ctv`.
 * >
 * > `ctv` is a pointer to the timestep of elements of type GG,
 * > so it must have size `(C.ntau_ + 1) * C.element_size_`
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * [int] given time-step
 * @param ctv
 * > [herm_matrix] left-mixing component \f$C^{\rceil}\f$ to which the result of the convolution is given.
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inverse temperature
 * @param h
 * > time step
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_tv(int n, std::complex<T> *ctv, GG &C, GG &A,
                             GG &Acc, GG &B, GG &Bcc,
                             integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;

    int k = I.get_k(); // order
    int sa = A.element_size();
    int sb = B.element_size();
    int sc = C.element_size();
    int ntau = A.ntau();
    int n1 = (n > k ? n : k);
    int size1 = C.size1();
    T dtau = beta / ntau;
    T weight;
    cplx *ctemp1, *ctv1, *btv, *atemp;
    int j, m, l;

    // check consistency:
    assert(sa * sb * sc != 0);
    assert(Acc.element_size() == sa);
    assert(Bcc.element_size() == sb);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(C.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);

    // CONTRIBUTION FROM Atv * Bmat:
    //     very similar to computing the matsubara convolution
    ctemp1 = new cplx[sc];
    for (m = 0; m <= ntau; m++) {
        matsubara_integral_2<T, SIZE1>(size1, m, ntau, ctemp1, A.tvptr(n, 0),
                                       B.matptr(0), I, B.sig());
        for (l = 0; l < sc; l++)
            ctv[m * sc + l] = dtau * ctemp1[l];
    }
    delete[] ctemp1;

    // CONTRIBUTION FROM Aret * Btv:
    //     loop over lines j
    atemp = new cplx[sa];
    for (j = 0; j <= n1; j++) { // j <= n1,  n1 = max(n, k)
        weight = I.gregory_weights(n, j);
        if (n < j) { // j > n  ==>  n < k ;  j <= max(n, k)
            element_conj<T, SIZE1>(size1, atemp, Acc.retptr(j, n));
            element_smul<T, SIZE1>(size1, atemp, -1);
            // atemp = dt (-Acc(j,n)*)
        } else { // j <= n
            element_set<T, SIZE1>(size1, atemp, A.retptr(n, j));
            // atemp = dt A(n,j)
        }
        element_smul<T, SIZE1>(size1, atemp, h); // here enters h
        btv = B.tvptr(j, 0);
        ctv1 = ctv;
        // ctv1 = ctv1 + weight atemp . btv
        if (weight != 1) {
            for (m = 0; m <= ntau; m++) {
                element_incr<T, SIZE1>(size1, ctv1, weight, atemp, btv);
                btv += sb;
                ctv1 += sc;
            }
        } else {
            for (m = 0; m <= ntau; m++) {
                element_incr<T, SIZE1>(size1, ctv1, atemp, btv);
                btv += sb;
                ctv1 += sc;
            }
        }
    }
    delete[] atemp;
}

/** \brief <b> Left-Mixing convolution at given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate
 * > \f{eqnarray*}{
 *     C &=& A*B\\
 *     C^{\rceil}(t,\tau') &=& \int_0^\beta ds  A^{\rceil}(t,s) B^M(s - \tau')
 *                       + \int_0^t ds A^R(t,s) B^{\rceil}(s,\tau') \f}
 * >
 * > At time-step `t = n h` for all \f$\tau'\f$.
 * > The result is written at a given time step into 'C'.
 * > The objects A,B and C are of type GG.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * [int] given time-step
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inverse temperature
 * @param h
 * > time step
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_tv(int n, GG &C, GG &A, GG &Acc, GG &B, GG &Bcc,
                             integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int ntau, m, sc, size1 = C.size1();
    cplx *ctv;
    ntau = C.ntau();
    sc = C.element_size();

    assert(sc > 0 && ntau > 0);
    ctv = new cplx[(ntau + 1) * sc];
    convolution_timestep_tv<T, GG, SIZE1>(n, ctv, C, A, Acc, B, Bcc, I, beta,
                                          h);
    for (m = 0; m <= ntau; m++)
        element_set<T, SIZE1>(size1, C.tvptr(n, m), ctv + m * sc);
    delete[] ctv;
}

/** \brief <b> Calculation of \f$C = A^{\rceil}*B^{\lceil}\f$ at a given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate \f$C = A^{\rceil}*B^{\lceil}\f$
 * > at a given time-step `t = n h`.
 * >
 * > Note: 'cles' is a pointer to the timestep of elements of type GG.
 * > It must have size   (max(n,k)+1) * C.element_size_. The result is added to 'cles'.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * [int] given time-step
 * @param j1
 * [int] lower integration limit
 * @param j2
 * [int] upper integration limit
 * @param cles
 * > [std::complex] part of the lesser component \f$C^{<}\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void
convolution_timestep_les_tvvt(int n, int j1, int j2, std::complex<T> *cles,
                              GG &C, GG &A, GG &Acc, GG &B, GG &Bcc,
                              integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1;
    int sa, sb, sc, ntau, j, m, l, n1, sig, size1 = C.size1();
    T weight, dtau;
    cplx *atv, *btv, *btemp, *cles1, idtau;

    ntau = A.ntau();
    sa = A.element_size();
    sb = B.element_size();
    sc = C.element_size();
    dtau = beta / ntau;
    sig = A.sig();
    n1 = (n < k ? k : n);
    // check consistency:
    assert(sa * sb * sc != 0);
    assert(ntau >= k);
    assert(Acc.element_size() == sa);
    assert(Bcc.element_size() == sb);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(C.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);
    assert(sig == Acc.sig());
    assert(sig == B.sig());
    assert(sig == Bcc.sig());
    assert(0 <= j1 && j1 <= n1);
    assert(j1 <= j2 && j2 <= n1);

    // contribution from Atv*Bvt = Atv(jh,tau) * Bcc^tv(nh,beta-tau)^* *
    // (-Bose/Fermi)
    btemp = new cplx[(ntau + 1) * sb];
    idtau = cplx(0, -dtau);
    for (m = 0; m <= ntau; m++)
        element_conj<T, SIZE1>(size1, btemp + m * sb, Bcc.tvptr(n, ntau - m));
    for (l = 0; l < (ntau + 1) * sb; l++)
        btemp[l] *= idtau * (-(T)sig);
    for (j = j1; j <= j2; j++) {
        btv = btemp;
        atv = A.tvptr(j, 0);
        cles1 = cles + j * sc;
        if (ntau < k2 - 1) {
            for (m = 0; m <= ntau; m++) {
                weight = I.gregory_weights(ntau, m);
                element_incr<T, SIZE1>(size1, cles1, weight, atv, btv);
                atv += sa;
                btv += sb;
            }
        } else {
            for (m = 0; m <= k; m++) {
                weight = I.gregory_omega(m);
                element_incr<T, SIZE1>(size1, cles1, weight, atv, btv);
                atv += sa;
                btv += sb;
            }
            for (m = k1; m < ntau - k; m++) {
                element_incr<T, SIZE1>(size1, cles1, atv, btv);
                atv += sa;
                btv += sb;
            }
            for (m = ntau - k; m <= ntau; m++) {
                weight = I.gregory_omega(ntau - m);
                element_incr<T, SIZE1>(size1, cles1, weight, atv, btv);
                atv += sa;
                btv += sb;
            }
        }
    }
    delete[] btemp;
    return;
}
/** \brief <b> Calculation of \f$C = A^{\rceil}*B^{\lceil}\f$ at a given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Calls a routine to calculate \f$C = A^{\rceil}*B^{\lceil}\f$
 * > at a given time-step `t = n h`.
 * >
 * > Note: 'cles' is a pointer to the timestep of elements of type GG.
 * > It must have size   (max(n,k)+1) * C.element_size_. The result is added to 'cles'.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * [int] given time-step
 * @param cles
 * > [std::complex] part of the lesser component \f$C^{<}\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_tvvt(int n, std::complex<T> *cles, GG &C, GG &A,
                                   GG &Acc, GG &B, GG &Bcc,
                                   integration::Integrator<T> &I, T beta,
                                   T h) {
    int k = I.get_k();
    int n1 = (n < k ? k : n);
    return convolution_timestep_les_tvvt<T, GG, SIZE1>(
        n, 0, n1, cles, C, A, Acc, B, Bcc, I, beta, h);
}
/** \brief <b> Calculation of \f$C = A^{<}*B^{A}\f$ at a given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate \f$C = A^{<}*B^{A}\f$
 * > at a given time-step `t = n h`.
 * >
 * > Note: 'cles' is a pointer to the timestep of elements of type GG.
 * > It must have size   (max(n,k)+1) * C.element_size_. The result is added to 'cles'.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * > [int] given time-step
 * @param j1
 * > [int] lower integration limit
 * @param j2
 * > [int] upper integration limit
 * @param cles
 * > [std::complex] part of the lesser component \f$C^{<}\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void
convolution_timestep_les_lesadv(int n, int j1, int j2, std::complex<T> *cles,
                                GG &C, GG &A, GG &Acc, GG &B, GG &Bcc,
                                integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k();
    int sa, sb, sc, ntau, j, m, l, n1, size1 = C.size1();
    T weight, dtau;
    cplx idtau, *ales, *badv;

    ntau = A.ntau();
    sa = A.element_size();
    sb = B.element_size();
    sc = C.element_size();
    dtau = beta / ntau;
    n1 = (n < k ? k : n); // j1 = 0, j2 = n1.
    // check consistency:
    assert(sa * sb * sc != 0);
    assert(ntau >= k);
    assert(Acc.element_size() == sa);
    assert(Bcc.element_size() == sb);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(C.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);
    assert(0 <= j1 && j1 <= n1);
    assert(j1 <= j2 && j2 <= n1);

    // contribution from Ales(j,m)*Badv(m,n)
    ales = new cplx[sa];
    badv = new cplx[(n1 + 1) * sb];
    for (m = 0; m <= n1; m++) {
        weight = h * I.gregory_weights(n, m);
        if (m <= n) {
            element_conj<T, SIZE1>(size1, badv + m * sb, Bcc.retptr(n, m));
            for (l = 0; l < sb; l++)
                badv[m * sb + l] *= weight;
        } else {
            for (l = 0; l < sb; l++)
                badv[m * sb + l] = -weight * B.retptr(m, n)[l];
        }
    }
    // Cles(t',t) += \int_0^{t'} ds -Ales*(s,t') * Bret*(t,s)
    for (j = j1; j <= j2; j++) {
        for (m = 0; m < j; m++) { // inner loop over cache-optimal `s`.
            element_minusconj<T, SIZE1>(size1, ales, Acc.lesptr(m, j));
            element_incr<T, SIZE1>(size1, cles + j * sc, ales, badv + m * sb);
        }
    }
    // Cles(t',t) += \int_{t'}^t ds Ales(t',s) * Bret*(t,s)
    for (m = 0; m <= n1; ++m) {
        int jmax = std::min(j2, m);
        for (j = j1; j <= jmax; ++j) { // inner loop over cache-optimal `t'`.
            element_incr<T, SIZE1>(size1, cles + j * sc, A.lesptr(j, m),
                                   badv + m * sb);
        }
    }
    delete[] ales;
    delete[] badv;
    return;
}
/** \brief <b> Calculation of \f$C = A^{<}*B^{A}\f$ at a given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Calls a routine to calculate \f$C = A^{<}*B^{A}\f$
 * > at a given time-step `t = n h`.
 * >
 * > Note: 'cles' is a pointer to the timestep of elements of type GG.
 * > It must have size   (max(n,k)+1) * C.element_size_. The result is added to 'cles'.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * > [int] given time-step
 * @param cles
 * > [std::complex] part of the lesser component \f$C^{<}\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_lesadv(int n, std::complex<T> *cles, GG &C,
                                     GG &A, GG &Acc, GG &B, GG &Bcc,
                                     integration::Integrator<T> &I, T beta,
                                     T h) {
    int k = I.get_k();
    int n1 = (n < k ? k : n);
    return convolution_timestep_les_lesadv<T, GG, SIZE1>(
        n, 0, n1, cles, C, A, Acc, B, Bcc, I, beta, h);
}
/** \brief <b> Calculation of \f$C = A^{R}*B^{<}\f$ at a given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here, we calculate \f$C = A^{R}*B^{<}\f$
 * > at a given time-step `t = n h`.
 * >
 * > Note: 'cles' is a pointer to the timestep of elements of type GG.
 * > It must have size   (max(n,k)+1) * C.element_size_. The result is added to 'cles'.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * [int] given time-step
 * @param cles
 * > [std::complex] part of the lesser component \f$C^{<}\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_retles(int n, std::complex<T> *cles, GG &C,
                                     GG &A, GG &Acc, GG &B, GG &Bcc,
                                     integration::Integrator<T> &I, T beta,
                                     T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1, size1 = C.size1();
    int sa, sb, sc, ntau, j, m, l, n1;
    T weight, dtau;
    cplx *aret, *atemp, *btemp, *cles1, *bles, idtau;

    ntau = A.ntau();
    sa = A.element_size();
    sb = B.element_size();
    sc = C.element_size();
    dtau = beta / ntau;
    n1 = (n < k ? k : n);
    // check consistency:
    assert(sa * sb * sc != 0);
    assert(ntau >= k);
    assert(Acc.element_size() == sa);
    assert(Bcc.element_size() == sb);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(C.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);

    // contribution from Aret*Bles
    btemp = new cplx[(n1 + 1) * sb];
    atemp = new cplx[sa];
    for (m = 0; m <= n1; m++) { // btemp(m) --> B^<(m,n)
        if (m <= n) {
            for (l = 0; l < sb; l++)
                btemp[m * sb + l] = h * B.lesptr(m, n)[l];
        } else {
            element_conj<T, SIZE1>(size1, btemp + m * sb, Bcc.lesptr(n, m));
            for (l = 0; l < sb; l++)
                btemp[m * sb + l] *= -h;
        }
    }
    for (j = 0; j <= n; j++) { // compute -> cles1(j,n)
        cles1 = cles + j * sc;
        // CONTRINBUTION  FROM A_RET
        if (j >= k2 - 1) {
            aret = A.retptr(j, 0);
            bles = btemp;
            for (m = 0; m <= k; m++) {
                weight = I.gregory_omega(m);
                element_incr<T, SIZE1>(size1, cles1, weight, aret, bles);
                bles += sb;
                aret += sa;
            }
            for (m = k1; m < j - k; m++) {
                element_incr<T, SIZE1>(size1, cles1, aret, bles);
                bles += sb;
                aret += sa;
            }
            for (m = j - k; m <= j; m++) {
                weight = I.gregory_omega(j - m);
                element_incr<T, SIZE1>(size1, cles1, weight, aret, bles);
                bles += sb;
                aret += sa;
            }
        } else {
            aret = A.retptr(j, 0);
            bles = btemp;
            for (m = 0; m <= j; m++) {
                weight = I.gregory_weights(j, m);
                element_incr<T, SIZE1>(size1, cles1, weight, aret, bles);
                bles += sb;
                aret += sa;
            }
        }
        // CONTRIBUTION FROM ACC_RET
        if (j < k) {
            for (m = j + 1; m <= k; m++) {
                element_conj<T, SIZE1>(size1, atemp, Acc.retptr(m, j));
                weight = -I.gregory_weights(j, m);
                bles = btemp + m * sb;
                element_incr<T, SIZE1>(size1, cles1, weight, atemp, bles);
            }
        }
    }
    delete[] btemp;
    delete[] atemp;
    return;
}
/** \brief <b> Lesser convolution at a given time-step of two matrices. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate
 * >   \f{eqnarray*}{
 *     C^{<} &=& A^{\rceil}*B^{\lceil} + A^<*B^{A} + A^R*B^<\\
 *     C^{<}(t,t') &=& -i \int_0^\beta ds  A^{\rceil}(t,s) B^{\lceil}(s,t')+ \int_0^t ds  A^{<}(t,s) B^A(s,t')
 *                       + \int_0^t ds  A^{R}(t,s) B^<(s,t') \f}
 * > at a given time-step `t = n h`.
 * > The result is written into a contor Green's function `C`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * > [int] given time-step
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_les(int n, GG &C, GG &A, GG &Acc, GG &B, GG &Bcc,
                              integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    cplx *cles;
    int m, sc, n1, k = I.get_k(), size1 = C.size1();
    sc = C.element_size();
    assert(sc > 0);
    n1 = (k > n ? k : n);
    cles = new cplx[(n1 + 1) * sc];
    for (m = 0; m < sc * (n1 + 1); m++)
        cles[m] = 0;
    convolution_timestep_les_tvvt<T, GG, SIZE1>(n, cles, C, A, Acc, B, Bcc, I,
                                                beta, h);
    convolution_timestep_les_lesadv<T, GG, SIZE1>(n, cles, C, A, Acc, B, Bcc,
                                                  I, beta, h);
    convolution_timestep_les_retles<T, GG, SIZE1>(n, cles, C, A, Acc, B, Bcc,
                                                  I, beta, h);
    for (m = 0; m <= n; m++)
        element_set<T, SIZE1>(size1, C.lesptr(m, n), cles + m * sc);
    delete[] cles;
    return;
}
/** \brief <b> Returns convolution of two matrices at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*B of the objects with class 'herm_matrix'
* > at a given time step 't=nh'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] complex conjugate to A
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] complex conjugate to B
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A,
                          herm_matrix<T> &Acc, herm_matrix<T> &B,
                          herm_matrix<T> &Bcc, integration::Integrator<T> &I,
                          T beta, T h) {
    int size1 = C.size1(), ntau = C.ntau(), k = I.k(), n1 = (n < k ? k : n);
    if (n == -1) {
        convolution_matsubara(C, A, B, I, beta);
        return;
    }
    assert(n >= 0);
    assert(A.size1() == size1);
    assert(Acc.size1() == size1);
    assert(B.size1() == size1);
    assert(Bcc.size1() == size1);
    assert(A.ntau() == ntau);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);
    assert(C.nt() >= n);
    if (size1 == 1) {
        convolution_timestep_ret<T, herm_matrix<T>, 1>(n, C, A, Acc, B, Bcc,
                                                       I, h);
        convolution_timestep_tv<T, herm_matrix<T>, 1>(n, C, A, Acc, B, Bcc, I,
                                                      beta, h);
        convolution_timestep_les<T, herm_matrix<T>, 1>(n, C, A, Acc, B, Bcc,
                                                       I, beta, h);
    } else {
        convolution_timestep_ret<T, herm_matrix<T>, LARGESIZE>(n, C, A, Acc,
                                                               B, Bcc, I, h);
        convolution_timestep_tv<T, herm_matrix<T>, LARGESIZE>(
            n, C, A, Acc, B, Bcc, I, beta, h);
        convolution_timestep_les<T, herm_matrix<T>, LARGESIZE>(
            n, C, A, Acc, B, Bcc, I, beta, h);
    }
}


/** \brief <b> Returns convolution of two matrices at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*B of the objects with class 'herm_matrix'
* > at a given time step 't=nh'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] complex conjugate to A
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] complex conjugate to B
* @param beta
* > inversed temperature
* @param h
* > time step interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A,
                          herm_matrix<T> &Acc, herm_matrix<T> &B,
                          herm_matrix<T> &Bcc,
                          T beta, T h, int SolveOrder) {
    int size1 = C.size1(), ntau = C.ntau(), n1 = (n < SolveOrder ? SolveOrder : n);
    if (n == -1) {
        convolution_matsubara(C, A, B, integration::I<T>(SolveOrder), beta);
        return;
    }
    assert(n >= 0);
    assert(A.size1() == size1);
    assert(Acc.size1() == size1);
    assert(B.size1() == size1);
    assert(Bcc.size1() == size1);
    assert(A.ntau() == ntau);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);
    assert(C.nt() >= n);
    if (size1 == 1) {
        convolution_timestep_ret<T, herm_matrix<T>, 1>(n, C, A, Acc, B, Bcc,
                                                       integration::I<T>(SolveOrder), h);
        convolution_timestep_tv<T, herm_matrix<T>, 1>(n, C, A, Acc, B, Bcc, integration::I<T>(SolveOrder),
                                                      beta, h);
        convolution_timestep_les<T, herm_matrix<T>, 1>(n, C, A, Acc, B, Bcc,
                                                       integration::I<T>(SolveOrder), beta, h);
    } else {
        convolution_timestep_ret<T, herm_matrix<T>, LARGESIZE>(n, C, A, Acc,
                                                               B, Bcc, integration::I<T>(SolveOrder), h);
        convolution_timestep_tv<T, herm_matrix<T>, LARGESIZE>(
            n, C, A, Acc, B, Bcc, integration::I<T>(SolveOrder), beta, h);
        convolution_timestep_les<T, herm_matrix<T>, LARGESIZE>(
            n, C, A, Acc, B, Bcc, integration::I<T>(SolveOrder), beta, h);
    }
}


/** \brief <b> Returns convolution of two hermitian matrices at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calls convolution routine to compute contour convolution C=A*B of two hermitian matrices
* > at a given time step 't=nh'.
* > The objects A,B and C are of the class 'herm_matrix'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function (A=Acc)
* @param B
* > [herm_matrix] contour Green's function (B=Bcc)
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A,
                          herm_matrix<T> &B, integration::Integrator<T> &I,
                          T beta, T h) {
    convolution_timestep<T>(n, C, A, A, B, B, I, beta, h);
}


/** \brief <b> Returns convolution of two hermitian matrices at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calls convolution routine to compute contour convolution C=A*B of two hermitian matrices
* > at a given time step 't=nh'.
* > The objects A,B and C are of the class 'herm_matrix'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function (A=Acc)
* @param B
* > [herm_matrix] contour Green's function (B=Bcc)
* @param beta
* > inversed temperature
* @param h
* > time step interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A,
                          herm_matrix<T> &B,
                          T beta, T h, int SolveOrder) {
    convolution_timestep<T>(n, C, A, A, B, B, beta, h, SolveOrder);
}

/** \brief <b> Returns the result of the contour convolution of two 'herm_matrix' objects</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calls convolution routines to compute contour convolution C=A*B.
* > The objects A,B and C are of the class 'herm_matrix'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param C
* > [herm_matrix] Matrix to which the result of the convolution is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] complex conjugate to A
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] complex conjugate to B
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                 herm_matrix<T> &B, herm_matrix<T> &Bcc,
                 integration::Integrator<T> &I, T beta, T h) {
    int tstp;
    convolution_matsubara(C, A, B, I, beta);
    for (tstp = 0; tstp <= C.nt(); tstp++)
        convolution_timestep<T>(tstp, C, A, Acc, B, Bcc, I, beta, h);
}

/// @private
/** \brief <b> Performs the Matsubara convolution of two matrices and a function, i.e. C=A*fxB</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Performs Matsubara convolution \f$C^M=A^M*fxB^M\f$, otherwise like contour convolution
* > (needs some functions defined in cntr_convolution.hpp). Here, we compute
* > \f$ C^M(\tau) = \int_0^\beta dx A^M(\tau-x) F(x)B^M(x)\f$
* > The objects A,B and C are of the general class 'GG'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param C
* > [GG] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [GG] contour Green's function
* @param *f0
* > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
* @param B
* > [GG] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
*/
template <typename T, class GG, int SIZE1>
void convolution_matsubara_dispatch(GG &C, GG &A, std::complex<T> *f0, GG &B,
                                    integration::Integrator<T> &I, T beta) {
    int ntau, l, m, size1 = C.size1(), sb = B.element_size();
    std::complex<T> *cmat, *bmat;
    T dtau;
    ntau = A.ntau();
    bmat = new std::complex<T>[sb * (ntau + 1)];
    for (m = 0; m <= ntau; m++)
        element_mult<T, SIZE1>(size1, bmat + m * sb, f0, B.matptr(m));
    for (m = 0; m <= ntau; m++) { // compute cmat(m*dtau) = int_0^beta dx amat(tau-x) b(x)
        matsubara_integral_1<T, SIZE1>(size1, m, ntau, C.matptr(m), A.matptr(0), bmat, I,
                                       A.sig());
    }
    delete[] bmat;
    // multiply by dtau:
    dtau = beta / ntau;
    cmat = C.matptr(0);
    m = (ntau + 1) * C.element_size();
    for (l = 0; l < m; l++)
        cmat[l] *= dtau;
    return;
}
/** \brief <b> Returns the result of the Matsubara convolution of two matrices and a function</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calls Matsubara convolution routines to compute Matsubara convolution C=A*fxB.
* > The objects A,B and C are of the general class 'GG'.
* > 'f0' is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param C
* > [GG] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [GG] contour Green's function
* @param *f0
* > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
* @param B
* > [GG] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
*/
template <typename T, class GG>
void convolution_matsubara(GG &C, GG &A, std::complex<T> *f0, GG &B,
                           integration::Integrator<T> &I, T beta) {
    int size1 = A.size1();
    assert(B.ntau() == A.ntau());
    assert(C.ntau() == A.ntau());
    assert(B.size1() == size1);
    assert(C.size1() == size1);
    if (size1 == 1)
        convolution_matsubara_dispatch<T, GG, 1>(C, A, f0, B, I, beta);
    else
        convolution_matsubara_dispatch<T, GG, LARGESIZE>(C, A, f0, B, I, beta);
}
/** \brief <b> Retarded convolution of two matrices and a function at given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate
 * > \f{eqnarray*}{
 *     C &=& A*ft*B\\
 *     C^{R}(t,t') &=&  \int_{t'}^t d\bar{t} A^{R}(t,\bar{t}) ft B^{R}(\bar{t},t')\f}
 * >
 * > At time-step `t = n h` for all \f$\tau'\f$. 'ft' is a pointer to \f$F(t)\f$ on the real axis.
 * > The result is written into `C`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * [int] given time-step
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_ret(int n, GG &C, GG &A, GG &Acc, std::complex<T> *ft, GG &B,
                              GG &Bcc, integration::Integrator<T> &I, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k();
    int sf, sa, sb, sc, j, m, j1, n1, l, size1 = C.size1();
    cplx *aret, *cret, *bret, *btemp, *atemp, *result;
    T weight;

    // duplicated arguments
    sa = A.element_size();
    sb = B.element_size();
    sc = C.element_size();
    sf = size1 * size1;
    n1 = (n < k ? k : n);
    atemp = new cplx[sa];
    aret = new cplx[sa];
    btemp = new cplx[sb];
    // such that G=Sigma*G can be called without creating a mess,
    // data are first written in a temporary variable and then written to C at the end
    result = new cplx[(n + 1) * sc];
    for (l = 0; l < (n + 1) * sc; l++)
        result[l] = 0;
    // check consistency:
    assert(sa * sb * sc != 0);
    assert(Acc.element_size() == sa);
    assert(Bcc.element_size() == sb);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);
    assert(C.nt() >= n);

    if (n >= k) {
        // CONTRIBUTION FROM BRET: loop over lines of Bret
        for (m = 0; m <= n; m++) { // contribution to integral from Bret(m,j)
            element_mult<T, SIZE1>(size1, aret, A.retptr(n, m), ft + m * sf);
            // aret = A.retptr(n,m);  // without f
            for (l = 0; l < sa; l++)
                atemp[l] = aret[l] * h; // here enters h
            // the triangle j <= m
            bret = B.retptr(m, 0);
            cret = result;
            // in the following sector the weights are 1
            if (m < n - k) {
                for (j = 0; j < m - k; j++) {
                    element_incr<T, SIZE1>(size1, cret, atemp, bret); // cret += aret*bret
                    bret += sb;
                    cret += sc;
                }
            } else { // m>=n-k
                weight = I.gregory_omega(n - m);
                for (j = 0; j < m - k; j++) {
                    element_incr<T, SIZE1>(size1, cret, weight, atemp,
                                           bret); // cret += aret*bret
                    bret += sb;
                    cret += sc;
                }
            }
            // contribution from the stripe m-j <= k, with different weight
            j1 = m - k;
            if (j1 < 0)
                j1 = 0;
            for (j = j1; j <= m; j++) {
                weight = I.gregory_weights(n - j, n - m);
                element_incr<T, SIZE1>(size1, cret, weight, atemp,
                                       bret); // cret += aret*bret
                bret += sb;
                cret += sc;
            }
        }
        // CONTRIBUTION FROM BRET^CONJ:
        for (m = n - k; m < n; m++) {
            element_mult<T, SIZE1>(size1, aret, A.retptr(n, m), ft + m * sf);
            // without f ... aret = A.retptr(n,m);
            for (l = 0; l < sa; l++)
                atemp[l] = aret[l] * h; // here enters h
            for (j = m + 1; j <= n; j++) {
                weight = I.gregory_weights(n - j, n - m);
                element_conj<T, SIZE1>(size1, btemp, Bcc.retptr(j, m));
                // mind the minus sign: B(m,j) continued to -Bcc(j,m)*
                element_incr<T, SIZE1>(size1, result + j * sc, -weight, atemp, btemp);
            }
        }
    } else { // n < k
        for (j = 0; j <= n; j++) {
            cret = result + j * sc;
            for (m = 0; m <= k; m++) {
                weight = I.poly_integration(j, n, m) * h;
                if (m >= j) {
                    for (l = 0; l < sb; l++)
                        btemp[l] = B.retptr(m, j)[l];
                } else {
                    element_conj<T, SIZE1>(size1, btemp, Bcc.retptr(j, m));
                    weight *= -1;
                }
                if (n >= m) {
                    for (l = 0; l < sa; l++)
                        atemp[l] = A.retptr(n, m)[l];
                } else {
                    element_conj<T, SIZE1>(size1, atemp, Acc.retptr(m, n));
                    weight *= -1;
                }
                element_mult<T, SIZE1>(size1, aret, atemp, ft + m * sf);
                element_incr<T, SIZE1>(size1, cret, weight, aret, btemp);
                // without f element_incr<T,SIZE1>(size1,cret,weight,atemp,btemp);
            }
        }
    }

    cret = C.retptr(n, 0);
    for (l = 0; l < (n + 1) * sc; l++)
        cret[l] = result[l];
    delete[] result;
    delete[] atemp;
    delete[] aret;
    delete[] btemp;
    return;
}
/** \brief <b> Left-Mixing convolution of two matrix and a function at given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate
 * > \f{eqnarray*}{
 *     C &=& A*B\\
 *     C^{\rceil}(t,\tau') &=& \int_0^\beta ds  A^{\rceil}(t,s) f0(s)B^M(s - \tau')
 *                       + \int_0^t ds A^R(t,s) ft(s)B^{\rceil}(s,\tau') \f}
 * >
 * > At time-step `t = n h` for all \f$\tau'\f$.
 * > The result is written into `ctv`. 'f0' is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis, and
 * > 'ft' is a pointer to \f$F(t)\f$ on the real axis.
 * >
 * > Note: `ctv` is a pointer to the timestep of elements of type GG,
 * > so it must have size `(C.ntau_ + 1) * C.element_size_`
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * [int] given time-step
 * @param ctv
 * > [herm_matrix] left-mixing component \f$C^{\rceil}\f$ to which the result of the convolution is given.
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *f0
 * > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inverse temperature
 * @param h
 * > time step
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_tv(int n, std::complex<T> *ctv, GG &C, GG &A, GG &Acc,
                             std::complex<T> *f0, std::complex<T> *ft, GG &B, GG &Bcc,
                             integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1;
    int sf, sa, sb, sc, ntau, j, m, l, n1, size1 = C.size1();
    T weight, dtau;
    cplx *ctemp1, *ctv1, *btv, *atemp, *bmat, *atemp1;

    ntau = A.ntau();
    sa = A.element_size();
    sb = B.element_size();
    sc = C.element_size();
    sf = size1 * size1;
    dtau = beta / ntau;
    n1 = (n > k ? n : k);

    // check consistency:
    assert(sa * sb * sc != 0);
    assert(ntau >= k2);
    assert(Acc.element_size() == sa);
    assert(Bcc.element_size() == sb);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(C.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);

    // CONTRIBUTION FROM Atv * Bmat : very similar to computing the matsubara convolution
    ctemp1 = new cplx[sc];
    bmat = new cplx[(ntau + 1) * sb];
    for (m = 0; m <= ntau; m++)
        element_mult<T, SIZE1>(size1, bmat + m * sb, f0, B.matptr(m));
    for (m = 0; m <= ntau; m++) {
        matsubara_integral_2<T, SIZE1>(size1, m, ntau, ctemp1, A.tvptr(n, 0), bmat, I,
                                       B.sig());
        for (l = 0; l < sc; l++)
            ctv[m * sc + l] = dtau * ctemp1[l];
    }
    delete[] bmat;
    delete[] ctemp1;
    // CONTRIBUTION FROM Aret * Btv: loop over lines j
    atemp = new cplx[sa];
    atemp1 = new cplx[sa];
    n1 = (n > k ? n : k);
    for (j = 0; j <= n1; j++) {
        weight = I.gregory_weights(n, j);
        if (n < j) {
            element_conj<T, SIZE1>(size1, atemp, Acc.retptr(j, n));
            element_smul<T, SIZE1>(size1, atemp, -1);
        } else {
            element_set<T, SIZE1>(size1, atemp, A.retptr(n, j));
        }
        element_smul<T, SIZE1>(size1, atemp, h); // here enters h
        element_mult<T, SIZE1>(size1, atemp1, atemp, ft + sf * j);
        btv = B.tvptr(j, 0);
        ctv1 = ctv;
        if (weight != 1) {
            for (m = 0; m <= ntau; m++) {
                element_incr<T, SIZE1>(size1, ctv1, weight, atemp1, btv);
                btv += sb;
                ctv1 += sc;
            }
        } else {
            for (m = 0; m <= ntau; m++) {
                element_incr<T, SIZE1>(size1, ctv1, atemp1, btv);
                btv += sb;
                ctv1 += sc;
            }
        }
    }
    delete[] atemp;
    delete[] atemp1;
    return;
}
/** \brief <b> Left-Mixing convolution of two matrices and a function at given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 Here we calculate
 * > \f{eqnarray*}{
 *     C &=& A*B\\
 *     C^{\rceil}(t,\tau') &=& \int_0^\beta ds  A^{\rceil}(t,s) f0(s)B^M(s - \tau')
 *                       + \int_0^t ds A^R(t,s) ft(s)B^{\rceil}(s,\tau') \f}
 * >
 * > at time-step `t = n h` for all \f$\tau'\f$. 'f0' is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis, and
 * > 'ft' is a pointer to \f$F(t)\f$ on the real axis.
 * > The result is written at a given time step into 'C'. The objects A,B and C are of type GG.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * [int] given time-step
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *f0
 * > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inverse temperature
 * @param h
 * > time step
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_tv(int n, GG &C, GG &A, GG &Acc, std::complex<T> *f0,
                             std::complex<T> *ft, GG &B, GG &Bcc,
                             integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int ntau, m, sc, size1 = C.size1();
    cplx *ctv;
    ntau = C.ntau();
    sc = C.element_size();

    assert(sc > 0 && ntau > 0);
    ctv = new cplx[(ntau + 1) * sc];
    convolution_timestep_tv<T, GG, SIZE1>(n, ctv, C, A, Acc, f0, ft, B, Bcc, I, beta, h);
    for (m = 0; m <= ntau; m++)
        element_set<T, SIZE1>(size1, C.tvptr(n, m), ctv + m * sc);
    delete[] ctv;
}

/** \brief <b> Calculation of \f$C = A^{\rceil}*f0*B^{\lceil}\f$ at a given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate \f$C = A^{\rceil}*f0*B^{\lceil}\f$
 * > at a given time-step `t = n h`. 'f0' is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis.
 * >
 * > Note: 'cles' is a pointer to the timestep of elements of type GG.
 * > It must have size   (max(n,k)+1) * C.element_size_. The result is added to 'cles'.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * [int] given time-step
 * @param cles
 * > [std::complex] part of the lesser component \f$C^{<}\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *f0
 * > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_tvvt(int n, std::complex<T> *cles, GG &C, GG &A, GG &Acc,
                                   std::complex<T> *f0, GG &B, GG &Bcc,
                                   integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1;
    int sa, sb, sc, ntau, j, m, l, n1, sig, size1 = C.size1();
    T weight, dtau;
    cplx *atv, *btv, *btemp, *btemp1, *cles1, idtau;

    ntau = A.ntau();
    sa = A.element_size();
    sb = B.element_size();
    sc = C.element_size();
    dtau = beta / ntau;
    sig = A.sig();
    n1 = (n < k ? k : n);
    // check consistency:
    assert(sa * sb * sc != 0);
    assert(ntau >= k);
    assert(Acc.element_size() == sa);
    assert(Bcc.element_size() == sb);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(C.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);
    assert(sig == Acc.sig());
    assert(sig == B.sig());
    assert(sig == Bcc.sig());

    // contribution from Atv*Bvt = Atv(jh,tau) * Bcc^tv(nh,beta-tau)^* * (-Bose/Fermi)
    btemp = new cplx[(ntau + 1) * sb];
    btemp1 = new cplx[sb];
    idtau = cplx(0, -dtau);
    for (m = 0; m <= ntau; m++) {
        element_conj<T, SIZE1>(size1, btemp1, Bcc.tvptr(n, ntau - m));
        element_mult<T, SIZE1>(size1, btemp + m * sb, f0, btemp1);
    }
    delete[] btemp1;
    for (l = 0; l < (ntau + 1) * sb; l++)
        btemp[l] *= idtau * (-(T)sig);
    for (j = 0; j <= n1; j++) {
        btv = btemp;
        atv = A.tvptr(j, 0);
        cles1 = cles + j * sc;
        if (ntau < k2 - 1) {
            for (m = 0; m <= ntau; m++) {
                weight = I.gregory_weights(ntau, m);
                element_incr<T, SIZE1>(size1, cles1, weight, atv, btv);
                atv += sa;
                btv += sb;
            }
        } else {
            for (m = 0; m <= k; m++) {
                weight = I.gregory_omega(m);
                element_incr<T, SIZE1>(size1, cles1, weight, atv, btv);
                atv += sa;
                btv += sb;
            }
            for (m = k1; m < ntau - k; m++) {
                element_incr<T, SIZE1>(size1, cles1, atv, btv);
                atv += sa;
                btv += sb;
            }
            for (m = ntau - k; m <= ntau; m++) {
                weight = I.gregory_omega(ntau - m);
                element_incr<T, SIZE1>(size1, cles1, weight, atv, btv);
                atv += sa;
                btv += sb;
            }
        }
    }
    delete[] btemp;
    return;
}
/** \brief <b> Calculation of \f$C = A^{<}*ft*B^{A}\f$ at a given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate \f$C = A^{<}*ft*B^{A}\f$
 * > at a given time-step `t = n h`. 'ft' is a pointer to a function \f$F(t)\f$ on the real axis.
 * >
 * > Note: 'cles' is a pointer to the timestep of elements of type GG.
 * > It must have size   (max(n,k)+1) * C.element_size_. The result is added to 'cles'.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * > [int] given time-step
 * @param cles
 * > [std::complex] part of the lesser component \f$C^{<}\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_lesadv(int n, std::complex<T> *cles, GG &C, GG &A, GG &Acc,
                                     std::complex<T> *ft, GG &B, GG &Bcc,
                                     integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k();
    int sf, sa, sb, sc, ntau, j, m, l, n1, size1 = C.size1();
    T weight, dtau;
    cplx *atemp, *btemp, idtau, *ales, *badv, *btemp1;

    ntau = A.ntau();
    sa = A.element_size();
    sb = B.element_size();
    sc = C.element_size();
    sf = size1 * size1;
    dtau = beta / ntau;
    n1 = (n < k ? k : n);
    // check consistency:
    assert(sa * sb * sc != 0);
    assert(ntau >= k);
    assert(Acc.element_size() == sa);
    assert(Bcc.element_size() == sb);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(C.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);

    // contribution from Ales(j,m)*Badv(m,n)
    atemp = new cplx[(n1 + 1) * sa];
    btemp = new cplx[(n1 + 1) * sb];
    btemp1 = new cplx[sb];
    for (m = 0; m <= n1; m++) {
        weight = h * I.gregory_weights(n, m);
        if (m <= n) {
            element_conj<T, SIZE1>(size1, btemp1, Bcc.retptr(n, m));
            for (l = 0; l < sb; l++)
                btemp1[l] *= weight;
        } else {
            for (l = 0; l < sb; l++)
                btemp1[l] = -weight * B.retptr(m, n)[l];
        }
        element_mult<T, SIZE1>(size1, btemp + m * sb, ft + m * sf, btemp1);
    }
    delete[] btemp1;
    for (j = 0; j <= n1; j++) {
        for (m = 0; m <= n1; m++) {
            if (m < j) {
                element_conj<T, SIZE1>(size1, atemp + m * sa, Acc.lesptr(m, j));
                for (l = 0; l < sa; l++)
                    atemp[m * sa + l] *= -1;
            } else {
                for (l = 0; l < sa; l++)
                    atemp[m * sa + l] = A.lesptr(j, m)[l];
            }
        }
        ales = atemp;
        badv = btemp;
        for (m = 0; m <= n1; m++) {
            element_incr<T, SIZE1>(size1, cles + j * sc, ales, badv);
            ales += sa;
            badv += sb;
        }
    }
    delete[] atemp;
    delete[] btemp;
    return;
}
/** \brief <b> Calculation of \f$C = A^{R}*ft*B^{<}\f$ at a given time-step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate \f$C = A^{R}*ft*B^{<}\f$
 * > at a given time-step `t = n h`. 'ft' is a pointer to a function \f$F(t)\f$ on the real axis.
 * >
 * > Note: 'cles' is a pointer to the timestep of elements of type GG.
 * > It must have size   (max(n,k)+1) * C.element_size_. The result is added to 'cles'.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * [int] given time-step
 * @param cles
 * > [std::complex] part of the lesser component \f$C^{<}\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_retles(int n, std::complex<T> *cles, GG &C, GG &A, GG &Acc,
                                     std::complex<T> *ft, GG &B, GG &Bcc,
                                     integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1, size1 = C.size1();
    int sf, sa, sb, sc, ntau, j, m, l, n1;
    T weight, dtau;
    cplx *aret, *atemp, *btemp, *btemp1, *cles1, *bles, idtau;

    ntau = A.ntau();
    sa = A.element_size();
    sb = B.element_size();
    sc = C.element_size();
    sf = size1 * size1;
    dtau = beta / ntau;
    n1 = (n < k ? k : n);
    // check consistency:
    assert(sa * sb * sc != 0);
    assert(ntau >= k);
    assert(Acc.element_size() == sa);
    assert(Bcc.element_size() == sb);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(C.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);

    // contribution from Aret*Bles
    btemp = new cplx[(n1 + 1) * sb];
    btemp1 = new cplx[sb];
    atemp = new cplx[sa];
    for (m = 0; m <= n1; m++) { // btemp(m) --> B^<(m,n)
        if (m <= n) {
            for (l = 0; l < sb; l++)
                btemp1[l] = h * B.lesptr(m, n)[l];
        } else {
            element_conj<T, SIZE1>(size1, btemp1, Bcc.lesptr(n, m));
            for (l = 0; l < sb; l++)
                btemp1[l] *= -h;
        }
        element_mult<T, SIZE1>(size1, btemp + m * sb, ft + m * sf, btemp1);
    }
    delete[] btemp1;
    for (j = 0; j <= n; j++) { // compute -> cles1(j,n)
        cles1 = cles + j * sc;
        // CONTRINBUTION  FROM A_RET
        if (j >= k2 - 1) {
            aret = A.retptr(j, 0);
            bles = btemp;
            for (m = 0; m <= k; m++) {
                weight = I.gregory_omega(m);
                element_incr<T, SIZE1>(size1, cles1, weight, aret, bles);
                bles += sb;
                aret += sa;
            }
            for (m = k1; m < j - k; m++) {
                element_incr<T, SIZE1>(size1, cles1, aret, bles);
                bles += sb;
                aret += sa;
            }
            for (m = j - k; m <= j; m++) {
                weight = I.gregory_omega(j - m);
                element_incr<T, SIZE1>(size1, cles1, weight, aret, bles);
                bles += sb;
                aret += sa;
            }
        } else {
            aret = A.retptr(j, 0);
            bles = btemp;
            for (m = 0; m <= j; m++) {
                weight = I.gregory_weights(j, m);
                element_incr<T, SIZE1>(size1, cles1, weight, aret, bles);
                bles += sb;
                aret += sa;
            }
        }
        // CONTRIBUTION FROM ACC_RET
        if (j < k) {
            for (m = j + 1; m <= k; m++) {
                element_conj<T, SIZE1>(size1, atemp, Acc.retptr(m, j));
                weight = -I.gregory_weights(j, m);
                bles = btemp + m * sb;
                element_incr<T, SIZE1>(size1, cles1, weight, atemp, bles);
            }
        }
    }
    delete[] btemp;
    delete[] atemp;
    return;
}
/** \brief <b> Lesser convolution at a given time-step of two matrices and a function. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Here we calculate
 * >   \f{eqnarray*}{
 *     C^{<} &=& A^{\rceil}*f0*B^{\lceil} + A^<*ft*B^{A} + A^R*ft*B^<\\
 *     C^{<}(t,t') &=& -i \int_0^\beta ds  A^{\rceil}(t,s) f0 B^{\lceil}(s,t')+ \int_0^t ds  A^{<}(t,s) ft B^A(s,t')
 *                       + \int_0^t ds  A^{R}(t,s) ft B^<(s,t') \f}
 * > at a given time-step `t = n h`.
 * > 'f0' is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis, and
 * > 'ft' is a pointer to \f$F(t)\f$ on the real axis. The result is written into a contour Green's function `C`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * > [int] given time-step
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *f0
 * > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_les(int n, GG &C, GG &A, GG &Acc, std::complex<T> *f0,
                              std::complex<T> *ft, GG &B, GG &Bcc,
                              integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    cplx *cles;
    int m, sc, n1, k = I.get_k(), size1 = C.size1();
    sc = C.element_size();
    assert(sc > 0);
    n1 = (k > n ? k : n);
    cles = new cplx[(n1 + 1) * sc];
    for (m = 0; m < sc * (n1 + 1); m++)
        cles[m] = 0;
    convolution_timestep_les_tvvt<T, GG, SIZE1>(n, cles, C, A, Acc, f0, B, Bcc, I, beta, h);
    convolution_timestep_les_lesadv<T, GG, SIZE1>(n, cles, C, A, Acc, ft, B, Bcc, I, beta,
                                                  h);
    convolution_timestep_les_retles<T, GG, SIZE1>(n, cles, C, A, Acc, ft, B, Bcc, I, beta,
                                                  h);
    for (m = 0; m <= n; m++)
        element_set<T, SIZE1>(size1, C.lesptr(m, n), cles + m * sc);
    delete[] cles;
    return;
}
/** \brief <b> Returns convolution of two matrices and a function at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*fxB of the objects 'A' and 'B' with class 'herm_matrix' and a function 'f'
* > at a given time step 't=nh'. If 'n=-1', one uses 'f0' (which is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis),
* >  otherwise 'ft' (which is a pointer to a function \f$F(t)\f$ on the real axis). Works for a scalar and square matrices.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] complex conjugate to A
* @param *f0
* > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
* @param *ft
* > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] complex conjugate to B
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                          std::complex<T> *f0, std::complex<T> *ft, herm_matrix<T> &B,
                          herm_matrix<T> &Bcc, integration::Integrator<T> &I, T beta, T h) {
    int size1 = C.size1(), ntau = C.ntau(), k = I.k(), n1 = (n < k ? k : n);
    if (n == -1) {
        convolution_matsubara(C, A, f0, B, I, beta);
        return;
    }
    assert(n >= 0);
    assert(A.size1() == size1);
    assert(Acc.size1() == size1);
    assert(B.size1() == size1);
    assert(Bcc.size1() == size1);
    assert(A.ntau() == ntau);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);
    assert(C.nt() >= n);
    if (size1 == 1) {
        convolution_timestep_ret<T, herm_matrix<T>, 1>(n, C, A, Acc, ft, B, Bcc, I, h);
        convolution_timestep_tv<T, herm_matrix<T>, 1>(n, C, A, Acc, f0, ft, B, Bcc, I, beta,
                                                      h);
        convolution_timestep_les<T, herm_matrix<T>, 1>(n, C, A, Acc, f0, ft, B, Bcc, I, beta,
                                                       h);
    } else {
        convolution_timestep_ret<T, herm_matrix<T>, LARGESIZE>(n, C, A, Acc, ft, B, Bcc, I,
                                                               h);
        convolution_timestep_tv<T, herm_matrix<T>, LARGESIZE>(n, C, A, Acc, f0, ft, B, Bcc,
                                                              I, beta, h);
        convolution_timestep_les<T, herm_matrix<T>, LARGESIZE>(n, C, A, Acc, f0, ft, B, Bcc,
                                                               I, beta, h);
    }
}
/** \brief <b> Returns convolution of two hermitian matrices and a function at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*fxB of the hermitian matrices 'A' and 'B' and a function 'f'
* > at a given time step 't=nh'. If 'n=-1', one uses 'f0' (which is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis),
* >  otherwise 'ft' (which is a pointer to \f$F(t)\f$ on the real axis). Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function (A=Acc)
* @param *f0
* > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
* @param *ft
* > [std::complex] Pointer to \f$F(t)\f$ on the real axis
* @param B
* > [herm_matrix] contour Green's function (B=Bcc)
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, std::complex<T> *f0,
                          std::complex<T> *ft, herm_matrix<T> &B,
                          integration::Integrator<T> &I, T beta, T h) {
    convolution_timestep<T>(n, C, A, A, f0, ft, B, B, I, beta, h);
}

/** \brief <b> Returns the result of the contour convolution of two matrices and a function object</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calls convolution routines to compute contour convolution C=A*fxB.
* > The objects A,B and C are of the class 'herm_matrix'.
* > On Matsubara axis one calls 'f0' (which is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis),
* >  otherwise 'ft' (which is a pointer to \f$F(t)\f$ on the real axis). Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param C
* > [herm_matrix] Matrix to which the result of the convolution is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] complex conjugate to A
* @param *f0
* > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (complex function)
* @param *ft
* > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] complex conjugate to B
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time interval
*/
template <typename T>
void convolution(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                 std::complex<T> *f0, std::complex<T> *ft, herm_matrix<T> &B,
                 herm_matrix<T> &Bcc, integration::Integrator<T> &I, T beta, T h) {
    int tstp;
    convolution_matsubara(C, A, f0, B, I, beta);
    for (tstp = 0; tstp <= C.nt(); tstp++)
        convolution_timestep<T>(tstp, C, A, Acc, f0, ft, B, Bcc, I, beta, h);
}

/** \brief <b> Returns convolution of two 'herm_matrix' objects and a contour function at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution \f$C=A*FxB \f$ of the 'herm_matrix' object 'A' and 'B' and a time-dependent contour function \f$F(t)\f$
* > at a given time step 't=nh'. If 'n=-1', one performs Matsubara convolution with a function \f$F(-1)\f$.
* > Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] complex conjugate to A
* @param ft
* > [function] contour function F(t).
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] complex conjugate to B
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time interval
*/
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                          function<T> &ft, herm_matrix<T> &B, herm_matrix<T> &Bcc,
                          integration::Integrator<T> &I, T beta, T h) {
    int size1 = C.size1(), ntau = C.ntau(), k = I.k(), n1 = (n < k ? k : n);
    assert(ft.size1() == size1 && ft.nt() >= -1);
    if (n == -1) {
        convolution_matsubara(C, A, ft.ptr(-1), B, I, beta);
        return;
    }
    assert(n >= 0);
    assert(A.size1() == size1);
    assert(Acc.size1() == size1);
    assert(B.size1() == size1);
    assert(Bcc.size1() == size1);
    assert(A.ntau() == ntau);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);
    assert(ft.nt() >= n1);
    assert(C.nt() >= n);
    if (size1 == 1) {
        convolution_timestep_ret<T, herm_matrix<T>, 1>(n, C, A, Acc, ft.ptr(0), B, Bcc, I,
                                                       h);
        convolution_timestep_tv<T, herm_matrix<T>, 1>(n, C, A, Acc, ft.ptr(-1), ft.ptr(0), B,
                                                      Bcc, I, beta, h);
        convolution_timestep_les<T, herm_matrix<T>, 1>(n, C, A, Acc, ft.ptr(-1), ft.ptr(0),
                                                       B, Bcc, I, beta, h);
    } else {
        convolution_timestep_ret<T, herm_matrix<T>, LARGESIZE>(n, C, A, Acc, ft.ptr(0), B,
                                                               Bcc, I, h);
        convolution_timestep_tv<T, herm_matrix<T>, LARGESIZE>(n, C, A, Acc, ft.ptr(-1),
                                                              ft.ptr(0), B, Bcc, I, beta, h);
        convolution_timestep_les<T, herm_matrix<T>, LARGESIZE>(
            n, C, A, Acc, ft.ptr(-1), ft.ptr(0), B, Bcc, I, beta, h);
    }
}


/** \brief <b> Returns convolution of two 'herm_matrix' objects and a contour function at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution \f$C=A*FxB \f$ of the 'herm_matrix' object 'A' and 'B' and a time-dependent contour function \f$F(t)\f$
* > at a given time step 't=nh'. If 'n=-1', one performs Matsubara convolution with a function \f$F(-1)\f$.
* > Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] complex conjugate to A
* @param ft
* > [function] contour function F(t).
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] complex conjugate to B
* @param beta
* > inversed temperature
* @param h
* > time interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                          function<T> &ft, herm_matrix<T> &B, herm_matrix<T> &Bcc,
                          T beta, T h, int SolveOrder) {
    int size1 = C.size1(), ntau = C.ntau(), n1 = (n < SolveOrder ? SolveOrder : n);
    assert(ft.size1() == size1 && ft.nt() >= -1);
    if (n == -1) {
        convolution_matsubara(C, A, ft.ptr(-1), B, integration::I<T>(SolveOrder), beta);
        return;
    }
    assert(n >= 0);
    assert(A.size1() == size1);
    assert(Acc.size1() == size1);
    assert(B.size1() == size1);
    assert(Bcc.size1() == size1);
    assert(A.ntau() == ntau);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);
    assert(ft.nt() >= n1);
    assert(C.nt() >= n);
    if (size1 == 1) {
        convolution_timestep_ret<T, herm_matrix<T>, 1>(n, C, A, Acc, ft.ptr(0), B, Bcc, integration::I<T>(SolveOrder),
                                                       h);
        convolution_timestep_tv<T, herm_matrix<T>, 1>(n, C, A, Acc, ft.ptr(-1), ft.ptr(0), B,
                                                      Bcc, integration::I<T>(SolveOrder), beta, h);
        convolution_timestep_les<T, herm_matrix<T>, 1>(n, C, A, Acc, ft.ptr(-1), ft.ptr(0),
                                                       B, Bcc, integration::I<T>(SolveOrder), beta, h);
    } else {
        convolution_timestep_ret<T, herm_matrix<T>, LARGESIZE>(n, C, A, Acc, ft.ptr(0), B,
                                                               Bcc, integration::I<T>(SolveOrder), h);
        convolution_timestep_tv<T, herm_matrix<T>, LARGESIZE>(n, C, A, Acc, ft.ptr(-1),
                                                              ft.ptr(0), B, Bcc, integration::I<T>(SolveOrder), beta, h);
        convolution_timestep_les<T, herm_matrix<T>, LARGESIZE>(
            n, C, A, Acc, ft.ptr(-1), ft.ptr(0), B, Bcc, integration::I<T>(SolveOrder), beta, h);
    }
}



/** \brief <b> Returns convolution of two hermitian matrices and a contour function at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*fxB of the hermitian matrices 'A' and 'B' and a time-dependent contour function 'F(t)'
* > at a given time step 't=nh'. If 'n=-1', one performs Matsubara convolution with a function 'F(-1)', otherwise with 'F(t)'
* > Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function (A=Acc)
* @param ft
* > [function] the contour function F(t).
* @param B
* > [herm_matrix] contour Green's function (B=Bcc)
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, function<T> &ft,
                          herm_matrix<T> &B, integration::Integrator<T> &I, T beta, T h) {
    convolution_timestep<T>(n, C, A, A, ft, B, B, I, beta, h);
}

/** \brief <b> Returns convolution of two hermitian matrices and a contour function at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*fxB of the hermitian matrices 'A' and 'B' and a time-dependent contour function 'F(t)'
* > at a given time step 't=nh'. If 'n=-1', one performs Matsubara convolution with a function 'F(-1)', otherwise with 'F(t)'
* > Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function (A=Acc)
* @param ft
* > [function] the contour function F(t).
* @param B
* > [herm_matrix] contour Green's function (B=Bcc)
* @param beta
* > inversed temperature
* @param h
* > time step interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T>
void convolution_timestep(int n, herm_matrix<T> &C, herm_matrix<T> &A, function<T> &ft,
                          herm_matrix<T> &B, T beta, T h, int SolveOrder) {
    convolution_timestep<T>(n, C, A, A, ft, B, B, beta, h, SolveOrder);
}

/** \brief <b> Returns the result of the contour convolution of two matrices and a function object</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calls convolution routines to compute contour convolution C=A*fxB.
* > The objects A,B and C are of the class 'herm_matrix'.
* > If 'n=-1', one performs Matsubara convolution with a function 'F(-1)', otherwise with 'F(t)'.
* > Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param C
* > [herm_matrix] Matrix to which the result of the convolution is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] complex conjugate to A
* @param ft
* > [std::complex] complex function
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] complex conjugate to B
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time interval
*/
template <typename T>
void convolution(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc, function<T> &ft,
                 herm_matrix<T> &B, herm_matrix<T> &Bcc, integration::Integrator<T> &I,
                 T beta, T h) {
    int tstp;
    convolution_matsubara(C, A, ft.ptr(-1), B, I, beta);
    for (tstp = 0; tstp <= C.nt(); tstp++)
        convolution_timestep<T>(tstp, C, A, Acc, ft, B, Bcc, I, beta, h);
}


/** \brief <b> Returns the result of the contour convolution of two matrices and a function object</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calls convolution routines to compute contour convolution C=A*fxB.
* > The objects A,B and C are of the class 'herm_matrix'.
* > If 'n=-1', one performs Matsubara convolution with a function 'F(-1)', otherwise with 'F(t)'.
* > Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param C
* > [herm_matrix] Matrix to which the result of the convolution is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] complex conjugate to A
* @param ft
* > [std::complex] complex function
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] complex conjugate to B
* @param beta
* > inversed temperature
* @param h
* > time interval
* @param SolveOrder
* > [int] integrator order

*/
template <typename T>
void convolution(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc, function<T> &ft,
                 herm_matrix<T> &B, herm_matrix<T> &Bcc,
                 T beta, T h, int SolveOrder) {
    int tstp;
    convolution_matsubara(C, A, ft.ptr(-1), B, integration::I<T>(SolveOrder), beta);
    for (tstp = 0; tstp <= C.nt(); tstp++)
        convolution_timestep<T>(tstp, C, A, Acc, ft, B, Bcc, beta, h, SolveOrder);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// same funcions, but compute only at one timepoint

/// @private
/** \brief <b> Performs the Matsubara convolution of two matrices and a function at one timepoint, i.e. C=A*fxB</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Same function, as 'convolution_matsubara_dispatch', but performs Matsubara convolution \f$C^M=A^M*fxB^M\f$ only at one timepoint
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param m1
* [int] time step
* @param cc
* > [std::complex] Matsubara component \f$C^M\f$
* @param sizec
* > [int] size of the matrix 'C'
* @param A
* > [GG] contour Green's function
* @param *f0
* > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
* @param B
* > [GG] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
*/
template <typename T, class GG, int SIZE1>
void convolution_matsubara_tau_dispatch(int m1, std::complex<T> *cc, int sizec, GG &A,
                                        std::complex<T> *f0, GG &B,
                                        integration::Integrator<T> &I, T beta) {
    int ntau, l, m, size1 = sizec, sb = B.element_size();
    std::complex<T> *bmat;
    T dtau;
    ntau = A.ntau();
    bmat = new std::complex<T>[sb * (ntau + 1)];
    if (f0 != NULL) {
        for (m = 0; m <= ntau; m++)
            element_mult<T, SIZE1>(size1, bmat + m * sb, f0, B.matptr(m));
    } else {
        for (m = 0; m <= ntau; m++)
            element_set<T, SIZE1>(size1, bmat + m * sb, B.matptr(m));
    }
    // compute cmat(m*dtau) = int_0^beta dx amat(tau-x) b(x)
    matsubara_integral_1<T, SIZE1>(size1, m1, ntau, cc, A.matptr(0), bmat, I, A.sig());
    delete[] bmat;
    // multiply by dtau:
    dtau = beta / ntau;
    for (l = 0; l < sizec * sizec; l++)
        cc[l] *= dtau;
    return;
}
/** \brief <b> Calculation of \f$C = A^{<}*ft*B^{A}\f$ at one timepoint. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Same function, as 'convolution_timestep_les_lesadv', but performs calculations only at one timepoint
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param j
 * [int] given time-step
 * @param n
 * [int] given time-step
 * @param cc
 * > [std::complex] part of the lesser component \f$C^{<}\f$
 * @param sizec
 * > [int] size of the matrix 'C'
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_jn_lesadv(int j, int n, std::complex<T> *cc, int sizec, GG &A,
                                        GG &Acc, std::complex<T> *ft, GG &B, GG &Bcc,
                                        integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k();
    int sf, sa, sb, sc, ntau, m, l, n1, size1 = sizec;
    T weight, dtau;
    cplx *atemp, *btemp, idtau, *ales, *badv, *btemp1;

    ntau = A.ntau();
    sa = A.element_size();
    sb = B.element_size();
    sc = sizec * sizec;
    sf = size1 * size1;
    dtau = beta / ntau;
    n1 = (n < k ? k : n);

    atemp = new cplx[(n1 + 1) * sa];
    btemp = new cplx[(n1 + 1) * sb];
    btemp1 = new cplx[sb];
    for (m = 0; m <= n1; m++) {
        weight = h * I.gregory_weights(n, m);
        if (m <= n) {
            element_conj<T, SIZE1>(size1, btemp1, Bcc.retptr(n, m));
            for (l = 0; l < sb; l++)
                btemp1[l] *= weight;
        } else {
            for (l = 0; l < sb; l++)
                btemp1[l] = -weight * B.retptr(m, n)[l];
        }
        if (ft != NULL) {
            element_mult<T, SIZE1>(size1, btemp + m * sb, ft + m * sf, btemp1);
        } else {
            element_set<T, SIZE1>(size1, btemp + m * sb, btemp1);
        }
    }
    delete[] btemp1;
    // for(j=0;j<=n1;j++){
    {
        for (m = 0; m <= n1; m++) {
            if (m < j) {
                element_conj<T, SIZE1>(size1, atemp + m * sa, Acc.lesptr(m, j));
                for (l = 0; l < sa; l++)
                    atemp[m * sa + l] *= -1;
            } else {
                for (l = 0; l < sa; l++)
                    atemp[m * sa + l] = A.lesptr(j, m)[l];
            }
        }
        ales = atemp;
        badv = btemp;
        for (m = 0; m <= n1; m++) {
            element_incr<T, SIZE1>(size1, cc, ales, badv);
            ales += sa;
            badv += sb;
        }
    }
    delete[] atemp;
    delete[] btemp;
    return;
}
/** \brief <b> Calculation of \f$C = A^{\rceil}*f0*B^{\lceil}\f$ at one timepoint. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Same function, as 'convolution_timestep_les_tvvt', but performs calculations only at one timepoint
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param j
 * [int] given time-step
 * @param n
 * [int] given time-step
 * @param cc
 * > [std::complex] part of the lesser component \f$C^{<}\f$
 * @param sizec
 * > [int] size of the matrix 'C'
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *f0
 * > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_jn_tvvt(int j, int n, std::complex<T> *cc, int sizec, GG &A,
                                      GG &Acc, std::complex<T> *f0, GG &B, GG &Bcc,
                                      integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1;
    int sa, sb, sc, ntau, m, l, n1, sig, size1 = sizec;
    T weight, dtau;
    cplx *atv, *btv, *btemp, *btemp1, *cles1, idtau;

    ntau = A.ntau();
    sa = A.element_size();
    sb = B.element_size();
    sc = sizec * sizec;
    dtau = beta / ntau;
    sig = A.sig();
    n1 = (n < k ? k : n);
    // contribution from Atv*Bvt = Atv(jh,tau) * Bcc^tv(nh,beta-tau)^* * (-Bose/Fermi)
    btemp = new cplx[(ntau + 1) * sb];
    btemp1 = new cplx[sb];
    idtau = cplx(0, -dtau);
    for (m = 0; m <= ntau; m++) {
        element_conj<T, SIZE1>(size1, btemp1, Bcc.tvptr(n, ntau - m));
        if (f0 != NULL) {
            element_mult<T, SIZE1>(size1, btemp + m * sb, f0, btemp1);
        } else {
            element_set<T, SIZE1>(size1, btemp + m * sb, btemp1);
        }
    }
    delete[] btemp1;
    for (l = 0; l < (ntau + 1) * sb; l++)
        btemp[l] *= idtau * (-(T)sig);
    // for(j=0;j<=n1;j++){
    {
        btv = btemp;
        atv = A.tvptr(j, 0);
        cles1 = cc;
        if (ntau < k2 - 1) {
            for (m = 0; m <= ntau; m++) {
                weight = I.gregory_weights(ntau, m);
                element_incr<T, SIZE1>(size1, cles1, weight, atv, btv);
                atv += sa;
                btv += sb;
            }
        } else {
            for (m = 0; m <= k; m++) {
                weight = I.gregory_omega(m);
                element_incr<T, SIZE1>(size1, cles1, weight, atv, btv);
                atv += sa;
                btv += sb;
            }
            for (m = k1; m < ntau - k; m++) {
                element_incr<T, SIZE1>(size1, cles1, atv, btv);
                atv += sa;
                btv += sb;
            }
            for (m = ntau - k; m <= ntau; m++) {
                weight = I.gregory_omega(ntau - m);
                element_incr<T, SIZE1>(size1, cles1, weight, atv, btv);
                atv += sa;
                btv += sb;
            }
        }
    }
    delete[] btemp;
    return;
}
/** \brief <b> Calculation of \f$C = A^{R}*ft*B^{<}\f$ at one timepoint. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Same function, as 'convolution_timestep_les_retles', but performs calculations only at one timepoint
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param j
 * > [int] given time-step
 * @param n
 * > [int] given time-step
 * @param cc
 * > [std::complex] part of the lesser component \f$C^{<}\f$
 * @param sizec
 * > [int] size of the matrix 'C'
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_jn_retles(int j, int n, std::complex<T> *cc, int sizec, GG &A,
                                        GG &Acc, std::complex<T> *ft, GG &B, GG &Bcc,
                                        integration::Integrator<T> &I, T beta, T h) {
    typedef std::complex<T> cplx;
    int k = I.get_k(), k1 = k + 1, k2 = 2 * k1, size1 = sizec;
    int sf, sa, sb, sc, ntau, m, l, n1;
    T weight, dtau;
    cplx *aret, *atemp, *btemp, *btemp1, *cles1, *bles, idtau;

    ntau = A.ntau();
    sa = A.element_size();
    sb = B.element_size();
    sc = sizec * sizec;
    sf = size1 * size1;
    dtau = beta / ntau;
    n1 = (n < k ? k : n);
    // contribution from Aret*Bles
    btemp = new cplx[(n1 + 1) * sb];
    btemp1 = new cplx[sb];
    atemp = new cplx[sa];
    for (m = 0; m <= n1; m++) { // btemp(m) --> B^<(m,n)
        if (m <= n) {
            for (l = 0; l < sb; l++)
                btemp1[l] = h * B.lesptr(m, n)[l];
        } else {
            element_conj<T, SIZE1>(size1, btemp1, Bcc.lesptr(n, m));
            for (l = 0; l < sb; l++)
                btemp1[l] *= -h;
        }
        if (ft != NULL) {
            element_mult<T, SIZE1>(size1, btemp + m * sb, ft + m * sf, btemp1);
        } else {
            element_set<T, SIZE1>(size1, btemp + m * sb, btemp1);
        }
    }
    delete[] btemp1;
    // for(j=0;j<=n;j++){ // compute -> cles1(j,n)
    {
        cles1 = cc;
        // CONTRINBUTION  FROM A_RET
        if (j >= k2 - 1) {
            aret = A.retptr(j, 0);
            bles = btemp;
            for (m = 0; m <= k; m++) {
                weight = I.gregory_omega(m);
                element_incr<T, SIZE1>(size1, cles1, weight, aret, bles);
                bles += sb;
                aret += sa;
            }
            for (m = k1; m < j - k; m++) {
                element_incr<T, SIZE1>(size1, cles1, aret, bles);
                bles += sb;
                aret += sa;
            }
            for (m = j - k; m <= j; m++) {
                weight = I.gregory_omega(j - m);
                element_incr<T, SIZE1>(size1, cles1, weight, aret, bles);
                bles += sb;
                aret += sa;
            }
        } else {
            aret = A.retptr(j, 0);
            bles = btemp;
            for (m = 0; m <= j; m++) {
                weight = I.gregory_weights(j, m);
                element_incr<T, SIZE1>(size1, cles1, weight, aret, bles);
                bles += sb;
                aret += sa;
            }
        }
        // CONTRIBUTION FROM ACC_RET
        if (j < k) {
            for (m = j + 1; m <= k; m++) {
                element_conj<T, SIZE1>(size1, atemp, Acc.retptr(m, j));
                weight = -I.gregory_weights(j, m);
                bles = btemp + m * sb;
                element_incr<T, SIZE1>(size1, cles1, weight, atemp, bles);
            }
        }
    }
    delete[] btemp;
    delete[] atemp;
    return;
}
/** \brief <b> Lesser convolution  at one timepoint of two matrices and a function. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Same function, as 'convolution_timestep_les', but performs calculations only at one timepoint
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param j
 * > [int] given time-step
 * @param n
 * [int] given time-step
 * > @param cc
 * > [std::complex] part of the lesser component \f$C^{<}\f$
 * @param sizec
 * > [int] size of the matrix 'C'
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *f0
 * > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void convolution_timestep_les_jn(int j, int n, std::complex<T> *cc, int sizec, GG &A,
                                 GG &Acc, std::complex<T> *f0, std::complex<T> *ft, GG &B,
                                 GG &Bcc, integration::Integrator<T> &I, T beta, T h) {
    element_set_zero<T, LARGESIZE>(sizec, cc);
    convolution_timestep_les_jn_tvvt<T, GG, SIZE1>(j, n, cc, sizec, A, Acc, f0, B, Bcc, I,
                                                   beta, h);
    convolution_timestep_les_jn_lesadv<T, GG, SIZE1>(j, n, cc, sizec, A, Acc, ft, B, Bcc, I,
                                                     beta, h);
    convolution_timestep_les_jn_retles<T, GG, SIZE1>(j, n, cc, sizec, A, Acc, ft, B, Bcc, I,
                                                     beta, h);
    return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// density matrix
// note: !! works for both (GG=herm_matrix and GG=herm_pseudo
// (because for density matrix, restricted convolution is the same as full convolution)

/** \brief <b> Returns the result of the contour convolution and a contour function for a density matrix at a given time-step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calls convolution routines to compute contour convolution \f$\rho=-i (A*FxB)^<\f$
* > of the objects 'A' and 'B' and a time-dependent contour function \f$F(t)\f$ at a given time step 't=nh'.
* > If 'n=-1', one performs the Matsubara convolution.
* > The objects 'A' and 'B' are of the class type 'GG'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param rho
* > [std::complex] complex function to which the result of the convolution is given
* @param A
* > [GG] contour Green's function
* @param Acc
* > [GG] complex conjugate to A
* @param ft
* > [function] contour function F(t).
* @param B
* > [GG] contour Green's function
* @param Bcc
* > [GG] complex conjugate to B
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time interval
*/
template <typename T, class GG>
void convolution_density_matrix(int n, std::complex<T> *rho, GG &A, GG &Acc, function<T> &ft,
                                GG &B, GG &Bcc, integration::Integrator<T> &I, T beta, T h) {
    int size1 = A.size1(), ntau = A.ntau(), k = I.k(), n1 = (n < k && n >= 0 ? k : n);
    assert(ft.size1() == size1 && ft.nt() >= n1);
    assert(A.size1() == size1);
    assert(Acc.size1() == size1);
    assert(B.size1() == size1);
    assert(Bcc.size1() == size1);
    assert(A.ntau() == ntau);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);
    assert(ft.nt() >= n1);
    if (n == -1) {
        if (size1 == 1) {
            // note:
            convolution_matsubara_tau_dispatch<T, GG, 1>(ntau, rho, size1, A, ft.ptr(-1), B,
                                                         I, beta);
        } else {
            convolution_matsubara_tau_dispatch<T, GG, LARGESIZE>(ntau, rho, size1, A,
                                                                 ft.ptr(-1), B, I, beta);
        }
        element_smul<T, LARGESIZE>(size1, rho, -1.0);
    } else {
        if (size1 == 1) {
            convolution_timestep_les_jn<T, GG, 1>(n, n, rho, size1, A, Acc, ft.ptr(-1),
                                                  ft.ptr(0), B, Bcc, I, beta, h);
        } else {
            convolution_timestep_les_jn<T, GG, LARGESIZE>(
                n, n, rho, size1, A, Acc, ft.ptr(-1), ft.ptr(0), B, Bcc, I, beta, h);
        }
        element_smul<T, LARGESIZE>(size1, rho, std::complex<T>(0, -1.0));
    }
}

/** \brief <b> Returns the result of the contour convolution and a contour function for a density matrix at a given time-step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calls convolution routines to compute contour convolution \f$\rho=-i (A*FxB)^<\f$
* > of the objects 'A' and 'B' and a time-dependent contour function \f$F(t)\f$ at a given time step 't=nh'.
* > If 'n=-1', one performs the Matsubara convolution.
* > The objects 'A' and 'B' are of the class type 'GG'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param rho
* > [cdmatrix] complex function to which the result of the convolution is given
* @param A
* > [GG] contour Green's function
* @param Acc
* > [GG] complex conjugate to A
* @param ft
* > [function] contour function F(t).
* @param B
* > [GG] contour Green's function
* @param Bcc
* > [GG] complex conjugate to B
* @param beta
* > inversed temperature
* @param h
* > time interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T, class GG>
void convolution_density_matrix(int n, cdmatrix &rho, GG &A, GG &Acc, function<T> &ft,
                                GG &B, GG &Bcc, T beta, T h, int SolveOrder) {
    int size1 = A.size1();
    int element_size = size1;
    std::complex<T> *rho_ptr;
    rho_ptr = new std::complex<T>[element_size];

    convolution_density_matrix<T, GG>(n, rho_ptr, A, Acc, ft, B, Bcc, integration::I<T>(SolveOrder), beta, h);
    map_ptr2matrix<T>(size1, size1, rho_ptr, rho);

    delete[] rho_ptr;
}

/** \brief <b> Returns the result of the contour convolution and a contour function for a density matrix</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calls convolution routine to compute contour convolution \f$\rho=-i (A\ast F\ast B)^<\f$
* > of the objects 'A' and 'B' and a time-dependent contour function \f$F(t)\f$.
* > The objects 'A' and 'B' are of the class type 'GG'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > [int] time step
* @param rho
* > [std::complex] complex function to which the result of the convolution is given
* @param A
* > [GG] contour Green's function
* @param ft
* > [function] contour function F(t).
* @param B
* > [GG] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time interval
*/
template <typename T, class GG>
void convolution_density_matrix(int tstp, std::complex<T> *rho, GG &A, function<T> &ft,
                                GG &B, integration::Integrator<T> &I, T beta, T h) {
    convolution_density_matrix<T, GG>(tstp, rho, A, A, ft, B, B, I, beta, h);
}


/** \brief <b> Returns the result of the contour convolution and a contour function for a density matrix at a given time-step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calls convolution routines to compute contour convolution \f$\rho=-i (A\ast F\ast B)^<\f$
* > of the objects 'A' and 'B' and a time-dependent contour function \f$F(t)\f$ at a given time step 't=nh'.
* > If 'n=-1', one performs the Matsubara convolution.
* > The objects 'A' and 'B' are of the class type 'GG'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param rho
* > [cdmatrix] complex function to which the result of the convolution is given
* @param A
* > [GG] contour Green's function
* @param ft
* > [function] contour function F(t).
* @param B
* > [GG] contour Green's function
* @param beta
* > inversed temperature
* @param h
* > time interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T, class GG>
void convolution_density_matrix(int n, cdmatrix &rho, GG &A, function<T> &ft,
                                GG &B, T beta, T h, int SolveOrder) {
    int size1 = A.size1();
    int element_size = size1;
    std::complex<T> *rho_ptr;
    rho_ptr = new std::complex<T>[element_size];

    convolution_density_matrix<T, GG>(n, rho_ptr, A, A, ft, B, B, integration::I<T>(SolveOrder), beta, h);
    map_ptr2matrix<T>(size1, size1, rho_ptr, rho);

    delete[] rho_ptr;
}


/** \brief <b> Returns the result of the contour convolution for a density matrix at a given time-step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calls convolution routines to compute contour convolution \f$\rho=-i (A*B)^<\f$ at a given time step 't=nh'.
* > The objects A and B are of the class type 'GG'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param rho
* > [std::complex] Matrix to which the result of the convolution is given
* @param A
* > [GG] contour Green's function
* @param Acc
* > [GG] complex conjugate to A
* @param B
* > [GG] contour Green's function
* @param Bcc
* > [GG] complex conjugate to B
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time interval
*/
template <typename T, class GG>
void convolution_density_matrix(int n, std::complex<T> *rho, GG &A, GG &Acc, GG &B, GG &Bcc,
                                integration::Integrator<T> &I, T beta, T h) {
    int size1 = A.size1(), ntau = A.ntau(), k = I.k(), n1 = (n < k && n >= 0 ? k : n);
    assert(A.size1() == size1);
    assert(Acc.size1() == size1);
    assert(B.size1() == size1);
    assert(Bcc.size1() == size1);
    assert(A.ntau() == ntau);
    assert(Acc.ntau() == ntau);
    assert(B.ntau() == ntau);
    assert(Bcc.ntau() == ntau);
    assert(A.nt() >= n1);
    assert(Acc.nt() >= n1);
    assert(B.nt() >= n1);
    assert(Bcc.nt() >= n1);
    if (n == -1) {
        if (size1 == 1) {
            // note:
            convolution_matsubara_tau_dispatch<T, GG, 1>(ntau, rho, size1, A, NULL, B, I,
                                                         beta);
        } else {
            convolution_matsubara_tau_dispatch<T, GG, LARGESIZE>(ntau, rho, size1, A, NULL,
                                                                 B, I, beta);
        }
        element_smul<T, LARGESIZE>(size1, rho, -1.0);
    } else {
        if (size1 == 1) {
            convolution_timestep_les_jn<T, GG, 1>(n, n, rho, size1, A, Acc, NULL, NULL, B,
                                                  Bcc, I, beta, h);
        } else {
            convolution_timestep_les_jn<T, GG, LARGESIZE>(n, n, rho, size1, A, Acc, NULL,
                                                          NULL, B, Bcc, I, beta, h);
        }
        element_smul<T, LARGESIZE>(size1, rho, std::complex<T>(0, -1.0));
    }
}


/** \brief <b> Returns the result of the contour convolution for a density matrix at a given time-step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calls convolution routines to compute contour convolution \f$\rho=-i (A*B)^<\f$ at a given time step 't=nh'.
* > The objects A and B are of the class type 'GG'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > [int] number of the time step ('t=nh')
* @param rho
* > [std::complex] Matrix to which the result of the convolution is given
* @param A
* > [GG] contour Green's function
* @param Acc
* > [GG] complex conjugate to A
* @param B
* > [GG] contour Green's function
* @param Bcc
* > [GG] complex conjugate to B
* @param beta
* > inversed temperature
* @param h
* > time interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T, class GG>
void convolution_density_matrix(int n, cdmatrix &rho, GG &A, GG &Acc, GG &B, GG &Bcc, 
    T beta, T h, int SolveOrder) {

    int size1 = A.size1();
    int element_size = size1;
    std::complex<T> *rho_ptr;
    rho_ptr = new std::complex<T>[element_size];

    convolution_density_matrix<T, GG>(n, rho_ptr, A, Acc, B, Bcc, integration::I<T>(SolveOrder), beta, h);
    map_ptr2matrix<T>(size1, size1, rho_ptr, rho);

    delete[] rho_ptr;

}


/** \brief <b> Returns the result of the contour convolution for a density matrix</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calls convolution routines to compute contour convolution \f$\rho=-i (A*B)^<\f$.
* > The objects A and B are of the class type 'GG'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* [int] given time-step
* @param rho
* > [std::complex] Matrix to which the result of the convolution is given
* @param A
* > [GG] contour Green's function
* @param B
* > [GG] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time interval
*/
template <typename T, class GG>
void convolution_density_matrix(int tstp, std::complex<T> *rho, GG &A, GG &B,
                                integration::Integrator<T> &I, T beta, T h) {
    convolution_density_matrix<T, GG>(tstp, rho, A, A, B, B, I, beta, h);
}



/** \brief <b> Returns the result of the contour convolution for a density matrix</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calls convolution routines to compute contour convolution \f$\rho=-i (A*B)^<\f$.
* > The objects A and B are of the class type 'GG'. Works for a scalar and square matrices.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* [int] given time-step
* @param rho
* > [std::complex] Matrix to which the result of the convolution is given
* @param A
* > [GG] contour Green's function
* @param B
* > [GG] contour Green's function
* @param beta
* > inversed temperature
* @param h
* > time interval
* @param SolveOrder
* > [int] integrator ordder
*/
template <typename T, class GG>
void convolution_density_matrix(int tstp, cdmatrix &rho, GG &A, GG &B,
                                T beta, T h, int SolveOrder) {
    int size1 = A.size1();
    int element_size = size1;
    std::complex<T> *rho_ptr;
    rho_ptr = new std::complex<T>[element_size];

    convolution_density_matrix<T, GG>(tstp, rho_ptr, A, A, B, B, integration::I<T>(SolveOrder), beta, h);
    map_ptr2matrix<T>(size1, size1, rho_ptr, rho);

    delete[] rho_ptr;
}

/// @private
template <typename T, class GG>
void convolution_les_timediag(int tstp, cdmatrix &Cles, GG &A, GG &B,
			      integration::Integrator<T> &I, T beta, T h){
  int size1=A.size1();
  std::complex<T> *Cles_ptr = new std::complex<T>[size1*size1];
  assert(Cles.rows() == size1);
  assert(Cles.cols() == size1);

  convolution_density_matrix<T, GG>(tstp, Cles_ptr, A, A, B, B, I, beta, h);

  for(int n1=0; n1 < size1; n1++){
    for(int n2=0; n2 < size1; n2++){
      Cles(n1,n2) = std::complex<T>(0, 1.0) * Cles_ptr[n1 * size1 + n2];
    }
  }

  delete Cles_ptr;

}

/// @private
template <typename T, class GG>
void convolution_density_matrix(int tstp, cdmatrix &Cles, GG &A, GG &B,
                  integration::Integrator<T> &I, T beta, T h){
  int size1=A.size1();
  std::complex<T> *Cles_ptr = new std::complex<T>[size1*size1];
  assert(Cles.rows() == size1);
  assert(Cles.cols() == size1);

  convolution_density_matrix<T, GG>(tstp, Cles_ptr, A, A, B, B, I, beta, h);
  for(int n1=0; n1 < size1; n1++){
    for(int n2=0; n2 < size1; n2++){
      Cles(n1,n2) = std::complex<T>(0, 1.0) * Cles_ptr[n1 * size1 + n2];
    }
  }
  delete Cles_ptr;

}

/* /////////////////////////////////////////////////////////////////////////////////////////
// INCREMENETAL
///////////////////////////////////////////////////////////////////////////////////////// */

/** \brief <b> Adds a Matsubara convolution of two matrices and a function with a given weight to a matrix </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Performs the operation \f$C^M \rightarrow C^M + \alpha A^M*f0*B^M\f$, where 'C', 'A',and 'B' are objects of the type 'GG',
 * > 'f0' is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis, and \f$\alpha\f$ is a complex weight.
 * >
 * > Note: parameter 'mask' allows to exclute some time steps from calculations. It can be useful for OMP calculations.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param mask
 * > [std::vector] parameter allowing to exclude some time steps
 * @param alpha
 * > [CPLX] The weight in front of \f$A*f0*B\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param *f0
 * > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
 * @param B
 * > [GG] contour Green's function
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 */
#define CPLX std::complex<T>
template <typename T, class GG, int SIZE1>
void incr_convolution_mat(std::vector<bool> &mask, CPLX alpha, GG &C, GG &A, CPLX *f0, GG &B,
                          integration::Integrator<T> &I, T beta) {
    int ntau = A.ntau();
    int size1 = A.size1();
    int sc = size1 * size1;
    bool func = (f0 == NULL ? false : true);
    CPLX adtau = alpha * beta * (1.0 / ntau);
    {
        // CONVOLUTION OF MATSUBARA GREENFUNCTIONS
        int m, sfb = sc;
        CPLX *ctemp = new CPLX[sc], *bmat1 = 0;
        CPLX *btemp;
        if (func) {
            bmat1 = new CPLX[(ntau + 1) * sfb];
            for (m = 0; m <= ntau; m++) {
                element_mult<T, SIZE1>(size1, bmat1 + m * sfb, f0, B.matptr(m));
            }
            btemp = bmat1;
        } else {
            btemp = B.matptr(0);
        }
        // t1=omp_get_wtime();
        for (m = 0; m <= ntau; m++) {
            if (mask[m]) {
                // compute cmat(m*dtau) = int_0^beta dx amat(tau-x) f0 b(x)
                matsubara_integral_1<T, SIZE1>(size1, m, ntau, ctemp, A.matptr(0), btemp, I,
                                               A.sig());
                element_incr<T, SIZE1>(size1, C.matptr(m), adtau, ctemp);
            }
        }
        // t2=omp_get_wtime();
        if (func)
            delete[] bmat1;
        delete[] ctemp;
    }
    // std::cout << "tid " << omp_get_thread_num() << " parallel time " << t2-t1  <<
    // std::endl;
    return;
}
/** \brief <b> Adds a retarded convolution of two matrices and a function with a given weight to a matrix at a given time step </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Performs the operation \f$C^R \rightarrow C^R + \alpha (A*ft*B)^R\f$, where 'C', 'A',and 'B' are objects of the type 'GG',
 * > 'ft' is a pointer to \f$F(t)\f$ on the real axis, and \f$\alpha\f$ is a complex weight. The operation is performed at given time step `tstp`.
 * >
 * > Note: parameter 'mask' allows to exclute some time steps from calculations. It can be useful for OMP calculations.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > [int] given time-step
 * @param mask
 * > [std::vector] parameter allowing to exclude some time steps
 * @param alpha
 * > [CPLX] The weight in front of \f$A*ft*B\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void incr_convolution_ret(int tstp, std::vector<bool> &mask, CPLX alpha, GG &C, GG &A,
                          GG &Acc, CPLX *ft, GG &B, GG &Bcc, integration::Integrator<T> &I,
                          T h) {
    int size1 = A.size1();
    int sc = size1 * size1;
    bool func = (ft == NULL ? false : true);
    CPLX adt = alpha * h;
    int SolveOrder = I.get_k();
    int n1 = (tstp >= SolveOrder ? tstp : SolveOrder);
    {
        // CONVOLUTION OF RET SECTION
        int j, n, l;
        int saf = size1 * size1;
        CPLX *ctemp0 = new CPLX[sc];
        CPLX *aret;
        CPLX *btmp = new CPLX[sc];
        CPLX *aret1 = 0;
        T wt;
        // aret[j]=Aret(tstp,j)*f(j)   j=0 ... n1
        {
            if (func) {
                CPLX *atemp = new CPLX[sc];
                aret1 = new CPLX[(n1 + 1) * saf];
                aret = aret1;
                for (j = 0; j <= tstp; j++) {
                    element_mult<T, SIZE1>(size1, aret1 + saf * j, A.retptr(tstp, j),
                                           ft + j * sc);
                }
                for (j = tstp + 1; j <= n1; j++) {
                    element_minusconj<T, SIZE1>(size1, atemp, Acc.retptr(j, tstp));
                    element_mult<T, SIZE1>(size1, aret1 + saf * j, atemp, ft + j * sc);
                }
                delete[] atemp;
            } else {
                if (n1 == tstp) {
                    aret = A.retptr(tstp, 0);
                } else {
                    aret1 = new CPLX[(n1 + 1) * saf];
                    aret = aret1;
                    for (j = 0; j <= tstp; j++)
                        element_set<T, SIZE1>(size1, aret1 + sc * j, A.retptr(tstp, j));
                    for (j = tstp + 1; j <= n1; j++)
                        element_minusconj<T, SIZE1>(size1, aret1 + sc * j,
                                                    Acc.retptr(j, tstp));
                }
            }
        }
        for (n = 0; n <= tstp; n++) {
            if (mask[n]) {
                // int_n^tstp dj Aret(tstp,j)Bret(j,n)  ==>  Cret(tstp,n)
                int n2 = tstp - n;
                element_set_zero<T, SIZE1>(size1, ctemp0);
                if (tstp < SolveOrder) {
                    for (j = 0; j <= SolveOrder; j++) {
                        wt = I.poly_integration(n, tstp, j);
                        if (j >= n) {
                            element_incr<T, SIZE1>(size1, ctemp0, wt, aret + j * saf,
                                                   B.retptr(j, n));
                        } else {
                            element_minusconj<T, SIZE1>(size1, btmp, Bcc.retptr(n, j));
                            element_incr<T, SIZE1>(size1, ctemp0, wt, aret + j * saf, btmp);
                        }
                    }
                } else if (n2 < SolveOrder) {
                    for (l = 0; l <= SolveOrder; l++) {
                        j = tstp - l;
                        wt = I.gregory_weights(n2, l);
                        if (j >= n) {
                            element_incr<T, SIZE1>(size1, ctemp0, wt, aret + j * saf,
                                                   B.retptr(j, n));
                        } else {
                            element_minusconj<T, SIZE1>(size1, btmp, Bcc.retptr(n, j));
                            element_incr<T, SIZE1>(size1, ctemp0, wt, aret + j * saf, btmp);
                        }
                    }
                } else if (n2 <= 2 * SolveOrder + 2) {
                    for (l = 0; l <= n2; l++) {
                        j = n + l;
                        wt = I.gregory_weights(n2, l);
                        element_incr<T, SIZE1>(size1, ctemp0, wt, aret + j * saf,
                                               B.retptr(j, n));
                    }
                } else {
                    for (j = n; j <= n + SolveOrder; j++) {
                        wt = I.gregory_omega(j - n);
                        element_incr<T, SIZE1>(size1, ctemp0, wt, aret + j * saf,
                                               B.retptr(j, n));
                    }
                    for (j = n + SolveOrder + 1; j < tstp - SolveOrder; j++) {
                        element_incr<T, SIZE1>(size1, ctemp0, aret + j * saf,
                                               B.retptr(j, n));
                    }
                    for (l = 0; l <= SolveOrder; l++) {
                        j = tstp - l;
                        wt = I.gregory_omega(l);
                        element_incr<T, SIZE1>(size1, ctemp0, wt, aret + j * saf,
                                               B.retptr(j, n));
                    }
                }
                element_incr<T, SIZE1>(size1, C.retptr(tstp, n), adt, ctemp0);
            }
        }
        delete[] ctemp0;
        delete[] btmp;
        if (aret1 != 0)
            delete[] aret1;
    }
    return;
}
/** \brief <b> Adds a left-mixing convolution of two matrices and a function with a given weight to a matrix at a given time step </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Performs the operation \f$C^{\rceil} \rightarrow C^{\rceil} + \alpha (A*ft*B)^{\rceil}\f$, where 'C', 'A',and 'B' are objects of the type 'GG',
 * > 'f0' is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis,'ft' is a pointer to \f$F(t)\f$ on the real axis, and \f$\alpha\f$ is a complex weight. The operation is performed at given time step `tstp`.
 * >
 * > Note: parameter 'mask' allows to exclute some time steps from calculations. It can be useful for OMP calculations.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > [int] given time-step
 * @param mask
 * > [std::vector] parameter allowing to exclude some time steps
 * @param alpha
 * > [CPLX] The weight in front of \f$A*ft*B\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *f0
 * > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void incr_convolution_tv(int tstp, std::vector<bool> &mask, CPLX alpha, GG &C, GG &A,
                         GG &Acc, CPLX *f0, CPLX *ft, GG &B, GG &Bcc,
                         integration::Integrator<T> &I, T beta, T h) {
    int ntau = A.ntau();
    int size1 = A.size1();
    int sc = size1 * size1;
    bool func = (ft == NULL ? false : true);
    CPLX adt = alpha * h;
    CPLX adtau = alpha * beta * (1.0 / ntau);
    int SolveOrder = I.get_k();
    int n1 = (tstp >= SolveOrder ? tstp : SolveOrder);
    {
        // CONVOLUTION OF TV SECTION
        int j, m, n, saf = size1 * size1, sfb = size1 * size1;
        CPLX *ctemp1 = new CPLX[sc];
        CPLX *ctemp2 = new CPLX[sc];
        CPLX *bmat;
        CPLX *aret;
        CPLX *bmat1 = 0;
        CPLX *aret1 = 0;
        {
            // aret[j]=Aret(tstp,j)*f(j)   j=0 ... n1
            if (func) {
                CPLX *atemp = new CPLX[sc];
                aret1 = new CPLX[(n1 + 1) * saf];
                aret = aret1;
                for (j = 0; j <= tstp; j++) {
                    element_mult<T, SIZE1>(size1, aret1 + saf * j, A.retptr(tstp, j),
                                           ft + j * sc);
                }
                for (j = tstp + 1; j <= n1; j++) {
                    element_minusconj<T, SIZE1>(size1, atemp, Acc.retptr(j, tstp));
                    element_mult<T, SIZE1>(size1, aret1 + saf * j, atemp, ft + j * sc);
                }
                delete[] atemp;
            } else {
                if (n1 == tstp) {
                    aret = A.retptr(tstp, 0);
                } else {
                    aret1 = new CPLX[(n1 + 1) * saf];
                    aret = aret1;
                    for (j = 0; j <= tstp; j++)
                        element_set<T, SIZE1>(size1, aret1 + sc * j, A.retptr(tstp, j));
                    for (j = tstp + 1; j <= n1; j++)
                        element_minusconj<T, SIZE1>(size1, aret1 + sc * j,
                                                    Acc.retptr(j, tstp));
                }
            }
            // bmat[m]=f0*Bmat(m)
            if (func) {
                bmat1 = new CPLX[(ntau + 1) * sfb];
                bmat = bmat1;
                for (m = 0; m <= ntau; m++)
                    element_mult<T, SIZE1>(size1, bmat1 + m * sfb, f0, B.matptr(m));
            } else {
                bmat = B.matptr(0);
            }
        }
        for (m = 0; m <= ntau; m++) {
            if (mask[m]) {
                // CONTRIBUTION FROM Atv * Bmat : very similar to computing the matsubara
                // convolution
                matsubara_integral_2<T, SIZE1>(size1, m, ntau, ctemp1, A.tvptr(tstp, 0),
                                               bmat, I, B.sig());
                element_smul<T, SIZE1>(size1, ctemp1, adtau);
                // CONTRIBUTION FROM Aret * Btv
                element_set_zero<T, SIZE1>(size1, ctemp2);
                if (tstp <= 2 * SolveOrder + 2) {
                    for (n = 0; n <= n1; n++) {
                        element_incr<T, SIZE1>(size1, ctemp2, I.gregory_weights(tstp, n),
                                               aret + n * saf, B.tvptr(n, m));
                    }
                } else {
                    for (n = 0; n <= SolveOrder; n++) {
                        element_incr<T, SIZE1>(size1, ctemp2, I.gregory_omega(n),
                                               aret + n * saf, B.tvptr(n, m));
                    }
                    for (n = SolveOrder + 1; n < tstp - SolveOrder; n++) {
                        element_incr<T, SIZE1>(size1, ctemp2, aret + n * saf, B.tvptr(n, m));
                    }
                    for (n = tstp - SolveOrder; n <= tstp; n++) {
                        element_incr<T, SIZE1>(size1, ctemp2, I.gregory_omega(tstp - n),
                                               aret + n * saf, B.tvptr(n, m));
                    }
                }
                element_smul<T, SIZE1>(size1, ctemp2, adt);
                element_incr<T, SIZE1>(size1, C.tvptr(tstp, m), ctemp1);
                element_incr<T, SIZE1>(size1, C.tvptr(tstp, m), ctemp2);
            }
        }
        if (bmat1 != 0)
            delete[] bmat1;
        if (aret1 != 0)
            delete[] aret1;
        delete[] ctemp1;
        delete[] ctemp2;
    }
    return;
}
/** \brief <b> Adds a lesser convolution of two matrices and a function with a given weight to a matrix at a given time step </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Performs the operation \f$C^{<} \rightarrow C^{<} + \alpha (A*ft*B)^{<}\f$, where 'C', 'A',and 'B' are objects of the type 'GG',
 * > 'f0' is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis,'ft' is a pointer to \f$F(t)\f$ on the real axis, and \f$\alpha\f$ is a complex weight. The operation is performed at given time step `tstp`.
 * >
 * > Note: parameter 'mask' allows to exclute some time steps from calculations. It can be useful for OMP calculations.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > [int] given time-step
 * @param mask
 * > [std::vector] parameter allowing to exclude some time steps
 * @param alpha
 * > [CPLX] The weight in front of \f$A*ft*B\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *f0
 * > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void incr_convolution_les(int tstp, std::vector<bool> &mask, CPLX alpha, GG &C, GG &A,
                          GG &Acc, CPLX *f0, CPLX *ft, GG &B, GG &Bcc,
                          integration::Integrator<T> &I, T beta, T h) {
    int ntau = A.ntau();
    int size1 = A.size1();
    int sc = size1 * size1;
    bool func = (ft == NULL ? false : true);
    CPLX adt = alpha * h;
    CPLX adtau = alpha * beta * (1.0 / ntau);
    int SolveOrder = I.get_k();
    int n1 = (tstp >= SolveOrder ? tstp : SolveOrder);
    {
        // CONVOLUTION OF LES SECTION
        int m, j, n, sfb = size1 * size1;
        CPLX *badv = new CPLX[(n1 + 1) * sfb];  // badv[j]=f(j)*Badv(j,tstp) j=0 ... n1
        CPLX *bles = new CPLX[(n1 + 1) * sfb];  // bles[j]=f(j)*Bles(j,tstp)   j=0 ... n1
        CPLX *bvt = new CPLX[(ntau + 1) * sfb]; // bvt[m]=f0*Bvt(m,tstp)    m=0 ... ntau
        CPLX *ctemp1 = new CPLX[sc];
        CPLX *ctemp2 = new CPLX[sc];
        CPLX *ctemp3 = new CPLX[sc];
        CPLX *atemp = new CPLX[sc];
        {
            CPLX *btemp = new CPLX[sc];
            if (func) {
                for (j = 0; j <= tstp; j++) {
                    element_conj<T, SIZE1>(size1, btemp, Bcc.retptr(tstp, j));
                    element_mult<T, SIZE1>(size1, badv + sfb * j, ft + sc * j, btemp);
                    element_mult<T, SIZE1>(size1, bles + sfb * j, ft + sc * j,
                                           B.lesptr(j, tstp));
                }
                for (j = tstp + 1; j <= n1; j++) {
                    element_mult<T, SIZE1>(size1, badv + sfb * j, ft + sc * j,
                                           B.retptr(j, tstp));
                    element_smul<T, SIZE1>(size1, badv + sfb * j, CPLX(-1, 0));
                    element_minusconj<T, SIZE1>(size1, btemp, Bcc.lesptr(tstp, j));
                    element_mult<T, SIZE1>(size1, bles + sfb * j, ft + sc * j, btemp);
                }
                if (B.sig() == -1) {
                    for (m = 0; m <= ntau; m++) {
                        element_conj<T, SIZE1>(size1, btemp, Bcc.tvptr(tstp, ntau - m));
                        element_mult<T, SIZE1>(size1, bvt + sfb * m, f0, btemp);
                    }
                } else {
                    // BOSE !!
                    for (m = 0; m <= ntau; m++) {
                        element_minusconj<T, SIZE1>(size1, btemp, Bcc.tvptr(tstp, ntau - m));
                        element_mult<T, SIZE1>(size1, bvt + sfb * m, f0, btemp);
                    }
                }
            } else {
                for (j = 0; j <= tstp; j++) {
                    element_conj<T, SIZE1>(size1, badv + sc * j, Bcc.retptr(tstp, j));
                    element_set<T, SIZE1>(size1, bles + sc * j, B.lesptr(j, tstp));
                }
                for (j = tstp + 1; j <= n1; j++) {
                    element_set<T, SIZE1>(size1, badv + sc * j, B.retptr(j, tstp));
                    element_smul<T, SIZE1>(size1, badv + sc * j, CPLX(-1, 0));
                    element_minusconj<T, SIZE1>(size1, bles + sc * j, Bcc.lesptr(tstp, j));
                }
                if (B.sig() == -1) {
                    for (m = 0; m <= ntau; m++) {
                        element_conj<T, SIZE1>(size1, bvt + sc * m,
                                               Bcc.tvptr(tstp, ntau - m));
                    }
                } else {
                    // BOSE !!
                    for (m = 0; m <= ntau; m++) {
                        element_minusconj<T, SIZE1>(size1, bvt + sc * m,
                                                    Bcc.tvptr(tstp, ntau - m));
                    }
                }
            }
            delete[] btemp;
        }
        for (n = 0; n <= tstp; n++) {
            if (mask[n]) {
                T wt;
                /* int_0^n dj Aret(n,j)Bles(j,tstp)  ==>  Cles(n,tstp)  */
                {
                    int nup = (n > SolveOrder ? n : SolveOrder);
                    element_set_zero<T, SIZE1>(size1, ctemp1);
                    if (n <= 2 * SolveOrder + 2) {
                        for (j = 0; j <= nup; j++) {
                            wt = I.gregory_weights(n, j);
                            if (j <= n) {
                                element_incr<T, SIZE1>(size1, ctemp1, wt, A.retptr(n, j),
                                                       bles + j * sfb);
                            } else {
                                element_minusconj<T, SIZE1>(size1, atemp, Acc.retptr(j, n));
                                element_incr<T, SIZE1>(size1, ctemp1, wt, atemp,
                                                       bles + j * sfb);
                            }
                        }
                    } else {
                        for (j = 0; j <= SolveOrder; j++) {
                            wt = I.gregory_omega(j);
                            element_incr<T, SIZE1>(size1, ctemp1, wt, A.retptr(n, j),
                                                   bles + j * sfb);
                        }
                        for (j = SolveOrder + 1; j < n - SolveOrder; j++) {
                            element_incr<T, SIZE1>(size1, ctemp1, A.retptr(n, j),
                                                   bles + j * sfb);
                        }
                        for (j = n - SolveOrder; j <= n; j++) {
                            wt = I.gregory_omega(n - j);
                            element_incr<T, SIZE1>(size1, ctemp1, wt, A.retptr(n, j),
                                                   bles + j * sfb);
                        }
                    }
                }
                /* int_0^beta dj Atv(n,j)Bvt(j,tstp) ==>  Cles(n,tstp) */
                {
                    element_set_zero<T, SIZE1>(size1, ctemp2);
                    for (m = 0; m <= SolveOrder; m++) {
                        wt = I.gregory_omega(m);
                        element_incr<T, SIZE1>(size1, ctemp2, wt, A.tvptr(n, m),
                                               bvt + sfb * m);
                    }
                    for (m = SolveOrder + 1; m < ntau - SolveOrder; m++) {
                        element_incr<T, SIZE1>(size1, ctemp2, A.tvptr(n, m), bvt + sfb * m);
                    }
                    for (m = ntau - SolveOrder; m <= ntau; m++) {
                        wt = I.gregory_omega(ntau - m);
                        element_incr<T, SIZE1>(size1, ctemp2, wt, A.tvptr(n, m),
                                               bvt + sfb * m);
                    }
                }
                /* int_0^tstp dj Ales(n,j)Badv(j,tstp)  ==>  Cles(n,tstp) */
                {
                    element_set_zero<T, SIZE1>(size1, ctemp3);
                    if (tstp <= 2 * SolveOrder + 2) {
                        for (j = 0; j <= n; j++) {
                            wt = I.gregory_weights(tstp, j);
                            element_minusconj<T, SIZE1>(size1, atemp, Acc.lesptr(j, n));
                            element_incr<T, SIZE1>(size1, ctemp3, wt, atemp, badv + sfb * j);
                        }
                        for (j = n + 1; j <= n1; j++) {
                            wt = I.gregory_weights(tstp, j);
                            element_incr<T, SIZE1>(size1, ctemp3, wt, A.lesptr(n, j),
                                                   badv + sfb * j);
                        }
                    } else {
                        for (j = 0; j <= SolveOrder; j++) {
                            wt = I.gregory_omega(j);
                            if (j < n) {
                                element_minusconj<T, SIZE1>(size1, atemp, Acc.lesptr(j, n));
                                element_incr<T, SIZE1>(size1, ctemp3, wt, atemp,
                                                       badv + sfb * j);
                            } else {
                                element_incr<T, SIZE1>(size1, ctemp3, wt, A.lesptr(n, j),
                                                       badv + sfb * j);
                            }
                        }
                        if (n <= SolveOrder) {
                            for (j = SolveOrder + 1; j < tstp - SolveOrder; j++) {
                                element_incr<T, SIZE1>(size1, ctemp3, A.lesptr(n, j),
                                                       badv + sfb * j);
                            }
                        } else if (n < tstp - SolveOrder) {
                            for (j = SolveOrder + 1; j < n; j++) {
                                element_minusconj<T, SIZE1>(size1, atemp, Acc.lesptr(j, n));
                                element_incr<T, SIZE1>(size1, ctemp3, atemp, badv + sfb * j);
                            }
                            for (j = n; j < tstp - SolveOrder; j++) {
                                element_incr<T, SIZE1>(size1, ctemp3, A.lesptr(n, j),
                                                       badv + sfb * j);
                            }
                        } else {
                            for (j = SolveOrder + 1; j < tstp - SolveOrder; j++) {
                                element_minusconj<T, SIZE1>(size1, atemp, Acc.lesptr(j, n));
                                element_incr<T, SIZE1>(size1, ctemp3, atemp, badv + sfb * j);
                            }
                        }
                        for (j = tstp - SolveOrder; j <= tstp; j++) {
                            wt = I.gregory_omega(tstp - j);
                            if (j < n) {
                                element_minusconj<T, SIZE1>(size1, atemp, Acc.lesptr(j, n));
                                element_incr<T, SIZE1>(size1, ctemp3, wt, atemp,
                                                       badv + sfb * j);
                            } else {
                                element_incr<T, SIZE1>(size1, ctemp3, wt, A.lesptr(n, j),
                                                       badv + sfb * j);
                            }
                        }
                    }
                }
                element_incr<T, SIZE1>(size1, C.lesptr(n, tstp), adt, ctemp1);
                element_incr<T, SIZE1>(size1, C.lesptr(n, tstp), adt, ctemp3);
                element_incr<T, SIZE1>(size1, C.lesptr(n, tstp), adtau * CPLX(0.0, -1.0),
                                       ctemp2);
            }
        }
        delete[] ctemp1;
        delete[] ctemp2;
        delete[] ctemp3;
        delete[] atemp;
        delete[] bles;
        delete[] badv;
        delete[] bvt;
    }
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// SINGLE-PROCESSOR ROUTINES ... function calls like for old version
// with function object

/** \brief <b> Adds a convolution of two matrices and a function with a given weight to a matrix at a given time step </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Performs the operation \f$C \rightarrow C + \alpha A*ft*B\f$, where 'C', 'A',and 'B' are objects of the type 'GG',
 * > 'f0' is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis,'ft' is a pointer to \f$F(t)\f$ on the real axis, and \f$\alpha\f$ is a complex weight. The operation is performed at given time step `tstp`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > [int] given time-step
 * @param alpha
 * > [CPLX] The weight in front of \f$A*f0*B\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *f0
 * > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void incr_convolution(int tstp, CPLX alpha, GG &C, GG &A, GG &Acc, CPLX *f0, CPLX *ft, GG &B,
                      GG &Bcc, integration::Integrator<T> &I, T beta, T h) {
    int ntau = A.ntau();
    // this function is still not on top level, so no asserts!
    if (tstp == -1) {
        std::vector<bool> mask(ntau + 1, true);
        incr_convolution_mat<T, GG, SIZE1>(mask, alpha, C, A, f0, B, I, beta);
    } else if (tstp >= 0) {
        std::vector<bool> mask_ret(tstp + 1, true);
        std::vector<bool> mask_tv(ntau + 1, true);
        std::vector<bool> mask_les(tstp + 1, true);
        incr_convolution_ret<T, GG, SIZE1>(tstp, mask_ret, alpha, C, A, Acc, ft, B, Bcc, I,
                                           h);
        incr_convolution_tv<T, GG, SIZE1>(tstp, mask_tv, alpha, C, A, Acc, f0, ft, B, Bcc, I,
                                          beta, h);
        incr_convolution_les<T, GG, SIZE1>(tstp, mask_les, alpha, C, A, Acc, f0, ft, B, Bcc,
                                           I, beta, h);
    }
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// NEW VERSIONS
/// @private
template <typename T>
void convolution_timestep_new(int tstp, herm_matrix<T> &C, herm_matrix<T> &A,
                              herm_matrix<T> &Acc, function<T> &ft, herm_matrix<T> &B,
                              herm_matrix<T> &Bcc, integration::Integrator<T> &I, T beta,
                              T h) {
    int SolveOrder = I.k();
    int ntmin = (tstp == -1 || tstp > SolveOrder ? tstp : SolveOrder);
    if (tstp < -1)
        return;
    int size1 = A.size1();
    std::complex<T> *fttemp;

    assert(C.size1() == size1);
    assert(ft.size1() == size1);
    assert(B.size1() == size1);
    assert(Bcc.size1() == size1);
    assert(Acc.size1() == size1);
    assert(SolveOrder > 0 && SolveOrder <= 5);
    assert(SolveOrder <= C.ntau());
    assert(ntmin <= C.nt());
    assert(ntmin <= A.nt());
    assert(ntmin <= Acc.nt());
    assert(ntmin <= B.nt());
    assert(ntmin <= Bcc.nt());
    assert(ntmin <= ft.nt());
    assert(C.ntau() == A.ntau());
    assert(C.ntau() == Acc.ntau());
    assert(C.ntau() == B.ntau());
    assert(C.ntau() == Bcc.ntau());

    C.set_timestep_zero(tstp);
    fttemp = (tstp == -1 ? ft.ptr(-1) : ft.ptr(0));

    switch (size1) {
    case 1:
        incr_convolution<T, herm_matrix<T>, 1>(tstp, CPLX(1, 0), C, A, Acc, ft.ptr(-1),
                                               fttemp, B, Bcc, integration::I<T>(SolveOrder), beta,
                                               h);
        break;
    case 2:
        incr_convolution<T, herm_matrix<T>, 2>(tstp, CPLX(1, 0), C, A, Acc, ft.ptr(-1),
                                               fttemp, B, Bcc, integration::I<T>(SolveOrder), beta,
                                               h);
        break;
    case 3:
        incr_convolution<T, herm_matrix<T>, 3>(tstp, CPLX(1, 0), C, A, Acc, ft.ptr(-1),
                                               fttemp, B, Bcc, integration::I<T>(SolveOrder), beta,
                                               h);
        break;
    case 4:
        incr_convolution<T, herm_matrix<T>, 4>(tstp, CPLX(1, 0), C, A, Acc, ft.ptr(-1),
                                               fttemp, B, Bcc, integration::I<T>(SolveOrder), beta,
                                               h);
        break;
    case 5:
        incr_convolution<T, herm_matrix<T>, 5>(tstp, CPLX(1, 0), C, A, Acc, ft.ptr(-1),
                                               fttemp, B, Bcc, integration::I<T>(SolveOrder), beta,
                                               h);
        break;
    case 6:
        incr_convolution<T, herm_matrix<T>, 6>(tstp, CPLX(1, 0), C, A, Acc, ft.ptr(-1),
                                               fttemp, B, Bcc, integration::I<T>(SolveOrder), beta,
                                               h);
        break;
    case 8:
        incr_convolution<T, herm_matrix<T>, 8>(tstp, CPLX(1, 0), C, A, Acc, ft.ptr(-1),
                                               fttemp, B, Bcc, integration::I<T>(SolveOrder), beta,
                                               h);
        break;
    default:
        incr_convolution<T, herm_matrix<T>, LARGESIZE>(tstp, CPLX(1, 0), C, A, Acc,
                                                       ft.ptr(-1), fttemp, B, Bcc,
                                                       integration::I<T>(SolveOrder), beta, h);
        break;
    }
}
/// @private
template <typename T>
void convolution_timestep_new(int n, herm_matrix<T> &C, herm_matrix<T> &A, function<T> &ft,
                              herm_matrix<T> &B, integration::Integrator<T> &I, T beta,
                              T h) {
    convolution_timestep_new<T>(n, C, A, A, ft, B, B, I, beta, h);
}
/// @private
template <typename T>
void convolution_matsubara_new(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &B,
                               integration::Integrator<T> &I, T beta) {
    convolution_timestep_new<T>(-1, C, A, A, B, B, I, beta, 0.0);
}
/// @private
template <typename T>
void convolution_matsubara_new(herm_matrix<T> &C, herm_matrix<T> &A, function<T> &ft,
                               herm_matrix<T> &B, integration::Integrator<T> &I, T beta) {
    convolution_timestep_new<T>(-1, C, A, A, ft, B, B, I, beta, 0.0);
}
/// @private
template <typename T>
void convolution_new(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                     function<T> &ft, herm_matrix<T> &B, herm_matrix<T> &Bcc,
                     integration::Integrator<T> &I, T beta, T h) {
    int tstp;
    for (tstp = -1; tstp <= C.nt(); tstp++)
        convolution_timestep_new<T>(tstp, C, A, Acc, ft, B, Bcc, I, beta, h);
}
/// @private
template <typename T>
void convolution_timestep_new(int tstp, herm_matrix<T> &C, herm_matrix<T> &A,
                              herm_matrix<T> &Acc, herm_matrix<T> &B, herm_matrix<T> &Bcc,
                              integration::Integrator<T> &I, T beta, T h) {
    int SolveOrder = I.k();
    int ntmin = (tstp == -1 || tstp > SolveOrder ? tstp : SolveOrder);
    if (tstp < -1)
        return;
    int size1 = A.size1();

    assert(C.size1() == size1);
    assert(B.size1() == size1);
    assert(Bcc.size1() == size1);
    assert(Acc.size1() == size1);
    assert(SolveOrder > 0 && SolveOrder <= 5);
    assert(SolveOrder <= C.ntau());
    assert(ntmin <= C.nt());
    assert(ntmin <= A.nt());
    assert(ntmin <= Acc.nt());
    assert(ntmin <= B.nt());
    assert(ntmin <= Bcc.nt());
    assert(C.ntau() == A.ntau());
    assert(C.ntau() == Acc.ntau());
    assert(C.ntau() == B.ntau());
    assert(C.ntau() == Bcc.ntau());

    C.set_timestep_zero(tstp);
    switch (size1) {
    case 1:
        incr_convolution<T, herm_matrix<T>, 1>(tstp, CPLX(1, 0), C, A, Acc, NULL, NULL, B,
                                               Bcc, integration::I<T>(SolveOrder), beta, h);
        break;
    case 2:
        incr_convolution<T, herm_matrix<T>, 2>(tstp, CPLX(1, 0), C, A, Acc, NULL, NULL, B,
                                               Bcc, integration::I<T>(SolveOrder), beta, h);
        break;
    case 3:
        incr_convolution<T, herm_matrix<T>, 3>(tstp, CPLX(1, 0), C, A, Acc, NULL, NULL, B,
                                               Bcc, integration::I<T>(SolveOrder), beta, h);
        break;
    case 4:
        incr_convolution<T, herm_matrix<T>, 4>(tstp, CPLX(1, 0), C, A, Acc, NULL, NULL, B,
                                               Bcc, integration::I<T>(SolveOrder), beta, h);
        break;
    case 5:
        incr_convolution<T, herm_matrix<T>, 5>(tstp, CPLX(1, 0), C, A, Acc, NULL, NULL, B,
                                               Bcc, integration::I<T>(SolveOrder), beta, h);
        break;
    case 6:
        incr_convolution<T, herm_matrix<T>, 6>(tstp, CPLX(1, 0), C, A, Acc, NULL, NULL, B,
                                               Bcc, integration::I<T>(SolveOrder), beta, h);
        break;
    case 8:
        incr_convolution<T, herm_matrix<T>, 8>(tstp, CPLX(1, 0), C, A, Acc, NULL, NULL, B,
                                               Bcc, integration::I<T>(SolveOrder), beta, h);
        break;
    default:
        incr_convolution<T, herm_matrix<T>, LARGESIZE>(
            tstp, CPLX(1, 0), C, A, Acc, NULL, NULL, B, Bcc, integration::I<T>(SolveOrder), beta, h);
        break;
    }
}
/// @private
template <typename T>
void convolution_timestep_new(int n, herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &B,
                              integration::Integrator<T> &I, T beta, T h) {
    convolution_timestep_new<T>(n, C, A, A, B, B, I, beta, h);
}
/// @private
template <typename T>
void convolution_new(herm_matrix<T> &C, herm_matrix<T> &A, herm_matrix<T> &Acc,
                     herm_matrix<T> &B, herm_matrix<T> &Bcc, integration::Integrator<T> &I,
                     T beta, T h) {
    int tstp;
    for (tstp = -1; tstp <= C.nt(); tstp++)
        convolution_timestep_new<T>(tstp, C, A, Acc, B, Bcc, I, beta, h);
}
#undef CPLX
//////////////////////////////////////////////////////////////////////////////////////////////////////
// OPEN-MP paralellized routines
#if CNTR_USE_OMP == 1

#define CPLX std::complex<T>
/** \brief <b> Adds a convolution of two matrices and a function with a given weight to a matrix at a given time step </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Performs the operation \f$C \rightarrow C + \alpha A*ft*B\f$, where 'C', 'A',and 'B' are objects of the type 'GG',
 * > 'f0' is a pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis,'ft' is a pointer to \f$F(t)\f$ on the real axis, 
 * > and \f$\alpha\f$ is a complex weight. The operation is performed at given time step `tstp`. `openMP` parallelized version.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param omp_num_threads
 * > [int] number of `openMP` threads
 * @param tstp
 * > [int] given time-step
 * @param alpha
 * > [CPLX] The weight in front of \f$A*f0*B\f$
 * @param C
 * > [GG] contour Green's function
 * @param A
 * > [GG] contour Green's function
 * @param Acc
 * > [GG] complex conjugate to A
 * @param *f0
 * > [std::complex] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix)
 * @param *ft
 * > [std::complex] Pointer to \f$F(t)\f$ on the real axis (complex function).
 * @param B
 * > [GG] contour Green's function
 * @param Bcc
 * > [GG] complex conjugate to B
 * @param I
* > [Integrator] integrator class
 * @param beta
 * > inversed temperature
 * @param h
 * > time step interval
 */
template <typename T, class GG, int SIZE1>
void incr_convolution_omp(int omp_num_threads, int tstp, CPLX alpha, GG &C, GG &A, GG &Acc,
                          CPLX *f0, CPLX *ft, GG &B, GG &Bcc, integration::Integrator<T> &I,
                          T beta, T h) {
#pragma omp parallel num_threads(omp_num_threads)
    {
        int nomp = omp_get_num_threads();
        int tid = omp_get_thread_num();
        int ntau = A.ntau(), i;
        // this function is still not on top level, so no asserts!
        if (tstp == -1) {
            std::vector<bool> mask(ntau + 1, false);
            for (i = 0; i <= ntau; i++)
                if (i % nomp == tid)
                    mask[i] = true;
            incr_convolution_mat<T, GG, SIZE1>(mask, alpha, C, A, f0, B, I, beta);
        } else if (tstp >= 0) {
            std::vector<bool> mask_ret(tstp + 1, false);
            std::vector<bool> mask_tv(ntau + 1, false);
            std::vector<bool> mask_les(tstp + 1, false);
            for (i = 0; i <= tstp; i++)
                if (i % nomp == tid)
                    mask_ret[i] = true;
            for (i = 0; i <= ntau; i++)
                if (i % nomp == tid)
                    mask_tv[i] = true;
            for (i = 0; i <= tstp; i++)
                if (mask_ret[i])
                    mask_les[tstp - i] = true;
            incr_convolution_ret<T, GG, SIZE1>(tstp, mask_ret, alpha, C, A, Acc, ft, B, Bcc,
                                               I, h);
            incr_convolution_tv<T, GG, SIZE1>(tstp, mask_tv, alpha, C, A, Acc, f0, ft, B,
                                              Bcc, I, beta, h);
            incr_convolution_les<T, GG, SIZE1>(tstp, mask_les, alpha, C, A, Acc, f0, ft, B,
                                               Bcc, I, beta, h);
        }
    }
    return;
}
/** \brief <b> Returns convolution \f$C = A\ast f\ast B\f$ at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*f*B of the objects with class 'herm_matrix'
* > at a given time step 't=nh'. Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param n
* > [int] number of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] hermitian conjugate of A
* @param ft
* > [herm_matrix] function \f$f\f$
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] hermitian conjugate of B
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &C,
                              herm_matrix<T> &A, herm_matrix<T> &Acc, function<T> &ft,
                              herm_matrix<T> &B, herm_matrix<T> &Bcc,
                              integration::Integrator<T> &I, T beta, T h) {
    int SolveOrder = I.k();
    int ntmin = (tstp == -1 || tstp > SolveOrder ? tstp : SolveOrder);
    if (tstp < -1)
        return;
    int size1 = A.size1();
    std::complex<T> *fttemp;
    assert(C.size1()==size1);
    assert(ft.size1()==size1);
    assert(B.size1()==size1);
    assert(Bcc.size1()==size1);
    assert(Acc.size1()==size1);
    assert(SolveOrder>=0 && SolveOrder <=5);
    assert(C.ntau()>=SolveOrder);
    assert(C.nt()>=ntmin);
    assert(A.nt()>=ntmin);
    assert(Acc.nt()>=ntmin);
    assert(B.nt()>=ntmin);
    assert(Bcc.nt()>=ntmin);
    assert(ft.nt()>=ntmin);
    assert(C.ntau()==A.ntau());
    assert(C.ntau()==Acc.ntau());
    assert(C.ntau()==B.ntau());
    assert(C.ntau()==Bcc.ntau());
    C.set_timestep_zero(tstp);
    fttemp = (tstp == -1 ? ft.ptr(-1) : ft.ptr(0));

    switch (size1) {
    case 1:
        incr_convolution_omp<T, herm_matrix<T>, 1>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 2:
        incr_convolution_omp<T, herm_matrix<T>, 2>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 3:
        incr_convolution_omp<T, herm_matrix<T>, 3>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 4:
        incr_convolution_omp<T, herm_matrix<T>, 4>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 5:
        incr_convolution_omp<T, herm_matrix<T>, 5>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 6:
        incr_convolution_omp<T, herm_matrix<T>, 6>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 8:
        incr_convolution_omp<T, herm_matrix<T>, 8>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    default:
        incr_convolution_omp<T, herm_matrix<T>, LARGESIZE>(
            omp_num_threads, tstp, CPLX(1, 0), C, A, Acc, ft.ptr(-1), fttemp, B, Bcc,
            integration::I<T>(SolveOrder), beta, h);
        break;
    }
}


/** \brief <b> Returns convolution \f$C = A\ast f\ast B\f$ at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*f*B of the objects with class 'herm_matrix'
* > at a given time step 't=nh'. Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param n
* > [int] number of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] hermitian conjugate of A
* @param ft
* > [herm_matrix] function \f$f\f$
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] hermitian conjugate of B
* @param beta
* > inversed temperature
* @param h
* > time step interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T>
void convolution_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &C,
                              herm_matrix<T> &A, herm_matrix<T> &Acc, function<T> &ft,
                              herm_matrix<T> &B, herm_matrix<T> &Bcc,
                              T beta, T h, int SolveOrder) {
    int ntmin = (tstp == -1 || tstp > SolveOrder ? tstp : SolveOrder);
    if (tstp < -1)
        return;
    int size1 = A.size1();
    std::complex<T> *fttemp;
    assert(C.size1()==size1);
    assert(ft.size1()==size1);
    assert(B.size1()==size1);
    assert(Bcc.size1()==size1);
    assert(Acc.size1()==size1);
    assert(SolveOrder>=0 && SolveOrder <=5);
    assert(C.ntau()>=SolveOrder);
    assert(C.nt()>=ntmin);
    assert(A.nt()>=ntmin);
    assert(Acc.nt()>=ntmin);
    assert(B.nt()>=ntmin);
    assert(Bcc.nt()>=ntmin);
    assert(ft.nt()>=ntmin);
    assert(C.ntau()==A.ntau());
    assert(C.ntau()==Acc.ntau());
    assert(C.ntau()==B.ntau());
    assert(C.ntau()==Bcc.ntau());
    C.set_timestep_zero(tstp);
    fttemp = (tstp == -1 ? ft.ptr(-1) : ft.ptr(0));

    switch (size1) {
    case 1:
        incr_convolution_omp<T, herm_matrix<T>, 1>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 2:
        incr_convolution_omp<T, herm_matrix<T>, 2>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 3:
        incr_convolution_omp<T, herm_matrix<T>, 3>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 4:
        incr_convolution_omp<T, herm_matrix<T>, 4>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 5:
        incr_convolution_omp<T, herm_matrix<T>, 5>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 6:
        incr_convolution_omp<T, herm_matrix<T>, 6>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 8:
        incr_convolution_omp<T, herm_matrix<T>, 8>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, ft.ptr(-1), fttemp, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    default:
        incr_convolution_omp<T, herm_matrix<T>, LARGESIZE>(
            omp_num_threads, tstp, CPLX(1, 0), C, A, Acc, ft.ptr(-1), fttemp, B, Bcc,
            integration::I<T>(SolveOrder), beta, h);
        break;
    }
}


/** \brief <b> Returns convolution \f$C = A\ast f\ast B\f$ at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*f*B of the objects with class 'herm_matrix'
* > at a given time step 't=nh'. Here we assume that A and B are hermitian. 
* > Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param tstp
* > [int] index of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param ft
* > [function] function \f$f\f$
* @param B
* > [herm_matrix] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &C,
                              herm_matrix<T> &A, function<T> &ft, herm_matrix<T> &B,
                              integration::Integrator<T> &I, T beta, T h) {
    convolution_timestep_omp<T>(omp_num_threads, tstp, C, A, A, ft, B, B, I, beta, h);
}


/** \brief <b> Returns convolution \f$C = A\ast f\ast B\f$ at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*f*B of the objects with class 'herm_matrix'
* > at a given time step 't=nh'. Here we assume that A and B are hermitian. 
* > Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param tstp
* > [int] index of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param ft
* > [function] function \f$f\f$
* @param B
* > [herm_matrix] contour Green's function
* @param beta
* > inversed temperature
* @param h
* > time step interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T>
void convolution_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &C,
                              herm_matrix<T> &A, function<T> &ft, herm_matrix<T> &B,
                              T beta, T h, int SolveOrder) {
    convolution_timestep_omp<T>(omp_num_threads, tstp, C, A, A, ft, B, B, beta, h, SolveOrder);
}


/** \brief <b> Returns convolution \f$C = A\ast B\f$ for the Matsubara component.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution \f$ C=A\ast B\f$ of the objects with class 'herm_matrix'
* > for the Matsubara component \f$C^\mathrm{M}(\tau)\f$. Here we assume that A and B are hermitian. 
* > Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param B
* > [herm_matrix] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_matsubara_omp(int omp_num_threads, herm_matrix<T> &C, herm_matrix<T> &A,
                               herm_matrix<T> &B, integration::Integrator<T> &I, T beta) {
    convolution_timestep_omp<T>(omp_num_threads, -1, C, A, A, B, B, I, beta, 0.0);
}
/** \brief <b> Returns convolution \f$C = A\ast f\ast B\f$ for the Matsubara component.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution \f$ C=A\ast f \ast B\f$ of the objects with class 'herm_matrix'
* > for the Matsubara component \f$C^\mathrm{M}(\tau)\f$. Here we assume that A and B are hermitian. 
* > Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param ft
* > [function] function \f$f\f$
* @param B
* > [herm_matrix] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_matsubara_omp(int omp_num_threads, herm_matrix<T> &C, herm_matrix<T> &A,
                               function<T> &ft, herm_matrix<T> &B,
                               integration::Integrator<T> &I, T beta) {
    convolution_timestep_omp<T>(omp_num_threads, -1, C, A, A, ft, B, B, I, beta, 0.0);
}
/** \brief <b> Returns convolution \f$C = A\ast f \ast B\f$ for all timesteps.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution \f$C = A\ast f \ast B\f$  of the objects with class 'herm_matrix'
* > for all time steps.
* > Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] hermitian conjugate of A
* @param ft
* > [function] function \f$f\f$
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] hermitian conjugate of B
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_omp(int omp_num_threads, herm_matrix<T> &C, herm_matrix<T> &A,
                     herm_matrix<T> &Acc, function<T> &ft, herm_matrix<T> &B,
                     herm_matrix<T> &Bcc, integration::Integrator<T> &I, T beta, T h) {
    int tstp;
    for (tstp = -1; tstp <= C.nt(); tstp++)
        convolution_timestep_omp<T>(omp_num_threads, tstp, C, A, Acc, ft, B, Bcc, I, beta,
                                    h);
}



/** \brief <b> Returns convolution \f$C = A\ast f \ast B\f$ for all timesteps.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution \f$C = A\ast f \ast B\f$  of the objects with class 'herm_matrix'
* > for all time steps.
* > Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] hermitian conjugate of A
* @param ft
* > [function] function \f$f\f$
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] hermitian conjugate of B
* @param beta
* > inversed temperature
* @param h
* > time step interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T>
void convolution_omp(int omp_num_threads, herm_matrix<T> &C, herm_matrix<T> &A,
                     herm_matrix<T> &Acc, function<T> &ft, herm_matrix<T> &B,
                     herm_matrix<T> &Bcc, T beta, T h, int SolveOrder) {
    int tstp;
    for (tstp = -1; tstp <= C.nt(); tstp++)
        convolution_timestep_omp<T>(omp_num_threads, tstp, C, A, Acc, ft, B, Bcc, beta,
                                    h, SolveOrder);
}



/** \brief <b> Returns convolution \f$C = A\ast B\f$ at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*B of the objects with class 'herm_matrix'
* > at a given time step 't=nh'.
* > Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param tstp
* > [int] index of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] hermitian conjugate of A
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] hermitian conjugate of B
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &C,
                              herm_matrix<T> &A, herm_matrix<T> &Acc, herm_matrix<T> &B,
                              herm_matrix<T> &Bcc, integration::Integrator<T> &I, T beta,
                              T h) {
    int SolveOrder = I.k();
    int ntmin = (tstp == -1 || tstp > SolveOrder ? tstp : SolveOrder);
    if (tstp < -1)
        return;
    int size1 = A.size1();
    assert(C.size1()==size1);
    assert(B.size1()==size1);
    assert(Bcc.size1()==size1);
    assert(Acc.size1()==size1);
    assert(SolveOrder>=0 && SolveOrder <=5);
    assert(C.ntau()>=SolveOrder);
    assert(C.nt()>=ntmin);
    assert(A.nt()>=ntmin);
    assert(Acc.nt()>=ntmin);
    assert(B.nt()>=ntmin);
    assert(Bcc.nt()>=ntmin);
    assert(C.ntau()==A.ntau());
    assert(C.ntau()==Acc.ntau());
    assert(C.ntau()==B.ntau());
    assert(C.ntau()==Bcc.ntau());
    C.set_timestep_zero(tstp);

    switch (size1) {
    case 1:
        incr_convolution_omp<T, herm_matrix<T>, 1>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 2:
        incr_convolution_omp<T, herm_matrix<T>, 2>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 3:
        incr_convolution_omp<T, herm_matrix<T>, 3>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 4:
        incr_convolution_omp<T, herm_matrix<T>, 4>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 5:
        incr_convolution_omp<T, herm_matrix<T>, 5>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 6:
        incr_convolution_omp<T, herm_matrix<T>, 6>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 8:
        incr_convolution_omp<T, herm_matrix<T>, 8>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    default:
        incr_convolution_omp<T, herm_matrix<T>, LARGESIZE>(omp_num_threads, tstp, CPLX(1, 0),
                                                           C, A, Acc, NULL, NULL, B, Bcc,
                                                           integration::I<T>(SolveOrder), beta, h);
        break;
    }
}


/** \brief <b> Returns convolution \f$C = A\ast B\f$ at a given time step</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*B of the objects with class 'herm_matrix'
* > at a given time step 't=nh'.
* > Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param tstp
* > [int] index of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] hermitian conjugate of A
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] hermitian conjugate of B
* @param beta
* > inversed temperature
* @param h
* > time step interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T>
void convolution_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &C,
                              herm_matrix<T> &A, herm_matrix<T> &Acc, herm_matrix<T> &B,
                              herm_matrix<T> &Bcc, T beta,
                              T h, int SolveOrder) {
    int ntmin = (tstp == -1 || tstp > SolveOrder ? tstp : SolveOrder);
    if (tstp < -1)
        return;
    int size1 = A.size1();
    assert(C.size1()==size1);
    assert(B.size1()==size1);
    assert(Bcc.size1()==size1);
    assert(Acc.size1()==size1);
    assert(SolveOrder>=0 && SolveOrder <=5);
    assert(C.ntau()>=SolveOrder);
    assert(C.nt()>=ntmin);
    assert(A.nt()>=ntmin);
    assert(Acc.nt()>=ntmin);
    assert(B.nt()>=ntmin);
    assert(Bcc.nt()>=ntmin);
    assert(C.ntau()==A.ntau());
    assert(C.ntau()==Acc.ntau());
    assert(C.ntau()==B.ntau());
    assert(C.ntau()==Bcc.ntau());
    C.set_timestep_zero(tstp);

    switch (size1) {
    case 1:
        incr_convolution_omp<T, herm_matrix<T>, 1>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 2:
        incr_convolution_omp<T, herm_matrix<T>, 2>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 3:
        incr_convolution_omp<T, herm_matrix<T>, 3>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 4:
        incr_convolution_omp<T, herm_matrix<T>, 4>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 5:
        incr_convolution_omp<T, herm_matrix<T>, 5>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 6:
        incr_convolution_omp<T, herm_matrix<T>, 6>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    case 8:
        incr_convolution_omp<T, herm_matrix<T>, 8>(omp_num_threads, tstp, CPLX(1, 0), C, A,
                                                   Acc, NULL, NULL, B, Bcc,
                                                   integration::I<T>(SolveOrder), beta, h);
        break;
    default:
        incr_convolution_omp<T, herm_matrix<T>, LARGESIZE>(omp_num_threads, tstp, CPLX(1, 0),
                                                           C, A, Acc, NULL, NULL, B, Bcc,
                                                           integration::I<T>(SolveOrder), beta, h);
        break;
    }
}




/** \brief <b> Returns convolution \f$C = A\ast f\ast B\f$ at a given time step.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*f*B of the objects with class 'herm_matrix'
* > at a given time step 't=nh'. Here we assume that A and B are hermitian. 
* > Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param tstp
* > [int] index of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param ft
* > [function] function \f$f\f$
* @param B
* > [herm_matrix] contour Green's function
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &C,
                              herm_matrix<T> &A, herm_matrix<T> &B,
                              integration::Integrator<T> &I, T beta, T h) {
    convolution_timestep_omp<T>(omp_num_threads, tstp, C, A, A, B, B, I, beta, h);
}

/** \brief <b> Returns convolution \f$C = A\ast f\ast B\f$ at a given time step.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*f*B of the objects with class 'herm_matrix'
* > at a given time step 't=nh'. Here we assume that A and B are hermitian. 
* > Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param tstp
* > [int] index of the time step ('t=nh')
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param ft
* > [function] function \f$f\f$
* @param B
* > [herm_matrix] contour Green's function
* @param beta
* > inversed temperature
* @param h
* > time step interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T>
void convolution_timestep_omp(int omp_num_threads, int tstp, herm_matrix<T> &C,
                              herm_matrix<T> &A, herm_matrix<T> &B,
                              T beta, T h, int SolveOrder) {
    convolution_timestep_omp<T>(omp_num_threads, tstp, C, A, A, B, B, beta, h, SolveOrder);
}

/** \brief <b> Returns convolution \f$C = A\ast B\f$ for all timesteps.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*B of the objects with class 'herm_matrix'
* > for all time steps.
* > Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] hermitian conjugate of A
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] hermitian conjugate of B
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time step interval
*/
template <typename T>
void convolution_omp(int omp_num_threads, herm_matrix<T> &C, herm_matrix<T> &A,
                     herm_matrix<T> &Acc, herm_matrix<T> &B, herm_matrix<T> &Bcc,
                     integration::Integrator<T> &I, T beta, T h) {
    int tstp;
    for (tstp = -1; tstp <= C.nt(); tstp++)
        convolution_timestep_omp<T>(omp_num_threads, tstp, C, A, Acc, B, Bcc, I, beta, h);
}

/** \brief <b> Returns convolution \f$C = A\ast B\f$ for all timesteps.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Computes contour convolution C=A*B of the objects with class 'herm_matrix'
* > for all time steps.
* > Works for a scalar and square matrices. `openMP` parallelized version.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param omp_num_threads
* > [int] number of `openMP` threads
* @param C
* > [herm_matrix] Matrix to which the result of the convolution on Matsubara axis is given
* @param A
* > [herm_matrix] contour Green's function
* @param Acc
* > [herm_matrix] hermitian conjugate of A
* @param B
* > [herm_matrix] contour Green's function
* @param Bcc
* > [herm_matrix] hermitian conjugate of B
* @param beta
* > inversed temperature
* @param h
* > time step interval
* @param SolveOrder
* > [int] integrator order
*/
template <typename T>
void convolution_omp(int omp_num_threads, herm_matrix<T> &C, herm_matrix<T> &A,
                     herm_matrix<T> &Acc, herm_matrix<T> &B, herm_matrix<T> &Bcc,
                     T beta, T h, int SolveOrder) {
    int tstp;
    for (tstp = -1; tstp <= C.nt(); tstp++)
        convolution_timestep_omp<T>(omp_num_threads, tstp, C, A, Acc, B, Bcc, beta, h, SolveOrder);
}

#undef CPLX

#endif

} // namespace cntr

#endif  // CNTR_CONVOLUTION_IMPL_H
