#ifndef CNTR_UTILITIES_IMPL_H
#define CNTR_UTILITIES_IMPL_H

#include "cntr_utilities_decl.hpp"
//#include "cntr_exception.hpp"
#include "cntr_elements.hpp"
#include "cntr_herm_matrix_decl.hpp"
#include "cntr_herm_matrix_timestep_view_decl.hpp"
#include "cntr_herm_matrix_timestep_view_impl.hpp"
#include "cntr_herm_matrix_timestep_decl.hpp"
#include "cntr_herm_pseudo_decl.hpp"
#include "cntr_function_decl.hpp"
#include "cntr_convolution_decl.hpp"

namespace cntr {

/// @private
//  K-th order polynomila extrapolate of realtime functions (G^ret, G^vt,
//  G^les)
//  to time t=(n+1), using information at times t=(n-j)  [j=0...k]
template <typename T, class GG, int SIZE1>
void extrapolate_timestep_dispatch(int n, GG &G,
                                   integration::Integrator<T> &I) {
    typedef std::complex<T> cplx;
    T *p1;
    cplx z1, *gtemp, *gtemp1;
    int m, l, j, ntau, nt, n1, k, sg, size1 = G.size1();
    ntau = G.ntau();
    nt = G.nt();
    n1 = n + 1;
    k = I.k();
    sg = G.element_size();
    gtemp = new cplx[sg];
    gtemp1 = new cplx[sg];
    if (n1 > nt || n < k) {
        std::cerr << " k= " << k << " n= " << n << " nt= " << nt << std::endl;
        std::cerr << "extrapolate_tstep: n out of range" << std::endl;
        abort();
    }
    if (k == 0) {
        for (m = 0; m <= ntau; m++)
            element_set<T, SIZE1>(size1, G.tvptr(n1, m), G.tvptr(n, m));
        for (j = 0; j <= n; j++)
            element_set<T, SIZE1>(size1, G.retptr(n1, n1 - j),
                                  G.retptr(n, n - j));
        element_set<T, SIZE1>(size1, G.retptr(n1, 0), G.retptr(n, 0));
        for (j = 1; j <= n1; j++)
            element_set<T, SIZE1>(size1, G.lesptr(j, n1), G.lesptr(j - 1, n));
        element_set<T, SIZE1>(size1, G.lesptr(0, n1), G.lesptr(0, n));
    } else {
        p1 = new T[k + 1];
        for (m = 0; m <= k; m++) {
            p1[m] = 0.0;
            for (l = 0; l <= k; l++)
                p1[m] += (1 - 2 * (l % 2)) * I.poly_interpolation(l, m);
        }
        // vt-component
        for (m = 0; m <= ntau; m++) {
            element_set_zero<T, SIZE1>(size1, gtemp);
            for (l = 0; l <= k; l++)
                element_incr<T, SIZE1>(size1, gtemp, p1[l],
                                       G.tvptr(n - l, m));
            element_set<T, SIZE1>(size1, G.tvptr(n1, m), gtemp);
        }
        // ret-component
        for (j = 0; j <= n - k; j++) {
            element_set_zero<T, SIZE1>(size1, gtemp);
            for (l = 0; l <= k; l++)
                element_incr<T, SIZE1>(size1, gtemp, p1[l],
                                       G.retptr(n - l, n - l - j));
            element_set<T, SIZE1>(size1, G.retptr(n1, n1 - j), gtemp);
        }
        for (j = 0; j <= k; j++) {
            element_set_zero<T, SIZE1>(size1, gtemp);
            for (l = 0; l <= k; l++) {
                if (n - l <= j) {
                    element_set<T, SIZE1>(size1, gtemp1, G.retptr(j, n - l));
                    element_conj<T, SIZE1>(size1, gtemp1);
                    element_smul<T, SIZE1>(size1, gtemp1, -1);
                } else {
                    element_set<T, SIZE1>(size1, gtemp1, G.retptr(n - l, j));
                }
                element_incr<T, SIZE1>(size1, gtemp, p1[l], gtemp1);
            }
            element_set<T, SIZE1>(size1, G.retptr(n1, j), gtemp);
        }
        // les-component
        for (j = 0; j <= k; j++) {
            element_set_zero<T, SIZE1>(size1, gtemp);
            for (l = 0; l <= k; l++) {
                if (n - l < j) {
                    element_set<T, SIZE1>(size1, gtemp1, G.lesptr(n - l, j));
                    element_conj<T, SIZE1>(size1, gtemp1);
                    element_smul<T, SIZE1>(size1, gtemp1, -1);
                } else {
                    element_set<T, SIZE1>(size1, gtemp1, G.lesptr(j, n - l));
                }
                element_incr<T, SIZE1>(size1, gtemp, p1[l], gtemp1);
            }
            element_set<T, SIZE1>(size1, G.lesptr(j, n1), gtemp);
        }
        for (j = k + 1; j <= n1; j++) {
            element_set_zero<T, SIZE1>(size1, gtemp);
            for (l = 0; l <= k; l++)
                element_incr<T, SIZE1>(size1, gtemp, p1[l],
                                       G.lesptr(j - l - 1, n - l));
            element_set<T, SIZE1>(size1, G.lesptr(j, n1), gtemp);
        }
        delete[] p1;
    }
    delete[] gtemp;
    delete[] gtemp1;
}
/// @private
/** \brief <b>  k-th order polynomial extrapolation to t=n+1 of the retarded, lesser and left-mixing components of herm_matrix. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 *
 * > k-th order polynomial extrapolation to t=n+1 of the retarded, lesser and left-mixing components of herm_matrix.
 * > Information at times t=n,n-1,,,,n-k is used.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * > t=n+1 data is obtained by the extrapolation.
 * @param G
 * > herm_matrix to be extrapolated.
 * @param I
 * > Class of 'Integrator', which includes the information of a k-th order polynomial interpolation/extrapolation.
 */
template <typename T>
void extrapolate_timestep(int n, herm_matrix<T> &G,
                          integration::Integrator<T> &I) {
    if (G.size1() == 1)
        extrapolate_timestep_dispatch<T, herm_matrix<T>, 1>(n, G, I);
    else
        extrapolate_timestep_dispatch<T, herm_matrix<T>, LARGESIZE>(n, G, I);
}


/** \brief <b>  k-th order polynomial extrapolation to t=n+1 of the retarded, lesser and left-mixing components of herm_matrix. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 *
 * > k-th order polynomial extrapolation to t=n+1 of the retarded, lesser and left-mixing components of herm_matrix.
 * > Information at times t=n,n-1,,,,n-k is used.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * > t=n+1 data is obtained by the extrapolation.
 * @param G
 * > herm_matrix to be extrapolated.
 * @param SolveOrder
 * > Extrapolation order
 */
template <typename T>
void extrapolate_timestep(int n, herm_matrix<T> &G, int SolveOrder) {
    if (G.size1() == 1)
        extrapolate_timestep_dispatch<T, herm_matrix<T>, 1>(n, G, integration::I<T>(SolveOrder));
    else
        extrapolate_timestep_dispatch<T, herm_matrix<T>, LARGESIZE>(n, G, integration::I<T>(SolveOrder));
}

/// @private
/// k-th order polynomial extrapolation to t=n+1 of the retarded, lesser and left-mixing components of herm_pseudo. </b>
template <typename T>
void extrapolate_timestep(int n, herm_pseudo<T> &G,
                          integration::Integrator<T> &I) {
    if (G.size1() == 1)
        extrapolate_timestep_dispatch<T, herm_pseudo<T>, 1>(n, G, I);
    else
        extrapolate_timestep_dispatch<T, herm_pseudo<T>, LARGESIZE>(n, G, I);
}


/// @private
//  k-th order polynomial extrapolate of contour function
//  to time t=(n+1), using information at times t=(n-j)  [j=0...k]
template <typename T, int SIZE1>
void extrapolate_timestep_dispatch(int n, cntr::function<T> &f,
                                   integration::Integrator<T> &I) {

    typedef std::complex<T> cplx;
    int size1 = f.size1_;
    int size2 = f.size2_;
    int nt = f.nt_;
    int kt = I.k();
    assert(n+1<=nt);
    assert(n>=kt);
    cdmatrix value(size1,size2),output(size1,size2);
    value.setZero();
    output.setZero();
    if (kt == 0) {
        f.get_value(n,value);
        output=value;
        f.set_value(n+1,output);
    } else {
        std::vector<double> weight(kt+1);
        // Set up weights for extrapolation from the interpolated ones
        for(int m=0; m<=kt; m++){
            weight[m]=0.0;
            for(int l=0; l <=kt;l++){
                weight[m] += (1 - 2 * (l % 2))*I.poly_interpolation(l, m);
            }
        }
        for(int l=0;l<= kt;l++){
            f.get_value(n-l,value);
            // std::cout << "extra2d " << value << " " <<  weight[l] << std::endl;
            output+=weight[l]*value;
            // std::cout << "extra2e " << output  << std::endl;
        }
        // std::cout << "extra2f " << output  << std::endl;
        f.set_value(n+1,output);
        // std::cout << "extra2g " << output  << std::endl;
    }
}
/// @private
/** \brief <b>  k-th order polynomial extrapolation to t=n+1 of the contour function. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > k-th order polynomial extrapolation to t=n+1 of the contour function.
 * > Information at times t=n,n-1,,,,n-k is used.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * > t=n+1 data is obtained by the extrapolation.
 * @param f
 * > the contour function to be extrapolated.
 * @param I
 * > Class of 'Integrator', which includes the information of a k-th order polynomial interpolatoion/extrapolation.
 */
template <typename T>
void extrapolate_timestep(int n, function<T> &f,
                          integration::Integrator<T> &I) {
    if (f.size1() == 1)
        extrapolate_timestep_dispatch<T, 1>(n, f, I);
    else
        extrapolate_timestep_dispatch<T, LARGESIZE>(n, f, I);
}


/** \brief <b>  k-th order polynomial extrapolation to t=n+1 of the contour function. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > k-th order polynomial extrapolation to t=n+1 of the contour function.
 * > Information at times t=n,n-1,,,,n-k is used.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param n
 * > t=n+1 data is obtained by the extrapolation.
 * @param f
 * > the contour function to be extrapolated.
 * @param ExtrapolationOrder
 * > Order of extrapolation
 */
template <typename T>
void extrapolate_timestep(int n, function<T> &f, int ExtrapolationOrder) {
    if (f.size1() == 1)
        extrapolate_timestep_dispatch<T, 1>(n, f, integration::I<T>(ExtrapolationOrder));
    else
        extrapolate_timestep_dispatch<T, LARGESIZE>(n, f, integration::I<T>(ExtrapolationOrder));
}


/// @private
//  k-th order polynomial interpolation of contour function
//  to arbitrary point on time interval [tstp-kt,tstp],
//  using information at times t=(n-j)  [j=0...k]
//  TODO: test this more throughly
//  TODO: can we do something smart for edge points ?
template <typename T, int SIZE1>
cdmatrix interpolation_dispatch(int tstp, double tinter, cntr::function<T> &f,
                                integration::Integrator<T> &I) {
    typedef std::complex<T> cplx;
    int kt=I.k();
    int size=f.size1();
    cdmatrix output(size,size),tmp(size,size);
    output.setZero();
    double t0=float(tstp-kt);
    double t1;

    // The k-th order approximation polynomial through function
    // values {(x_l,f(x_l)),l=0...k} is given by
    // P(x) = sum_{p,l=0}^{k} x^p*I.interpolation(p,l)*f_l

    for(int l=0;l<=kt;l++){
        t1=1.0;
        double weight = I.poly_interpolation(0,l);
        f.get_value(tstp-kt+l,tmp);
        for(int p=1;p<=kt;p++){
            t1*=tinter-t0;
            weight +=t1 * I.poly_interpolation(p,l);
        }

        output+=tmp*weight;

    }
    return output;
}
/// @private
/** \brief <b>  k-th order polynomial interpolation of the contour function </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > k-th order polynomial interpolation of the contour function to arbitrary point 'tinter' on time interval [tstp-k,tstp].
 * > One uses information at times t=(n-j)  [j=0...k]
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > [int] Time step.
 * @param tinter
 * > [double] interpolation time point
 * @param f
 * > [function] the contour function to be interpolated
 * @param I
 * > [Integrator] Class of 'Integrator', which includes the information of a k-th order polynomial interpolation/extrapolation.
 */
template <typename T>
cdmatrix interpolation(int tstp,double tinter,function<T> &f,
                          integration::Integrator<T> &I) {

    int kt=I.k();
    assert(tstp>=kt);
    assert(tinter<=tstp);
    assert(tinter>=(tstp-kt));

    if (f.size1() == 1)
        return interpolation_dispatch<T, 1>(tstp,tinter, f, I);
    else
        return interpolation_dispatch<T, LARGESIZE>(tstp,tinter, f, I);
}

/** \brief <b>  k-th order polynomial interpolation of the contour function </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > k-th order polynomial interpolation of the contour function to arbitrary point 'tinter' on time interval [tstp-k,tstp].
 * > One uses information at times t=(n-j)  [j=0...k]
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > [int] Time step.
 * @param tinter
 * > [double] interpolation time point
 * @param f
 * > [function] the contour function to be interpolated
 * @param InterpolationOrder
 * > [int] Order of interpolation
 */
template <typename T>
cdmatrix interpolation(int tstp,double tinter,function<T> &f,
                          int InterpolationOrder) {

    assert(tstp>=InterpolationOrder);
    assert(tinter<=tstp);
    assert(tinter>=(tstp-InterpolationOrder));

    if (f.size1() == 1)
        return interpolation_dispatch<T, 1>(tstp,tinter, f, integration::I<T>(InterpolationOrder));
    else
        return interpolation_dispatch<T, LARGESIZE>(tstp,tinter, f, integration::I<T>(InterpolationOrder));
}

/** \brief <b>  Set t=0 components of the two-time contour object from the Matsubara component. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Set t=0 components of herm_matrix (retarded,lesser,left-mixing components) from the Matsubara component.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param G
 * > herm_matrix to be modified.
 */
template <typename T>
void set_t0_from_mat(herm_matrix<T> &G) {
    int m, ntau = G.ntau(), size1 = G.size1();
    for (m = 0; m <= ntau; m++) {
        element_set<T, LARGESIZE>(size1, G.tvptr(0, m), G.matptr(ntau - m));
        element_smul<T, LARGESIZE>(size1, G.tvptr(0, m),
                                   std::complex<T>(0, G.sig()));
    }
    element_set<T, LARGESIZE>(size1, G.lesptr(0, 0), G.tvptr(0, 0));
    element_set<T, LARGESIZE>(size1, G.retptr(0, 0), G.matptr(0));
    element_smul<T, LARGESIZE>(size1, G.retptr(0, 0), std::complex<T>(0, 1));
    element_incr<T, LARGESIZE>(size1, G.retptr(0, 0), -1, G.lesptr(0, 0));
}


// Fix first kt points from mat
/** \brief <b>  Set components of herm_matrix at t=0,1,..kt from the Matsubara component. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Set components of herm_matrix at t=0,1,..kt (retarded,lesser,left-mixing components) from the Matsubara component.
 *
 * \note This function put contant values for t=0,1,..kt for each components.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param G
 * > herm_matrix to be modified
 * @param kt
 * > the time up to which the values are set.
 */
template <typename T>
void set_tk_from_mat(herm_matrix<T> &G,int kt){
   int m,n,ntau=G.ntau(),size1=G.size1();
   for(n=0;n<=kt;n++){
    for(m=0;m<=ntau;m++){
        element_set<T,LARGESIZE>(size1,G.tvptr(n,m),G.matptr(ntau-m));
        element_smul<T,LARGESIZE>(size1,G.tvptr(n,m),std::complex<T>(0,G.sig()));
    }
    element_set<T,LARGESIZE>(size1,G.lesptr(n,n),G.tvptr(0,0));
    element_set<T,LARGESIZE>(size1,G.retptr(n,n),G.matptr(0));
    element_smul<T,LARGESIZE>(size1,G.retptr(n,n),std::complex<T>(0,1));
    element_incr<T,LARGESIZE>(size1,G.retptr(n,n),-1,G.lesptr(0,0));
   }
}

/// @private
template <typename T>
void set_t0_from_mat(herm_pseudo<T> &G) {
    int m, ntau = G.ntau(), size1 = G.size1();
    for (m = 0; m <= ntau; m++) {
        element_set<T, LARGESIZE>(size1, G.tvptr(0, m), G.matptr(ntau - m));
        element_smul<T, LARGESIZE>(size1, G.tvptr(0, m),
                                   std::complex<T>(0, G.sig()));
    }
    element_set<T, LARGESIZE>(size1, G.lesptr(0, 0), G.tvptr(0, 0));
    element_set<T, LARGESIZE>(size1, G.retptr(0, 0), G.matptr(0));
    element_smul<T, LARGESIZE>(size1, G.retptr(0, 0), std::complex<T>(0, 1));
    // element_incr<T,LARGESIZE>(size1,G.retptr(0,0),-1,G.lesptr(0,0));
}

/** \brief <b> Returns the result of the correlation energy at a given time-step from the time-diagonal convolution</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calculates the correlation energy given by \f$ -i/2\f$ Tr\f$[G*\Sigma]^<\f$ at a given time step 'tstp'.
* > The objects 'G' and 'Sigma' are of the class type 'herm_matrix'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > [int] time step
* @param G
* > [herm_matrix] contour Green's function
* @param Sigma
* > [herm_matrix] self-energy
* @param I
* > [Integrator] integrator class
* @param beta
* > inversed temperature
* @param h
* > time interval
*/
template <typename T>
T correlation_energy(int tstp, herm_matrix<T> &G, herm_matrix<T> &Sigma,
		     integration::Integrator<T> &I, T beta, T h){
  int size1=G.size1();
  std::complex<T> trGxSGM (0.0,0.0);
  std::complex<T> *GxSGM = new std::complex<T>[size1*size1];

  convolution_density_matrix<T,herm_matrix<T>>(tstp, GxSGM, Sigma, G, I, beta, h);
  //trGxSGM = (0.0,0.0);
  for(int i=0; i< size1; i++){
    trGxSGM += GxSGM[i*size1 + i];
  }
  delete[] GxSGM;
  return 0.5*trGxSGM.real();
}

/** \brief <b> Returns the result of the correlation energy at a given time-step from the time-diagonal convolution</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >  Calculates the correlation energy given by \f$ -i/2\f$ Tr\f$[G*\Sigma]^<\f$ at a given time step 'tstp'.
* > The objects 'G' and 'Sigma' are of the class type 'herm_matrix'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > [int] time step
* @param G
* > [herm_matrix] contour Green's function
* @param Sigma
* > [herm_matrix] self-energy
* @param beta
* > inversed temperature
* @param h
* > time interval
* @param SolveOrder
* > [int] integration order
*/
template <typename T>
T correlation_energy(int tstp, herm_matrix<T> &G, herm_matrix<T> &Sigma,
              T beta, T h, int SolveOrder){
  int size1=G.size1();
  std::complex<T> trGxSGM (0.0,0.0);
  std::complex<T> *GxSGM = new std::complex<T>[size1*size1];

  convolution_density_matrix<T,herm_matrix<T>>(tstp, GxSGM, Sigma, G, integration::I<T>(SolveOrder), beta, h);
  //trGxSGM = (0.0,0.0);
  for(int i=0; i< size1; i++){
    trGxSGM += GxSGM[i*size1 + i];
  }
  return 0.5*trGxSGM.real();

  delete[] GxSGM;
}


/// @private
template <typename T, class GG, int SIZE1>
T distance_norm2_ret_dispatch(int tstp, GG &g1, GG &g2) {
    int size1 = g1.size1(), i;
    T err = 0.0;
    std::complex<T> *temp = new std::complex<T>[size1 * size1];
    assert(tstp>=-1 && g1.nt()>=tstp && g2.nt()>=tstp && g1.size1() == g2.size1());

    if (tstp == -1)
        return 0.0;
    for (i = 0; i <= tstp; i++) {
        element_set<T, SIZE1>(size1, temp, g1.retptr(tstp, i));
        element_incr<T, SIZE1>(size1, temp, -1.0, g2.retptr(tstp, i));
        err += element_norm2<T, SIZE1>(size1, temp);
    }
    delete[] temp;
    return err;
}
/// @private
template <typename T, class GG, int SIZE1>
T distance_norm2_tv_dispatch(int tstp, GG &g1, GG &g2) {
    int size1 = g1.size1(), ntau = g1.ntau(), i;
    T err = 0.0;
    std::complex<T> *temp = new std::complex<T>[size1 * size1];
    assert(tstp>=-1 && g1.nt()>=tstp && g2.nt() >= tstp);
    assert(g1.ntau()==g2.ntau() && g1.size1()==g2.size1());
    if (tstp == -1)
        return 0.0;
    for (i = 0; i <= ntau; i++) {
        element_set<T, SIZE1>(size1, temp, g1.tvptr(tstp, i));
        element_incr<T, SIZE1>(size1, temp, -1.0, g2.tvptr(tstp, i));
        err += element_norm2<T, SIZE1>(size1, temp);
    }
    delete[] temp;
    return err;
}
/// @private
template <typename T, class GG, int SIZE1>
T distance_norm2_les_dispatch(int tstp, GG &g1, GG &g2) {
    int size1 = g1.size1(), i;
    T err = 0.0;
    std::complex<T> *temp = new std::complex<T>[size1 * size1];

    assert(tstp>=-1 && g1.nt() >= tstp && g2.nt() >= tstp && g1.size1() == g2.size1());
    if (tstp == -1)
        return 0.0;
    for (i = 0; i <= tstp; i++) {
        element_set<T, SIZE1>(size1, temp, g1.lesptr(i, tstp));
        element_incr<T, SIZE1>(size1, temp, -1.0, g2.lesptr(i, tstp));
        err += element_norm2<T, SIZE1>(size1, temp);
    }
    delete[] temp;
    return err;
}

/// @private
template <typename T, class GG, int SIZE1>
T distance_norm2_dispatch(int tstp, GG &g1, GG &g2) {
    int size1 = g1.size1(), ntau = g1.ntau(), i;
    T err = 0.0;
    std::complex<T> *temp = new std::complex<T>[size1 * size1];

    assert(tstp >= -1 && g1.nt() >= tstp && g2.nt() >= tstp);
    assert(g1.ntau() == g2.ntau());
    assert(g1.size1() == g2.size1());

    if (tstp == -1) {
        for (i = 0; i <= ntau; i++) {
            element_set<T, SIZE1>(size1, temp, g1.matptr(i));
            element_incr<T, SIZE1>(size1, temp, -1.0, g2.matptr(i));
            err += element_norm2<T, SIZE1>(size1, temp);
        }
    } else {
        for (i = 0; i <= tstp; i++) {
            element_set<T, SIZE1>(size1, temp, g1.retptr(tstp, i));
            element_incr<T, SIZE1>(size1, temp, -1.0, g2.retptr(tstp, i));
            err += element_norm2<T, SIZE1>(size1, temp);
        }
        for (i = 0; i <= ntau; i++) {
            element_set<T, SIZE1>(size1, temp, g1.tvptr(tstp, i));
            element_incr<T, SIZE1>(size1, temp, -1.0, g2.tvptr(tstp, i));
            err += element_norm2<T, SIZE1>(size1, temp);
        }
        for (i = 0; i <= tstp; i++) {
            element_set<T, SIZE1>(size1, temp, g1.lesptr(i, tstp));
            element_incr<T, SIZE1>(size1, temp, -1.0, g2.lesptr(i, tstp));
            err += element_norm2<T, SIZE1>(size1, temp);
        }
    }
    delete[] temp;
    return err;
}








/// @private
template <typename T, int SIZE1>
T distance_norm2_ret_dispatch(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep_view<T> &g2) {
    int size1 = g1.size1();
    int s1 = size1*size1;
    T err = 0.0;
    std::complex<T> *temp = new std::complex<T>[size1 * size1];

    if (tstp == -1)
        return 0.0;
    for (int i = 0; i <= tstp; i++) {
        element_set<T, SIZE1>(size1, temp, g1.ret_ + i * s1);
        element_incr<T, SIZE1>(size1, temp, -1.0, g2.ret_ + i * s1);
        err += element_norm2<T, SIZE1>(size1, temp);
    }
    delete[] temp;
    return err;
}
/// @private
template <typename T, int SIZE1>
T distance_norm2_tv_dispatch(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep_view<T> &g2) {
    int size1 = g1.size1(), ntau = g1.ntau();
    int s1 = size1*size1;
    T err = 0.0;
    std::complex<T> *temp = new std::complex<T>[size1 * size1];

    if (tstp == -1)
        return 0.0;
    for (int i = 0; i <= ntau; i++) {
        element_set<T, SIZE1>(size1, temp, g1.tv_ + i * s1);
        element_incr<T, SIZE1>(size1, temp, -1.0, g2.tv_ + i * s1);
        err += element_norm2<T, SIZE1>(size1, temp);
    }
    delete[] temp;
    return err;
}
/// @private
template <typename T, int SIZE1>
T distance_norm2_les_dispatch(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep_view<T> &g2) {
    int size1 = g1.size1();
    int s1 = size1*size1;
    T err = 0.0;
    std::complex<T> *temp = new std::complex<T>[size1 * size1];

    if (tstp == -1)
        return 0.0;
    for (int i = 0; i <= tstp; i++) {
        element_set<T, SIZE1>(size1, temp, g1.les_ + i * s1);
        element_incr<T, SIZE1>(size1, temp, -1.0, g2.les_ + i * s1);
        err += element_norm2<T, SIZE1>(size1, temp);
    }
    delete[] temp;
    return err;
}

/// @private
template <typename T, int SIZE1>
T distance_norm2_dispatch(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep_view<T> &g2) {
    int size1 = g1.size1(), ntau = g1.ntau();
    int s1 = size1*size1;
    T err = 0.0;
    std::complex<T> *temp = new std::complex<T>[size1 * size1];

    if (tstp == -1) {
        for (int i = 0; i <= ntau; i++) {
            element_set<T, SIZE1>(size1, temp, g1.mat_ + i * s1);
            element_incr<T, SIZE1>(size1, temp, -1.0, g2.mat_ + i * s1);
            err += element_norm2<T, SIZE1>(size1, temp);
        }
    } else {
        for (int i = 0; i <= tstp; i++) {
            element_set<T, SIZE1>(size1, temp, g1.ret_ + i * s1);
            element_incr<T, SIZE1>(size1, temp, -1.0, g2.ret_ + i * s1);
            err += element_norm2<T, SIZE1>(size1, temp);
        }
        for (int i = 0; i <= ntau; i++) {
            element_set<T, SIZE1>(size1, temp, g1.tv_ + i * s1);
            element_incr<T, SIZE1>(size1, temp, -1.0, g2.tv_ + i * s1);
            err += element_norm2<T, SIZE1>(size1, temp);
        }
        for (int i = 0; i <= tstp; i++) {
            element_set<T, SIZE1>(size1, temp, g1.les_ + i * s1);
            element_incr<T, SIZE1>(size1, temp, -1.0, g2.les_ + i * s1);
            err += element_norm2<T, SIZE1>(size1, temp);
        }
    }
    delete[] temp;
    return err;
}


/** \brief <b>  Evaluate the Euclidean norm between the retarded components of two `herm_matrix` at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the retarded components of two `herm_matrix` (\f$g_1,g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix
 * @param g2
 * > herm_matrix
 */
template <typename T>
T distance_norm2_ret(int tstp, herm_matrix<T> &g1, herm_matrix<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.nt() >= tstp);
    assert(g2.nt() >= tstp);
    if (g1.size1() == 1)
        return distance_norm2_ret_dispatch<T, herm_matrix<T>, 1>(tstp, g1,
                                                                 g2);
    else
        return distance_norm2_ret_dispatch<T, herm_matrix<T>, LARGESIZE>(
            tstp, g1, g2);
}

/** \brief <b>  Evaluate the Euclidean norm between the retarded components a `herm_matrix_timestep` and `herm_matrix`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the retarded components of a `herm_matrix_timestep` \f$g_1\f$ and a `herm_matrix` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep
 * @param g2
 * > herm_matrix
 */
template <typename T>
T distance_norm2_ret(int tstp, herm_matrix_timestep<T> &g1, herm_matrix<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.nt() >= tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_ret_dispatch<T, 1>(tstp, g1_view, g2_view);
    else
        return distance_norm2_ret_dispatch<T, LARGESIZE>(tstp, g1_view, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between the retarded components a `herm_matrix` and `herm_matrix_timestep`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the retarded components of a `herm_matrix` \f$g_1\f$ and a `herm_matrix_timestep` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix
 * @param g2
 * > herm_matrix_timestep
 */
template <typename T>
T distance_norm2_ret(int tstp, herm_matrix<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.nt() >= tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_ret_dispatch<T, 1>(tstp, g1_view, g2_view);
    else
        return distance_norm2_ret_dispatch<T, LARGESIZE>(tstp, g1_view, g2_view);
}


/** \brief <b>  Evaluate the Euclidean norm between the retarded components a `herm_matrix_timestep` and `herm_matrix_timestep`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the retarded components of a `herm_matrix_timestep` \f$g_1\f$ and a `herm_matrix_timestep` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep
 * @param g2
 * > herm_matrix_timestep
 */
template <typename T>
T distance_norm2_ret(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_ret_dispatch<T, 1>(tstp, g1_view, g2_view);
    else
        return distance_norm2_ret_dispatch<T, LARGESIZE>(tstp, g1_view, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between the retarded components a `herm_matrix_timestep_view` and `herm_matrix`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the retarded components of a `herm_matrix_timestep_view` \f$g_1\f$ and a `herm_matrix` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep_view
 * @param g2
 * > herm_matrix
 */
template <typename T>
T distance_norm2_ret(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.nt() >= tstp);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_ret_dispatch<T, 1>(tstp, g1, g2_view);
    else
        return distance_norm2_ret_dispatch<T, LARGESIZE>(tstp, g1, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between the retarded components a `herm_matrix` and `herm_matrix_timestep_view`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the retarded components of a `herm_matrix` \f$g_1\f$ and a `herm_matrix_timestep_view` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix
 * @param g2
 * > herm_matrix_timestep_view
 */
template <typename T>
T distance_norm2_ret(int tstp, herm_matrix<T> &g1, herm_matrix_timestep_view<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.nt() >= tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    if (g1.size1() == 1)
        return distance_norm2_ret_dispatch<T, 1>(tstp, g1_view, g2);
    else
        return distance_norm2_ret_dispatch<T, LARGESIZE>(tstp, g1_view, g2);
}


/** \brief <b>  Evaluate the Euclidean norm between the retarded components a `herm_matrix_timestep` and `herm_matrix_timestep_view`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the retarded components of a `herm_matrix_timestep` \f$g_1\f$ and a `herm_matrix_timestep_view` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep
 * @param g2
 * > herm_matrix_timestep_view
 */
template <typename T>
T distance_norm2_ret(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep_view<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    if (g1.size1() == 1)
        return distance_norm2_ret_dispatch<T, 1>(tstp, g1_view, g2);
    else
        return distance_norm2_ret_dispatch<T, LARGESIZE>(tstp, g1_view, g2);
}

/** \brief <b>  Evaluate the Euclidean norm between the retarded components a `herm_matrix_timestep_view` and `herm_matrix_timestep`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the retarded components of a `herm_matrix_timestep_view` \f$g_1\f$ and a `herm_matrix_timestep` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep_view
 * @param g2
 * > herm_matrix_timestep
 */
template <typename T>
T distance_norm2_ret(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_ret_dispatch<T, 1>(tstp, g1, g2_view);
    else
        return distance_norm2_ret_dispatch<T, LARGESIZE>(tstp, g1, g2_view);
}


/** \brief <b>  Evaluate the Euclidean norm between the retarded components a `herm_matrix_timestep_view` and `herm_matrix_timestep_view`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the retarded components of a `herm_matrix_timestep_view` \f$g_1\f$ and a `herm_matrix_timestep_view` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep_view
 * @param g2
 * > herm_matrix_timestep_view
 */
template <typename T>
T distance_norm2_ret(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep_view<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    if (g1.size1() == 1)
        return distance_norm2_ret_dispatch<T, 1>(tstp, g1, g2);
    else
        return distance_norm2_ret_dispatch<T, LARGESIZE>(tstp, g1, g2);
}



/** \brief <b>  Evaluate the Euclidean norm between the left-mixing components of two herm_matrces at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the left-mixing components of two herm_matrces (\f$g_1,g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix
 * @param g2
 * > herm_matrix
 */
template <typename T>
T distance_norm2_tv(int tstp, herm_matrix<T> &g1, herm_matrix<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.nt() >= tstp);
    assert(g2.nt() >= tstp);
    if (g1.size1() == 1)
        return distance_norm2_tv_dispatch<T, herm_matrix<T>, 1>(tstp, g1, g2);
    else
        return distance_norm2_tv_dispatch<T, herm_matrix<T>, LARGESIZE>(
            tstp, g1, g2);
}

/** \brief <b>  Evaluate the Euclidean norm between the left-mixing components a `herm_matrix_timestep` and `herm_matrix`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the left-mixing components of a `herm_matrix_timestep` \f$g_1\f$ and a `herm_matrix` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep
 * @param g2
 * > herm_matrix
 */
template <typename T>
T distance_norm2_tv(int tstp, herm_matrix_timestep<T> &g1, herm_matrix<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.nt() >= tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_tv_dispatch<T, 1>(tstp, g1_view, g2_view);
    else
        return distance_norm2_tv_dispatch<T, LARGESIZE>(tstp, g1_view, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between the left-mixing components a `herm_matrix` and `herm_matrix_timestep`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the left-mixing components of a `herm_matrix` \f$g_1\f$ and a `herm_matrix_timestep` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix
 * @param g2
 * > herm_matrix_timestep
 */
template <typename T>
T distance_norm2_tv(int tstp, herm_matrix<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.nt() >= tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_tv_dispatch<T, 1>(tstp, g1_view, g2_view);
    else
        return distance_norm2_tv_dispatch<T, LARGESIZE>(tstp, g1_view, g2_view);
}


/** \brief <b>  Evaluate the Euclidean norm between the left-mixing components a `herm_matrix_timestep` and `herm_matrix_timestep`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the left-mixing components of a `herm_matrix_timestep` \f$g_1\f$ and a `herm_matrix_timestep` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep
 * @param g2
 * > herm_matrix_timestep
 */
template <typename T>
T distance_norm2_tv(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_tv_dispatch<T, 1>(tstp, g1_view, g2_view);
    else
        return distance_norm2_tv_dispatch<T, LARGESIZE>(tstp, g1_view, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between the left-mixing components a `herm_matrix_timestep_view` and `herm_matrix`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the left-mixing components of a `herm_matrix_timestep_view` \f$g_1\f$ and a `herm_matrix` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep_view
 * @param g2
 * > herm_matrix
 */
template <typename T>
T distance_norm2_tv(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.nt() >= tstp);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_tv_dispatch<T, 1>(tstp, g1, g2_view);
    else
        return distance_norm2_tv_dispatch<T, LARGESIZE>(tstp, g1, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between the left-mixing components a `herm_matrix` and `herm_matrix_timestep_view`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the left-mixing components of a `herm_matrix` \f$g_1\f$ and a `herm_matrix_timestep_view` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix
 * @param g2
 * > herm_matrix_timestep_view
 */
template <typename T>
T distance_norm2_tv(int tstp, herm_matrix<T> &g1, herm_matrix_timestep_view<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.nt() >= tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    if (g1.size1() == 1)
        return distance_norm2_tv_dispatch<T, 1>(tstp, g1_view, g2);
    else
        return distance_norm2_tv_dispatch<T, LARGESIZE>(tstp, g1_view, g2);
}


/** \brief <b>  Evaluate the Euclidean norm between the left-mixing components a `herm_matrix_timestep` and `herm_matrix_timestep_view`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the left-mixing components of a `herm_matrix_timestep` \f$g_1\f$ and a `herm_matrix_timestep_view` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep
 * @param g2
 * > herm_matrix_timestep_view
 */
template <typename T>
T distance_norm2_tv(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep_view<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    if (g1.size1() == 1)
        return distance_norm2_tv_dispatch<T, 1>(tstp, g1_view, g2);
    else
        return distance_norm2_tv_dispatch<T, LARGESIZE>(tstp, g1_view, g2);
}

/** \brief <b>  Evaluate the Euclidean norm between the left-mixing components a `herm_matrix_timestep_view` and `herm_matrix_timestep`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the left-mixing components of a `herm_matrix_timestep_view` \f$g_1\f$ and a `herm_matrix_timestep` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep_view
 * @param g2
 * > herm_matrix_timestep
 */
template <typename T>
T distance_norm2_tv(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_tv_dispatch<T, 1>(tstp, g1, g2_view);
    else
        return distance_norm2_tv_dispatch<T, LARGESIZE>(tstp, g1, g2_view);
}


/** \brief <b>  Evaluate the Euclidean norm between the left-mixing components a `herm_matrix_timestep_view` and `herm_matrix_timestep_view`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the left-mixing components of a `herm_matrix_timestep_view` \f$g_1\f$ and a `herm_matrix_timestep_view` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep_view
 * @param g2
 * > herm_matrix_timestep_view
 */
template <typename T>
T distance_norm2_tv(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep_view<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    if (g1.size1() == 1)
        return distance_norm2_tv_dispatch<T, 1>(tstp, g1, g2);
    else
        return distance_norm2_tv_dispatch<T, LARGESIZE>(tstp, g1, g2);
}


/** \brief <b>  Evaluate the Euclidean norm between the lesser components of two herm_matrces at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the lesser components of two herm_matrces (\f$g_1,g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix
 * @param g2
 * > herm_matrix
 */
template <typename T>
T distance_norm2_les(int tstp, herm_matrix<T> &g1, herm_matrix<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.nt() >= tstp);
    assert(g2.nt() >= tstp);
    if (g1.size1() == 1)
        return distance_norm2_les_dispatch<T, herm_matrix<T>, 1>(tstp, g1,
                                                                 g2);
    else
        return distance_norm2_les_dispatch<T, herm_matrix<T>, LARGESIZE>(
            tstp, g1, g2);
}



/** \brief <b>  Evaluate the Euclidean norm between the lesser components a `herm_matrix_timestep` and `herm_matrix`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the lesser components of a `herm_matrix_timestep` \f$g_1\f$ and a `herm_matrix` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep
 * @param g2
 * > herm_matrix
 */
template <typename T>
T distance_norm2_les(int tstp, herm_matrix_timestep<T> &g1, herm_matrix<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.nt() >= tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_les_dispatch<T, 1>(tstp, g1_view, g2_view);
    else
        return distance_norm2_les_dispatch<T, LARGESIZE>(tstp, g1_view, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between the lesser components a `herm_matrix` and `herm_matrix_timestep`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the lesser components of a `herm_matrix` \f$g_1\f$ and a `herm_matrix_timestep` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix
 * @param g2
 * > herm_matrix_timestep
 */
template <typename T>
T distance_norm2_les(int tstp, herm_matrix<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.nt() >= tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_les_dispatch<T, 1>(tstp, g1_view, g2_view);
    else
        return distance_norm2_les_dispatch<T, LARGESIZE>(tstp, g1_view, g2_view);
}


/** \brief <b>  Evaluate the Euclidean norm between the lesser components a `herm_matrix_timestep` and `herm_matrix_timestep`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the lesser components of a `herm_matrix_timestep` \f$g_1\f$ and a `herm_matrix_timestep` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep
 * @param g2
 * > herm_matrix_timestep
 */
template <typename T>
T distance_norm2_les(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_les_dispatch<T, 1>(tstp, g1_view, g2_view);
    else
        return distance_norm2_les_dispatch<T, LARGESIZE>(tstp, g1_view, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between the lesser components a `herm_matrix_timestep_view` and `herm_matrix`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the lesser components of a `herm_matrix_timestep_view` \f$g_1\f$ and a `herm_matrix` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep_view
 * @param g2
 * > herm_matrix
 */
template <typename T>
T distance_norm2_les(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.nt() >= tstp);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_les_dispatch<T, 1>(tstp, g1, g2_view);
    else
        return distance_norm2_les_dispatch<T, LARGESIZE>(tstp, g1, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between the lesser components a `herm_matrix` and `herm_matrix_timestep_view`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the lesser components of a `herm_matrix` \f$g_1\f$ and a `herm_matrix_timestep_view` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix
 * @param g2
 * > herm_matrix_timestep_view
 */
template <typename T>
T distance_norm2_les(int tstp, herm_matrix<T> &g1, herm_matrix_timestep_view<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.nt() >= tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    if (g1.size1() == 1)
        return distance_norm2_les_dispatch<T, 1>(tstp, g1_view, g2);
    else
        return distance_norm2_les_dispatch<T, LARGESIZE>(tstp, g1_view, g2);
}


/** \brief <b>  Evaluate the Euclidean norm between the lesser components a `herm_matrix_timestep` and `herm_matrix_timestep_view`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the lesser components of a `herm_matrix_timestep` \f$g_1\f$ and a `herm_matrix_timestep_view` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep
 * @param g2
 * > herm_matrix_timestep_view
 */
template <typename T>
T distance_norm2_les(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep_view<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    if (g1.size1() == 1)
        return distance_norm2_les_dispatch<T, 1>(tstp, g1_view, g2);
    else
        return distance_norm2_les_dispatch<T, LARGESIZE>(tstp, g1_view, g2);
}

/** \brief <b>  Evaluate the Euclidean norm between the lesser components a `herm_matrix_timestep_view` and `herm_matrix_timestep`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the lesser components of a `herm_matrix_timestep_view` \f$g_1\f$ and a `herm_matrix_timestep` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep_view
 * @param g2
 * > herm_matrix_timestep
 */
template <typename T>
T distance_norm2_les(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_les_dispatch<T, 1>(tstp, g1, g2_view);
    else
        return distance_norm2_les_dispatch<T, LARGESIZE>(tstp, g1, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between the lesser components a `herm_matrix_timestep_view` and `herm_matrix_timestep_view`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between the lesser components of a `herm_matrix_timestep_view` \f$g_1\f$ and a `herm_matrix_timestep_view` (\f$g_2\f$) at a given time step (tstp).
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep_view
 * @param g2
 * > herm_matrix_timestep_view
 */
template <typename T>
T distance_norm2_les(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep_view<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    if (g1.size1() == 1)
        return distance_norm2_les_dispatch<T, 1>(tstp, g1, g2);
    else
        return distance_norm2_les_dispatch<T, LARGESIZE>(tstp, g1, g2);
}



/** \brief <b>  Evaluate the Euclidean norm between two `herm_matrix` at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between two`herm_matrix` (\f$g_1,g_2\f$) at a given time step (tstp).
 * > To evaluate the norm, the elements of retarded, lesser and left-mixing components at the time step is used.
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > [int]time step
 * @param g1
 * > herm_matrix
 * @param g2
 * > herm_matrix
 */
template <typename T>
T distance_norm2(int tstp, herm_matrix<T> &g1, herm_matrix<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.nt() >= tstp);
    assert(g2.nt() >= tstp);
    if (g1.size1() == 1)
        return distance_norm2_dispatch<T, herm_matrix<T>, 1>(tstp, g1, g2);
    else
        return distance_norm2_dispatch<T, herm_matrix<T>, LARGESIZE>(tstp, g1,
                                                                     g2);
}


/** \brief <b>  Evaluate the Euclidean norm between  a `herm_matrix_timestep` and `herm_matrix`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between  of a `herm_matrix_timestep` \f$g_1\f$ and a `herm_matrix` (\f$g_2\f$) at a given time step (tstp).
 * > To evaluate the norm, the elements of retarded, lesser and left-mixing components at the time step is used. The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep
 * @param g2
 * > herm_matrix
 */
template <typename T>
T distance_norm2(int tstp, herm_matrix_timestep<T> &g1, herm_matrix<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.nt() >= tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_dispatch<T, 1>(tstp, g1_view, g2_view);
    else
        return distance_norm2_dispatch<T, LARGESIZE>(tstp, g1_view, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between  a `herm_matrix` and `herm_matrix_timestep`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between  of a `herm_matrix` \f$g_1\f$ and a `herm_matrix_timestep` (\f$g_2\f$) at a given time step (tstp).
 * > To evaluate the norm, the elements of retarded, lesser and left-mixing components at the time step is used. The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix
 * @param g2
 * > herm_matrix_timestep
 */
template <typename T>
T distance_norm2(int tstp, herm_matrix<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.nt() >= tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_dispatch<T, 1>(tstp, g1_view, g2_view);
    else
        return distance_norm2_dispatch<T, LARGESIZE>(tstp, g1_view, g2_view);
}


/** \brief <b>  Evaluate the Euclidean norm between  a `herm_matrix_timestep` and `herm_matrix_timestep`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between  of a `herm_matrix_timestep` \f$g_1\f$ and a `herm_matrix_timestep` (\f$g_2\f$) at a given time step (tstp).
 * > To evaluate the norm, the elements of retarded, lesser and left-mixing components at the time step is used. The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep
 * @param g2
 * > herm_matrix_timestep
 */
template <typename T>
T distance_norm2(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_dispatch<T, 1>(tstp, g1_view, g2_view);
    else
        return distance_norm2_dispatch<T, LARGESIZE>(tstp, g1_view, g2_view);
}


/// @private
template <typename T>
T distance_norm2(herm_matrix_timestep<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == g2.tstp());
    herm_matrix_timestep_view<T> g1_view(g1.tstp(), g1);
    herm_matrix_timestep_view<T> g2_view(g1.tstp(), g2);
    if (g1.size1() == 1)
        return distance_norm2_dispatch<T, 1>(g1.tstp(), g1_view, g2_view);
    else
        return distance_norm2_dispatch<T, LARGESIZE>(g1.tstp(), g1_view, g2_view);
}

/// @private
template <typename T>
T distance_norm2(herm_matrix_timestep_view<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == g2.tstp());
    herm_matrix_timestep_view<T> g2_view(g1.tstp(), g2);
    if (g1.size1() == 1)
        return distance_norm2_dispatch<T, 1>(g1.tstp(), g1, g2_view);
    else
        return distance_norm2_dispatch<T, LARGESIZE>(g1.tstp(), g1, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between  a `herm_matrix_timestep_view` and `herm_matrix`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between  of a `herm_matrix_timestep_view` \f$g_1\f$ and a `herm_matrix` (\f$g_2\f$) at a given time step (tstp).
 * > To evaluate the norm, the elements of retarded, lesser and left-mixing components at the time step is used. The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep_view
 * @param g2
 * > herm_matrix
 */
template <typename T>
T distance_norm2(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.nt() >= tstp);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_dispatch<T, 1>(tstp, g1, g2_view);
    else
        return distance_norm2_dispatch<T, LARGESIZE>(tstp, g1, g2_view);
}

/** \brief <b>  Evaluate the Euclidean norm between  a `herm_matrix` and `herm_matrix_timestep_view`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between  of a `herm_matrix` \f$g_1\f$ and a `herm_matrix_timestep_view` (\f$g_2\f$) at a given time step (tstp).
 * > To evaluate the norm, the elements of retarded, lesser and left-mixing components at the time step is used. The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix
 * @param g2
 * > herm_matrix_timestep_view
 */
template <typename T>
T distance_norm2(int tstp, herm_matrix<T> &g1, herm_matrix_timestep_view<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.nt() >= tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    if (g1.size1() == 1)
        return distance_norm2_dispatch<T, 1>(tstp, g1_view, g2);
    else
        return distance_norm2_dispatch<T, LARGESIZE>(tstp, g1_view, g2);
}


/** \brief <b>  Evaluate the Euclidean norm between  a `herm_matrix_timestep` and `herm_matrix_timestep_view`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between  of a `herm_matrix_timestep` \f$g_1\f$ and a `herm_matrix_timestep_view` (\f$g_2\f$) at a given time step (tstp).
 * > To evaluate the norm, the elements of retarded, lesser and left-mixing components at the time step is used. The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep
 * @param g2
 * > herm_matrix_timestep_view
 */
template <typename T>
T distance_norm2(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep_view<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g1_view(tstp, g1);
    if (g1.size1() == 1)
        return distance_norm2_dispatch<T, 1>(tstp, g1_view, g2);
    else
        return distance_norm2_dispatch<T, LARGESIZE>(tstp, g1_view, g2);
}

/** \brief <b>  Evaluate the Euclidean norm between  a `herm_matrix_timestep_view` and `herm_matrix_timestep`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between  of a `herm_matrix_timestep_view` \f$g_1\f$ and a `herm_matrix_timestep` (\f$g_2\f$) at a given time step (tstp).
 * > To evaluate the norm, the elements of retarded, lesser and left-mixing components at the time step is used. The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep_view
 * @param g2
 * > herm_matrix_timestep
 */
template <typename T>
T distance_norm2(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    herm_matrix_timestep_view<T> g2_view(tstp, g2);
    if (g1.size1() == 1)
        return distance_norm2_dispatch<T, 1>(tstp, g1, g2_view);
    else
        return distance_norm2_dispatch<T, LARGESIZE>(tstp, g1, g2_view);
}


/** \brief <b>  Evaluate the Euclidean norm between  a `herm_matrix_timestep_view` and `herm_matrix_timestep_view`  at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between  of a `herm_matrix_timestep_view` \f$g_1\f$ and a `herm_matrix_timestep_view` (\f$g_2\f$) at a given time step (tstp).
 * > To evaluate the norm, the elements of retarded, lesser and left-mixing components at the time step is used. The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > herm_matrix_timestep_view
 * @param g2
 * > herm_matrix_timestep_view
 */
template <typename T>
T distance_norm2(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep_view<T> &g2) {
    assert(g1.size1() == g2.size1());
    assert(g1.ntau() == g2.ntau());
    assert(g1.tstp() == tstp);
    assert(g2.tstp() == tstp);
    if (g1.size1() == 1)
        return distance_norm2_dispatch<T, 1>(tstp, g1, g2);
    else
        return distance_norm2_dispatch<T, LARGESIZE>(tstp, g1, g2);
}


// template <typename T>
// T distance_norm2(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix<T> &g2) {
//     herm_matrix_timestep<T> tmp(tstp,g2.ntau(),g2.size1(),g2.size2(),g2.sig());
//     g1.get_data(tmp);
//     return distance_norm2(tstp,tmp,g2);
// }

// template <typename T>
// T distance_norm2(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep<T> &g2) {
//     herm_matrix_timestep<T> tmp(tstp,g2.ntau(),g2.size1(),g2.size2(),g2.sig());
//     g1.get_data(tmp);
//     return distance_norm2(tmp,g2);
// }

template <typename T>
T distance_norm2(int tstp, herm_pseudo<T> &g1, herm_pseudo<T> &g2) {
    if (g1.size1() == 1)
        return distance_norm2_dispatch<T, herm_pseudo<T>, 1>(tstp, g1, g2);
    else
        return distance_norm2_dispatch<T, herm_pseudo<T>, LARGESIZE>(tstp, g1,
                                                                     g2);
}

/** \brief <b>  Evaluate the Euclidean norm between 'herm_matrix_timestep' and 'herm_matrix' at a given time step. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the Euclidean norm between 'herm_matrix_timestep' (\f$g_1 \f$) and 'herm_matrix' (\f$ g_2\f$)
 * > at a given time step 'tstp' using the standard 'norm'-routine.
 * > To evaluate the norm, the elements of retarded, lesser and left-mixing components at the time step is used.
 * > The norm is not normalized per elements, but it is the summention of all the elements.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > time step
 * @param g1
 * > object of the class 'herm_matrix_timestep'
 * @param g2
 * > object of the class 'herm_matrix'
 */
template <typename T>
T distance_norm2_eigen(int tstp,herm_matrix_timestep<T> &g1,herm_matrix<T> &g2){
   int size1=g2.size1(),size2=g2.size2(),ntau=g2.ntau();
   T err=0.0;

   assert(g1.tstp_==tstp && g1.tstp_<=g2.nt());
   assert(g1.ntau_==g2.ntau() && g1.size1_==g2.size1() && g1.size2_==g2.size2() );
   cdmatrix matg1(size1,size2);
   cdmatrix matg2(size1,size2);

   if(tstp==-1){
	 for(int i=0;i<=ntau;i++){
	    g1.get_mat(i,matg1);
	    g2.get_mat(i,matg2);
	    err+=(matg1-matg2).norm();
	 }
   }else{
   	 // Ret
	 for(int i=0;i<=tstp;i++){
	 	g1.get_ret_tstp_t(i,matg1);
	    g2.get_ret(tstp,i,matg2);
	    err+=(matg1-matg2).norm();
	 }
	 // tv
	 for(int i=0;i<=ntau;i++){
	 	g1.get_tv(i,matg1);
	    g2.get_tv(tstp,i,matg2);
	    err+=(matg1-matg2).norm();
	 }
	 // Les
	 for(int i=0;i<=tstp;i++){
	  	g1.get_les_t_tstp(i,matg1);
	    g2.get_les(i,tstp,matg2);
	    err+=(matg1-matg2).norm();
	 }
   }
}

/** \brief <b> Evaluate the memory necessary for a herm_matrix. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 *  > Evaluate the memory necessary for a herm_matrix.
 *  > Square matrix is assumed.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 * @param nt
 * > maximum number of the time steps
 * @param ntau
 * > number of time grids on the Matsubara axis
 * @param size
 * > size of the colums and the rows of the square Matrix
 */
template <typename T>
size_t mem_herm_matrix(int nt, int ntau, int size) {
    size_t mem = ((nt + 1) * (nt + 2) + (nt + 2) * (ntau + 1)) * size * size *
                 sizeof(std::complex<T>);
    return mem;
}

/** \brief <b> Evaluate the memory necessary for a herm_matrix_timestep. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 * > Evaluate the memory necessary for a herm_matrix_timestep. Squre matrix is assumed.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 * @param tstp
 * > time step
 * @param ntau
 * > number of time grids on the Matsubara axis
 * @param size
 * > size of the colums and the rows of the square Matrix
 */
template <typename T>
size_t mem_herm_matrix_timestep(int tstp, int ntau, int size) {
    size_t mem = ((tstp + 2) + (ntau + 1)) * size * size *
                 sizeof(std::complex<T>);
    return mem;
}

/** \brief <b> Evaluate the memory necessary for a contour function. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 *  > Evaluate the memory necessary for a contour function. Squre matrix is assumed.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 * @param nt
 * > maximum number of the time steps
 * @param size
 * > size of the colums and the rows of the square Matrix
 */
template <typename T>
size_t mem_function(int nt, int size) {
    size_t mem = (nt + 2) * size * size * sizeof(std::complex<T>);
    return mem;
}

/** \brief <b> Force the Matsubara component of herm_matrix to be a hermitian matrix at each \f$ \tau \f$. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 *  > Force the Matsubara component of herm_matrix to be a hermitian matrix at each \f$ \tau \f$ by \f$ \hat{G}^M(\tau)=\frac{1}{2} \bigl\{\hat{G}^M(\tau)+[\hat{G}^{M}(\tau) ]^{\dagger} \bigl\}.\f$
 *
 * <!-- ARGUMENTS
 *      ========= -->
 * @param G
 * > herm_matrix to be modulated
 */
template <typename T>
void force_matsubara_hermitian(herm_matrix<T> &G) {
  int ntau = G.ntau(), m, s1 = G.size1(), p1, p2;
  cdmatrix Gmat,Gmat_herm;
  for  (m = 0; m <= ntau; m++) {
      G.get_mat(m,Gmat);
      Gmat_herm = 0.5*(Gmat + Gmat.adjoint());
      G.set_mat(m,Gmat_herm);
    }
}

/** \brief <b> Force the Matsubara component of a two-time contour object to be a hermitian matrix at each \f$ \tau \f$. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *  \par Purpose
 * <!-- ========= -->
 *  > Force the Matsubara component of herm_matrix to be a two-time contour object  at each \f$ \tau \f$ by \f$ \hat{G}^M(\tau)=\frac{1}{2} \bigl\{\hat{G}^M(\tau)+[\hat{G}^{M}(\tau) ]^{\dagger} \bigl\}.\f$
 *
 * <!-- ARGUMENTS
 *      ========= -->
 * @param G
 * > [GG] a two-time contour object to be modified
 */
template <class GG>
void force_matsubara_hermitian(GG &G) {
    int ntau = G.ntau(), m, s1 = G.size1(), p1, p2;
    cdmatrix Gmat,Gmat_herm;
    for  (m = 0; m <= ntau; m++) {
      G.get_mat(m,Gmat);
      Gmat_herm = 0.5*(Gmat + Gmat.adjoint());
      G.set_mat(m,Gmat_herm);
    }
}

} // namespace cntr

#endif  // CNTR_UTILITIES_IMPL_H
