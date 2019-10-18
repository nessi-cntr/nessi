#ifndef CNTR_MATSUBARA_IMPL_H
#define CNTR_MATSUBARA_IMPL_H

#include "cntr_matsubara_decl.hpp"
#include "fourier.hpp"
#include "cntr_elements.hpp"

namespace cntr {

// ----------------------------------------------------------------------

// Helper functions for imaginary time and matsubara frequencies 
/// @private   
template <typename F, typename I> F get_tau(const I tau_idx, const F beta, const I ntau) {
    return tau_idx * (beta / ntau);
}

/// @private
template <typename F, typename I> F get_omega(const I m, const F beta, const I sig) {

    // matsub_one = 1 if sig = -1, 0 if sig = 1
    // (shift is integer division by 2)
    const I matsub_one = ((1 - sig) >> 1);
    // const I matsub_one = (sig == -1 ? 1.0 : 0.0);

    // return fourier::pi/beta * (2*m + matsub_one);
    return PI / beta * (2 * m + matsub_one);
}

// ----------------------------------------------------------------------
/// @private
/** \brief <b> Returns the first-order tail correction term \f$c_0(\tau)\f$, such that
    \f$\tilde{C}^\mathrm{M}(\tau) = C^\mathrm{M}(\tau) + c_0(\tau)\f$ becomes continuous. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > The Matsubara function \f$ C^\mathrm{M}(\tau) \f$ is generally discontinuous at \f$\tau=0\f$ and
* > \f$\tau=\beta\f$; the modified function \f$\tilde{C}^\mathrm{M}(\tau) = C^\mathrm{M}(\tau) + c_0(\tau)\f$, 
* > where \f$c_0(\tau)\f$ is the tail correction, is then a continuous function. Thus, the Fourier series
* > coefficients \f$ C^\mathrm{M}(i\omega_m) \f$ tend to zero as \f$(i\omega)^{-2}\f$ 
* > instead of \f$(i\omega)^{-1}\f$. 
* >
* > For fermions, one finds \f$c_0(\tau) = -a_0/2\f$, while for bosons \f$c_0(\tau)= a_0(\tau/\beta-1/2)\f$.
* > The coefficient \f$a_0\f$ is given by the jump at \f$\tau=0\f$: 
* > \f$a_0=-(C^\mathrm{M}(0)+C^\mathrm{M}(\beta))\f$ for fermions, \f$a_0=C^\mathrm{M}(\beta)-C^\mathrm{M}(0)\f$
* > for bosons.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param xmat
* > [complex] on return, pointer storing the correction \f$c_0(\tau)\f$
* @param coeff
* > [complex] coefficient \f$a_0\f$ (square matrix in general)
* @param beta
* > [double] inverse temperature
* @param sg
* > [int] element size
* @param ntau
* > [int] number of points on Matsubara track
* @param sig
* > [int] `sig=-1` for fermions, `sig=+1` for bosons
* @param size1
* > [int] matrix dimension
*/
template <typename T, int SIZE1>
void set_first_order_tail(std::complex<T> *xmat, std::complex<T> *coeff, T beta, int sg,
                          int ntau, int sig, int size1) {

    // First order tail correction tau-polynomial for Bosons & Fermions
    // Fermions : g_0(\tau) = -1/2
    // Bosons   : g_0(\tau) = -1/2 + \tau/\beta

    std::complex<T> *z1 = new std::complex<T>[sg];

    for (int r = 0; r <= ntau; r++) {
        // -1/2 factor
        element_set<T, SIZE1>(size1, xmat + r * sg, coeff);
        element_smul<T, SIZE1>(size1, xmat + r * sg, -0.5);

        if (sig == 1) {
            // tau/beta factor
            double tau = r * (beta / ntau);
            element_set<T, SIZE1>(size1, z1, coeff);
            element_smul<T, SIZE1>(size1, z1, tau / beta);
            element_incr<T, SIZE1>(size1, xmat + r * sg, z1);
        }
    }

    delete[] z1;
}

// ----------------------------------------------------------------------

/*-------------------------------------------------------
compute (WITHOUT FFT)
I=sum_{r=0}^m ftau[r] exp(i omega_n tau_r)
for n=0...m
tau_r=r*dtau,tau=beta/m
omega_n= (2n+1)pi/beta, (sgn=-1)
omega_n= 2n*pi/beta,    (sgn=+1)

I=sum_{r=0}^m ftau[r] exp(i pi (2n+1)r/m)
-------------------------------------------------------*/


/** \brief <b> Computes the Fourier series coefficients of a Matsubara function by plain DFT. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Computes \f$ I = \sum_{r=0}^m f(\tau_r) e^{i \omega_n \tau_r}\f$ for \f$n=0,\dot,m\f$, assuming
* > \f$\tau_r = r \beta/m\f$. The Matsubara frequencies are given by \f$\omega_n= (2n+1)\pi/\beta\f$ 
* > for fermions, while \f$\omega_n= 2n*\pi/\beta\f$ for bosons.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param mdft
* > [complex] on return, pointer storing the Fourier coeffients in element representation 
* @param G
* > [GG] Matsubara function to be Fourier transformed
* @param sig
* > [int] `sig=-1` for fermions, `sig=+1` for bosons
*/
template <typename T, class GG, int SIZE1>
void matsubara_dft(std::complex<T> *mdft, GG &G, int sig) {
    typedef std::complex<T> cplx;
    int ntau, r, m, sg, size1 = G.size1();
    double arg, one;
    cplx *z, *z1, expfac;
    sg = G.element_size();
    ntau = G.ntau();
    z = new cplx[sg];
    z1 = new cplx[sg];
    one = (sig == -1 ? 1.0 : 0.0);
    for (m = 0; m <= ntau; m++) {
        element_set_zero<T, SIZE1>(size1, z);
        for (r = 0; r <= ntau; r++) {
            arg = ((2 * m + one) * r * PI) / ntau;
            expfac = cplx(cos(arg), sin(arg));
            element_set<T, SIZE1>(size1, z1, G.matptr(r));
            element_smul<T, SIZE1>(size1, z1, expfac);
            element_incr<T, SIZE1>(size1, z, z1);
        }
        element_set<T, SIZE1>(size1, mdft + m * sg, z);
    }
    delete[] z;
    delete[] z1;
}

/** \brief <b> Computes the Fourier series coefficients of a Matsubara 
    function by cubically corrected DFT. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Evaluates the Fourier integral 
* > \f$G^\mathrm{M}(i \omega_m) = \int^\beta_0 d\tau\,e^{i \omega_m \tau} G^\mathrm{M}(\tau) \f$ 
* > by linearly or cubically corrected DFT.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param result
* > [complex] on return, pointer storing the Fourier coeffient \f$G^\mathrm{M}(i \omega_m)\f$
* @param m
* > [int] index of Matsubara frequency \f$\omega_m\f$
* @param G
* > [GG] Matsubara function to be Fourier transformed
* @param mdft
* > [complex] pointer storing the Fourier coeffients computed by `matsubara_dft` in element representation 
* @param sig
* > [int] `sig=-1` for fermions, `sig=+1` for bosons
* @param beta
* > [double] inverse temperature
* @param order
* > [int] `order=1` / `order=3` for linearly/cubically corrected DFT
*/
template <typename T, class GG, int SIZE1>
void matsubara_ft(std::complex<T> *result, int m, GG &G, std::complex<T> *mdft, int sig,
                  T beta, int order) {
    typedef std::complex<T> cplx;
    double corfac, arg, dtau, one;
    int ntau, m1, r, l, sg, size1 = G.size1();
    cplx *z1, *z2, *dft, bcorr[4];

    one = (sig == -1 ? 1.0 : 0.0);
    ntau = G.ntau();
    sg = G.element_size();
    dtau = beta / ntau;
    z1 = new cplx[sg];
    z2 = new cplx[sg];
    dft = new cplx[sg];

    arg = (2 * m + one) * PI / ntau; /*arg=omn*dtau, omn=(2m+1)Pi/beta*/
    if (order == 0) {
        corfac = 1.0;
        bcorr[0] = 0.0;
    } else if (order == 1) {
        fourier::get_dftcorr_linear(arg, &corfac, bcorr);
    } else if (order == 3) {
        fourier::get_dftcorr_cubic(arg, &corfac, bcorr);
    } else {
        std::cerr << "matsubara_ft: unknown order" << std::endl;
        abort();
    }
    m1 = m; /*shift m1 to interval [0,ntau-1]*/
    while (m1 < 0)
        m1 += ntau;
    while (m1 >= ntau)
        m1 -= ntau;
    for (l = 0; l < sg; l++)
        dft[l] = mdft[m1 * sg + l] * corfac;
    element_set_zero<T, SIZE1>(size1, z1);
    for (r = 0; r <= order; r++) {
        element_set<T, SIZE1>(size1, z2, G.matptr(r));
        element_smul<T, SIZE1>(size1, z2, bcorr[r]);
        element_incr<T, SIZE1>(size1, z1, z2);
        element_set<T, SIZE1>(size1, z2, G.matptr(ntau - r));
        // element_smul<T,SIZE1>(size1,z2,-cplx(bcorr[r].real(),-bcorr[r].imag()));
        element_smul<T, SIZE1>(size1, z2,
                               (1. * sig) * cplx(bcorr[r].real(), -bcorr[r].imag()));
        element_incr<T, SIZE1>(size1, z1, z2);
    }
    for (l = 0; l < sg; l++)
        result[l] = dtau * (z1[l] + dft[l]);

    delete[] z1;
    delete[] z2;
    delete[] dft;
}

// ----------------------------------------------------------------------

}  // namespace cntr

#endif  // CNTR_MATSUBARA_IMPL_H
