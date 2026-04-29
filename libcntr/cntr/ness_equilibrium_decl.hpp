#ifndef NESS_EQUILIBRIUM_DECL_H
#define NESS_EQUILIBRIUM_DECL_H

#include "ness_global_settings.hpp"
#include "cntr/cntr.hpp"

using namespace std;
using namespace cntr;

namespace ness {

// Has to be here in decl file because they are template functions
/// @private
template <class dos_function>
void green_equilibrium_ret_ness(GF &G,dos_function &dos,double h,int limit,int nn,double mu)
{
    typedef std::complex<double> cplx;
    int nt=G.ngrid_,l,size1=G.size1_;
    int sign=G.sign_;
    double t;
    cplx res,err,cplx_i=cplx(0,1);
    fourier::adft_func adft;
    dos_wrapper<dos_function> dos1(dos,sign,mu);
    dos1.mu_=mu;
    dos1.x_=ret;
    adft.sample(0.0,dos.lo_,dos.hi_,dos1,nn,limit);
    for(int idx1 = 0; idx1 < size1; idx1++)
    {
        for(l=0;l<nt;l++){ // l = t-t'
            t=h*l;
            adft.dft(-t,res,err);
            res *= std::complex<double>(0,-1.0);
            (*G.p_ret(l, idx1, idx1))=res;
        }
    }
}

/// @private
template <class dos_function>
void green_equilibrium_les_ness(GF &G,dos_function &dos,double beta,double h,int limit,int nn,double mu)
{
    typedef std::complex<double> cplx;
    int nt=G.ngrid_,l,size1=G.size1_;
    int sign=G.sign_;
    double t;
    cplx res,err,cplx_i=cplx(0,1);
    fourier::adft_func adft;
    dos_wrapper<dos_function> dos1(dos,sign,mu);
    dos1.mu_=mu;
    dos1.x_=les;
    dos1.beta_=beta;
    adft.sample(0.0,dos.lo_-0.0,dos.hi_-0.0,dos1,nn,limit);
    for(int idx1 = 0; idx1 < size1; idx1++)
    {
        for(l=0;l<nt;l++){ // l = t-t'
            t=h*l;
            adft.dft(t,res,err);
            res *= std::complex<double>(0,-1.0);
            (*G.p_les(l, idx1, idx1))=res;
        }
    }
}

/** \brief <b> Equilibrium steady state propagator for the given density of states  </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Calculate the steady state equilibrium `GF` \f$ G \f$ for the given density of states (own type `dos_function`)
 * > via \f$G(t-t') = -i \int d\omega A(\omega+\mu) \exp(i\omega (t'-t)) [ \Theta(t-t') - \Theta(t'-t)\exp(-\beta*\omega)]\f$,
 * for all times and the retarded and lesser component. The calculation uses an adaptive Fourier transform algorithm.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param G
 *  The output Greens function GF set to the equilibrium free propagator
 * @param dos
 *  density of states
 * @param beta
 *  inverse temperature
 * @param h
 *  timestep
 * @param mu
 *  chemical potential
 * @param limit
 *  max number of intervals in Fourier transform (default: 100)
 * @param nn
 *  number of points in each interval of the Fourier transform (default: 20)
 */
template <class dos_function>
void green_equilibrium_ness(GF &G, dos_function &dos,double beta,double h,int limit,int nn,double mu)
{
    green_equilibrium_ret_ness(G,dos,h,limit,nn,mu);
    green_equilibrium_les_ness(G,dos,beta,h,limit,nn,mu);
}

/// @private
template <class dos_function>
void green_equilibrium_ret_ness(GF_pair &G,dos_function &dos,double h,int limit,int nn,double mu)
{
    typedef std::complex<double> cplx;
    int nt=G.time.ngrid_,l,size1=G.time.size1_;
    int sign=G.time.sign_;
    double t;
    cplx res,err,cplx_i=cplx(0,1);
    fourier::adft_func adft;
    dos_wrapper<dos_function> dos1(dos,sign,mu);
    dos1.mu_=mu;
    dos1.x_=ret;
    adft.sample(0.0,dos.lo_,dos.hi_,dos1,nn,limit);
    for(int idx1 = 0; idx1 < size1; idx1++)
    {
        for(l=0;l<nt;l++){ // l = t-t'
            t=h*l;
            adft.dft(-t,res,err);
            res *= std::complex<double>(0,-1.0);
            (*G.time.p_ret(l, idx1, idx1))=res;
        }
    }
}

/// @private
template <class dos_function>
void green_equilibrium_les_ness(GF_pair &G,dos_function &dos,double beta,double h,int limit,int nn,double mu)
{
    typedef std::complex<double> cplx;
    int nt=G.time.ngrid_,l,size1=G.time.size1_;
    int sign=G.time.sign_;
    double t;
    cplx res,err,cplx_i=cplx(0,1);
    fourier::adft_func adft;
    dos_wrapper<dos_function> dos1(dos,sign,mu);
    dos1.mu_=mu;
    dos1.x_=les;
    dos1.beta_=beta;
    adft.sample(0.0,dos.lo_-0.0,dos.hi_-0.0,dos1,nn,limit);
    for(int idx1 = 0; idx1 < size1; idx1++)
    {
        for(l=0;l<nt;l++){ // l = t-t'
            t=h*l;
            adft.dft(-t,res,err);
            res *= std::complex<double>(0,-1.0);
            (*G.time.p_les(l, idx1, idx1))=res;
        }
    }
    // flip to negatie times.
}

/** \brief <b> Equilibrium steady state propagator for the given density of states  </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Calculate the steady state equilibrium `GF_pair` \f$ G \f$ (in time domain) for the given density of states (own type `dos_function`)
 * > via \f$G(t-t') = -i \int d\omega A(\omega+\mu) \exp(i\omega (t'-t)) [ \Theta(t-t') - \Theta(t'-t)\exp(-\beta*\omega)]\f$,
 * for all times and the retarded and lesser component. The calculation uses an adaptive Fourier transform algorithm.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param G
 *  The output Greens function pair GF_pair set to the equilibrium free propagator (time component)
 * @param dos
 *  density of states
 * @param beta
 *  inverse temperature
 * @param h
 *  timestep
 * @param mu
 *  chemical potential
 * @param limit
 *  max number of intervals in Fourier transform (default: 100)
 * @param nn
 *  number of points in each interval of the Fourier transform (default: 20)
 */
template <class dos_function>
void green_equilibrium_ness(GF_pair &G, dos_function &dos,double beta,double h,int limit,int nn,double mu)
{
    green_equilibrium_ret_ness(G,dos,h,limit,nn,mu);
    green_equilibrium_les_ness(G,dos,beta,h,limit,nn,mu);
}


/// @private
/** \brief <b> Class `gauss` represent the gaussian density of states </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > This class contains the data structures for representing the density of states
 * > for a gaussian with mean `mu_=0` and a standard deviation of `sigma_=1`.
 *
 */
class gauss {
public:
    /** \brief <b> Higher edge of the density of states </b> */
    double hi_;
    /** \brief <b> Lower edge of the density of states </b> */
    double lo_;
    /** \brief <b> Mean of gaussian </b> */
    double mu_;
    /** \brief <b> Standard deviation of gaussian </b> */
    double sigma_;
    gauss();
    gauss(double mu, double sigma);
    double operator()(double x);
};


/// @private
/** \brief <b> Class `ohmic` represents a symmetric ohmic bath \f$ x^2 \cdot exp(-x/x_c) \f$ density of states </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > This class contains the data structures for representing the density of states
 * > for an ohmic bath and includes a free parameter: the cutoff energy  \f$ \omega_c \f$
 * > and is adapted so the spectrum is antisymmetric (boson).
 *
 */
class ohmic_sym{
public:
    /** \brief <b> Higher edge of the density of states </b> */
    double hi_;
    /** \brief <b> Lower edge of the density of states </b> */
    double lo_;
    /** \brief <b> Cutoff </b> */
    double omegac_;
    ohmic_sym();
    ohmic_sym(double omegac);
    double operator()(double x);
};


} // namespace ness

#endif  // NESS_EQUILIBRIUM_DECL_H
