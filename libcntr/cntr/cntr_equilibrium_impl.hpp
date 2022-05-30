#ifndef CNTR_EQUILIBRIUM_IMPL_H
#define CNTR_EQUILIBRIUM_IMPL_H

#include "cntr_global_settings.hpp"
#include "cntr_equilibrium_decl.hpp"
#include "fourier.hpp"
#include "integration.hpp"
#include "cntr_elements.hpp"
//#include "cntr_exception.hpp"
#include "cntr_herm_matrix_decl.hpp"
#include "cntr_herm_matrix_timestep_decl.hpp"
#include "cntr_herm_pseudo_decl.hpp"
#include "cntr_function_decl.hpp"
#include "cntr_utilities_decl.hpp"

namespace cntr {

/*####################################################################################
#
#   The computation of the equilibrium propagator from the density of states
#
######################################################################################*/

#define EXPMAX 100
/** \brief <b> Evaluates the Fermi distribution function. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Evaluates the Fermi distribution function \f$ n_\mathrm{F}(\omega)\f$ at given inverse temperature
* > \f$\beta\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param beta
* > [T] inverse temperature \f$\beta\f$
* @param omega
* > [T] energy
*/
template <typename T>
T fermi(T beta, T omega) {
    T arg = omega * beta;
    if (fabs(arg) > EXPMAX) {
        return (arg > 0.0 ? 0.0 : 1.0);
    } else {
        return 1.0 / (1.0 + exp(arg));
    }
}
/** \brief <b> Evaluates the Fermi distribution function for multiple energies. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Evaluates the Fermi distribution function \f$ n_\mathrm{F}(\omega)\f$ at given inverse temperature
* > \f$\beta\f$ for a vector of energies.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param beta
* > [T] inverse temperature \f$\beta\f$
* @param omega
* > [dvector] vector of energies
*/
template<typename T>
dvector fermi(T beta,dvector &omega){
   int size=omega.size();
   dvector tmp(size);
   for(int i=0;i<size;i++){
      tmp(i)=fermi(beta,omega(i));
   }
   return tmp;
}
/// @private
template <typename T>
T fermi_exp(T beta, T tau, T omega) {
    if (omega < 0) {
        return exp(omega * tau) *
               fermi(beta, omega); // exp(w*t)/(1+exp(b*w)) always OK for w<0
    } else {
        return exp((tau - beta) * omega) *
               fermi(beta, -omega); // exp((t-b)*w)/(1+exp(-w*b))
    }
}
/// @private
template<typename T>
dvector fermi_exp(T beta,T tau,dvector &omega){
   int size=omega.size();
   dvector tmp(size);

   for(int i=0;i<size;i++){
      tmp(i)=fermi_exp(beta,tau,omega(i));
   }
   return tmp;
}
/// @private
template<typename T>
cdmatrix diag_prop(T time,dvector &omega){
   int size=omega.size();
   cdmatrix tmp(size,size);
   for(int i=0;i<size;i++){
      tmp(i,i)=std::complex<T>(cos(omega(i)*time),sin(omega(i)*time));
   }
   return tmp;
}

/** \brief <b> Evaluates the Bose distribution function. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Evaluates the Bose distribution function \f$ n_\mathrm{B}(\omega)\f$ at given inverse temperature
* > \f$\beta\f$. The calculation is numerically stable for any \f$\omega\ne 0 \f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param beta
* > [T] inverse temperature \f$\beta\f$
* @param omega
* > [T] energy
*/
template <typename T>
T bose(T beta, T omega) {
    T arg = omega * beta;
    if (arg < 0)
        return (-1.0 - bose<T>(beta, -omega));
    if (fabs(arg) > EXPMAX) {
        return 0.0;
    } else if (arg < 1e-10) {
        return 1.0 / arg;
    } else {
        return 1.0 / (exp(arg) - 1.0);
    }
}

/** \brief <b> Evaluates the Bose distribution function. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > Evaluates the Bose distribution function \f$ n_\mathrm{B}(\omega)\f$ at given inverse temperature
* > \f$\beta\f$ for a vector of energies. The calculation is numerically stable for any \f$\omega\ne 0 \f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param beta
* > [T] inverse temperature \f$\beta\f$
* @param omega
* > [dvector] vector of energies
*/
template<typename T>
dvector bose(T beta,dvector &omega){
   int size=omega.size();
   dvector tmp(size);
   for(int i=0;i<size;i++){
      tmp(i)=bose(beta,omega(i));
   }
   return tmp;
}


// exp(tau*w)b(w) ... actually assuming positive tau
/// @private
template <typename T>
T bose_exp(T beta, T tau, T omega) {
    if (omega < 0)
        return exp(tau * omega) * bose<T>(beta, omega);
    else
        return -exp((tau - beta) * omega) * bose<T>(beta, -omega);
}
/// @private
template<typename T>
dvector bose_exp(T beta,T tau,dvector &omega){
   int size=omega.size();
   dvector tmp(size);

   for(int i=0;i<size;i++){
      tmp(i)=bose_exp(beta,tau,omega(i));
   }
   return tmp;
}


#undef EXPMAX

#define EXPMAX 100
/// @private
template <typename T>
double green_cntr_full_equilibrium_fermi(T beta, T omega) {
    T arg = omega * beta;
    if (fabs(arg) > EXPMAX) {
        return (arg > 0.0 ? 0.0 : 1.0);
    } else {
        return 1.0 / (1.0 + exp(arg));
    }
}
/// @private
template<typename T>
double distribution_eq(T beta,T omega,int sign)
{
  if(sign==-1){
    return green_cntr_full_equilibrium_fermi(beta,omega);
  }else if(sign==1){
    return bose(beta,omega);
  }
}
/// @private
template<typename T>
double distribution_exp_eq(T beta,T tau,T omega,int sign)
{
  if(sign==-1){
    return fermi_exp(beta,tau,omega);
  }else if(sign==1){
    return bose_exp(beta,tau,omega);
  }
}


#undef EXPMAX

/// @private
enum flavor { ret, adv, mat, tv, vt, les, gtr };

template < class dos > struct dos_wrapper{
    flavor x_;
    double beta_;
    double tau_;
    int sign_;
    double mu_;
    dos A_;
    // init with a dos object:
    dos_wrapper(const dos& f,int sign,double mu=0.0){
      beta_=0;
      tau_=0;
      A_=f;
      mu_=mu;
      sign_=sign;
    }
    double operator()(double omega){
      double a=A_(omega),fm;
      if(x_==ret){
          fm=1;
      }else if(x_==les){
        fm=sign_*distribution_eq(beta_,omega-mu_,sign_);
      }else if (x_==tv){
        fm=sign_*distribution_exp_eq(beta_,tau_,omega-mu_,sign_);
      }else if(x_==mat){
        if(omega-mu_>0){
           fm=sign_*exp(-tau_*(omega-mu_))*distribution_eq(beta_,-(omega-mu_),sign_);
        }else{
           fm=(-1.0)*exp((beta_-tau_)*(omega-mu_))*distribution_eq(beta_,(omega-mu_),sign_);
        }
      }else{
        fm=0;
      }
      return a*fm;
    }
};


/*####################################################################################
#
#   Here is the definition of some dos functions:
#
#   The minimal functionality of the dos_function class is an
#   operator()(double) which returns the values of the dos at some
#   previously set parameters (such as the bandwidth etc.) and
#   members lo_, hi_
#
######################################################################################*/

/** \brief <b> Class `bethe` represent the semicircular density of states </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > This class contains the data structures for representing the density of states
 * > for a Bethe semicircular with the bandwidth 4V.
 *
 */
class bethedos {
  public:
    /** \brief <b> Higher edge of the density of states </b> */
    double hi_;
    /** \brief <b> Lower edge of the density of states </b> */
    double lo_;
    /** \brief <b> Hopping integral and 4V corresponds to the bandwidth </b> */
    double V_;

    double Eshift_;
    bethedos() {
        Eshift_ = 0.0;
        V_ = 1;
        lo_ = -2;
        hi_ = 2;
    }
    bethedos(double Eshift) {
        Eshift_ = Eshift;
        V_ = 1;
        lo_ = -2 + Eshift_;
        hi_ = 2 + Eshift_;
    }
    double operator()(double x) {
        double xv = x - Eshift_;
        double arg = 4.0 * V_ * V_ - xv * xv;
        double num = V_ * V_ * 3.14159265358979323846 * 2;
        return (arg < 0 ? 0.0 : sqrt(arg) / num);
    }
};
/// @private
/** \brief <b> Class 'ohmic` represent a ohmic bath x^2*exp(-x/x_c) density of states </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > This class contains the data structures for representing the density of states
 * > for an ohmic bath and includes a free parameter: the cutoff energy  \f$ \omega_c \f$
 *
 */
class ohmic{
  public:
  	/** \brief <b> Higher edge of the density of states </b> */
    double hi_;
    /** \brief <b> Lower edge of the density of states </b> */
    double lo_;
    /** \brief <b> Cutoff </b> */
    double omegac_;
    ohmic() {
        omegac_ = 1;
        lo_ = 0.01;
        hi_ = omegac_*20;
    }
    ohmic(double omegac) {
        omegac_ = omegac;
        lo_=0.01;
        hi_=omegac*20; // value is 10^-6
    }
    double operator()(double x) {
        return x*x*std::exp(-x/omegac_)/(2.0*omegac_*omegac_*omegac_);
    }
};


/** \brief <b> Class 'smooth_box` represent a smooth box-like density of states </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > This class contains the data structures for representing the density of states
 * > for a smooth box, which is modeled by two overlaping Fermi functions and the
 * > smoothness parameter \f$ \nu \f$ corresponds to the inverse temperature.
 *
 */
class smooth_box {
  public:
    /** \brief <b> Internal value determining the higher edge of the Fourier transform region  \f$ </b> */
    double hi_;
    /** \brief <b> Internal value determining the lower edge of the Fourier transform region  \f$ </b> */
    double lo_;
    /** \brief <b> Lower border of the box \f$ </b> */
    double A_;
    /** \brief <b> Higher border of the box \f$ </b> */
    double B_;
    /** \brief <b> Smoothness parameter of the box \f$ </b> */
    double nu_;
    /** \brief <b> Normalization of density of states \f$ </b> */
    double N_;
    smooth_box() {
        A_ = -1.0;
        B_ = 1.0;
        nu_ = 10.0;
        N_ = -1.0;
        lo_ = A_ - 20.0 / nu_;
        hi_ = B_ + 20.0 / nu_;
    }
    smooth_box(double a, double b, double nu) {
        A_ = a;
        B_ = b;
        nu_ = nu;
        N_ = -1.0;
        if (nu_ > 100 || nu_ <= 0.1 || A_ >= B_) {
            std::cerr << "nu=" << nu_ << " in smooth_box " << std::endl;
            exit(0);
        }
        lo_ = A_ - 20.0 / nu_;
        hi_ = B_ + 20.0 / nu_;
    }
    void set_norm(void) {
        typedef std::complex<double> cplx;
        cplx res, err;
        N_ = 1.0;
        fourier::adft_func adft;
        adft.sample(0.0, lo_, hi_, *this, 20, 100);
        adft.dft(0.0, res, err);
        N_ = res.real();
    }
    double operator()(double x) {
        if (N_ == -1)
            set_norm();
        double arg1 = fermi<double>(nu_, x - B_);
        double arg2 = fermi<double>(nu_, A_ - x);
        return arg1 * arg2 / N_;
    }
};


/// @private
template <typename T,class dos_function>
void green_equilibrium_ret(herm_matrix<T> &G,dos_function &dos,double h,int limit,int nn,double mu)
{
  typedef std::complex<double> cplx;
  int nt=G.nt(),i,l,size1=G.size1();
  int sign=G.sig();
  double t;
  cplx res,err,cplx_i=cplx(0,1);
  fourier::adft_func adft;
  dos_wrapper<dos_function> dos1(dos,sign,mu);
  dos1.mu_=mu;
  dos1.x_=ret;
  adft.sample(0.0,dos.lo_,dos.hi_,dos1,nn,limit);
  for(l=0;l<=nt;l++){ // l = t-t'
   t=h*l;
   adft.dft(-t,res,err);
   res *= std::complex<double>(0,-1.0);
   for(i=l;i<=nt;i++) element_set<T,LARGESIZE>(size1,G.retptr(i,i-l),(std::complex<T>)(res));
  }
}

/// @private
template <typename T,class dos_function>
void green_equilibrium_mat(herm_matrix<T> &G,dos_function &dos,double beta,int limit,int nn,double mu)
{
  typedef std::complex<double> cplx;
  int ntau=G.ntau(),m,size1=G.size1();
  int sign=G.sig();
  double dtau;
  cplx res,err;
  fourier::adft_func adft;
  dos_wrapper<dos_function> dos1(dos,sign,mu);
  dos1.beta_=beta;
  dtau=beta/ntau;
  dos1.mu_=mu;
  dos1.x_=mat;
  for(m=0;m<=ntau;m++){
    dos1.tau_=m*dtau;
    adft.sample(0.0,dos.lo_,dos.hi_,dos1,nn,limit);
    adft.dft(0.0,res,err);
    element_set<T,LARGESIZE>(size1,G.matptr(m),(std::complex<T>)(res));
  }
}
/// @private
template <typename T,class dos_function>
void green_equilibrium_tv(herm_matrix<T> &G,dos_function &dos,double beta,double h,int limit,int nn,double mu)
{
  typedef std::complex<double> cplx;
  int ntau=G.ntau(),nt=G.nt(),n,m,size1=G.size1();
  double dtau;
  int sign=G.sig();
  cplx res,err,cplx_i=cplx(0,1);
  fourier::adft_func adft;
  dos_wrapper<dos_function> dos1(dos,sign,mu);
  dos1.beta_=beta;
  dtau=beta/ntau;
  dos1.x_=tv;
  dos1.mu_=mu;
  for(m=0;m<=ntau;m++){
    dos1.tau_=m*dtau;
    adft.sample(0.0,dos.lo_,dos.hi_,dos1,nn,limit);
    for(n=0;n<=nt;n++){
      adft.dft(-n*h,res,err);
      res *= std::complex<double>(0,-1.0);
      element_set<T,LARGESIZE>(size1,G.tvptr(n,m),(std::complex<T>)(res));
    }
  }
}
/// @private
template <typename T,class dos_function>
void green_equilibrium_les(herm_matrix<T> &G,dos_function &dos,double beta,double h,int limit,int nn,double mu)
{
  typedef std::complex<double> cplx;
  int nt=G.nt(),i,l,size1=G.size1();
  double t;
  int sign=G.sig();
  cplx res,err,cplx_i=cplx(0,1);
  fourier::adft_func adft;
  dos_wrapper<dos_function> dos1(dos,sign,mu);
  dos1.mu_=mu;
  dos1.x_=les;
  dos1.beta_=beta;
  adft.sample(0.0,dos.lo_-0.0,dos.hi_-0.0,dos1,nn,limit);
  for(l=0;l<=nt;l++){ // l = t'-t
   t=h*l;
   adft.dft(t,res,err);
   res *= std::complex<double>(0,-1.0);
   for(i=l;i<=nt;i++) element_set<T,LARGESIZE>(size1,G.lesptr(i-l,i),(std::complex<T>)res);
  }
}

/// @private
template <typename T,class dos_function>
void green_equilibrium(herm_matrix<T> &G,dos_function &dos,double beta,double h,int limit,int nn,double mu)
{
  green_equilibrium_mat(G,dos,beta,limit,nn,mu);
  green_equilibrium_ret(G,dos,h,limit,nn,mu);
  green_equilibrium_tv(G,dos,beta,h,limit,nn,mu);
  green_equilibrium_les(G,dos,beta,h,limit,nn,mu);
}

/** \brief <b> Equilibrium  propagator for the given density of states  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calculate the equilibrium propagator G for the given density of states
* > via \f$G(t,t') = -i \int d\omega A(\omega+\mu) exp(i\omega (t'-t)) [ \Theta(t,t') - \Theta(t',t)exp(-\beta*\omega)]\f$
* <!-- ARGUMENTS
*      ========= -->
*
* @param G
* > The output Greens function set to the equilibrium free propagator
* @param dos
* > density of states
* @param beta
* > inverse temperature
* @param h
* > timestep
* @param mu
* > chemical potential
* @param limit
* > max number of intervals in Fourier transform (default: 100)
* @param nn
* > number of points in each interval of the Fourier transform (default: 20)
*/
template <typename T,class dos_function>
void green_equilibrium(herm_matrix<T> &G,dos_function &dos,double beta,double h,double mu,int limit,int nn)
{
  green_equilibrium_mat(G,dos,beta,limit,nn,mu);
  green_equilibrium_ret(G,dos,h,limit,nn,mu);
  green_equilibrium_tv(G,dos,beta,h,limit,nn,mu);
  green_equilibrium_les(G,dos,beta,h,limit,nn,mu);
}


/** \brief <b> Equilibrium propagator for Bethe semicircular density of states. Matsubara only.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calculate the Matsubara Green's function \f$G\f$ for the Bethe density of states.
* <!-- ARGUMENTS
*      ========= -->
*
* @param G
* > The output Greens function set to the equilibrium free propagator
* @param beta
* > inverse temperature
* @param limit
* > max number of intervals in Fourier transform
* @param nn
* > number of points in each interval of the Fourier transform
* @param mu
* > chemical potential
*/
template <typename T>
void green_equilibrium_mat_bethe(herm_matrix<T> &G, double beta, int limit,
                                 int nn,double mu,double Eshift) {
    bethedos dos(Eshift);
    green_equilibrium_mat(G, dos, beta, limit, nn,mu);
}


/** \brief <b> Equilibrium  propagator for Bethe semicircular density of states  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calculate the equilibrium propagator G for the Bethe density of states
* <!-- ARGUMENTS
*      ========= -->
*
* @param G
* > The output Greens function set to the equilibrium free propagator
* @param mu
* > chemical potential
* @param beta
* > inverse temperature
* @param h
* > timestep
* @param limit
* > max number of intervals in Fourier transform
* @param nn
* > number of points in each interval of the Fourier transform
*/

template <typename T>
void green_equilibrium_bethe(herm_matrix<T> &G, double beta, double h,
                             int limit, int nn, double mu, double Eshift) {
    bethedos dos(Eshift);
    green_equilibrium(G, dos, beta, h, limit, nn,mu);
}





/** \brief <b> Interpolation for the second order propagator </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Interpolates H(t-dt/2) from H(t+dt), H(t), H(t-dt)
*   where H(t+dt) is obtained from the extrapolation if not yet known
*   H is original function
*   If fixham=true we assume that the hamiltonian is known for all times and there is no extrapolation
*
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > timestep at which propagator is calculated
* @param H
* > Time dependent Hamiltonian
* @param Hinter
* > Interpolated  Hamiltonian
* @param kt
* > Order of integrator used for extrapolation and interpolation
* @param fixHam
* > Hamiltonian is known for all times and no extrapolation is needed for the predictor/corrector
*/
template<typename T>
void interpolate_CF2(int tstp,cntr::function<T> &H,cdmatrix &Hinter,int kt,bool fixHam=false){
    int size=H.size1_;
    cdmatrix tmp1(size,size),tmp2(size,size);
    int ktextrap=(tstp<=kt ? tstp : kt);
    int n1=(tstp<=kt ? kt : tstp); //Once you are at kt the next point needs to be included
    // Extrapolate to get H(t + dt)
    if(!fixHam && tstp>kt){
        extrapolate_timestep(tstp-1,H,integration::I<double>(ktextrap));
    }
    // Interpolate to get H(t + dt/2)
    H.get_value(tstp-1,tmp1);
    H.get_value(tstp,tmp2);

    Hinter=interpolation(n1,(tstp-0.5),H,integration::I<double>(kt)); //ktextrap+1 since we added one more point
}


/** \brief <b> Interpolation for the forth order propagator </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Interpolates  \f$H(t - dt + dt*c1)\f$ and \f$H(t - dt + dt*c2)\f$, where
*   \f$c1 = 1/2 - \sqrt(3)/6\f$
*   \f$c2 = 1/2 + \sqrt(3)/6\f$
*   from \f$H(t+dt)\f$, \f$H(t)\f$, \f$H(t-dt)\f$
*   where \f$H(t+dt)\f$ is obtained from the extrapolation if not yet known.
*   If fixham=true we assume that the hamiltonian is known for all times and there is no extrapolation
*
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > timestep at which propagator is calculated
* @param H
* > Time dependent Hamiltonian
* @param Hinte1
* > Interpolated  Hamiltonian at c1
* @param Hinte2
* > Interpolated  Hamiltonian at c2
* @param kt
* > Order of integrator used for extrapolation and interpolation
* @param fixHam
* > Hamiltonian is known for all times and no extrapolation is needed for the predictor/corrector
*/
template<typename T>
void interpolate_CF4(int tstp,cntr::function<T> &H,cdmatrix &Hinte1,cdmatrix &Hinte2,int kt,bool fixHam=false){
    int size=H.size1_;
    cdmatrix tmp(size,size);
    int ktextrap=(tstp<=kt ? tstp : kt);
    int n1=(tstp<=kt ? kt : tstp);
    // Extrapolate to get H(t)
    if(!fixHam && tstp>kt){
        extrapolate_timestep(tstp-1,H,integration::I<double>(ktextrap));
    }
    Hinte1=interpolation(n1,(tstp-0.5-sqrt(3)/6.0),H,integration::I<double>(kt)); //ktextrap+1 since we added one more point
    Hinte2=interpolation(n1,(tstp-0.5+sqrt(3)/6.0),H,integration::I<double>(kt)); //ktextrap+1 since we added one more point

}

/** \brief <b> Propagator for time-dependent free Hamiltonian </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calculate the free propagator at time tstp from time dependent free Hamiltonian using high-order commutator-free exponential time-propagation,
*   see https://doi.org/10.1016/j.jcp.2011.04.006 for the description.
*   Currently implemented versions are the second order using one exponential CF2:1 (order=2) and fourth order using two exponentials CF4:2 (order=4),
*   see also article for more details.
*
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > timestep at which propagator is calculated
* @param U
* > propagator as an output of the function
* @param H
* > Time dependent Hamiltonian
* @param dt
* > time step
* @param order
* > Order of approximation for commutator-free exponential, currently implemented orders = 2,4
* @param kt
* > Order of integrator used for extrapolation and interpolation
* @param fixHam
* > Hamiltonian is known for all times and no extrapolation is needed for the predictor/corrector
*/
template<typename T>
void propagator_exp(int tstp,cntr::function<T> &U,cntr::function<T> &H,double dt,int order,int kt,bool fixHam=false){
  int nt=U.nt_;
  int size=U.size1_;
  cdmatrix prop(size,size);

  assert(tstp<=nt);
  assert(order==2 || order==4);
  assert(size==H.size1_);

  if(tstp==-1 || tstp==0){
    prop.setIdentity();
    U.set_value(tstp,prop);
  }else{
    if(order==2){
        cdmatrix arg(size,size);
        // Get H(t+dt/2)-> Extrapolate and interpolate
        interpolate_CF2(tstp,H,arg,kt,fixHam);
        arg=arg*std::complex<double>(0.0,-1.0)*dt;
        U.get_value(tstp-1,prop);
        prop=arg.exp()*prop;
        U.set_value(tstp,prop);
    }else if(order==4){
       cdmatrix H1(size,size),H2(size,size);
       cdmatrix arg1(size,size),arg2(size,size);
       // Get H(t+dt*c1) and H(t+dt*c2) -> Extrapolate and interpolate
       interpolate_CF4(tstp,H,H1,H2,kt,fixHam);
       double a1=(3.0-2.0*sqrt(3.0))/12.0;
       double a2=(3.0+2.0*sqrt(3.0))/12.0;
       arg1=std::complex<double>(0.0,-1.0)*dt*(a1*H1+a2*H2);
       arg2=std::complex<double>(0.0,-1.0)*dt*(a2*H1+a1*H2);
       U.get_value(tstp-1,prop);
       prop=arg1.exp()*arg2.exp()*prop;
       U.set_value(tstp,prop);
    }
  }
}

/// @private
template<typename T,int SIZE>
void green_from_H_dispatch(herm_matrix<T> &G,T mu,cntr::function<T> &eps,T beta,T h,int kt,int order,bool fixHam=false){
  std::complex<T> iu = std::complex<T>(0.0, 1.0);
  int nt=G.nt(),ntau=G.ntau();
  int size=G.size1();
  int sign=G.sig();
  double tau,t,dtau=beta/ntau;
  cdmatrix H0(size,size),H1(size,size);
  eps.get_value(-1,H0);
  H1=mu*cdmatrix::Identity(size,size)-H0;
  cntr::function<T> Ut(nt,size);
  cdmatrix evec0(size,size),value(size,size);
  dvector eval0(size),eval0m(size);
  Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(H1);
  evec0=eigensolver.eigenvectors();
  eval0=eigensolver.eigenvalues();
  eval0m=(-1.0)*eval0;


  for(int m=0;m<=ntau;m++){
    tau=m*dtau;
    if(sign==-1){
      value=(-1.0)*evec0*fermi_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
    }else if(sign==1){
      value=evec0*bose_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
    }
    G.set_mat(m,value);
  }

  if(nt >= 0){
    cdmatrix idm(size,size);
    idm = MatrixXcd::Identity(size,size);
    Ut.set_value(-1,idm);
    Ut.set_value(0,idm);
    for(int tstp=1;tstp<=nt;tstp++){
      propagator_exp(tstp,Ut,eps,h,order,kt,fixHam);
    }
    for(int tstp=0;tstp<=nt;tstp++){
      cdmatrix tmp;
      Ut.get_value(tstp,tmp);
      tmp = tmp * std::complex<T>(cos(mu * h * tstp),sin(mu * h * tstp));
      Ut.set_value(tstp,tmp);
    }
    //!! Propagator here checked !!

    cdmatrix expp(size,size);
    for(int m=0;m<=ntau;m++){
      tau=m*dtau;
      for(int n=0;n<=nt;n++){
        Ut.get_value(n,expp);
        if(sign==-1){
          value=iu*expp*evec0*fermi_exp(beta,tau,eval0m).asDiagonal()*evec0.adjoint();
        }else if(sign==1){
          value=-iu*expp*evec0*bose_exp(beta,tau,eval0m).asDiagonal()*evec0.adjoint();
        }
        G.set_tv(n,m,value);
      }
    }
    if(sign==-1){
      value=evec0*fermi(beta,eval0m).asDiagonal()*evec0.adjoint();
    }else if(sign==1){
      value=-1.0*evec0*bose(beta,eval0m).asDiagonal()*evec0.adjoint();
    }
    cdmatrix exppt1(size,size);
    cdmatrix exppt2(size,size);
    for(int m=0;m<=nt;m++){
      for(int n=0;n<=m;n++){
	      cdmatrix tmp(size,size);
	      Ut.get_value(m,exppt1);
	      Ut.get_value(n,exppt2);
	      tmp = -iu*exppt1*exppt2.adjoint();
	      G.set_ret(m,n,tmp);
	      tmp=iu*exppt2*value*exppt1.adjoint();
	      G.set_les(n,m,tmp);
      }
    }
  }
}

/// @private
template<typename T,int SIZE>
void green_from_H_const_dispatch(herm_matrix<T> &G,T mu,cdmatrix &H0,T beta,T h){
  std::complex<T> iu = std::complex<T>(0.0, 1.0);
  int nt=G.nt(),ntau=G.ntau();
  int size=G.size1();
  int sign=G.sig();
  double tau,t,dtau=beta/ntau;
  cdmatrix idm(size,size);
  cdmatrix Udt(size,size);
  cdmatrix IHdt(size,size);
  cdmatrix Hmu(size,size);
  cdmatrix evec0(size,size),value(size,size);
  dvector eval0(size),eval0m(size);

  idm = MatrixXcd::Identity(size,size);
  Hmu = -H0 + mu * idm;

  Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(Hmu);
  evec0=eigensolver.eigenvectors();
  eval0=eigensolver.eigenvalues();
  eval0m=(-1.0)*eval0;

  for(int m=0;m<=ntau;m++){
    tau=m*dtau;
    if(sign==-1){
      value=(-1.0)*evec0*fermi_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
    }else if(sign==1){
      value=(1.0)*evec0*bose_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
    }
    G.set_mat(m,value);
  }

  if(nt >=0 ){
    // IHdt = -iu * h * H0;
    IHdt = iu * h * Hmu;
    Udt = IHdt.exp();

    cntr::function<T> Ut(nt,size);
    cdmatrix Un(size,size);
    Ut.set_value(-1,idm);
    Ut.set_value(0,idm);
    for(int n=1;n<=nt;n++){
      Ut.get_value(n-1,Un);
      Un = Un * Udt;
      Ut.set_value(n,Un);
    }

    cdmatrix expp(size,size);
    for(int m=0;m<=ntau;m++){
      tau=m*dtau;
      for(int n=0;n<=nt;n++){
        Ut.get_value(n,expp);
        if(sign==-1){
          value=iu*expp*evec0*fermi_exp(beta,tau,eval0m).asDiagonal()*evec0.adjoint();
        }else if(sign==1){
          value=-iu*expp*evec0*bose_exp(beta,tau,eval0m).asDiagonal()*evec0.adjoint();
        }
        G.set_tv(n,m,value);
      }
    }

    if(sign==-1){
      value=evec0*fermi(beta,eval0m).asDiagonal()*evec0.adjoint();
    }else if(sign==1){
      value=-1.0*evec0*bose(beta,eval0m).asDiagonal()*evec0.adjoint();
    }
    cdmatrix exppt1(size,size);
    cdmatrix exppt2(size,size);
    for(int m=0;m<=nt;m++){
      for(int n=0;n<=m;n++){
	      cdmatrix tmp(size,size);
	      Ut.get_value(m,exppt1);
	      Ut.get_value(n,exppt2);
        tmp = -iu*exppt1*exppt2.adjoint();
	      G.set_ret(m,n,tmp);
        tmp = iu*exppt2*value*exppt1.adjoint();
	      G.set_les(n,m,tmp);
      }
    }
  }
}

/// @private
template<typename T,int SIZE>
void green_from_H_dispatch(herm_matrix_timestep<T> &G,T mu,cntr::function<T> &eps,T beta,T h,int kt,int order,bool fixHam=false){
  std::complex<T> iu = std::complex<T>(0.0, 1.0);
  int tstp=G.tstp(),ntau=G.ntau();
  int size=G.size1();
  int sign=G.sig();
  cdmatrix H0(size,size),H1(size,size);
  eps.get_value(-1,H0);
  H1=mu*cdmatrix::Identity(size,size)-H0;
  cntr::function<T> Ut(tstp,size);
  cdmatrix evec0(size,size),value(size,size);
  dvector eval0(size),eval0m(size);
  Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(H1);
  evec0=eigensolver.eigenvectors();
  eval0=eigensolver.eigenvalues();
  eval0m=(-1.0)*eval0;

  if(tstp==-1){
    for(int m=0;m<=ntau;m++){
        double tau,t,dtau=beta/ntau;
        tau=m*dtau;
        if(sign==-1){
          value=(-1.0)*evec0*fermi_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
        }else if(sign==1){
          value=evec0*bose_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
        }
        G.set_mat(m,value);
    }
  }else{
    cdmatrix expp(size,size);
    cdmatrix exppt1(size,size);
    cdmatrix exppt2(size,size);
    cdmatrix idm(size,size);
    idm = MatrixXcd::Identity(size,size);
    Ut.set_value(-1,idm);
    Ut.set_value(0,idm);
    for(int t=1;t<=tstp;t++){
        propagator_exp(t,Ut,eps,h,order,kt,fixHam);
    }

    for(int n=0; n<=tstp;n++){
      cdmatrix tmp;
      Ut.get_value(tstp,tmp);
      tmp = tmp * std::complex<T>(cos(mu * h * n),sin(mu * h * n));
      Ut.set_value(tstp,tmp);
    }


    // TV
    for(int m=0;m<=ntau;m++){
        double tau,dtau=beta/ntau;
        tau=m*dtau;
        Ut.get_value(tstp,expp);
        if(sign==-1){
          value=iu*expp*evec0*fermi_exp(beta,tau,eval0m).asDiagonal()*evec0.adjoint();
        }else if(sign==1){
          value=-iu*expp*evec0*bose_exp(beta,tau,eval0m).asDiagonal()*evec0.adjoint();
        }
        G.set_tv(m,value);
    }
    // RET
    for(int n=0;n<=tstp;n++){
       cdmatrix tmp(size,size);
       Ut.get_value(tstp,exppt1);
       Ut.get_value(n,exppt2);
       tmp = -iu*exppt1*exppt2.adjoint();
       G.set_ret(n,tmp);
    }
    // LES
    if(sign==-1){
      value=evec0*fermi(beta,eval0m).asDiagonal()*evec0.adjoint();
    }else if(sign==1){
      value=-1.0*evec0*bose(beta,eval0m).asDiagonal()*evec0.adjoint();
    }
    for(int n=0;n<=tstp;n++){
       cdmatrix tmp(size,size);
       Ut.get_value(tstp,exppt1);
       Ut.get_value(n,exppt2);
       tmp = iu*exppt2*value*exppt1.adjoint();
       G.set_les(n,tmp);
    }
  }
}

/// @private
template<typename T,int SIZE>
void green_from_H_const_dispatch(herm_matrix_timestep<T> &G,T mu,cdmatrix &H0,T beta,T h){
  std::complex<T> iu = std::complex<T>(0.0, 1.0);
  int tstp=G.tstp(),ntau=G.ntau();
  int size=G.size1();
  int sign=G.sig();
  cdmatrix evec0(size,size),value(size,size);
  dvector eval0(size),eval0m(size);
  cdmatrix H1;
  H1 = mu*cdmatrix::Identity(size,size)-H0;

  Eigen::SelfAdjointEigenSolver<cdmatrix> eigensolver(H1);
  evec0=eigensolver.eigenvectors();
  eval0=eigensolver.eigenvalues();
  eval0m=(-1.0)*eval0;

  if(tstp==-1){
    for(int m=0;m<=ntau;m++){
        double tau,t,dtau=beta/ntau;
        tau=m*dtau;
        if(sign==-1){
          value=(-1.0)*evec0*fermi_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
        }else if(sign==1){
          value=evec0*bose_exp(beta,tau,eval0).asDiagonal()*evec0.adjoint();
        }
        G.set_mat(m,value);
    }
  }else{
    cdmatrix Udt(size,size);
    cdmatrix IHdt(size,size);
    cdmatrix expp(size,size);
    cdmatrix exppt1(size,size);
    cdmatrix exppt2(size,size);

    // IHdt = -iu * h * H0;
    IHdt = iu * h * H1;
    Udt = IHdt.exp();

    cntr::function<T> Ut(tstp,size);
    cdmatrix idm(size,size);
    idm = MatrixXcd::Identity(size,size);
    Ut.set_value(-1,idm);
    Ut.set_value(0,idm);
    for(int n=1;n<=tstp;n++){
      cdmatrix Un(size,size);
      Ut.get_value(n-1,Un);
      Un = Un * Udt;
      Ut.set_value(n,Un);
    }
    // TV
    for(int m=0;m<=ntau;m++){
        double tau,dtau=beta/ntau;
        tau=m*dtau;
        Ut.get_value(tstp,expp);
        if(sign==-1){
          value=iu*expp*evec0*fermi_exp(beta,tau,eval0m).asDiagonal()*evec0.adjoint();
        }else if(sign==1){
          value=-iu*expp*evec0*bose_exp(beta,tau,eval0m).asDiagonal()*evec0.adjoint();
        }
        G.set_tv(m,value);
    }
    // RET
    for(int n=0;n<=tstp;n++){
       cdmatrix tmp(size,size);
       Ut.get_value(tstp,exppt1);
       Ut.get_value(n,exppt2);
       tmp = -iu * (exppt1 * exppt2.adjoint());
       G.set_ret(n,tmp);
    }
    // LES
    if(sign==-1){
      value=evec0*fermi(beta,eval0m).asDiagonal()*evec0.adjoint();
    }else if(sign==1){
      value=-1.0*evec0*bose(beta,eval0m).asDiagonal()*evec0.adjoint();
    }
    for(int n=0;n<=tstp;n++){
       cdmatrix tmp(size,size);
       Ut.get_value(tstp,exppt1);
       Ut.get_value(n,exppt2);
       tmp = iu * (exppt2 * value * exppt1.adjoint());
       G.set_les(n,tmp);
    }
  }
}

/** \brief <b> Propagator for time-independent free Hamiltonian  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calculate the free propagator G from fixed quadratic Hamiltonian using high-order commutator-free exponential time-propagation,
* >  see https://doi.org/10.1016/j.jcp.2011.04.006 for the description.
* >  Currently implemented versions are the second order using one exponential CF2:1 (order=2) and fourth order using two exponentials CF4:2 (order=4),
* >  see also article for more details.
* <!-- ARGUMENTS
*      ========= -->
*
* @param G
* > The output Greens function set to time dependent free propagator
* @param mu
* > chemical potential
* @param eps
* > time-independent representation of quadratical hamiltonian
* @param beta
* > inverse temperature
* @param h
* > timestep
*/
template<typename T>
void green_from_H(herm_matrix<T> &G,T mu,cdmatrix &eps,T beta,T h){
  assert(G.size1()==eps.rows());
  assert(eps.rows()==eps.cols());

  int size=G.size1();
  if(size==1) green_from_H_const_dispatch<T,1>(G,mu,eps,beta,h);
  else green_from_H_const_dispatch<T,LARGESIZE>(G,mu,eps,beta,h);
}

/// @private
template<typename T>
void green_from_H(herm_matrix<T> &G,T mu,cntr::function<T> &eps,T beta,T h,bool fixHam,int SolveOrder,int cf_order){
  assert(G.size1()==eps.size2_);
  assert(eps.size1_==eps.size2_);
  assert(SolveOrder <= MAX_SOLVE_ORDER);
  assert(cf_order == 2 || cf_order == 4);

  int size=G.size1();
  if(size==1) green_from_H_dispatch<T,1>(G,mu,eps,beta,h,SolveOrder,cf_order,fixHam);
  else green_from_H_dispatch<T,LARGESIZE>(G,mu,eps,beta,h,SolveOrder,cf_order,fixHam);
}


/** \brief <b> Propagator for time-dependent free Hamiltonian  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calculate the free propagator G from time dependent quadratic Hamiltonian using high-order commutator-free exponential time-propagation,
* >  see https://doi.org/10.1016/j.jcp.2011.04.006 for the description.
* >  Currently implemented versions are the second order using one exponential CF2:1 (order=2) and fourth order using two exponentials CF4:2 (order=4),
* >  see also article for more details.
* <!-- ARGUMENTS
*      ========= -->
*
* @param G
* > The output Greens function set to time dependent free propagator
* @param mu
* > chemical potential
* @param eps
* > time dependent representation of quadratical hamiltonian
* @param beta
* > inverse temperature
* @param h
* > timestep
* @param SolveOrder
* > Order of integrator used for extrapolation and interpolation
* @param cf_order
* > Order of approximation for commutator-free exponential, currently implemented orders = 2,4
*/

template<typename T>
void green_from_H(herm_matrix<T> &G,T mu,cntr::function<T> &eps,T beta,T h,int SolveOrder,int cf_order){
  assert(G.size1()==eps.size2_);
  assert(eps.size1_==eps.size2_);
  assert(SolveOrder <= MAX_SOLVE_ORDER);
  assert(cf_order == 2 || cf_order == 4);

  int size=G.size1();
  if(size==1) green_from_H_dispatch<T,1>(G,mu,eps,beta,h,SolveOrder,cf_order,true);
  else green_from_H_dispatch<T,LARGESIZE>(G,mu,eps,beta,h,SolveOrder,cf_order,true);
}



/// @private
template<typename T>
void green_from_H(herm_matrix_timestep<T> &G,T mu,cdmatrix &eps,T beta,T h){
  assert(G.size1()==eps.rows());
  assert(eps.rows()==eps.cols());

  int size=G.size1();
  if(size==1) green_from_H_const_dispatch<T,1>(G,mu,eps,beta,h);
  else green_from_H_const_dispatch<T,LARGESIZE>(G,mu,eps,beta,h);
}


/** \brief <b> Propagator for time-independent free Hamiltonian  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calculate the free propagator G from fixed quadratic Hamiltonian using high-order commutator-free exponential time-propagation,
* >  see https://doi.org/10.1016/j.jcp.2011.04.006 for the description.
* >  Currently implemented versions are the second order using one exponential CF2:1 (order=2) and fourth order using two exponentials CF4:2 (order=4),
* >  see also article for more details.
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > given time step
* @param G
* > The output: a timestep of the Greens function set to time dependent free propagator
* @param mu
* > chemical potential
* @param eps
* > time-independent representation of quadratical hamiltonian
* @param beta
* > inverse temperature
* @param h
* > timestep
*/
template<typename T>
void green_from_H(int tstp, herm_matrix_timestep<T> &G,T mu,cdmatrix &eps,T beta,T h){
  assert(tstp == G.tstp());
  assert(G.size1()==eps.rows());
  assert(eps.rows()==eps.cols());

  int size=G.size1();
  if(size==1) green_from_H_const_dispatch<T,1>(G,mu,eps,beta,h);
  else green_from_H_const_dispatch<T,LARGESIZE>(G,mu,eps,beta,h);
}


/** \brief <b> Propagator for time-independent free Hamiltonian  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calculate the free propagator G from fixed quadratic Hamiltonian using high-order commutator-free exponential time-propagation,
* >  see https://doi.org/10.1016/j.jcp.2011.04.006 for the description.
* >  Currently implemented versions are the second order using one exponential CF2:1 (order=2) and fourth order using two exponentials CF4:2 (order=4),
* >  see also article for more details.
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > the index of the time step
* @param G
* > The output: a time step of the Greens function set to time dependent free propagator
* @param mu
* > chemical potential
* @param eps
* > time-independent representation of quadratical hamiltonian
* @param beta
* > inverse temperature
* @param h
* > timestep interval
*/
template<typename T>
void green_from_H(int tstp, herm_matrix<T> &G,T mu,cdmatrix &eps,T beta,T h){
  assert(G.size1()==eps.rows());
  assert(eps.rows()==eps.cols());
  herm_matrix_timestep<T> Gstep(tstp,G.ntau(),G.size1(),G.sig());

  int size=G.size1();
  if(size==1) green_from_H_const_dispatch<T,1>(Gstep,mu,eps,beta,h);
  else green_from_H_const_dispatch<T,LARGESIZE>(Gstep,mu,eps,beta,h);

  G.set_timestep(tstp, Gstep);
}


/// @private
template<typename T>
void green_from_H(herm_matrix_timestep<T> &G,T mu,cntr::function<T> &eps,T beta,T h,bool fixHam,int SolveOrder,int cf_order){
  assert(G.size1()==eps.size2_);
  assert(eps.size1_==eps.size2_);
  assert(SolveOrder <= MAX_SOLVE_ORDER);
  assert(cf_order == 2 || cf_order == 4);

  int size=G.size1();
  if(size==1) green_from_H_dispatch<T,1>(G,mu,eps,beta,h,SolveOrder,cf_order,fixHam);
  else green_from_H_dispatch<T,LARGESIZE>(G,mu,eps,beta,h,SolveOrder,cf_order,fixHam);
}

/** \brief <b> Propagator for time-dependent free Hamiltonian  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calculate the free propagator G from time dependent quadratic Hamiltonian using high-order commutator-free exponential time-propagation,
* >  see https://doi.org/10.1016/j.jcp.2011.04.006 for the description.
* >  Currently implemented versions are the second order using one exponential CF2:1 (order=2) and fourth order using two exponentials CF4:2 (order=4),
* >  see also article for more details.
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > the index of the time step
* @param G
* > The output: a timestep of the Greens function set to time dependent free propagator
* @param mu
* > chemical potential
* @param eps
* > time dependent representation of quadratical hamiltonian
* @param beta
* > inverse temperature
* @param h
* > time step interval
* @param fixHam
* > If True Hamiltonian is known for all times and no extrapolation is needed for the predictor/corrector
* @param SolveOrder
* > Order of integrator used for extrapolation and interpolation
* @param cf_order
* > Order of approximation for commutator-free exponential, currently implemented orders = 2,4
*/
template<typename T>
void green_from_H(int tstp, herm_matrix_timestep<T> &G,T mu,cntr::function<T> &eps,T beta,T h,bool fixHam,int SolveOrder,int cf_order){
  assert(tstp == G.tstp());
  assert(G.size1()==eps.size2_);
  assert(eps.size1_==eps.size2_);
  assert(SolveOrder <= MAX_SOLVE_ORDER);
  assert(cf_order == 2 || cf_order == 4);

  int size=G.size1();
  if(size==1) green_from_H_dispatch<T,1>(G,mu,eps,beta,h,SolveOrder,cf_order,fixHam);
  else green_from_H_dispatch<T,LARGESIZE>(G,mu,eps,beta,h,SolveOrder,cf_order,fixHam);
}

/** \brief <b> Propagator for time-dependent free Hamiltonian  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calculate the free propagator G from time dependent quadratic Hamiltonian using high-order commutator-free exponential time-propagation,
* >  see https://doi.org/10.1016/j.jcp.2011.04.006 for the description.
* >  Currently implemented versions are the second order using one exponential CF2:1 (order=2) and fourth order using two exponentials CF4:2 (order=4),
* >  see also article for more details.
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > the index of the time step
* @param G
* > The output: a timestep of the Greens function set to time dependent free propagator
* @param mu
* > chemical potential
* @param eps
* > time dependent representation of quadratical hamiltonian
* @param beta
* > inverse temperature
* @param h
* > time step interval
* @param fixHam
* > If True Hamiltonian is known for all times and no extrapolation is needed for the predictor/corrector
* @param SolveOrder
* > Order of integrator used for extrapolation and interpolation
* @param cf_order
* > Order of approximation for commutator-free exponential, currently implemented orders = 2,4
*/
template<typename T>
void green_from_H(int tstp, herm_matrix<T> &G,T mu,cntr::function<T> &eps,T beta,T h,bool fixHam,int SolveOrder,int cf_order){
  assert(tstp <= G.nt());
  assert(G.size1()==eps.size2_);
  assert(eps.size1_==eps.size2_);
  assert(SolveOrder <= MAX_SOLVE_ORDER);
  assert(cf_order == 2 || cf_order == 4);
  herm_matrix_timestep<T> Gstep(tstp,G.ntau(),G.size1(),G.sig());

  int size=G.size1();
  if(size==1) green_from_H_dispatch<T,1>(Gstep,mu,eps,beta,h,SolveOrder,cf_order,fixHam);
  else green_from_H_dispatch<T,LARGESIZE>(Gstep,mu,eps,beta,h,SolveOrder,cf_order,fixHam);

  G.set_timestep(tstp, Gstep);
}


///////////////////////////////////////////////////////////////////////////////////////
// BOSONIC GREENS FUNCTION:
// G(t,t') = - ii * <TC b(t) bdag(t') >
// with H = w * b^dag b
/*template<typename T>
void green_single_pole_bose(herm_matrix<T> &G,T *eps,T beta,T h){
  int nt=G.nt(),size1=G.size1();
  int ntau=G.ntau();
  int n,m,i0=0;
  double tau,fm,eps1;
  std::complex<T> x,ii=std::complex<T>(0,1.0);
  std::vector<std::complex<T> > expp;
  // expp[n] = -exp(ii*w*n)
  expp.resize(nt+1);
  G.clear();
  double dtau=beta/ntau;
  for(int i=0;i<size1;i++){
    i0=i*size1+i;
    eps1=-eps[i];
    for(m=0;m<=nt;m++){
      expp[m]=std::complex<T>(cos(eps1*m*h),sin(eps1*m*h));
    }
    fm=bose<T>(beta,-eps1);
    for(int m=0;m<=ntau;m++){
      tau=m*dtau;
      G.matptr(m)[i0]=bose_exp<T>(beta,tau,eps1);
      x=-ii*bose_exp<T>(beta,tau,-eps1);
      for(n=0;n<=nt;n++) G.tvptr(n,m)[i0]=x*expp[n];
    }
    x=bose<T>(beta,-eps1);
    for(m=0;m<=nt;m++){
      for(n=0;n<=m;n++){
        G.retptr(m,n)[i0]=-ii*expp[m-n];
        G.lesptr(n,m)[i0]=-ii*x*conj(expp[m-n]);
      }
    }
  }
}
*/

/// @private
template<typename T>
void green_single_pole_XX_timestep(int tstp,int ntau,int size1,std::complex<T> *ret,std::complex<T> *tv,std::complex<T> *les,T w,T beta,T h) {

	double bw = bose(beta,w);
	double c,s,ewv,emwv;
	int s2=size1*size1;
	std::complex<T> d_les, d_gtr, d_ret, d_tv;
	std::complex<T> i = std::complex<T>(0.0, 1.0);
	double mat_axisunit = beta/ntau;
	assert(w*beta>0);

	for(int tp=0;tp<=tstp;tp++){
		c = cos(w*(tstp-tp)*h);
		s = sin(w*(tstp-tp)*h);
		d_les = -i *  ( (c-i*s) + 2.0 * c * bw );	// D_les(tp, t)
		d_gtr = -i *  ( (c-i*s) + 2.0 * c * bw );	// D_gtr(t, tp)
		d_ret = d_gtr + conj(d_les);				// D_ret(t, tp)
		les[tp*s2]=0.5*d_les;
		ret[tp*s2]=0.5*d_ret;
	}
	for(int v=0;v<=ntau;v++) {
		c = cos(w*tstp*h);
		s = sin(w*tstp*h);
		ewv  = exp(w*v*mat_axisunit);
		emwv = exp(-w*v*mat_axisunit);
		d_tv = -i * ( (c+i*s)*emwv + ((c+i*s)*emwv + (c-i*s)*ewv)*bw ); // ok
		tv[v*s2]=0.5*d_tv;
	}
}
/// @private
template<typename T>
void green_single_pole_XX_mat(int ntau,int size1,std::complex<T> *mat,T w,T beta) {
	double ewv,emwv;
	int s2=size1*size1;
	std::complex<T> d_mat;
	std::complex<T> i = std::complex<T>(0.0, 1.0);
	double mat_axisunit = beta/ntau;
	assert(w*beta>0);

	for(int v=0;v<=ntau;v++) {
		ewv=bose_exp(beta,v*mat_axisunit,w); // exp(tau*w)b(w)
		emwv=-bose_exp(beta,v*mat_axisunit,-w); // -exp(-tau*w)b(-w)
		d_mat = -0.5*(ewv+emwv);
		mat[v*s2]=d_mat;
	}
}
/// @private
template<typename T>
void green_single_pole_XX_timestep(herm_matrix_timestep<T> &D0,T w,T beta,T h) {
	int tstp=D0.tstp_;
	if(tstp==-1) green_single_pole_XX_mat(D0.ntau(),D0.size1(),D0.matptr(0),w,beta);
	else green_single_pole_XX_timestep(tstp,D0.ntau(),D0.size1(),D0.retptr(0),D0.tvptr(0),D0.lesptr(0),w,beta,h);
}
/// @private
template<typename T>
void green_single_pole_XX_timestep(int tstp, herm_matrix_timestep<T> &D0,T w,T beta,T h) {
  assert(tstp == D0.tstp());
  if(tstp==-1) green_single_pole_XX_mat(D0.ntau(),D0.size1(),D0.matptr(0),w,beta);
  else green_single_pole_XX_timestep(tstp,D0.ntau(),D0.size1(),D0.retptr(0),D0.tvptr(0),D0.lesptr(0),w,beta,h);
}
/// @private
template<typename T>
void green_single_pole_XX_timestep(int tstp,herm_matrix<T> &D0,T w,T beta,T h) {
	if(tstp==-1) green_single_pole_XX_mat(D0.ntau(),D0.size1(),D0.matptr(0),w,beta);
	else green_single_pole_XX_timestep(tstp,D0.ntau(),D0.size1(),D0.retptr(tstp,0),D0.tvptr(tstp,0),D0.lesptr(0,tstp),w,beta,h);
}


/** \brief <b> Propagator for the displacement operator \f$X=(b+b^{\dagger})/\sqrt(2)\f$ in the thermal equilibrium  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Calculate the free displacement \f$X=(b+b^{\dagger})/\sqrt(2)\f$ propagator GD0 in equilibrium
* <!-- ARGUMENTS
*      ========= -->
*
* @param D0
* > The output Greens function set to the equilibrium free propagator
* @param w
* > bosonic frequency
* @param beta
* > inverse temperature
* @param h
* > timestep
*/

template<typename T>
void green_single_pole_XX(herm_matrix<T> &D0,T w,T beta,T h) {
	int nt=D0.nt(),t;
	if(nt>=-1) green_single_pole_XX_mat(D0.ntau(),D0.size1(),D0.matptr(0),w,beta);
	for(t=0;t<=nt;t++){
	    green_single_pole_XX_timestep(t,D0.ntau(),D0.size1(),D0.retptr(t,0),D0.tvptr(t,0),D0.lesptr(0,t),w,beta,h);
	}
}


} // namespace cntr

#endif  // CNTR_EQUILIBRIUM_IMPL_H
