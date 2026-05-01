#ifndef NESS2_EQUILIBRIUM_DECL_H
#define NESS2_EQUILIBRIUM_DECL_H

#include "cntr/ness2_global_settings.hpp"
#include "cntr/ness2_herm_matrix_ness_decl.hpp"
#include "cntr/cntr.hpp"

namespace ness2 {

  //@private
  template <typename DosFunction>
  void green_equilibrium_ret_ness_adaptive(herm_matrix_ness &G, DosFunction &dos, double h, int nn,int limit);

  //@private
  template <typename Integrand>
  void compute_lesser(fft_array &G, int nt,Integrand &integrand, double lo, double hi, double h, int nn, int limit);
  

/**
 * \brief <b>Equilibrium steady-state Green's function systems</b>
 *
 * Constructs a retarded and lesser Green's function in the time domain for a system
 * in thermal equilibrium, using a specified density of states and numerical
 * Fourier integration method.
 *
 * \tparam DOS
 *   A density of states class providing:
 *   - `double lo_` — lower frequency bound
 *   - `double hi_` — upper frequency bound
 *   - `double operator()(double omega)` — density of states  \f$ A(\omega) \f$
 *
 * \param bose_fermi
 * = BOSON (+1) or FERMION (-1) sign
 * \param G
 *   Output: Green's function object (`herm_matrix_ness`) to be filled in time domain
 * \param dos
 *   Density of states object used to define the spectral function
 * \param beta
 *   Inverse temperature \f$ \beta \f$
 * \param h
 *   Time step size
 * \param mu
 *   Chemical potential \f$ \mu \f$
 * \param method
 *   Integration method to use for the Fourier transform. One of:
 *   - `FFT_TRAPEZ`: Simple piecewise-constant quadrature (trapezoidal rule)
 *   - `FFT_ADAPTIVE`: Adaptive quadrature as in libcntr (uses `nn` and `limit`)
 * \param nn
 *   Number of sub-sample points per interval (used only if `method == FFT_ADAPTIVE`)
 * \param limit
 *   Maximum number of intervals in adaptive Fourier transform (used only if `method == FFT_ADAPTIVE`)
 *
 * \details
 * Computes the time-domain Green's functions:
 * \f[
 * G^R(t) = -i \int_{\mathrm{lo}}^{\mathrm{hi}} d\omega\, A(\omega) e^{-i \omega t}
 * \f]
 * \f[
 * G^<(t) = - \xi i \int_{\mathrm{lo}}^{\mathrm{hi}} d\omega\, A(\omega) f_F(\omega - \mu) e^{-i \omega t},
 * \quad f_F(x) = \frac{1}{e^{\beta x} - \xi}
 * \f]
 * where $\xi$ = BOSE or FERMI, and \f$ A(\omega)\f$ is computed from the density of states.
 *
 * The integration bounds are taken from `dos.lo_` to `dos.hi_`.
 *
 * For `FFT_ADAPTIVE`, a frequency-adaptive scheme is used, refining the interval until convergence.
 * For `FFT_TRAPEZ` the FFT grid defined by Nft_ and h is used, nn and limit are ignored
 */
  template <typename DOS>
  void green_equilibrium_ness(int bosefermi,herm_matrix_ness &G, DOS &dos, double beta,double mu,
			      double h,fft_integral_method method, int limit=100,int nn=20);
  
  /** \brief <b> Class `bethe` represent the semicircular density of states </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * > This class contains the data structures for representing the density of states
   * > for a Bethe semicircular with the bandwidth 4V, centered around E0
   *
   */
  class bethedos {
  public:
    /** \brief <b> Higher edge of the density of states </b> */
    double hi_;
    /** \brief <b> Lower edge of the density of states </b> */
    double lo_;
    /** \brief <b>  4V corresponds to the bandwidth </b> */
    double V_;
    double E0_;
    bethedos();
    /** \brief <b> DOS between a and b, i.e. center = (a+b)/2, width=(b-a)</b> */
    bethedos(double a, double b);
    double operator()(double x);
  };
  
  /** \brief <b> Class `gauss` represent the gaussian density of states </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * > This class contains the data structures for representing the density of states
   * > for a gaussian  with mean 0 and a standard deviation of 1.
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
  /** \brief <b> Class `ohmic` represents a symmetric ohmic bath x^2*exp(-x/x_c) density of states </b>
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
  
  
  // implementation:
  //@private
  template <typename DOS>
  void green_equilibrium_ness_fermi(herm_matrix_ness &G, DOS &dos, double beta,double mu, 
				    double h,fft_integral_method method, int limit, int nn) {
    if(method==FFT_ADAPTIVE){
      // Retarded Green’s function
      green_equilibrium_ret_ness_adaptive(G,dos,h,nn,limit);      
      // Lesser Green’s function
      auto fermi_integrand = [&](double omega) {
	double arg = beta * (omega - mu);
	double fm;
	if (arg > 100) fm=0.0;
	else if (arg < -100) fm=1.0;
	else fm=1.0/(1.0 + std::exp(arg));      
	return -dos(omega) *fm;
      };
      // Compute lesser Green’s function
      compute_lesser(G.les_,G.Nft_/2,fermi_integrand,dos.lo_,dos.hi_,h,nn,limit);
    }else{
      // METHOD=TRAPEZ
      int Nft=G.Nft_;
      int size1=G.size1_;
      grid_info grid(Nft,h);
      cdmatrix Id = cdmatrix::Identity(size1, size1);
      int nfreq=Nft/2;
      G.ret_.set_zero(fft_domain::freq);
      for(int w = -nfreq+1; w < nfreq; ++w){
	double omega=grid.freq_at(w);
	double aw=dos(omega);
	double fm;
	cdmatrix tmp=Id*cdouble(0,-M_PI*aw);
	// set only Im(Gret), as only Re(Gret) is not used in backtransform
	G.set_ret(w,tmp,fft_domain::freq);
	double arg = beta * (omega - mu);
	if (arg > 100) fm=0.0;
	else if (arg < -100) fm=1.0;
	else fm=1.0/(1.0 + std::exp(arg));      
	G.set_les(w,-2*fm*tmp,fft_domain::freq);
      }
      // this transforms back the imaginary part of G(w) only
      G.integral_transform_to_time(h,FFT_TRAPEZ);
    }
    return;
  }

  //@private
  template <typename DOS>
  void green_equilibrium_ness_bose(herm_matrix_ness &G, DOS &dos, double beta, double mu,
				   double h,fft_integral_method method,  int limit,int nn) {
    auto bose_integrand = [&](double omega) {
      // A(w)*bose(w-mu), assuming that A(w)=0 at w=mu
      double arg = beta * (omega - mu);
      if (arg > 100) return 0.0;
      if (arg < -100) return -dos(omega);
      if (std::abs(arg) < 1e-8) {
	// Use small-argument limit: A(w) ~ (w-mu) assumed:
	double delta = 1e-6;
	double slope = (dos(delta + mu) - dos(mu)) / delta;
	return slope / beta;
      }
      return dos(omega) / (std::exp(arg) - 1.0);
    };          
    if(method==FFT_ADAPTIVE){
      green_equilibrium_ret_ness_adaptive(G,dos,h,nn,limit);
      compute_lesser(G.les_,G.Nft_/2, bose_integrand,dos.lo_,dos.hi_, h, nn, limit);
    }else{
      //method==FFT_TRAPEZ
      int Nft=G.Nft_;
      int size1=G.size1_;
      grid_info grid(Nft,h);
      cdmatrix Id = cdmatrix::Identity(size1, size1);
      int nfreq=Nft/2;
      G.ret_.set_zero(fft_domain::freq);
      for(int w = -nfreq+1; w < nfreq; ++w){
	double omega=grid.freq_at(w);
	double aw=dos(omega);
	cdouble fm=cdouble(0,-bose_integrand(omega)*2*M_PI);
	cdmatrix tmp=Id*cdouble(0,-M_PI*aw);
	G.set_ret(w,tmp,fft_domain::freq);
	tmp=Id*fm;	  
	G.set_les(w,tmp,fft_domain::freq);
      }
      // this transforms back the imaginary part of G(w) only
      G.integral_transform_to_time(h,FFT_TRAPEZ);
    }
  }

  // implementation:
  template <typename DOS>
  void green_equilibrium_ness(int bosefermi,herm_matrix_ness &G, DOS &dos, double beta,  double mu,
			      double h,fft_integral_method method, int limit, int nn) {
    if(bosefermi==-1){
      return green_equilibrium_ness_fermi(G,dos,beta,mu,h,method,limit,nn);
    }else{
      return green_equilibrium_ness_bose(G,dos,beta,mu,h,method,limit,nn);
    }    
  }
  

  //@private
  template <typename DosFunction>
  void green_equilibrium_ret_ness_adaptive(herm_matrix_ness &G, DosFunction &dos, double h, int nn,int limit) {
    int nt = G.Nft_/2;
    int size1 = G.retarded().size1_;
    int Ntot = 2 * nt;    
    fourier::adft_func adft;    
    auto integrand = [&](double omega) {
      return dos(omega);
    };    
    adft.sample(0.0, dos.lo_, dos.hi_, integrand, nn, limit);    
    // Compute positive times t = n h
    for (int n = 0; n < nt; ++n) {
      cplx res, err;
      double tval = h * n;
      adft.dft(-tval, res, err);
      res *= cplx(0, -1.0);  // apply prefactor -i      
      Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(size1, size1);
      for (int i = 0; i < size1; ++i) {
	M(i, i) = res;
      }      
      G.retarded().set_element(n, M, fft_domain::time);
    }    
    // Set negatice times to 0:
    Eigen::MatrixXcd zeroMat = Eigen::MatrixXcd::Zero(size1, size1);
    for (int n = nt; n < 2*nt; ++n) {
      G.retarded().set_element(n, zeroMat, fft_domain::time);
    }
  }

  //@private
  template <typename Integrand>
  void compute_lesser(fft_array &G, int nt, Integrand &integrand, double lo, double hi, double h, int nn, int limit) {
    int Nft = G.Nft_;
    int size = G.size1_;
    assert(size > 0);
    assert(nt <= Nft/2);

    G.set_zero(fft_domain::time);
    fourier::adft_func adft;
    adft.sample(10*h, lo, hi, integrand, nn, limit);
    
    // Prepare the full matrix at each time
    for (int n = 0; n < nt; ++n) {
      double t = n * h;
      std::complex<double> res, err;
      adft.dft(-t, res, err);
      res *= std::complex<double>(0, -1.0);
      
      Eigen::MatrixXcd mat(size, size);
      mat.setZero();
      mat.diagonal().setConstant(res);
      
      G.set_element(n, mat, fft_domain::time);
      
      if (n > 0) {
	Eigen::MatrixXcd mat_neg = -mat.adjoint();
	G.set_element(-n, mat_neg, fft_domain::time);
      }
    }
    
    // Set midpoint (Nft / 2) to zero matrix
    Eigen::MatrixXcd zero_mat = Eigen::MatrixXcd::Zero(size, size);
    G.set_element(Nft / 2, zero_mat, fft_domain::time);
  }
  


} // namespace ness2


#endif  // NESS_EQUILIBRIUM_DECL_H
