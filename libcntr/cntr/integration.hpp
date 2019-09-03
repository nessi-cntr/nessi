/*********************************************************
 *  Martin Eckstein, 2010
 *  Gregory integration rules, 
 *  and backward differntiation formulae
 *********************************************************/
#ifndef INTEGRATION
#define INTEGRATION

#include <cmath>
#include <cassert>
#include <iostream>
#include <complex>
#include <vector>
#include <stdlib.h>

namespace integration{

#define GREGORY_KMAX 8
  void read_poly_interpolation(int k,long double *P); // P must be of size k1*k1
  void read_poly_differentiation(int k,long double *w); // w must be of size  k1*k1 
  void read_poly_integration(int k,long double *w); // w must be of size  k1*k1*k1 
  void read_poly_interpolation(int k,double *P); // P must be of size k1*k1
  void read_poly_differentiation(int k,double *P); // P must be of size k1*k1
  void read_poly_integration(int k,double *P); // P must be of size k1*k1*k1
  void read_poly_interpolation(int k,float *P); // P must be of size k1*k1
  void read_poly_differentiation(int k,float *P); // P must be of size k1*k1
  void read_poly_integration(int k,float *P); // P must be of size k1*k1*k1
  void read_gregory_weights(int k,double *w);
  void read_bd_weights(int k,double *w);
  void read_rcorr(int k,double *w);

  /* #######################################################################################
     #
     #  Integrator class
     # 
     #  The integrator contains all kinds of weights for integration and differntiation
     #  of a function at equally spaced abscissae, up to a certain order k of accuray
     #  
     #  !!!!  OUT OF RANGE ERRORS ARE NOT CATCHED FOR THE ARRAYS BELOW !!!!
     #
     #  (1) GREGORY WEIGHTS: I.gregory_weight and I.gregory_omega
     #
     #      The Gregory integration formula of order k is defined by
     #
     #         int_0^m f(x) dx  =  sum_{l=0}^{ max(k,m) } f_l I.gregory_weights(m,l)
     #
     #      for m>2*k, the weights have a simpler form:
     #        gregory_weights(m,j)= I.gregory_omega(j)  j<=k
     #        gregory_weights(m,j)= 1 k < j < m-k
     #        gregory_weights(m,j)= I.gregory_omega(m-j)  j>=m-k
     #      I.gregory_omega returns for 0 <= j <= k
     #
     #
     #  (2) RCORR:  I.rcorr
     #  
     #      These are weights for the special boundary term, which is needed for
     #      the convolution  of matsubara Greenfunctions, which cannot be 
     #      contuinuously be continued to x<0 and x>beta: 
     #
     #         For 0 < m < k:     
     #         R=int_0^m dx a(m-x) b(x) = sum_{j,l=0}^k I.rcorr(m,j,l) a(j)*b(l)
     #
     #      The formula should be used if both a(x) and b(x) are known only for 
     #      x>0. The integral is computed by approximating a(m-x) in the 
     #      interval [m-k,m] and b(x) in the interval [0,k], and integrating the 
     #      product of the approxiating polynomials. 
     #
     #
     #  (3) BACKWARD DIFFERENTIATION:  I.bd_weight
     #
     #      Backward differentiation formula  !!!! OF ORDER K+1 !!!!
     #
     #         d/dx f(x=x0) = sum_{l=0}^{k+1} I.bd_weights(l)*f(x0-l)
     #
     #      The weights are related to the poly_differentiation weights (see 4a),
     #      I.bd_weights(l)=I._poly_differentiation(k+1,k+1-l)
     #      but (...k+1  are not stored)
     #
     #  (4a) POLYNOMIAL INTERPOLATION: I.poly_interpolation
     #
     #      Provided for completeness only: 
     #      The k-th order approximation polynomial trough function 
     #      values {(x_l,f(x_l)),l=0...k} is given by
     #
     #        P(x) = sum_{alpha,l=0}^{k} x^alpha*I.interpolation(alpha,l)*f_l
     #
     #  (4b) POLYNOMIAL DIFFERENTIATION: I.poly_differentiation
     # 
     #      To compute the derivative of the kth order polynomial
     #      through points (x_l,f_l), l=0...k:
     #         
     #         d/dx f(x=x_i) = sum_{l=0}^k  I.poly_differentiation(i,l)*f(l)
     #
     #
     #  (4c) POLYNOMIAL INTEGRATION:  I.poly_integration
     #
     #      These are weights to compute the integral
     #
     #        For 0 <= i,j <= k 
     #        int_i^j dx f(x) = sum_{l=0}^k  I.poly_integration(i,j,l)*f(l)
     #  
     #      This is needed for starting formulas of the volterra Integration. 
     #      The integral is ontained by computing the approximating polynomial of
     #      f(x) in [0,k] and integrating it.
     #
     ####################################################################################### */



  /** \brief <b> Class `Integrator` contains all kinds of weights for integration and differentiation
   * of a function at equally spaced abscissae, up to a certain order k of accuray. </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   *  Class `Integrator` contains all kinds of weights for integration and differentiation
   *  of a function at equally spaced abscissae, up to a certain order k of accuray. In particular,
   *  the class contains coefficients for
   *   - Gregory quadrature
   *   - polynomial quadrature
   *   - polynomial interpolation
   *   - polynomial differentiation
   *   - coefficients for the backwards differencing formula (BDF)
   *
   */
  template <typename T=double> class Integrator {
public:
  /** \brief <b> Initializes the `Integrator` class for a given order k. </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * Initializes the `Integrator` class for given order \f$k\f$ between 1 and 8. If \f$k>5\f$,
   * Gregory integration up to order \f$k\f$ is used, while all other operations are limited to 
   * maximum order 5.
   *
   * <!-- ARGUMENTS
   *      ========= -->
   *
   * @param k
   * > The order of the weights.
   */
  Integrator(int k=0){
  int k1=k+1;
  if(k<0 || k > GREGORY_KMAX ){ std::cout << "Integrator: k out of range " << std::endl; abort();}
  poly_interpolation_= new T [k1*k1];
  poly_differentiation_= new T [k1*k1];
  poly_integration_= new T [k1*k1*k1];
  bd_weights_= new T [k1+1]; 
  gregory_weights_= new T [ 4*k1*k1 ];
  gregory_omega_=gregory_weights_ + ((2*k+1)*2*k1); /*pointer to last line*/
  if(k>1) rcorr_ = new T [ (k-1)*k1*k1] ; else rcorr_=0;
  k_=k;
  integration::read_poly_interpolation(k_,poly_interpolation_); 
  integration::read_poly_differentiation(k_,poly_differentiation_); 
  integration::read_poly_integration(k_,poly_integration_); 
  integration::read_bd_weights(k+1,bd_weights_);
  integration::read_gregory_weights(k,gregory_weights_); 
  integration::read_rcorr(k_,rcorr_);
}	 
  Integrator &operator=(const Integrator& I){
  int k1;
  if( this->k_==I.k_ ) return *this;
  delete [] poly_interpolation_;
  delete [] poly_differentiation_;
  delete [] poly_integration_;
  delete [] bd_weights_;
  delete [] gregory_weights_;
  delete [] rcorr_;
  k_=I.k_;
  k1=k_+1;
  poly_interpolation_= new T [k1*k1];
  poly_differentiation_= new T [k1*k1];
  poly_integration_= new T [k1*k1*k1];
  bd_weights_= new T [k1+1]; 
  gregory_weights_= new T [ 4*k1*k1 ];
  gregory_omega_=gregory_weights_ + ((2*k_+1)*2*k1); /*pointer to last line*/
  if(k_>1) rcorr_ = new T [ (k_-1)*k1*k1] ; else rcorr_=0;
  integration::read_poly_interpolation(k_,poly_interpolation_); 
  integration::read_poly_differentiation(k_,poly_differentiation_); 
  integration::read_poly_integration(k_,poly_integration_); 
  integration::read_bd_weights(k_+1,bd_weights_);
  integration::read_gregory_weights(k_,gregory_weights_); 
  integration::read_rcorr(k_,rcorr_);
  return *this;
}
  ~Integrator(){
  delete [] poly_interpolation_;
  delete [] poly_differentiation_;
  delete [] poly_integration_;
  delete [] bd_weights_;
  delete [] gregory_weights_;
  delete [] rcorr_;
}
  Integrator(const Integrator &I){
  int k1;
  k_=I.k_;
  k1=k_+1;
  poly_interpolation_= new T [k1*k1];
  poly_differentiation_= new T [k1*k1];
  poly_integration_= new T [k1*k1*k1];
  bd_weights_= new T [k1+1]; 
  gregory_weights_= new T [ 4*k1*k1 ];
  gregory_omega_=gregory_weights_ + ((2*k_+1)*2*k1); /*pointer to last line*/
  if(k_>1) rcorr_ = new T [ (k_-1)*k1*k1] ; else rcorr_=0;
  integration::read_poly_interpolation(k_,poly_interpolation_); 
  integration::read_poly_differentiation(k_,poly_differentiation_); 
  integration::read_poly_integration(k_,poly_integration_); 
  integration::read_bd_weights(k_+1,bd_weights_);
  integration::read_gregory_weights(k_,gregory_weights_); 
  integration::read_rcorr(k_,rcorr_);
}
   /** \brief <b> Returns the order \f$k\f$ of the integrator class.  </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * Returns the order \f$k\f$ of the integrator class.
   */
  int get_k(void) {return k_;}
   /** \brief <b> Returns the order \f$k\f$ of the integrator class.  </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * Returns the order \f$k\f$ of the integrator class.
   */
  int k(void) {return k_;}
  
  /** \brief <b> Returns the the weight needed for polynomial interpolation.  </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * The \f$k\f$-th order approximation polynomial \f$P(x)\f$ through function values
   * \f$\{(x_l,f(x_l),l=0,\dots,k)\}\f$ with \f$x_{l+1}-x_{l}=h\f$ is given by
   *
   * \f{align*}{
   *   P(x) = \sum^k_{\alpha=0}\sum^k_{l=0} c_{\alpha,l}\, \left(\frac{x}{h}\right)^\alpha  f(x_l)
   * \f}
   * 
   * `poly_interpolation` returns the coeffient \f$c_{\alpha,l}\f$ for given 
   * indicies \f$\alpha,l\f$.
   *
   * <!-- ARGUMENTS
   *      ========= -->
   * @param alpha
   * > First index of coefficients.
   * @param l
   * > Second index of coefficients.
   */
  T poly_interpolation(int alpha,int l) {return  poly_interpolation_[alpha*(k_+1)+l];}


   /** \brief <b> Returns the the weight needed for polynomial differentiation.  </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * Assuming the function \f$f(x)\f$ is given through points \f$\{(x_l,f(x_l),l=0,\dots,k)\}\f$\
   * with \f$x_{l+1}-x_{l}=h\f$,
   * polynomial differentiation is defined by
   *
   * \f{align*}{
   *   \frac{d}{dx}f(x)\big|_{x_i} = \frac{1}{h}\sum^k_{l=0} c_{i,l}\, f(x_l)
   * \f} 
   * 
   * `poly_differentiation` returns the coefficient \f$c_{i,l}\f$ for given 
   * indicies \f$i,l\f$.
   *
   * <!-- ARGUMENTS
   *      ========= -->
   * @param i
   * > First index of coefficients.
   * @param l
   * > Second index of coefficients.
   */
  T poly_differentiation(int i,int l) {return  poly_differentiation_[i*(k_+1)+l];}

  /** \brief <b> Returns the the weight needed for polynomial integration.  </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * Assuming the function \f$f(x)\f$ is given through points \f$\{(x_l,f(x_l),l=0,\dots,k)\}\f$,
   * polynomial interpolation is defined by
   *
   * \f{align*}{
   *   \int^{x_j}_{x_i}dx \, f(x) = h \sum^k_{l=0} w_{i,j,l}\, f(x_l)
   * \f} 
   *
   * This is needed for computing integrals over small intervals with less grid points
   * than needed for Gregory quadrature. Polynomial integration is used for the 
   * polynomial collocation method for solving the Volterra integral equations.
   *
   * ` poly_integration` returns the weight \f$w_{i,j,l}\f$ for given indices.
   *
   * <!-- ARGUMENTS
   *      ========= -->
   * @param i
   * > Index of lower integral bound.
   * @param j
   * > Index of upper integral bound.
   * @param l
   * > Index of weight.
   */
  T poly_integration(int i,int j,int l) {return  poly_integration_[i*(k_+1)*(k_+1) + j*(k_+1) + l];}

  
  /** \brief <b> Returns the the backwards differencing coefficients.  </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * Backward differentiation formula of order \f$k+1\f$ is given by
   *  
   * \f{align*}{
   *   \frac{d}{d x}f(x)\big|_{x_i} = \frac{1}{h}\sum^{k+1}_{l=0} c_l\, f(x_{i-l}) 
   * \f} 
   *
   * The weights are related to the `poly_differentiation` weights.
   *
   * `bd_weights` returns the coefficient \f$c_l\f$ for given \f$l\f$.
   *
   * <!-- ARGUMENTS
   *      ========= -->
   * @param l
   * > Index of weight.
   */
  T bd_weights(int l) {return bd_weights_[l];}

  /** \brief <b> Returns the Gregory weights for integration.  </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * The Gregory integration formula of order \f$k\f$ is defined by
   *
   * \f{align*}{
   *   \int^{x_m}_0 dx \, f(x) = \sum^{\mathrm{max}(k,m)}_{l=0} w_{m,l}\, f(x_l)
   * \f}
   * 
   * For \f$m > 2 k\f$, the weights have a simpler form:
   *
   * \f{align*}{
   *  w_{m,j} = \begin{cases} 
   *  \omega_j & : j \le k \\
   *   1 & : k < j < m-k \\
   *  \omega_{m-j} & : j \ge m-k
   *  \end{cases}
   * \f} 
   *
   * Here, \f$\omega_j\f$ are tabulated weights for \f$0 \le j \le k\f$.
   * 
   * `gregory_weights` returns \f$w_{n,j}\f$ for given \f$n, j\f$.
   * <!-- ARGUMENTS
   *      ========= -->
   *
   * @param n
   * > Integration runs up to \f$x_n\f$.
   * @param j
   * > Index of the Gregory weight.
   */
  T gregory_weights(int n,int j){
  int k1=k_+1,k2=2*k1;
  assert( j>=0 && ( j<=k_ || j<=n) );
  if( n>=k2 ){
  if( j > k_ && j < n-k_ ) { return 1;}
  else if ( j<k1 ) { return gregory_omega_[j];}
  else { return gregory_omega_[n-j]; }
}else{ return gregory_weights_[n*k2 + j];}
}

  /** \brief <b> Returns the Gregory weights at the integral boundaries (see `gregory_weights`).  </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * Returns the weight \f$\omega_j\f$ needed close to the integral boundaries (see `gregory weights`).
   * <!-- ARGUMENTS
   *      ========= -->
   * @param j
   * > Index of the Gregory weight.
   */
  T gregory_omega(int j) {return gregory_omega_[j];}


  /** \brief <b> Returns the special quadrature weights for computing integrals on the Matsubara axis.  </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * For calculation integrals on the Matsubara axis of the type
   * \f{align*}{
   *  R(\tau_m) = \int^{\tau_m}_0 d\tau\, A(\tau_m-\tau)B(\tau) \ ,
   * \f} 
   * Gregory quadrature can only be used for \f$m \ge k\f$. For \f$m < k\f$,
   * special weights are needed, since the Matsubara Green's functions can 
   * not contuinuously be continued to \f$\tau <0\f$ and \f$\tau >\beta\f$.
   * Using two-dimensional polynomial interpolation, one can represent the 
   * convolution integral by
   * \f{align*}{
   *  R(\tau_m) = \Delta\tau\sum^k_{j,l=0} w_{m,j,l}\, A(\tau_j) B(\tau_l) \ .
   * \f} 
   * `rcorr` returns the quadrature weights \f$w_{m,j,l}\f$ for given indices.
   *
   * <!-- ARGUMENTS
   *      ========= -->
   *
   * @param m
   * > The index of the upper integration bound.
   * @param j
   * > The second indes of the weights.
   * @param l
   * > The second indes of the weights.
   */
  T rcorr(int m,int j,int l){ 
  int k1=k_+1;
  if(k_<2) return 0; 
  assert(m<k_ && j<=k_ && l<=k_);
  return *(rcorr_ + (m-1)*k1*k1 + j*k1 +l);
}
  
  // The following is a particulary useless function,
  // more or less provided as an example
  // compute the integral int_j^n dx of a function f(x) - 
  // in the class M, M=0 must make sense :

  /** \brief <b> Integrates a function using Gregory quadrature.  </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * Computes the integral
   *
   * \f{align*}{
   *   \int^{x_n}_0 dx \, f(x) = h\sum^{\mathrm{max}(k,m)}_{l=0} w_{n,l}\, f(x_l)
   * \f}
   * 
   * The function \f$f\f$ represents any class `M` for which `M=0` is defined. 
   * 
   * <!-- ARGUMENTS
   *      ========= -->
   *
   * @param f
   * > A standard vector of class `M` with at least lenght `n+1`.
   * @param n
   * > The index of the upper integral bound.
   */
  template <class M> M integrate(const std::vector<M> &f,int n){
  int k1=k_+1,i,k2=2*k1,n1=n-k_;
  T *w;
  M sum=f[0];
  sum=0;

  n1=(n<=k_ ? k_ : n);   
  assert( int(f.size()) >= n1); 
   
  if( n>=k2 ){
  n1=n-k_;
  for(i=0;i<=k_;i++) sum += f[i]*gregory_omega_[i];
  for(i=k1;i<n1;i++) sum += f[i];
  for(i=n1;i<=n;i++) sum += f[i]*gregory_omega_[n-i];
}else if(n>0){
  w=gregory_weights_ + n*k2;
  for(i=0;i<=n1;i++) sum += f[i]*w[i];
}
  return sum;
}
  
private:
  int k_;
  T *poly_interpolation_;
  T *poly_differentiation_;
  T *poly_integration_;
  T *bd_weights_;
  T *gregory_weights_;
  T *gregory_omega_;
  T *rcorr_;
};

  template<typename T> Integrator<T> &I(int k){
  static Integrator<T> Isave[6];
  static int init=0;
  if(!init){
  init=1;
  for(int k=1;k<=5;k++) Isave[k]=Integrator<T>(k);
}
  return Isave[k%6];
}

#ifndef CNTR_NO_EXTERN_TEMPLATES
  extern template class Integrator<double>;
  extern template Integrator<double> &I<double>(int k);
#endif

} //namespace

#endif


