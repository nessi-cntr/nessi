#ifndef _FOURIER_H
#define _FOURIER_H


#include <cmath>
#include <cassert>
#include <iostream>
#include <complex>
#include <vector>
#include <stdlib.h>



namespace fourier {

void get_dftcorr_linear(double th,double *corfac,std::complex<double> *endcor);
void get_dftcorr_cubic(double th,double *corfac,std::complex<double> *endcor);

void complex_dftcor_cubic(double w, double delta,double a,double b, 
  std::complex<double> *endpts,std::complex<double> *endcor,double *corfac);
void dft_cplx(double w,int n,double a,double b,std::complex<double> *f,
  std::complex<double> &res, std::complex<double> &err);


#define PI 3.14159265358979323846
#define ADFT_MINPTS 16

  /** \brief <b> Class `adft_func` contains functions for computing accurate Fourier transforms
    by adaptive splitting of intervals.  </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   *  Class `adft_func` contains functions for computing accurate Fourier transform
   *  by adaptive splitting of intervals. For computing the integral
   *  \f$I(\omega) = \int^b_a dt\, e^{i \omega t} f(t)\f$, cubically corrected DFT (see dft_cplx)
   *  is used. From the error estimate of dft_cplx, the integration interval \f$[a,b]\f$ can be split further
   *  until the error is below a given threshold for all interval.
   */
class adft_func{
 public:
   typedef  std::complex<double> cplx;   
   adft_func(){
      data_= 0;
      size_=0;
      tot_=0;
      nmax_=0;
      n_=0;
      a_=0;
      b_=0;
      len_=0;
      f_=0;
   }
   /// @private
   /** \brief <b> Initializes the `adft_func` class </b>
    *
    * <!-- ====== DOCUMENTATION ====== -->
    *
    *  \par Purpose
    * <!-- ========= -->
    *
    * Initializes the `adft_func` class for a given number of points `size` and max. numer of intervals `nmax`.
    *
    * <!-- ARGUMENTS
    *      ========= -->
    *
    * @param size
    * > number of points
    * @param nmax
    * > max. numer of intervals
    */    
   adft_func(int size,int nmax){
      assert(size>0 && nmax>0);
      data_=new cplx [size];
      size_=size;
      tot_=0;
      nmax_=nmax;
      n_=0;
      a_=new double [nmax];
      b_=new double [nmax];
      len_=new int [nmax];
      f_=new cplx* [nmax];
   }
   /// @private
   ~adft_func(){
      delete [] a_;
      delete [] b_;
      delete [] len_;
      delete [] data_;
      delete [] f_;
   }
   /// @private
   adft_func &operator=(const adft_func &a){
     std::cerr << "= not implemented for adft" << std::endl;
     abort();
   }
   /// @private
   adft_func(const adft_func &a){
     std::cerr << " copy constructor not implemented for for adft" << std::endl;
     abort();
   }
    /// @private
   /** \brief <b> Resizes the `adft_func` class. </b>
    *
    * <!-- ====== DOCUMENTATION ====== -->
    *
    *  \par Purpose
    * <!-- ========= -->
    *
    * Resizes the `adft_func` class for a given number of points `size` and max. numer of intervals `nmax`.
    *
    * <!-- ARGUMENTS
    *      ========= -->
    *
    * @param size
    * > number of points
    * @param nmax
    * > max. numer of intervals
    */  
   void resize(int size,int nmax){  // a resize always cancels all entries !!!
     if(size<0 || nmax<0) {std::cerr << "adft: size<0 or nmax<0" << std::endl;abort();}
     delete [] data_;
     delete  [] a_;
     delete  [] b_;
     delete  [] len_;
     delete [] f_;
     data_=new cplx [size];
     size_=size;
     tot_=0;
     nmax_=nmax;
     n_=0;
     a_=new double [nmax];
     b_=new double [nmax];
     len_=new int [nmax];   
     f_=new cplx* [nmax];
   }
   /// @private
   /** \brief <b> Returns the function values on the interval points. </b>
    *
    * <!-- ====== DOCUMENTATION ====== -->
    *
    *  \par Purpose
    * <!-- ========= -->
    *
    * Returns a pointer to the function values in the interval `i`.
    *
    * <!-- ARGUMENTS
    *      ========= -->
    *
    * @param i
    * > index of interval
    */  
   cplx *interval_data(int i){
     if(i>=n_)  {std::cerr << "interval i="<<i<<" >n="<< n_<<" in adft" << std::endl;abort();}
     return f_[i]; 
   }
   /// @private
   /** \brief <b> Initializes the internal interval data and function values. </b>
    *
    * <!-- ====== DOCUMENTATION ====== -->
    *
    *  \par Purpose
    * <!-- ========= -->
    *
    * Initializes the internal interval data and function values on the internal grid.
    *
    * <!-- ARGUMENTS
    *      ========= -->
    *
    * @param a
    * > lower interval bound
    * @param b
    * > upper interval bound
    * @param n
    * > number of points in [a,b]
    * @param fz
    * > [template class function] scalar complex function
    */  
   template < class function>
   void init(double a,double b,int n,function &fz){
     int i;
     cplx *f;
     double x,dx;
     
     if(a>=b){std::cerr << "adft_func: a>=b" << std::endl;abort();}
     if(n<1){std::cerr << "adft_func: n<1" << std::endl;abort();}
     if(nmax_<1 || size_<n) resize(n+1,1);
     n_=1;
     tot_=n+1;     
     len_[0]=n;
     a_[0]=a;
     b_[0]=b;
     f=data_;
     dx=(b-a)/((double) n);
     for(i=0;i<=n;i++){
       x=a+i*dx;
       f[i]=fz(x);
     }
     f_[0]=data_;     
   }   
   /// @private
   /** \brief <b> Split interval interv by doubling points. </b>
    *
    * <!-- ====== DOCUMENTATION ====== -->
    *
    *  \par Purpose
    * <!-- ========= -->
    *
    * Split interval interv by doubling points.
    *
    * <!-- ARGUMENTS
    *      ========= -->
    *
    * @param interv
    * > index of interval
    * @param fz
    * > [template class function] scalar complex function
    */  
   template < class function>
   void split(int interv,function &fz){
     // split interval interv by doubling points
     int j,j0,j1,jm,len0,n0=n_,i1=interv+1,tot0;
     double a0,b0,mid,x,dx;
     cplx *f0,*f1;
     
     
     if(interv>=n0 || interv<0) {std::cerr << "adft_cplx_split: interv out of range" << std::endl; abort();}
     n0=n_;
     tot0=tot_;
     a0=a_[interv];
     b0=b_[interv];
     len0=len_[interv];
     
     n_++;
     tot_ += len_[interv] + 1	; 
     if(n_ > nmax_ ||  tot_ > size_ ) {std::cerr << "adft_cplx_split: something too small" << std::endl; abort();}
  
     /*shift all intervals after interv*/
     for(j=n0-1;j>interv;j--){
        a_[j+1]=a_[j];
        b_[j+1]=b_[j];
        len_[j+1]=len_[j];
        f_[j+1]=f_[j];
     }
     /*new interval-edges*/
     mid=(a0+b0)/2.0;
     b_[interv]=mid;
     b_[i1]=b0;
     a_[i1]=mid;
     a_[interv]=a0;
     len_[i1]=len0;
     /*reshuffle data and recalculate*/
     f0=interval_data(interv);
     f1=data_+tot0;
     f_[i1]=f1;     
     dx=(b0-a0)/(2.0*len0);
     if(len0%2==0){
       f1[len0]=f0[len0];
       jm=len0/2;
       j1=len0;
       for(j0=len0-1;j0>=jm;j0--){
         j1--;
         x=mid+dx*j1;
         f1[j1]=fz(x);
         j1--;
         f1[j1]=f0[j0];
       }
       j1=len0;
       for(j0=jm;j0>0;j0--){
         f0[j1]=f0[j0];
         j1--;
         x=a0+dx*j1;
         f0[j1]=fz(x);
         j1--;
       }
     }else{
       jm=(len0+1)/2;
       j1=len0;
       for(j0=len0;j0>=jm;j0--){
         f1[j1]=f0[j0];
         j1--;
         x=mid+dx*j1;
         f1[j1]=fz(x);
         j1--;
        }
       j1=len0;
       for(j0=jm-1;j0>=0;j0--){
         x=a0+dx*j1;
         f0[j1]=fz(x);
         j1--;
         f0[j1]=f0[j0];
         j1--;
       }
     }
   }
   /** \brief <b> Computes the cubically corrected DFT. </b>
    *
    * <!-- ====== DOCUMENTATION ====== -->
    *
    *  \par Purpose
    * <!-- ========= -->
    *
    * Computes the cubically corrected DFT by integrating over the individual intervals.
    *
    * <!-- ARGUMENTS
    *      ========= -->
    *
    * @param w
    * > frequency \f$\omega\f$ in the Fourier integral
    * @param result
    * > value of the Fourier integral
    * @param err
    * > estimate of the error
    */  
   void dft(double w,cplx &result,cplx &err){
      int i,n_i;
      double a,b;
      cplx res_i,err_i,*f_i;
      if(n_<=0) {std::cerr << "adft_cplx: no interval in f" << std::endl; abort();}
      result=0;
      err=0;
  
      for(i=0;i<n_;i++){
         f_i=interval_data(i);
         a=a_[i];
         b=b_[i];
         n_i=len_[i];
         fourier::dft_cplx(w,n_i,a,b,f_i,res_i,err_i);
         result += res_i;
         err += err_i;
      }
   }
   // template <class T>
   // void sample(double w,double a,double b,T (*fz)(double),int n,int limit){

   /** \brief <b> Fourier samples the given function to determine the individual intervals. </b>
    *
    * <!-- ====== DOCUMENTATION ====== -->
    *
    *  \par Purpose
    * <!-- ========= -->
    *
    * Computes the Fourier integral for fixed frequency \f$\omega\f$. The interval \f$[a,b]\f$
    * is iteratively split until the error of the Fourier integral over each subinterval
    * is smaller or equal than the smallest error.
    *
    * <!-- ARGUMENTS
    *      ========= -->
    *
    * @param w
    * > frequency \f$\omega\f$ in the Fourier integral
    * @param a
    * > lower interval bound
    * @param b
    * > upper interval bound
    * @param fz
    * > [template class function] scalar complex function
    * @param n
    * > number of points in each interval
     * @param limit
    * > max number of intervals 
    */ 
   template <class function>
   void sample(double w,double a,double b,function &fz,int n,int limit){
       // n:      number of points in each interval
       // limit:  max number of intervals 
       int i,n_cur,isplit;
       double abserr,abserr_i;
       cplx z,res_i,err_i,res_l,err_l,res_r,err_r;
       std::vector<cplx> results(limit),errors(limit);
    
       if(n<ADFT_MINPTS || n%2!=0){ std::cerr << "adft_cplx_sample: n<16 || n%2!=0" << std::endl; abort();}
       if(limit<1){std::cerr << "adft_cplx_sample: limit<1" << std::endl; abort();}
       if(a>=b){ std::cerr << "adft_cplx_sample: a>=b" << std::endl;abort();}
       
       if(limit != nmax_ || (n+1)*limit != size_ ) resize((n+1)*limit,limit);  
       // Sample f in first interval
       init(a,b,n,fz);
       if(limit==1) return;
       // split first interval
       split(0,fz);
       // do the first fourier transform
       fourier::dft_cplx(w,n,a_[0],b_[0],interval_data(0),results[0],errors[0]);
       fourier::dft_cplx(w,n,a_[1],b_[1],interval_data(1),results[1],errors[1]);
       while(n_+1<limit){
          // look for largest error in res-table
          n_cur=n_;
          abserr=-1.0;
          isplit=-1;
          for(i=0;i<n_cur;i++){
             z=errors[i];
             abserr_i=sqrt(z.real()*z.real()+z.imag()*z.imag());
             if(abserr_i>abserr){
               isplit=i;
               abserr=abserr_i;
             }
          }
          // split interval isplit and integrate partial intervals

          //std::cout << "isplit= "<< isplit << "n= "<< n_ << std::endl; 
//          std::cout << "n " << n_ << "i " << isplit << "[ " << a_[isplit] << "," << b_[isplit] << "]" << std::endl; 
          split(isplit,fz);

          fourier::dft_cplx(w,n,a_[isplit],b_[isplit],interval_data(isplit),res_l,err_l);
          fourier::dft_cplx(w,n,a_[isplit+1],b_[isplit+1],interval_data(isplit+1),res_r,err_r);
          // update res-table
          for(i=n_cur-1;i>isplit;i--){
             results[i+1]=results[i];
             errors[i+1]=errors[i];
          }
          results[isplit+1]=res_r;
          errors[isplit+1]=err_r;
          results[isplit]=res_l;
          errors[isplit]=err_l;
       }
   }        
 private:
   /// @private
   /** \brief <b> function data </b> */
   cplx *data_; 
    /// @private
   /** \brief <b> size of data block </b> */
   int size_;    
   /// @private
   /** \brief <b> space currently used </b> */
   int tot_;    
   /// @private
   /** \brief <b> size of data block </b> */
   int nmax_;    
   /// @private
   /** \brief <b> number of intervals </b> */
   int n_;      
   /// @private
   /** \brief <b> a[0..n-1], left endpoints </b> */
   double *a_;  
   /// @private
   /** \brief <b> a[0..n-1], right endpoints </b> */
   double *b_; 
   /// @private
   /** \brief <b> len[0..n-1] number of data points in interval </b> */
   int *len_;    
   cplx **f_;
};



} // namespace


#endif      // fourier

