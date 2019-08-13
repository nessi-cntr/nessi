/*--------------------------------------------------
       fourier.c
       Martin Eckstein, Dec. 7, 2006
-----------------------------------------------------
PURPOSE:   fourier transforms of complex functions         */
#include "./fourier.hpp"
#include <cmath>

namespace fourier{

typedef std::complex<double> cplx;

/** \brief <b> Returns correction factor \f$W(\theta)\f$ and boundary
	 correction terms \f$\alpha_j(\theta)\f$ needed for cubically corrected DFT. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > For a function \f$f(t)\f$, represented on an equidistant grid \f$t_n= n h\f$, \f$n=0,\dots,N\f$,
* > we would like to compute the Fourier integral \f$I(\omega) = \int^b_a dt\, e^{i \omega t}f(t)\f$.
* > Using piecewie cubic interpolation of \f$f(t)\f$, the \f$I(\omega)\f$ can be computed by 
* > by cubically corrected discrete Fourier transformation (DFT). The algorithm is explained in
* > 
* > W. H. Press, S. A. Teukolosky, W. T. Vetterling, B. P. Flannery, Numerical Recipes 3rd Edition: 
* > The Art of Scientific Computing, Cambridge University Press, 2007, chapter 13.9
* 
* > `get_dftcorr_cubic` computes the required correction factor \f$W(\theta)\f$ and the boundary
* > terms \f$\alpha_i(\theta)\f$, \f$i=0,\dots,3\f$, where \f$\theta = \omega h\f$.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param th
* > [double] \f$\theta = \omega h\f$ as explained above.
* @param corfac
* > [double] on return, correction factor \f$W(\theta)\f$
* @param endcor
* > [complex] on return, boundary term \f$\alpha_j(\theta)\f$, stored in 'endcor[j]'
*/
void get_dftcorr_cubic(double th,double *corfac,cplx *endcor)
{
	double ai[4],ar[4],t;
	double t2,t4,t6;
	double cth,ctth,spth2,sth,sth4i,stth,th2,th4,tmth2,tth4i;
        int j;
	
	if (fabs(th) < 5.0e-2) {
		t=th;
		t2=t*t;
		t4=t2*t2;
		t6=t4*t2;
		*corfac=1.0-(11.0/720.0)*t4+(23.0/15120.0)*t6;
		ar[0]=(-2.0/3.0)+t2/45.0+(103.0/15120.0)*t4-(169.0/226800.0)*t6;
		ar[1]=(7.0/24.0)-(7.0/180.0)*t2+(5.0/3456.0)*t4-(7.0/259200.0)*t6;
		ar[2]=(-1.0/6.0)+t2/45.0-(5.0/6048.0)*t4+t6/64800.0;
		ar[3]=(1.0/24.0)-t2/180.0+(5.0/24192.0)*t4-t6/259200.0;
		ai[0]=t*(2.0/45.0+(2.0/105.0)*t2-(8.0/2835.0)*t4+(86.0/467775.0)*t6);
		ai[1]=t*(7.0/72.0-t2/168.0+(11.0/72576.0)*t4-(13.0/5987520.0)*t6);
		ai[2]=t*(-7.0/90.0+t2/210.0-(11.0/90720.0)*t4+(13.0/7484400.0)*t6);
		ai[3]=t*(7.0/360.0-t2/840.0+(11.0/362880.0)*t4-(13.0/29937600.0)*t6);
	} else {
		cth=cos(th);
		sth=sin(th);
		ctth=cth*cth-sth*sth;
		stth=2.0e0*sth*cth;
		th2=th*th;
		th4=th2*th2;
		tmth2=3.0e0-th2;
		spth2=6.0e0+th2;
		sth4i=1.0/(6.0e0*th4);
		tth4i=2.0e0*sth4i;
		*corfac=tth4i*spth2*(3.0e0-4.0e0*cth+ctth);
		ar[0]=sth4i*(-42.0e0+5.0e0*th2+spth2*(8.0e0*cth-ctth));
		ai[0]=sth4i*(th*(-12.0e0+6.0e0*th2)+spth2*stth);
		ar[1]=sth4i*(14.0e0*tmth2-7.0e0*spth2*cth);
		ai[1]=sth4i*(30.0e0*th-5.0e0*spth2*sth);
		ar[2]=tth4i*(-4.0e0*tmth2+2.0e0*spth2*cth);
		ai[2]=tth4i*(-12.0e0*th+2.0e0*spth2*sth);
		ar[3]=sth4i*(2.0e0*tmth2-spth2*cth);
		ai[3]=sth4i*(6.0e0*th-spth2*sth);
	}
	for(j=0;j<=3;j++) endcor[j]=cplx(ar[j],ai[j]);
}
/** \brief <b> Returns correction factor \f$W(\theta)\f$ and boundary
	 correction term \f$\alpha_0(\theta)\f$ needed for linearly corrected DFT. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > For a function \f$f(t)\f$, represented on an equidistant grid \f$t_n= n h\f$, \f$n=0,\dots,N\f$,
* > we would like to compute the Fourier integral \f$I(\omega) = \int^b_a dt\, e^{i \omega t}f(t)\f$.
* > Using linear interpolation of \f$f(t)\f$, the \f$I(\omega)\f$ can be computed by 
* > by linearlt corrected discrete Fourier transformation (DFT). The algorithm is explained in
* > 
* > W. H. Press, S. A. Teukolosky, W. T. Vetterling, B. P. Flannery, Numerical Recipes 3rd Edition: 
* > The Art of Scientific Computing, Cambridge University Press, 2007, chapter 13.9
* >
* > `get_dftcorr_linear` computes the required correction factor \f$W(\theta)\f$ and the boundary
* > term \f$\alpha_0(\theta)\f$, where \f$\theta = \omega h\f$.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param th
* > [double] \f$\theta = \omega h\f$ as explained above.
* @param corfac
* > [double] on return, correction factor \f$W(\theta)\f$
* @param endcor
* > [complex] on return, boundary term \f$\alpha_0(\theta)\f$
*/
void get_dftcorr_linear(double th,double *corfac,cplx *endcor)
{
	double ai,ar;
	double th2,th4,th6,cth,sth;
	
	if (fabs(th) < 5.0e-2) {
		th2=th*th;
		th4=th2*th2;
		th6=th4*th2;
		*corfac=1.0-(1.0/12.0)*th2+(1.0/360.0)*th4-(1.0/20160.0)*th6;
		ar=-0.5+th2/24.0-(1.0/720.0)*th4+(1.0/40320.0)*th6;
		ai=th*(1.0/6.0-(1.0/120.0)*th2+(1.0/5040.0)*th4-(1.0/362880.0)*th6);
	} else {
		cth=cos(th);
		sth=sin(th);
		th2=th*th;
		*corfac=2.0*(1.0-cth)/th2;
		ar=(cth-1.0)/th2;
		ai=(th-sth)/th2;
	}
	endcor[0]=cplx(ar,ai);
	endcor[1]=0;
}


/** \brief <b> Returns correction factor \f$W(\theta)\f$ and evaluate boundary
	 correction terms needed for cubically corrected DFT. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > For a function \f$f(t)\f$, represented on an equidistant grid \f$t_n= n h\f$, \f$n=0,\dots,N\f$,
* > we would like to compute the Fourier integral \f$I(\omega) = \int^b_a dt\, e^{i \omega t}f(t)\f$.
* > Using piecewie cubic interpolation of \f$f(t)\f$, the \f$I(\omega)\f$ can be computed by 
* > by cubically corrected discrete Fourier transformation (DFT). The algorithm is explained in
* > 
* > W. H. Press, S. A. Teukolosky, W. T. Vetterling, B. P. Flannery, Numerical Recipes 3rd Edition: 
* > The Art of Scientific Computing, Cambridge University Press, 2007, chapter 13.9
* >
* > `complex_dftcor_cubic` computes the required correction factor \f$W(\theta)\f$ and the boundary
* > terms \f$\alpha_j(\theta)f_j\f$, \f$j=0,\dots,3\f$, and 
* > \f$ e^{i\omega(b-a)}\alpha^*_j(\theta)f_{N-j} \f$ where \f$\theta = \omega h\f$. 
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param w
* > [double] frequency \f$\omega\f$
* @param delta
* > [double] grid spacing \f$h\f$
* @param a
* > [double] starting point of interval
* @param b
* > [double] end point of interval
* @param endpts
* > [complex] 'endpts[j]', 'j=0,1,2,3', corresponds to \f$f_j\f$, while 'endpts[7-j]', 'j=0,1,2,3', 
* >	 corresponds to \f$f_{N-j}\f$
* @param endcor
* > [complex] on return, boundary correction terms, stored in the order as 'endpts'
* @param corfac
* > [double] on return, correction factor \f$W(\theta)\f$
*/
void complex_dftcor_cubic(double w, double delta,double a,double b, cplx *endpts,cplx *endcor,double *corfac)
{
	double ai[4],ar[4],arg,t;
	double t2,t4,t6;
	double cth,ctth,spth2,sth,sth4i,stth,th,th2,th4,tmth2,tth4i;
        cplx expfac,cr,cl,cor,al;
	int j;
	
	th=w*delta;
	if (fabs(th) < 5.0e-2) {
		t=th;
		t2=t*t;
		t4=t2*t2;
		t6=t4*t2;
		*corfac=1.0-(11.0/720.0)*t4+(23.0/15120.0)*t6;
		ar[0]=(-2.0/3.0)+t2/45.0+(103.0/15120.0)*t4-(169.0/226800.0)*t6;
		ar[1]=(7.0/24.0)-(7.0/180.0)*t2+(5.0/3456.0)*t4-(7.0/259200.0)*t6;
		ar[2]=(-1.0/6.0)+t2/45.0-(5.0/6048.0)*t4+t6/64800.0;
		ar[3]=(1.0/24.0)-t2/180.0+(5.0/24192.0)*t4-t6/259200.0;
		ai[0]=t*(2.0/45.0+(2.0/105.0)*t2-(8.0/2835.0)*t4+(86.0/467775.0)*t6);
		ai[1]=t*(7.0/72.0-t2/168.0+(11.0/72576.0)*t4-(13.0/5987520.0)*t6);
		ai[2]=t*(-7.0/90.0+t2/210.0-(11.0/90720.0)*t4+(13.0/7484400.0)*t6);
		ai[3]=t*(7.0/360.0-t2/840.0+(11.0/362880.0)*t4-(13.0/29937600.0)*t6);
	} else {
		cth=cos(th);
		sth=sin(th);
		ctth=cth*cth-sth*sth;
		stth=2.0e0*sth*cth;
		th2=th*th;
		th4=th2*th2;
		tmth2=3.0e0-th2;
		spth2=6.0e0+th2;
		sth4i=1.0/(6.0e0*th4);
		tth4i=2.0e0*sth4i;
		*corfac=tth4i*spth2*(3.0e0-4.0e0*cth+ctth);
		ar[0]=sth4i*(-42.0e0+5.0e0*th2+spth2*(8.0e0*cth-ctth));
		ai[0]=sth4i*(th*(-12.0e0+6.0e0*th2)+spth2*stth);
		ar[1]=sth4i*(14.0e0*tmth2-7.0e0*spth2*cth);
		ai[1]=sth4i*(30.0e0*th-5.0e0*spth2*sth);
		ar[2]=tth4i*(-4.0e0*tmth2+2.0e0*spth2*cth);
		ai[2]=tth4i*(-12.0e0*th+2.0e0*spth2*sth);
		ar[3]=sth4i*(2.0e0*tmth2-spth2*cth);
		ai[3]=sth4i*(6.0e0*th-spth2*sth);
	}
	cl=0;
	for(j=0;j<=3;j++){
	   al=cplx(ar[j],ai[j]);
	   cor=endpts[j]*al;
	   cl+=cor;
	}
	cr=0;
	for(j=0;j<=3;j++){
	   al=cplx(ar[j],-ai[j]);
	   cor=endpts[7-j]*al;
	   cr += cor;
	}
	arg=w*(b-a);
	expfac=cplx(cos(arg),sin(arg));
	cr *= expfac;
	*endcor = cl+cr;
}

/** \brief <b> Computes the Fourier integral \f$I(\omega) = \int^b_a dt\, e^{i \omega t}f(t)\f$ using
	cubically corrected DFT.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*   \par Purpose
* <!-- ========= -->
*
* > For a function \f$f(t)\f$, represented on an equidistant grid \f$t_j= j h\f$, \f$j=0,\dots,n\f$,
* > we would like to compute the Fourier integral \f$I(\omega) = \int^b_a dt\, e^{i \omega t}f(t)\f$.
* > Using piecewie cubic interpolation of \f$f(t)\f$, the \f$I(\omega)\f$ can be computed by 
* > by cubically corrected discrete Fourier transformation (DFT). The algorithm is explained in
* > 
* > W. H. Press, S. A. Teukolosky, W. T. Vetterling, B. P. Flannery, Numerical Recipes 3rd Edition: 
* > The Art of Scientific Computing, Cambridge University Press, 2007, chapter 13.9
* >
* > `dft_cplx` computes the integral  \f$I(\omega)\f$, using the cubic correction factor and boundary
* > correction terms computed by `complex_dftcor_cubic`. Furthermore, an approximate error is returned, 
* > given the difference using `n` and `(n+1)/2` points.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param w
* > [double] frequency \f$\omega\f$
* @param n
* > [int] number of time points \f$t_j=j h, j=0,\dots,n\f$
* @param a
* > [double] starting point of interval
* @param b
* > [double] end point of interval
* @param f
* > [complex] pointer containing the function values 
* @param res
* > [complex] Fourier integral \f$I(\omega)\f$
* @param err
* > [complex] difference using 'n' and '(N+1)/2' points.
*/
void dft_cplx(double w,int n,double a,double b,cplx *f,cplx &res,cplx &err)
{
  /*same as fourier transform complex_ft_cubic, but
  an approximate error is returned, given the difference between 
  complex_ft_cubic using m and (m+1)/2 points*/
  double theta,delta,len,corfac,corfac2,arg;
  cplx endpoints[8],endcor,endcor2,res2,expfac,z;
  int j;
  
  len=b-a;
  if(n<=16  || (n/2)*2!=n  || len<=0.0) std::cerr << "n=" << n << " len= " << len << std::endl;
  if(n<=16  || (n/2)*2!=n  || len<=0.0) {std::cerr << "complex_ft_cubic_err: wrong input"  << std::endl;abort();}
  delta=len/((double) n);
  theta=w*delta;
  /*endcorrections using all points*/
  /*store endpoints*/
  for(j=0;j<=3;j++){
    endpoints[j]=f[j];
    endpoints[7-j]=f[n-j];
  }
  /*get corrections*/
  complex_dftcor_cubic(w,delta,a,b,endpoints,&endcor,&corfac);
  /*endcorrections using every second points*/
  /*store endpoints*/
  for(j=0;j<=3;j++){
    endpoints[j]=f[2*j];
    endpoints[7-j]=f[n-2*j];
  }
  /*get corrections*/
  complex_dftcor_cubic(w,2.0*delta,a,b,endpoints,&endcor2,&corfac2);
  /*compute discrete fourier transform 
  (naive version, without FFT)*/
  /*first using every second point*/
  res2=0;
  for(j=0;j<=n;j+=2){
    arg=((double) j)*theta;
    expfac=cplx(cos(arg),sin(arg));
    res2 += expfac*f[j];
  }
  /*using remaining points*/
  res=res2;
  for(j=1;j<n;j+=2){
    arg=((double) j)*theta;
    expfac=cplx(cos(arg),sin(arg));
    res += expfac*f[j];
  }
  /*apply corrections*/
  res *= corfac;
  res += endcor;
  arg=w*a;
  expfac=cplx(delta*cos(arg),delta*sin(arg));
  res *= expfac;
  
  res2 *= corfac2;
  res2 += endcor2;
  expfac *=2;
  res2 *=  expfac;
  err=res-res2;
}

} //namespace
