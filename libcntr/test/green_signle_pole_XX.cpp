#include "catch.hpp"
#include "cntr.hpp"
#include "vector"
#include <iostream>

#define GREEN cntr::herm_matrix<double> 
#define CPLX std::complex<double>

using namespace std;

//------------------------------------
const complex<double> xj(0.0,1.0);
double f_boson(double w0,double beta){
  
  return 1.0/(exp(beta*w0)-1.0);
}

//---------------------------------------------------------------------------
//  Explicit expression for phonon propagator -i<X(t)X(t')> with X=a+a^\dagger
//---------------------------------------------------------------------------
double D0_ph_matsu_t(double w0, double beta, double tau){
  
  double value;
  value=-(1.0+f_boson(w0,beta))*exp(-w0*tau)-f_boson(w0,beta)*exp(w0*tau);
  
  return value;
}

complex<double> D0_ph_ret(double w0, double t1,double t2){
  
  complex<double> value;
  value=(-2.0)*sin(w0*(t1-t2))*( t1>=t2 ? 1.0 : 0.0);
  
  return value;
}

complex<double> D0_ph_less(double w0, double beta, double t1,double t2){
  
  complex<double> value;
  
  value=-xj*f_boson(w0,beta)*exp(-xj*w0*(t1-t2))
    -xj*(1.0+f_boson(w0,beta))*exp(xj*w0*(t1-t2));
  
  return value;
}

complex<double> D0_ph_gtr(double w0, double beta, double t1,double t2){
  
  complex<double> value;
  
  value=-xj*(1.0+f_boson(w0,beta))*exp(-xj*w0*(t1-t2))
    -xj*f_boson(w0,beta)*exp(xj*w0*(t1-t2));
  
  return value;
}

complex<double> D0_ph_left_mix(double w0, double beta, double t1,double tau2){
  
  complex<double> value;
  
  value=-xj*f_boson(w0,beta)*exp(tau2*w0)*exp(-xj*w0*t1)
    -xj*(1.0+f_boson(w0,beta))*exp(-tau2*w0)*exp(xj*w0*t1);
  
  return value;
}

//------------------------------------
TEST_CASE("green_signel_pole_XX","[green_signel_pole_X]"){
 
  int nt,ntau,tstp;
  double w0,beta,dt,dtau;
  GREEN D0,D0_ref;
  double err=0.0;
  double eps=1e-7;
 
  nt=100;
  ntau=100;
  w0=1.4;
  beta=20.0;
  dt=0.01;
  dtau=beta/(double)ntau;
  
  D0 = GREEN(nt,ntau,1,1);
  cntr::green_single_pole_XX(D0,w0,beta,dt);

  //set reference value
  D0_ref = GREEN(nt,ntau,1,1);
  cdmatrix tmp(1,1);
  //matsubara
  for(int i=0;i<=ntau;i++){
    double tau=dtau*(double)i;
    tmp(0,0)=D0_ph_matsu_t(w0, beta, tau)/2.0;
    D0_ref.set_mat(i,tmp);
  }
  //left_mixing
  for(int i2=0;i2<=ntau;i2++){
    double tau2=dtau*(double)i2;
    for(int i1=0;i1<=nt;i1++){
      double t1=dt*(double)i1;
      tmp(0,0)=D0_ph_left_mix(w0, beta, t1,tau2)/2.0;
      D0_ref.set_tv(i1,i2,tmp);
    }
  }
  //lesser
  for(int i2=0;i2<=nt;i2++){
    double t2=dt*(double)i2;
    for(int i1=0;i1<=nt;i1++){
      double t1=dt*(double)i1;
      if(i1<=i2){
	tmp(0,0)=D0_ph_less(w0, beta, t1, t2)/2.0;
	D0_ref.set_les(i1,i2,tmp);
      }
    }
  }

  //retareded
  for(int i2=0;i2<=nt;i2++){
    double t2=dt*(double)i2;
    for(int i1=0;i1<=nt;i1++){
      double t1=dt*(double)i1;
      if(i1>=i2){
	tmp(0,0)=D0_ph_ret(w0, t1, t2)/2.0;
	D0_ref.set_ret(i1,i2,tmp);
      }
    }
  }
  
  SECTION("green_signel_pole_XX"){
    err=0.0;
    for(tstp=-1;tstp<=nt;tstp++){
      err+=cntr::distance_norm2(tstp,D0_ref,D0);
    }
    REQUIRE(err<eps);
  }
}
