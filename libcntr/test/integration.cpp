#include "catch.hpp"

#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <sys/stat.h>

#include "cntr.hpp"

#define CPLX std::complex<double>
using namespace std;

TEST_CASE("integrator class", "[integrator]") {
  int k=5;
  double h = 0.1;
  double eps = 1.0e-4;
  
  SECTION("poly_interpolation"){
    integration::Integrator<double> I(k);
    dvector ft(k+1);
    
    for(int i=0;i<=k;i++) ft(i) = cos(h*i);

    double err=0.0;
    for(int i=0;i<k;i++){
      double t = (i+0.5)*h;
      double finter=0.0;

      for(int l=0;l<=k;l++){
        double t1=1.0;
        double weight = I.poly_interpolation(0,l);
        for(int n=1;n<=k;n++){
	  t1*=i+0.5;
	  weight +=t1 * I.poly_interpolation(n,l);
        }
        
        finter+=ft(l)*weight;
      }     
      err += fabs(cos(t)-finter);
    }
        
    REQUIRE(err < eps);
  }

  SECTION("poly_differentiation"){
    integration::Integrator<double> I(k);
    dvector ft(k+1);
    
    for(int i=0;i<=k;i++) ft(i) = cos(h*i);

    double err=0.0;
    for(int i=0;i<=k;i++){
      double df_exact = -sin(h*i);
      double df_approx = 0.0;
      for(int l=0;l<=k;l++) df_approx += (1.0/h) * I.poly_differentiation(i,l)*ft(l);
      err += fabs(df_exact-df_approx);
    }
    REQUIRE(err < eps);
  }
  
  SECTION("poly_integration"){
    integration::Integrator<double> I(k);
    dvector ft(k+1);
    
    for(int i=0;i<=k;i++) ft(i) = cos(h*i);

    double err=0.0;
    for(int i=1;i<=k;i++){
      for(int j=0;j<i;j++){
	double I_approx=0.0;
	for(int l=0;l<=k;l++) I_approx += h*I.poly_integration(i,j,l)*ft(l);
	double I_exact = sin(h*j) - sin(h*i);
	err += fabs(I_exact-I_approx);
      }
    }

    REQUIRE(err < eps);
    
  }

  SECTION("bd_weights"){
    integration::Integrator<double> I(k);
    double t1=0.5;

    double df_approx=0.0;
    for(int l=0;l<=k+1;l++) df_approx += I.bd_weights(l) * cos(t1-l*h)/h;
    double df_exact = -sin(t1);
      
    double err = fabs(df_exact-df_approx);
    REQUIRE(err < eps);
  }

  SECTION("rcorr"){
    integration::Integrator<double> I(k);
    int ntau=100;
    double beta=5.0, eps1=-0.1, eps2=0.2;
    double dtau=beta/ntau;
    dvector A(ntau+1), B(ntau+1);
    dvector R_exact(k);

    R_exact(1) = 0.0137659;
    R_exact(2) = 0.0274638;
    R_exact(3) = 0.0410948;
    R_exact(4) = 0.0546598;

    for(int m=0;m<=ntau;m++){
      A(m) = -cntr::fermi_exp(beta,m*dtau,-eps1);
      B(m) = -cntr::fermi_exp(beta,m*dtau,-eps2);
    }

    double err=0.0;
    for(int m=1;m<k;m++){
    
      double R_approx=0.0;
      for(int j=0;j<=k;j++){
	for(int l=0;l<=k;l++){
	  R_approx += dtau*I.rcorr(m,j,l)*A(j)*B(l);
	}
      }

      err += fabs(R_exact(m)-R_approx);
    }

    REQUIRE(err < eps);
  }

  SECTION("integrate"){
    integration::Integrator<double> I(k);
    int nt=100;
    double tmax=10.0;
    double h=tmax/nt;
    vector<double> ft(nt+1);
    
    for(int i=0;i<=nt;i++) ft[i] = cos(h*i);

    double err=0.0;
    for(int i=k+1;i<=nt;i++){
      double I_approx = h*I.integrate(ft,i);
      double I_exact = sin(h*i);
      err += fabs(I_exact-I_approx);
    }
    
    REQUIRE(err < eps);
    
  }
  
}
