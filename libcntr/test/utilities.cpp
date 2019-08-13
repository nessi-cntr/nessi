#include "catch.hpp"
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

#include "cntr.hpp"

#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define CPLX std::complex<double>
#define CFUNC cntr::function<double>
using namespace std;

/*///////////////////////////////////////////////////////////////////////////////////////


  This is a minimal test program to test cntr_utilities.

///////////////////////////////////////////////////////////////////////////////////////*/

//--------------------------------------------------------------------------------
void LagrangeWeights(int p, double h, double x, vector<double> &w){
  for(int i=0; i<=p; i++){
    w[i] = 1.0;
    for(int j=0; j<=p; j++){
      if(i != j){
	w[i] *= (x-h*j)/(i*h-h*j);
      }
    }
  }
}
//--------------------------------------------------------------------------------
void extrapolate_tv(int tstp, int p, double h, GREEN &G){
  int size1=G.size1();
  vector<double> wlagr(p+1);
  LagrangeWeights(p, h, (p+1)*h, wlagr);
  cdmatrix G_tv(size1,size1),G_tv_new(size1,size1);

  LagrangeWeights(p, h, (p+1)*h, wlagr);

  for(int itau=0; itau<=G.ntau(); itau++){
    G_tv_new.setZero();
    for(int k=0; k<=p; k++){
      G.get_tv(tstp-k,itau,G_tv);
      G_tv_new += wlagr[p-k]*G_tv;
    }
    G.set_tv(tstp+1,itau,G_tv_new);
  }

}
//--------------------------------------------------------------------------------
void extrapolate_les(int tstp, int p, double h, GREEN &G){
  int size1=G.size1();
  vector<double> wlagr(p+1);
  LagrangeWeights(p, h, (p+1)*h, wlagr);
  cdmatrix G_les(size1,size1),G_les_new(size1,size1);

  LagrangeWeights(p, h, (p+1)*h, wlagr);

  for(int j=0; j<=tstp; j++){
    G_les_new.setZero();
    for(int k=0; k<=p; k++){
      G.get_les(j,tstp-k,G_les);
      G_les_new += wlagr[p-k]*G_les;
    }
    G.set_les(j,tstp+1,G_les_new);
  }

  G_les_new.setZero();
  for(int k=0; k<=p; k++){
    G.get_les(tstp-k,tstp+1,G_les);
    G_les_new += wlagr[p-k]*G_les;
  }
  G.set_les(tstp+1,tstp+1,G_les_new);

}
//--------------------------------------------------------------------------------
void extrapolate_diag(int tstp, int p, double h, GREEN &G){
  int size1=G.size1();
  vector<double> wlagr(p+1);
  LagrangeWeights(p, h, (p+1)*h, wlagr);
  cdmatrix G_les(size1,size1),G_les_new(size1,size1);

  LagrangeWeights(p, h, (p+1)*h, wlagr);

  G_les_new.setZero();
  for(int k=0; k<=p; k++){
    G.get_les(tstp-k,tstp-k,G_les);
    G_les_new += wlagr[p-k]*G_les;
  }
  G.set_les(tstp+1,tstp+1,G_les_new);

}
//--------------------------------------------------------------------------------
void extrapolate_ret(int tstp, int p, double h, GREEN &G){
  int size1=G.size1();
  vector<double> wlagr(p+1);
  LagrangeWeights(p, h, (p+1)*h, wlagr);
  cdmatrix G_gtr(size1,size1),G_les(size1,size1),G_ret_new(size1,size1);

  LagrangeWeights(p, h, (p+1)*h, wlagr);

  for(int j=0; j<=tstp; j++){
    G_ret_new.setZero();
    for(int k=0; k<=p; k++){
      G.get_gtr(tstp-k,j,G_gtr);
      G.get_les(tstp-k,j,G_les);
      G_ret_new += wlagr[p-k]*(G_gtr-G_les);
    }
    G.set_ret(tstp+1,j,G_ret_new);
  }

}
//--------------------------------------------------------------------------------
void extrapolate_function(int tstp, int p, double h, CFUNC &F){
  int size1=F.size1();
  vector<double> wlagr(p+1);
  cdmatrix f1(size1,size1),fk(size1,size1);

  LagrangeWeights(p, h, (p+1)*h, wlagr);

  f1.setZero();
  for(int k=0; k<=p; k++){
    F.get_value(tstp-k,fk);
    f1 += wlagr[p-k]*fk;
  }
  F.set_value(tstp+1,f1);

}
//--------------------------------------------------------------------------------
cdmatrix interpolate_function(int tstp, int p, double h, double t, CFUNC &F){
  int size1=F.size1();
  vector<double> wlagr(p+1);
  LagrangeWeights(p, h, t - h*(tstp-p), wlagr);
  cdmatrix f1(size1,size1),fk(size1,size1);

  f1.setZero();
  for(int k=0; k<=p; k++){
    F.get_value(tstp-k,fk);
    f1 += wlagr[p-k]*fk;
  }
  return f1;

}
//--------------------------------------------------------------------------------
double dist_green(int tstp, GREEN &G1, GREEN &G2){
  int size=G1.size1();
  cdmatrix G1_mat(size,size),G2_mat(size,size);
  cdmatrix G1_tv(size,size),G2_tv(size,size);
  cdmatrix G1_les(size,size),G2_les(size,size);
  cdmatrix G1_ret(size,size),G2_ret(size,size);
  double dist=0.0;

  if(tstp == -1) {
    for(int m=0; m<=G1.ntau(); m++){
      G1.get_mat(m,G1_mat);
      G2.get_mat(m,G2_mat);
      dist += (G1_mat-G2_mat).norm();
    }
  } else {
    for(int m=0; m<=G1.ntau(); m++){
      G1.get_tv(tstp,m,G1_tv);
      G2.get_tv(tstp,m,G2_tv);
      dist += (G1_tv-G2_tv).norm();
    }
    for(int j=0; j<=tstp; j++){
      G1.get_les(j,tstp,G1_les);
      G2.get_les(j,tstp,G2_les);
      dist += (G1_les-G2_les).norm();
      G1.get_ret(tstp,j,G1_ret);
      G2.get_ret(tstp,j,G2_ret);
      dist += (G1_ret-G2_ret).norm();
    }
  }
  return dist;
}
//--------------------------------------------------------------------------------
double dist_green_step(GREEN_TSTP &G1, GREEN_TSTP &G2){
  int size=G1.size1();
  cdmatrix G1_mat(size,size),G2_mat(size,size);
  cdmatrix G1_tv(size,size),G2_tv(size,size);
  cdmatrix G1_les(size,size),G2_les(size,size);
  cdmatrix G1_ret(size,size),G2_ret(size,size);
  double dist=0.0;

  if(G1.tstp() == -1) {
    for(int m=0; m<=G1.ntau(); m++){
      G1.get_mat(m,G1_mat);
      G2.get_mat(m,G2_mat);
      dist += (G1_mat-G2_mat).norm();
    }
  } else {
    for(int m=0; m<=G1.ntau(); m++){
      G1.get_tv(m,G1_tv);
      G2.get_tv(m,G2_tv);
      dist += (G1_tv-G2_tv).norm();
    }
    for(int j=0; j<=G1.tstp(); j++){
      G1.get_les_t_tstp(j,G1_les);
      G2.get_les_t_tstp(j,G2_les);
      dist += (G1_les-G2_les).norm();
      G1.get_ret_tstp_t(j,G1_ret);
      G2.get_ret_tstp_t(j,G2_ret);
      dist += (G1_ret-G2_ret).norm();
    }
  }
  return dist;
}
//--------------------------------------------------------------------------------

TEST_CASE("utilities","[utilities]"){
  const int nt=20, ntau=100, kt=5;
  const double beta=5.0, h=0.05, eps=1.0e-6;
  std::complex<double> I(0.0,1.0);
  int tstp;

  SECTION("extrapolate timestep"){
    GREEN G1(nt,ntau,1,-1);
    GREEN G2(nt,ntau,1,-1);
    GREEN G3(nt,ntau,1,-1);
    cdmatrix ham1(1,1);
    double err=0.0;

    ham1(0,0) = 1.0;
    cntr::green_from_H(G1,0.0,ham1,beta,h);
    cntr::green_from_H(G2,0.0,ham1,beta,h);
    cntr::green_from_H(G3,0.0,ham1,beta,h);

    for(tstp=kt; tstp<nt; tstp++){
      G2.set_timestep(tstp,G1);
      G3.set_timestep(tstp,G1);
      cntr::extrapolate_timestep(tstp, G2, integration::I<double>(kt));
      extrapolate_tv(tstp, kt, h, G3);
      extrapolate_les(tstp, kt, h, G3);
      extrapolate_ret(tstp, kt, h, G3);
      //extrapolate_diag(tstp, kt, h, G3);

      cdmatrix G2_tv,G3_tv;
      G2.get_tv(tstp+1,ntau,G2_tv);
      G3.get_tv(tstp+1,ntau,G3_tv);

      err += cntr::distance_norm2_les(tstp+1,G2,G3);

    }

    REQUIRE(err<eps);
  }

  SECTION("Extrapolate function"){
    const int size=2;
    const double w1=-1.0,w2=0.8,w3=0.7;
    CFUNC func1(nt,size), func2(nt,size), func3(nt,size);
    cdmatrix f(size,size);
    cdmatrix fref(size,size);
    double err=0.0;

    for(tstp=0;tstp<=nt;tstp++){
      f(0,0) = cos(w1*tstp*h);
      f(1,0) = I*0.2*sin(w3*tstp*h);
      f(0,1) = -I*0.2*sin(w3*tstp*h);
      f(1,1) = cos(w2*tstp*h);
      func1.set_value(tstp,f);
      func2.set_value(tstp,f);
      func3.set_value(tstp,f);
    }

    for(tstp=kt; tstp<nt; tstp++){
      cntr::extrapolate_timestep(tstp, func1, integration::I<double>(kt));
      func1.get_value(tstp+1,f);

      extrapolate_function(tstp, kt, h, func3);
      func3.get_value(tstp+1,fref);

      err += (f - fref).norm();
      //std::cout<< "err(" << tstp << ") = " << err << std::endl;

      func2.get_value(tstp+1,f);
      func1.set_value(tstp+1,f);
      func3.set_value(tstp+1,f);
    }

    REQUIRE(err<eps);

  }

  SECTION("Interpolate function"){
    const int size=2;
    const double w1=-1.0,w2=0.8,w3=0.7;
    double tinter;
    CFUNC func1(nt,size);
    cdmatrix f(size,size);
    cdmatrix fref(size,size);
    double err=0.0;

    for(tstp=0;tstp<=nt;tstp++){
      f(0,0) = cos(w1*tstp*h);
      f(1,0) = I*0.2*sin(w3*tstp*h);
      f(0,1) = -I*0.2*sin(w3*tstp*h);
      f(1,1) = cos(w2*tstp*h);
      func1.set_value(tstp,f);
    }

    for(tstp=kt;tstp<nt;tstp++){
      tinter = tstp + 0.5;
      f = cntr::interpolation(tstp+1, tinter, func1, integration::I<double>(kt));
      fref = interpolate_function(tstp+1, kt, h, h*tinter, func1);

      err += (f - fref).norm();
      //std::cout<< "err(" << tstp << ") = " << err << std::endl;
    }

    REQUIRE(err<eps);

  }

  SECTION("set_t0_from_mat"){
    const int size=2;
    const double w1=-1.0,w2=0.8,w3=0.7;
    GREEN G1(nt,ntau,size,FERMION);
    GREEN G2(nt,ntau,size,FERMION);
    cdmatrix ham1(size,size);
    double err=0.0;

    ham1(0,0) = w1;
    ham1(1,0) = I*w3;
    ham1(0,1) = -I*w3;
    ham1(1,1) = w2;

    cntr::green_from_H(G1, 0.0, ham1, beta, h);
    G2.set_timestep(-1,G1);

    set_t0_from_mat(G2);

    err = cntr::distance_norm2(0,G1,G2);
    //std::cout << "err = " << err << std::endl;
    REQUIRE(err<eps);
  }

  SECTION("set_tk_from_mat"){
    const int size=2;
    const double w1=-1.0,w2=0.8,w3=0.7;
    GREEN G1(nt,ntau,size,FERMION);
    GREEN G2(nt,ntau,size,FERMION);
    cdmatrix ham1(size,size);
    double err=0.0;

    ham1(0,0) = w1;
    ham1(1,0) = I*w3;
    ham1(0,1) = -I*w3;
    ham1(1,1) = w2;

    cntr::green_from_H(G1, 0.0, ham1, beta, h);
    G2.set_timestep(-1,G1);

    set_tk_from_mat(G2, kt);

    cdmatrix Gles1(size,size),Gret1(size,size),Gtv1(size,size);
    cdmatrix Gles2(size,size),Gret2(size,size),Gtv2(size,size);
    for(tstp=0; tstp<=kt; tstp++){
      for(int m=0; m<=ntau; m++){
        G2.get_tv(tstp,m,Gtv2);
        G1.get_tv(0,m,Gtv1);
        err += (Gtv1-Gtv2).norm();
      }
      G2.get_les(tstp,tstp,Gles2);
      G1.get_les(tstp,tstp,Gles1);
      err += (Gles1-Gles2).norm();
      G2.get_ret(tstp,tstp,Gret2);
      G1.get_ret(tstp,tstp,Gret1);
      err += (Gret1-Gret2).norm();

      //std::cout << "err = " << err << std::endl;
    }

    REQUIRE(err<eps);
  }

  SECTION("distance_norm2 (herm_matrix)"){
    const int size=2;
    const double w1=-1.0,w2=0.8,w3=0.7;
    const double nu1=-2.0,nu2=1.6,nu3=1.4;
    GREEN G1(nt,ntau,size,FERMION);
    GREEN G2(nt,ntau,size,FERMION);
    cdmatrix ham1(size,size);
    cdmatrix ham2(size,size);
    double dist1,dist2;
    double err=0.0;

    ham1(0,0) = w1;
    ham1(1,0) = I*w3;
    ham1(0,1) = -I*w3;
    ham1(1,1) = w2;

    ham2(0,0) = nu1;
    ham2(1,0) = I*nu3;
    ham2(0,1) = -I*nu3;
    ham2(1,1) = nu2;

    cntr::green_from_H(G1, 0.0, ham1, beta, h);
    cntr::green_from_H(G2, 0.0, ham2, beta, h);

    for(tstp=-1; tstp<=nt; tstp++){
      dist1 = cntr::distance_norm2(tstp, G1, G2);
      dist2 = dist_green(tstp, G1, G2);
      //std::cout << dist1 << " | " << dist2 << std::endl;
      err += fabs(dist1-dist2);
    }
    REQUIRE(err<eps);


  }

  SECTION("distance_norm2 (herm_matrix_timestep)"){
    const int size=2;
    const double w1=-1.0,w2=0.8,w3=0.7;
    const double nu1=-2.0,nu2=1.6,nu3=1.4;
    GREEN G1(nt,ntau,size,FERMION);
    GREEN G2(nt,ntau,size,FERMION);
    GREEN_TSTP G1_step(nt,ntau,size,FERMION);
    GREEN_TSTP G2_step(nt,ntau,size,FERMION);
    cdmatrix ham1(size,size);
    cdmatrix ham2(size,size);
    double dist1,dist2;
    double err=0.0;

    ham1(0,0) = w1;
    ham1(1,0) = I*w3;
    ham1(0,1) = -I*w3;
    ham1(1,1) = w2;

    ham2(0,0) = nu1;
    ham2(1,0) = I*nu3;
    ham2(0,1) = -I*nu3;
    ham2(1,1) = nu2;

    cntr::green_from_H(G1, 0.0, ham1, beta, h);
    cntr::green_from_H(G2, 0.0, ham2, beta, h);

    for(tstp=-1; tstp<=nt; tstp++){
      G1.get_timestep(tstp,G1_step);
      G2.get_timestep(tstp,G2_step);
      dist1 = cntr::distance_norm2(G1_step, G2_step);
      dist2 = dist_green_step(G1_step, G2_step);
      //std::cout << dist1 << " | " << dist2 << std::endl;
      err += fabs(dist1-dist2);
    }
    REQUIRE(err<eps);
  }

  SECTION("force_matsubara_hermitian"){
    const int size=2;
    const double w1=-1.0,w2=0.8,w3=0.7;
    const std::complex<double> wz=(1.0,-0.3);
    GREEN G1(-1,ntau,size,FERMION);
    GREEN G2(-1,ntau,size,FERMION);
    cdmatrix ham1(size,size),ham2(size,size);
    cdmatrix G1_mat(size,size),G2_mat(size,size);
    double err=0.0;

    ham1(0,0) = w1;
    ham1(1,0) = I*w3;
    ham1(0,1) = -I*w3;
    ham1(1,1) = w2;

    ham2(0,0) = w2;
    ham2(1,0) = w3;
    ham2(0,1) = w3;
    ham2(1,1) = w1;

    cntr::green_from_H(G1, 0.0, ham1, beta, h);
    cntr::green_from_H(G2, 0.0, ham2, beta, h);

    G1.incr_timestep(-1,G2,I);

    cntr::force_matsubara_hermitian(G1);

    for(int m=0; m<=ntau; m++){
      G1.get_mat(m,G1_mat);
      err += (G1_mat - G1_mat.adjoint()).norm();
      //      std::cout << m << "  " << G1_mat(0,1) << " " << G1_mat(1,0)  << std::endl;
    }
    REQUIRE(err<eps);

  }

}
