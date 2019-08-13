#include <stdio.h>
#include <string.h>
#include "catch.hpp"
#include "cntr.hpp"
#include <cmath>


// TODO:
// Set matrixelement
// Incr
// Right/left multiply
// Interaction with herm_matrix_timestep

// TODO:
// read/write (slices) for hdf5

#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define cfunction cntr::function<double>
#define CPLX std::complex<double>
using namespace std;

#if CNTR_USE_HDF5 == 1
// Test get_set for herm_matrix,herm_matrix_timestep,cfunction
double setget(GREEN &A,cdmatrix &a){
  double toterr=0.0;
  // Matsubara
  for(int m=0;m<A.ntau();m++){
    cdmatrix tmp;
    A.set_mat(m,a);
    A.get_mat(m,tmp);
    toterr+=(tmp-a).norm();
    // tv
    for(int n=0;n<A.nt();n++){
      A.set_tv(n,m,a);
      A.get_tv(n,m,a);
      toterr+=(tmp-a).norm();
    }
  }

  for(int m=0;m<A.nt();m++){
    for(int n=0;n<m;n++){
      cdmatrix ret,les;
      // Ret
      A.set_ret(m,n,a);
      A.get_ret(m,n,ret);
      toterr+=(ret-a).norm();
      // Les
      A.set_les(n,m,a);
      A.get_les(n,m,les);
      toterr+=(les-a).norm();
    }
  }
  return toterr;
}


void exact_rightmultiply(double beta,double dt,GREEN &G){
  int ntau=G.ntau();
  int nt=G.nt();
  int size=G.size1();
  assert(size==2);

  double dtau=beta/ntau,tau;
  std::complex<double> I(0.0,1.0);

  // Matsubara
  cdmatrix mat(size,size),tv(size,size);
  for(int m=0;m<=ntau;m++){
    tau=m*dtau;
    mat(0,0)=std::complex<double>(-1.7071067811865475,-0.17677669529663675)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-0.2928932188134524,+0.17677669529663687)*exp(2.0*tau);
    mat(0,1)=std::complex<double>(-0.4267766952966371,-1.0606601717798207)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-0.0732233047033631,+1.0606601717798212)*exp(2.0*tau);
    mat(1,0)=std::complex<double>(-0.07322330470336319,0.7071067811865475)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-0.4267766952966368,-0.7071067811865475)*exp(2.0*tau);
    mat(1,1)=std::complex<double>(-0.43933982822017864,0.17677669529663687)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-2.560660171779821,-0.17677669529663687)*exp(2.0*tau);
    mat=mat/(1.0+exp(2.0*beta));
    // std::cout << "Mat " << tau << " " << mat << std::endl;

    G.set_mat(m,mat);
    for(int n=0;n<=nt;n++){
      double t1=n*dt;
      tv(0,0)=std::complex<double>(0.17677669529663687,0.2928932188134524)*exp(2.0*I*t1 + 2.0*(beta-tau)) - std::complex<double>(0.1767766952966371,-1.707106781186548)*exp(-2.0*I*t1 + tau*2.0);
      tv(0,1)=std::complex<double>(1.0606601717798212,0.07322330470336319)*exp(2.0*I*t1 + 2.0*(beta-tau)) - std::complex<double>(1.0606601717798216,-0.426776695296637)*exp(-2.0*I*t1 + tau*2.0);
      tv(1,0)=-1.0*std::complex<double>(0.7071067811865475,-0.4267766952966369)*exp(2.0*I*t1 + 2.0*(beta-tau)) + std::complex<double>(0.7071067811865475,0.07322330470336313)*exp(-2.0*I*t1 + tau*2.0);
      tv(1,1)=-1.0*std::complex<double>(0.17677669529663675,-2.5606601717798214)*exp(2.0*I*t1 + 2.0*(beta-tau)) + std::complex<double>(0.1767766952966369,0.4393398282201787)*exp(-2.0*I*t1 + tau*2.0);
      tv=tv/(1.0+exp(2.0*beta));
      G.set_tv(n,m,tv);
    }
  }
  // Les + ret
  cdmatrix les(size,size),ret(size,size);
  for(int m=0;m<=nt;m++){
    for(int n=0;n<=m;n++){
      double t1=m*dt;
      double t2=n*dt;

      // Ret
      ret(0,0)=exp(-2.0*I*(t2+t1))*cos(t2)*(exp(4.0*I*t1)*std::complex<double>(-0.17677669529663687,-0.2928932188134524)+exp(4.0*I*t2)*std::complex<double>(0.1767766952966371,-1.707106781186548));
      ret(0,1)=exp(-2.0*I*(t2+t1))*cos(t2)*(exp(4.0*I*t1)*std::complex<double>(-1.0606601717798212,-0.07322330470336319)+exp(4.0*I*t2)*std::complex<double>(1.0606601717798216,-0.426776695296637));
      ret(1,0)=exp(-2.0*I*(t2+t1))*cos(t2)*(exp(4.0*I*t1)*std::complex<double>(0.7071067811865475,-0.4267766952966369)-exp(4.0*I*t2)*std::complex<double>(0.7071067811865475,0.07322330470336313));
      ret(1,1)=exp(-2.0*I*(t2+t1))*cos(t2)*(exp(4.0*I*t1)*std::complex<double>(0.17677669529663675,-2.5606601717798214)-exp(4.0*I*t2)*std::complex<double>(0.1767766952966369,0.4393398282201787));
      G.set_ret(m,n,ret);

      // Les
      t1=n*dt;
      t2=m*dt;
      les(0,0)=exp(-2.0*I*(t2+t1))*cos(t2)*(std::complex<double>(-0.1767766952966371,1.707106781186548)*exp(4.0*I*t2)+std::complex<double>(0.17677669529663687,0.2928932188134524)*exp(4.0*I*t1+2.0*beta));
      les(0,1)=exp(-2.0*I*(t2+t1))*cos(t2)*(std::complex<double>(-1.0606601717798216,0.426776695296637)*exp(4.0*I*t2)+std::complex<double>(1.0606601717798212,0.07322330470336319)*exp(4.0*I*t1+2.0*beta));
      les(1,0)=cos(t2)*(std::complex<double>(0.7071067811865475,0.07322330470336313)*exp(2.0*I*(t2-t1))+std::complex<double>(-0.7071067811865475,0.4267766952966369)*exp(2.0*I*(t1-t2)+2.0*beta));
      les(1,1)=cos(t2)*(std::complex<double>(0.1767766952966369,0.4393398282201787)*exp(2.0*I*(t2-t1))+std::complex<double>(-0.17677669529663675,2.5606601717798214)*exp(2.0*I*(t1-t2)+2.0*beta));
      les=les/(1.0+exp(2.0*beta));
      G.set_les(n,m,les);
    }
  }
}


void exact_leftmultiply(double beta,double dt,GREEN &G){
  int ntau=G.ntau();
  int nt=G.nt();
  int size=G.size1();
  assert(size==2);

  double dtau=beta/ntau,tau;
  std::complex<double> I(0.0,1.0);

  // Matsubara
  cdmatrix mat(size,size),tv(size,size);
  for(int m=0;m<=ntau;m++){
    tau=m*dtau;
    mat(0,0)=std::complex<double>(-1.7071067811865475,0.17677669529663675)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-0.2928932188134524,-0.17677669529663687)*exp(2.0*tau);
    mat(0,1)=std::complex<double>(-0.07322330470336319,-0.7071067811865475)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-0.4267766952966368,0.7071067811865475)*exp(2.0*tau);
    mat(1,0)=std::complex<double>(-0.4267766952966371,1.0606601717798207)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-0.0732233047033631,-1.0606601717798212)*exp(2.0*tau);
    mat(1,1)=std::complex<double>(-0.43933982822017864,-0.17677669529663687)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-2.560660171779821,0.17677669529663687)*exp(2.0*tau);
    mat=mat/(1.0+exp(2.0*beta));
    // std::cout << "Mat " << tau << " " << mat << std::endl;

    G.set_mat(m,mat);
    for(int n=0;n<=nt;n++){
      double t1=n*dt;
      tv(0,0)=cos(t1)*(std::complex<double>(-0.17677669529663687,0.2928932188134524)*exp(2.0*I*t1 + 2.0*(beta-tau)) + std::complex<double>(0.1767766952966371,1.707106781186548)*exp(-2.0*I*t1 + tau*2.0));
      tv(0,1)=cos(t1)*(std::complex<double>(0.7071067811865475,0.4267766952966369)*exp(2.0*I*t1 + 2.0*(beta-tau)) + std::complex<double>(-0.7071067811865475,0.07322330470336313)*exp(-2.0*I*t1 + tau*2.0));
      tv(1,0)=cos(t1)*(std::complex<double>(-1.0606601717798212,0.07322330470336319)*exp(2.0*I*t1 + 2.0*(beta-tau)) + std::complex<double>(1.0606601717798216,0.426776695296637)*exp(-2.0*I*t1 + tau*2.0));
      tv(1,1)=cos(t1)*(std::complex<double>(0.17677669529663675,2.5606601717798214)*exp(2.0*I*t1 + 2.0*(beta-tau)) + std::complex<double>(-0.1767766952966369,0.4393398282201787)*exp(-2.0*I*t1 + tau*2.0));
      tv=tv/(1.0+exp(2.0*beta));
      G.set_tv(n,m,tv);
    }
  }
  // Les + ret
  cdmatrix les(size,size),ret(size,size);
  for(int m=0;m<=nt;m++){
    for(int n=0;n<=m;n++){
      double t1=m*dt;
      double t2=n*dt;

      // Ret
      ret(0,0)=exp(-2.0*I*(t2+t1))*cos(t1)*(exp(4.0*I*t1)*std::complex<double>(0.17677669529663687,-0.2928932188134524)-exp(4.0*I*t2)*std::complex<double>(0.1767766952966371,1.707106781186548));
      ret(0,1)=exp(-2.0*I*(t2+t1))*cos(t1)*(exp(4.0*I*t1)*std::complex<double>(-0.7071067811865475,-0.4267766952966369)+exp(4.0*I*t2)*std::complex<double>(0.7071067811865475,-0.07322330470336313));
      ret(1,0)=exp(-2.0*I*(t2+t1))*cos(t1)*(exp(4.0*I*t1)*std::complex<double>(1.0606601717798212,-0.07322330470336319)-exp(4.0*I*t2)*std::complex<double>(1.0606601717798216,0.426776695296637));
      ret(1,1)=exp(-2.0*I*(t2+t1))*cos(t1)*(exp(4.0*I*t1)*std::complex<double>(-0.17677669529663675,-2.5606601717798214)+exp(4.0*I*t2)*std::complex<double>(0.1767766952966369,-0.4393398282201787));
      G.set_ret(m,n,ret);

      // Les
      t1=n*dt;
      t2=m*dt;
      les(0,0)=cos(t1)*(std::complex<double>(0.1767766952966371,1.707106781186548)*exp(2.0*I*(t2-t1))-std::complex<double>(0.17677669529663687,-0.2928932188134524)*exp(2.0*I*(t1-t2)+2.0*beta));
      les(0,1)=cos(t1)*exp(-2.0*I*(t2+t1))*(std::complex<double>(-0.7071067811865475,0.07322330470336325)*exp(4.0*I*t2)+std::complex<double>(0.7071067811865476,0.4267766952966368)*exp(4.0*I*t1+2.0*beta));
      les(1,0)=cos(t1)*(std::complex<double>(1.0606601717798219,0.4267766952966371)*exp(2.0*I*(t2-t1))+std::complex<double>(-1.0606601717798212,0.07322330470336325)*exp(2.0*I*(t1-t2)+2.0*beta));
      les(1,1)=cos(t1)*exp(-2.0*I*(t2+t1))*(std::complex<double>(-0.1767766952966369,0.4393398282201787)*exp(4.0*I*t2)+std::complex<double>(0.17677669529663675,2.5606601717798214)*exp(4.0*I*t1+2.0*beta));
      les=les/(1.0+exp(2.0*beta));
      G.set_les(n,m,les);
    }
  }
}



TEST_CASE("Herm matrix","[Herm_matrix]"){

  std::complex<double> I(0.0,1.0);

  int nt=100;
  int ntau=50;
  int size=2;
  int size2=5;
  double eps=1e-7;
  double beta=5.0;
  double h=0.01;
  double tmax=h*nt;

  double eps1=-0.4,eps2=0.6,lam1=0.1;
  double mu=0.0;

  cdmatrix h0(2,2);
  cdmatrix h1(2,2);
  cdmatrix h2(2,2);
  GREEN G1(nt,ntau,size,-1);
  GREEN G2(nt,ntau,size,-1);
  GREEN G3(nt,ntau,size,-1);
  GREEN G4(nt,ntau,size,-1);


  {
    GREEN g=GREEN(nt,ntau,size,-1);
    //REQUIRE(g.nt()==nt);

    // Square Green's function
    GREEN A=GREEN(nt,ntau,size,-1);
    cdmatrix a(size,size);
    a.setZero();
    a(0,0)=sqrt(2.0);
    a(0,1)=sqrt(2.0)*std::complex<double>(0.0,1.0);
    a(1,0)=sqrt(2.0)*std::complex<double>(0.0,-1.0);
    a(1,1)=-sqrt(2.0);
    cntr::green_from_H(A,0.0,a,beta,h);

    // Rectangular Green's function
    GREEN B=GREEN(nt,ntau,size,size2,-1);
    cdmatrix b(size,size2);
    b.setZero();

    for (int i=0;i<size;i++){
      for (int j=0;j<size2;j++){
        b(i,j)=i+j*2;
      }
    }

    // Function
    cfunction funcC(nt,size);
    cdmatrix c(size,size);
    c.setZero();
    double t;
    for (int tstp=-1;tstp<nt;tstp++){

      if(tstp==-1){
        t=0.0;
      }else{
        t=tstp*h;
      }
      c(0,0)=2.0*cos(t);
      c(0,1)=0.5*cos(t);
      c(1,0)=0.5*cos(t);
      c(1,1)=3.0*cos(t);
      funcC.set_value(tstp,c);
    }

    SECTION ("Set/get for herm_matrix"){
      // Writing/reading from normal green's function
      REQUIRE(setget(A,a)<eps);
      // Writing/reading from rectangular Green's function
      REQUIRE(setget(B,b)<eps);
    }

    SECTION ("Left/right multiply"){
      GREEN Ar,Al;
      //Left/Right multiply
      Ar=A;
      Al=A;
      for(int tstp=-1;tstp<=nt;tstp++){
        Ar.right_multiply(tstp,funcC);
        Al.left_multiply(tstp,funcC);
      }

      double err=0.0;
      GREEN exactR(nt,ntau,size,-1);
      exact_rightmultiply(beta,h,exactR);
      for(int tstp=-1;tstp<nt;tstp++){
        err+=cntr::distance_norm2(tstp,Ar,exactR);
      }
      REQUIRE(err<eps);

      err=0.0;
      GREEN exactL(nt,ntau,size,-1);
      exact_leftmultiply(beta,h,exactL);
      for(int tstp=-1;tstp<nt;tstp++){
        err+=cntr::distance_norm2(tstp,Al,exactL);
      }
      REQUIRE(err<eps);

      // REQUIRE(simple_readwrite(A,a)<eps);
      // // Writing/reading from rectangular Green's function
      // REQUIRE(simple_readwrite(B,b)<eps);
    }

  }
  ///////////////////////////////////
  //Herm Matrix Algebra
  ////////////////////////////////////
  {
    double eps3=0.435,eps4=0.5676,lam2=0.1566;
    double wr=0.3;
    std::complex<double> wz;

    wz = 1.0 -0.3*I;

    h1(0,0) = eps1;
    h1(1,1) = eps2;
    h1(0,1) = I*lam1;
    h1(1,0) = -I*lam1;

    h2(0,0) = eps3;
    h2(1,1) = eps4;
    h2(0,1) = I*lam2;
    h2(1,0) = -I*lam2;

    cntr::green_from_H(G1,mu,h1,beta,h);
    cntr::green_from_H(G2,mu,h2,beta,h);

    SECTION("herm_matrix algebra: incr_timestep"){
      double err=0.0;
      cdmatrix mat1(2,2),tv1(2,2),les1(2,2),ret1(2,2);
      cdmatrix mat2(2,2),tv2(2,2),les2(2,2),ret2(2,2);
      cdmatrix mat3(2,2),tv3(2,2),les3(2,2),ret3(2,2);

      for(int q=0; q<=ntau; q++){
        G1.get_mat(q,mat1);
        G2.get_mat(q,mat2);
        mat3 = mat1 + wz*mat2;
        G3.set_mat(q,mat3);
      }

      for(int i=0; i<=nt; i++){
        for(int j=0; j<=i; j++){
          G1.get_les(j,i,les1);
          G2.get_les(j,i,les2);
          les3 = les1 + wz*les2;
          G3.set_les(j,i,les3);

          G1.get_ret(i,j,ret1);
          G2.get_ret(i,j,ret2);
          ret3 = ret1 + wz*ret2;
          G3.set_ret(i,j,ret3);
        }
        for(int q=0; q<=ntau; q++){
          G1.get_tv(i,q,tv1);
          G2.get_tv(i,q,tv2);
          tv3 = tv1 + wz*tv2;
          G3.set_tv(i,q,tv3);
        }
      }

      cntr::green_from_H(G4,mu,h1,beta,h);
      err = 0.0;
      G4.incr_timestep(G2,wz);
      for(int tstp=-1; tstp<=nt; tstp++){
        err += cntr::distance_norm2(tstp,G3,G4);
      }
      REQUIRE(err<eps);

      cntr::green_from_H(G4,mu,h1,beta,h);
      err = 0.0;
      for(int tstp=-1; tstp<=nt; tstp++){
        GREEN_TSTP A(tstp,ntau,size);
        G2.get_timestep(tstp,A);
        G4.incr_timestep(tstp,A,wz);
        err += cntr::distance_norm2(tstp,G3,G4);
      }
      REQUIRE(err<eps);

      cntr::green_from_H(G4,mu,h1,beta,h);
      err = 0.0;
      for(int tstp=-1; tstp<=nt; tstp++){
        G4.incr_timestep(tstp,G2,wz);
        err += cntr::distance_norm2(tstp,G3,G4);
      }
      REQUIRE(err<eps);

    }

    SECTION("herm_matrix algebra: smul (complex weight)"){
      double err=0.0;
      cdmatrix mat1(2,2),tv1(2,2),les1(2,2),ret1(2,2);
      cdmatrix mat3(2,2),tv3(2,2),les3(2,2),ret3(2,2);

      for(int q=0; q<=ntau; q++){
        G1.get_mat(q,mat1);
        mat3 = wz*mat1;
        G3.set_mat(q,mat3);
      }

      for(int i=0; i<=nt; i++){
        for(int j=0; j<=i; j++){
          G1.get_les(j,i,les1);
          les3 = wz*les1;
          G3.set_les(j,i,les3);

          G1.get_ret(i,j,ret1);
          ret3 = wz*ret1;
          G3.set_ret(i,j,ret3);
        }
        for(int q=0; q<=ntau; q++){
          G1.get_tv(i,q,tv1);
          tv3 = wz*tv1;
          G3.set_tv(i,q,tv3);
        }
      }

      err = 0.0;
      cntr::green_from_H(G4,mu,h1,beta,h);
      for(int tstp=-1; tstp<=nt; tstp++){
        G4.smul(tstp,wz);
        err += cntr::distance_norm2(tstp,G3,G4);
      }
      REQUIRE(err<eps);
    }

    SECTION("herm_matrix algebra: smul (real weight)"){
      double err=0.0;
      cdmatrix mat1(2,2),tv1(2,2),les1(2,2),ret1(2,2);
      cdmatrix mat3(2,2),tv3(2,2),les3(2,2),ret3(2,2);

      for(int q=0; q<=ntau; q++){
        G1.get_mat(q,mat1);
        mat3 = wr*mat1;
        G3.set_mat(q,mat3);
      }

      for(int i=0; i<=nt; i++){
        for(int j=0; j<=i; j++){
          G1.get_les(j,i,les1);
          les3 = wr*les1;
          G3.set_les(j,i,les3);

          G1.get_ret(i,j,ret1);
          ret3 = wr*ret1;
          G3.set_ret(i,j,ret3);
        }
        for(int q=0; q<=ntau; q++){
          G1.get_tv(i,q,tv1);
          tv3 = wr*tv1;
          G3.set_tv(i,q,tv3);
        }
      }

      err = 0.0;
      cntr::green_from_H(G4,mu,h1,beta,h);
      for(int tstp=-1; tstp<=nt; tstp++){
        G4.smul(tstp,wr);
        err += cntr::distance_norm2(tstp,G3,G4);
      }
      REQUIRE(err<eps);
    }

  }
  ///////////////////////////////////
  //Herm Matrix set/get timestep
  ////////////////////////////////////
  SECTION("set_get_timestep"){

    h0(0,0) = eps1;
    h0(1,1) = eps2;
    h0(0,1) = I*lam1;
    h0(1,0) = -I*lam1;

    cntr::green_from_H(G1,mu,h0,beta,h);


    double err=0.0;
    for(int tstp=-1; tstp<=nt; tstp++){
      GREEN_TSTP A(tstp,ntau,size);
      G1.get_timestep(tstp,A);
      G2.set_timestep(tstp,A);
      err += cntr::distance_norm2(tstp,G1,G2);
    }
    REQUIRE(err<eps);

  }

  ///////////////////////////////////
  //Herm Matrix readwrite
  ////////////////////////////////////

  {

    int prec=12;

    h0(0,0) = eps1;
    h0(1,1) = eps2;
    h0(0,1) = I*lam1;
    h0(1,0) = -I*lam1;

    GREEN G5(nt,ntau,size,-1);
    GREEN G6(nt,ntau,size,-1);

    cntr::green_from_H(G1,mu,h0,beta,h);

    SECTION("read/write"){
      double err=0.0;
      G1.print_to_file("herm_matrix_plain.dat",prec);
      G2.read_from_file("herm_matrix_plain.dat");
      for(int tstp=-1; tstp<=nt; tstp++){
        err += cntr::distance_norm2(tstp,G1,G2);
      }
      REQUIRE(err<eps);
    }

    SECTION("read/write up to given time step"){
      int nt1=50;
      double err=0.0;
      G1.print_to_file("herm_matrix_plain.dat",prec);
      G2.read_from_file(nt1, "herm_matrix_plain.dat");
      for(int tstp=-1; tstp<=nt1; tstp++){
        err += cntr::distance_norm2(tstp,G1,G2);
      }
      REQUIRE(err<eps);
    }

  }


  ///////////////////////////////////
  //Herm Matrix hdf5 for bethe
  ////////////////////////////////////
  /*
  {
  int size_tmp = 1;

  GREEN G_tmp(nt,ntau,size_tmp,-1);
  GREEN G1_tmp(nt,ntau,size_tmp,-1);

  GREEN_TSTP Gtstp;
  cntr::herm_matrix_timestep_view<double> Gview;
  cfunction f, f1;

  //cdmatrix a_tmp(size_tmp,size_tmp);
  //a_tmp(0,0)=eps1;
  //cntr::green_from_H(G_tmp,0.0,a_tmp,beta,h,5,4,true);

  SECTION("hdf5 for bethe :write/read"){
  // Read entire Green's function from file
  cntr::green_equilibrium_bethe(G_tmp, 1.0, 0.1);
  G_tmp.write_to_hdf5("gbethe.h5", "G");
  G1_tmp.read_from_hdf5("gbethe.h5", "G");

  {
  double err = 0.;
  for (int tstp = -1; tstp <= nt; tstp++)
  err += cntr::distance_norm2(tstp, G_tmp, G1_tmp);
  REQUIRE(err < eps);
}
}

SECTION("hdf5 for bethe :write/read upto read_nt"){
// Read from t=0 to read_nt from file
cntr::green_equilibrium_bethe(G_tmp, 1.0, 0.1);
G1_tmp.clear();
int read_nt = 2;
G1_tmp.read_from_hdf5(read_nt, "gbethe.h5", "G");

{
double err = 0.;
for (int tstp = -1; tstp <= read_nt; tstp++)
err += cntr::distance_norm2(tstp, G_tmp, G1_tmp);
REQUIRE(err < eps);
}
}

SECTION("hdf5 for bethe :write/read for a single time step"){
// Write and read a single time step
cntr::green_equilibrium_bethe(G_tmp, 1.0, 0.1);
int select_timestep = 4;
Gview = cntr::herm_matrix_timestep_view<double>(select_timestep, G_tmp);
Gview.write_to_hdf5("gbethe_tstp4.h5", "G");
Gtstp.read_from_hdf5("gbethe_tstp4.h5", "G");
G1_tmp.set_timestep(Gtstp.tstp_, Gtstp);

{
double err = cntr::distance_norm2(select_timestep, G_tmp, G1_tmp);
REQUIRE(err < eps);
}
}

SECTION("hdf5 for bethe:write/read for a certain time slice"){
// Write individual time slices out to disk
// and read only time slice no 8
cntr::green_equilibrium_bethe(G_tmp, 1.0, 0.1);
G_tmp.write_to_hdf5_slices("g_slices.h5", "G", 2);

hid_t file_id = read_hdf5_file("g_slices.h5");
hid_t group_id = open_group(file_id, "G");
hid_t sub_group_id = open_group(group_id, "t8");
Gview = cntr::herm_matrix_timestep_view<double>(8, G1_tmp);
Gview.read_from_hdf5(sub_group_id);
close_group(sub_group_id);
close_group(group_id);
close_hdf5_file(file_id);

{
double err = cntr::distance_norm2(8, G_tmp, G1_tmp);
REQUIRE(err < eps);
}
}


{
// Write and read some contour functions to hdf5 file

f = cntr::function<double>(nt, size);
for (int tstp = -1; tstp <= nt; tstp++) {
f.ptr(tstp)[0] = 100 * tstp + 0;
f.ptr(tstp)[1] = 100 * tstp + 1;
f.ptr(tstp)[2] = 100 * tstp + 2;
f.ptr(tstp)[3] = 100 * tstp + 3;
}

f.write_to_hdf5("func.h5", "f");
f1.read_from_hdf5("func.h5", "f");

SECTION("hdf5:write/read contour function"){
double err = 0.;
for (int tstp = -1; tstp <= nt; tstp++)
err += std::abs(f.ptr(tstp)[1] - f1.ptr(tstp)[1]);
REQUIRE(err < eps);
}

f1 = cntr::function<double>(nt, size);

int max_timestep = 4;
f1.read_from_hdf5(max_timestep, "func.h5", "f");

SECTION("hdf5:write/read contour function upto max_timestep "){
double err = 0.;
for (int tstp = -1; tstp <= max_timestep; tstp++)
err += std::abs(f.ptr(tstp)[1] - f1.ptr(tstp)[1]);
REQUIRE(err < eps);
}
}

}
*/

///////////////////////////////////
//Herm Matrix hdf5 readwrite
////////////////////////////////////

{

  int prec=10;

  h0(0,0) = eps1;
  h0(1,1) = eps2;
  h0(0,1) = I*lam1;
  h0(1,0) = -I*lam1;


  cntr::green_from_H(G1,mu,h0,beta,h);

  SECTION("read/write (hdf5): group"){
    double err=0.0;
    G1.write_to_hdf5("herm_matrix_hdf5.h5","testgroup");
    G2.read_from_hdf5("herm_matrix_hdf5.h5","testgroup");
    for(int tstp=-1; tstp<=nt; tstp++){
      err += cntr::distance_norm2(tstp,G1,G2);
    }
    REQUIRE(err<eps);
  }

  SECTION("read/write (hdf5): up to given time step"){
    int nt1=50;
    double err=0.0;
    G1.write_to_hdf5("herm_matrix_hdf5.h5","testgroup");
    G2.read_from_hdf5(nt1, "herm_matrix_hdf5.h5","testgroup");
    for(int tstp=-1; tstp<=nt1; tstp++){
      err += cntr::distance_norm2(tstp,G1,G2);
    }
    REQUIRE(err<eps);
  }


  SECTION("read/write (slices)"){
    int output_every=10;
    double err=0.0;

    G1.write_to_hdf5_slices("herm_matrix_hdf5_slices.h5","testgroup",output_every);

    for(int j=0; j<output_every; j++){
      int tstp=j*output_every;
      char subgroup[20];
      char groupname[40];

      std::sprintf(subgroup, "t%d", tstp);
      strcpy(groupname,"testgroup/");
      strcat(groupname,subgroup);

      GREEN_TSTP A(tstp,ntau,size);
      A.read_from_hdf5("herm_matrix_hdf5_slices.h5",groupname);
      G2.set_timestep(tstp,A);
      err += cntr::distance_norm2(tstp,G1,G2);
    }

    REQUIRE(err<eps);
  }
}


///////////////////////////////////
//Herm Matrix submatirx
////////////////////////////////////

{

  GREEN A,B,Asub,B1;

  // Test diagonal
  A=GREEN(nt,ntau,size2,-1);
  Asub=GREEN(nt,ntau,size,-1);
  B=GREEN(nt,ntau,size,-1);
  B1=GREEN(nt,ntau,3,-1);

  cdmatrix a(size2,size2);
  a.setZero();
  cdmatrix b(2,2);

  for (int i=0;i<size;i++){
    for (int j=0;j<size;j++){
      a(i,j)=i+j;
    }
  }
  for (int i=size2-3;i<size2;i++){
    for (int j=size2-3;j<size2;j++){
      a(i,j)=i+j+1;
    }
  }
  cntr::green_from_H(A,0.0,a,beta,h);

  std::vector<int> i1{0,0,1,1};
  std::vector<int> i2{0,1,0,1};
  std::vector<int> j1{0,0,1,1};
  std::vector<int> j2{0,1,0,1};
  Asub.set_submatrix(j1,j2,A,i1,i2);


  cdmatrix b1(3,3);
  b1.setZero();
  std::vector<int> i11{2,2,2,3,3,3,4,4,4};
  std::vector<int> i12{2,3,4,2,3,4,2,3,4};
  std::vector<int> j11{0,0,0,1,1,1,2,2,2};
  std::vector<int> j12{0,1,2,0,1,2,0,1,2};
  for (int i=0;i<i11.size();i++){
    b1(j11[i],j12[i])=a(i11[i],i12[i]);
  }

  cntr::green_from_H(B1,0.0,b1,beta,h);

  SECTION ("Green of submatrix = Submatrix of Green"){
    double err=0;
    // Take elements of first submatrix
    b.setZero();
    for (int i=0;i<i1.size();i++){
      b(j1[i],j2[i])=a(i1[i],i2[i]);
    }
    cntr::green_from_H(B,0.0,b,beta,h);
    for(int tstp=-1;tstp<nt;tstp++){
      err+=cntr::distance_norm2(tstp,B,Asub);
    }
    REQUIRE(err<eps);
  }

  SECTION ("Green of submatrix = Submatrix of Green, Test 2"){
    double err=0;
    GREEN Asub1=GREEN(nt,ntau,3,-1);
    Asub1.set_submatrix(j11,j12,A,i11,i12);

    for(int tstp=-1;tstp<nt;tstp++){
      err+=cntr::distance_norm2(tstp,B1,Asub1);
    }
    REQUIRE(err<eps);
  }

}


///////////////////////////////////
//Herm Matrix get matsubara
////////////////////////////////////
SECTION("get matminus for fermion"){

  h0(0,0)=sqrt(2.0);
  h0(0,1)=sqrt(2.0)*complex<double>(0.0,1.0);
  h0(1,0)=sqrt(2.0)*complex<double>(0.0,-1.0);
  h0(1,1)=-sqrt(2.0);

  double err=0.0;
  cntr::green_from_H(G1,0.0,h0,beta,h);

  for(int itau=0;itau<=ntau;itau++){
    cdmatrix A(size,size);
    cdmatrix B(size,size);

    G1.get_matminus(itau,A);
    G1.get_mat(ntau-itau,B);

    err+=(A+B).norm();
  }

  REQUIRE(err<eps);
}
}

#endif
