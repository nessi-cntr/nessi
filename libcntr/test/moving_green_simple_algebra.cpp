#include "catch.hpp"
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <cstdio>  // For std::remove

#include "cntr.hpp"

using namespace std;

/*///////////////////////////////////////////////////////////////////////////////////////

  Testing data exchange between _moving and full Greenfunctions

///////////////////////////////////////////////////////////////////////////////////////*/

//------------------------------------------------------------------------------

TEST_CASE("moving Greens functions: simple algebra","[moving Greens functions: simple algebra]"){
  const int Nst=2;
  const double eps1 = -1.0;
  const double eps2 = 1.0;
  const double lam = 0.2;
  const double mu = -0.1;
  const double beta = 1.0;
  const int ntau = 100;
  const int nt = 100;
  const double dt=0.02;
  const int SolveOrder=5;
  const double tol_coarse=1.0e-3;
  const double tol_tight=1.0e-6;
  int tstp;
  cdmatrix h2x2(Nst,Nst);
  double err;
  int tc=nt;

  GREEN G1(nt,ntau,Nst,FERMION);
  GREEN G2(nt,ntau,Nst,FERMION);
  CFUNC J1(nt,Nst), J2(nt,Nst);
  std::complex<double> I(0.0,1.0);
  // fill some Greensfunctions 
  h2x2(0,0) = eps1;
  h2x2(1,1) = eps2;
  h2x2(0,1) = I*lam;
  h2x2(1,0) = -I*lam;
  cntr::green_from_H(G1,mu,h2x2,beta,dt);
  cntr::green_from_H(G2,mu,h2x2,beta,dt);
  for(int t1=0;t1<=nt;t1++){
    cdmatrix tmp(Nst,Nst);
    tmp(0,0)=cos(t1*0.1);
    tmp(0,1)=-sin(t1*0.1);
    tmp(1,0)=sin(t1*0.1);
    tmp(1,1)=cos(t1*0.2);
    J1.set_value(t1,tmp);
    tmp(0,0)=cos(t1*0.2);
    tmp(0,1)=I*sin(t1*0.2);
    tmp(1,0)=I*sin(t1*0.2);
    tmp(1,1)=cos(t1*0.3);
    J2.set_value(t1,tmp);
  }
  // G(t,t') ->  J1(t)*G(t,t')*J2(t')
  // Gcc(t,t') ->  J2(t)^cc * G(t,t') * J1(t')^cc
  for(int t1=0;t1<=nt;t1++){
    G1.left_multiply(t1,J1,1.0);
    G1.right_multiply_hermconj(t1,J1,1.0);
    G2.left_multiply(t1,J2,1.0);
    G2.right_multiply_hermconj(t1,J2,1.0);
  }
  GTRUNC G1_t(tc,Nst,FERMION);
  GTRUNC G2_t(tc,Nst,FERMION);
  G1_t.set_from_G_backward(G1,G1,tc);
  G2_t=G1_t;
  SECTION("set matrixelement"){
    GTRUNC_TSTP G1x1_t(tc,1,FERMION);
    int i=5;
    G2_t.set_matrixelement(0,1,0,G1_t,i,0,1);
    err=0.0;
    // do the comparison
    for(int t2=0;t2<=tc;t2++){
      cdmatrix m1,m2,m3,dm;
      G1_t.get_les(0,t2,m1);
      G1_t.get_les(i,t2,m3);
      m1(1,0)=m3(0,1);
      G2_t.get_les(0,t2,m2);
      dm=m1-m2;
      err += dm.norm();
      G1_t.get_ret(0,t2,m1);
      G1_t.get_ret(i,t2,m3);
      m1(1,0)=m3(0,1);
      G2_t.get_ret(0,t2,m2);
      dm=m1-m2;
      err += dm.norm();
    }
    //cout << "err=" << err/nt << endl;
    REQUIRE(err/nt<tol_tight);
  }  
  SECTION("incr and smul"){
    int i=5;
    GTRUNC_TSTP tG_t(G2_t,i); // store  step 0 of G2
    cdouble fac(4,2),fac1(2,2);
    G2_t.incr_timestep(i,G1_t,i,fac1);
    G2_t.smul(i,fac);
    err=0.0;
    // do the comparison
    for(int t2=0;t2<=tc;t2++){
      cdmatrix m1,m2,m3,m4,dm;
      G1_t.get_les(i,t2,m1);
      tG_t.get_les(t2,m2);
      m4=fac*(m2+fac1*m1);
      G2_t.get_les(i,t2,m3);
      dm=m3-m4;
      err += dm.norm();
      G1_t.get_ret(i,t2,m1);
      tG_t.get_ret(t2,m2);
      m4=fac*(m2+fac1*m1);
      G2_t.get_ret(i,t2,m3);
      //cout << m3<<endl;
      dm=m3-m4;
      err += dm.norm();
    }
    //cout << "err=" << err/nt << endl;
    REQUIRE(err/nt<tol_tight);
  }  
  SECTION(" function multiplication"){
    int i=5;
    GTRUNC_TSTP tG_t; // store  step 0 of G2
    CTRUNC J_t(tc,Nst);
    cdmatrix m1,m2,m3,m4,dm;
    
    J_t.set_from_function_backward(J2,tc);
    tG_t=GTRUNC_TSTP(G2_t,0);
    G2_t.left_multiply(J_t); // G to J*G
    err=0.0;
    // do the comparison
    for(int t2=0;t2<=tc;t2++){
      tG_t.get_les(t2,m1);
      J_t.get_value(0,m2);
      m3=m2*m1;
      G2_t.get_les(0,t2,m4);
      dm=m3-m4;
      err += dm.norm();
      tG_t.get_ret(t2,m1);
      J_t.get_value(0,m2);
      m3=m2*m1;
      G2_t.get_ret(0,t2,m4);
      dm=m3-m4;
      err += dm.norm();
    }
    //cout << "err=" << err/nt << endl;
    REQUIRE(err/nt<tol_tight);

    J_t.set_from_function_backward(J2,tc);
    tG_t=GTRUNC_TSTP(G2_t,0);
    G2_t.left_multiply_hermconj(J_t); // G to J*G
    err=0.0;
    // do the comparison
    for(int t2=0;t2<=tc;t2++){
      tG_t.get_les(t2,m1);
      J_t.get_value(0,m2);
      m3=m2.adjoint()*m1;
      G2_t.get_les(0,t2,m4);
      dm=m3-m4;
      err += dm.norm();
      tG_t.get_ret(t2,m1);
      J_t.get_value(0,m2);
      m3=m2.adjoint()*m1;
      G2_t.get_ret(0,t2,m4);
      dm=m3-m4;
      err += dm.norm();
    }
    //cout << "err=" << err/nt << endl;
    REQUIRE(err/nt<tol_tight);

    J_t.set_from_function_backward(J2,tc);
    tG_t=GTRUNC_TSTP(G2_t,0);
    G2_t.right_multiply(J_t); // G to J*G
    err=0.0;
    // do the comparison
    for(int t2=0;t2<=tc;t2++){
      tG_t.get_les(t2,m1);
      J_t.get_value(t2,m2);
      m3=m1*m2;
      G2_t.get_les(0,t2,m4);
      dm=m3-m4;
      err += dm.norm();
      tG_t.get_ret(t2,m1);
      J_t.get_value(t2,m2);
      m3=m1*m2;
      G2_t.get_ret(0,t2,m4);
      dm=m3-m4;
      err += dm.norm();
    }
    //cout << "err=" << err/nt << endl;
    REQUIRE(err/nt<tol_tight);

    J_t.set_from_function_backward(J2,tc);
    tG_t=GTRUNC_TSTP(G2_t,0);
    G2_t.right_multiply_hermconj(J_t); // G to J*G
    err=0.0;
    // do the comparison
    for(int t2=0;t2<=tc;t2++){
      tG_t.get_les(t2,m1);
      J_t.get_value(t2,m2);
      m3=m1*m2.adjoint();
      G2_t.get_les(0,t2,m4);
      dm=m3-m4;
      err += dm.norm();
      tG_t.get_ret(t2,m1);
      J_t.get_value(t2,m2);
      m3=m1*m2.adjoint();
      G2_t.get_ret(0,t2,m4);
      dm=m3-m4;
      err += dm.norm();
    }
    //cout << "err=" << err/nt << endl;
    REQUIRE(err/nt<tol_tight);    
  }

  
}


  

  
