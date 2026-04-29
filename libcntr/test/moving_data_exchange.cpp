#include "catch.hpp"
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

#include "cntr.hpp"

using namespace std;

/*///////////////////////////////////////////////////////////////////////////////////////

  Testing data exchange between _moving and full Greenfunctions

///////////////////////////////////////////////////////////////////////////////////////*/

//------------------------------------------------------------------------------

TEST_CASE("moving Greens functions: data exchange with herm_matrix","[moving Greens functions: data exchange with herm_matrix]"){
  const int Nst=2;
  const double eps1 = -1.0;
  const double eps2 = 1.0;
  const double lam = 0.2;
  const double mu = -0.1;
  const double beta = 1.0;
  const int ntau = 100;
  const int nt = 200;
  const double dt=0.02;
  const int SolveOrder=5;
  const double tol_coarse=1.0e-3;
  const double tol_tight=1.0e-6;
  int tstp = 0;
  cdmatrix h2x2(Nst,Nst);
  double err;
  int tc=100;

  GREEN G(nt,ntau,Nst,FERMION);
  GREEN Gcc(nt,ntau,Nst,FERMION);
  CFUNC J1(nt,Nst),J2(nt,Nst);
  GTRUNC G_t(tc,Nst,FERMION);
  GTRUNC_TSTP tG_t(tc,Nst,FERMION);
  std::complex<double> I(0.0,1.0);
  // fill some Greensfunctions 
  h2x2(0,0) = eps1;
  h2x2(1,1) = eps2;
  h2x2(0,1) = I*lam;
  h2x2(1,0) = -I*lam;
  cntr::green_from_H(G,mu,h2x2,beta,dt);
  Gcc=G;
  for(int t1=0;t1<=nt;t1++){
    cdmatrix tmp(Nst,Nst);
    tmp(0,0)=cos(t1*0.1);
    tmp(0,1)=-sin(t1*0.1);
    tmp(1,0)=sin(t1*0.1);
    tmp(1,1)=cos(t1*0.2);
    J1.set_value(t1,tmp);
    tmp(0,0)=cos(t1*0.1);
    tmp(0,1)=I*sin(t1*0.1);
    tmp(1,0)=I*sin(t1*0.1);
    tmp(1,1)=cos(t1*0.2);
    J2.set_value(t1,tmp);
  }
  // G(t,t') ->  J1(t)*G(t,t')*J2(t')
  // Gcc(t,t') ->  J2(t)^cc * G(t,t') * J1(t')^cc
  for(int t1=0;t1<=nt;t1++){
    G.left_multiply(t1,J1,1.0);
    G.right_multiply(t1,J2,1.0);
    Gcc.left_multiply_hermconj(t1,J2,1.0);
    Gcc.right_multiply_hermconj(t1,J1,1.0);
  }
  SECTION("reading one timestep"){
    // test is done at three timesteps tstp=tc-10,tc,tc+10
    // also we test reading to timestep variables and Greensfunctions
    err=0;
    for(tstp=0;tstp<=nt;tstp++){
      int smax=(tstp>=tc ? tc : tstp);
      // reading timestep tstp of G  to timestep i of G_t or to teimstep 0 of tG_t
      int i=5;
      G_t.set_timestep(i,G,Gcc,tstp);
      tG_t.set_timestep(G,Gcc,tstp);
      // do the comparison
      for(int t2=0;t2<=tc;t2++){
	cdmatrix m1(2,2),m2(2,2),m3(2,2),dm;
	G_t.get_les(i,t2,m1);
	tG_t.get_les(t2,m3);
	if(tstp-t2>=0){
	  get_les(tstp,tstp-t2,m2,G,Gcc);
	}else{
	  m2.setZero();
	}
	dm=m1-m2;
	err += dm.norm();
	dm=m3-m2;
	err += dm.norm();

	G_t.get_ret(i,t2,m1);
	tG_t.get_ret(t2,m3);
	if(tstp-t2>=0){
	  get_ret(tstp,tstp-t2,m2,G,Gcc);
	}else{
	  m2.setZero();
	}
	dm=m1-m2;
	err += dm.norm();
	dm=m3-m2;
	err += dm.norm();
      }
    }
    //cout << "err=" << err/nt << endl;
    REQUIRE(err/nt<tol_tight);
  }  
}

