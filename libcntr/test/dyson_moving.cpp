#include "catch.hpp"
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

#include "cntr.hpp"

using namespace std;

/*///////////////////////////////////////////////////////////////////////////////////////

  This is a minimal test program to test the solution of the trunctated Kadanoff-Baym 
  equations, by comparing them to the full KBE in a case where the self-energy is 
  actually short range.

  note that the two only agree up to the typical time propagatioon error

///////////////////////////////////////////////////////////////////////////////////////*/

//------------------------------------------------------------------------------

TEST_CASE("truncated Kadanoff-Baym equations (Dyson)","[truncated Kadanoff-Baym equations]"){
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
  int tstp;
  cdmatrix h2x2(Nst,Nst);
  double err;
  int tc=100;

  CFUNC hfunc(nt,Nst);
  GREEN G(nt,ntau,Nst,FERMION),Sigma(nt,ntau,Nst,FERMION);
  CTRUNC hfunc_t(tc,Nst);
  GTRUNC G_t(tc,Nst,FERMION),Sigma_t(tc,Nst,FERMION);
  std::complex<double> I(0.0,1.0);

  // for sigma, take a single orbital G, and multiply with 
  // free 2x2 Hamiltonian
  h2x2(0,0) = eps1;
  h2x2(1,1) = eps2;
  h2x2(0,1) = I*lam;
  h2x2(1,0) = -I*lam;
  cntr::green_from_H(Sigma,mu,h2x2,beta,dt);
  // multiply Sigma with some time-dependent functions exp(-(5*|t-t'|/tc))**2)
  Sigma.set_timestep_zero(-1); // matsubara set to zero
  for(int t1=0;t1<=nt;t1++){
    cdmatrix tmp(Nst,Nst);
    for(int t2=0;t2<=t1;t2++){
      double fac=(t2<=tc ? exp(-(25.0*t2*t2)/(tc*tc)) : 0.0);
      Sigma.get_ret(t1,t1-t2,tmp);
      tmp *= fac*cos(0.2*dt*t1)*cos(0.2*dt*(t1-t2)); // break ant time transl. symmetry
      Sigma.set_ret(t1,t1-t2,tmp);
      Sigma.get_les(t1-t2,t1,tmp);
      tmp *= fac*cos(0.2*dt*(t1-t2))*cos(0.2*dt*t1); 
      Sigma.set_les(t1-t2,t1,tmp);
    }
    // set tv part to zero
    tmp.setZero();
    for (int tau=0;tau<=ntau;tau++){
      Sigma.set_tv(t1,tau,tmp);
    }    
  }
  // solve the full dyson equastion:
  hfunc = CFUNC(nt,Nst);
  hfunc.set_value(-1,h2x2); // use same H
  for(int t1=0;t1<=nt;t1++){
     cdmatrix tmp;
     tmp=h2x2*cos((0.1*dt*t1));
     hfunc.get_value(t1,tmp);
  }
  cntr::dyson(G,mu,hfunc,Sigma,beta,dt,SolveOrder,CNTR_MAT_FIXPOINT);

  // truncated time evolution:

  // init moving G's from G's:
  G_t.set_from_G_backward(G,G,tc);
  Sigma_t.set_from_G_backward(Sigma,Sigma,tc);
  hfunc_t.set_from_function_backward(hfunc,tc);         
  // truncatedtime-stepping
  err=0.0;
  for(tstp=tc+1;tstp<=nt;tstp++){
    cdmatrix tmp;
    GTRUNC_TSTP gtmp(tc,Nst,FERMION);
    Sigma_t.forward();
    // read in data from Sigma at leading timestep      
    Sigma_t.set_timestep(0,Sigma,Sigma,tstp);
    // set matrix element of hfunc:
    hfunc_t.forward();
    hfunc.get_value(tstp,tmp);
    hfunc_t.set_value(0,tmp);
    G_t.forward();
    cntr::dyson_timestep(G_t,Sigma_t,hfunc_t,mu,SolveOrder,dt);
    // compare:
    gtmp.set_timestep(G,G,tstp);
    double err1=cntr::distance_norm2(G_t,0,gtmp);
    //cout << "tstp= " << tstp << " err= " << err1 << endl;
    err+=err1;
  }
  //cout << "|Dyson trunc - Dyson | / nt = " << err/nt  << endl;
  REQUIRE(err/nt<tol_coarse);
}

