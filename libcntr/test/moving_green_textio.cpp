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

  Testing data 

///////////////////////////////////////////////////////////////////////////////////////*/

//------------------------------------------------------------------------------

TEST_CASE("moving Greens functions: text-file I/O","[moving Greens functions: text-file I/O]"){
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

  GREEN G(nt,ntau,Nst,FERMION);
  CFUNC J1(nt,Nst);
  GTRUNC G_t(tc,Nst,FERMION);
  GTRUNC G_t_1;
  std::complex<double> I(0.0,1.0);
  // fill some Greensfunctions 
  h2x2(0,0) = eps1;
  h2x2(1,1) = eps2;
  h2x2(0,1) = I*lam;
  h2x2(1,0) = -I*lam;
  cntr::green_from_H(G,mu,h2x2,beta,dt);
  for(int t1=0;t1<=nt;t1++){
    cdmatrix tmp(Nst,Nst);
    tmp(0,0)=cos(t1*0.1);
    tmp(0,1)=-sin(t1*0.1);
    tmp(1,0)=sin(t1*0.1);
    tmp(1,1)=cos(t1*0.2);
    J1.set_value(t1,tmp);
  }
  // G(t,t') ->  J1(t)*G(t,t')*J2(t')
  // Gcc(t,t') ->  J2(t)^cc * G(t,t') * J1(t')^cc
  for(int t1=0;t1<=nt;t1++){
    G.left_multiply(tstp,J1);
    G.right_multiply_hermconj(tstp,J1);
  }
    G_t.set_from_G_backward(G,G,tc);
  SECTION("hdf5"){
    // test is done at three timesteps tstp=tc-10,tc,tc+10
    // also we test reading to timestep variables and Greensfunctions
    err=0;
    std::string filename =  "test_Cyiwek_data.txt";
    G_t.print_to_file(filename.c_str(),16); // file now has group a/ with a/ret/, a/les/, ...
    G_t_1.read_from_file(filename.c_str()); // file now has group a/ with a/ret/, a/les/, ...
    err=0;
    for(int t1=0;t1<=tc;t1++){
      err += cntr::distance_norm2(G_t,t1,G_t_1,t1);
    }
    std::remove(filename.c_str());
    //cout << err << endl;		   
    REQUIRE(err/nt<tol_tight);
  }  
}


  

  
