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

  Testing hdf5 read and write for moving Greenfunctions

///////////////////////////////////////////////////////////////////////////////////////*/

//------------------------------------------------------------------------------

TEST_CASE("moving Greens functions: hdf5","[moving Greens functions: hdf5]"){
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
  cdmatrix h2x2(Nst,Nst);
  double err;
  int tc=nt;

  GREEN G(nt,ntau,Nst,FERMION);
  CFUNC J1(nt,Nst);
  GTRUNC G_t(tc,Nst,FERMION);
  GTRUNC G_t_1(tc,Nst,FERMION);
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
    G.left_multiply(t1,J1,1.0);
    G.right_multiply_hermconj(t1,J1,1.0);
  }
  G_t.set_from_G_backward(G,G,tc);
  
std::string filename =  "test_Cyiwek_data.h5";
std::string groupname =  "G";

SECTION("hdf5"){
err=0;
hid_t file_id = open_hdf5_file(filename);
G_t.write_to_hdf5(file_id, "G"); // file now has group a/ with a/ret/, a/les/, ...
close_hdf5_file(file_id);

file_id = read_hdf5_file(filename.c_str());
G_t_1.read_from_hdf5(file_id, "G");
close_hdf5_file(file_id);
for(int t1=0;t1<=tc;t1++){
  err += cntr::distance_norm2(G_t,t1,G_t_1,t1);
}
std::remove(filename.c_str());
//cout << err << endl;
REQUIRE(err/nt<tol_tight);
}
}


  

  
