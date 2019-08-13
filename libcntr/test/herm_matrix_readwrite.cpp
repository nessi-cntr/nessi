#include "catch.hpp"
#include "cntr.hpp"
#include <cmath>

using namespace std;
#define GREEN cntr::herm_matrix<double>
#define CFUNC cntr::function<double>

TEST_CASE("read/write (plain text)","[herm_matrix_read_write]"){
    
    

  int prec=10;
  int size=2;
  int nt=100, ntau=50;
  double eps=1e-6;
  double h=0.01, mu=0.0, beta=10.0;
  double tmax=h*nt;
  double eps1=-0.4,eps2=0.6,lam=0.1;
  std::complex<double> I(0.0,1.0);
  cdmatrix h0(2,2);
  GREEN G1(nt,ntau,size,-1);
  GREEN G2(nt,ntau,size,-1);
  
  h0(0,0) = eps1;
  h0(1,1) = eps2;
  h0(0,1) = I*lam;
  h0(1,0) = -I*lam;

  cntr::green_from_H(G1,mu,h0,beta,h,5,4,true);

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
