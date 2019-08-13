#include "catch.hpp"
#include "cntr.hpp"


using namespace std;
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>


TEST_CASE("set/get timestep","[herm_matrix_set_get_timestep]"){
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

  SECTION("set_get_timestep"){
    double err=0.0;
    for(int tstp=-1; tstp<=nt; tstp++){
      GREEN_TSTP A(tstp,ntau,size);
      G1.get_timestep(tstp,A);
      G2.set_timestep(tstp,A);     
      err += cntr::distance_norm2(tstp,G1,G2);
    }
    REQUIRE(err<eps);
  }
 
}

