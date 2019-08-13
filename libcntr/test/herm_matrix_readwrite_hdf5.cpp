#include <stdio.h>
#include <string.h>
#include "catch.hpp"
#include "cntr.hpp"
#include <cmath>


using namespace std;
#define GREEN cntr::herm_matrix<double>
#define CFUNC cntr::function<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>

#if CNTR_USE_HDF5 == 1

TEST_CASE("read/write (hdf5)","[herm_matrix_read_write_hdf5]"){
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

  SECTION("read/write (group)"){
    double err=0.0;
    G1.write_to_hdf5("herm_matrix_hdf5.h5","testgroup");
    G2.read_from_hdf5("herm_matrix_hdf5.h5","testgroup");
    for(int tstp=-1; tstp<=nt; tstp++){
      err += cntr::distance_norm2(tstp,G1,G2);
    }
    REQUIRE(err<eps);
  }
  
  SECTION("read/write up to given time step"){
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

#endif
