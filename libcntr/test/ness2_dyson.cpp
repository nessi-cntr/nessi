#include "catch.hpp"
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

#include "cntr.hpp"
#include "ness2.hpp"

using namespace std;
using namespace ness2;

/*///////////////////////////////////////////////////////////////////////////////////////

  This is a minimal test program to test the solution of the Kadanoff-Baym equations in the steady state.

///////////////////////////////////////////////////////////////////////////////////////*/

//------------------------------------------------------------------------------

TEST_CASE("Steady-state Kadanoff-Baym equations (Dyson)","[Steady-state Kadanoff-Baym equations]"){
    double h= 0.01;
    int Nft = 1024;
    int ntau = 10;
    int size = 1;
    int SolveOrder=5;
    double tol = 1e-3;
    double beta = 10.0;
    double mu = 0.0;
    cdmatrix ham(1,1);
    double eps = 0.0;
    ham(0,0) = eps;

    gauss dos(0.0, 2.0);

    // === Set up sigma in real time and ness ===
    cntr::herm_matrix<double> S_cntr(Nft,ntau,size,FERMION);
    ness2::herm_matrix_ness S_ness(Nft, size);
    cntr::green_equilibrium(S_cntr, dos, beta, h, mu, 100, 20);
    S_ness.read_from_cntr(S_cntr,Nft-1);

    // steady-state solution
    ness2::herm_matrix_ness G_ness(Nft,1);
    ness2::dyson(G_ness, mu, ham, S_ness, h, ness2::FFT_TRAPEZ);
    
    // solve the full dyson equation:
    cntr::function<double> hfunc(Nft,size);
    cntr::herm_matrix<double> G_cntr(Nft,ntau,size,FERMION);
    for(int t=0;t<=Nft;t++){
        hfunc.set_value(t,ham);
    }
    cntr::dyson(G_cntr,mu,hfunc,S_cntr,beta,h,SolveOrder,CNTR_MAT_FIXPOINT);
    ness2::herm_matrix_ness G_ness2(Nft,1);
    ness2::cntr2ness(G_ness2,G_cntr,Nft-1);
      
    double err=0.0;
    err += ness2::distance_norm2(G_ness,G_ness2,ness2::fft_domain::time);
    REQUIRE(err/Nft<tol);
}
