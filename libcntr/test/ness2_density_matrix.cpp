#include "catch.hpp"
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

#include "cntr.hpp"
#include "ness2.hpp"

using namespace std;

/*///////////////////////////////////////////////////////////////////////////////////////

  Test the resampling of steady-state Green's functions.

///////////////////////////////////////////////////////////////////////////////////////*/

//------------------------------------------------------------------------------

TEST_CASE("Steady-state density_matrix check with real time","[Steady-state density_matrix]"){
    double h = 0.01;
    int Nt = 1024;
    int Ntau = 100;
    int size = 1;
    double tol = 1e-6;
    double beta = 10.0;
    double mu = 0.0;
    double err=0.0;

    // set up Green's fcts
    cdmatrix Mness(size,size);
    cdmatrix Mcntr(size,size);
    ness2::herm_matrix_ness Gness(Nt,1);
    cntr::herm_matrix<double> Gcntr(Nt,Ntau,size,-1);
    
    ness2::gauss dos(0.0, 2.0);

    ness2::green_equilibrium_ness(-1, Gness, dos, beta, mu, h, ness2::FFT_TRAPEZ);
    cntr::green_equilibrium(Gcntr, dos, beta, h, mu, 100, 20);
    
    // calculate density matrix
    Gcntr.density_matrix(0,Mcntr);
    ness2::density_matrix(Mness, -1, Gness);

    REQUIRE(std::abs(Mcntr(0,0).real()-Mness(0,0).real())<tol);
}

TEST_CASE("Steady-state convolution_density_matrix check with real time","[Steady-state convolution_density_matrix]"){
    double h = 0.01;
    int Nt = 1024;
    int Ntau = 100;
    int size = 1;
    double tol = 1e-3;
    double beta = 10.0;
    double mu = 0.0;
    double err=0.0;
    int SolveOrder = 5;

    // set up Green's fcts
    cdmatrix Mness(size,size);
    cdmatrix Mcntr(size,size);
    ness2::herm_matrix_ness Aness(Nt,1);
    cntr::herm_matrix<double> Acntr(Nt,Ntau,size,-1);
    ness2::herm_matrix_ness Bness(Nt,1);
    cntr::herm_matrix<double> Bcntr(Nt,Ntau,size,-1);
    
    ness2::gauss dosa(0.0, 2.0);
    ness2::gauss dosb(0.0, 4.0);

    ness2::green_equilibrium_ness(-1, Aness, dosa, beta, mu, h, ness2::FFT_TRAPEZ);
    cntr::green_equilibrium(Acntr, dosa, beta, h, mu, 100, 20);
    ness2::green_equilibrium_ness(-1, Bness, dosb, beta, mu, h, ness2::FFT_TRAPEZ);
    cntr::green_equilibrium(Bcntr, dosb, beta, h, mu, 100, 20);
    
    // calculate convolution density matrix with A,B
    cntr::convolution_density_matrix(0,Mcntr,Acntr,Bcntr,beta,h,SolveOrder);
    ness2::convolution_density_matrix(Mness, -1, Aness, Bness, h);

    REQUIRE(std::abs(Mcntr(0,0).real()-Mness(0,0).real())<tol);
}
