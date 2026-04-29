#include "catch.hpp"
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

#include "cntr.hpp"

using namespace std;
using namespace ness2;

/*///////////////////////////////////////////////////////////////////////////////////////

  Test the resampling of steady-state Green's functions.

///////////////////////////////////////////////////////////////////////////////////////*/

//------------------------------------------------------------------------------

TEST_CASE("Steady-state resampling","[Steady-state resampling]"){
    double h1= 0.01;
    double h2 = 0.1;
    int Nt1 = 1000;
    int Nt2 = 100;
    int size = 1;
    double tol = 1e-6;
    double beta = 10.0;
    double mu = 0.0;
    double err=0.0;

    herm_matrix_ness G1(Nt1,1);
    herm_matrix_ness G2(Nt2,1);
    herm_matrix_ness Gdo(Nt1,1);
    herm_matrix_ness Gup(Nt2,1);
    
    gauss dos(0.0, 2.0);

    green_equilibrium_ness(-1, G1, dos, beta, mu, h1, FFT_TRAPEZ);
    green_equilibrium_ness(-1, G2, dos, beta, mu, h2, FFT_TRAPEZ);
    
    // Downsampling
    Gdo = downsample(G1,10);
    REQUIRE(distance_norm2(Gdo,G2,fft_domain::time)/Nt2<tol);
    
    // Upsampling
    Gup = upsample(h2,G2,10);
    REQUIRE(distance_norm2(Gup,G1,fft_domain::time)/Nt1<tol);
    
    // Rwsampling
    Gdo = resample(h1,G1,10);
    REQUIRE(distance_norm2(Gdo,G2,fft_domain::time)/Nt2<tol);
    Gup = resample(h2,G2,-10);
    REQUIRE(distance_norm2(Gup,G1,fft_domain::time)/Nt1<tol);
}
