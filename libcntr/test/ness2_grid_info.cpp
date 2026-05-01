#include "catch.hpp"
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include "ness2.hpp"

using namespace ness2;

/*///////////////////////////////////////////////////////////////////////////////////////

  This is a minimal test program to test the grid info helper class in the steady state.

///////////////////////////////////////////////////////////////////////////////////////*/

//------------------------------------------------------------------------------

TEST_CASE("Steady-state FFT grid info","[Steady-state grid info]"){
    const int Nt = 1000;
    const double h1=0.1;
    const double h2=0.01;
    const double dw2=0.1;
    const double tol=1.0e-10;
    double terr = 0.0;
    double werr = 0.0;
    double t1,t2,w1,w2;

    // check initialization of grid info object
    grid_info grid(Nt, h1);
    REQUIRE((std::abs(grid.h_-h1))<tol);
    REQUIRE((std::abs(grid.Nft_-Nt))<tol);
    REQUIRE((std::abs(grid.dw_- M_PI / ((Nt/2) * h1)))<tol);
    
    // check the time_at, freq_at method
    for (int i = 0; i < Nt; ++i) {
        t1 = grid.time_at(i);
        t2 = (i<(Nt/2)) ? i*grid.h_ : (i-Nt)*grid.h_;
        w1 = grid.freq_at(i);
        w2 = (i<(Nt/2)) ? i*grid.dw_ : (i-Nt)*grid.dw_;
        terr += std::abs(t2-t1);
        werr += std::abs(w2-w1);
    }
    REQUIRE(terr/Nt<tol);
    REQUIRE(werr/Nt<tol);
    
    // check the time_grid, freq_grid method
    terr = 0.0;
    werr = 0.0;
    std::vector<double> tgrid = grid.time_grid();
    std::vector<double> wgrid = grid.freq_grid();
    for (int i = 0; i < Nt; ++i) {
        t1 = tgrid[i];
        t2 = (i<(Nt/2)) ? i*grid.h_ : (i-Nt)*grid.h_;
        w1 = wgrid[i];
        w2 = (i<(Nt/2)) ? i*grid.dw_ : (i-Nt)*grid.dw_;
        terr += std::abs(t2-t1);
        werr += std::abs(w2-w1);
    }
    REQUIRE(terr/Nt<tol);
    REQUIRE(werr/Nt<tol);
    
    
    // check the set_h method
    grid.set_h(h2);
    terr = 0.0;
    werr = 0.0;
    for (int i = 0; i < Nt; ++i) {
        t1 = grid.time_at(i);
        t2 = (i<(Nt/2)) ? i*grid.h_ : (i-Nt)*grid.h_;
        w1 = grid.freq_at(i);
        w2 = (i<(Nt/2)) ? i*grid.dw_ : (i-Nt)*grid.dw_;
        terr += std::abs(t2-t1);
        werr += std::abs(w2-w1);
    }
    REQUIRE(terr/Nt<tol);
    REQUIRE(werr/Nt<tol);
    
    // check the set_dw method
    grid.set_dw(dw2);
    terr = 0.0;
    werr = 0.0;
    for (int i = 0; i < Nt; ++i) {
        t1 = grid.time_at(i);
        t2 = (i<(Nt/2)) ? i*grid.h_ : (i-Nt)*grid.h_;
        w1 = grid.freq_at(i);
        w2 = (i<(Nt/2)) ? i*grid.dw_ : (i-Nt)*grid.dw_;
        terr += std::abs(t2-t1);
        werr += std::abs(w2-w1);
    }
    REQUIRE(terr/Nt<tol);
    REQUIRE(werr/Nt<tol);
}
