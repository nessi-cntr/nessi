#include "catch.hpp"
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

#include "cntr.hpp"

#define GREEN cntr::herm_matrix<double>
#define CPLX std::complex<double>
#define CFUNC cntr::function<double>
using namespace std;

/*///////////////////////////////////////////////////////////////////////////////////////

  This is a minimal test program to test the solution of the Kadanoff-Baym equations.

///////////////////////////////////////////////////////////////////////////////////////*/

//------------------------------------------------------------------------------

TEST_CASE("Dyson, new interface","[Dyson, new interface]"){
  const int fermion = -1;
  const int Nst=2;
  const double eps1 = -1.0;
  const double eps2 = 1.0;
  const double lam = 0.2;
  const double mu = 0.0;
  const double beta = 20.0;
  const int SolveOrder = 5;
  const int Ntau = 400;
  const int Nt = 200;
  const double dt=0.05;
  const double tol_coarse=1.0e-3;
  const double tol_tight=1.0e-6;
  int tstp;
  cdmatrix h2x2(Nst,Nst), h1x1(1,1);
  CFUNC hfunc;
  GREEN G_approx;
  GREEN G2x2,G_exact,Sigma;
  double err;
  std::complex<double> I(0.0,1.0);

  // free 2x2 Hamiltonian
  h2x2(0,0) = eps1;
  h2x2(1,1) = eps2;
  h2x2(0,1) = I*lam;
  h2x2(1,0) = -I*lam;

  // free 1x1 Hamiltonian
  h1x1(0,0) = eps1;
  hfunc = CFUNC(Nt,1);
  hfunc.set_constant(h1x1);

  // exact Green's function for 2x2 Hamiltonian
  G2x2 = GREEN(Nt,Ntau,Nst,fermion);
  cntr::green_from_H(G2x2,mu,h2x2,beta,dt);
  // exact Green's function for embedded 1x1 problem
  G_exact = GREEN(Nt,Ntau,1,fermion);
  for(tstp=-1; tstp<=Nt; tstp++){
    G_exact.set_matrixelement(tstp, 0, 0, G2x2, 0, 0);
  }

  // embedding self-energy
  Sigma = GREEN(Nt,Ntau,1,fermion);
  cdmatrix h22(1,1);
  h22(0,0) = eps2;
  cntr::green_from_H(Sigma,mu,h22,beta,dt);
  for(tstp=-1; tstp<=Nt; tstp++){
    Sigma.smul(tstp,lam*lam);
  }

  G_approx = GREEN(Nt,Ntau,1,fermion);

  // solve equilibrium
  cntr::dyson_mat(G_approx, mu, hfunc, Sigma, beta, SolveOrder, CNTR_MAT_CG);

  SECTION("Integro-differential form"){
    cntr::dyson_start(G_approx,mu,hfunc,Sigma,beta,dt,SolveOrder);

    for(tstp=SolveOrder+1;tstp<=Nt;tstp++){
      cntr::dyson_timestep(tstp,G_approx,mu,hfunc,Sigma,beta,dt,SolveOrder);
    }

    err=0.0;
    for(tstp=0; tstp<=Nt; tstp++){
      err += cntr::distance_norm2(tstp,G_exact,G_approx);
    }
    //cout << "Error [Dyson] : " << err << endl;
    REQUIRE(err/Nt<tol_coarse);
  }

  SECTION("Integral form"){
    GREEN G0(Nt,Ntau,1,fermion);
    GREEN G0xSGM(Nt,Ntau,1,fermion);
    GREEN SGMxG0(Nt,Ntau,1,fermion);

    cntr::green_from_H(G0,mu,h1x1,beta,dt);

    for(tstp=0; tstp<=SolveOrder;tstp++){
      cntr::convolution_timestep(tstp,G0xSGM,G0,Sigma,beta,dt);
      cntr::convolution_timestep(tstp,SGMxG0,Sigma,G0,beta,dt);
      G0xSGM.smul(tstp,-1);
      SGMxG0.smul(tstp,-1);
    }

    cntr::vie2_start(G_approx,G0xSGM,SGMxG0,G0,beta,dt);

    for(tstp=SolveOrder+1;tstp<=Nt;tstp++){
      cntr::convolution_timestep(tstp,G0xSGM,G0,Sigma,beta,dt);
      cntr::convolution_timestep(tstp,SGMxG0,Sigma,G0,beta,dt);
      G0xSGM.smul(tstp,-1);
      SGMxG0.smul(tstp,-1);
      cntr::vie2_timestep(tstp,G_approx,G0xSGM,SGMxG0,G0,beta,dt);
    }

    err=0.0;
    for(tstp=0; tstp<=Nt; tstp++){
      err += cntr::distance_norm2(tstp,G_exact,G_approx);
    }
    //cout << "Error [VIE2] : " << err << endl;
    REQUIRE(err/Nt<tol_tight);
  }

#if CNTR_USE_OMP==1
  SECTION("Integral form (omp)"){
    GREEN G0(Nt,Ntau,1,fermion);
    GREEN G0xSGM(Nt,Ntau,1,fermion);
    GREEN SGMxG0(Nt,Ntau,1,fermion);
    int nthreads;

    #pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }

    cntr::green_from_H(G0,mu,h1x1,beta,dt);

    for(tstp=0;tstp<=Nt;tstp++){
      cntr::convolution_timestep_omp(nthreads,tstp,G0xSGM,G0,Sigma,beta,dt,SolveOrder);
      cntr::convolution_timestep_omp(nthreads,tstp,SGMxG0,Sigma,G0,beta,dt,SolveOrder);
      G0xSGM.smul(tstp,-1);
      SGMxG0.smul(tstp,-1);
      cntr::vie2_timestep_omp(nthreads,tstp,G_approx,G0xSGM,SGMxG0,G0,beta,dt,SolveOrder);
    }

    err=0.0;
    for(tstp=0; tstp<=Nt; tstp++){
      err += cntr::distance_norm2(tstp,G_exact,G_approx);
    }
    //cout << "Error [VIE2 (omp)] : " << err << endl;
    REQUIRE(err/Nt<tol_tight);
  }
#endif // CNTR_USE_OMP==1
}
