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

  This is a minimal test program to test the solution of the Matsubara-Dyson equation.

///////////////////////////////////////////////////////////////////////////////////////*/

//------------------------------------------------------------------------------
void get_1x1(GREEN &G,GREEN &G11){
  complex<double> g11;
  cdmatrix g;

  for(int n=0; n<=G.ntau(); n++){
    G.get_mat(n,g);
    g11 = g(0,0);
    G11.set_mat(n,g11);
  }
}
//------------------------------------------------------------------------------

TEST_CASE("Matsubara Dyson equation","[Matsubara Dyson equation]"){
  const int fermion = -1;
  const int Nst=2;
  const double eps1 = -1.0;
  const double eps2 = 1.0;
  const double lam = 0.2;
  const double mu = 0.0;
  const double beta = 20.0;
  const int Ntau = 400;
  const int SolverOrder=5;
  const double tol_coarse=1.0e-3;
  const double tol_tight=1.0e-7;

  cdmatrix h2x2(Nst,Nst), h1x1(1,1);
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

  // exact Green's function for 2x2 Hamiltonian
  G2x2 = GREEN(-1,Ntau,Nst,fermion);
  cntr::green_from_H(G2x2,mu,h2x2,beta,1.0);
  // exact Green's function for embedded 1x1 problem
  G_exact = GREEN(-1,Ntau,1,fermion);
  get_1x1(G2x2,G_exact);
  // embedding self-energy
  Sigma = GREEN(-1,Ntau,1,-1);
  cdmatrix h22(1,1);
  h22(0,0) = eps2;
  cntr::green_from_H(Sigma,mu,h22,beta,1.0);
  Sigma.smul(-1,lam*lam);

  G_approx = GREEN(-1,Ntau,1,-1);

  SECTION("Dyson"){
    CFUNC H(-1,1);
    H.set_value(-1,h1x1);
    cntr::dyson_mat(G_approx, Sigma, mu, H, integration::I<double>(SolverOrder), beta, CNTR_MAT_FOURIER);
    err = cntr::distance_norm2(-1,G_exact,G_approx);
    // cout << "Error [Dyson] : " << err << endl;
    REQUIRE(err<tol_coarse);
  }

  SECTION("Fixpoint"){
    CFUNC H(-1,1);
    H.set_value(-1,h1x1);
    cntr::dyson_mat(G_approx, Sigma, mu, H, integration::I<double>(SolverOrder),
		    beta, CNTR_MAT_FIXPOINT);
    err = cntr::distance_norm2(-1,G_exact,G_approx);
    // cout << "Error [Fixpoint] : " << err << endl;
    REQUIRE(err<tol_tight);
  }

  SECTION("Conjugate gradient"){
    CFUNC H(-1,1);
    H.set_value(-1,h1x1);
    cntr::dyson_mat(G_approx, Sigma, mu, H, integration::I<double>(SolverOrder),
  		    beta, CNTR_MAT_CG);
    err = cntr::distance_norm2(-1,G_exact,G_approx);
    //cout << "Error [Steepest descent] : " << err << endl;
    REQUIRE(err<tol_tight);
  }

}
