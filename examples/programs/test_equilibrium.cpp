#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

// contour library headers
#include "cntr/cntr.hpp"
#include "cntr/utils/read_inputfile.hpp"

// local headers to include
#include "formats.hpp"

using namespace std;

//==============================================================================
//         main program
//==============================================================================
int main(int argc,char *argv[]){
  const double eps1 = -1.0;
  const double eps2 = 1.0;
  const double lam = 0.5;
  const double mu = 0.0;
  const double beta = 20.0;
  const double h = 1.0;
  //..................................................
  //                input
  //..................................................
  int Ntau,SolveOrder;
  char* flin;
  //..................................................
  //                internal
  //..................................................
  double err_fourier,err_fixpoint;
  //.................................................
  try{
    //============================================================================
    //                          (II) READ INPUT
    //============================================================================
    {
      if(argc<2) throw("COMMAND LINE ARGUMENT MISSING");

      // solver parameters
      find_param(argv[1],"__Ntau=",Ntau);
      find_param(argv[1],"__SolveOrder=",SolveOrder);
    }

    {
      // free 2x2 Hamiltonian
      cdmatrix eps_2x2(2,2);
      eps_2x2(0,0) = eps1;
      eps_2x2(1,1) = eps2;
      eps_2x2(0,1) = II*lam;
      eps_2x2(1,0) = -II*lam;

      // free 1x1 Hamiltonian
      CFUNC eps_11_func(-1,1);
      eps_11_func.set_constant(eps1*MatrixXcd::Identity(1,1));

      // exact Green's function for 2x2 Hamiltonian
      GREEN G2x2(-1,Ntau,2,FERMION);
      cntr::green_from_H(G2x2,mu,eps_2x2,beta,h);

      // exact Green's function for embedded 1x1 problem
      GREEN G_exact(-1,Ntau,1,FERMION);
      G_exact.set_matrixelement(-1,0,0,G2x2,0,0);

      // embedding self-energy
      GREEN Sigma(-1,Ntau,1,FERMION);
      cdmatrix eps_22=eps2*MatrixXcd::Identity(1,1);
      cntr::green_from_H(Sigma, mu, eps_22, beta, h);
      Sigma.smul(-1,lam*lam);

      // approximate Greens function
      GREEN G_approx(-1,Ntau,1,FERMION);

      // solve using Fourier method
      cntr::dyson_mat(G_approx, mu, eps_11_func, Sigma, beta, SolveOrder, CNTR_MAT_FOURIER);
      err_fourier = cntr::distance_norm2(-1,G_exact,G_approx) / Ntau;

      // solve using fix-point iteration
      cntr::dyson_mat(G_approx, mu, eps_11_func, Sigma, beta, SolveOrder, CNTR_MAT_FIXPOINT);
      err_fixpoint = cntr::distance_norm2(-1,G_exact,G_approx) / Ntau;

    }

    {
      FILE *pErrDysonEq = fopen("out/test_equilibrium.dat","a");
      fprintf(pErrDysonEq,"%.5g ", err_fourier);
      fprintf(pErrDysonEq,"%.5g ", err_fixpoint);
      fprintf(pErrDysonEq,"\n");

      fclose(pErrDysonEq);
    }

  } // try
  catch(char *message){
    cerr << "exception\n**** " << message << " ****" << endl;
    cerr << " No input file found. Exiting ... " << endl;
  }
  catch(...){
    cerr << " No input file found. Exiting ... " << endl;
  }
  return 0;
}
