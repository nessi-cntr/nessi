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
//    contains
//==============================================================================

//------------------------------------------------------------------------------
void GenKernel(int tstp, GREEN &G0, GREEN &Sigma, GREEN &F, GREEN &Fcc,
  const double beta, const double h, const int SolveOrder){

    cntr::convolution_timestep(tstp, F, G0, Sigma, beta, h, SolveOrder);
    cntr::convolution_timestep(tstp, Fcc, Sigma, G0, beta, h, SolveOrder);
    F.smul(tstp,-1);
    Fcc.smul(tstp,-1);
  }
//------------------------------------------------------------------------------

//==============================================================================
//         main program
//==============================================================================
int main(int argc,char *argv[]){
  const double eps1 = -1.0;
  const double eps2 = 1.0;
  const double lam = 0.5;
  const double mu = 0.0;
  const double beta = 20.0;
  //..................................................
  //                input
  //..................................................
  int Nt,Ntau,SolveOrder;
  double Tmax;
  char* flin;
  //..................................................
  //                internal
  //..................................................
  int tstp;
  double h;
  double err_dyson,err_vie2;
  CFUNC eps_11_func;
  GREEN Sigma;
  GREEN G_exact;
  std::complex<double> I(0.0,1.0);

  try{
    //============================================================================
    //                          (II) READ INPUT
    //============================================================================
    {
      if(argc<2) throw("COMMAND LINE ARGUMENT MISSING");

      flin=argv[1];

      // solver parameters
      find_param(flin,"__Nt=",Nt);
      find_param(flin,"__Ntau=",Ntau);
      find_param(flin,"__Tmax=",Tmax);
      find_param(flin,"__SolveOrder=",SolveOrder);

      h = Tmax/Nt;
    }

    {
      // free 2x2 Hamiltonian
      cdmatrix eps_2x2(2,2);
      eps_2x2(0,0) = eps1;
      eps_2x2(1,1) = eps2;
      eps_2x2(0,1) = I*lam;
      eps_2x2(1,0) = -I*lam;

      // free 1x1 Hamiltonian
      eps_11_func = CFUNC(Nt,1);
      eps_11_func.set_constant(eps1*MatrixXcd::Identity(1,1));

      // exact Green's function for 2x2 Hamiltonian
      GREEN G2x2(Nt,Ntau,2,FERMION);
      cntr::green_from_H(G2x2,mu,eps_2x2,beta,h);

      // exact Green's function for embedded 1x1 problem
      G_exact = GREEN(Nt,Ntau,1,FERMION);
      for(tstp=-1; tstp<=Nt; tstp++) {
        G_exact.set_matrixelement(tstp,0,0,G2x2,0,0);
      }

      // embedding self-energy
      Sigma = GREEN(Nt,Ntau,1,FERMION);
      cdmatrix eps_22=eps2*MatrixXcd::Identity(1,1);
      cntr::green_from_H(Sigma, mu, eps_22, beta, h);
      for(tstp=-1; tstp<=Nt; tstp++) {
        Sigma.smul(tstp,lam*lam);
      }

    }
    //=============================================
    //  Solve in integro-differential form (Dyson)
    //=============================================
    {
      GREEN G_approx(Nt,Ntau,1,FERMION);

      // equilibrium
      cntr::dyson_mat(G_approx, mu, eps_11_func, Sigma, beta, SolveOrder);

      // start
      cntr::dyson_start(G_approx, mu, eps_11_func, Sigma, beta, h, SolveOrder);

      // time stepping
      for (tstp=SolveOrder+1; tstp<=Nt; tstp++) {
        cntr::dyson_timestep(tstp, G_approx, mu, eps_11_func, Sigma, beta, h, SolveOrder);
      }

      // check error
      err_dyson=0.0;
      for(tstp=0; tstp<=Nt; tstp++){
        err_dyson += 2.0 * cntr::distance_norm2_les(tstp, G_exact, G_approx) / (Nt*Nt);
        err_dyson += 2.0 * cntr::distance_norm2_ret(tstp, G_exact, G_approx) / (Nt*Nt);
        err_dyson += cntr::distance_norm2_tv(tstp, G_exact, G_approx) / (Nt*Ntau);
      }
    }
    //=============================================
    //      Solve in integral form (VIE2)
    //=============================================
    {
      // noninteracting 1x1 Greens function (Sigma=0)
      GREEN G0(Nt,Ntau,1,FERMION);
      cdmatrix eps_11=eps1*MatrixXcd::Identity(1,1);
      cntr::green_from_H(G0, mu, eps_11, beta, h);

      GREEN G_approx(Nt,Ntau,1,FERMION);
      GREEN F(Nt,Ntau,1,FERMION);
      GREEN Fcc(Nt,Ntau,1,FERMION);

      // equilibrium
      GenKernel(-1, G0, Sigma, F, Fcc, beta, h, SolveOrder);
      cntr::vie2_mat(G_approx, F, Fcc, G0, beta, SolveOrder);

      // start
      for(tstp=0; tstp <= SolveOrder; tstp++){
        GenKernel(tstp, G0, Sigma, F, Fcc, beta, h, SolveOrder);
      }
      cntr::vie2_start(G_approx, F, Fcc, G0, beta, h, SolveOrder);

      // time stepping
      for (tstp=SolveOrder+1; tstp<=Nt; tstp++) {
        GenKernel(tstp, G0, Sigma, F, Fcc, beta, h, SolveOrder);
        cntr::vie2_timestep(tstp, G_approx, F, Fcc, G0, beta, h, SolveOrder);
      }

      // check error
      err_vie2=0.0;
      for(tstp=0; tstp<=Nt; tstp++){
        err_vie2 += 2.0 * cntr::distance_norm2_les(tstp, G_exact, G_approx) / (Nt*Nt);
        err_vie2 += 2.0 * cntr::distance_norm2_ret(tstp, G_exact, G_approx) / (Nt*Nt);
        err_vie2 += cntr::distance_norm2_tv(tstp, G_exact, G_approx) / (Nt*Ntau);
      }
    }

    //==============================
    //    save error to file
    //==============================
    {
      FILE *pErrDysonEq = fopen("out/test_nonequilibrium.dat","a");
      fprintf(pErrDysonEq,"%.5g ", err_dyson);
      fprintf(pErrDysonEq,"%.5g ", err_vie2);
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
