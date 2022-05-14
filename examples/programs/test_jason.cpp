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
  int size=2;
  double h;
  double err_dyson,err_vie2;
  std::complex<double> I(0.0,1.0);

  //============================================================================
  //                          (II) READ INPUT
  //============================================================================
  
  if(argc<2) throw("COMMAND LINE ARGUMENT MISSING");

  flin=argv[1];

  // solver parameters
  find_param(flin,"__Nt=",Nt);
  find_param(flin,"__Ntau=",Ntau);
  // find_param(flin,"__h=",h);
  find_param(flin,"__SolveOrder=",SolveOrder);
  
  // free 2x2 Hamiltonian
  cdmatrix eps_2x2(size,size);
  eps_2x2(0,0) = eps1;
  eps_2x2(1,1) = eps2;
  eps_2x2(0,1) = lam;
  eps_2x2(1,0) = lam;

  CFUNC eps_2x2_func(Nt,size);
  eps_2x2_func.set_constant(eps_2x2);

  // embedding self-energy
  GREEN Sigma(Nt,Ntau,size,FERMION);
  GREEN G(Nt,Ntau,size,FERMION);
  cntr::green_from_H(Sigma, mu, eps_2x2, beta, h);
    
  //=============================================
  //  Solve in integro-differential form (Dyson)
  //=============================================
    
  // equilibrium
  cntr::dyson_mat(G, mu, eps_2x2_func, Sigma, beta, SolveOrder);

  // // start
  // cntr::dyson_start(G, mu, eps_2x2_func, Sigma, beta, h, SolveOrder);

  // // time stepping
  // for (tstp=SolveOrder+1; tstp<=Nt; tstp++) {
  //   cntr::dyson_timestep(tstp, G, mu, eps_2x2_func, Sigma, beta, h, SolveOrder);
  // }

  G.print_to_file("jason.txt",10);
  return 0;
}
