#ifndef HUBB_H
#define HUBB_H

#include <sys/stat.h>
#include <complex>
#include "cntr/cntr.hpp"

using namespace cntr;

namespace hubb {
  void Polarization(int tstp, GREEN &G, GREEN &Pol);
  void Sigma_MF(int tstp, GREEN &G, CFUNC &U, cdmatrix &Sigma_mf);
  void Ham_MF(int tstp, GREEN &G, CFUNC &U, cdmatrix &h0, CFUNC &hmf);
  void Sigma_2B(int tstp, GREEN &G, CFUNC &U, GREEN &Sigma);
  void Sigma_GW(int tstp, GREEN &G, CFUNC &U, GREEN &Chi, GREEN &Sigma);
  void Sigma_TPP(int tstp, GREEN &G, CFUNC &U, GREEN &TPP, GREEN &Sigma);
 
  void GenChi(double h, double beta, GREEN &Pol, CFUNC &U,
	    GREEN &PxU, GREEN &UxP, GREEN &Chi, int SolveOrder=MAX_SOLVE_ORDER);
  void GenChi(int tstp, double h, double beta, GREEN &Pol, CFUNC &U,
	    GREEN &PxU, GREEN &UxP, GREEN &Chi, int SolveOrder=MAX_SOLVE_ORDER);


  void GenTPP(double h, double beta, GREEN &G, GREEN &Phi, CFUNC &U,
	      GREEN &UxPhi, GREEN &PhixU, GREEN &TPP, int SolveOrder=MAX_SOLVE_ORDER);
  void GenTPP(int tstp, double h, double beta, GREEN &G, GREEN &Phi, CFUNC &U,
	      GREEN &UxPhi, GREEN &PhixU, GREEN &TPP, int SolveOrder=MAX_SOLVE_ORDER);
}

#endif
