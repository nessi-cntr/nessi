#ifndef HOLS_H
#define HOLS_H

#include <sys/stat.h>
#include <complex>
#include "cntr/cntr.hpp"
#include <vector>

using namespace cntr;
using namespace std;

#define CFUNC cntr::function<double>
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>

namespace Hols {
    
    void Sigma_Mig(int tstp, GREEN &G, GREEN &Sigma, GREEN &D0, GREEN &D,
                      GREEN &Pi, GREEN &D0_Pi, GREEN &Pi_D0, CFUNC &g_el_ph, double beta,
                      double h, int SolverOrder, int MAT_METHOD=CNTR_MAT_FIXPOINT);
    void Sigma_Mig(GREEN &G, GREEN &Sigma, GREEN &D0, GREEN &D,
                      GREEN &Pi, GREEN &D0_Pi, GREEN &Pi_D0, CFUNC &g_el_ph, double beta,
                      double h, int SolverOrder);
    
    void Sigma_Mig_sc(int tstp, GREEN &G, GREEN &Sigma, GREEN &D0, GREEN &D,
                      GREEN &Pi, GREEN &D0_Pi, GREEN &Pi_D0, CFUNC &g_el_ph, double beta,
                      double h, int SolverOrder, int MAT_METHOD=CNTR_MAT_FIXPOINT);
    void Sigma_Mig_sc(GREEN &G, GREEN &Sigma, GREEN &D0, GREEN &D,
                      GREEN &Pi, GREEN &D0_Pi, GREEN &Pi_D0, CFUNC &g_el_ph, double beta,
                      double h, int SolverOrder);
    
    void Sigma_uMig(int tstp, GREEN &G, GREEN &D0, CFUNC &g_el_ph, GREEN &Sigma);
    void Sigma_uMig(GREEN &G, GREEN &D0, CFUNC &g_el_ph, GREEN &Sigma,int SolverOrder);
    
    void Sigma_uMig_sc(int tstp, GREEN &G, GREEN &D0, CFUNC &g_el_ph, GREEN &Sigma);
    void Sigma_uMig_sc(GREEN &G, GREEN &D0, CFUNC &g_el_ph, GREEN &Sigma,int SolverOrder);
    
    
    void step_D(int tstp,GREEN &D0, GREEN &D, GREEN &Pi, GREEN &D0_Pi,
      GREEN &Pi_D0, double beta, double h, int SolverOrder,
      int MAT_METHOD=CNTR_MAT_FIXPOINT);
    void start_D(GREEN &D0, GREEN &D, GREEN &Pi, GREEN &D0_Pi,
        GREEN &Pi_D0, double beta, double h, int SolverOrder);
    
    void get_phonon_displace(int tstp,CFUNC &X_ph, CFUNC &n_tot,CFUNC &g_el_ph,GREEN &D0, double w0, int SolverOrder,double h);
    void get_phonon_displace(CFUNC &X_ph, CFUNC &n_tot,CFUNC &g_el_ph,GREEN &D0, double w0, int SolverOrder,double h);
}

#endif
