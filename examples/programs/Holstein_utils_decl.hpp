#ifndef HOLS_UTILS
#define HOLS_UTILS

#include <sys/stat.h>
#include <complex>
#include "cntr/cntr.hpp"

using namespace cntr;
using namespace std;

#define CFUNC cntr::function<double>
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define CPLX complex<double>

namespace Hols {
    
    void evaluate_phonon_energy_qu(CFUNC & Eph_corr_t, GREEN & D, GREEN & Pi, int SolverOrder, double beta, double h, double w0);
    void get_phonon_momentum(CFUNC &P_ph, CFUNC &n_tot,CFUNC &g_el_ph, GREEN &D0, double w0, int SolverOrder, double h);
    void evaluate_phonon_energy_cl(CFUNC & Eph_t, CFUNC &X_ph, CFUNC &P_ph ,double w0);
}

#endif
