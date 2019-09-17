#ifndef CNTR_DYSON_OMP_EXTERN_TEMPLATES_H
#define CNTR_DYSON_OMP_EXTERN_TEMPLATES_H

#include "cntr_dyson_omp_decl.hpp"

namespace cntr {

#if CNTR_USE_OMP == 1
extern template void dyson_timestep_omp<double>(int omp_num_threads, int n, herm_matrix<double> &G,
	double mu, function<double> &H, herm_matrix<double> &Sigma, integration::Integrator<double> &I,
	double beta, double h);
extern template void dyson_timestep_omp<double>(int omp_num_threads, int n, herm_matrix<double> &G,
	double mu, function<double> &H, herm_matrix<double> &Sigma,
	double beta, double h, int SolveOrder);
#endif // CNTR_USE_OMP

}  // namespace cntr

#endif  // CNTR_DYSON_OMP_EXTERN_TEMPLATES_H
