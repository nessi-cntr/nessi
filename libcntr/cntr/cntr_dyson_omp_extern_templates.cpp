#include "cntr_dyson_omp_extern_templates.hpp"
#include "cntr_dyson_omp_impl.hpp"

namespace cntr {

#if CNTR_USE_OMP == 1
template void dyson_timestep_omp<double>(int omp_num_threads, int n, herm_matrix<double> &G, 
	double mu, function<double> &H, herm_matrix<double> &Sigma, integration::Integrator<double> &I,
	double beta, double h);
template void dyson_timestep_omp<double>(int omp_num_threads, int n, herm_matrix<double> &G, 
	double mu, function<double> &H, herm_matrix<double> &Sigma,
	double beta, double h, int SolveOrder);
#endif // CNTR_USE_OMP

}  // namespace cntr
