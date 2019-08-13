#ifndef CNTR_PSEUDO_CONVOLUTION_EXTERN_TEMPLATES_H
#define CNTR_PSEUDO_CONVOLUTION_EXTERN_TEMPLATES_H

#include "cntr_pseudo_convolution_decl.hpp"

/*
namespace cntr {
#if CNTR_USE_OMP == 1
extern template void
convolution_timestep_omp<double>(int omp_num_threads, int tstp, herm_matrix<double> &C,
                                 herm_matrix<double> &A, herm_matrix<double> &B,
                                 integration::Integrator<double> &I, double beta, double h);
extern template void
convolution_omp<double>(int omp_num_threads, herm_matrix<double> &C, herm_matrix<double> &A,
                        herm_matrix<double> &Acc, herm_matrix<double> &B,
                        herm_matrix<double> &Bcc, integration::Integrator<double> &I,
                        double beta, double h);
#endif // CNTR_USE_OMP

}  // namespace cntr
*/
#endif  // CNTR_PSEUDO_CONVOLUTION_EXTERN_TEMPLATES_H
