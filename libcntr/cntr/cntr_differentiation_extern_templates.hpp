#ifndef CNTR_DIFFERENTIATION_EXTERN_TEMPLATES_H
#define CNTR_DIFFERENTIATION_EXTERN_TEMPLATES_H

#include "cntr_differentiation_decl.hpp"

namespace cntr {

extern template void deriv1_matsubara<double>(herm_matrix<double> &dA,herm_matrix<double> &A, integration::Integrator<double> &I, double beta);
extern template void deriv1_element<double>(int tstp1,int tstp2,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, integration::Integrator<double> &I,double h);
extern template void deriv2_element<double>(int tstp1,int tstp2,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, integration::Integrator<double> &I,double h);
extern template void deriv1_tv<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A, integration::Integrator<double> &I,double h);
extern template void deriv2_tv<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A, integration::Integrator<double> &I,double beta);
extern template void deriv1_timestep<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, integration::Integrator<double> &I, double beta,double h);
extern template void deriv2_timestep<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, integration::Integrator<double> &I, double beta,double h);

extern template void deriv1_matsubara<double>(herm_matrix<double> &dA,herm_matrix<double> &A,double beta, int SolveOrder);
extern template void deriv1_element<double>(int tstp1,int tstp2,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, double h, int SolveOrder);
extern template void deriv2_element<double>(int tstp1,int tstp2,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, double h, int SolveOrder);
extern template void deriv1_tv<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A, double h, int SolveOrder);
extern template void deriv2_tv<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A, double beta, int SolveOrder);
extern template void deriv1_timestep<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, double beta,double h, int SolveOrder);
extern template void deriv2_timestep<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, double beta,double h, int SolveOrder);


}  // namespace cntr

#endif  // CNTR_DIFFERENTIATION_EXTERN_TEMPLATES_H
