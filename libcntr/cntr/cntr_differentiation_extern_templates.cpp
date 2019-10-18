#include "cntr_differentiation_extern_templates.hpp"
#include "cntr_differentiation_impl.hpp"

namespace cntr {

//
// DIFFERENCIATION /////////////////////////////////////////
template void deriv1_matsubara<double>(herm_matrix<double> &dA,herm_matrix<double> &A, integration::Integrator<double> &I, double beta);
template void deriv1_element<double>(int tstp1,int tstp2,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, integration::Integrator<double> &I,double h);
template void deriv2_element<double>(int tstp1,int tstp2,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, integration::Integrator<double> &I,double h);
template void deriv1_tv<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A, integration::Integrator<double> &I,double h);
template void deriv2_tv<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A, integration::Integrator<double> &I,double beta);
template void deriv1_timestep<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, integration::Integrator<double> &I, double beta,double h);
template void deriv2_timestep<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, integration::Integrator<double> &I, double beta,double h);


template void deriv1_matsubara<double>(herm_matrix<double> &dA,herm_matrix<double> &A, double beta, int SolveOrder);
template void deriv1_element<double>(int tstp1,int tstp2,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc,double h, int SolveOrder);
template void deriv2_element<double>(int tstp1,int tstp2,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc,double h, int SolveOrder);
template void deriv1_tv<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A, double h, int SolveOrder);
template void deriv2_tv<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A, double beta, int SolveOrder);
template void deriv1_timestep<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, double beta,double h, int SolveOrder);
template void deriv2_timestep<double>(int tstp,herm_matrix<double> &dA,herm_matrix<double> &A,herm_matrix<double> &Acc, double beta,double h, int SolveOrder);


}  // namespace cntr
