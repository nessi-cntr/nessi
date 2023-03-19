#ifndef CNTR_EQUILIBRIUM_EXTERN_TEMPLATES_H
#define CNTR_EQUILIBRIUM_EXTERN_TEMPLATES_H

#include "cntr_global_settings.hpp"
#include "cntr_equilibrium_decl.hpp"

namespace cntr {

  extern template double fermi<double>(double beta,double omega);
  extern template double fermi_exp<double>(double beta,double tau,double omega);
  extern template dvector fermi<double>(double beta,dvector &omega);
  extern template dvector fermi_exp<double>(double beta,double tau,dvector &omega);
  extern template cdmatrix diag_prop<double>(double time,dvector &omega);
  extern template void green_equilibrium_mat_bethe<double>(herm_matrix<double> &G,double beta,int limit,int nn,double mu,double Eshift);
  extern template void green_equilibrium_bethe<double>(herm_matrix<double> &G,double beta,double h,int limit,int nn,double mu,double Eshift);


  // prefered interfaces
  extern template void green_from_H<double>(herm_matrix<double> &G,double mu,cdmatrix &eps,double beta,double h);
  extern template void green_from_H<double>(int tstp, herm_matrix_timestep<double> &G,double mu,cdmatrix &eps,double beta,double h);
  extern template void green_from_H<double>(int tstp, herm_matrix<double> &G,double mu,cdmatrix &eps,double beta,double h);
  extern template void green_from_H<double>(herm_matrix<double> &G,double mu,cntr::function<double> &eps,
    double beta,double h,int SolveOrder,int cf_order);
  extern template void green_from_H<double>(int tstp, herm_matrix_timestep<double> &G,double mu,cntr::function<double> &eps,
           double beta,double h,bool fixHam,int SolveOrder,int cf_order);
  extern template void green_from_H<double>(int tstp, herm_matrix<double> &G,double mu,cntr::function<double> &eps,
           double beta,double h,bool fixHam,int SolveOrder,int cf_order);

  // legacy interfaces 
  extern template void green_from_H<double>(herm_matrix<double> &G,double mu,cntr::function<double> &eps,
    double beta,double h,bool fixHam,int SolveOrder,int cf_order);
  extern template void green_from_H<double>(herm_matrix_timestep<double> &G,double mu,cntr::function<double> &eps,
			     double beta,double h,bool fixHam,int SolveOrder,int cf_order);



  extern template void green_single_pole_XX_timestep(herm_matrix_timestep<double> &D0,
						     double w, double beta, double h);

  extern template void green_single_pole_XX_timestep(int tstp, herm_matrix_timestep<double> &D0,
                 double w, double beta, double h);
  extern template void green_single_pole_XX_timestep(int tstp, herm_matrix<double> &D0,
						     double w, double beta,double h);
  extern template void green_single_pole_XX(herm_matrix<double> &D0, double w,
					    double beta, double h);
  
}  // namespace cntr

#endif  // CNTR_EQUILIBRIUM_EXTERN_TEMPLATES_H
