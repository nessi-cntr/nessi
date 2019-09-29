
#include "cntr_global_settings.hpp"
#include "cntr_equilibrium_extern_templates.hpp"
#include "cntr_equilibrium_impl.hpp"

namespace cntr {

/// @private
  template double fermi<double>(double beta,double omega);
/// @private
  template double fermi_exp<double>(double beta,double tau,double omega);
/// @private
  template dvector fermi<double>(double beta,dvector &omega);
/// @private
  template dvector fermi_exp<double>(double beta,double tau,dvector &omega);
  /// @private
  template cdmatrix diag_prop<double>(double time,dvector &omega);

  // prefered interfaces
  /// @private
  template void green_from_H<double>(herm_matrix<double> &G,double mu,cdmatrix &eps,double beta,double h);
  /// @private
  template void green_from_H<double>(int tstp, herm_matrix_timestep<double> &G,double mu,cdmatrix &eps,double beta,double h);
  /// @private
  template void green_from_H<double>(int tstp, herm_matrix<double> &G,double mu,cdmatrix &eps,double beta,double h);
  /// @private
  template void green_from_H<double>(herm_matrix<double> &G,double mu,cntr::function<double> &eps,
    double beta,double h,int SolveOrder,int cf_order);
  /// @private
  template void green_from_H<double>(int tstp, herm_matrix_timestep<double> &G,double mu,cntr::function<double> &eps,
           double beta,double h,bool fixHam,int SolveOrder,int cf_order);
  /// @private
  template void green_from_H<double>(int tstp, herm_matrix<double> &G,double mu,cntr::function<double> &eps,
           double beta,double h,bool fixHam,int SolveOrder,int cf_order);

  // legacy interfaces 
  /// @private
  template void green_from_H<double>(herm_matrix<double> &G,double mu,cntr::function<double> &eps,
    double beta,double h,bool fixHam,int SolveOrder,int cf_order);
  /// @private
  template void green_from_H<double>(herm_matrix_timestep<double> &G,double mu,cntr::function<double> &eps,
           double beta,double h,bool fixHam,int SolveOrder,int cf_order);




  // template void green_single_pole_bose(herm_matrix<double> &G, double *w, double beta, double h);
  /// @private
  template void green_single_pole_XX_timestep(herm_matrix_timestep<double> &D0,
					      double w, double beta, double h);
  /// @private
  template void green_single_pole_XX_timestep(int tstp, herm_matrix<double> &D0,
					      double w, double beta,double h);
  /// @private
  template void green_single_pole_XX(herm_matrix<double> &D0, double w,
				     double beta, double h);

}  // namespace cntr
