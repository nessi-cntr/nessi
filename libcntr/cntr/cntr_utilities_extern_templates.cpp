#include "cntr_utilities_extern_templates.hpp"
#include "cntr_utilities_impl.hpp"

namespace cntr {

  ///@private
  template void extrapolate_timestep<double>(int n,herm_matrix<double> &G,integration::Integrator<double> &I);
  ///@private
  template void extrapolate_timestep<double>(int n,herm_pseudo<double> &G,integration::Integrator<double> &I);

  ///@private
  template void extrapolate_timestep(int n, function<double> &f,integration::Integrator<double> &I);
  ///@private
  template cdmatrix interpolation(int tstp,double tinter,function<double> &f,integration::Integrator<double> &I);
  
  ///@private
  template void set_t0_from_mat<double>(herm_matrix<double> &G);
  ///@private
  template void set_t0_from_mat<double>(herm_pseudo<double> &G);
  ///@private
  template void set_tk_from_mat<double>(herm_matrix<double> &G,int kt);

  // correlation energy from time-diagonal convolution
  /// @private
  template
  double correlation_energy(int tstp, herm_matrix<double> &G, herm_matrix<double> &Sigma,
		       integration::Integrator<double> &I, double beta, double h);
  
  /// @private
  template double distance_norm2<double>(int tstp,herm_matrix<double> &g1,herm_matrix<double> &g2);
  /// @private
  template double distance_norm2_ret<double>(int tstp,herm_matrix<double> &g1,herm_matrix<double> &g2);
  /// @private
  template double distance_norm2_tv<double>(int tstp,herm_matrix<double> &g1,herm_matrix<double> &g2);
  /// @private
  template double distance_norm2_les<double>(int tstp,herm_matrix<double> &g1,herm_matrix<double> &g2);
  /// @private
  template double distance_norm2<double>(int tstp,herm_matrix_timestep<double> &g1,herm_matrix<double> &g2);
  /// @private
  template double distance_norm2<double>(int tstp,herm_matrix_timestep_view<double> &g1,herm_matrix<double> &g2);
  /// @private
  template double distance_norm2<double>(herm_matrix_timestep_view<double> &g1,herm_matrix_timestep<double> &g2);
  /// @private
  template double distance_norm2<double>(int tstp, herm_matrix_timestep_view<double> &g1,herm_matrix_timestep<double> &g2);
  /// @private
  template double distance_norm2<double>(herm_matrix_timestep<double> &g1,herm_matrix_timestep<double> &g2);
  /// @private
  template double distance_norm2<double>(int tstp, herm_matrix_timestep<double> &g1,herm_matrix_timestep<double> &g2);
  /// @private
  template double distance_norm2_eigen<double>(int tstp, herm_matrix_timestep<double> &g1,herm_matrix<double> &g2);
  /// @private
  template size_t mem_herm_matrix<double>(int nt,int ntau,int size);
  /// @private
  template size_t mem_herm_matrix_timestep<double>(int tstp,int ntau,int size);
  /// @private
  template size_t mem_function<double>(int nt,int size);

  /// @private
  template void force_matsubara_hermitian(herm_matrix<double> &G);

  
}  // namespace cntr
