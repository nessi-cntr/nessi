#ifndef CNTR_UTILITIES_DECL_H
#define CNTR_UTILITIES_DECL_H

#include "cntr_global_settings.hpp"
#include "integration.hpp"

namespace cntr {

  template <typename T> class function;
  template <typename T> class herm_matrix;
  template <typename T> class herm_matrix_timestep;
  template <typename T> class herm_matrix_timestep_view;
  template <typename T> class herm_pseudo;

  /* /////////////////////////////////////////////////////////////////////////////////////////
  // poly extrapolation/interpolation to tstp n+1, using n,...,n-k, with k = integration order
  ///////////////////////////////////////////////////////////////////////////////////////// */

  /// @private
  template <typename T>
  void extrapolate_timestep(int n,herm_matrix<T> &G,integration::Integrator<T> &I);
    template <typename T>
  void extrapolate_timestep(int n, herm_matrix<T> &G, int SolveOrder=MAX_SOLVE_ORDER);
  template <typename T>
  /// @private
  void extrapolate_timestep(int n,herm_pseudo<T> &G,integration::Integrator<T> &I);

  /// @private
  template <typename T>
  void extrapolate_timestep(int n, function<T> &f,integration::Integrator<T> &I);

  template <typename T>
  void extrapolate_timestep(int n, function<T> &f, int ExtrapolationOrder=MAX_SOLVE_ORDER);

  /// @private
  template <typename T>
  cdmatrix interpolation(int tstp,double tinter,function<T> &f,integration::Integrator<T> &I);

  template <typename T>
  cdmatrix interpolation(int tstp,double tinter,function<T> &f,int InterpolationOrder=MAX_SOLVE_ORDER);

  template <typename T> void set_t0_from_mat(herm_matrix<T> &G);
  /// @private
  template <typename T> void set_t0_from_mat(herm_pseudo<T> &G);

  template <typename T> void set_tk_from_mat(herm_matrix<T> &G,int SolveOrder);

  /* /////////////////////////////////////////////////////////////////////////////////////////
  // correlation energy from time-diagonal convolution
  ///////////////////////////////////////////////////////////////////////////////////////// */

  template <typename T>
  T correlation_energy(int tstp, herm_matrix<T> &G, herm_matrix<T> &Sigma, T beta, T h, 
    int SolveOrder=MAX_SOLVE_ORDER);

  /* /////////////////////////////////////////////////////////////////////////////////////////
  // compare to GF's on timestep tstp
  ///////////////////////////////////////////////////////////////////////////////////////// */

  template <typename T> T distance_norm2_ret(int tstp, herm_matrix<T> &g1,herm_matrix<T> &g2);
  template <typename T> T distance_norm2_ret(int tstp, herm_matrix_timestep<T> &g1, herm_matrix<T> &g2);
  template <typename T> T distance_norm2_ret(int tstp, herm_matrix<T> &g1, herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2_ret(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2_ret(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix<T> &g2);
  template <typename T> T distance_norm2_ret(int tstp, herm_matrix<T> &g1, herm_matrix_timestep_view<T> &g2);
  template <typename T> T distance_norm2_ret(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep_view<T> &g2);
  template <typename T> T distance_norm2_ret(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2_ret(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep_view<T> &g2);

  template <typename T> T distance_norm2_tv(int tstp, herm_matrix<T> &g1,herm_matrix<T> &g2);
  template <typename T> T distance_norm2_tv(int tstp, herm_matrix_timestep<T> &g1, herm_matrix<T> &g2);
  template <typename T> T distance_norm2_tv(int tstp, herm_matrix<T> &g1, herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2_tv(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2_tv(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix<T> &g2);
  template <typename T> T distance_norm2_tv(int tstp, herm_matrix<T> &g1, herm_matrix_timestep_view<T> &g2);
  template <typename T> T distance_norm2_tv(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep_view<T> &g2);
  template <typename T> T distance_norm2_tv(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2_tv(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep_view<T> &g2);

  template <typename T> T distance_norm2_les(int tstp, herm_matrix<T> &g1,herm_matrix<T> &g2);
  template <typename T> T distance_norm2_les(int tstp, herm_matrix_timestep<T> &g1, herm_matrix<T> &g2);
  template <typename T> T distance_norm2_les(int tstp, herm_matrix<T> &g1, herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2_les(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2_les(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix<T> &g2);
  template <typename T> T distance_norm2_les(int tstp, herm_matrix<T> &g1, herm_matrix_timestep_view<T> &g2);
  template <typename T> T distance_norm2_les(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep_view<T> &g2);
  template <typename T> T distance_norm2_les(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2_les(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep_view<T> &g2);

  template <typename T> T distance_norm2(int tstp, herm_matrix<T> &g1,herm_matrix<T> &g2);
  template <typename T> T distance_norm2(int tstp, herm_matrix_timestep<T> &g1, herm_matrix<T> &g2);
  template <typename T> T distance_norm2(int tstp, herm_matrix<T> &g1, herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix<T> &g2);
  template <typename T> T distance_norm2(int tstp, herm_matrix<T> &g1, herm_matrix_timestep_view<T> &g2);
  template <typename T> T distance_norm2(int tstp, herm_matrix_timestep<T> &g1, herm_matrix_timestep_view<T> &g2);
  template <typename T> T distance_norm2(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2(int tstp, herm_matrix_timestep_view<T> &g1, herm_matrix_timestep_view<T> &g2);

  // legacy interfaces
  template <typename T> T distance_norm2(herm_matrix_timestep_view<T> &g1,herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2(herm_matrix_timestep<T> &g1,herm_matrix_timestep<T> &g2);
  template <typename T> T distance_norm2_eigen(int tstp,herm_matrix_timestep<T> &g1,herm_matrix<T> &g2);

  /* /////////////////////////////////////////////////////////////////////////////////////////
  // evaluate memory for various objects:
  ///////////////////////////////////////////////////////////////////////////////////////// */

  template <typename T> size_t mem_herm_matrix(int nt,int ntau,int size);

  template <typename T> size_t mem_herm_matrix_timestep(int tstp,int ntau,int size);

  template <typename T> size_t mem_function(int nt,int size);

  template <typename T> void force_matsubara_hermitian(herm_matrix<T> &G);

  template <class GG> void force_matsubara_hermitian(GG &G);

}  // namespace cntr

#endif  // CNTR_UTILITIES_DECL_H
