#ifndef CNTR_GETSET_EXTERN_TEMPLATES_H
#define CNTR_GETSET_EXTERN_TEMPLATES_H

#include "cntr_global_settings.hpp"
#include "cntr_getset_decl.hpp"

namespace cntr {
  extern template void get_mat<double>(const int m, std::complex<double> &G_mat, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  extern template void get_mat<double>(const int m, cdmatrix &G_mat, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  extern template void get_mat<double>(const int m, std::complex<double> &G_mat, herm_matrix<double> &G);
  extern template void get_mat<double>(const int m, cdmatrix &G_mat, herm_matrix<double> &G);
  extern template void get_les<double>(const int i, const int j, std::complex<double> &G_les, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  extern template void get_les<double>(const int i, const int j, cdmatrix &G_les, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  extern template void get_ret<double>(const int i, const int j, std::complex<double> &G_ret, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  extern template void get_ret<double>(const int i, const int j, cdmatrix &G_ret, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  extern template void get_tv<double>(const int i, const int m, std::complex<double> &G_tv, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  extern template void get_tv<double>(const int i, const int m, cdmatrix &G_ret, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  extern template void get_vt<double>(const int m, const int i, std::complex<double> &G_vt, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  extern template void get_vt<double>(const int m, const int i, cdmatrix &G_vt, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  extern template void get_gtr<double>(const int i, const int j, std::complex<double> &G_gtr, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  extern template void get_gtr<double>(const int i, const int j, cdmatrix &G_gtr, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  extern template void get_les<double>(const int i, const int j, std::complex<double> &G_les, herm_matrix<double> &G);
  extern template void get_les<double>(const int i, const int j, cdmatrix &G_les, herm_matrix<double> &G);
  extern template void get_ret<double>(const int i, const int j, std::complex<double> &G_ret, herm_matrix<double> &G);
  extern template void get_ret<double>(const int i, const int j, cdmatrix &G_ret, herm_matrix<double> &G);
  extern template void get_tv<double>(const int i, const int m, std::complex<double> &G_tv, herm_matrix<double> &G);
  extern template void get_tv<double>(const int i, const int m, cdmatrix &G_ret, herm_matrix<double> &G);
  extern template void get_vt<double>(const int m, const int i, std::complex<double> &G_vt, herm_matrix<double> &G);
  extern template void get_vt<double>(const int m, const int i, cdmatrix &G_vt, herm_matrix<double> &G);
  extern template void get_gtr<double>(const int i, const int j, std::complex<double> &G_gtr, herm_matrix<double> &G);
  extern template void get_gtr<double>(const int i, const int j, cdmatrix &G_gtr, herm_matrix<double> &G);

  extern template void get_mat<double>(const int m, std::complex<double> &G_mat, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  extern template void get_mat<double>(const int m, cdmatrix &G_mat, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  extern template void get_mat<double>(const int m, std::complex<double> &G_mat, herm_matrix_timestep<double> &G);
  extern template void get_mat<double>(const int m, cdmatrix &G_mat, herm_matrix_timestep<double> &G);
  extern template void get_les<double>(const int i, const int j, std::complex<double> &G_les, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  extern template void get_les<double>(const int i, const int j, cdmatrix &G_les, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  extern template void get_ret<double>(const int i, const int j, std::complex<double> &G_ret, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  extern template void get_ret<double>(const int i, const int j, cdmatrix &G_ret, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  extern template void get_tv<double>(const int i, const int m, std::complex<double> &G_tv, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  extern template void get_tv<double>(const int i, const int m, cdmatrix &G_ret, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  extern template void get_vt<double>(const int m, const int i, std::complex<double> &G_vt, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  extern template void get_vt<double>(const int m, const int i, cdmatrix &G_vt, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  extern template void get_gtr<double>(const int i, const int j, std::complex<double> &G_gtr, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  extern template void get_gtr<double>(const int i, const int j, cdmatrix &G_gtr, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  extern template void get_les<double>(const int i, const int j, std::complex<double> &G_les, herm_matrix_timestep<double> &G);
  extern template void get_les<double>(const int i, const int j, cdmatrix &G_les, herm_matrix_timestep<double> &G);
  extern template void get_ret<double>(const int i, const int j, std::complex<double> &G_ret, herm_matrix_timestep<double> &G);
  extern template void get_ret<double>(const int i, const int j, cdmatrix &G_ret, herm_matrix_timestep<double> &G);
  extern template void get_tv<double>(const int i, const int m, std::complex<double> &G_tv, herm_matrix_timestep<double> &G);
  extern template void get_tv<double>(const int i, const int m, cdmatrix &G_ret, herm_matrix_timestep<double> &G);
  extern template void get_vt<double>(const int m, const int i, std::complex<double> &G_vt, herm_matrix_timestep<double> &G);
  extern template void get_vt<double>(const int m, const int i, cdmatrix &G_vt, herm_matrix_timestep<double> &G);
  extern template void get_gtr<double>(const int i, const int j, std::complex<double> &G_gtr, herm_matrix_timestep<double> &G);
  extern template void get_gtr<double>(const int i, const int j, cdmatrix &G_gtr, herm_matrix_timestep<double> &G);

  extern template void get_mat<double>(const int m, std::complex<double> &G_mat, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  extern template void get_mat<double>(const int m, cdmatrix &G_mat, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  extern template void get_mat<double>(const int m, std::complex<double> &G_mat, herm_matrix_timestep_view<double> &G);
  extern template void get_mat<double>(const int m, cdmatrix &G_mat, herm_matrix_timestep_view<double> &G);
  extern template void get_les<double>(const int i, const int j, std::complex<double> &G_les, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  extern template void get_les<double>(const int i, const int j, cdmatrix &G_les, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  extern template void get_ret<double>(const int i, const int j, std::complex<double> &G_ret, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  extern template void get_ret<double>(const int i, const int j, cdmatrix &G_ret, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  extern template void get_tv<double>(const int i, const int m, std::complex<double> &G_tv, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  extern template void get_tv<double>(const int i, const int m, cdmatrix &G_ret, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  extern template void get_vt<double>(const int m, const int i, std::complex<double> &G_vt, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  extern template void get_vt<double>(const int m, const int i, cdmatrix &G_vt, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  extern template void get_gtr<double>(const int i, const int j, std::complex<double> &G_gtr, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  extern template void get_gtr<double>(const int i, const int j, cdmatrix &G_gtr, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  extern template void get_les<double>(const int i, const int j, std::complex<double> &G_les, herm_matrix_timestep_view<double> &G);
  extern template void get_les<double>(const int i, const int j, cdmatrix &G_les, herm_matrix_timestep_view<double> &G);
  extern template void get_ret<double>(const int i, const int j, std::complex<double> &G_ret, herm_matrix_timestep_view<double> &G);
  extern template void get_ret<double>(const int i, const int j, cdmatrix &G_ret, herm_matrix_timestep_view<double> &G);
  extern template void get_tv<double>(const int i, const int m, std::complex<double> &G_tv, herm_matrix_timestep_view<double> &G);
  extern template void get_tv<double>(const int i, const int m, cdmatrix &G_ret, herm_matrix_timestep_view<double> &G);
  extern template void get_vt<double>(const int m, const int i, std::complex<double> &G_vt, herm_matrix_timestep_view<double> &G);
  extern template void get_vt<double>(const int m, const int i, cdmatrix &G_vt, herm_matrix_timestep_view<double> &G);
  extern template void get_gtr<double>(const int i, const int j, std::complex<double> &G_gtr, herm_matrix_timestep_view<double> &G);
  extern template void get_gtr<double>(const int i, const int j, cdmatrix &G_gtr, herm_matrix_timestep_view<double> &G);

} // namespace cntr

#endif  // CNTR_GETSET_DECL_H