
#include "cntr_global_settings.hpp"
#include "cntr_getset_extern_templates.hpp"
#include "cntr_getset_impl.hpp"

namespace cntr {
  template void get_mat<double>(const int m, std::complex<double> &G_mat, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  template void get_mat<double>(const int m, cdmatrix &G_mat, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  template void get_mat<double>(const int m, std::complex<double> &G_mat, herm_matrix<double> &G);
  template void get_mat<double>(const int m, cdmatrix &G_mat, herm_matrix<double> &G);
  template void get_les<double>(const int i, const int j, std::complex<double> &G_les, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  template void get_les<double>(const int i, const int j, cdmatrix &G_les, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  template void get_ret<double>(const int i, const int j, std::complex<double> &G_ret, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  template void get_ret<double>(const int i, const int j, cdmatrix &G_ret, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  template void get_tv<double>(const int i, const int m, std::complex<double> &G_tv, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  template void get_tv<double>(const int i, const int m, cdmatrix &G_ret, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  template void get_vt<double>(const int m, const int i, std::complex<double> &G_vt, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  template void get_vt<double>(const int m, const int i, cdmatrix &G_vt, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  template void get_gtr<double>(const int i, const int j, std::complex<double> &G_gtr, herm_matrix<double> &G, herm_matrix<double> &Gcc);
  template void get_gtr<double>(const int i, const int j, cdmatrix &G_gtr, herm_matrix<double> &G, herm_matrix<double> &Gcc);

  template void get_mat<double>(const int m, std::complex<double> &G_mat, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  template void get_mat<double>(const int m, cdmatrix &G_mat, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  template void get_mat<double>(const int m, std::complex<double> &G_mat, herm_matrix_timestep<double> &G);
  template void get_mat<double>(const int m, cdmatrix &G_mat, herm_matrix_timestep<double> &G);
  template void get_les<double>(const int i, const int j, std::complex<double> &G_les, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  template void get_les<double>(const int i, const int j, cdmatrix &G_les, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  template void get_ret<double>(const int i, const int j, std::complex<double> &G_ret, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  template void get_ret<double>(const int i, const int j, cdmatrix &G_ret, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  template void get_tv<double>(const int i, const int m, std::complex<double> &G_tv, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  template void get_tv<double>(const int i, const int m, cdmatrix &G_ret, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  template void get_vt<double>(const int m, const int i, std::complex<double> &G_vt, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  template void get_vt<double>(const int m, const int i, cdmatrix &G_vt, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  template void get_gtr<double>(const int i, const int j, std::complex<double> &G_gtr, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);
  template void get_gtr<double>(const int i, const int j, cdmatrix &G_gtr, herm_matrix_timestep<double> &G, herm_matrix_timestep<double> &Gcc);

  template void get_mat<double>(const int m, std::complex<double> &G_mat, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  template void get_mat<double>(const int m, cdmatrix &G_mat, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  template void get_mat<double>(const int m, std::complex<double> &G_mat, herm_matrix_timestep_view<double> &G);
  template void get_mat<double>(const int m, cdmatrix &G_mat, herm_matrix_timestep_view<double> &G);
  template void get_les<double>(const int i, const int j, std::complex<double> &G_les, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  template void get_les<double>(const int i, const int j, cdmatrix &G_les, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  template void get_ret<double>(const int i, const int j, std::complex<double> &G_ret, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  template void get_ret<double>(const int i, const int j, cdmatrix &G_ret, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  template void get_tv<double>(const int i, const int m, std::complex<double> &G_tv, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  template void get_tv<double>(const int i, const int m, cdmatrix &G_ret, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  template void get_vt<double>(const int m, const int i, std::complex<double> &G_vt, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  template void get_vt<double>(const int m, const int i, cdmatrix &G_vt, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  template void get_gtr<double>(const int i, const int j, std::complex<double> &G_gtr, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);
  template void get_gtr<double>(const int i, const int j, cdmatrix &G_gtr, herm_matrix_timestep_view<double> &G, herm_matrix_timestep_view<double> &Gcc);

} // namespace cntr

