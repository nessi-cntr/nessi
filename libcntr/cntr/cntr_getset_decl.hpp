#ifndef CNTR_GETSET_DECL_H
#define CNTR_GETSET_DECL_H

#include "cntr_global_settings.hpp"
#include "eigen_typedef.h"

namespace cntr {

  template <typename T> class herm_matrix;
  template <typename T> class herm_matrix_timestep;
  template <typename T> class herm_matrix_timestep_view;

  template<typename T> void get_mat(const int m, std::complex<T> &G_mat, herm_matrix<T> &G, herm_matrix<T> &Gcc);
  template<typename T> void get_mat(const int m, cdmatrix &G_mat, herm_matrix<T> &G, herm_matrix<T> &Gcc);
  template<typename T> void get_mat(const int m, std::complex<T> &G_mat, herm_matrix<T> &G);
  template<typename T> void get_mat(const int m, cdmatrix &G_mat, herm_matrix<T> &G);
  template<typename T> void get_les(const int i, const int j, std::complex<T> &G_les, herm_matrix<T> &G, herm_matrix<T> &Gcc);
  template<typename T> void get_les(const int i, const int j, cdmatrix &G_les, herm_matrix<T> &G, herm_matrix<T> &Gcc);
  template<typename T> void get_ret(const int i, const int j, std::complex<T> &G_ret, herm_matrix<T> &G, herm_matrix<T> &Gcc);
  template<typename T> void get_ret(const int i, const int j, cdmatrix &G_ret, herm_matrix<T> &G, herm_matrix<T> &Gcc);
  template<typename T> void get_tv(const int i, const int m, std::complex<T> &G_tv, herm_matrix<T> &G, herm_matrix<T> &Gcc);
  template<typename T> void get_tv(const int i, const int m, cdmatrix &G_ret, herm_matrix<T> &G, herm_matrix<T> &Gcc);
  template<typename T> void get_vt(const int m, const int i, std::complex<T> &G_vt, herm_matrix<T> &G, herm_matrix<T> &Gcc);
  template<typename T> void get_vt(const int m, const int i, cdmatrix &G_vt, herm_matrix<T> &G, herm_matrix<T> &Gcc);
  template<typename T> void get_gtr(const int i, const int j, std::complex<T> &G_gtr, herm_matrix<T> &G, herm_matrix<T> &Gcc);
  template<typename T> void get_gtr(const int i, const int j, cdmatrix &G_gtr, herm_matrix<T> &G, herm_matrix<T> &Gcc);

  template<typename T> void get_les(const int i, const int j, std::complex<T> &G_les, herm_matrix<T> &G);
  template<typename T> void get_les(const int i, const int j, cdmatrix &G_les, herm_matrix<T> &G);
  template<typename T> void get_ret(const int i, const int j, std::complex<T> &G_ret, herm_matrix<T> &G);
  template<typename T> void get_ret(const int i, const int j, cdmatrix &G_ret, herm_matrix<T> &G);
  template<typename T> void get_tv(const int i, const int m, std::complex<T> &G_tv, herm_matrix<T> &G);
  template<typename T> void get_tv(const int i, const int m, cdmatrix &G_ret, herm_matrix<T> &G);
  template<typename T> void get_vt(const int m, const int i, std::complex<T> &G_vt, herm_matrix<T> &G);
  template<typename T> void get_vt(const int m, const int i, cdmatrix &G_vt, herm_matrix<T> &G);
  template<typename T> void get_gtr(const int i, const int j, std::complex<T> &G_gtr, herm_matrix<T> &G);
  template<typename T> void get_gtr(const int i, const int j, cdmatrix &G_gtr, herm_matrix<T> &G);



  template<typename T> void get_mat(const int m, std::complex<T> &G_mat, herm_matrix_timestep<T> &G, herm_matrix_timestep<T> &Gcc);
  template<typename T> void get_mat(const int m, cdmatrix &G_mat, herm_matrix_timestep<T> &G, herm_matrix_timestep<T> &Gcc);
  template<typename T> void get_mat(const int m, std::complex<T> &G_mat, herm_matrix_timestep<T> &G);
  template<typename T> void get_mat(const int m, cdmatrix &G_mat, herm_matrix_timestep<T> &G);
  template<typename T> void get_les(const int i, const int j, std::complex<T> &G_les, herm_matrix_timestep<T> &G, herm_matrix_timestep<T> &Gcc);
  template<typename T> void get_les(const int i, const int j, cdmatrix &G_les, herm_matrix_timestep<T> &G, herm_matrix_timestep<T> &Gcc);
  template<typename T> void get_ret(const int i, const int j, std::complex<T> &G_ret, herm_matrix_timestep<T> &G, herm_matrix_timestep<T> &Gcc);
  template<typename T> void get_ret(const int i, const int j, cdmatrix &G_ret, herm_matrix_timestep<T> &G, herm_matrix_timestep<T> &Gcc);
  template<typename T> void get_tv(const int i, const int m, std::complex<T> &G_tv, herm_matrix_timestep<T> &G, herm_matrix_timestep<T> &Gcc);
  template<typename T> void get_tv(const int i, const int m, cdmatrix &G_ret, herm_matrix_timestep<T> &G, herm_matrix_timestep<T> &Gcc);
  template<typename T> void get_vt(const int m, const int i, std::complex<T> &G_vt, herm_matrix_timestep<T> &G, herm_matrix_timestep<T> &Gcc);
  template<typename T> void get_vt(const int m, const int i, cdmatrix &G_vt, herm_matrix_timestep<T> &G, herm_matrix_timestep<T> &Gcc);
  template<typename T> void get_gtr(const int i, const int j, std::complex<T> &G_gtr, herm_matrix_timestep<T> &G, herm_matrix_timestep<T> &Gcc);
  template<typename T> void get_gtr(const int i, const int j, cdmatrix &G_gtr, herm_matrix_timestep<T> &G, herm_matrix_timestep<T> &Gcc);

  template<typename T> void get_les(const int i, const int j, std::complex<T> &G_les, herm_matrix_timestep<T> &G);
  template<typename T> void get_les(const int i, const int j, cdmatrix &G_les, herm_matrix_timestep<T> &G);
  template<typename T> void get_ret(const int i, const int j, std::complex<T> &G_ret, herm_matrix_timestep<T> &G);
  template<typename T> void get_ret(const int i, const int j, cdmatrix &G_ret, herm_matrix_timestep<T> &G);
  template<typename T> void get_tv(const int i, const int m, std::complex<T> &G_tv, herm_matrix_timestep<T> &G);
  template<typename T> void get_tv(const int i, const int m, cdmatrix &G_ret, herm_matrix_timestep<T> &G);
  template<typename T> void get_vt(const int m, const int i, std::complex<T> &G_vt, herm_matrix_timestep<T> &G);
  template<typename T> void get_vt(const int m, const int i, cdmatrix &G_vt, herm_matrix_timestep<T> &G);
  template<typename T> void get_gtr(const int i, const int j, std::complex<T> &G_gtr, herm_matrix_timestep<T> &G);
  template<typename T> void get_gtr(const int i, const int j, cdmatrix &G_gtr, herm_matrix_timestep<T> &G);



  template<typename T> void get_mat(const int m, std::complex<T> &G_mat, herm_matrix_timestep_view<T> &G, herm_matrix_timestep_view<T> &Gcc);
  template<typename T> void get_mat(const int m, cdmatrix &G_mat, herm_matrix_timestep_view<T> &G, herm_matrix_timestep_view<T> &Gcc);
  template<typename T> void get_mat(const int m, std::complex<T> &G_mat, herm_matrix_timestep_view<T> &G);
  template<typename T> void get_mat(const int m, cdmatrix &G_mat, herm_matrix_timestep_view<T> &G);
  template<typename T> void get_les(const int i, const int j, std::complex<T> &G_les, herm_matrix_timestep_view<T> &G, herm_matrix_timestep_view<T> &Gcc);
  template<typename T> void get_les(const int i, const int j, cdmatrix &G_les, herm_matrix_timestep_view<T> &G, herm_matrix_timestep_view<T> &Gcc);
  template<typename T> void get_ret(const int i, const int j, std::complex<T> &G_ret, herm_matrix_timestep_view<T> &G, herm_matrix_timestep_view<T> &Gcc);
  template<typename T> void get_ret(const int i, const int j, cdmatrix &G_ret, herm_matrix_timestep_view<T> &G, herm_matrix_timestep_view<T> &Gcc);
  template<typename T> void get_tv(const int i, const int m, std::complex<T> &G_tv, herm_matrix_timestep_view<T> &G, herm_matrix_timestep_view<T> &Gcc);
  template<typename T> void get_tv(const int i, const int m, cdmatrix &G_ret, herm_matrix_timestep_view<T> &G, herm_matrix_timestep_view<T> &Gcc);
  template<typename T> void get_vt(const int m, const int i, std::complex<T> &G_vt, herm_matrix_timestep_view<T> &G, herm_matrix_timestep_view<T> &Gcc);
  template<typename T> void get_vt(const int m, const int i, cdmatrix &G_vt, herm_matrix_timestep_view<T> &G, herm_matrix_timestep_view<T> &Gcc);
  template<typename T> void get_gtr(const int i, const int j, std::complex<T> &G_gtr, herm_matrix_timestep_view<T> &G, herm_matrix_timestep_view<T> &Gcc);
  template<typename T> void get_gtr(const int i, const int j, cdmatrix &G_gtr, herm_matrix_timestep_view<T> &G, herm_matrix_timestep_view<T> &Gcc);

  template<typename T> void get_les(const int i, const int j, std::complex<T> &G_les, herm_matrix_timestep_view<T> &G);
  template<typename T> void get_les(const int i, const int j, cdmatrix &G_les, herm_matrix_timestep_view<T> &G);
  template<typename T> void get_ret(const int i, const int j, std::complex<T> &G_ret, herm_matrix_timestep_view<T> &G);
  template<typename T> void get_ret(const int i, const int j, cdmatrix &G_ret, herm_matrix_timestep_view<T> &G);
  template<typename T> void get_tv(const int i, const int m, std::complex<T> &G_tv, herm_matrix_timestep_view<T> &G);
  template<typename T> void get_tv(const int i, const int m, cdmatrix &G_ret, herm_matrix_timestep_view<T> &G);
  template<typename T> void get_vt(const int m, const int i, std::complex<T> &G_vt, herm_matrix_timestep_view<T> &G);
  template<typename T> void get_vt(const int m, const int i, cdmatrix &G_vt, herm_matrix_timestep_view<T> &G);
  template<typename T> void get_gtr(const int i, const int j, std::complex<T> &G_gtr, herm_matrix_timestep_view<T> &G);
  template<typename T> void get_gtr(const int i, const int j, cdmatrix &G_gtr, herm_matrix_timestep_view<T> &G);


} // namespace cntr

#endif  // CNTR_GETSET_DECL_H
