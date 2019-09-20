#ifndef CNTR_HERM_MATRIX_TIMESTEP_EXTERN_TEMPLATES_H
#define CNTR_HERM_MATRIX_TIMESTEP_EXTERN_TEMPLATES_H

#include <Eigen/Core>
#include "cntr_herm_matrix_timestep_decl.hpp"

#define GREEN_TSTP cntr::herm_matrix_timestep<double>

namespace cntr {

  // preferred interfaces
  extern template void herm_matrix_timestep<double>::set_ret<Eigen::MatrixXcd>(int i, int j,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::set_les<Eigen::MatrixXcd>(int i, int j,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::set_tv<Eigen::MatrixXcd>(int i, int j,Eigen::MatrixXcd &M); 
  extern template void herm_matrix_timestep<double>::set_mat<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M);

  extern template void herm_matrix_timestep<double>::get_les<Eigen::MatrixXcd>(int i, int j,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::get_ret<Eigen::MatrixXcd>(int i, int j,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::get_tv<Eigen::MatrixXcd>(int i, int mtau ,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::get_mat<Eigen::MatrixXcd>(int i, Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::get_matminus<Eigen::MatrixXcd>(int i, Eigen::MatrixXcd &M);

  extern template void herm_matrix_timestep<double>::density_matrix<Eigen::MatrixXcd>(int tstp, Eigen::MatrixXcd &rho);

  // legacy interfaces

  extern template void herm_matrix_timestep<double>::set_ret<Eigen::MatrixXcd>(int j,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::set_les<Eigen::MatrixXcd>(int j,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::set_tv<Eigen::MatrixXcd>(int j,Eigen::MatrixXcd &M); 

  extern template void herm_matrix_timestep<double>::get_les_t_tstp<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::get_les_tstp_t<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::get_ret_tstp_t<Eigen::MatrixXcd>(int j,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::get_ret_t_tstp<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::get_tv<Eigen::MatrixXcd>(int j,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::get_vt<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M,int sig);
  extern template void herm_matrix_timestep<double>::get_matminus<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M,int sig);
  extern template void herm_matrix_timestep<double>::get_gtr_tstp_t<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M);
  extern template void herm_matrix_timestep<double>::get_gtr_t_tstp<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M);
  
  extern template void herm_matrix_timestep<double>::density_matrix<Eigen::MatrixXcd>(Eigen::MatrixXcd &rho);

}  // namespace cntr

#endif  // CNTR_HERM_MATRIX_TIMESTEP_EXTERN_TEMPLATES_H
