#ifndef CNTR_HERM_MATRIX_EXTERN_TEMPLATES_H
#define CNTR_HERM_MATRIX_EXTERN_TEMPLATES_H

#include "cntr_herm_matrix_decl.hpp"
#include <Eigen/Core>

#define GREEN cntr::herm_matrix<double>

namespace cntr {

extern template class herm_matrix<double>;

extern template void herm_matrix<double>::get_les<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M) const;
extern template void herm_matrix<double>::get_gtr<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M) const;
extern template void herm_matrix<double>::get_ret<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M) const;
extern template void herm_matrix<double>::get_vt<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M) const;
extern template void herm_matrix<double>::get_tv<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M) const;
extern template void herm_matrix<double>::get_mat<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M) const;
extern template void herm_matrix<double>::get_matminus<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M) const;
extern template void herm_matrix<double>::set_les<Eigen::MatrixXcd>(int i, int j, Eigen::MatrixXcd &M);
extern template void herm_matrix<double>::set_ret<Eigen::MatrixXcd>(int i, int j, Eigen::MatrixXcd &M);
extern template void herm_matrix<double>::set_tv<Eigen::MatrixXcd>(int i, int j, Eigen::MatrixXcd &M);
extern template void herm_matrix<double>::set_mat<Eigen::MatrixXcd>(int i, Eigen::MatrixXcd &M);

}  // namespace cntr

#endif  // CNTR_HERM_MATRIX_EXTERN_TEMPLATES_H
