#include "cntr_herm_matrix_extern_templates.hpp"
#include "cntr_herm_matrix_impl.hpp"

namespace cntr {

template class herm_matrix<double>;

template void herm_matrix<double>::get_les<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M) const;
template void herm_matrix<double>::get_gtr<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M) const;
template void herm_matrix<double>::get_ret<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M) const;
template void herm_matrix<double>::get_vt<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M) const;
template void herm_matrix<double>::get_tv<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M) const;
template void herm_matrix<double>::get_mat<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M) const;
template void herm_matrix<double>::get_matminus<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M) const;
template void herm_matrix<double>::set_les<Eigen::MatrixXcd>(int i, int j, Eigen::MatrixXcd &M);
template void herm_matrix<double>::set_ret<Eigen::MatrixXcd>(int i, int j, Eigen::MatrixXcd &M);
template void herm_matrix<double>::set_tv<Eigen::MatrixXcd>(int i, int j, Eigen::MatrixXcd &M);
template void herm_matrix<double>::set_mat<Eigen::MatrixXcd>(int i, Eigen::MatrixXcd &M);

}  // namespace cntr
