#include "cntr_herm_pseudo_extern_templates.hpp"
#include "cntr_herm_pseudo_impl.hpp"

namespace cntr {

template class herm_pseudo<double>;

template void herm_pseudo<double>::get_les<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M);
template void herm_pseudo<double>::get_gtr<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M);
template void herm_pseudo<double>::get_ret<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M);
template void herm_pseudo<double>::get_vt<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M);
template void herm_pseudo<double>::get_tv<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M);
template void herm_pseudo<double>::get_mat<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M);
template void herm_pseudo<double>::get_matminus<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M);
template void herm_pseudo<double>::set_les<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M);
template void herm_pseudo<double>::set_ret<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M);
template void herm_pseudo<double>::set_tv<Eigen::MatrixXcd>(int i,int j,Eigen::MatrixXcd &M);
template void herm_pseudo<double>::set_mat<Eigen::MatrixXcd>(int i,Eigen::MatrixXcd &M);
template void herm_pseudo<double>::density_matrix<Eigen::MatrixXcd>(int tstp,Eigen::MatrixXcd &M);

}  // namespace cntr
