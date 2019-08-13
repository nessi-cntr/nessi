
#include "cntr_function_extern_templates.hpp"
#include "cntr_function_impl.hpp"

namespace cntr {

  template class function<double>;

  template void function<double>::set_value<Eigen::MatrixXcd>(int tstp,Eigen::MatrixXcd &M);
  template void function<double>::set_value<Eigen::MatrixXcd>(int tstp,cplx x);
  template void function<double>::get_value<Eigen::MatrixXcd>(int tstp,Eigen::MatrixXcd &M) const;
  template void function<double>::set_matrixelement<Eigen::MatrixXcd>(int tstp,int i1,int i2,Eigen::MatrixXcd &M,int j1,int j2);
  //template void function<double>::set_matrixelement(int i1,int i2,function<double> &g,int j1,int j2);

}  // namespace cntr
