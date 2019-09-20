#ifndef CNTR_FUNCTION_EXTERN_TEMPLATES_H
#define CNTR_FUNCTION_EXTERN_TEMPLATES_H

#include <Eigen/Core>
#include "cntr_function_decl.hpp"

#define CFUNC cntr::function<double>

namespace cntr {

  extern template class function<double>;

  extern template void function<double>::set_value<Eigen::MatrixXcd>(int tstp,Eigen::MatrixXcd &M);
  extern template void function<double>::set_value<Eigen::MatrixXcd>(int tstp,cplx x);
  extern template void function<double>::get_value<Eigen::MatrixXcd>(int tstp,Eigen::MatrixXcd &M) const;
  extern template void function<double>::set_matrixelement<Eigen::MatrixXcd>(int tstp,int i1,int i2,Eigen::MatrixXcd &M,int j1,int j2);
  //extern template void function<double>::set_matrixelement(int i1,int i2,function<double> &g,int j1,int j2);

}  // namespace cntr

#endif  // CNTR_FUNCTION_EXTERN_TEMPLATES_H
