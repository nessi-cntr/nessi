#ifndef _CNTR_EIGEN_MAP_H_
#define _CNTR_EIGEN_MAP_H_

#include <eigen3/Eigen/Core>

namespace cntr {

// All the data in libcntr is stored in row-major order in unaligned, dense
// memory.
// Set Size1, Size2 to Eigen::Dynamic if the size is not known at
// compile-time.
/// @private
template <int Size1, int Size2, typename Scalar>
struct ElementMap {
    typedef Eigen::OuterStride<Size1> Stride;
    typedef Eigen::Map<Eigen::Matrix<Scalar, Size1, Size2, Eigen::RowMajor>,
                       Eigen::Unaligned, Stride> Map;
    typedef Map type;
    static Map map(int size1, int size2, Scalar *data) {
        return Map(data, size1, size2, Stride(size1));
    }
};

/// @private
template <int Size1, int Size2, typename Scalar>
typename ElementMap<Size1, Size2, Scalar>::type
element_map(int size1, int size2, Scalar *data) {
    return ElementMap<Size1, Size2, Scalar>::map(size1, size2, data);
}

/// @private
template <int Size1, typename Scalar>
typename ElementMap<Size1, Size1, Scalar>::type element_map(int size1,
                                                            Scalar *data) {
    return ElementMap<Size1, Size1, Scalar>::map(size1, size1, data);
}

/// @private  
/** \brief <b> Maps a raw pointer to a complex eigen3 matrix. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* Maps a raw pointer to a complex eigen3 matrix. Involves a copy.
* <!-- ARGUMENTS
*      ========= -->
*
* @param size1
* > Rows of the matrix
* @param size2
* > Cols of the matrix
* @param ptr
* > Pointer storing the matrix.
* @param M
* > eigen3 matrix (output)
*/
template <typename T>
void map_ptr2matrix(int size1, int size2, std::complex<T> *ptr, cdmatrix &M){
  assert(size1 == size2);
  // M.resize(size1, size2);
    // for(int i=0;i<size1;i++) for(int j=0;j<size1;j++) M(i,j)=ptr[i*size1+j];

  switch (size1){
    case 1:
      M = element_map<1, 1>(size1, size2, ptr);
      break;
    case 2:
      M = element_map<2, 2>(size1, size2, ptr);
      break;
    case 3:
      M = element_map<3, 3>(size1, size2, ptr);
      break;
    case 4:
      M = element_map<4, 4>(size1, size2, ptr);
      break;
    case 5:
      M = element_map<5, 5>(size1, size2, ptr);
      break;
    case 6:
      M = element_map<6, 6>(size1, size2, ptr);
      break;
    case 7:
      M = element_map<7, 7>(size1, size2, ptr);
      break;
    case 8:
      M = element_map<8, 8>(size1, size2, ptr);
      break;
    default:
      M = cntr::element_map<-1, -1>(size1, size2, ptr);
  }
}


}  // namespace cntr

#endif  // _CNTR_EIGEN_MAP_H_
