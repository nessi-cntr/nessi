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

}  // namespace cntr

#endif  // _CNTR_EIGEN_MAP_H_
