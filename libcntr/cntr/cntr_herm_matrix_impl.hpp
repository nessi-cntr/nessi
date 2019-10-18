#ifndef CNTR_HERM_MATRIX_IMPL_H
#define CNTR_HERM_MATRIX_IMPL_H

#include "cntr_herm_matrix_decl.hpp"
//#include "cntr_exception.hpp"
#include "cntr_elements.hpp"
#include "cntr_function_decl.hpp"
#include "cntr_herm_matrix_timestep_decl.hpp"
#include "cntr_herm_matrix_timestep_view_impl.hpp"

namespace cntr {

/* #######################################################################################
#
#   CONSTRUCTION/DESTRUCTION
#
########################################################################################*/
template <typename T>
herm_matrix<T>::herm_matrix() {
    les_ = 0;
    tv_ = 0;
    ret_ = 0;
    mat_ = 0;
    ntau_ = 0;
    nt_ = 0;
    size1_ = 0;
    size2_ = 0;
    element_size_ = 0;
    sig_ = -1;
}
template <typename T>
herm_matrix<T>::~herm_matrix() {
    delete[] les_;
    delete[] ret_;
    delete[] tv_;
    delete[] mat_;
}

/** \brief <b> Initializes the `herm_matrix` class for a square-matrix two-time contour function.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Initializes the `herm_matrix` class for a square-matrix two-time contour function.
* If `nt = 0`, memory will be allocated for the Matsubara component only.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nt
* > Number of time steps
* @param ntau
* > Number of points on Matsubara axis
* @param size1
* > Matrix rank of the contour function
* @param sig
* > Set `sig = -1` for fermions or `sig = +1` for bosons.
*/
template <typename T>
herm_matrix<T>::herm_matrix(int nt, int ntau, int size1, int sig) {
    assert(size1 >= 0 && nt >= -1 && sig * sig == 1 && ntau >= 0);
    nt_ = nt;
    ntau_ = ntau;
    sig_ = sig;
    size1_ = size1;
    size2_ = size1;
    element_size_ = size1 * size1;
    if (size1 > 0) {
        mat_ = new cplx[(ntau_ + 1) * element_size_];
    } else {
        mat_ = 0;
    }
    if (nt >= 0 && size1 > 0) {
        les_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        ret_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        tv_ = new cplx[(nt_ + 1) * (ntau_ + 1) * element_size_];
    } else {
        les_ = 0;
        tv_ = 0;
        ret_ = 0;
    }
}

/** \brief <b> Initializes the `herm_matrix` class for a general matrix two-time contour function.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Initializes the `herm_matrix` class for a general (square or non-square)
* matrix two-time contour function.
* If `nt = 0`, memory will be allocated for the Matsubara component only.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nt
* > Number of time steps
* @param ntau
* > Number of points on Matsubara axis
* @param size1
* > Number of matrix rows
* @param size2
* > Number of matrix columns
* @param sig
* > Set `sig = -1` for fermions or `sig = +1` for bosons.
*/
template <typename T> herm_matrix<T>::herm_matrix(int nt,int ntau,int size1,int size2,int sig){
   assert(size1>=0 && size2>=0 && nt>=-1 && sig*sig==1 && ntau>=0);
   nt_=nt;
   ntau_=ntau;
   sig_=sig;
   size1_=size1;
   size2_=size2;
   element_size_=size1*size2;
   if(size1>0){
	   mat_ = new cplx [(ntau_+1)*element_size_];
	   memset(mat_, 0, sizeof(cplx)*(ntau_+1)*element_size_);
   }else{
	   mat_=0;
   }
   if(nt>=0 && size1>0){
	   les_ = new cplx [((nt_+1)*(nt_+2))/2*element_size_];
	   ret_ = new cplx [((nt_+1)*(nt_+2))/2*element_size_];
	   tv_ = new cplx [(nt_+1)*(ntau_+1)*element_size_];
	   memset(les_, 0, sizeof(cplx)*((nt_+1)*(nt_+2))/2*element_size_);
	   memset(ret_, 0, sizeof(cplx)*((nt_+1)*(nt_+2))/2*element_size_);
	   memset(tv_, 0, sizeof(cplx)*(nt_+1)*(ntau_+1)*element_size_);
   }else{
	   les_=0;
	   tv_=0;
	   ret_=0;
   }
}

/** \brief <b> Initializes the `herm_matrix` class with the same layout as a given `herm_matrix`.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Initializes the `herm_matrix` class with the same number of time steps `nt`,
* number of points on the imaginary branch `ntau`, matrix rank `size1` and
* bosonic/fermionic symmetry `sig`. Works for scalar or square-matrix contour objects
* only.
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix` according to which the class should be initialized
*/
template <typename T>
herm_matrix<T>::herm_matrix(const herm_matrix &g) {
    nt_ = g.nt_;
    ntau_ = g.ntau_;
    sig_ = g.sig_;
    size1_ = g.size1_;
    size2_ = g.size1_;
    element_size_ = size1_ * size1_;
    if (size1_ > 0) {
        mat_ = new cplx[(ntau_ + 1) * element_size_];
        memcpy(mat_, g.mat_, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        mat_ = 0;
    }
    if (nt_ >= 0 && size1_ > 0) {
        les_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        ret_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        tv_ = new cplx[(nt_ + 1) * (ntau_ + 1) * element_size_];
        memcpy(les_, g.les_,
               sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
        memcpy(ret_, g.ret_,
               sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
        memcpy(tv_, g.tv_,
               sizeof(cplx) * (nt_ + 1) * (ntau_ + 1) * element_size_);
    } else {
        les_ = 0;
        ret_ = 0;
        tv_ = 0;
    }
}
template <typename T>
herm_matrix<T> &herm_matrix<T>::operator=(const herm_matrix &g) {
    if (this == &g)
        return *this;
    sig_ = g.sig_;
    if (nt_ != g.nt_ || ntau_ != g.ntau_ || size1_ != g.size1_) {
        delete[] les_;
        delete[] ret_;
        delete[] tv_;
        delete[] mat_;
        nt_ = g.nt_;
        ntau_ = g.ntau_;
        size1_ = g.size1_;
        size2_ = g.size1_;
        element_size_ = size1_ * size1_;
        if (size1_ > 0) {
            mat_ = new cplx[(ntau_ + 1) * element_size_];
        } else {
            mat_ = 0;
        }
        if (size1_ > 0 && nt_ >= 0) {
            les_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
            ret_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
            tv_ = new cplx[(nt_ + 1) * (ntau_ + 1) * element_size_];
        } else {
            les_ = 0;
            ret_ = 0;
            tv_ = 0;
        }
    }
    if (size1_ > 0) {
        memcpy(mat_, g.mat_, sizeof(cplx) * (ntau_ + 1) * element_size_);
        if (nt_ >= 0) {
            memcpy(les_, g.les_, sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 *
                                     element_size_);
            memcpy(ret_, g.ret_, sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 *
                                     element_size_);
            memcpy(tv_, g.tv_,
                   sizeof(cplx) * (nt_ + 1) * (ntau_ + 1) * element_size_);
        }
    }
    return *this;
}
#if __cplusplus >= 201103L
template <typename T>
herm_matrix<T>::herm_matrix(herm_matrix &&g) noexcept
    : les_(g.les_),
      ret_(g.ret_),
      tv_(g.tv_),
      mat_(g.mat_),
      nt_(g.nt_),
      ntau_(g.ntau_),
      size1_(g.size1_),
      size2_(g.size2_),
      element_size_(g.element_size_),
      sig_(g.sig_) {
    g.les_ = nullptr;
    g.tv_ = nullptr;
    g.ret_ = nullptr;
    g.mat_ = nullptr;
    g.ntau_ = 0;
    g.nt_ = 0;
    g.size1_ = 0;
    g.size2_ = 0;
    g.element_size_ = 0;
}
template <typename T>
herm_matrix<T> &herm_matrix<T>::operator=(herm_matrix &&g) noexcept {
    if (&g == this)
        return *this;

    les_ = g.les_;
    ret_ = g.ret_;
    tv_ = g.tv_;
    mat_ = g.mat_;
    nt_ = g.nt_;
    ntau_ = g.ntau_;
    size1_ = g.size1_;
    size2_ = g.size2_;
    element_size_ = g.element_size_;
    sig_ = g.sig_;

    g.les_ = nullptr;
    g.tv_ = nullptr;
    g.ret_ = nullptr;
    g.mat_ = nullptr;
    g.ntau_ = 0;
    g.nt_ = 0;
    g.size1_ = 0;
    g.size2_ = 0;
    g.element_size_ = 0;

    return *this;
}
#endif
/* #######################################################################################
#
#   RESIZE
#
########################################################################################*/

/** \brief <b> Discards the `herm_matrix` object and resize with respect to the number of
 * time points `nt`, points on the Matsubara branch `ntau` and the square-matrix size
 * `size1`.  </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 *
 * > Discard `herm_matrix` class and then resize with respect to number of time steps `nt`,
 * > points on the Matsubara branch `ntau` and the maxtix size `size1`.
 * > If `nt = 0`, real-time components are deallocated (i.e. only the Matsubara
 * > component is kept in memory). Internal routine; see top-level interface `resize`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param nt
 * > New number of time steps.
 * @param ntau
 * > New number of points on Matsubara axis
 * @param size1
 * > New matrix rank of the contour function
 */
template <typename T>
void herm_matrix<T>::resize_discard(int nt, int ntau, int size1) {
    assert(ntau >= 0 && nt >= -1 && size1 >= 0);
    delete[] les_;
    delete[] ret_;
    delete[] tv_;
    delete[] mat_;
    nt_ = nt;
    ntau_ = ntau;
    size1_ = size1;
    size2_ = size1;
    element_size_ = size1 * size1;
    if (size1 > 0) {
        mat_ = new cplx[(ntau_ + 1) * element_size_];
    } else {
        mat_ = 0;
    }
    if (nt_ >= 0 && size1_ > 0) {
        les_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        ret_ = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        tv_ = new cplx[(nt_ + 1) * (ntau_ + 1) * element_size_];
    } else {
        les_ = 0;
        tv_ = 0;
        ret_ = 0;
    }
}

/** \brief <b> Resizes `herm_matrix` object with respect to the number
 * of time points `nt`. Internal routine; see top-level interface `resize`.  </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 *
 * > Resizes `herm_matrix` class with respect to number of time steps `nt`. If `nt >= 0`
 * > real-time components are resized, otherwise deallocated (i.e. only the Matsubara
 * > component is kept in memory).
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param nt
 * > New number of time steps.
 *
 */
template <typename T>
void herm_matrix<T>::resize_nt(int nt) {
    int nt1 = (nt_ > nt ? nt : nt_);
    cplx *ret, *les, *tv;
    assert(nt >= -1);
    nt_ = nt;
    if (size1_ == 0)
        return;
    if (nt_ >= 0) {
        les = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        ret = new cplx[((nt_ + 1) * (nt_ + 2)) / 2 * element_size_];
        tv = new cplx[(nt_ + 1) * (ntau_ + 1) * element_size_];
        if (nt1 >= 0) {
            memcpy(les, les_, sizeof(cplx) * ((nt1 + 1) * (nt1 + 2)) / 2 *
                                  element_size_);
            memcpy(ret, ret_, sizeof(cplx) * ((nt1 + 1) * (nt1 + 2)) / 2 *
                                  element_size_);
            memcpy(tv, tv_,
                   sizeof(cplx) * (nt1 + 1) * (ntau_ + 1) * element_size_);
        }
    } else {
        les = 0;
        ret = 0;
        tv = 0;
    }
    delete[] les_;
    delete[] ret_;
    delete[] tv_;
    les_ = les;
    ret_ = ret;
    tv_ = tv;
}
/** \brief <b> Resizes `herm_matrix` object with respect to the number of
 * time points `nt`, points on the Matsubara branch `ntau` or the matrix size
 * `size1`.  </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 *
 * > Resizes `herm_matrix` class with respect to number of time steps `nt`. If `nt >= 0`
 * > real-time components are resized, otherwise deallocated (i.e. only the Matsubara
 * > component is kept in memory).
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param nt
 * > New number of time steps.
 * @param ntau
 * > Number of points on Matsubara axis
 * @param size1
 * > size of the matrix
 */
template <typename T>
void herm_matrix<T>::resize(int nt, int ntau, int size1) {
    // std::cout  << "herm_matrix<T>::" << __FUNCTION__ << " " << nt << " " <<
    // ntau << " " << size1 << std::endl;
    if (ntau == ntau_ && size1_ == size1)
        resize_nt(nt);
    else
        resize_discard(nt, ntau, size1);
}
/** \brief <b> Sets all values to zero. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 *
 * > Sets all values of all components to zero.
 *
 */
template <typename T>
void herm_matrix<T>::clear(void) {
    if (size1_ == 0)
        return;
    memset(mat_, 0, sizeof(cplx) * (ntau_ + 1) * element_size_);
    if (nt_ >= 0) {
        memset(les_, 0,
               sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
        memset(ret_, 0,
               sizeof(cplx) * ((nt_ + 1) * (nt_ + 2)) / 2 * element_size_);
        memset(tv_, 0,
               sizeof(cplx) * (nt_ + 1) * (ntau_ + 1) * element_size_);
    }
}
/* #######################################################################################
#
#   RAW POINTERS TO ELEMENTS
#
########################################################################################*/
/// @private
template <typename T>
int herm_matrix<T>::les_offset(int t, int t1) const {
    assert(t >= 0 && t1 >= 0 && t <= t1 && t1 <= nt_);
    return ((t1 * (t1 + 1)) / 2 + t) * element_size_;
}
/// @private
template <typename T>
int herm_matrix<T>::ret_offset(int t, int t1) const {
    assert(t >= 0 && t1 >= 0 && t <= nt_ && t1 <= t);
    return ((t * (t + 1)) / 2 + t1) * element_size_;
}
/// @private
template <typename T>
int herm_matrix<T>::tv_offset(int t, int tau) const {
    assert(t >= 0 && tau >= 0 && t <= nt_ && tau <= ntau_);
    return (t * (ntau_ + 1) + tau) * element_size_;
}
/// @private
template <typename T>
int herm_matrix<T>::mat_offset(int tau) const {
    assert(tau >= 0 && tau <= ntau_);
    return tau * element_size_;
}
/// @private
template <typename T>
inline std::complex<T> *herm_matrix<T>::lesptr(int t, int t1) {
    return les_ + les_offset(t, t1);
}
/// @private
template <typename T>
inline std::complex<T> *herm_matrix<T>::retptr(int t, int t1) {
    return ret_ + ret_offset(t, t1);
}
/// @private
template <typename T>
inline std::complex<T> *herm_matrix<T>::tvptr(int t, int tau) {
    return tv_ + tv_offset(t, tau);
}
/// @private
template <typename T>
inline std::complex<T> *herm_matrix<T>::matptr(int tau) {
    return mat_ + mat_offset(tau);
}
/// @private
template <typename T>
inline const std::complex<T> *herm_matrix<T>::lesptr(int t, int t1) const {
    return les_ + les_offset(t, t1);
}
/// @private
template <typename T>
inline const std::complex<T> *herm_matrix<T>::retptr(int t, int t1) const {
    return ret_ + ret_offset(t, t1);
}
/// @private
template <typename T>
inline const std::complex<T> *herm_matrix<T>::tvptr(int t, int tau) const {
    return tv_ + tv_offset(t, tau);
}
/// @private
template <typename T>
inline const std::complex<T> *herm_matrix<T>::matptr(int tau) const {
    return mat_ + mat_offset(tau);
}
/* #######################################################################################
#
#   READING ELEMENTS TO ANY MATRIX TYPE
#   OR TO COMPLEX NUMBERS (then only the (0,0) element is addressed for dim>0)
#   note: these are not efficient, in particular the conjugation
#
########################################################################################*/
#define herm_matrix_READ_ELEMENT                                             \
    {                                                                        \
        int r, s;                                              \
        M.resize(size1_, size2_);                                                  \
        for (r = 0; r < size1_; r++)                                            \
            for (s = 0; s < size2_; s++)                                        \
                M(r, s) = x[r * size2_ + s];                                    \
    }
// TODO: Please think about this function, what happens for rectangular
#define herm_matrix_READ_ELEMENT_MINUS_CONJ                                  \
    {                                                                        \
        cplx w;                                                              \
        int r, s, dim = size1_;                                              \
        M.resize(dim, dim);                                                  \
        for (r = 0; r < dim; r++)                                            \
            for (s = 0; s < dim; s++) {                                      \
                w = x[s * dim + r];                                          \
                M(r, s) = std::complex<T>(-w.real(), w.imag());              \
            }                                                                \
    }
/// @private
/** \brief <b> Returns the lesser component at given times.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the lesser component \f$ C^>(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_i\f$ to a given matrix class. Hermitian symmetry is used for the
* > case \f$ t_i > t_j \f$ to express \f$ C^>(t_i,t_j)  = - [C^>(t_j,t_i)]^\ddagger\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param j
* > Index of time \f$ t_j\f$ .
* @param M
* > Matrix to which the lesser component is given.
*/
  template <typename T> template <class Matrix>
  void herm_matrix<T>::get_les(int i, int j, Matrix &M) const {
    assert(i <= nt_ && j <= nt_);
    const cplx *x;
    if (i <= j) {
      x = lesptr(i, j);
      herm_matrix_READ_ELEMENT
	} else {
      x = lesptr(j, i);
      herm_matrix_READ_ELEMENT_MINUS_CONJ
	}
  }
/// @private
/** \brief <b> Returns the retarded component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the retarded component \f$ C^\mathrm{R}(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_j\f$ to a given matrix class. \f$ C^\mathrm{R}(t_i,t_j) =0 \f$ for
* > \f$ t_i < t_j \f$; however, \f$ - [C^\mathrm{R}(t_j,t_i)]^\ddagger\f$ is
* > returned in this case.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param j
* > Index of time \f$ t_j\f$ .
* @param M
* > Matrix to which the retarded component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix<T>::get_ret(int i, int j, Matrix &M) const {
  assert(i <= nt_ && j <= nt_);
    const cplx *x;
    if (i >= j) {
        x = retptr(i, j);
        herm_matrix_READ_ELEMENT
    } else {
        x = retptr(j, i);
        herm_matrix_READ_ELEMENT_MINUS_CONJ
    }
}
/// @private
/** \brief <b> Returns the left-mixing component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the left-mixing component \f$ C^\rceil(t_i,\tau_j) \f$ at given times \f$ t_i\f$
* > and imaginary time \f$ \tau_j\f$ to a given matrix class.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param j
* > Index of time \f$ \tau_j\f$ .
* @param M
* > Matrix to which the left-mixing component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix<T>::get_tv(int i, int j, Matrix &M) const {
  assert(i <= nt_ && j <= ntau_);
    const cplx *x = tvptr(i, j);
    herm_matrix_READ_ELEMENT
}
/// @private
/** \brief <b> Returns the right-mixing component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the right-mixing component \f$ C^\lceil(\tau_i,t_j) \f$ at given
* > imaginary \f$ \tau_i\f$ and real \f$ t_j\f$ time. The right-mixing
* > component is expressed by the left-mixing component.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of imaginary time \f$ \tau_i\f$ .
* @param j
* > Index of time \f$ t_j\f$ .
* @param M
* > Matrix to which the right-mixing component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix<T>::get_vt(int i, int j, Matrix &M) const {
  assert(i <= ntau_ && j <= nt_);
    const cplx *x = tvptr(j, ntau_ - i);
    herm_matrix_READ_ELEMENT_MINUS_CONJ if (sig_ == -1) M = -M;
}
/// @private
/** \brief <b> Returns the Matsubara component at given imaginary time. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the Matsubara component \f$ C^\mathrm{M}(\tau_i) \f$ at given
* > imaginary time \f$ \tau_i\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of imaginary time \f$ \tau_i\f$ .
* @param M
* > Matrix to which the Matsubara component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix<T>::get_mat(int i, Matrix &M) const {
  assert(i <= ntau_);
    const cplx *x = matptr(i);
    herm_matrix_READ_ELEMENT
}
/// @private
/** \brief <b> Returns the Matsubara component for the negative
* of a given imaginary time. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the Matsubara component \f$ C^\mathrm{M}(-\tau_i) \f$ at given
* > imaginary time \f$ \tau_i\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] Index of imaginary time \f$ \tau_i\f$ .
* @param M
* > [Matrix] Matrix to which the Matsubara component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix<T>::get_matminus(int i, Matrix &M) const {
  assert(i <= ntau_);
    const cplx *x = matptr(ntau_ - i);
    herm_matrix_READ_ELEMENT if (sig_ == -1) M = -M;
}
#undef herm_matrix_READ_ELEMENT
#undef herm_matrix_READ_ELEMENT_MINUS_CONJ

/// @private
/** \brief <b> Returns the greater component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the greater component \f$ C^>(t_i,t_j) \f$ at given
* > times \f$t_i\f$ \f$t_j\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param j
* > Index of time \f$ t_j\f$ .
* @param M
* > Matrix to which the Matsubara component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix<T>::get_gtr(int i, int j, Matrix &M) const {
  assert(i <= nt_ && j <= nt_);
    Matrix M1;
    get_ret(i, j, M);
    get_les(i, j, M1);
    M += M1;
}
/// @private
/** \brief <b> Returns the retarded component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the scalar-valued retarded component of a scaler type \f$ C^\mathrm{R}(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_i\f$. \f$ C^\mathrm{R}(t_j,t_i) =0 \f$ for
* \f$ t_j > t_i \f$; however, \f$ - [C^\mathrm{R}(t_j,t_i)]^*\f$ is
* returned in this case.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param j
* > Index of time \f$ t_j\f$ .
* @param x
* > The value of the retarded component.
*/
template <typename T>
inline void herm_matrix<T>::get_ret(int i, int j, cplx &x) const {
  assert(i <= nt_ && j <= nt_);
    if (i >= j)
        x = *retptr(i, j);
    else {
        x = *retptr(j, i);
        x = -std::conj(x);
    }
}
/// @private
/** \brief <b> Returns the lesser component at given times. </b>
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the scalar-valued lesser component of a scaler type \f$ C^>(t_i,t_j) \f$ at given times \f$ t_i\f$
* > and \f$ t_i\f$. Hermitian symmetry is used for the
* > case \f$ t_i > t_j \f$ to express \f$ C^>(t_i,t_j)  = - [C^>(t_j,t_i)]^\ddagger\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param j
* > Index of time \f$ t_j\f$ .
* @param x
* > The value of the lesser component.
*/
template <typename T>
inline void herm_matrix<T>::get_les(int i, int j, cplx &x) const {
  assert(i <= nt_ && j <= nt_);
    if (i <= j)
        x = *lesptr(i, j);
    else {
        x = *lesptr(j, i);
        x = -std::conj(x);
    }
}
/// @private
/** \brief <b> Returns the left-mixing component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the scalar-valued left-mixing component \f$ C^\rceil(t_i,\tau_j) \f$ at given times \f$ t_i\f$
* > and imaginary time \f$ \tau_j\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param j
* > Index of time \f$ \tau_j\f$ .
* @param x
* > The value of the left-mixing component.
*/
template <typename T>
inline void herm_matrix<T>::get_tv(int i, int j, cplx &x) const {
  assert(i <= nt_ && j <= ntau_);
    x = *tvptr(i, j);
}
/// @private
/** \brief <b> Returns the right-mixing component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the scalar-valued right-mixing component \f$ C^\lceil(\tau_i,t_j) \f$ at given
* > imaginary \f$ \tau_i\f$ and real \f$ t_j\f$ time. The right-mixing
* > component is expressed by the left-mixing component.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of imaginary time \f$ \tau_i\f$ .
* @param j
* > Index of time \f$ t_j\f$ .
* @param x
* > The value of the right-mixing component.
*/
template <typename T>
inline void herm_matrix<T>::get_vt(int i, int j, cplx &x) const {
  assert(i <= ntau_ && j <= nt_);
    x = *tvptr(j, ntau_ - i);
    if (sig_ == -1)
        x = std::conj(x);
    else
        x = -std::conj(x);
}
/// @private
/** \brief <b> Returns the Matsubara component at given imaginary time. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the scalar-valued Matsubara component \f$ C^\mathrm{M}(\tau_i) \f$ at given
* > imaginary time \f$ \tau_i\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of imaginary time \f$ \tau_i\f$ .
* @param x
* > The value of the Matsubara component.
*/
template <typename T>
inline void herm_matrix<T>::get_mat(int i, cplx &x) const {
  assert(i <= ntau_);
    x = *matptr(i);
}
/// @private
/** \brief <b> Returns the Matsubara component for the negative
* of a given imaginary time. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the scalar-valued Matsubara component \f$ C^\mathrm{M}(-\tau_i) \f$ at given
* > imaginary time \f$ \tau_i\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of imaginary time \f$ \tau_i\f$ .
* @param x
* > The value of the Matsubara component.
*/
template <typename T>
inline void herm_matrix<T>::get_matminus(int i, cplx &x) const {
  assert(i <= ntau_);
    x = *matptr(ntau_ - i);
    if (sig_ == -1)
        x = -x;
}
/// @private
/** \brief <b> Returns the greater component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the scalar-valued greater component \f$ C^>(t_i,t_j) \f$ at given
* > times \f$t_i\f$ \f$t_j\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param j
* > Index of time \f$ t_j\f$ .
* @param x
* > The value of the greater component.
*/
template <typename T>
inline void herm_matrix<T>::get_gtr(int i, int j, cplx &x) const {
  assert(i <= nt_ && j <= nt_);
    cplx x1;
    get_ret(i, j, x);
    get_les(i, j, x1);
    x += x1;
}
/// @private
/** \brief <b> Returns the density matrix at given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Returns the scalar-valued density matrix (occupation, that is) at
* given time step `tstp`. Setting `tstp = -1` returns the equilibrium density matrix
* \f$ \rho = -C^\mathrm{M}(\beta) \f$, while `tstp >= 0` returns
* \f$ \rho(t) = i \eta C^<(t,t) \f$.
* The return value is formally complex.
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > The time step at which the density matrix is returned.
*/
template <typename T>
std::complex<T> herm_matrix<T>::density_matrix(int tstp) const {
  assert(tstp >= -1 && tstp <= nt_);
    cplx x1;
    if (tstp == -1) {
        get_mat(ntau_, x1);
        return -x1;
    } else {
        get_les(tstp, tstp, x1);
        return std::complex<T>(0.0, sig_) * x1;
    }
}

/** \brief <b> Returns the density matrix at given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Returns the matrix-valued density matrix at
* given time step `tstp`. Setting `tstp = -1` returns the equilibrium density matrix
* \f$ \rho = -C^\mathrm{M}(\beta) \f$, while `tstp >= 0` returns
* \f$ \rho(t) = i \eta C^<(t,t) \f$.
* Works for square-matrices only.
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > The time step at which the density matrix is returned.
* @param M
* > The density matrix at time step `tstp`.
*/
template <typename T>
template <class Matrix>
void herm_matrix<T>::density_matrix(int tstp, Matrix &M) const {
  assert(M.rows() == size1_ && M.cols() == size2_);
  assert(tstp >= -1 && tstp <= nt_);
    if (tstp == -1) {
        get_mat(ntau_, M);
        M *= (-1.0);
    } else {
        get_les(tstp, tstp, M);
        M *= std::complex<T>(0.0, 1.0 * sig_);
    }
}

/* #######################################################################################
#
#   WRITING ELEMENTS FROM ANY MATRIX TYPE
#   OR FROM COMPLEX NUMBERS (then only the (0,0) element is addressed for
dim>0)
#
########################################################################################*/
#define herm_matrix_SET_ELEMENT_MATRIX                                       \
    {                                                                        \
        int r, s;                                                             \
        for (r = 0; r < size1_; r++)                                            \
            for (s = 0; s < size2_; s++)                                        \
                x[r * size2_ + s] = M(r, s);                                    \
    }

/// @private
/** \brief <b> Copies a given matrix to the retarded component of at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* > Copies a given matrix to the retarded component of at given times. Matrix type is assumed.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$.
* @param j
* > Index of time \f$ t_j\f$.
* @param M
* > Matrix which is copied to the retarded component is given
*/
template <typename T>
template <class Matrix>
void herm_matrix<T>::set_ret(int i, int j, Matrix &M) {
  assert(i <= nt_ && j <= nt_);
    cplx *x = retptr(i, j);
    herm_matrix_SET_ELEMENT_MATRIX
}
/// @private
/** \brief <b> Copies a given matrix to the lesser component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*  > Copies a given matrix to the lesser component of at given times. Matrix type is assumed.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$.
* @param j
* > Index of time \f$ t_j\f$.
* @param M
* > Matrix which is copied to the retarded component is given
*/
template <typename T>
template <class Matrix>
void herm_matrix<T>::set_les(int i, int j, Matrix &M) {
  assert(i <= nt_ && j <= nt_);
    cplx *x = lesptr(i, j);
    herm_matrix_SET_ELEMENT_MATRIX
}
/// @private
/** \brief <b> Copies a given matrix to the left-mixing component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*  > Copies a given matrix to the left-mixing component of at given times. Matrix type is assumed.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$.
* @param j
* > Index of time \f$ \tau_j\f$.
* @param M
* > Matrix which is copied to the retarded component is given
*/
template <typename T>
template <class Matrix>
void herm_matrix<T>::set_tv(int i, int j, Matrix &M) {
  assert(i <= nt_ && j <= ntau_);
    cplx *x = tvptr(i, j);
    herm_matrix_SET_ELEMENT_MATRIX
}
/// @private
/** \brief <b> Copies a given matrix to the Matsubara component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*  > Copies a given matrix to the Matsubara component of at given times. Matrix type is assumed.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ \tau_i\f$.
* @param M
* > Matrix which is copied to the retarded component is given
*/
template <typename T>
template <class Matrix>
void herm_matrix<T>::set_mat(int i, Matrix &M) {
  assert(i <= ntau_);
    cplx *x = matptr(i);
    herm_matrix_SET_ELEMENT_MATRIX
}


/** \brief <b> Hermitianizes the Matsubara component at given imaginary time. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* > Forces the Matsubara component at \f$ \tau_i \f$ to be hermitian by
* > transforming \f$ C^\mathrm{M}(\tau_i) \rightarriw [C^\mathrm{M}(\tau_i) + (C^\mathrm{M}(\tau_i))^\ddagger]/2\f$.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of imaginary time \f$ \tau_i\f$.
*/
template<typename T>
void herm_matrix<T>::set_mat_herm(int i){
  cdmatrix tmp(size1_,size2_);
  this->get_mat(i,tmp);
  tmp=(tmp+tmp.adjoint())/2.0;
  this->set_mat(i,tmp);
}
/// @private
/** \brief <b> Hermitianizes the Matsubara component at for all imaginary times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* > Forces the Matsubara component to be hermitian by
* > transforming \f$C^\mathrm{M}(\tau_i) \rightarriw [C^\mathrm{M}(\tau_i) + (C^\mathrm{M}(\tau_i))^\ddagger]/2\f$
* > for all `i=0,...,ntau`.
*/
template<typename T>
void herm_matrix<T>::set_mat_herm(void){
  cdmatrix tmp(size1_,size2_);
  for(int i=0;i<ntau_;i++){
    this->get_mat(i,tmp);
    tmp=(tmp+tmp.adjoint())/2.0;
    this->set_mat(i,tmp);
  }
}
/// @private
/** \brief <b> Copies a given complex number to the lesser component of at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* > Copies a given matrix to the lesser component of at given times. Matrix type is assumed.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$.
* @param j
* > Index of time \f$ t_j\f$.
* @param x
* > Complex number which is copied to the retarded component is given
*/
template <typename T>
inline void herm_matrix<T>::set_les(int i, int j, cplx x) {
  assert(i <= nt_ && j <= nt_);
    *lesptr(i, j) = x;
}

/// @private
/** \brief <b> Copies a given complex number to the retarded component of at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* > Copies a given matrix to the retarded component of at given times. Matrix type is assumed.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$.
* @param j
* > Index of time \f$ t_j\f$.
* @param x
* > Complex number which is copied to the retarded component is given
*/
template <typename T>
inline void herm_matrix<T>::set_ret(int i, int j, cplx x) {
  assert(i <= nt_ && j <= nt_);
    *retptr(i, j) = x;
}
/// @private
/** \brief <b> Copies a given complex number to the left-mixing component of at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* > Copies a given matrix to the left-mixing component of at given times. Matrix type is assumed.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$.
* @param j
* > Index of time \f$ \tau_j\f$.
* @param x
* > Complex number which is copied to the retarded component is given
*/
template <typename T>
inline void herm_matrix<T>::set_tv(int i, int j, cplx x) {
  assert(i <= nt_ && j <= ntau_);
    *tvptr(i, j) = x;
}
/// @private
/** \brief <b> Copies a given complex number to the Matsubara component of at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
* > Copies a given matrix to the Matsubara component of at given times. Matrix type is assumed.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ \tau_i\f$.
* @param x
* > Complex number which is copied to the retarded component is given
*/
template <typename T>
inline void herm_matrix<T>::set_mat(int i, cplx x) {
  assert(i <= ntau_ );
    *matptr(i) = x;
}



/* #######################################################################################
#
#   INPUT/OUTPUT FROM/TO FILES
#
########################################################################################*/

/** \brief <b> Saves the `herm_matrix` to file in text format. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Saves the `herm_matrix` to file in plain text format, which can
 * > be read by the Python tools or by the corresponding `read_from_file`
 * > routine. Warning: this output format generates very large files.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param file
 * > The output file name.
 * @param precision
 * > The number of digits to be stored.
 */
template <typename T>
void herm_matrix<T>::print_to_file(const char *file, int precision) const {
    int i, j, l, sg = element_size_;
    std::ofstream out;
    out.open(file, std::ios::out);
    out.precision(precision);
    out << "# " << nt_ << " " << ntau_ << " " << size1_ << " "
        << " " << sig_ << std::endl;
    for (j = 0; j <= ntau_; j++) {
        out << "mat: " << j;
        for (l = 0; l < sg; l++)
            out << " " << matptr(j)[l].real() << " " << matptr(j)[l].imag();
        out << std::endl;
    }
    out << std::endl;
    if (nt_ >= 0) {
        for (i = 0; i <= nt_; i++) {
            for (j = 0; j <= i; j++) {
                out << "ret: " << i << " " << j;
                for (l = 0; l < sg; l++)
                    out << " " << retptr(i, j)[l].real() << " "
                        << retptr(i, j)[l].imag();
                out << std::endl;
            }
            out << std::endl;
        }
        out << std::endl;
        for (i = 0; i <= nt_; i++) {
            for (j = 0; j <= ntau_; j++) {
                out << "tv: " << i << " " << j;
                for (l = 0; l < sg; l++)
                    out << " " << tvptr(i, j)[l].real() << " "
                        << tvptr(i, j)[l].imag();
                out << std::endl;
            }
            out << std::endl;
        }
        out << std::endl;
        for (j = 0; j <= nt_; j++) {
            for (i = 0; i <= j; i++) {
                out << "les: " << i << " " << j;
                for (l = 0; l < sg; l++)
                    out << " " << lesptr(i, j)[l].real() << " "
                        << lesptr(i, j)[l].imag();
                out << std::endl;
            }
            out << std::endl;
        }
        out << std::endl;
    }
    out.close();
}
/** \brief <b> Reads the `herm_matrix` from file in text format. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Reads the `herm_matrix` to file in plain text format, which happens
 * > been created with the corresponding routine `herm_matrix::print_to_file`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param file
 * > The input file name.
 */
template <typename T>
void herm_matrix<T>::read_from_file(const char *file) {
    int i, n, m, j, l, size1, sg, sig;
    double real, imag;
    std::string s;
    std::ifstream out;
    out.open(file, std::ios::in);
    if (!(out >> s >> n >> m >> size1 >> sig)) {
        std::cerr << "read G from file " << file << " error in file"
                  << std::endl;
        abort();
    }
    if (n > nt_ || m != ntau_ || size1 != size1_)
        resize(n, m, size1);
    sig_ = sig;
    sg = element_size_;
    for (j = 0; j <= ntau_; j++) {
        out >> s >> s;
        for (l = 0; l < sg; l++) {
            if (!(out >> real >> imag)) {
                std::cerr << "read G from file " << file << " error at mat ("
                          << j << ")" << std::endl;
                abort();
            }
            matptr(j)[l] = std::complex<T>(real, imag);
        }
    }
    if (n >= 0) {
        for (i = 0; i <= n; i++) {
            for (j = 0; j <= i; j++) {
                out >> s >> s >> s;
                for (l = 0; l < sg; l++) {
                    if (!(out >> real >> imag)) {
                        std::cerr << "read G from file " << file
                                  << " error at ret (" << i << "," << j << ")"
                                  << std::endl;
                        abort();
                    }
                    retptr(i, j)[l] = std::complex<T>(real, imag);
                }
            }
        }
        for (i = 0; i <= n; i++) {
            for (j = 0; j <= ntau_; j++) {
                out >> s >> s >> s;
                for (l = 0; l < sg; l++) {
                    if (!(out >> real >> imag)) {
                        std::cerr << "read G from file " << file
                                  << " error at tv (" << i << "," << j << ")"
                                  << std::endl;
                        abort();
                    }
                    tvptr(i, j)[l] = std::complex<T>(real, imag);
                }
            }
        }
        for (j = 0; j <= n; j++) {
            for (i = 0; i <= j; i++) {
                out >> s >> s >> s;
                for (l = 0; l < sg; l++) {
                    if (!(out >> real >> imag)) {
                        std::cerr << "read G from file " << file
                                  << " error at les (" << i << "," << j << ")"
                                  << std::endl;
                        abort();
                    }
                    lesptr(i, j)[l] = std::complex<T>(real, imag);
                }
            }
        }
    }
    out.close();
}
/** \brief <b> Reads the `herm_matrix` up to a given number of time steps
 *  from file in text format. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Reads the `herm_matrix` to file in plain text format, which happens
 * > been created with the corresponding routine `herm_matrix:print_to_file`.
 * > The `herm_matrix` is read up to a given number of time steps `nt1`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param nt1
 * > Read up to this time step.
 * @param file
 * > The input file name.
 */
template <typename T>
void herm_matrix<T>::read_from_file(int nt1, const char *file) {
    int i, n, m, j, l, size1, sg, sig;
    double real, imag;
    std::string s;
    std::ifstream out;
    out.open(file, std::ios::in);

    if (!(out >> s >> n >> m >> size1 >> sig)) {
        std::cerr << "read G from file " << file << " error in file"
                  << std::endl;
        abort();
    }

    assert(nt1 <= nt_ && "nt1 > nt_");
    assert(nt1 <= n && "nt1 > n");
    assert(m == ntau_ && "m /= ntau_");
    assert(size1 == size1_ && "size1 /= size1_");

    sg = size1 * size1;
    if (nt1 >= -1) {
        for (j = 0; j <= ntau_; j++) {
            out >> s >> s;
            for (l = 0; l < sg; l++) {
                if (!(out >> real >> imag)) {
                    std::cerr << "read G from file " << file
                              << " error at mat (" << j << ")" << std::endl;
                    abort();
                }
                matptr(j)[l] = std::complex<T>(real, imag);
            }
        }
    }
    if (n >= 0) {
        for (i = 0; i <= n; i++) {
            for (j = 0; j <= i; j++) {
                out >> s >> s >> s;
                for (l = 0; l < sg; l++) {
                    if (!(out >> real >> imag)) {
                        std::cerr << "read G from file " << file
                                  << " error at ret (" << i << "," << j << ")"
                                  << std::endl;
                        abort();
                    }
                    if (i <= nt1)
                        retptr(i, j)[l] = std::complex<T>(real, imag);
                }
            }
        }
        for (i = 0; i <= n; i++) {
            for (j = 0; j <= ntau_; j++) {
                out >> s >> s >> s;
                for (l = 0; l < sg; l++) {
                    if (!(out >> real >> imag)) {
                        std::cerr << "read G from file " << file
                                  << " error at tv (" << i << "," << j << ")"
                                  << std::endl;
                        abort();
                    }
                    if (i <= nt1)
                        tvptr(i, j)[l] = std::complex<T>(real, imag);
                }
            }
        }
        for (j = 0; j <= n; j++) {
            for (i = 0; i <= j; i++) {
                out >> s >> s >> s;
                for (l = 0; l < sg; l++) {
                    if (!(out >> real >> imag)) {
                        std::cerr << "read G from file " << file
                                  << " error at les (" << i << "," << j << ")"
                                  << std::endl;
                        abort();
                    }
                    if (j <= nt1)
                        lesptr(i, j)[l] = std::complex<T>(real, imag);
                }
            }
        }
    }
    out.close();
}

#if CNTR_USE_HDF5 == 1

/** \brief <b> Stores `herm_matrix` to a given group in HDF5 format. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Stores the `herm_matrix`, including the matrix size, number of
 * > points on real and imaginary axis and fermionic/bosonic character,
 * > to a given HDF5 group.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param group_id
 * > The HDF5 group handle under which the `herm_matrix` is stored.
 */
template <typename T>
void herm_matrix<T>::write_to_hdf5(hid_t group_id) {
    store_int_attribute_to_hid(group_id, std::string("ntau"), ntau_);
    store_int_attribute_to_hid(group_id, std::string("nt"), nt_);
    store_int_attribute_to_hid(group_id, std::string("sig"), sig_);
    store_int_attribute_to_hid(group_id, std::string("size1"), size1_);
    store_int_attribute_to_hid(group_id, std::string("size2"), size2_);
    store_int_attribute_to_hid(group_id, std::string("element_size"),
                               element_size_);
    hsize_t len_shape = 3, shape[3];
    shape[1] = size1_;
    shape[2] = size2_;
    if (nt_ > -2) {
        shape[0] = ntau_ + 1;
        store_cplx_array_to_hid(group_id, std::string("mat"), matptr(0),
                                shape, len_shape);
    }
    if (nt_ > -1) {
        shape[0] = ((nt_ + 1) * (nt_ + 2)) / 2;
        // CHECK: implement store_cplx_array_to_hid with template typename T
        store_cplx_array_to_hid(group_id, std::string("ret"), retptr(0, 0),
                                shape, len_shape);
        store_cplx_array_to_hid(group_id, std::string("les"), lesptr(0, 0),
                                shape, len_shape);
        shape[0] = (nt_ + 1) * (ntau_ + 1);
        store_cplx_array_to_hid(group_id, std::string("tv"), tvptr(0, 0),
                                shape, len_shape);
    }
}
/** \brief <b> Stores `herm_matrix` to a given group in HDF5 format. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Stores the `herm_matrix`, including the matrix size, number of
 * > points on real and imaginary axis and fermionic/bosonic character,
 * > to a given HDF5 group with given groupname.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param group_id
 * > [hid_t] The HDF5 group handle under which the `herm_matrix` is stored.
 * @param groupname
 * > [char*] The name of the HDF5 group.
 */
template <typename T>
void herm_matrix<T>::write_to_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = create_group(group_id, groupname);
    this->write_to_hdf5(sub_group_id);
    close_group(sub_group_id);
}
/** \brief <b> Stores `herm_matrix` to a given file in HDF5 format. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Stores the `herm_matrix`, including the matrix size, number of
 * > points on real and imaginary axis and fermionic/bosonic character,
 * > to a given HDF5 file with given groupname.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param filename
 * > [char*] The name of the file to which the `herm_matrix` is stored.
 * @param groupname
 * > [char*] The name of the HDF5 group.
 */
template <typename T>
void herm_matrix<T>::write_to_hdf5(const char *filename,
                                   const char *groupname) {
    hid_t file_id = open_hdf5_file(filename);
    this->write_to_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
}
/** \brief <b> Reads `herm_matrix` from a given HDF5 group handle. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
*  \par Purpose
 * <!-- ========= -->
 *
 * > Reads the `herm_matrix`, including the matrix size, number of
 * > points on real and imaginary axis and fermionic/bosonic character,
 * > from a given HDF5 group handle.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param group_id
 * > [hid_t] The HDF5 group handle from which the `herm_matrix` is read.
 */
template <typename T>
void herm_matrix<T>::read_from_hdf5(hid_t group_id) {
    // -- Read dimensions
    int nt = read_primitive_type<int>(group_id, "nt");
    int ntau = read_primitive_type<int>(group_id, "ntau");
    int sig = read_primitive_type<int>(group_id, "sig");
    int size1 = read_primitive_type<int>(group_id, "size1");
    // RESIZE G
    this->resize(nt, ntau, size1);
    sig_ = sig;
    if (nt > -2) {
        hsize_t mat_size = (ntau + 1) * element_size_;
        read_primitive_type_array(group_id, "mat", mat_size, matptr(0));
    }
    if (nt > -1) {
        hsize_t ret_size = ((nt + 1) * (nt + 2)) / 2 * element_size_;
        hsize_t les_size = ((nt + 1) * (nt + 2)) / 2 * element_size_;
        hsize_t tv_size = ((nt + 1) * (ntau + 1)) * element_size_;
        read_primitive_type_array(group_id, "ret", ret_size, retptr(0, 0));
        read_primitive_type_array(group_id, "les", les_size, lesptr(0, 0));
        read_primitive_type_array(group_id, "tv", tv_size, tvptr(0, 0));
    }
}
/** \brief <b> Reads `herm_matrix` from a given HDF5 group handle and given group name. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Reads the `herm_matrix`, including the matrix size, number of
 * > points on real and imaginary axis and fermionic/bosonic character,
 * > from a given HDF5 group handle with given group name.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param group_id
 * > The HDF5 group handle from which the `herm_matrix` is read.
 * @param groupname
 * > The group name from which the `herm_matrix` is read,
 */
template <typename T>
void herm_matrix<T>::read_from_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(sub_group_id);
    close_group(sub_group_id);
}
/** \brief <b> Reads `herm_matrix` from a given HDF5 file and given group name. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Reads the `herm_matrix`, including the matrix size, number of
 * > points on real and imaginary axis and fermionic/bosonic character,
 * > from a given HDF5 file with given group name.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param filename
 * > The name of the HDF5 file.
 * @param groupname
 * > The group name from which the `herm_matrix` is read,
 */
template <typename T>
void herm_matrix<T>::read_from_hdf5(const char *filename,
                                    const char *groupname) {
    hid_t file_id = read_hdf5_file(filename);
    this->read_from_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
}
/** \brief <b> Reads the `herm_matrix` up to a given number of time steps from an HDF5 group. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Reads the `herm_matrix`, including the matrix size, number of
 * > points on real and imaginary axis and fermionic/bosonic character,
 * > from a given HDF5 group handle.
 * > The `herm_matrix` is read up to a given number of time steps `nt1`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param nt1
 * > Read up to this time step.
 * @param group_id
 * > The HDF5 group handle.
 */
template <typename T>
void herm_matrix<T>::read_from_hdf5(int nt1, hid_t group_id) {
    herm_matrix<T> gtmp;
    gtmp.read_from_hdf5(group_id);

    assert(nt1 >= -1 && nt1 <= gtmp.nt() && "nt1 >= -1 && nt1 <= gtmp.nt()");
    assert(nt1 >= -1 && nt1 <= nt_ && "nt1 >= -1 && nt1 <= nt_");
    assert(gtmp.size1() == size1_ && "gtmp.size1() == size1_");
    assert(gtmp.element_size() == element_size_ && "gtmp.element_size() == element_size_");
    assert(gtmp.ntau() == ntau_ && "gtmp.ntau() == ntau_");

    for (int n = -1; n <= nt1; n++)
        this->set_timestep(n, gtmp);
}
/** \brief <b> Reads the `herm_matrix` up to a given number of time steps from an HDF5 group. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Reads the `herm_matrix`, including the matrix size, number of
 * > points on real and imaginary axis and fermionic/bosonic character,
 * > from a given HDF5 group handle and group name.
 * > The `herm_matrix` is read up to a given number of time steps `nt1`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param nt1
 * > Read up to this time step.
 * @param group_id
 * > The HDF5 group handle.
 * @param groupname
 * > The HDF5 group name.
 */
template <typename T>
void herm_matrix<T>::read_from_hdf5(int nt1, hid_t group_id,
                                    const char *groupname) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(nt1, sub_group_id);
    close_group(sub_group_id);
}
/** \brief <b> Reads the `herm_matrix` up to a given number of time steps from an HDF5 file. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Reads the `herm_matrix`, including the matrix size, number of
 * > points on real and imaginary axis and fermionic/bosonic character,
 * > from a given HDF5 file and group name.
 * > The `herm_matrix` is read up to a given number of time steps `nt1`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param nt1
 * > Read up to this time step.
 * @param filename
 * > The name of the HDF5 file.
 * @param groupname
 *> The HDF5 group name.
 */
template <typename T>
void herm_matrix<T>::read_from_hdf5(int nt1, const char *filename,
                                    const char *groupname) {
    hid_t file_id = read_hdf5_file(filename);
    this->read_from_hdf5(nt1, file_id, groupname);
    close_hdf5_file(file_id);
}
/** \brief <b> Stores time slices of `herm_matrix` with a given interval to a HDF5 group handle. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Converts the `herm_matrix` every `dt` time steps to a time slice and stores them to
 * > the given HDF5 group handle.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param group_id
 * > The HDF5 group handle.
 * @param dt
 * > Store the slices every `dt` time steps.
 */
template <typename T>
void herm_matrix<T>::write_to_hdf5_slices(hid_t group_id, int dt) {
  assert(dt >= 1);

    char groupname[20];
    hid_t subgroup_id;
    for (int tstp = -1; tstp <= nt_; tstp++) {
        if (tstp == -1 || tstp % dt == 0) {
            herm_matrix_timestep_view<T> tmp(tstp, *this);
            std::sprintf(groupname, "t%d", tstp);
            subgroup_id = create_group(group_id, std::string(groupname));
            tmp.write_to_hdf5(subgroup_id);
        }
    }
}
 /** \brief <b> Stores time slices of `herm_matrix` with a given interval to a HDF5 group handle and group name. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Converts the `herm_matrix` every `dt` time steps to a time slice and stores them to
 * > the given HDF5 group handle with given group name.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param group_id
 * > The HDF5 group handle.
 * @param groupname
 * > The name of the HDF5 group.
 * @param dt
 * > Store the slices every `dt` time steps.
 */
template <typename T>
void herm_matrix<T>::write_to_hdf5_slices(hid_t group_id,
                                          const char *groupname, int dt) {
    hid_t sub_group_id = create_group(group_id, groupname);
    this->write_to_hdf5_slices(sub_group_id, dt);
    close_group(sub_group_id);
}
 /** \brief <b> Stores time slices of `herm_matrix` with a given interval to a HDF5 file under a given and group name. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Converts the `herm_matrix` every `dt` time steps to a time slice and stores them to
 * > the given HDF5 file under a given group name.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param filename
 * > The HDF5 output file.
 * @param groupname
 * > The name of the HDF5 group.
 * @param dt
 * > Store the slices every `dt` time steps.
 */
template <typename T>
void herm_matrix<T>::write_to_hdf5_slices(const char *filename,
                                          const char *groupname, int dt) {
    hid_t file_id = open_hdf5_file(filename);
    this->write_to_hdf5_slices(file_id, groupname, dt);
    close_hdf5_file(file_id);
}
/** \brief <b> Stores greater and lesser components of `herm_matrix` average-relative time representation to a given HDF5 group handle. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Stores greater and lesser components of `herm_matrix` \f$C(t,t^\prime)\f$ in Wigner time representation to a given HDF5 group handle. The
 * > representation is defined by \f$ \widetilde{C}^\gtrless}(T,t) = C^\gtrless(T+t/2,T-t/2)\f$. Every `dt` time steps are stored.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param group_id
 * > The HDF5 group handle.
 * @param dt
 * > Store the slices as function of the relative time every `dt` time steps.
 */
template <typename T>
void herm_matrix<T>::write_to_hdf5_tavtrel(hid_t group_id, int dt) {
  assert(dt >= 1);

    typedef std::complex<double> complex;
    char name[100];
    if (nt_ < 0)
        return;
    complex *ggtr = new complex[(nt_ + 1) * element_size_]; // temp storage
    complex *gles = new complex[(nt_ + 1) * element_size_]; // temp storage
    store_int_attribute_to_hid(group_id, std::string("nt"), nt_);
    store_int_attribute_to_hid(group_id, std::string("size1"), size1_);
    store_int_attribute_to_hid(group_id, std::string("size2"), size2_);
    store_int_attribute_to_hid(group_id, std::string("element_size"),
                               element_size_);
    hid_t group_id_les = create_group(group_id, std::string("les"));
    hid_t group_id_gtr = create_group(group_id, std::string("gtr"));
    hsize_t len_shape = 3, shape[3];
    shape[1] = size1_;
    shape[2] = size2_;
    for (int tav = 0; tav <= nt_; tav += dt) {
        int len = (tav <= nt_ - tav ? tav : nt_ - tav);
        std::sprintf(name, "%d", tav);
        shape[0] = len + 1;
        // read data: trel2 is trel/2
        for (int trel2 = 0; trel2 <= len; trel2++) {
            cdmatrix ggtrtt, glestt;
            this->get_les(tav + trel2, tav - trel2, glestt);
            this->get_gtr(tav + trel2, tav - trel2, ggtrtt);
            for (int i = 0; i < element_size_; i++) {
                gles[trel2 * element_size_ + i] =
                    glestt(i / size2_, i % size2_);
                ggtr[trel2 * element_size_ + i] =
                    ggtrtt(i / size2_, i % size2_);
            }
        }
        store_cplx_array_to_hid(group_id_les, std::string(name), gles, shape,
                                len_shape);
        store_cplx_array_to_hid(group_id_gtr, std::string(name), ggtr, shape,
                                len_shape);
    }
    delete[] ggtr;
    delete[] gles;
}
 /** \brief <b> Stores greater and lesser components of `herm_matrix` average-relative time representation to a given HDF5 group handle and group name. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Stores greater and lesser components of `herm_matrix` \f$C(t,t^\prime)\f$ in Wigner time representation to a given HDF5 group handle with
 * > a given group name. The representation is defined by
 * >  \f$ \widetilde{C}^\gtrless}(T,t) = C^\gtrless(T+t/2,T-t/2)\f$. Every `dt` time steps are stored.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param group_id
 * > [hid_t] The HDF5 group handle.
 * @param groupname
 * > [char*] The name of the HDF5 group.
 * @param dt
 * > [int] Store the slices as function of the relative time every `dt` time steps.
 */
template <typename T>
void herm_matrix<T>::write_to_hdf5_tavtrel(hid_t group_id,
                                           const char *groupname, int dt) {
    hid_t sub_group_id = create_group(group_id, groupname);
    this->write_to_hdf5_tavtrel(sub_group_id, dt);
    close_group(sub_group_id);
}
 /** \brief <b> Stores greater and lesser components of `herm_matrix` average-relative time representation to a given HDF5 file under given group name. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Stores greater and lesser components of `herm_matrix` \f$C(t,t^\prime)\f$ in Wigner time representation to a given HDF5 file with
 * > a given group name. The representation is defined by
 * >  \f$ \widetilde{C}^\gtrless}(T,t) = C^\gtrless(T+t/2,T-t/2)\f$. Every `dt` time steps are stored.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param filename
 * > The name of the output HDF5 file.
 * @param groupname
 * > The name of the HDF5 group.
 * @param dt
 * > Store the slices as function of the relative time every `dt` time steps.
 */
template <typename T>
void herm_matrix<T>::write_to_hdf5_tavtrel(const char *filename,
                                           const char *groupname, int dt) {
    hid_t file_id = open_hdf5_file(filename);
    this->write_to_hdf5_tavtrel(file_id, groupname, dt);
    close_hdf5_file(file_id);
}

#endif

/* #######################################################################################
#
#   SIMPLE OPERATIONS ON TIMESTEPS
#   NOTE: tstp IS A PHYSICAL TIME, tstp=-1 is the matsubara branch
#
########################################################################################*/
/** \brief <b> Sets all components at time step `tstp` to zero. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets all components of the `herm_matrix` at time step `tstp` to zero. If
 * > `tstp = -1`, only the Matsubara component will be set to zero.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > The time step at which the components are set to zero.
 *
 */
template <typename T>
void herm_matrix<T>::set_timestep_zero(int tstp) {
    assert(tstp >= -1 && tstp <= nt_ && "tstp >= -1 && tstp <= nt_");
    if (tstp == -1) {
        memset(matptr(0), 0, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memset(retptr(tstp, 0), 0, sizeof(cplx) * (tstp + 1) * element_size_);
        memset(tvptr(tstp, 0), 0, sizeof(cplx) * (ntau_ + 1) * element_size_);
        memset(lesptr(0, tstp), 0, sizeof(cplx) * (tstp + 1) * element_size_);
    }
}
/** \brief <b> Sets all components at time step `tstp` to the components of
 *  a given `herm_matrix`. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets all components of the `herm_matrix` at time step `tstp` to
 * > the components of given `herm_matrix`. If `tstp = -1`, only the
 * > Matsubara component will be copied.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > The time step at which the components are set.
 *
 * @param g1
 * > The `herm_matrix` from which the time step is copied.
 *
 */
template <typename T>
void herm_matrix<T>::set_timestep(int tstp, herm_matrix<T> &g1) {
    assert(tstp >= -1 && tstp <= nt_ && tstp <= g1.nt() && "tstp >= -1 && tstp <= nt_ && tstp <= g1.nt()");
    assert(g1.size1() == size1_ && "g1.size1() == size1_");
    assert(g1.ntau() == ntau_ && "g1.ntau() == ntau_");
    if (tstp == -1) {
        memcpy(mat_, g1.mat_, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memcpy(retptr(tstp, 0), g1.retptr(tstp, 0),
               sizeof(cplx) * (tstp + 1) * element_size_);
        memcpy(tvptr(tstp, 0), g1.tvptr(tstp, 0),
               sizeof(cplx) * (ntau_ + 1) * element_size_);
        memcpy(lesptr(0, tstp), g1.lesptr(0, tstp),
               sizeof(cplx) * (tstp + 1) * element_size_);
    }
}

/** \brief <b> Sets all components at time step `tstp` to the components of
 *  a given `herm_matrix_timestep`. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets all components of the `herm_matrix` at time step `tstp` to
 * > the components of given `herm_matrix_timestep`. If `tstp = -1`, only the
 * > Matsubara component will be copied.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > [int] The time step at which the components are set.
 *
 * @param timestep
 * > [herm_matrix_timestep] The `herm_matrix_timestep` from which the time step is copied.
 *
 */
template <typename T>
void herm_matrix<T>::set_timestep(int tstp,
                                  herm_matrix_timestep<T> &timestep) {
    cplx *x = timestep.data_;
    assert(tstp >= -1 && tstp <= nt_ && "(tstp >= -1 && tstp <= nt_");
    assert(timestep.tstp_ == tstp && timestep.ntau_ == ntau_ &&
           timestep.size1_ == size1_ &&
	   "timestep.tstp_ == tstp && timestep.ntau_ == ntau_ && timestep.size1_ == size1_");
    if (tstp == -1) {
        memcpy(mat_, x, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memcpy(retptr(tstp, 0), x, sizeof(cplx) * (tstp + 1) * element_size_);
        memcpy(tvptr(tstp, 0), x + (tstp + 1) * element_size_,
               sizeof(cplx) * (ntau_ + 1) * element_size_);
        memcpy(lesptr(0, tstp), x + (tstp + 1 + ntau_ + 1) * element_size_,
               sizeof(cplx) * (tstp + 1) * element_size_);
    }
}

template <typename T>
void herm_matrix<T>::get_timestep(int tstp,
                                  herm_matrix_timestep<T> &timestep) const {
    int len = (2 * (tstp + 1) + ntau_ + 1) * element_size_;
    cplx *x;
    assert(tstp >= -1 && tstp <= nt_ && "tstp >= -1 && tstp <= nt_");
    if (timestep.total_size_ < len)
        timestep.resize(tstp, ntau_, size1_);
    x = timestep.data_;
    timestep.tstp_ = tstp;
    timestep.ntau_ = ntau_;
    timestep.size1_ = size1_;
    if (tstp == -1) {
        memcpy(x, mat_, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memcpy(x, retptr(tstp, 0), sizeof(cplx) * (tstp + 1) * element_size_);
        memcpy(x + (tstp + 1) * element_size_, tvptr(tstp, 0),
               sizeof(cplx) * (ntau_ + 1) * element_size_);
        memcpy(x + (tstp + 1 + ntau_ + 1) * element_size_, lesptr(0, tstp),
               sizeof(cplx) * (tstp + 1) * element_size_);
    }
}

template <typename T>
void herm_matrix<T>::get_timestep(int tstp, herm_matrix<T> &g1) const {

    assert(tstp >= -1 && tstp <= nt_ && tstp <= g1.nt() && "tstp >= -1 && tstp <= nt_ && tstp <= g1.nt()");
    assert(g1.size1() == size1_ && "g1.size1() == size1_");
    assert(g1.ntau() == ntau_ && "g1.ntau() == ntau_");
    if (tstp == -1) {
        memcpy(g1.mat_, mat_, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memcpy(g1.retptr(tstp, 0), retptr(tstp, 0),
               sizeof(cplx) * (tstp + 1) * element_size_);
        memcpy(g1.tvptr(tstp, 0), tvptr(tstp, 0),
               sizeof(cplx) * (ntau_ + 1) * element_size_);
        memcpy(g1.lesptr(0, tstp), lesptr(0, tstp),
               sizeof(cplx) * (tstp + 1) * element_size_);
    }
}


/** \brief <b> Sets given matrix elements of all components at time step `tstp` to the components of
 *  a given `herm_matrix_timestep` of a scaler type. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets the (sub-) matrix elements `(i1,i2)` of all components at time step `tstp` to the components of
 * >  a given `herm_matrix_timestep` of scaler type.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > The time step at which the components are set.
 * @param i1
 * > First index for submatrix.
 * @param i2
 * > Second index for submatrix.
 *
 * @param g
 * > The `herm_matrix_timestep` from which the time step is copied.
 *
 */
template <typename T>
void herm_matrix<T>::set_matrixelement(int tstp, int i1, int i2,
                                       herm_matrix_timestep<T> &g) {
    int i, sij = i1 * size2_ + i2;
    assert(tstp == g.tstp_ && tstp <= nt_);
    assert(0 <= i1 && i1 < size1_ && 0 <= i2 && i2 < size2_ && g.size1() == 1);
    if (tstp == -1) {
        for (i = 0; i <= ntau_; i++)
            matptr(i)[sij] = *(g.matptr(i));
    } else {
        for (i = 0; i <= tstp; i++)
            retptr(tstp, i)[sij] = *(g.retptr(i));
        for (i = 0; i <= ntau_; i++)
            tvptr(tstp, i)[sij] = *(g.tvptr(i));
        for (i = 0; i <= tstp; i++)
            lesptr(i, tstp)[sij] = *(g.lesptr(i));
    }
}
/** \brief <b> Sets given matrix elements of all components at time step `tstp` to the components of
 *  a given `herm_matrix` of scaler type. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets the (sub-) matrix elements `(i1,i2)` of all components at time step `tstp` to the components at time step `tstp`  of
 * >  a given `herm_matrix`.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > The time step at which the components are set.
 * @param i1
 * > First index for submatrix.
 * @param i2
 * > Second index for submatrix.
 *
 * @param g
 * > The `herm_matrix` from which the time step is copied.
 *
 */
template <typename T>
void herm_matrix<T>::set_matrixelement(int tstp, int i1, int i2,
                                       herm_matrix<T> &g) {
    int i, sij = i1 * size2_ + i2;
    assert(tstp <= g.nt() && tstp <= nt_ && "tstp <= g.nt() && tstp <= nt_");
    assert(0 <= i1 && i1 < size1_ && 0 <= i2 && i2 < size2_ && g.size1() == 1
           && "0 <= i1 && i1 < size1_ && 0 <= i2 && i2 < size2_ && g.size1() == 1");
    if (tstp == -1) {
        for (i = 0; i <= ntau_; i++)
            matptr(i)[sij] = *(g.matptr(i));
    } else {
        for (i = 0; i <= tstp; i++)
            retptr(tstp, i)[sij] = *(g.retptr(tstp, i));
        for (i = 0; i <= ntau_; i++)
            tvptr(tstp, i)[sij] = *(g.tvptr(tstp, i));
        for (i = 0; i <= tstp; i++)
            lesptr(i, tstp)[sij] = *(g.lesptr(i, tstp));
    }
}

/** \brief <b> Sets given matrix elements of all components at time step `tstp` to
* given matrix elements of the corresponding components of a given `herm_matrix_timestep` of matrix type. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Sets the (sub-) matrix elements `(i1,i2)` of all components at time step `tstp` to the matrix elements `(j1,j2)`
* > of the corresponding componentof a given `herm_matrix_timestep`.
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > The time step at which the components are set.
* @param i1
* > First index for submatrix.
* @param i2
* > Second index for submatrix.
*
* @param g
* > The `herm_matrix_timestep` of a matrix type from which the time step is copied.
* @param j1
* > First index for submatrix of `g`.
* @param j2
* > Second index for submatrix of `g`.
*
*/

template <typename T>
void herm_matrix<T>::set_matrixelement(int tstp, int i1, int i2,
                                       herm_matrix_timestep<T> &g, int j1,
                                       int j2) {
    int i, sij = i1 * size2_ + i2, tij = j1 * g.size2() + j2;
    assert(tstp == g.tstp_ && tstp <= nt_ && "tstp == g.tstp_ && tstp <= nt_");
    assert(0 <= i1 && i1 < size1_ && 0 <= i2 && i2 < size2_ && "0 <= i1 && i1 < size1_ && 0 <= i2 && i2 < size2_");
    assert(0 <= j1 && j1 < g.size1() && 0 <= j2 && j2 < g.size1()
	   && "0 <= j1 && j1 < g.size1() && 0 <= j2 && j2 < g.size1()");
    if (tstp == -1) {
        for (i = 0; i <= ntau_; i++)
            matptr(i)[sij] = g.matptr(i)[tij];
    } else {
        for (i = 0; i <= tstp; i++)
            retptr(tstp, i)[sij] = g.retptr(i)[tij];
        for (i = 0; i <= ntau_; i++)
            tvptr(tstp, i)[sij] = g.tvptr(i)[tij];
        for (i = 0; i <= tstp; i++)
            lesptr(i, tstp)[sij] = g.lesptr(i)[tij];
    }
}

/** \brief <b> Sets given matrix elements of all components at time step `tstp` to
* given matrix elements the components of a given `herm_matrix` of matrix type. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Sets the (sub-) matrix elements `(i1,i2)` of all components at time step `tstp`
* > to matrix elements `(j1,j2)` of the corresponding componentat time step `tstp` of
* >  a given `herm_matrix` of matrix type.
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > The time step at which the components are set.
* @param i1
* > First index for submatrix.
* @param i2
* > Second index for submatrix.
*
* @param g
* > The `herm_matrix` from which the time step is copied.
* @param j1
* > First index for submatrix of `g`.
* @param j2
* > Second index for submatrix of `g`.
*
*/

template <typename T>
void herm_matrix<T>::set_matrixelement(int tstp, int i1, int i2,
                                       herm_matrix<T> &g, int j1, int j2) {
    int i, sij = i1 * size2_ + i2, tij = j1 * g.size2() + j2;
    assert(tstp <= g.nt_ && tstp <= nt_ && "tstp <= g.nt_ && tstp <= nt_");
    assert(0 <= i1 && i1 < size1_ && 0 <= i2 && i2 < size2_
	   && "0 <= i1 && i1 < size1_ && 0 <= i2 && i2 < size2_");
    assert(0 <= j1 && j1 < g.size1() && 0 <= j2 && j2 < g.size2()
	   && "0 <= j1 && j1 < g.size1() && 0 <= j2 && j2 < g.size2()");
    if (tstp == -1) {
        for (i = 0; i <= ntau_; i++)
            matptr(i)[sij] = g.matptr(i)[tij];
    } else {
        for (i = 0; i <= tstp; i++)
            retptr(tstp, i)[sij] = g.retptr(tstp, i)[tij];
        for (i = 0; i <= ntau_; i++)
            tvptr(tstp, i)[sij] = g.tvptr(tstp, i)[tij];
        for (i = 0; i <= tstp; i++)
            lesptr(i, tstp)[sij] = g.lesptr(i, tstp)[tij];
    }
}


/** \brief <b> Sets given matrix elements of all components to the for all the time components of
*  a given `herm_matrix` of scaler type for all the time. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Sets the (sub-) matrix elements `(i1,i2)` of all components  to the components of
* >  a given `herm_matrix` of scaler type for all the time.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i1
* > First index for submatrix.
* @param i2
* > Second index for submatrix.
*
* @param g
* > The `herm_matrix` from which the time steps are copied.
*
*/

template<typename T> void herm_matrix<T>::set_matrixelement(int i1,int i2,herm_matrix<T> &g){
	assert(nt_ == g.nt_);
	for(int tstp=-1;tstp<=nt_;tstp++){
		this->set_matrixelement(tstp,i1,i2,g);
	}
}

/** \brief <b> Sets given matrix elements of all components to
* given matrix elements of the correponding components of a given `herm_matrix` of matrix type for all the time. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Sets the (sub-) matrix elements `(i1,i2)` of all components
* > to matrix elements `(j1,j2)` of the corresponding componentof
* >  a given `herm_matrix` of matrix type for all the time.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i1
* > First index for submatrix.
* @param i2
* > Second index for submatrix.
*
* @param g
* > The `herm_matrix` from which the time steps are copied.
* @param j1
* > First index for submatrix of `g`.
* @param j2
* > Second index for submatrix of `g`.
*
*/

template<typename T> void herm_matrix<T>::set_matrixelement(int i1,int i2,herm_matrix<T> &g,int j1,int j2){
	assert(nt_ == g.nt_);
	for(int tstp=-1;tstp<=nt_;tstp++){
		this->set_matrixelement(tstp,i1,i2,g,j1,j2);
	}
}

// Set the submatrix of the Green'd function
// G_{j_1(k) j_2(k)}=G_{i_1(k) i_2(k) }, where k is an iterator over subspace


/** \brief <b> Sets a (sub-) matrix of this contour object at all time steps
* to a (sub-) matrix of a given `herm_matrix`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Sets the (sub-) matrix of a two-time contour object \f$C(t,t')\f$ \f$C(t,t')\f$ the (sub-)
* of a given `herm_matrix` \f$g\f$ according to
* \f$ C_{i_1(k) i_2(k)}(t,t') = g_{j_1(k) j_2(k)}(t,t') \f$ with \f$k\f$ denoting
* an iterator of a subspace for all times \f$t,t'\f$.
* <!-- ARGUMENTS
*      ========= -->
*
* @param i1
* > Vector of row indices of `C`.
* @param i2
* > Vector of column indices of `C`.
* @param g
* > The `herm_matrix` from which the (sub-) matrix is copied.
* @param j1
* > Vector of row indices of `g`.
* @param j2
* > Vector of row indices of `g`.
*/

template<typename T> void herm_matrix<T>::set_submatrix(std::vector<int> &i1,
  std::vector<int> &i2,herm_matrix<T> &g,std::vector<int> &j1,std::vector<int> &j2){
	assert(nt_ == g.nt_);
	assert(i1.size()==i2.size() && i1.size()==j1.size() && j1.size()==j2.size());
	assert(size1_*size2_==i1.size());
	for(int k1=0;k1<i1.size();k1++){
		this->set_matrixelement(i1[k1],i2[k1],g,j1[k1],j2[k1]);
	}
}

/** \brief <b> Sets a (sub-) matrix of this contour object at a given time step
* to a (sub-) matrix of a given `herm_matrix`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Sets the (sub-) matrix of a two-time contour object \f$C(t,t')\f$ \f$C(t,t')\f$ the (sub-)
* of a given `herm_matrix` \f$g\f$ according to
* \f$ C_{i_1(k) i_2(k)}(t,t') = g_{j_1(k) j_2(k)}(t,t') \f$ with \f$k\f$ denoting
* an iterator of a subspace (at a given time step).
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Time step `tstp`.
* @param i1
* > Vector of row indices of `C`.
* @param i2
* > Vector of column indices of `C`.
* @param g
* > The `herm_matrix` from which the (sub-) matrix is copied.
* @param j1
* > Vector of row indices of `g`.
* @param j2
* > Vector of row indices of `g`.
*/

template<typename T> void herm_matrix<T>::set_submatrix(int tstp, std::vector<int> &i1,
  std::vector<int> &i2,herm_matrix<T> &g,std::vector<int> &j1,std::vector<int> &j2){
    assert(tstp <= nt_);
    assert(nt_ == g.nt_);
    assert(i1.size()==i2.size() && i1.size()==j1.size() && j1.size()==j2.size());
    assert(size1_*size2_==i1.size());
    for(int k1=0;k1<i1.size();k1++){
        this->set_matrixelement(tstp,i1[k1],i2[k1],g,j1[k1],j2[k1]);
    }
}


#define HERM_MATRIX_INCR_TSTP                                                \
    if (alpha == cplx(1.0, 0.0)) {                                           \
        for (i = 0; i < len; i++)                                            \
            x0[i] += x[i];                                                   \
    } else {                                                                 \
        for (i = 0; i < len; i++)                                            \
            x0[i] += alpha * x[i];                                           \
    }



/** \brief <b> Adds a `herm_matrix_timestep` with given weight to the `herm_matrix`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Performs the operation \f$C \rightarrow C + \alpha A\f$, where \f$C\f$ is the
* > `herm_matrix`, \f$A\f$ is a time slice described by a `herm_matrix_timestep` and \f$\alpha\f$
* > is a complex weight. The operation is performed at given time step `tstp`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > [int] The time step where the `herm_matrix` is increased.
* @param timestep
* > [herm_matrix_timestep] The `herm_matrix_timestep` which is added to the `herm_matrix`.
* @param alpha
* > [complex<T>] The weight in front of `timestep`.
*/
template <typename T>
void herm_matrix<T>::incr_timestep(int tstp,
                                   herm_matrix_timestep<T> &timestep,
                                   std::complex<T> alpha) {
    int i, len;
    cplx *x, *x0;
    assert(tstp >= -1 && tstp <= nt_ && "tstp >= -1 && tstp <= nt_");
    assert(timestep.tstp_ == tstp && timestep.ntau_ == ntau_ && timestep.size1_ == size1_
      && "timestep.tstp_ == tstp && timestep.ntau_ == ntau_ && timestep.size1_ == size1_");
    if (tstp == -1) {
        len = (ntau_ + 1) * element_size_;
        x0 = matptr(0);
        x = timestep.data_;
        HERM_MATRIX_INCR_TSTP
    } else {
        len = (tstp + 1) * element_size_;
        x0 = retptr(tstp, 0);
        x = timestep.data_;
        HERM_MATRIX_INCR_TSTP
        len = (ntau_ + 1) * element_size_;
        x0 = tvptr(tstp, 0);
        x = timestep.data_ + (tstp + 1) * element_size_;
        HERM_MATRIX_INCR_TSTP
        len = (tstp + 1) * element_size_;
        x0 = lesptr(0, tstp);
        x = timestep.data_ + (tstp + 1 + ntau_ + 1) * element_size_;
        HERM_MATRIX_INCR_TSTP
    }
}
#undef HERM_MATRIX_INCR_TSTP

/** \brief <b> Adds a `herm_matrix` with given weight to the `herm_matrix` at given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Performs the operation \f$C \rightarrow C + \alpha g\f$, where \f$C\f$ is the
* `herm_matrix`, \f$g\f$ is a time slice taken from a `herm_matrix` and \f$\alpha\f$
* is a complex weight. The operation is performed at given time step `tstp`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > [int] The time step where the `herm_matrix` is increased.
* @param g
* > [herm_matrix] The `herm_matrix` \f$g\f$ which is added to the `herm_matrix` at the time step.
* @param alpha
* > [complex<T>] The weight in front of \f$g\f$.
*/
template <typename T>
void herm_matrix<T>::incr_timestep(int tstp, herm_matrix<T> &g,
                                   std::complex<T> alpha) {
    int m;
    cplx *x, *x0;
    assert(tstp >= -1 && tstp <= nt_ && "tstp >= -1 && tstp <= nt_");
    assert(g.nt_ >= tstp && g.ntau_ == ntau_ && g.size1_ == size1_
	   && "g.nt_ >= tstp && g.ntau_ == ntau_ && g.size1_ == size1_");
    if (tstp == -1) {
        x0 = matptr(0);
        x = g.matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_incr<T, LARGESIZE>(size1_, x0 + m * element_size_, alpha,
                                       x + m * element_size_);
        }
    } else {
        x0 = retptr(tstp, 0);
        x = g.retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_incr<T, LARGESIZE>(size1_, x0 + m * element_size_, alpha,
                                       x + m * element_size_);
        }
        x0 = tvptr(tstp, 0);
        x = g.tvptr(tstp, 0);
        for (m = 0; m <= ntau_; m++) {
            element_incr<T, LARGESIZE>(size1_, x0 + m * element_size_, alpha,
                                       x + m * element_size_);
        }
        x0 = lesptr(0, tstp);
        x = g.lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_incr<T, LARGESIZE>(size1_, x0 + m * element_size_, alpha,
                                       x + m * element_size_);
        }
    }
}
/** \brief <b> Adds a `herm_matrix` with given weight to the `herm_matrix` at all time steps. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Performs the operation \f$C \rightarrow C + \alpha g\f$, where \f$C\f$ is the
* `herm_matrix`, \f$g\f$ is a `herm_matrix` and \f$\alpha\f$
* is a complex weight.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > [herm_matrix] The `herm_matrix` \f$g\f$ which is added to the `herm_matrix`.
* @param alpha
* > [complex<T>] The weight in front of \f$g\f$`.
*/
template <typename T>
void herm_matrix<T>::incr_timestep(herm_matrix<T> &g, std::complex<T> alpha) {
    assert(g.nt_ >= nt_ && g.ntau_ == ntau_ && g.size1_ == size1_
	   && "g.nt_ >= nt_ && g.ntau_ == ntau_ && g.size1_ == size1_");
    for (int m = -1; m <= nt_; m++)
        this->incr_timestep(m, g, alpha);
}
/// @private
// G(t,t') ==> F(t)G(t,t')   ... ft+t*element_size_ points to F(t)

/** \brief <b> Left-multiplies the `herm_matrix` with contour function. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Performs the operation \f$C(t,t^\prime) \rightarrow w F(t)C(t,t^\prime)\f$, where \f$C(t,t^\prime)\f$ is the
* `herm_matrix`, \f$F(t)\f$ is a `function` given in pointer format and \f$w\f$ is a real weight. The operation is performed
* at given time step `tstp`. This is a lower-level routine. The interface where \f$F(t)\f$ is supplied as
* contour function is preferred.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > [int] The time step at which \f$F(t)\f$ and \f$C(t,t^\prime)\f$ are multiplied.
* @param *f0
* > [complex<T>] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix).
* @param *ft
* > [complex<T>] Pointer to \f$F(t)\f$ on the real axis. It is assumed that ft+t*element_size_ points to \f$ F(t) \f$.
* @param weight
* > [T] The weight as above.
*/
template <typename T>
void herm_matrix<T>::left_multiply(int tstp, std::complex<T> *f0,
                                   std::complex<T> *ft, T weight) {
    int m;
    cplx *xtemp, *ftemp, *x0;
    xtemp = new cplx[element_size_];
    assert(tstp >= -1 && tstp <= nt_ && "tstp >= -1 && tstp <= nt_");
    if (tstp == -1) {
        x0 = matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, f0,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    } else {
        ftemp = ft + tstp * element_size_;
        x0 = retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, ftemp,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = tvptr(tstp, 0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, ftemp,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, ft + m * element_size_,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    }
    delete[] xtemp;
}

/** \brief <b> Left-multiplies the `herm_matrix` with contour function. Preferred interface. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Performs the operation \f$C(t,t^\prime) \rightarrow w F(t)C(t,t^\prime)\f$, where \f$C(t,t^\prime)\f$ is the
* `herm_matrix`, \f$F(t)\f$ is a contour `function` and \f$w\f$ is a real weight. The operation is performed
* at given time step `tstp`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > [int] The time step at which \f$F(t)\f$ and \f$C(t,t^\prime)\f$ are multiplied.
* @param ft
* > [function] The contour function.
* @param weight
* > [T] The weight as above.
*/
template <typename T>
void herm_matrix<T>::left_multiply(int tstp, function<T> &ft, T weight) {
    switch (size1_) {
    case 1:
        left_multiply_(tstp, ft, weight, std::integral_constant<int, 1>{});
        break;
    case 2:
        left_multiply_(tstp, ft, weight, std::integral_constant<int, 2>{});
        break;
    case 3:
        left_multiply_(tstp, ft, weight, std::integral_constant<int, 3>{});
        break;
    case 4:
        left_multiply_(tstp, ft, weight, std::integral_constant<int, 4>{});
        break;
    case 5:
        left_multiply_(tstp, ft, weight, std::integral_constant<int, 5>{});
        break;
    case 8:
        left_multiply_(tstp, ft, weight, std::integral_constant<int, 8>{});
        break;
    default:
        left_multiply_(tstp, ft, weight, std::integral_constant<int, LARGESIZE>{});
        break;
    }
}
/// @private
/** \brief <b> Left-multiplies the `herm_matrix` with contour function. Internal version, not intended
*  to be called directly. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Performs the operation \f$C(t,t^\prime) \rightarrow w F(t)C(t,t^\prime)\f$, where \f$C(t,t^\prime)\f$ is the
* `herm_matrix`, \f$F(t)\f$ is a contour `function` and \f$w\f$ is a real weight. The operation is performed
* at given time step `tstp`. This is an internal version with specializations with respect to the matrix size.
*
* \note parameter SIZE determines the matrix size specified by the top-level interface.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > [int] The time step at which \f$F(t)\f$ and \f$C(t,t^\prime)\f$ are multiplied.
* @param ft
* > [function] The contour function.
* @param weight
* > [T] The weight as above.
*/
template <typename T>
template <int SIZE>
void herm_matrix<T>::left_multiply_(int tstp, function<T> &ft, T weight,
                                    std::integral_constant<int, SIZE>) {
    int m;
    cplx *xtemp, *ftemp, *x0, *f0;
    xtemp = new cplx[element_size_];
    assert(tstp >= -1 && tstp <= nt_ && ft.nt() >= tstp &&
           ft.size1() == size1_ && ft.size2() == size2_
	   && "tstp >= -1 && tstp <= nt_ && ft.nt() >= tstp && ft.size1() == size1_ && ft.size2() == size2_");
    f0 = ft.ptr(-1);
    if (tstp == -1) {
        x0 = matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, SIZE>(size1_, xtemp, f0, x0 + m * element_size_);
            element_smul<T, SIZE>(size1_, xtemp, weight);
            element_set<T, SIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    } else {
        ftemp = ft.ptr(tstp);
        x0 = retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, SIZE>(size1_, xtemp, ftemp, x0 + m * element_size_);
            element_smul<T, SIZE>(size1_, xtemp, weight);
            element_set<T, SIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = tvptr(tstp, 0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, SIZE>(size1_, xtemp, ftemp, x0 + m * element_size_);
            element_smul<T, SIZE>(size1_, xtemp, weight);
            element_set<T, SIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, SIZE>(size1_, xtemp, ft.ptr(m), x0 + m * element_size_);
            element_smul<T, SIZE>(size1_, xtemp, weight);
            element_set<T, SIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    }
    delete[] xtemp;
}

/** \brief <b> Left-multiplies the `herm_matrix` with the hermitian conjugate of a contour function. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Performs the operation \f$C(t,t^\prime) \rightarrow w F^\ddagger(t)C(t,t^\prime)\f$, where \f$C(t,t^\prime)\f$ is the
* `herm_matrix`, \f$F(t)\f$ is a contour `function` and \f$w\f$ is a real weight. The operation is performed
* at given time step `tstp`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > [int] The time step at which \f$F^\ddagger(t)\f$ and \f$C(t,t^\prime)\f$ are multiplied.
* @param ft
* > [function] The contour function.
* @param weight
* > [T] The weight as above.
*/
template <typename T>
void herm_matrix<T>::left_multiply_hermconj(int tstp, function<T> &ft,
                                            T weight) {
    int m;
    cplx *xtemp, *ftemp, *x0, *f0, *fcc;
    xtemp = new cplx[element_size_];
    fcc = new cplx[element_size_];
    assert(tstp >= -1 && tstp <= nt_ && ft.nt() >= tstp &&
           ft.size1() == size1_ && ft.size2() == size2_ &&
	   "tstp >= -1 && tstp <= nt_ && ft.nt() >= tstp && ft.size1() == size1_ && ft.size2() == size2_");

    f0 = ft.ptr(-1);
    element_conj<T, LARGESIZE>(size1_, fcc, f0);
    if (tstp == -1) {
        x0 = matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, fcc,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    } else {
        ftemp = ft.ptr(tstp);
        element_conj<T, LARGESIZE>(size1_, fcc, ftemp);
        x0 = retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, fcc,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = tvptr(tstp, 0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, fcc,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_conj<T, LARGESIZE>(size1_, fcc, ft.ptr(m));
            element_mult<T, LARGESIZE>(size1_, xtemp, fcc,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    }
    delete[] xtemp;
    delete[] fcc;
}
/// @private
/** \brief <b> Right-multiplies the `herm_matrix` with contour function. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Performs the operation \f$C(t,t^\prime) \rightarrow w C(t,t^\prime) F(t^\prime)\f$, where \f$C(t,t^\prime)\f$ is the
* `herm_matrix`, \f$F(t)\f$ is a `function` given in pointer format and \f$w\f$ is a real weight. The operation is performed
* at given time step `tstp`. This is a lower-level routine. The interface where \f$F(t)\f$ is supplied as
* contour function is preferred.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > [int] The time step at which \f$F(t)\f$ and \f$C(t,t^\prime)\f$ are multiplied.
* @param *f0
* > [complex<T>] Pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix).
* @param *ft
* > [complex<T>] Pointer to \f$F(t)\f$ on the real axis. It is assumed that ft+t*element_size_ points to \f$ F(t) \f$.
* @param weight
* > [T] The weight as above.
*/
template <typename T>
void herm_matrix<T>::right_multiply(int tstp, std::complex<T> *f0,
                                    std::complex<T> *ft, T weight) {
    int m;
    cplx *xtemp, *ftemp, *x0;
    xtemp = new cplx[element_size_];
    assert(tstp >= -1 && tstp <= nt_ && "tstp >= -1 && tstp <= nt_");
    if (tstp == -1) {
        x0 = matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       f0);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    } else {
        x0 = retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       ft + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = tvptr(tstp, 0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       f0);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        ftemp = ft + tstp * element_size_;
        x0 = lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       ftemp);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    }
    delete[] xtemp;
}
/** \brief <b> Right-multiplies the `herm_matrix` with contour function. Preferred interface. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Performs the operation \f$C(t,t^\prime) \rightarrow w C(t,t^\prime)F(t^\prime)\f$, where \f$C(t,t^\prime)\f$ is the
* `herm_matrix`, \f$F(t)\f$ is a contour `function` and \f$w\f$ is a real weight. The operation is performed
* at given time step `tstp`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > [int] The time step at which \f$F(t)\f$ and \f$C(t,t^\prime)\f$ are multiplied.
* @param ft
* > [function] The contour function.
* @param weight
* > [T] The weight as above.
*/
template <typename T>
void herm_matrix<T>::right_multiply(int tstp, function<T> &ft, T weight) {
    switch (size1_) {
    case 1:
        right_multiply_(tstp, ft, weight, std::integral_constant<int, 1>{});
        break;
    case 2:
        right_multiply_(tstp, ft, weight, std::integral_constant<int, 2>{});
        break;
    case 3:
        right_multiply_(tstp, ft, weight, std::integral_constant<int, 3>{});
        break;
    case 4:
        right_multiply_(tstp, ft, weight, std::integral_constant<int, 4>{});
        break;
    case 5:
        right_multiply_(tstp, ft, weight, std::integral_constant<int, 5>{});
        break;
    case 8:
        right_multiply_(tstp, ft, weight, std::integral_constant<int, 8>{});
        break;
    default:
        right_multiply_(tstp, ft, weight, std::integral_constant<int, LARGESIZE>{});
        break;
    }
}
/// @private
/** \brief <b> Right-multiplies the `herm_matrix` with contour function.
* Internal version, not intended to be called directly. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Performs the operation \f$C(t,t^\prime) \rightarrow w C(t,t^\prime)F(t^\prime)\f$, where \f$C(t,t^\prime)\f$ is the
* `herm_matrix`, \f$F(t)\f$ is a contour `function` and \f$w\f$ is a real weight. The operation is performed
* at given time step `tstp`. This is an internal version with specializations with respect to the matrix size.
*
* \note Parameter SIZE determines the matrix size specified by the top-level interface.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > The time step at which \f$F(t)\f$ and \f$C(t,t^\prime)\f$ are multiplied.
* @param ft
* > The contour function.
* @param weight
* > The weight as above.
*/
template <typename T>
template <int SIZE>
void herm_matrix<T>::right_multiply_(int tstp, function<T> &ft, T weight,
                                     std::integral_constant<int, SIZE>) {
    int m;
    cplx *xtemp, *ftemp, *x0, *f0;
    xtemp = new cplx[element_size_];
    assert(ft.size1() == size1_ && "ft.size1() == size1_");
    assert(ft.nt() >= tstp && "ft.nt() >= tstp");
    assert(tstp <= nt_ && "tstp <= nt_");
    assert(tstp >= -1 && "tstp >= -1");
    f0 = ft.ptr(-1);
    if (tstp == -1) {
        x0 = matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, SIZE>(size1_, xtemp, x0 + m * element_size_, f0);
            element_smul<T, SIZE>(size1_, xtemp, weight);
            element_set<T, SIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    } else {
        x0 = retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, SIZE>(size1_, xtemp, x0 + m * element_size_, ft.ptr(m));
            element_smul<T, SIZE>(size1_, xtemp, weight);
            element_set<T, SIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = tvptr(tstp, 0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, SIZE>(size1_, xtemp, x0 + m * element_size_, f0);
            element_smul<T, SIZE>(size1_, xtemp, weight);
            element_set<T, SIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        ftemp = ft.ptr(tstp);
        x0 = lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, SIZE>(size1_, xtemp, x0 + m * element_size_, ftemp);
            element_smul<T, SIZE>(size1_, xtemp, weight);
            element_set<T, SIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    }
    delete[] xtemp;
}

/** \brief <b> Right-multiplies the `herm_matrix` with the
*   hermitian conjugate of a contour function. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Performs the operation \f$C(t,t^\prime) \rightarrow w C(t,t^\prime)F^\ddagger(t^\prime)\f$, where \f$C(t,t^\prime)\f$ is the
* `herm_matrix`, \f$F(t)\f$ is a contour `function` and \f$w\f$ is a real weight. The operation is performed
* at given time step `tstp`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > The time step at which \f$F^\ddagger(t)\f$ and \f$C(t,t^\prime)\f$ are multiplied.
* @param ft
* > The contour function.
* @param weight
* > The weight as above.
*/
template <typename T>
void herm_matrix<T>::right_multiply_hermconj(int tstp, function<T> &ft,
                                             T weight) {
    int m;
    cplx *xtemp, *ftemp, *x0, *f0, *fcc;
    xtemp = new cplx[element_size_];
    fcc = new cplx[element_size_];
    assert(ft.size1() == size1_ && "ft.size1() == size1_");
    assert(ft.nt() >= tstp && "ft.nt() >= tstp");
    assert(tstp <= nt_&& "ft.nt() >= tstp");
    assert(tstp >= -1&& "ft.nt() >= tstp");
    f0 = ft.ptr(-1);
    element_conj<T, LARGESIZE>(size1_, fcc, f0);
    if (tstp == -1) {
        x0 = matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       fcc);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    } else {
        x0 = retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_conj<T, LARGESIZE>(size1_, fcc, ft.ptr(m));
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       fcc);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = tvptr(tstp, 0);
        element_conj<T, LARGESIZE>(size1_, fcc, f0);
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       fcc);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        ftemp = ft.ptr(tstp);
        element_conj<T, LARGESIZE>(size1_, fcc, ftemp);
        x0 = lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       fcc);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    }
    delete[] fcc;
    delete[] xtemp;
}

/** \brief <b> Multiplies all component of the `herm_matrix` with a real scalar. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * Multiplies all component of the `herm_matrix` at a given time step
 * `tstp` with a real scalar `weight`.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > The time step at which the contour object is weighted with a factor.
 * @param weight
 * > The weight.
 */
template <typename T>
void herm_matrix<T>::smul(int tstp, T weight) {
    int m;
    cplx *x0;
    assert(tstp >= -1 && tstp <= nt_);
    if (tstp == -1) {
        x0 = matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
                                       weight);
        }
    } else {
        x0 = retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
                                       weight);
        }
        x0 = tvptr(tstp, 0);
        for (m = 0; m <= ntau_; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
                                       weight);
        }
        x0 = lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
                                       weight);
        }
    }
}

/** \brief <b> Multiplies all component of the `herm_matrix` with a complex scalar. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * Multiplies all component of the `herm_matrix` at a given time step
 * `tstp` with a complex scalar `weight`.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > [int] The time step at which the contour object is weighted with a factor.
 * @param weight
 * > [T] The weight.
 */

template <typename T>
void herm_matrix<T>::smul(int tstp, std::complex<T> weight) {
    int m;
    cplx *x0;
    assert(tstp >= -1 && tstp <= nt_ && "tstp >= -1 && tstp <= nt_");
    if (tstp == -1) {
        x0 = matptr(0);
        for (m = 0; m <= ntau_; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
                                       weight);
        }
    } else {
        x0 = retptr(tstp, 0);
        for (m = 0; m <= tstp; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
                                       weight);
        }
        x0 = tvptr(tstp, 0);
        for (m = 0; m <= ntau_; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
                                       weight);
        }
        x0 = lesptr(0, tstp);
        for (m = 0; m <= tstp; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
                                       weight);
        }
    }
}

/* #######################################################################################
#
#   MPI UTILS
#
########################################################################################*/
#if CNTR_USE_MPI == 1

/** \brief <b> MPI reduce for the `herm_matrix` at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > MPI reduce for the `herm_matrix` to the `root` at a given time step.
* > Works for scalar or square-matrix contour objects.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > time step
* @param root
* > Index of root
*/
template <typename T>
void herm_matrix<T>::Reduce_timestep(int tstp, int root) {
    assert(tstp <= nt_);

    herm_matrix_timestep<T> Gtemp;
    Gtemp.resize(tstp, ntau_, size1_);
    this->get_timestep(tstp, Gtemp);
    
    Gtemp.Reduce_timestep(tstp, root);

    this->set_timestep(tstp, Gtemp);
}


/** \brief <b> Broadcasts the `herm_matrix` at a given time step to all tasks. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Broadcasts the `herm_matrix` at a given time step `tstp` to all tasks.
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Time step which should be broadcasted.
* @param root
* > The task rank from which the `herm_matrix` should be broadcasted.
*/

template <typename T>
void herm_matrix<T>::Bcast_timestep(int tstp, int root) {
    int numtasks, taskid;
    herm_matrix_timestep<T> Gtemp;
    numtasks = MPI::COMM_WORLD.Get_size();
    taskid = MPI::COMM_WORLD.Get_rank();
    Gtemp.resize(tstp, ntau_, size1_);
    if (taskid == root)
        this->get_timestep(tstp, Gtemp);
    if (sizeof(T) == sizeof(double))
        MPI::COMM_WORLD.Bcast(Gtemp.data_, Gtemp.total_size_,
                              MPI::DOUBLE_COMPLEX, root);
    else
        MPI::COMM_WORLD.Bcast(Gtemp.data_, Gtemp.total_size_, MPI::COMPLEX,
                              root);
    if (taskid != root)
        this->set_timestep(tstp, Gtemp);
}

/** \brief <b> Sends the `herm_matrix` at a given time step to a specific task. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Sends the `herm_matrix` at a given time step `tstp` to a specific task with
* rank `dest`.
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Time step which should be send.
* @param dest
* > The task rank to which the `herm_matrix` should be send.
* @param tag
* > The MPI error flag.
*/
template <typename T>
void herm_matrix<T>::Send_timestep(int tstp, int dest, int tag) {
    int taskid = MPI::COMM_WORLD.Get_rank();
    if (!(taskid == dest)) {
        herm_matrix_timestep<T> Gtemp;
        Gtemp.resize(tstp, ntau_, size1_);
        this->get_timestep(tstp, Gtemp);
        if (sizeof(T) == sizeof(double))
            MPI::COMM_WORLD.Send(Gtemp.data_, Gtemp.total_size_,
                                 MPI::DOUBLE_COMPLEX, dest, tag);
        else
            MPI::COMM_WORLD.Send(Gtemp.data_, Gtemp.total_size_, MPI::COMPLEX,
                                 dest, tag);
    }
}

/** \brief <b> Recevies the `herm_matrix` at a given time step from a specific task. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * Receives the `herm_matrix` at a given time step `tstp` from a specific task with
 * rank `root`.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > Time step which should be received.
 * @param root
 * > The task rank from which the `herm_matrix` should be received.
 * @param tag
 * > The MPI error flag.
 */
template <typename T>
void herm_matrix<T>::Recv_timestep(int tstp, int root, int tag) {
    int taskid = MPI::COMM_WORLD.Get_rank();
    if (!(taskid == root)) {
        herm_matrix_timestep<T> Gtemp;
        Gtemp.resize(tstp, ntau_, size1_);
        if (sizeof(T) == sizeof(double))
            MPI::COMM_WORLD.Recv(Gtemp.data_, Gtemp.total_size_,
                                 MPI::DOUBLE_COMPLEX, root, tag);
        else
            MPI::COMM_WORLD.Recv(Gtemp.data_, Gtemp.total_size_, MPI::COMPLEX,
                                 root, tag);
        this->set_timestep(tstp, Gtemp);
    }
}
#endif

} // namespace cntr

#endif  // CNTR_HERM_MATRIX_IMPL_H
