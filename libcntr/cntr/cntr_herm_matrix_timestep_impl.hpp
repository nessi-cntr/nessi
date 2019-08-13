#ifndef CNTR_HERM_MATRIX_TIMESTEP_IMPL_H
#define CNTR_HERM_MATRIX_TIMESTEP_IMPL_H

#include "cntr_herm_matrix_timestep_decl.hpp"
//#include "cntr_exception.hpp"
#include "cntr_elements.hpp"
#include "cntr_function_decl.hpp"
#include "cntr_herm_matrix_decl.hpp"
#include "cntr_herm_matrix_timestep_view_impl.hpp"

namespace cntr {

  /* #######################################################################################
  #
  #   CONSTRUCTION/DESTRUCTION
  #
  ########################################################################################*/
template <typename T>
herm_matrix_timestep<T>::herm_matrix_timestep() {
    data_ = 0;
    ntau_ = 0;
    tstp_ = 0;
    size1_ = 0;
    size2_ = 0;
    element_size_ = 0;
    total_size_ = 0;
    sig_ = -1;
}
template <typename T>
herm_matrix_timestep<T>::~herm_matrix_timestep() {
    delete[] data_;
}

/** \brief <b> Initializes the `herm_matrix_timestep` class for a square-matrix for fermions. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `herm_matrix_timestep` class for a square-matrix. Parameter `sig_=-1`
* > is set for fermions. Obsolete, please use constructor with explicit 'sig'.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Time step
* @param ntau
* > Number of points on Matsubara axis
* @param size1
* > Matrix rank of the contour function. (size2=size1)
*/
template <typename T> herm_matrix_timestep<T>::herm_matrix_timestep(int tstp,int ntau,int size1){
   int len=((tstp+1)*2+(ntau+1))*size1*size1;
   assert(size1>=0 && tstp>=-1 && ntau>=0);
   if(len==0) data_=0;
   else{
      data_ = new cplx [len];
	  memset(data_, 0, sizeof(cplx)*len);
   }
   size1_=size1;
   size2_=size1;
   element_size_=size1*size1_;
   tstp_=tstp;
   ntau_=ntau;
   total_size_=len;
   sig_=-1;
}

/** \brief <b> Initializes the `herm_matrix_timestep` class for a general matrix. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `herm_matrix_timestep` class for a general matrix,
* > where the number of column and rows can be different.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Time step
* @param ntau
* > Number of points on Matsubara axis
* @param size1
* > Number of matrix rows
* @param size2
* > Number of matrix columns
* @param sig
* > Set `sig = -1` for fermions or `sig = +1` for bosons.
*/
template <typename T> herm_matrix_timestep<T>::herm_matrix_timestep(int tstp,int ntau,int size1,int size2,int sig){
   int len=((tstp+1)*2+(ntau+1))*size1*size2;
   assert(size1>=0 && tstp>=-1 && ntau>=0);
   if(len==0) data_=0;
   else{
      data_ = new cplx [len];
	  memset(data_, 0, sizeof(cplx)*len);
   }
   size1_=size1;
   size2_=size2;
   element_size_=size1*size2;
   tstp_=tstp;
   ntau_=ntau;
   total_size_=len;
   sig_=sig;
}

/** \brief <b> Initializes the `herm_matrix_timestep` class for a square-matrix for fermions/bosons. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `herm_matrix_timestep` class for a square-matrix.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Time step
* @param ntau
* > Number of points on Matsubara axis
* @param size1
* > Matrix rank of the contour function. (size2=size1)
@param sig
* > Set `sig = -1` for fermions or `sig = +1` for bosons.
*/
template <typename T> herm_matrix_timestep<T>::herm_matrix_timestep(int tstp,int ntau,int size1,int sig){
   int len=((tstp+1)*2+(ntau+1))*size1*size1;
   assert(size1>=0 && tstp>=-1 && ntau>=0 && sig*sig==1);
   if(len==0) data_=0;
   else{
      data_ = new cplx [len];
	  memset(data_, 0, sizeof(cplx)*len);
   }
   size1_=size1;
   size2_=size1;
   element_size_=size1*size1_;
   tstp_=tstp;
   ntau_=ntau;
   total_size_=len;
   sig_=sig;
}

/** \brief <b> Initializes the `herm_matrix_timestep` class with the same layout as a given `herm_matrix_timestep`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `herm_matrix_timestep` class with the same value of the time steps `tstp`,
* > number of points on the imaginary branch `ntau`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig`. Works for scalar or square-matrix contour objects only.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix` according to which the class should be initialized
*/
template <typename T>
herm_matrix_timestep<T>::herm_matrix_timestep(const herm_matrix_timestep &g) {
    tstp_ = g.tstp_;
    ntau_ = g.ntau_;
    size1_ = g.size1_;
    size2_ = g.size1_;
    element_size_ = size1_ * size1_;
    total_size_ = g.total_size_;
    if (total_size_ > 0) {
        data_ = new cplx[total_size_];
        memcpy(data_, g.data_, sizeof(cplx) * total_size_);
    } else {
        data_ = 0;
    }
    sig_ = g.sig_;
}
template <typename T>
herm_matrix_timestep<T> &herm_matrix_timestep<T>::
operator=(const herm_matrix_timestep &g) {
    if (this == &g)
        return *this;
    if (total_size_ != g.total_size_) {
        delete[] data_;
        total_size_ = g.total_size_;
        if (total_size_ > 0)
            data_ = new cplx[total_size_];
        else
            data_ = 0;
    }
    tstp_ = g.tstp_;
    ntau_ = g.ntau_;
    size1_ = g.size1_;
    size2_ = g.size1_;
    sig_ = g.sig_;
    element_size_ = size1_ * size1_;
    if (total_size_ > 0)
        memcpy(data_, g.data_, sizeof(cplx) * total_size_);
    return *this;
}
#if __cplusplus >= 201103L
template <typename T>
herm_matrix_timestep<T>::herm_matrix_timestep(
    herm_matrix_timestep &&g) noexcept : data_(g.data_),
                                         tstp_(g.tstp_),
                                         ntau_(g.ntau_),
                                         size1_(g.size1_),
                                         size2_(g.size2_),
                                         element_size_(g.element_size_),
                                         total_size_(g.total_size_),
                                         sig_(g.sig_) {
    g.data_ = nullptr;
    g.tstp_ = 0;
    g.ntau_ = 0;
    g.size1_ = 0;
    g.size2_ = 0;
    g.element_size_ = 0;
    g.total_size_ = 0;
}
template <typename T>
herm_matrix_timestep<T> &herm_matrix_timestep<T>::
operator=(herm_matrix_timestep &&g) noexcept {
    if (&g == this)
        return *this;

    data_ = g.data_;
    tstp_ = g.tstp_;
    ntau_ = g.ntau_;
    size1_ = g.size1_;
    size2_ = g.size2_;
    element_size_ = g.element_size_;
    total_size_ = g.total_size_;
    sig_ = g.sig_;

    g.data_ = nullptr;
    g.tstp_ = 0;
    g.ntau_ = 0;
    g.size1_ = 0;
    g.size2_ = 0;
    g.element_size_ = 0;
    g.total_size_ = 0;

    return *this;
}
#endif

/** \brief <b> Resizes `herm_matrix_timestep` object with respect to the number of points on the Matsubara branch and the matrix size at a given timestep </b>.
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Resizes `herm_matrix_timestep` class with respect to number
* > of points on the Matsubara branch `ntau` and the matrix size `size1`
* > at a given timestep `tstp`. Works for a square matrices.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > time steps
* @param ntau
* > number of points on Matsubara branch
* @param size1
* > size of the square matrix
*/
template <typename T>
void herm_matrix_timestep<T>::resize(int tstp, int ntau, int size1) {
    int len = ((tstp + 1) * 2 + (ntau + 1)) * size1 * size1;
    assert(ntau >= 0 && tstp >= -1 && size1 >= 0);
    delete[] data_;
    if (len == 0)
        data_ = 0;
    else {
        data_ = new cplx[len];
    }
    size1_ = size1;
    size2_ = size1;
    element_size_ = size1_ * size1_;
    tstp_ = tstp;
    ntau_ = ntau;
    total_size_ = len;
}

/** \brief <b> Sets all values to zero. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 *
 * > Sets all values to zero.
 *
 */
template <typename T>
void herm_matrix_timestep<T>::clear(void) {
    if (total_size_ > 0)
        memset(data_, 0, sizeof(cplx) * total_size_);
}
////////////////////////////////////////////////////////////////////////////////////////
// the following routines are not very "efficent" but sometimes simple to
// implement
#define herm_matrix_READ_ELEMENT                                             \
    {                                                                        \
        int r, s, dim = size1_;                                              \
        M.resize(dim, dim);                                                  \
        for (r = 0; r < dim; r++)                                            \
            for (s = 0; s < dim; s++)                                        \
                M(r, s) = x[r * dim + s];                                    \
    }
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

/** \brief <b> Returns the lesser component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the lesser component \f$ C^\mathrm{<}(t_i,tstp) \f$
* > at a given time \f$ t_i\f$ for particular time step `tstp`
* > to a given matrix class.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param M
* > Matrix to which the lesser component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix_timestep<T>::get_les_t_tstp(int i, Matrix &M) {
    cplx *x;
    x = lesptr(i);
    herm_matrix_READ_ELEMENT
}

/** \brief <b> Returns the lesser component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the lesser component \f$ C^\mathrm{<}(tstp, t_i) \f$
* > at a given time \f$ t_i\f$ for particular time step `tstp`
* > to a given matrix class.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param M
* > Matrix to which the lesser component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix_timestep<T>::get_les_tstp_t(int i, Matrix &M) {
    cplx *x;
    x = lesptr(i);
    herm_matrix_READ_ELEMENT_MINUS_CONJ
}

/** \brief <b> Returns the retarded component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the retarded component \f$ C^\mathrm{R}(tstp, t_j) \f$
* > at a given time \f$ t_j\f$ for particular time step `tstp`
* > to a given matrix class.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Index of time \f$ t_j\f$ .
* @param M
* > Matrix to which the retarded component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix_timestep<T>::get_ret_tstp_t(int j, Matrix &M) {
    cplx *x;
    x = retptr(j);
    herm_matrix_READ_ELEMENT
}

/** \brief <b> Returns the retarded component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the retarded component \f$ C^\mathrm{R}(t_i,tstp) \f$
* > at a given time \f$ t_i\f$ for particular time step `tstp`
* > to a given matrix class.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param M
* > Matrix to which the retarded component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix_timestep<T>::get_ret_t_tstp(int i, Matrix &M) {
    cplx *x;
    x = retptr(i);
    herm_matrix_READ_ELEMENT_MINUS_CONJ
}

/** \brief <b> Returns the left-mixing component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the left-mixing component \f$ C^\rceil(tstp,t_j) \f$
* > at a given time \f$ t_j\f$ for particular time step `tstp`
* > to a given matrix class.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Index of time \f$ t_j\f$ .
* @param M
* > Matrix to which the left-mixing component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix_timestep<T>::get_tv(int j, Matrix &M) {
    cplx *x = tvptr(j);
    herm_matrix_READ_ELEMENT
}

/** \brief <b> Returns the right-mixing component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the right-mixing component \f$ C^\lceil(t_i,t_j) \f$
* > at a given time \f$ t_i\f$ for particular time step `tstp`
* > to a given matrix class M. If 'sig' is negative, then M = -M
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param M
* > Matrix to which the right-mixing component is given.
* @param sig
* > Set `sig = -1` for fermions or `sig = +1` for bosons
*/
template <typename T>
template <class Matrix>
void herm_matrix_timestep<T>::get_vt(int i, Matrix &M, int sig) {
    cplx *x = tvptr(ntau_ - i);
    herm_matrix_READ_ELEMENT_MINUS_CONJ if (sig == -1) M = -M;
}

/** \brief <b> Returns the Matsubara component at given imaginary time.</b>
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
void herm_matrix_timestep<T>::get_mat(int i, Matrix &M) {
    cplx *x = matptr(i);
    herm_matrix_READ_ELEMENT
}

/** \brief <b> Returns the Matsubara component for the negative of a given imaginary time.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the Matsubara component \f$ C^\mathrm{M}(-\tau_i) \f$ at given
* > imaginary time \f$ \tau_i\f$ to a given matrix class M. If
* 'sig' is negative, return M =-M
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of imaginary time \f$ \tau_i\f$ .
* @param M
* > Matrix to which the Matsubara component is given.
* @param sig
* > Set `sig = -1` for fermions or `sig = +1` for bosons
*/
template <typename T>
template <class Matrix>
void herm_matrix_timestep<T>::get_matminus(int i, Matrix &M, int sig) {
    cplx *x = matptr(ntau_ - i);
    herm_matrix_READ_ELEMENT if (sig == -1) M = -M;
}

/** \brief <b> Returns the greater component at given times.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the greater component \f$ C^>(tstp,t_i) \f$ at given
* > time \f$t_i\f$ for particular time step `tstp`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$.
* @param M
* > Matrix to which the greater component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix_timestep<T>::get_gtr_tstp_t(int i, Matrix &M) {
    Matrix M1;
    get_ret_tstp_t(i, M);
    get_les_tstp_t(i, M1);
    M += M1;
}

/** \brief <b> Returns the greater component at given times.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the greater component \f$ C^>(t_i, tstp) \f$ at given
* > time \f$t_i\f$ for particular time step `tstp`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$.
* @param M
* > Matrix to which the greater component is given.
*/
template <typename T>
template <class Matrix>
void herm_matrix_timestep<T>::get_gtr_t_tstp(int i, Matrix &M) {
    Matrix M1;
    get_ret_t_tstp(i, M);
    get_les_t_tstp(i, M1);
    M += M1;
}

////////////////////////////////////////////////////////////////////////////////////////
// Note: same syntax like for herm_matrix (with dummy argument tstp).
// The following routines are not very "efficent" but sometimes simple to implement
////////////////////////////////////////////////////////////////////////////////////////

/** \brief <b> Returns the retarded component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the scalar-valued retarded component of a scalar type
* > \f$ C^\mathrm{R}(t_i,t_j) \f$ at a given times \f$ t_i\f$ and \f$ t_j\f$.
* > If \f$ t_i\f$ = `tstp`, \f$ C^\mathrm{R}(tstp,t_j) \f$ is returned.
* > Otherwise, one gets \f$ C^\mathrm{R}(t_i, tstp) \f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param j
* > Index of time \f$ t_j\f$ .
* @param x
* > complex number to which the retarded component is given.
*
* \note
* Same syntax like for herm_matrix (with dummy argument `tstp`)
*
*/
template <typename T>
inline void herm_matrix_timestep<T>::get_ret(int i, int j, cplx &x) {
  assert(i == tstp_ || j == tstp_);

    if (i == tstp_) {
        x = *retptr(j);
    } else {
        x = *retptr(i);
        x = -std::conj(x);
    }
}

/** \brief <b> Returns the lesser component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the scalar-valued lesser component of a scalar type
* > \f$ C^\mathrm{<}(t_i,t_j) \f$ at a given times \f$ t_i\f$ and \f$ t_j\f$.
* > If \f$ t_i\f$ = `tstp`, \f$ C^\mathrm{R}(tstp,t_j) \f$ is returned.
* > Otherwise, one gets \f$ C^\mathrm{R}(t_i, tstp) \f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param j
* > Index of time \f$ t_j\f$ .
* @param x
* > complex number to which the lesser component is given.
*
* \note
* Same syntax like for herm_matrix (with dummy argument `tstp`)
*
*/
template <typename T>
inline void herm_matrix_timestep<T>::get_les(int i, int j, cplx &x) {
  assert(i == tstp_ || j == tstp_);

    if (j == tstp_) {
        x = *lesptr(i);
    } else {
        x = *lesptr(j);
        x = -std::conj(x);
    }
}

/** \brief <b> Returns the left-mixing component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the scalar-valued left-mixing component
* > of a scalar type \f$ C^\rceil(tstp,t_j) \f$
* >  at a given time \f$ t_j\f$. Note, \f$ t_i\f$ is set to `tstp`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param j
* > Index of time \f$ t_j\f$ .
* @param x
* > complex number to which the left-mixing component is given.
*
* \note
* Same syntax like for herm_matrix (with dummy argument `tstp`)
*
*/
template <typename T>
inline void herm_matrix_timestep<T>::get_tv(int i, int j, cplx &x) {
  assert(i == tstp_);
     x = *tvptr(j);
}

/** \brief <b> Returns the right-mixing component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the scalar-valued right-mixing component
* > \f$ C^\lceil(t_i,t_j) \f$ at a given time \f$ \tau_i\f$. \f$ \tau_j\f$
* > is set to `tstp`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$ .
* @param j
* > Index of time \f$ t_j\f$ .
* @param x
* > complex number to which the right-mixing component is given.
*
* \note
* Same syntax like for herm_matrix (with dummy argument `tstp`)
*
*/
template <typename T>
inline void herm_matrix_timestep<T>::get_vt(int i, int j, cplx &x) {
  assert(j == tstp_);
     x = *tvptr(ntau_ - i);
    if (sig_ == -1)
        x = std::conj(x);
    else
        x = -std::conj(x);
}

/** \brief <b> Returns the matsubara component at a given imaginary time. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the scalar-valued Matsubara component \f$ C^\mathrm{M}(\tau_i) \f$
* > at given imaginary time \f$ \tau_i\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of imaginary time \f$ \tau_i\f$
* @param x
* > The value of the Matsubara component.
*/
template <typename T>
inline void herm_matrix_timestep<T>::get_mat(int i, cplx &x) {
  assert(tstp_ == -1);
    x = *matptr(i);
}

/** \brief <b> Returns the matsubara component for the negative of a given imaginary time. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the scalar-valued Matsubara component \f$ C^\mathrm{M}(-\tau_i) \f$ at given
* > imaginary time \f$ \tau_i\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of imaginary time \f$ \tau_i\f$
* @param x
* > The value of the Matsubara component.
*/
template <typename T>
inline void herm_matrix_timestep<T>::get_matminus(int i, cplx &x) {
  assert(tstp_ == -1);
    x = *matptr(ntau_ - i);
    if (sig_ == -1)
        x = -x;
}

/** \brief <b> Returns the greater component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the scalar-valued greater component \f$ C^>(t_i,t_j) \f$
* > at given times \f$t_i\f$ \f$t_j\f$.
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
inline void herm_matrix_timestep<T>::get_gtr(int i, int j, cplx &x) {
    cplx x1;
    get_ret(i, j, x);
    get_les(i, j, x1);
    x += x1;
}

/** \brief <b> Returns the density matrix at given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the scalar-valued density matrix (occupation, that is) at
* > given time step `tstp`. Setting `tstp = -1` returns the equilibrium density matrix
* > \f$ \rho = -C^\mathrm{M}(\beta) \f$, while `tstp >= 0` returns
* > \f$ \rho(t) = i \eta C^<(t,t) \f$.
* > The return value is formally complex.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > The time step at which the density matrix is returned.
*/
template <typename T>
std::complex<T> herm_matrix_timestep<T>::density_matrix(int tstp) {
  assert(tstp_ == tstp);
    cplx x1;
    if (tstp_ == -1) {
        get_mat(ntau_, x1);
        return -x1;
    } else {
        get_les(tstp_, tstp_, x1);
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
* > Returns the matrix-valued density matrix at
* > given time step `tstp`. Setting `tstp = -1` returns the equilibrium density matrix
* > \f$ \rho = -C^\mathrm{M}(\beta) \f$, while `tstp >= 0` returns \f$ \rho(t) = i \eta C^<(t,t) \f$.
* > Works for square-matrices only.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param M
* > The density matrix at time step `tstp`.
*/
template<typename T> template<class Matrix>
void herm_matrix_timestep<T>::density_matrix(Matrix &M){

   if(tstp_==-1){
     get_mat(ntau_,M);
     M *= (-1.0);
   }else{
     get_les_tstp_t(tstp_,M);
     M *= std::complex<T>(0.0,1.0*sig_);
   }
}


/* #######################################################################################
#
#   WRITING ELEMENTS FROM ANY MATRIX TYPE
#   OR FROM COMPLEX NUMBERS (then only the (0,0) element is addressed for dim>0)
#
########################################################################################*/
#define herm_matrix_SET_ELEMENT_MATRIX                                       \
    {                                                                        \
        int r, s;                                                             \
        for (r = 0; r < size1_; r++)                                            \
            for (s = 0; s < size2_; s++)                                        \
                x[r * size2_ + s] = M(r, s);                                    \
    }

/** \brief <b> Returns the retarded component into the reserved memory. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return the retarded component in `herm_matrix_timestep`
* > at a given time \f$ t_j\f$ into the reserved memory.
* > Works for scalar or general-matrix contour objects
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Index of time \f$ t_j\f$.
* @param M
* > Matrix in which the retarded component is given.
*/
template<typename T> template <class Matrix> void herm_matrix_timestep<T>::set_ret(int j,Matrix &M){
  cplx *x=retptr(j);
  herm_matrix_SET_ELEMENT_MATRIX
}

/** \brief <b> Returns the lesser component into the reserved memory. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return the lesser component in `herm_matrix_timestep`
* > at a given time \f$ t_j\f$ into the reserved memory.
* > Works for scalar or general-matrix contour objects
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Index of time \f$ t_j\f$.
* @param M
* > Matrix in which the lesser component is given.
*/
template<typename T> template <class Matrix> void herm_matrix_timestep<T>::set_les(int j,Matrix &M){
  cplx *x=lesptr(j);
  herm_matrix_SET_ELEMENT_MATRIX
}

/** \brief <b> Returns the left-mixing component into the reserved memory. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return the left-mixing component in `herm_matrix_timestep`
* > at a given time \f$ t_j\f$ into the reserved memory.
* > Works for scalar or general-matrix contour objects
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Index of time \f$ t_j\f$.
* @param M
* > Matrix in which the left-mixing component is given.
*/
template<typename T> template <class Matrix> void herm_matrix_timestep<T>::set_tv(int j,Matrix &M){
  cplx *x=tvptr(j);
  herm_matrix_SET_ELEMENT_MATRIX
}

/** \brief <b> Returns the Matsubara component into the reserved memory. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return the Matsubara component in `herm_matrix_timestep`
* > at a given time \f$ t_j\f$ into the reserved memory.
* > Works for scalar or general-matrix contour objects
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Index of time \f$ t_j\f$.
* @param M
* > Matrix in which the Matsubara component is given.
*/
template<typename T> template <class Matrix> void herm_matrix_timestep<T>::set_mat(int j,Matrix &M){
  cplx *x=matptr(j);
  herm_matrix_SET_ELEMENT_MATRIX
}




// G(t,t') ==> F(t)G(t,t')   ... ft+t*element_size_ points to F(t)

/** \brief <b> Left-multiplication of the `herm_matrix_timestep` with a contour function </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Left-multiplication of the `herm_matrix_timestep` with a time dependent contour function F(t)
* > i.e. it performs operation \f$G(t,t') \rightarrow w F(t)G(t,t')\f$
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param f0
* > pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix).
* @param ft
* > the contour function F(t)
* @param weight
* > some number (weight)
*/
template <typename T>
void herm_matrix_timestep<T>::left_multiply(std::complex<T> *f0,
                                            std::complex<T> *ft, T weight) {
    int m;
    cplx *x0, *xtemp, *ftemp;
    xtemp = new cplx[element_size_];
    if (tstp_ == -1) {
        x0 = data_;
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, f0,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    } else {
        ftemp = ft + tstp_ * element_size_;
        x0 = data_;
        for (m = 0; m <= tstp_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, ftemp,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = data_ + (tstp_ + 1) * element_size_;
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, ftemp,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = data_ + (tstp_ + 1 + ntau_ + 1) * element_size_;
        for (m = 0; m <= tstp_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, ft + m * element_size_,
                                       x0 + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    }
    delete[] xtemp;
}

/** \brief <b> Left-multiplication of the `herm_matrix_timestep` with a contour function </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Right-multiplication of the `herm_matrix_timestep` with a time dependent contour function F(t)
* > i.e. it performs operation \f$G(t,t') \rightarrow w F(t)G(t,t')\f$
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param ft
* > the contour function F(t)
* @param weight
* > some number (weight)
*/
template <typename T>
void herm_matrix_timestep<T>::left_multiply(function<T> &ft, T weight) {
  assert( ft.nt() >= tstp_);

    this->left_multiply(ft.ptr(-1), ft.ptr(0), weight);
}

// G(t,t') ==> G(t,t')F(t')

/** \brief <b> Right-multiplication of the `herm_matrix_timestep` with a contour function </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Right-multiplication of the `herm_matrix_timestep` with a time dependent contour function F(t)
* > i.e. it performs operation \f$G(t,t') \rightarrow w G(t,t')F(t')\f$
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param f0
* > pointer to \f$F(-\mathrm{i}\beta)\f$ on the Matsubara axis (this is just a constant matrix).
* @param ft
* > the contour function F(t)
* @param weight
* > some number (weight)
*/
template <typename T>
void herm_matrix_timestep<T>::right_multiply(std::complex<T> *f0,
                                             std::complex<T> *ft, T weight) {
    int m;
    cplx *xtemp, *ftemp, *x0;
    xtemp = new cplx[element_size_];
    if (tstp_ == -1) {
        x0 = data_;
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       f0);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    } else {
        x0 = data_;
        for (m = 0; m <= tstp_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       ft + m * element_size_);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        x0 = data_ + (tstp_ + 1) * element_size_;
        for (m = 0; m <= ntau_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       f0);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
        ftemp = ft + tstp_ * element_size_;
        x0 = data_ + (tstp_ + 1 + ntau_ + 1) * element_size_;
        for (m = 0; m <= tstp_; m++) {
            element_mult<T, LARGESIZE>(size1_, xtemp, x0 + m * element_size_,
                                       ftemp);
            element_smul<T, LARGESIZE>(size1_, xtemp, weight);
            element_set<T, LARGESIZE>(size1_, x0 + m * element_size_, xtemp);
        }
    }
    delete[] xtemp;
}

/** \brief <b> Right-multiplication of the `herm_matrix_timestep` with a contour function </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Right-multiplication of the `herm_matrix_timestep` with a time dependent contour function F(t)
* > i.e. it performs operation \f$G(t,t') \rightarrow w G(t,t')F(t')\f$
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param ft
* > the contour function F(t)
* @param weight
* > some number (weight)
*/
template <typename T>
void herm_matrix_timestep<T>::right_multiply(function<T> &ft, T weight) {
  assert(ft.nt() >= tstp_);

    this->right_multiply(ft.ptr(-1), ft.ptr(0), weight);
}

/** \brief <b> Increase the value of the  `herm_matrix_timestep` by `weight` \f$ * g(t) \f$ </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Increase the value of `herm_matrix_timestep` by a value of `weight`*\f$ g(t)\f$,
* > where \f$ g(t)\f$ is a `herm_matrix_timestep` and `weight` is a constant scalar number
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param g1
 * > The `herm_matrix_timestep` which is added
 * @param weight
 * > Scalar constant multiplication factor
 *
 */
template <typename T>
void herm_matrix_timestep<T>::incr(herm_matrix_timestep<T> &g1, T weight) {
    int m;
    assert(g1.size1_ == size1_ && g1.ntau_ == ntau_ && g1.tstp_ == tstp_);
    for (m = 0; m < total_size_; m++)
        data_[m] += weight * g1.data_[m];
}
#define HERM_MATRIX_INCR_TSTP                                                \
    if (alpha == cplx(1.0, 0.0)) {                                           \
        for (i = 0; i < len; i++)                                            \
            x0[i] += x[i];                                                   \
    } else {                                                                 \
        for (i = 0; i < len; i++)                                            \
            x0[i] += alpha * x[i];                                           \
    }

/** \brief <b> Increase the value of the `herm_matrix_timestep` by \f$\alpha g(t)\f$ </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Increase the value of `herm_matrix_timestep` by a value of \f$\alpha g(t) \f$, where \f$ g(t)\f$ is a `herm_matrix` and
* > \f$\alpha \f$ is a complex number. If \f$t>-1\f$ then `ret,les,tv` components are set, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param g
 * > The `herm_matrix` which is added
 * @param alpha
 * > Constant multiplication factor
 *
 */
template <typename T>
void herm_matrix_timestep<T>::incr(herm_matrix<T> &g, T alpha) {
    int i, len;
    cplx *x, *x0;
    assert(tstp_ <= g.nt() && ntau_ == g.ntau() && size1_ == g.size1());
    if (tstp_ == -1) {
        len = (ntau_ + 1) * element_size_;
        x = g.matptr(0);
        x0 = data_;
        HERM_MATRIX_INCR_TSTP
    } else {
        len = (tstp_ + 1) * element_size_;
        x = g.retptr(tstp_, 0);
        x0 = data_;
        HERM_MATRIX_INCR_TSTP
        len = (ntau_ + 1) * element_size_;
        x = g.tvptr(tstp_, 0);
        x0 = data_ + (tstp_ + 1) * element_size_;
        HERM_MATRIX_INCR_TSTP
        len = (tstp_ + 1) * element_size_;
        x = g.lesptr(0, tstp_);
        x0 = data_ + (tstp_ + 1 + ntau_ + 1) * element_size_;
        HERM_MATRIX_INCR_TSTP
    }
}
#undef HERM_MATRIX_INCR_TSTP

/** \brief <b> Multiply  `herm_matrix_timestep` with scalar `weight`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Multiply `herm_matrix_timestep` with a scalar.
* > If \f$t>-1\f$ then `ret,les,tv` components are set, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param weight
* > The `template argument` multiplication factor
*/

template <typename T>
void herm_matrix_timestep<T>::smul(T weight) {
    herm_matrix_timestep_view<T> g_view(*this);

    g_view.smul(weight);

}



/** \brief <b> Gets matrix elements of all components at time step `tstp` from the components of a given `herm_matrix` of a scaler type. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Gets the (sub-) matrix elements `(i1,i2)` of all components at time step `tstp` from the components of
 * >  a given `herm_matrix` of scaler type.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param i1
 * > First index for submatrix.
 * @param i2
 * > Second index for submatrix.
 * @param g
 * > The `herm_matrix` from which the time step is copied.
 *
 */
template <typename T>
void herm_matrix_timestep<T>::get_matrixelement(int i1, int i2,
                                                herm_matrix<T> &g) {
    int i, sij = i1 * g.size1() + i2;
    cplx *x;
    // std::cout << tstp_ << " " << g.nt() << " "  << i1 << " " << i2 << " "
    // << g.size1() << " " << g.size2() <<  std::endl;
    assert(tstp_ <= g.nt() && 0 <= i1 && i1 < g.size1() && 0 <= i2 &&
           i2 < g.size1() && size1_ == 1);
    if (tstp_ == -1) {
        x = data_;
        for (i = 0; i <= ntau_; i++)
            x[i] = *(g.matptr(i) + sij);
    } else {
        x = data_;
        for (i = 0; i <= tstp_; i++)
            x[i] = *(g.retptr(tstp_, i) + sij);
        x = data_ + (tstp_ + 1) * element_size_;
        for (i = 0; i <= ntau_; i++)
            x[i] = *(g.tvptr(tstp_, i) + sij);
        x = data_ + (tstp_ + 1 + ntau_ + 1) * element_size_;
        for (i = 0; i <= tstp_; i++)
            x[i] = *(g.lesptr(i, tstp_) + sij);
    }
}

/** \brief <b> Set the matrix element of \f$C_{[i1,i2]}\f$ for each component from  \f$ g_{[i1,i2]} \f$ </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set the matrix element \f$C_{[i1,i2]}\f$ for each component from  \f$ g_{[i1,i2]} \f$.
* > Use implementation of `herm_matrix_timestep_view`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i1
* > Row index of `herm_matrix_timestep`
* @param i2
* > Column index of `herm_matrix_timestep`
* @param g
* > The `herm_matrix_timestep` from which the matrix element element is given
* @param j1
* > Row index of `herm_matrix_timestep`
* @param j2
* > Column index of `herm_matrix_timestep`
*/
template <typename T>
void herm_matrix_timestep<T>::set_matrixelement(int i1, int i2,
                                                herm_matrix_timestep<T> &g,
                                                int j1, int j2) {
    herm_matrix_timestep_view<T> tmp(g);
    herm_matrix_timestep_view<T> tmp1(*this);
    tmp1.set_matrixelement(i1, i2, tmp, j1, j2);
}

/** \brief <b> Set the matrix element of \f$C_{[i1,i2]}\f$ for each component from \f$ g_{[i1,i2]} \f$ </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set the matrix element \f$C_{[i1,i2]}\f$ for each component from  \f$ g_{[i1,i2]} \f$,
* > which is given by `herm_matrix_timestep_view`.
* > Use implementation of `herm_matrix_timestep_view`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i1
* > Row index of `herm_matrix_timestep`
* @param i2
* > Column index of `herm_matrix_timestep`
* @param g
* > The `herm_matrix_timestep_view` from which the matrix element element is given
* @param j1
* > Row index of `herm_matrix_timestep`
* @param j2
* > Column index of `herm_matrix_timestep`
*/
template <typename T>
void herm_matrix_timestep<T>::set_matrixelement(
    int i1, int i2, herm_matrix_timestep_view<T> &g, int j1, int j2) {
    herm_matrix_timestep_view<T> tmp1(*this);
    tmp1.set_matrixelement(i1, i2, g, j1, j2);
}

/** \brief <b> Set the matrix element of \f$C_{[i1,i2]}\f$ for each component from \f$ g_{[i1,i2]} \f$ </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set the matrix element \f$C_{[i1,i2]}\f$ for each component from \f$ g_{[i1,i2]}\f$,
* > which is given by `herm_matrix`.
* > Use implementation of `herm_matrix_timestep_view` at each time step `tstp`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i1
* > Row index of `herm_matrix_timestep`
* @param i2
* > Column index of `herm_matrix_timestep`
* @param g
* > The `herm_matrix` from which the matrix element element is given
* @param j1
* > Row index of `herm_matrix_timestep`
* @param j2
* > Column index of `herm_matrix_timestep`
*/
template <typename T>
void herm_matrix_timestep<T>::set_matrixelement(int i1, int i2,
                                                herm_matrix<T> &g, int j1,
                                                int j2) {
    herm_matrix_timestep_view<T> tmp(tstp_, g);
    herm_matrix_timestep_view<T> tmp1(*this);
    tmp1.set_matrixelement(i1, i2, tmp, j1, j2);
}

#if CNTR_USE_MPI == 1

/** \brief <b> MPI reduce for the `herm_matrix_timestep` </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > MPI reduce for the `herm_matrix_timestep` to the `root`
* > Works for scalar or square-matrix contour objects.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param root
* > Index of root
*/
template <typename T>
void herm_matrix_timestep<T>::MPI_Reduce(int root) {
    int taskid;
    int len1 = 2 * ((tstp_ + 1) * 2 + (ntau_ + 1)) * size1_ * size1_;
    taskid = MPI::COMM_WORLD.Get_rank();
    if (sizeof(T) == sizeof(double)) {
        if (taskid == root) {
            MPI::COMM_WORLD.Reduce(MPI::IN_PLACE, (double *)this->data_, len1,
                                   MPI::DOUBLE, MPI::SUM, root);
        } else {
            MPI::COMM_WORLD.Reduce((double *)this->data_,
                                   (double *)this->data_, len1, MPI::DOUBLE,
                                   MPI::SUM, root);
        }
    } else {
        std::cerr << "herm_matrix_timestep<T>::MPI_Reduce only for double "
                  << std::endl;
        exit(0);
    }
}
#endif

/* #########################################################################
#
#   MPI UTILS
#
############################################################################ */


#if CNTR_USE_MPI == 1

/** \brief <b> Broadcasts the `herm_matrix` at a given time step to all tasks. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Broadcasts the `herm_matrix_timestep` at a given time step `tstp` to all tasks.
* > Works for a square matrices.
*
*<!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Time step which should be broadcasted.
* @param ntau
* > number of point on the Matsubara branch
* @param size1
* > size of the matrix
* @param root
* > The task rank from which the `herm_matrix` should be broadcasted.
*
* \note
* The green's function is resized before broadcast with respect to number of points on the Matsubara branch `ntau`
* and the matrix size `size1` at a given timestep `tstp`.
*/
template <typename T>
void herm_matrix_timestep<T>::Bcast_timestep(int tstp, int ntau, int size1,
                                             int root) {
    int numtasks = MPI::COMM_WORLD.Get_size();
    int taskid = MPI::COMM_WORLD.Get_rank();
    if (taskid != root)
        resize(tstp, ntau, size1);
    int len = (2 * (tstp_ + 1) + ntau_ + 1) * element_size_;
    // test effective on root:
    assert(tstp == tstp_);
    assert(ntau == ntau_);
    assert(size1 == size1_);
    if (sizeof(T) == sizeof(double))
        MPI::COMM_WORLD.Bcast(data_, len, MPI::DOUBLE_COMPLEX, root);
    else
        MPI::COMM_WORLD.Bcast(data_, len, MPI::COMPLEX, root);
}

/** \brief <b> Sends the `herm_matrix_timestep` at a given time step to a specific task. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Sends the `herm_matrix_timestep` at a given time step `tstp` to a specific task with rank `dest`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Time step which should be send.
* @param ntau
* > number of points on the Matsubara branch
* @param size1
* > size of the matrix
* @param dest
* > The task rank to which the `herm_matrix` should be send.
* @param tag
* > The MPI error flag.
*/
template <typename T>
void herm_matrix_timestep<T>::Send_timestep(int tstp, int ntau, int size1,
                                            int dest, int tag) {
    int taskid = MPI::COMM_WORLD.Get_rank();
    int len = (2 * (tstp_ + 1) + ntau_ + 1) * element_size_;
    if (!(taskid == dest)) {
      assert(tstp == tstp_);
      assert(ntau == ntau_);
      assert(size1 == size1_);

        // std::cout << "SEnd timestep tstp= " << tstp << " ntau=  " << ntau
        //<< " size1= " << size1 << " rank= " << taskid << " dest= " << dest<<
        //" tot= " << les << " tag " << tag << std::endl;
        if (sizeof(T) == sizeof(double))
            MPI::COMM_WORLD.Send(data_, len, MPI::DOUBLE_COMPLEX, dest, tag);
        else
            MPI::COMM_WORLD.Send(data_, len, MPI::COMPLEX, dest, tag);
    } else {
        // std::cout << "Sending timestep  DEST==ROOT tstp= " << tstp << "
        // ntau=  " << ntau
        //<< " size1= " << size1 << " root= " << taskid << " dest= " << dest
        //<< std::endl;
    }
}

/** \brief <b> Recevies the `herm_matrix` at a given time step from a specific task. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Receives the `herm_matrix_timestep` at a given time step `tstp` from a specific task with rank `root`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > Time step which should be received.
 * @param ntau
 * > number of points on the Matsubara branch
 * @param size1
 * > size of the matrix
 * @param root
 * > The task rank from which the `herm_matrix` should be received.
 * @param tag
 * > The MPI error flag.
 */
template <typename T>
void herm_matrix_timestep<T>::Recv_timestep(int tstp, int ntau, int size1,
                                            int root, int tag) {
    int taskid = MPI::COMM_WORLD.Get_rank();
    if (!(taskid == root)) {
        resize(tstp, ntau, size1);
        int len = (2 * (tstp_ + 1) + ntau_ + 1) * element_size_;
        // std::cout << "Rcv timestep tstp= " << tstp << " ntau=  " << ntau
        //<< " size1= " << size1 << " rank= " << taskid << " root= " << root<<
        //" tot= " << les << " tag " << tag << std::endl;
        if (sizeof(T) == sizeof(double))
            MPI::COMM_WORLD.Recv(data_, len, MPI::DOUBLE_COMPLEX, root, tag);
        else
            MPI::COMM_WORLD.Recv(data_, len, MPI::COMPLEX, root, tag);
        // std::cout << "FINISHED Rcv timestep tstp= " << tstp << " ntau=  "
        // << ntau
        //<< " size1= " << size1 << " rank= " << taskid << " root= " << root<<
        //" tot= " << total_size_ << " tag " << tag << std::endl;
    }
}
#endif

#undef herm_matrix_READ_ELEMENT
#undef herm_matrix_READ_ELEMENT_MINUS_CONJ

/* ##########################################################################
#
#   HDF5 I/O (using implementation in herm_matrix_timestep_view)
#
#############################################################################*/

#if CNTR_USE_HDF5 == 1

/** \brief <b> Write `herm_matrix_timestep` to hdf5 group given by `group_id` </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Write `herm_matrix_timestep` to the hdf5 group given by `group_id`.
* > Use implementation given in `herm_matrix_timestep_view`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param group_id
* > Group ID of the hdf5 file to write in
*/
template <typename T>
void herm_matrix_timestep<T>::write_to_hdf5(hid_t group_id) {
    herm_matrix_timestep_view<T> tmp(*this);
    tmp.write_to_hdf5(group_id);
}

/** \brief <b> Write `herm_matrix_timestep` to hdf5 group given by `group_id` and name `groupname` </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Write `herm_matrix_timestep` to the hdf5 group given by `group_id` and `groupname`.
* > Use implementation given in `herm_matrix_timestep_view`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param group_id
* > hdf5 group ID for group to write in
* @param groupname
* > hdf5 group name for group to write in
*/
template <typename T>
void herm_matrix_timestep<T>::write_to_hdf5(hid_t group_id,
                                            const char *groupname) {
    herm_matrix_timestep_view<T> tmp(*this);
    tmp.write_to_hdf5(group_id, groupname);
}

/** \brief <b> Write `herm_matrix_timestep` to hdf5 group given by name `groupname` and store in file `filename` </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Write `herm_matrix_timestep` to the hdf5 group given by `groupname` and store in file `filename`.
* > Use implementation given in `herm_matrix_timestep_view`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param filename
* > hdf5 filename for file to write in
* @param groupname
* > hdf5 group name for group to write in
*/
template <typename T>
void herm_matrix_timestep<T>::write_to_hdf5(const char *filename,
                                            const char *groupname) {
    herm_matrix_timestep_view<T> tmp(*this);
    tmp.write_to_hdf5(filename, groupname);
}

/** \brief <b> Read `herm_matrix_timestep` from hdf5 group given by `group_id`</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read `herm_matrix_timestep` from the hdf5 group given by `group_id`.
* > Note, the green's function is resized before implementing a routine from `herm_matrix_timestep_view`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
@param group_id
* > Group ID of the hdf5 file to read from
*/
template <typename T>
void herm_matrix_timestep<T>::read_from_hdf5(hid_t group_id) {
    // -- implementation is different from herm_matrix_timestep_view, because G can be resized:
    int tstp = read_primitive_type<int>(group_id, "tstp");
    int ntau = read_primitive_type<int>(group_id, "ntau");
    int sig = read_primitive_type<int>(group_id, "sig");
    int size1 = read_primitive_type<int>(group_id, "size1");
    // RESIZE G
    this->resize(tstp, ntau, size1);
    sig_ = sig;
    herm_matrix_timestep_view<T> tmp(*this);
    tmp.read_from_hdf5(group_id);
}

/** \brief <b> Read `herm_matrix_timestep` from hdf5 group given by `group_id` and group name `groupname`</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read `herm_matrix_timestep` from the hdf5 group given by `group_id` and group name `groupname`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
@param group_id
* > Group ID of the hdf5 file to read from
* @param groupname
* > hdf5 group name for group to read from
*/
template <typename T>
void herm_matrix_timestep<T>::read_from_hdf5(hid_t group_id,
                                             const char *groupname) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(sub_group_id);
    close_group(sub_group_id);
}

/** \brief <b> Read `herm_matrix_timestep` from hdf5 group given by group name `groupname` in file `filename`</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read `herm_matrix_timestep` from the hdf5 group given by group name `groupname` in file `filename`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
@param filename
* > hdf5 filename for file to read from
* @param groupname
* > hdf5 group name for group to read from
*/
template <typename T>
void herm_matrix_timestep<T>::read_from_hdf5(const char *filename,
                                             const char *groupname) {
    hid_t file_id = read_hdf5_file(filename);
    this->read_from_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
}
#endif

} // namespace cntr

#endif  // CNTR_HERM_MATRIX_TIMESTEP_IMPL_H
