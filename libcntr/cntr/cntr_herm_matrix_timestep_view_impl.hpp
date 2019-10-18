#ifndef CNTR_HERM_TIMESTEP_VIEW_IMPL_H
#define CNTR_HERM_TIMESTEP_VIEW_IMPL_H

#include "cntr_herm_matrix_timestep_view_decl.hpp"
#include "cntr_elements.hpp"
//#include "cntr_exception.hpp"

namespace cntr {

#define CPLX std::complex<T>
template <typename T>
/* #######################################################################################
#
#   CONSTRUCTION/DESTRUCTION
#
########################################################################################*/
herm_matrix_timestep_view<T>::herm_matrix_timestep_view() {
    ret_ = 0;
    les_ = 0;
    tv_ = 0;
    mat_ = 0;
    ntau_ = 0;
    tstp_ = -2;
    size1_ = 0;
    size2_ = 0;
    element_size_ = 0;
    sig_ = -1;
}
template <typename T>
herm_matrix_timestep_view<T>::~herm_matrix_timestep_view() {
    // donothing
}

/** \brief <b> Initializes the `herm_matrix_timestep_view` class with the same layout as a given `herm_matrix_timestep_view g`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and copy the `herm_matrix_timestep_view` class with the same number of time steps `nt`,
* > number of points on the imaginary branch `ntau`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig`. Works for scalar or square-matrix contour objects
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix_timestep_view` according to which the class should be initialized
*/

template <typename T>
herm_matrix_timestep_view<T>::herm_matrix_timestep_view(
    const herm_matrix_timestep_view &g) {
    tstp_ = g.tstp_;
    ntau_ = g.ntau_;
    size1_ = g.size1_;
    size2_ = g.size1_;
    element_size_ = size1_ * size1_;
    sig_ = g.sig_;
    // copy pointers
    ret_ = g.ret_;
    les_ = g.les_;
    tv_ = g.tv_;
    mat_ = g.mat_;
}



/** \brief <b> Initializes the `herm_matrix_timestep_view` class for a general matrix. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `herm_matrix_timestep_view` class for a general matrix,
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
template <typename T> herm_matrix_timestep_view<T>::herm_matrix_timestep_view(int tstp,int ntau,int size1,int size2,int sig){

   assert(size1>=0 && tstp>=-1 && ntau>=0);
   ret_ = 0;
   les_ = 0;
   tv_ = 0;
   mat_ = 0;
   size1_=size1;
   size2_=size2;
   element_size_=size1*size2;
   tstp_=tstp;
   ntau_=ntau;
   sig_=sig;

}



/** \brief <b> Copy assignment operator for `herm_matrix_timestep_view g` class. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and copy the `herm_matrix_timestep_view` class with the same number of time steps `nt`,
* > number of points on the imaginary branch `ntau`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig`. Works for scalar or square-matrix contour objects
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix_timestep_view` according to which the class should be initialized
*/

template <typename T>
herm_matrix_timestep_view<T> &herm_matrix_timestep_view<T>::
operator=(const herm_matrix_timestep_view &g) {
    if (this == &g)
        return *this;
    tstp_ = g.tstp_;
    ntau_ = g.ntau_;
    size1_ = g.size1_;
    size2_ = g.size1_;
    sig_ = g.sig_;
    element_size_ = size1_ * size1_;
    // copy pointers
    ret_ = g.ret_;
    les_ = g.les_;
    tv_ = g.tv_;
    mat_ = g.mat_;
    return *this;
}

/** \brief <b> Initializes the `herm_matrix_timestep_view` class with the same layout as a given `herm_matrix_timestep g`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and copy the `herm_matrix_timestep_view` class from the `herm_matrix_timestep`  with the same number of time steps `nt`,
* > number of points on the imaginary branch `ntau`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig`. No data, only pointers are copied.
* > Works for scalar or square-matrix contour objects
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix_timestep` according to which the class should be initialized
*/

template <typename T>
herm_matrix_timestep_view<T>::herm_matrix_timestep_view(
    herm_matrix_timestep<T> &g) {
    tstp_ = g.tstp_;
    ntau_ = g.ntau_;
    size1_ = g.size1_;
    size2_ = g.size2_;
    element_size_ = size1_ * size2_;
    sig_ = g.sig_;
    if (tstp_ == -1) {
        mat_ = g.matptr(0);
        ret_ = 0;
        les_ = 0;
        tv_ = 0;
    } else if (tstp_ >= 0) {
        mat_ = 0;
        ret_ = g.retptr(0);
        les_ = g.lesptr(0);
        tv_ = g.tvptr(0);
    }
}

/** \brief <b> Initializes the `herm_matrix_timestep_view` class with the same layout as a given `herm_matrix_timestep g`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and copy the `herm_matrix_timestep_view` class from the `herm_matrix_timestep`  with the same number of time steps `nt`,
* > number of points on the imaginary branch `ntau`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig`. No data, only pointers are copied. First argument `tstp` is redundant, but present due to the safety [and historical] reason.
* > Works for scalar or square-matrix contour objects
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Index of time
* @param g
* > The `herm_matrix_timestep` according to which the class should be initialized
*/

template <typename T>
herm_matrix_timestep_view<T>::herm_matrix_timestep_view(
    int tstp, herm_matrix_timestep<T> &g) {
    assert(tstp==g.tstp_);

    tstp_ = g.tstp_;
    ntau_ = g.ntau_;
    size1_ = g.size1_;
    size2_ = g.size2_;
    element_size_ = size1_ * size2_;
    sig_ = g.sig_;
    if (tstp_ == -1) {
        mat_ = g.matptr(0);
        ret_ = 0;
        les_ = 0;
        tv_ = 0;
    } else if (tstp_ >= 0) {
        mat_ = 0;
        ret_ = g.retptr(0);
        les_ = g.lesptr(0);
        tv_ = g.tvptr(0);
    }
}

/** \brief <b> Initializes the `herm_matrix_timestep_view` class with the same layout as a given `herm_matrix_timestep_view g`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and copy the `herm_matrix_timestep_view` class from the `herm_matrix_timestep`  with the same number of time steps `nt`,
* > number of points on the imaginary branch `ntau`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig`. No data, only pointers are copied. First argument `tstp` is redundant, but present due to the safety [and historical] reason.
* > Works for scalar or square-matrix contour objects
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Index of time
* @param g
* > The `herm_matrix_timestep` according to which the class should be initialized
*/
template <typename T>
herm_matrix_timestep_view<T>::herm_matrix_timestep_view(
    int tstp, herm_matrix_timestep_view<T> &g) {
    assert(tstp == g.tstp_);

    tstp_ = g.tstp_;
    ntau_ = g.ntau_;
    size1_ = g.size1_;
    size2_ = g.size1_;
    element_size_ = size1_ * size1_;
    sig_ = g.sig_;
    // copy pointers
    ret_ = g.ret_;
    les_ = g.les_;
    tv_ = g.tv_;
    mat_ = g.mat_;
}

/** \brief <b> Initializes the `herm_matrix_timestep_view` class with the same layout as a given `herm_matrix g`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and copy the `herm_matrix_timestep_view` class from the `herm_matrix`  with the same number of time steps `nt`,
* > number of points on the imaginary branch `ntau`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig` for time slice `tstp`. No data, only pointers are copied.
* > Works for scalar or square-matrix contour objects
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Index of time
* @param g
* > The `herm_matrix` according to which the class should be initialized
*/

template <typename T>
herm_matrix_timestep_view<T>::herm_matrix_timestep_view(int tstp,
                                                        herm_matrix<T> &g) {
    assert(tstp>=-1 && tstp <=g.nt());
    tstp_ = tstp;
    ntau_ = g.ntau();
    size1_ = g.size1();
    size2_ = g.size2();
    element_size_ = size1_ * size2_;
    sig_ = g.sig();
    if (tstp_ == -1) {
        mat_ = g.matptr(0);
        ret_ = 0;
        les_ = 0;
        tv_ = 0;
    } else if (tstp_ >= 0) {
        mat_ = 0;
        ret_ = g.retptr(tstp_, 0);
        les_ = g.lesptr(0, tstp_);
        tv_ = g.tvptr(tstp_, 0);
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

/// @private
/** \brief <b> Sets the retarded component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Sets \f$G^\mathrm{R}(t_i, t_j)\f$ to given matrix `M`.
* > Restricted to the domain of `herm_matrix`: \f$i \ge j\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$.
* @param j
* > Index of time \f$ t_j\f$.
* @param M
* > Matrix in which the retarded component is given.
*/
template<typename T> template <class Matrix> void herm_matrix_timestep_view<T>::set_ret(int i, int j,Matrix &M){
         assert(i == tstp_);
         assert(j <= i);
         cplx *x=retptr(j);
         herm_matrix_SET_ELEMENT_MATRIX
      }

/// @private
/** \brief <b> Sets the lesser component at given times. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Sets the lesser component \f$G^<(t_i,t_j)\f$ of `herm_matrix_timestep`
* > to a given complex matrix `M`. Restricted to the domain of `herm_matrix`: \f$j \ge i\f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$.
* @param j
* > Index of time \f$ t_j\f$.
* @param M
* > Matrix in which the lesser component is given.
*/
template<typename T> template <class Matrix> void herm_matrix_timestep_view<T>::set_les(int i, int j, Matrix &M){
      assert(j == tstp_);
      assert(i <= j);
      cplx *x=lesptr(i);
      herm_matrix_SET_ELEMENT_MATRIX
   }

/// @private
/** \brief <b> Sets the left-mixing component to a given matrix. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Sets the lesser component \f$G^\rceil(t_i,\tau_j)\f$ of `herm_matrix_timestep`
* > to a given complex matrix `M`. 
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of time \f$ t_i\f$.
* @param j
* > Index of imaginary time \f$ \tau_j\f$.
* @param M
* > Matrix in which the left-mixing component is given.
*/
template<typename T> template <class Matrix> void herm_matrix_timestep_view<T>::set_tv(int i, int j, Matrix &M){
   assert(i == tstp_);
   assert(j <= ntau_);
   cplx *x=tvptr(j);
   herm_matrix_SET_ELEMENT_MATRIX
}

/// @private
/** \brief <b> Sets the Matsubara component to a given matrix. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Sets the Matsubara component \f$G^\mathrm{M}(\tau_j)\f$ of `herm_matrix_timestep`
* > to a given complex matrix `M`. 
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > Index of imaginary time \f$ \tau_i\f$.
* @param M
* > Matrix in which the Matsubara component is given.
*/
template<typename T> template <class Matrix> void herm_matrix_timestep_view<T>::set_mat(int i,Matrix &M){
   assert(i <= ntau_ );
   cplx *x=matptr(i);
   herm_matrix_SET_ELEMENT_MATRIX
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

/// @private
/** \brief <b> Returns the lesser component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the lesser component \f$ C^\mathrm{<}(t_j,t_i) \f$
* > at a given time \f$ t_j\f$ with \f$ j <= tstp\f$ for particular time step `i=tstp`
* > to a given matrix class.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Index of time \f$ t_j\f$ .
* @param M
* > Matrix to which the lesser component is given.
*/
template <typename T> template <class Matrix> 
      void herm_matrix_timestep_view<T>::get_les(int i, int j, Matrix &M){
         assert(j == tstp_ && i <= tstp_);
         cplx *x;
         x = lesptr(i);
         herm_matrix_READ_ELEMENT
      }

/// @private
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
      void herm_matrix_timestep_view<T>::get_ret(int i, int j, Matrix &M) {
         assert(i == tstp_ && j <= tstp_);
         cplx *x;
         x = retptr(j);
         herm_matrix_READ_ELEMENT
      }

/// @private
/** \brief <b> Returns the left-mixing component at a given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the left-mixing component \f$ C^\rceil(tstp,\tau_j) \f$
* > at a given time \f$ \tau_j\f$ for particular time step `tstp`
* > to a given matrix class.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Index of time \f$ \tau_j\f$ .
* @param M
* > Matrix to which the left-mixing component is given.
*/
template <typename T>
template <class Matrix>
      void herm_matrix_timestep_view<T>::get_tv(int i, int j, Matrix &M) {
         assert(i == tstp_);
         cplx *x = tvptr(j);
         herm_matrix_READ_ELEMENT
      }

/// @private
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
      void herm_matrix_timestep_view<T>::get_mat(int i, Matrix &M) {
         cplx *x = matptr(i);
         herm_matrix_READ_ELEMENT
      }

/// @private
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
      void herm_matrix_timestep_view<T>::get_matminus(int i, Matrix &M) {
         cplx *x = matptr(ntau_ - i);
         herm_matrix_READ_ELEMENT if (sig_ == -1) M = -M;
      }


/** \brief <b> Sets all components at time step `tstp` to zero. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets all components of the `herm_matrix_timestep_view` to zero. If
 * > `tstp = -1`, only the Matsubara component will be set to zero.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > The time step at which the components are set to zero. 
 * > Dummy argument in release mode.
 *
 */
template <typename T>
void herm_matrix_timestep_view<T>::set_timestep_zero(int tstp) {
    assert(tstp == tstp_);
    assert(tstp >= -1);
    if (tstp == -1) {
        memset(matptr(0), 0, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memset(retptr(0), 0, sizeof(cplx) * (tstp + 1) * element_size_);
        memset(tvptr(0), 0, sizeof(cplx) * (ntau_ + 1) * element_size_);
        memset(lesptr(0), 0, sizeof(cplx) * (tstp + 1) * element_size_);
    }
}


/** \brief <b> Sets all components of `herm_matrix_timestep_view`  to the components of
 *  a given `herm_matrix` at time step `tstp`. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets all components of the `herm_matrix_timestep` at time step `tstp` to
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
void herm_matrix_timestep_view<T>::set_timestep(int tstp, herm_matrix<T> &g1) {
    assert(tstp == tstp_);
    assert(tstp >= -1 && tstp <= g1.nt() && "tstp >= -1 && tstp <= g1.nt()");
    assert(g1.size1() == size1_ && "g1.size1() == size1_");
    assert(g1.ntau() == ntau_ && "g1.ntau() == ntau_");
    if (tstp == -1) {
        memcpy(matptr(0), g1.matptr(0), sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memcpy(retptr(0), g1.retptr(tstp, 0),
               sizeof(cplx) * (tstp + 1) * element_size_);
        memcpy(tvptr(0), g1.tvptr(tstp, 0),
               sizeof(cplx) * (ntau_ + 1) * element_size_);
        memcpy(lesptr(0), g1.lesptr(0, tstp),
               sizeof(cplx) * (tstp + 1) * element_size_);
    }
}

/** \brief <b> Sets all components of `herm_matrix_timestep_view`  to the components of
 *  a given `herm_matrix_timestep` at time step `tstp`. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets all components of the `herm_matrix_timestep` at time step `tstp` to
 * > the components of given `herm_matrix_timestep`. If `tstp = -1`, only the
 * > Matsubara component will be copied.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > The time step at which the components are set.
 *
 * @param g1
 * > The `herm_matrix_timestep` from which the time step is copied.
 *
 */
template <typename T>
void herm_matrix_timestep_view<T>::set_timestep(int tstp, herm_matrix_timestep<T> &g1) {
    assert(tstp == tstp_);
    assert(tstp == g1.tstp());
    assert(tstp >= -1 && tstp <= g1.nt() && "tstp >= -1 && tstp <= g1.nt()");
    assert(g1.size1() == size1_ && "g1.size1() == size1_");
    assert(g1.ntau() == ntau_ && "g1.ntau() == ntau_");
    if (tstp == -1) {
        memcpy(matptr(0), g1.matptr(0), sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memcpy(retptr(0), g1.retptr(0),
               sizeof(cplx) * (tstp + 1) * element_size_);
        memcpy(tvptr(0), g1.tvptr(0),
               sizeof(cplx) * (ntau_ + 1) * element_size_);
        memcpy(lesptr(0), g1.lesptr(0),
               sizeof(cplx) * (tstp + 1) * element_size_);
    }
}



/** \brief <b> Left-multiplication of the `herm_matrix_timestep_view` with a contour function </b>
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
void herm_matrix_timestep_view<T>::left_multiply(int tstp, function<T> &ft, T weight) {
   assert(tstp == tstp_);
   assert( ft.nt() >= tstp_);

   this->left_multiply(ft.ptr(-1), ft.ptr(0), weight);
}

/// @private
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
void herm_matrix_timestep_view<T>::left_multiply(int tstp, std::complex<T> *f0,
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


/** \brief <b> Right-multiplication of the `herm_matrix_timestep_view` with a contour function </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Right-multiplication of the `herm_matrix_timestep_view` with a time dependent contour function F(t)
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
void herm_matrix_timestep_view<T>::right_multiply(int tstp, function<T> &ft, T weight) {
   assert(tstp == tstp_);
   assert( ft.nt() >= tstp_);

   this->right_multiply(ft.ptr(-1), ft.ptr(0), weight);
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
void herm_matrix_timestep_view<T>::right_multiply(int tstp, std::complex<T> *f0,
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






/** \brief <b> Left-multiplies the `herm_matrix_timestep_view` with the hermitian conjugate of a contour function. </b>
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
void herm_matrix_timestep_view<T>::left_multiply_hermconj(int tstp, function<T> &ft,
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






/** \brief <b> Reset the pointers of `herm_matrix_timestep_view` class to the pointer given by `data`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and reset the `herm_matrix_timestep_view` class, where data pointer are set to *data
* > If \f$ tstp==-1 \f$ the pointer of data is set to the first element of `mat` and for \f$ tstp>-1 \f$ to the first element of `ret`.
* > Number of time steps is set to `nt`, number of points on the imaginary branch to `ntau`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig` for time slice `tstp`. No data, only pointers are copied.
* > Works for scalar or square-matrix contour objects
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param *data
* > Pointer to the data where you want to reset herm_matrix_timestep_view
* @param tstp
* > Index of time
* @param ntau
* > Number of points on Matsubara axis
* @param size1
* > Matrix rank of the contour function
* @param sig
* > Set `sig = -1` for fermions or `sig = +1` for bosons.
*/

template <typename T>
void herm_matrix_timestep_view<T>::set_to_data(CPLX *data, int tstp, int ntau,
                                               int size, int sig) {
    tstp_ = tstp;
    ntau_ = ntau;
    size1_ = size;
    size2_ = size;
    element_size_ = size1_ * size2_;
    sig_ = sig;
    if (tstp == -1) {
        mat_ = data;
        ret_ = 0;
        les_ = 0;
        tv_ = 0;
    } else if (tstp >= 0) {
        mat_ = 0;
        ret_ = data;
        tv_ = data + (tstp_ + 1) * element_size_;
        les_ = data + (tstp_ + 1 + ntau_ + 1) * element_size_;
    }
}

/** \brief <b> Reset the pointers of `herm_matrix_timestep_view` class to the pointer given by `herm_matrix g` at timestep `tstp`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and reset the `herm_matrix_timestep_view` class, where data pointer are set to `g` at timestep `tstp`.
* > If \f$ tstp==-1 \f$ the pointer of data is set to the first element of `mat` and for \f$ tstp>-1 \f$ to the first element of `ret`.
* > Number of time steps is set to `nt`, number of points on the imaginary branch to `ntau`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig` for time slice `tstp`. No data, only pointers are copied.
* > Works for scalar or square-matrix contour objects
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Index of time
* @param g
* > The `herm_matrix` according to where the pointer of the class members should be reset
*/

template <typename T>
void herm_matrix_timestep_view<T>::set_to_data(int tstp, herm_matrix<T> &g) {
    assert(tstp>=-1 && tstp<=g.nt());
    tstp_ = tstp;
    ntau_ = g.ntau();
    size1_ = g.size1();
    size2_ = g.size2();
    element_size_ = size1_ * size2_;
    sig_ = g.sig();
    if (tstp == -1) {
        mat_ = g.matptr(0);
        ret_ = 0;
        les_ = 0;
        tv_ = 0;
    } else if (tstp >= 0) {
        mat_ = 0;
        ret_ = g.retptr(tstp, 0);
        tv_ = g.tvptr(tstp, 0);
        les_ = g.lesptr(0, tstp);
    }
}

/** \brief <b> Reset the pointers of `herm_matrix_timestep_view` class to the pointer given by `herm_matrix_timestep g` at timestep `tstp`</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and reset the `herm_matrix_timestep_view` class, where data pointer are set to `g` at timestep `tstp`.
* > If \f$ tstp==-1 \f$ the pointer of data is set to the first element of `mat` and for \f$ tstp>-1 \f$ to the first element of `ret`.
* > Number of time steps is set to `nt`, number of points on the imaginary branch to `ntau`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig` for time slice `tstp`. No data, only pointers are copied.
* > Works for scalar or square-matrix contour objects
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix_timestep` according to where the pointer of the class members should be reset
*/

template <typename T>
void herm_matrix_timestep_view<T>::set_to_data(herm_matrix_timestep<T> &g) {
    tstp_ = g.tstp_;
    ntau_ = g.ntau();
    size1_ = g.size1();
    size2_ = g.size2();
    element_size_ = size1_ * size2_;
    sig_ = g.sig();
    if (tstp_ == -1) {
        mat_ = g.matptr(0);
        ret_ = 0;
        les_ = 0;
        tv_ = 0;
    } else if (tstp_ >= 0) {
        mat_ = 0;
        ret_ = g.retptr(0);
        tv_ = g.tvptr(0);
        les_ = g.lesptr(0);
    }
}

/** \brief <b> Return the data</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return the data for each component in `herm_matrix_timestep_view` into the reserved memory
* > Works for scalar or square-matrix contour objects
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param *ret
* > Pointer to memory where ‘ret‘ component is copied
* @param *les
* > Pointer to memory where ‘les‘ component is copied
* @param *tv
* > Pointer to memory where ‘tv‘ component is copied
* @param *mat
* > Pointer to memory where ‘mat‘ component is copied
*/
template <typename T>
void herm_matrix_timestep_view<T>::get_data(CPLX *ret, CPLX *les, CPLX *tv,
                                            CPLX *mat) {
    if (tstp_ == -1) {
        memcpy(mat_, mat, sizeof(CPLX) * (ntau_ + 1) * element_size_);
        ret_ = 0;
        les_ = 0;
        tv_ = 0;
    } else if (tstp_ >= 0) {
        mat_ = 0;
        memcpy(ret_, ret, sizeof(CPLX) * (tstp_ + 1) * element_size_);
        memcpy(les_, les, sizeof(CPLX) * (tstp_ + 1) * element_size_);
        memcpy(tv_, tv, sizeof(CPLX) * (ntau_ + 1) * element_size_);
    }
}

/** \brief <b> Return the data into `herm_matrix_timestep_view`</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return the data for each component in `herm_matrix_timestep_view` into new `herm_matrix_timestep_view`
* > Works for scalar or square-matrix contour objects. I believe this function is obsolete.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix_timestep_view` according to which the class should be initialized
*/

template <typename T>
void herm_matrix_timestep_view<T>::get_data(herm_matrix_timestep_view<T> &g) {
    assert(tstp_==g.tstp_ && size1_ == g.size1() && size2_ == g.size2() && ntau_ == g.ntau());
    if (tstp_ == -1)
        get_data(NULL, NULL, NULL, g.mat_);
    else
        get_data(g.ret_, g.les_, g.tv_, NULL);
}


/** \brief <b> Return the data into object given by the `template argument`</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return the data for each component in `template argument` [ either herm_matrix_timestep_view, herm_matrix_timestep,herm_matrix] into new `herm_matrix_timestep_view`
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `template argument` according to which the data should be set
*/

template <typename T>
template <class GG>
void herm_matrix_timestep_view<T>::get_data(GG &g) {
    herm_matrix_timestep_view<T> tmp(tstp_, g);
    get_data(tmp);
}

/** \brief <b> Get timestep into herm_matrix_timestep </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Get timestep into herm_matrix_timestep
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > time step
* @param timestep
* > The `template argument` according to which the data should be set
*/

template <typename T>
void herm_matrix_timestep_view<T>::get_timestep(int tstp,
                                  herm_matrix_timestep<T> &timestep) {
    int len = (2 * (tstp + 1) + ntau_ + 1) * element_size_;
    cplx *x;
    assert(tstp_ ==  tstp);
    if (timestep.total_size_ < len)
        timestep.resize(tstp, ntau_, size1_);
    x = timestep.data_;
    timestep.tstp_ = tstp;
    timestep.ntau_ = ntau_;
    timestep.size1_ = size1_;
    if (tstp == -1) {
        memcpy(x, mat_, sizeof(cplx) * (ntau_ + 1) * element_size_);
    } else {
        memcpy(x, ret_, sizeof(cplx) * (tstp + 1) * element_size_);
        memcpy(x + (tstp + 1) * element_size_, tv_,
               sizeof(cplx) * (ntau_ + 1) * element_size_);
        memcpy(x + (tstp + 1 + ntau_ + 1) * element_size_, les_,
               sizeof(cplx) * (tstp + 1) * element_size_);
    }
}


/** \brief <b> Set the matrix element of \f$C_{[i1,i2]}\f$ for each component from  \f$ g_{[i1,i2]} \f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set the matrix element \f$C_{[i1,i2]}\f$ for each component from  \f$ g_{[i1,i2]} \f$.
* > If \f$t>-1\f$ then `ret,les,tv` components are set, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i1
* > Row index of `herm_matrix_timestep_view`
* @param i2
* > Column index of `herm_matrix_timestep_view`
* @param g
* > The `herm_matrix_timestep_view` from which the matrix element element is given
* @param j1
* > Row index of `herm_matrix_timestep_view`
* @param j2
* > Column index of `herm_matrix_timestep_view`
*/

template <typename T>
void herm_matrix_timestep_view<T>::set_matrixelement(
    int i1, int i2, herm_matrix_timestep_view<T> &g, int j1, int j2) {
    int i, sij = i1 * size2_ + i2, tij = j1 * g.size2() + j2;
    assert(tstp_==g.tstp_ && ntau_ == g.ntau_);

    assert(i1>=0 && i1 <=size1_ -1);
    assert(i2>=0 && i2 <=size2_ -1);
    assert(j1>=0 && j1 <=g.size1_-1);
    assert(j2>=0 && j2 <=g.size2_-1);

    if (tstp_ == -1) {
        for (i = 0; i <= ntau_; i++)
            matptr(i)[sij] = g.matptr(i)[tij];
    } else {
        for (i = 0; i <= tstp_; i++)
            retptr(i)[sij] = g.retptr(i)[tij];
        for (i = 0; i <= ntau_; i++)
            tvptr(i)[sij] = g.tvptr(i)[tij];
        for (i = 0; i <= tstp_; i++)
            lesptr(i)[sij] = g.lesptr(i)[tij];
    }
}


/** \brief <b> Set the matrix element of \f$C_{[i1,i2]}\f$ for each component from \f$ g_{[i1,i2]} \f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set the matrix element \f$C_{[i1,i2]}\f$ for each component from  \f$ g_{[i1,i2]} \f$.
* > If \f$t>-1\f$ then `ret,les,tv` components are set, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i1
* > Row index of `herm_matrix_timestep_view`
* @param i2
* > Column index of `herm_matrix_timestep_view`
* @param g
* > The `template argument` from which the matrix element is given
* @param j1
* > Row index of `herm_matrix_timestep_view`
* @param j2
* > Column index of `herm_matrix_timestep_view`
*/

template <typename T>
template <class GG>
void herm_matrix_timestep_view<T>::set_matrixelement(int i1, int i2, GG &g,
                                                     int j1, int j2) {
    herm_matrix_timestep_view<T> tmp(tstp_, g);
    set_matrixelement(i1, i2, tmp, j1, j2);
}


#define HERM_MATRIX_INCR_TSTP                                                \
    if (alpha == CPLX(1.0, 0.0)) {                                           \
        for (i = 0; i < len; i++)                                            \
            x0[i] += x[i];                                                   \
    } else {                                                                 \
        for (i = 0; i < len; i++)                                            \
            x0[i] += alpha * x[i];                                           \
    }


/** \brief <b> Increase the value of the  `herm_matrix_timestep_view` for  \f$\alpha g(t) \f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Increase the value of `herm_matrix_timestep_view` for a value of \f$\alpha g(t)\f$, where \f$ g(t)\f$
* > is a `herm_matrix_timestep_view` and \f$\alpha\f$ is a complex number
* > If \f$t>-1\f$ then `ret,les,tv` components are set, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix_timestep_view` which is added
* @param alpha
* > Constant multiplication factor
*/

template <typename T>
void
herm_matrix_timestep_view<T>::incr_timestep(herm_matrix_timestep_view<T> &g,
                                            CPLX alpha) {
    assert(tstp_==g.tstp_ && ntau_ ==g.ntau_ && size1_ == g.size1_ && size2_ == g.size2_);

    int i, len;
    CPLX *x, *x0;
    if (tstp_ == -1) {
        len = (ntau_ + 1) * element_size_;
        x0 = mat_;
        x = g.mat_;
        HERM_MATRIX_INCR_TSTP
    } else {
        len = (tstp_ + 1) * element_size_;
        x0 = ret_;
        x = g.ret_;
        HERM_MATRIX_INCR_TSTP
        x0 = les_;
        x = g.les_;
        HERM_MATRIX_INCR_TSTP
        len = (ntau_ + 1) * element_size_;
        x0 = tv_;
        x = g.tv_;
        HERM_MATRIX_INCR_TSTP
    }
}
#undef HERM_MATRIX_INCR_TSTP

/** \brief <b> Increase the value of the  `herm_matrix_timestep_view` for  \f$\alpha g(t) \f$.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Increase the value of `herm_matrix_timestep_view` for a value of \f$\alpha g(t) \f$, where \f$ g(t)\f$
* > is a `template argument` and \f$\alpha\f$ is a complex number
* > If \f$t>-1\f$ then `ret,les,tv` components are set, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `template argument` which is added
* @param alpha
* > Constant multiplication factor
*/


template <typename T>
template <class GG>
void herm_matrix_timestep_view<T>::incr_timestep(GG &g, CPLX alpha) {
    herm_matrix_timestep_view<T> tmp(tstp_, g);
    incr_timestep(tmp, alpha);
}
#define HERM_MATRIX_INCR_TSTP                                                \
    if (alpha == 1.0) {                                                      \
        for (i = 0; i < len; i++)                                            \
            x0[i] += x[i];                                                   \
    } else {                                                                 \
        for (i = 0; i < len; i++)                                            \
            x0[i] += alpha * x[i];                                           \
    }


/** \brief <b> Increase the value of the  `herm_matrix_timestep_view` for  \f$\alpha g(t) \f$.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Increase the value of `herm_matrix_timestep_view` for a value of \f$\alpha g(t)\f$, where \f$ g(t)\f$
* > is a `herm_matrix_timestep_view` and \f$\alpha \f$ is a `template argument`
* > If \f$t>-1\f$ then `ret,les,tv` components are set, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `template argument` which is added
* @param alpha
* > Constant `template argument` multiplication factor
*/

template <typename T>
void
herm_matrix_timestep_view<T>::incr_timestep(herm_matrix_timestep_view<T> &g,
                                            T alpha) {
    assert(tstp_==g.tstp_ && ntau_ == g.ntau_ && size1_ == g.size1_ && size2_ == g.size2_);

    int i, len;
    CPLX *x, *x0;
    if (tstp_ == -1) {
        len = (ntau_ + 1) * element_size_;
        x0 = mat_;
        x = g.mat_;
        HERM_MATRIX_INCR_TSTP
    } else {
        len = (tstp_ + 1) * element_size_;
        x0 = ret_;
        x = g.ret_;
        HERM_MATRIX_INCR_TSTP
        x0 = les_;
        x = g.les_;
        HERM_MATRIX_INCR_TSTP
        len = (ntau_ + 1) * element_size_;
        x0 = tv_;
        x = g.tv_;
        HERM_MATRIX_INCR_TSTP
    }
}
#undef HERM_MATRIX_INCR_TSTP

/** \brief <b> Increase the value of the  `template argument for  \f$\alpha g(t) \f$.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Increase the value of `template argument` for a value of \f$\alpha g(t)\f$, where \f$ g(t)\f$
* > is a `template argument`[ either herm_matrix_timestep_view, herm_matrix_timestep,herm_matrix]
* > and \f$\alpha\f$ is a `template argument`. If \f$t>-1\f$ then `ret,les,tv` components are set, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `template argument` added
* @param alpha
* > Constant `template argument` multiplication factor
*/

template <typename T>
template <class GG>
void herm_matrix_timestep_view<T>::incr_timestep(GG &g, T alpha) {
    herm_matrix_timestep_view<T> tmp(tstp_, g);
    incr_timestep(tmp, alpha);
}

/** \brief <b> Multiply  `herm_matrix_timestep_view` with scalar `weight`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Multiply `herm_matrix_timestep_view` with a scalar.
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
void herm_matrix_timestep_view<T>::smul(T weight) {
    int m;
    CPLX *x0;
    if (tstp_ == -1) {
        x0 = mat_;
        for (m = 0; m <= ntau_; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
                                       weight);
        }
    } else {
        x0 = ret_;
        for (m = 0; m <= tstp_; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
                                       weight);
        }
        x0 = tv_;
        for (m = 0; m <= ntau_; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
                                       weight);
        }
        x0 = les_;
        for (m = 0; m <= tstp_; m++) {
            element_smul<T, LARGESIZE>(size1_, x0 + m * element_size_,
                                       weight);
        }
    }
}
/////////////////

#if CNTR_USE_MPI == 1

/// @private
template <typename T>
void my_mpi_reduce(std::complex<T> *data, int len, int root) {
    std::cerr << __PRETTY_FUNCTION__ << ", LEN=" << len
               << " ... NOT DEFINED FOR THIS TYPE " << std::endl;
     exit(0);
 }


/// @private
template<>
inline void my_mpi_reduce<double>(std::complex<double> *data, int len, int root) {
    int tid = MPI::COMM_WORLD.Get_rank();
    int ntasks = MPI::COMM_WORLD.Get_size();
    assert(root>=0 && root <= ntasks -1);
    assert(len>=0);
    if (tid == root) {
        MPI::COMM_WORLD.Reduce(MPI::IN_PLACE, (double *)data, 2 * len,
                               MPI::DOUBLE, MPI::SUM, root);
    } else {
        MPI::COMM_WORLD.Reduce((double *)data, (double *)data, 2 * len,
                               MPI::DOUBLE, MPI::SUM, root);
    }
}

/// @private
/** \brief <b> MPI reduce for the `herm_matrix_timestep_view`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > MPI reduce for the `herm_matrix_timestep_view` to the `root`
* > If \f$t>-1\f$ then `ret,les,tv` components are set, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param root
* > Index of root
*/

template <typename T>
void herm_matrix_timestep_view<T>::MPI_Reduce(int root) {
    if (tstp_ == -1) {
        my_mpi_reduce<T>(mat_, (ntau_ + 1) * element_size_, root);
    } else {
        my_mpi_reduce<T>(les_, (tstp_ + 1) * element_size_, root);
        my_mpi_reduce<T>(ret_, (tstp_ + 1) * element_size_, root);
        my_mpi_reduce<T>(tv_, (ntau_ + 1) * element_size_, root);
    }
}

/** \brief <b> Broadcasts the `herm_matrix_timestep_view` at a given time step to all ranks. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Broadcasts the `herm_matrix_timestep_view` at a given time step `tstp` to all ranks.
* > Works for a square matrices.
*
*<!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Time step which should be broadcasted.
* @param root
* > The rank from which the `herm_matrix_timestep_view` should be broadcasted.
*
*/
template <typename T>
void herm_matrix_timestep_view<T>::Bcast_timestep(int tstp, int root){
   int numtasks = MPI::COMM_WORLD.Get_size();
   int taskid = MPI::COMM_WORLD.Get_rank();
   // test effective on root:
   assert(tstp == tstp_);
   if (sizeof(T) == sizeof(double)){
      if (tstp_ == -1){
         MPI::COMM_WORLD.Bcast(mat_, (ntau_ + 1) * element_size_, MPI::DOUBLE_COMPLEX, root);
      } else {
         MPI::COMM_WORLD.Bcast(les_, (tstp_ + 1) * element_size_, MPI::DOUBLE_COMPLEX, root);
         MPI::COMM_WORLD.Bcast(ret_, (tstp_ + 1) * element_size_, MPI::DOUBLE_COMPLEX, root);
         MPI::COMM_WORLD.Bcast(tv_, (ntau_ + 1) * element_size_, MPI::DOUBLE_COMPLEX, root);
      }
   } else { // assuming single precision
      if (tstp_ == -1){
         MPI::COMM_WORLD.Bcast(mat_, (ntau_ + 1) * element_size_, MPI::COMPLEX, root);
      } else {
         MPI::COMM_WORLD.Bcast(les_, (tstp_ + 1) * element_size_, MPI::COMPLEX, root);
         MPI::COMM_WORLD.Bcast(ret_, (tstp_ + 1) * element_size_, MPI::COMPLEX, root);
         MPI::COMM_WORLD.Bcast(tv_, (ntau_ + 1) * element_size_, MPI::COMPLEX, root);
      }
   }

}

/** \brief <b> Sends the `herm_matrix_timestep_view` at a given time step to a specific task. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Sends the `herm_matrix_timestep_view` at a given time step `tstp` to a specific task with rank `dest`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Time step which should be send.
* @param dest
* > The task rank to which the `herm_matrix_timestep_view` should be send.
* @param tag
* > The MPI error flag.
*/
template <typename T>
void herm_matrix_timestep_view<T>::Send_timestep(int tstp, int dest, int tag) {
   int taskid = MPI::COMM_WORLD.Get_rank();
   assert(tstp == tstp_);
   if (!(taskid == dest)) {
      if (sizeof(T) == sizeof(double)){
         if (tstp_ == -1){
            MPI::COMM_WORLD.Send(mat_, (ntau_ + 1) * element_size_, MPI::DOUBLE_COMPLEX, dest, tag);
         } else {
            MPI::COMM_WORLD.Send(les_, (tstp_ + 1) * element_size_, MPI::DOUBLE_COMPLEX, dest, tag);
            MPI::COMM_WORLD.Send(ret_, (tstp_ + 1) * element_size_, MPI::DOUBLE_COMPLEX, dest, tag);
            MPI::COMM_WORLD.Send(tv_, (ntau_ + 1) * element_size_, MPI::DOUBLE_COMPLEX, dest, tag);
         }
      }
      else {
         if (tstp_ == -1){
            MPI::COMM_WORLD.Send(mat_, (ntau_ + 1) * element_size_, MPI::COMPLEX, dest, tag);
         } else {
            MPI::COMM_WORLD.Send(les_, (tstp_ + 1) * element_size_, MPI::COMPLEX, dest, tag);
            MPI::COMM_WORLD.Send(ret_, (tstp_ + 1) * element_size_, MPI::COMPLEX, dest, tag);
            MPI::COMM_WORLD.Send(tv_, (ntau_ + 1) * element_size_, MPI::COMPLEX, dest, tag);
         }
      }
   } else {
      // do nothing
   }
}

/** \brief <b> Recevies the `herm_matrix_timestep_view` at a given time step from a specific task. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Receives the `herm_matrix_timestep_view` at a given time step `tstp` from a specific task with rank `root`.
*
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
void herm_matrix_timestep_view<T>::Recv_timestep(int tstp, int root, int tag) {
   int taskid = MPI::COMM_WORLD.Get_rank();
   assert(tstp == tstp_);
   if (!(taskid == root)) {
      if (sizeof(T) == sizeof(double))
         if (tstp_ == -1){
            MPI::COMM_WORLD.Recv(mat_, (ntau_ + 1) * element_size_, MPI::DOUBLE_COMPLEX, root, tag);
         } else {
            MPI::COMM_WORLD.Recv(les_, (tstp_ + 1) * element_size_, MPI::DOUBLE_COMPLEX, root, tag);
            MPI::COMM_WORLD.Recv(ret_, (tstp_ + 1) * element_size_, MPI::DOUBLE_COMPLEX, root, tag);
            MPI::COMM_WORLD.Recv(tv_, (ntau_ + 1) * element_size_, MPI::DOUBLE_COMPLEX, root, tag);
         }
      else{
         if (tstp_ == -1){
            MPI::COMM_WORLD.Recv(mat_, (ntau_ + 1) * element_size_, MPI::COMPLEX, root, tag);
         } else {
            MPI::COMM_WORLD.Recv(les_, (tstp_ + 1) * element_size_, MPI::COMPLEX, root, tag);
            MPI::COMM_WORLD.Recv(ret_, (tstp_ + 1) * element_size_, MPI::COMPLEX, root, tag);
            MPI::COMM_WORLD.Recv(tv_, (ntau_ + 1) * element_size_, MPI::COMPLEX, root, tag);
         }
      }
   }
}


#endif

/** \brief <b> Write `herm_matrix_timestep_view` to hdf5 group given by `group_id`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Write `herm_matrix_timestep_view` to the hdf5 group given by `group_id`
* > If \f$t>-1\f$ then `ret,les,tv` components are written, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param group_id
* > Group ID of the hdf5 file to write in
*/

#if CNTR_USE_HDF5 == 1
template <typename T>
void herm_matrix_timestep_view<T>::write_to_hdf5(hid_t group_id) {
    store_int_attribute_to_hid(group_id, std::string("ntau"), ntau_);
    store_int_attribute_to_hid(group_id, std::string("tstp"), tstp_);
    store_int_attribute_to_hid(group_id, std::string("sig"), sig_);
    store_int_attribute_to_hid(group_id, std::string("size1"), size1_);
    store_int_attribute_to_hid(group_id, std::string("size2"), size2_);
    store_int_attribute_to_hid(group_id, std::string("element_size"),
                               element_size_);
    hsize_t len_shape = 3, shape[3];
    shape[1] = size1_;
    shape[2] = size2_;
    if (tstp_ == -1) {
        shape[0] = ntau_ + 1;
        store_cplx_array_to_hid(group_id, std::string("mat"), matptr(0),
                                shape, len_shape);
    } else if (tstp_ > -1) {
        shape[0] = (tstp_ + 1);
        // CHECK: implement store_cplx_array_to_hid with template typename T
        store_cplx_array_to_hid(group_id, std::string("ret"), retptr(0),
                                shape, len_shape);
        store_cplx_array_to_hid(group_id, std::string("les"), lesptr(0),
                                shape, len_shape);
        shape[0] = (ntau_ + 1);
        store_cplx_array_to_hid(group_id, std::string("tv"), tvptr(0), shape,
                                len_shape);
    }
}

/** \brief <b> Write `herm_matrix_timestep_view` to hdf5 group given by `group_id` and name `groupname`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Write `herm_matrix_timestep_view` to the hdf5 group given by `group_id` and `groupname`
* > If \f$t>-1\f$ then `ret,les,tv` components are written, otherwise `mat`.
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
void herm_matrix_timestep_view<T>::write_to_hdf5(hid_t group_id,
                                                 const char *groupname) {
    hid_t sub_group_id = create_group(group_id, groupname);
    this->write_to_hdf5(sub_group_id);
    close_group(sub_group_id);
}

/** \brief <b> Write `herm_matrix_timestep_view` to hdf5 group `groupname` in file `filename`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Write `herm_matrix_timestep_view` to the hdf5 group `groupname` in file  `filename`
* > If \f$t>-1\f$ then `ret,les,tv` components are written, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param filename
* > hdf5 filename for file to write in
* @param groupname
* > hdf5 groupname for file to write in
*/

template <typename T>
void herm_matrix_timestep_view<T>::write_to_hdf5(const char *filename,
                                                 const char *groupname) {
    hid_t file_id = open_hdf5_file(filename);
    this->write_to_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
}

/** \brief <b> Read `herm_matrix_timestep_view` from hdf5 group given by group ID `group_id`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read `herm_matrix_timestep_view` from the hdf5 group given by group ID `group_id`
* > If \f$t>-1\f$ then `ret,les,tv` components are written, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param group_id
* > group ID for group to read from
*/

template <typename T>
void herm_matrix_timestep_view<T>::read_from_hdf5(hid_t group_id) {
    // -- Read data: NO RESIZE POSSIBLE
    int tstp = read_primitive_type<int>(group_id, "tstp");
    int ntau = read_primitive_type<int>(group_id, "ntau");
    int sig = read_primitive_type<int>(group_id, "sig");
    int size1 = read_primitive_type<int>(group_id, "size1");
    int size2 = read_primitive_type<int>(group_id, "size2");
    int element_size = read_primitive_type<int>(group_id, "element_size");

    assert(element_size == element_size_ && size1 == size1_ && size2 == size2_ && tstp == tstp_ && ntau ==ntau_);

    sig_ = sig;
    if (tstp == -1) {
        hsize_t mat_size = (ntau + 1) * element_size;
        read_primitive_type_array(group_id, "mat", mat_size, matptr(0));
    } else if (tstp_ > -1) {
        hsize_t ret_size = (tstp_ + 1) * element_size;
        read_primitive_type_array(group_id, "ret", ret_size, retptr(0));
        hsize_t les_size = (tstp_ + 1) * element_size;
        read_primitive_type_array(group_id, "les", les_size, lesptr(0));
        hsize_t tv_size = (ntau_ + 1) * element_size;
        read_primitive_type_array(group_id, "tv", tv_size, tvptr(0));
    }
}

/** \brief <b> Read `herm_matrix_timestep_view` from hdf5 group given by group ID `group_id` and group name `groupname`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read `herm_matrix_timestep_view` from the hdf5 group given by group ID `group_id` and group name `groupname`
* > If \f$t>-1\f$ then `ret,les,tv` components are written, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param group_id
* > group ID `group_id` for group to read from
* @param groupname
* > group name `groupname` for group to read from
*/

template <typename T>
void herm_matrix_timestep_view<T>::read_from_hdf5(hid_t group_id,
                                                  const char *groupname) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(sub_group_id);
    close_group(sub_group_id);
}

/** \brief <b> Read `herm_matrix_timestep_view` from hdf5 group `groupname` in file `filename`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read `herm_matrix_timestep_view` from the hdf5 group `groupname` in file  `filename`
* > If \f$t>-1\f$ then `ret,les,tv` components are written, otherwise `mat`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param filename
* > hdf5 filename for file to read from
* @param groupname
* > hdf5 groupname for file to read from
*/
template <typename T>
void herm_matrix_timestep_view<T>::read_from_hdf5(const char *filename,
                                                  const char *groupname) {
    hid_t file_id = read_hdf5_file(filename);
    this->read_from_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
}
#endif

/// @private
/** \brief <b> Return single particle density matrix from `herm_matrix_timestep_view`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return single particle density matrix (occupation) from `herm_matrix_timestep_view`
* > Works for scalar.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > timestep index
*/
template <typename T>
CPLX herm_matrix_timestep_view<T>::density_matrix(int tstp) {
    assert(tstp==tstp_);

    CPLX x1;
    if (tstp_ == -1) {
        x1 = *matptr(ntau_);
        return -x1;
    } else {
        x1 = *lesptr(tstp_);
        return CPLX(0.0, sig_) * x1;
    }
}


/** \brief <b> Return single particle density matrix from `herm_matrix_timestep_view`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return single particle density matrix (occupation) from `herm_matrix_timestep_view`
* > Works for scalar.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > timestep index
* @param M
* > density matrix, returned as complex number
*/
template <typename T>
void herm_matrix_timestep_view<T>::density_matrix(int tstp, std::complex<T> &M) {
    assert(tstp==tstp_);
    if (tstp_ == -1) {
        M = -(*matptr(ntau_));
    } else {
        M = CPLX(0.0, sig_)*(*lesptr(tstp_));
    }
}

/** \brief <b> Return single particle density matrix as Eigen matrix from `herm_matrix_timestep_view`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return single particle density matrix as Eigen matrix from `herm_matrix_timestep_view`
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > timestep index
* @param M [Matrix]
* > Eigen matrix representation of the density matrix
*/
template <typename T>
template <class Matrix>
void herm_matrix_timestep_view<T>::density_matrix(int tstp, Matrix &M) {
    assert(tstp == tstp_);

    CPLX *x;
    if (tstp == -1) {
        x = matptr(ntau_);
        herm_matrix_READ_ELEMENT M *= (-1.0);
    } else {
        x = lesptr(tstp);
        herm_matrix_READ_ELEMENT M *= CPLX(0.0, 1.0 * sig_);
    }
}

#undef herm_matrix_READ_ELEMENT
#undef herm_matrix_READ_ELEMENT_MINUS_CONJ
#undef CPLX

} // namespace cntr

#endif  // CNTR_HERM_TIMESTEP_VIEW_IMPL_H
