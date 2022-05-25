#ifndef CNTR_HERM_TIMESTEP_MOVING_VIEW_IMPL_H
#define CNTR_HERM_TIMESTEP_MOVING_VIEW_IMPL_H

#include "cntr_herm_matrix_timestep_moving_view_decl.hpp"
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
  herm_matrix_timestep_moving_view<T>::herm_matrix_timestep_moving_view() {
    ret_=0;
    les_=0;
    tc_=-1;
    t0_=0;
    size1_=0;
    size2_=0;
    element_size_=0;
    sig_= -1;
  }
  
  template <typename T>
  herm_matrix_timestep_moving_view<T>::~herm_matrix_timestep_moving_view() {
    // donothing
  }
/** \brief <b> Initializes the `herm_matrix_timestep_moving_view` class with the same layout as a given `herm_matrix_timestep_moving_view g`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and copy the `herm_matrix_timestep_moving_view` class with the same 
* > cutoff time steps `tc`, physical time `t0`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig`.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix_timestep_moving_view` according to which the class should be initialized
*/ 
  template <typename T>
  herm_matrix_timestep_moving_view<T>::herm_matrix_timestep_moving_view(const herm_matrix_timestep_moving_view &g) {
    tc_=g.tc_;
    t0_=g.t0_;
    size1_=g.size1_;
    size2_=g.size2_;
    element_size_=g.element_size_;
    sig_= g.sig_;
    // copy pointers upgrade to shared pointers??
    ret_=g.ret_;
    les_=g.les_;
  }
  
  template <typename T>
  herm_matrix_timestep_moving_view<T> &herm_matrix_timestep_moving_view<T>::operator=(const herm_matrix_timestep_moving_view &g){
    if(this == &g) return *this;
    tc_=g.tc_;
    t0_=g.t0_;
    size1_=g.size1_;
    size2_=g.size2_;
    element_size_=size1_*size2_;
    sig_= g.sig_;
    // copy pointers upgrade to shared pointers??
    ret_=g.ret_;
    les_=g.les_;
    return *this;
  }
/** \brief <b> Initializes the `herm_matrix_timestep_moving_view` class with the same layout as a given `herm_matrix_timestep_moving g`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and copy the `herm_matrix_timestep_moving_view` class with the same 
* > cutoff time steps `tc`, physical time `t0`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig`.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix_timestep_moving` according to which the class should be initialized
*/   
  template<typename T>
  herm_matrix_timestep_moving_view<T>::herm_matrix_timestep_moving_view(herm_matrix_timestep_moving<T> &g){
    tc_=g.tc();
    t0_=g.t0();
    size1_=g.size1();
    size2_=g.size2();
    element_size_=size1_*size2_;
    sig_= g.sig();
    ret_=g.retptr(0);
    les_=g.lesptr(0);
  }
/** \brief <b> Initializes the `herm_matrix_timestep_moving_view` class with the same layout as a given `herm_matrix_moving g`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and copy the `herm_matrix_timestep_moving_view` class with the same 
* > cutoff time steps `tc`, physical time `t0`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig`. References to the data stored in physical timestep 'tstp'.
*
*
* <!-- ARGUMENTS
*      ========= -->
** @param tstp
* > Time step 
* @param g
* > The `herm_matrix_moving` according to which the class should be initialized
*/ 
  template <typename T>
  herm_matrix_timestep_moving_view<T>::herm_matrix_timestep_moving_view(int i, herm_matrix_moving<T> &g){
    assert(tstp<=g.t0() && tstp>=(g.t0()-g.tc()) && g.t0()>0);
    tc_=g.tc();
    t0_=g.t0();
    size1_=g.size1();
    size2_=g.size2();
    element_size_=size1_*size2_;
    sig_= g.sig();
    //only valid for the wasteful square storage scheme currently implemented
    ret_=g.retptr(t0_-i,0);
    les_=g.lesptr(t0_-i,0);
  }
  /** \brief <b> Initializes the `herm_matrix_timestep_moving_view` class with the same layout as a given `herm_matrix_timestep_moving g`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and copy the `herm_matrix_timestep_moving_view` class with the same 
* > cutoff time steps `tc`, physical time `t0`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig`. References to the data stored in physical timestep 't0'.
*
*
* <!-- ARGUMENTS
*      ========= -->
* @param t0
* > Time step 
* @param g
* > The `herm_matrix_timestep_moving` according to which the class should be initialized
*/
  template <typename T>
  herm_matrix_timestep_moving_view<T>::herm_matrix_timestep_moving_view(int t0, herm_matrix_timestep_moving<T> &g){
    assert(t0==g.t0() && g.t0()>0);
    tc_=g.tc();
    t0_=t0;
    size1_=g.size1();
    size2_=g.size2();
    element_size_=size1_*size2_;
    sig_= g.sig();
    //only valid for the wasteful square storage scheme currently implemented
    ret_=g.retptr(0);
    les_=g.lesptr(0);
  }
/** \brief <b> Initializes the `herm_matrix_timestep_moving_view` class with the same layout as a given `herm_matrix_timestep_moving_view g`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and copy the `herm_matrix_timestep_moving_view` class with the same 
* > cutoff time steps `tc`, physical time `t0`, matrix rank `size1` and
* > bosonic/fermionic symmetry `sig`.
*
*
* <!-- ARGUMENTS
*      ========= -->
* @param t0
* > Time step 
* @param g
* > The `herm_matrix_timestep_moving_view` according to which the class should be initialized
*/
  template <typename T>
  herm_matrix_timestep_moving_view<T>::herm_matrix_timestep_moving_view(int t0, herm_matrix_timestep_moving_view<T> &g){
    assert(t0==g.t0() && g.t0()>0);
    tc_=g.tc_;
    t0_=t0;
    size1_=g.size1_;
    size2_=g.size2_;
    element_size_=size1_*size2_;
    sig_= g.sig_;
    //only valid for the wasteful square storage scheme currently implemented
    ret_=g.ret_;
    les_=g.les_;
  }
/** \brief <b> Initializes the `herm_matrix_timestep_moving_view` class for a general matrix. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `herm_matrix_timestep_moving_view` class for a general matrix,
* > where the number of column and rows can be different.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tc
* > Cutoff time
* @param t0
* > Physical time
* @param size1
* > Number of matrix rows
* @param size2
* > Number of matrix columns
* @param sig
* > Set `sig = -1` for fermions or `sig = +1` for bosons.
*/  
  template <typename T>
  herm_matrix_timestep_moving_view<T>::herm_matrix_timestep_moving_view(int tc, int t0, int size1, int size2, int sig) {
    assert(size1>=0 && size2>=0 && tc>=1 && t0>=tc);
    tc_=tc;
    t0_=t0;
    size1_=size1;
    size2_=size2;
    element_size_=size1_*size2_;
    sig_=sig;
    ret_=0;
    les_=0;
  }
  /** \brief <b> Reset the pointers of `herm_matrix_timestep_moving_view` class to the pointer given by `data`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and reset the `herm_matrix_timestep_view` class, where data pointer are set to *data
* > Number of cutoff time is set to `tc`, physical time is set to 't0', matrix rank `size1` and
* > bosonic/fermionic symmetry `sig`. No data, only pointers are copied.
* > Works for scalar or square-matrix contour objects
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param *data
* > Pointer to the data where you want to reset herm_matrix_timestep_moving_view
* @param tc
* > cutoff time
* @param t0
* > physical time
* @param size1
* > Matrix rank of the contour function
* @param sig
* > Set `sig = -1` for fermions or `sig = +1` for bosons.
*/
  template <typename T>
  void herm_matrix_timestep_moving_view<T>::set_to_data(std::complex<T> *data, int tc, int t0, int size, int sig){
    assert(t0>0 && tc>0 && size>0 && tc>0);
    tc_=tc;
    t0_=t0;
    size1_=size;
    size2_=size;
    element_size_=size1_*size2_;
    sig_=sig;
    long ndata1=(tc_+1)*element_size_;
    ret_=data;//
    les_=data+ndata1;//
  }
  /** \brief <b> Reset the pointers of `herm_matrix_timestep_moving_view` class to the pointer given by `herm_matrix_moving g` at timestep `tstp`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and reset the `herm_matrix_timestep_moving_view` class, where data pointer are set to `g` at timestep `tstp`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > Time step
* @param g
* > The `herm_matrix_moving` according to where the pointer of the class members should be reset
*/
  template <typename T>
  void herm_matrix_timestep_moving_view<T>::set_to_data(int tstp, herm_matrix_moving<T> &g){
    assert(tstp<=g.t0() && tstp>=(g.t0()-g.tc()) && g.t0()>0);
    tc_=g.tc_;
    t0_=tstp;
    size1_=g.size1_;
    size2_=g.size2_;
    element_size_=size1_*size2_;
    sig_= g.sig_;
    ret_=g.retptr(t0_-tstp,0);
    les_=g.lesptr(t0_-tstp,0);
  }
/** \brief <b> Reset the pointers of `herm_matrix_timestep_moving_view` class to the pointer given by `herm_matrix_timestep_moving g`</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes and reset the `herm_matrix_timestep_moving_view` class, where data pointer are set to `g`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix_timestep_moving` according to where the pointer of the class members should be reset
*/
  template <typename T>
  void herm_matrix_timestep_moving_view<T>::set_to_data(herm_matrix_timestep_moving<T> &g){
    tc_=g.tc_;
    t0_=g.t0_;
    size1_=g.size1_;
    size2_=g.size2_;
    element_size_=size1_*size2_;
    sig_= g.sig_;
    ret_=g.retptr(0);
    les_=g.lesptr(0);   
  }
/** \brief <b> Return the data</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return the data for each component in `herm_matrix_timestep_moving_view` into the 
* > reserved memory. Works for scalar or square-matrix contour objects
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param *ret
* > Pointer to memory where ‘ret‘ component is copied
* @param *les
* > Pointer to memory where ‘les‘ component is copied
*/  
  template <typename T>
  void  herm_matrix_timestep_moving_view<T>::get_data(CPLX *ret,CPLX *les){
    memcpy(ret_, ret, sizeof(CPLX) * (tc_ + 1) *element_size_);
    memcpy(les_, les, sizeof(CPLX) * (tc_ + 1) *element_size_); 
  }
/** \brief <b> Return the data into `herm_matrix_timestep_moving_view`</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return the data for each component in `herm_matrix_timestep_moving_view` 
* > into new `herm_matrix_timestep_moving_view`.
* > Works for scalar or square-matrix contour objects.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix_timestep_moving_view` according to which the class should be initialized
*/
  template <typename T>
  void  herm_matrix_timestep_moving_view<T>::get_data(herm_matrix_timestep_moving_view<T> &g){
    assert(t0_==g.t0() && tc_==g.tc() && size1_==g.size1() && size2_==g.size2());
    get_data(g.ret_,g.les_);   
      }
/** \brief <b> Return the data into object given by the `template argument`</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return the data for each component in `template argument` [ either 
* > herm_matrix_timestep_moving_view, herm_matrix_timestep_moving,herm_matrix_moving] 
* > into new `herm_matrix_timestep_moving_view`.
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
  void  herm_matrix_timestep_moving_view<T>::get_data(GG &g){
    herm_matrix_timestep_moving_view<T> tmp(t0_, g);
    get_data(tmp);   
  }
/** \brief <b> Return the data into object given by the `template argument`</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return the data for each component in `template argument` [ either 
* > herm_matrix_timestep_moving_view, herm_matrix_timestep_moving,herm_matrix_moving] 
* > into new `herm_matrix_timestep_moving_view`.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
* @param tstp
* > Physical time passed to the herm_matrix_timestep_moving_view
* @param g
* > The `template argument` according to which the data should be set
*/
  template <typename T>
  template <class GG>
  void  herm_matrix_timestep_moving_view<T>::get_data(int tstp ,GG &g){
    herm_matrix_timestep_moving_view<T> tmp(tstp, g);
    get_data(tmp);   
  }
/** \brief <b> Set the matrix element of \f$C_{[i1,i2]}\f$ for each component from  \f$ g_{[j1,j2]} \f$</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set the matrix element \f$C_{[i1,i2]}\f$ for each component from  \f$ g_{[j1,j2]} \f$.
* > Works for scalar or square-matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i1
* > Row index of `herm_matrix_timestep_moving_view`
* @param i2
* > Column index of `herm_matrix_timestep_moving_view`
* @param g
* > The `herm_matrix_timestep_moving_view` from which the matrix element element is given
* @param j1
* > Row index of `herm_matrix_timestep_moving_view`
* @param j2
* > Column index of `herm_matrix_timestep_moving_view`
*/  
  template <typename T>
  void herm_matrix_timestep_moving_view<T>::set_matrixelement(int i1, int i2, herm_matrix_timestep_moving_view<T> &g,
							      int j1, int j2){
    int i, sij = i1 * size2_ + i2, tij = j1 * g.size2() + j2;
    assert(tc_==g.tc_);
    assert(i1>=0 && i1 <=size1_ -1);
    assert(i2>=0 && i2 <=size2_ -1);
    assert(j1>=0 && j1 <=g.size1_-1);
    assert(j2>=0 && j2 <=g.size2_-1);
    for(i =0; i<=tc_;i++)
      retptr(i)[sij]=g.retptr(i)[tij];
    for(i =0; i<=tc_;i++)
      lesptr(i)[sij]=g.lesptr(i)[tij];
  }
  // template <typename T>
  // template <class GG>
  // void herm_matrix_timestep_moving_view<T>::set_matrixelement(int i1, int i2, GG &g, int j1, int j2){
  //   //due to missing constructors following this scheme only for GG==herm_matrix_moving<T>. Extend?
  //   herm_matrix_timestep_moving_view<T> tmp(t0_,g);
  //   get_data(tmp)
  //     }
#if CNTR_USE_MPI ==1
  /// @private
  template <typename T>
  void my_mpi_reduce_moving(std::complex<T> *data, int len, int root){
    std::cerr << __PRETTY_FUNCTION__ << ", LEN=" << len
	      << " ... NOT DEFINED FOR THIS TYPE " << std::endl;
    exit(0);
  }
  /// @private
  template<>
  inline void my_mpi_reduce_moving<double>(std::complex<double> *data, int len, int root) {
    int tid,ntasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &tid);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    assert(root>=0 && root <= ntasks -1);
    assert(len>=0);
    if (tid == root) {
      MPI_Reduce(MPI_IN_PLACE, (double *)data, 2 * len,
		 MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    } else {
      MPI_Reduce((double *)data, (double *)data, 2 * len,
		 MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    }
  }
  /// @private
  template <typename T>
  void herm_matrix_timestep_moving_view<T>::MPI_Reduce(int root) {
    my_mpi_reduce_moving<T>(les_, (tc_ + 1) * element_size_, root);
    my_mpi_reduce_moving<T>(ret_, (tc_ + 1) * element_size_, root);
  }  
#endif
/// @private
/** \brief <b> Return single particle density matrix from `herm_matrix_timestep_moving_view`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return single particle density matrix (occupation) from `herm_matrix_timestep_moving_view`
* > Works for scalar.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param t0
* > Physical time
*/
  template <typename T>
  std::complex<T> herm_matrix_timestep_moving_view<T>::density_matrix(int t0){
    assert(t0==t0_);

    CPLX x1;
    x1= *lesptr(0);//Check this
    return CPLX(0.0,sig_)*x1;
  }
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
        CPLX w;                                                              \
        int r, s, dim = size1_;                                              \
        M.resize(dim, dim);                                                  \
        for (r = 0; r < dim; r++)                                            \
            for (s = 0; s < dim; s++) {                                      \
                w = x[s * dim + r];                                          \
                M(r, s) = CPLX(-w.real(), w.imag());                         \
            }                                                                \
    }

/** \brief <b> Return single particle density matrix from `herm_matrix_timestep_moving_view`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Return single particle density matrix (occupation) from `herm_matrix_timestep_moving_view`
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param t0
* > Physical time
* @param M
* > density matrix, returned as complex number
*/  
  template <typename T>
  template <class Matrix>
  void herm_matrix_timestep_moving_view<T>::density_matrix(int t0, Matrix &M){
    CPLX *x;
    x=lesptr(0);
    herm_matrix_READ_ELEMENT M *=CPLX(0.0,1.0*sig_);   
  }
#undef herm_matrix_READ_ELEMENT
#undef herm_matrix_READ_ELEMENT_MINUS_CONJ
#undef CPLX

  
}// namespace cntr
#endif
