#ifndef CNTR_HERM_MATRIX_TIMESTEP_MOVING_IMPL_H
#define CNTR_HERM_MATRIX_TIMESTEP_MOVING_IMPL_H

#include "cntr_herm_matrix_timestep_decl.hpp"
#include "cntr_herm_matrix_timestep_moving_decl.hpp"
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
  herm_matrix_timestep_moving<T>::herm_matrix_timestep_moving(){
    data_=0;
    les_=0;
    ret_=0;
    tc_=-1;
    t0_=0;
    size1_=0;
    size2_=0;
    element_size_=0;
    sig_=-1;
  }
  
  template <typename T>
  herm_matrix_timestep_moving<T>::~herm_matrix_timestep_moving(){
    if (data_!=0) delete [] data_;
  }
/** \brief <b> Initializes the `herm_matrix_timestep_moving` class for a square-matrix for fermions. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `herm_matrix_timestep_moving` class for a square-matrix.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tc
* > Cutoff time
* @param t0
* > Current physical time step
* @param size1
* > Matrix rank of the contour function. (size2=size1)
* @param sig
* > Set `sig = -1` for fermions or `sig = +1` for bosons.
*/
  
  template <typename T>
  herm_matrix_timestep_moving<T>::herm_matrix_timestep_moving(int tc,int t0,int size1,int sig){
    assert(-1<=tc);
    assert(1==sig*sig);
    // CNTR_ASSERT(MOVING_HERM_ASSERT_LEVEL,!(tc==-1 && size1>0),__PRETTY_FUNCTION__)
    tc_=tc;
    t0_=t0;
    sig_=sig;
    size1_=size1;
    size2_=size1;
    element_size_=size1*size1;
    if(tc_==-1){
      data_=0;
      les_=0;
      ret_=0;
    }else{
      // here tc>=0 AND size>0
      long ndata1=(tc_+1)*element_size_;
      data_ = new cplx [2*ndata1];
      ret_ = data_;
      les_ = data_+ndata1;
      memset(data_, 0, 2*sizeof(cplx)*ndata1);
    }
  }
/** \brief <b> Initializes the `herm_matrix_timestep_moving` class with the same layout as a given `herm_matrix_timestep_moving`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `herm_matrix_timestep_moving` class with the same value of 
* > the cutoff time `tc`, physical time `t0`, colum rank `size1`, row rank `size2` and
* > bosonic/fermionic symmetry `sig`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix_timestep_moving` according to which the class should be initialized
*/
  template <typename T>
  herm_matrix_timestep_moving<T>::herm_matrix_timestep_moving(const herm_matrix_timestep_moving &g){
    tc_=g.tc();
    t0_=g.t0();
    sig_=g.sig();
    size1_=g.size1();
    size2_=g.size2();
    element_size_=g.element_size_;
    if(tc_==-1){
      data_=0;
      les_=0;
      ret_=0;
    }else{
      // here tc>0 AND size>0
      long ndata1=(tc_+1)*element_size_;
      data_ = new cplx [2*ndata1];
      ret_ = data_;
      les_ = data_+ndata1;
      memcpy(data_, g.data_, 2*sizeof(cplx)*ndata1);
    }
  }
/** \brief <b> Initializes the `herm_matrix_timestep_moving` class form a `herm_matrix_moving`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Initializes the `herm_matrix_timestep_moving` class with the same value of 
* > the cutoff time `tc`, physical time `t0`, colum rank `size1`, row rank `size2` and
* > bosonic/fermionic symmetry `sig`. Copies the data from the timestep \f$ t_0-n \f$ from the 
* > 'herm_matrix_moving'
*
* <!-- ARGUMENTS
*      ========= -->
* @param n
* > The time index of the time step copied from `herm_matrix_moving`
* @param g
* > The `herm_matrix_timestep_moving` according to which the class should be initialized
*/
  template <typename T>
  herm_matrix_timestep_moving<T>::herm_matrix_timestep_moving(int n,herm_matrix_moving<T> &g){
    int tc=g.tc();
    assert(0<=n && n<=tc);
    tc_=g.tc();
    t0_=g.t0();
    sig_=g.sig();
    size1_=g.size1();
    size2_=g.size2();
    element_size_=g.element_size();
    if(tc_==-1){
      data_=0;
      les_=0;
      ret_=0;
    }else{
      // here tc>0 AND size>0
      long ndata1=(tc_+1)*element_size_;
      data_ = new cplx [2*ndata1];
      ret_ = data_;
      les_ = data_+ndata1;
      memcpy(ret_, g.retptr(n,0), sizeof(cplx)*ndata1);
      memcpy(les_, g.lesptr(n,0), sizeof(cplx)*ndata1);
    }
  }
  
  template <typename T>
  herm_matrix_timestep_moving<T> &  herm_matrix_timestep_moving<T>::operator=(const  herm_matrix_timestep_moving &g){
    if(this==&g) return *this;
    sig_=g.sig_;
    if( tc_!=g.tc_ || size1_!=g.size1_ || size2_!=g.size2_){
      // reallocate
      if (data_!=0) delete [] data_;
      tc_=g.tc_;
      t0_=g.t0();
      size1_=g.size1_;
      size2_=g.size2_;
      element_size_=g.element_size_;
      if(tc_>=0){
	// here tc>0 AND size>0
	long ndata1=(tc_+1)*element_size_;
	data_ = new cplx [2*ndata1];
	ret_ = data_;
	les_ = data_+ndata1;
      }else{
	data_=0;
	les_=0;
	ret_=0;
      }
    }
    if(tc_>=0){
      memcpy(data_, g.data_, 2*sizeof(cplx)*(tc_+1)*element_size_);
    }
    return *this;
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
  void herm_matrix_timestep_moving<T>::clear(void){
    if(tc_==-1) return;
    memset(data_, 0, 2*sizeof(cplx)*(tc_+1)*element_size_);
  }

/** \brief <b> Resizes `herm_matrix_timestep_moving` object with respect to the given cutoff time and the matrix size </b>.
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Resizes `herm_matrix_timestep_moving` class with respect to the given cutoff time 'tc' 
* > and the matrix size `size1`. Works for a square matrices.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tc
* > Cutoff time
* @param size1
* > size of the square matrix
*/
  template <typename T>
  void  herm_matrix_timestep_moving<T>::resize(int tc,int size1){
    if( tc!=tc_ || size1!=size1_){
      // reallocate
      if (data_!=0) delete [] data_;
      tc_=tc;
      size1_=size1;
      size2_=size1;
      element_size_=size1*size1;
      if(tc_>0){
	long ndata1=(tc_+1)*element_size_;
	data_ = new cplx [2*ndata1];
	ret_ = data_;
	les_ = data_+ndata1;
      }else{
	data_=0;
	les_=0;
	ret_=0;
      }
    }
    clear();
  }
  // READING ELEMENTS TO ANY MATRIX TYPE OR TO COMPLEX NUMBERS
  // (then only the (0,0) element is addressed for dim>0)
#define herm_matrix_timestep_moving_READ_ELEMENT {int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++) M(r,s)=x[r*dim+s];}
#define herm_matrix_timestep_moving_READ_ELEMENT_MINUS_CONJ {cplx w;int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++){ w=x[s*dim+r];M(r,s)=std::complex<T>(-w.real(),w.imag());}}
/** \brief <b> Returns the lesser component at a given time index. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the lesser component \f$ C^\mathrm{<}(t_0,t_0-j) \f$
* > at a given time \f$ t_0-j\f$ with \f$ j <= t_c\f$ for particular time step \f$ t_0 \f$
* > to a given matrix class.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Time index
* @param M
* > Matrix to which the lesser component is given.
*/
  template<typename T> template <class Matrix> 
  void herm_matrix_timestep_moving<T>::get_les(int j,Matrix &M){
    assert(0<=j && j<=tc_);
    cplx *x;
    x=lesptr(j);
    herm_matrix_timestep_moving_READ_ELEMENT
      }
/** \brief <b> Returns the retarded component at a given time index. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the retarded component \f$ C^\mathrm{R}(t_0,t_0-j) \f$
* > at a given time \f$ t_0-j\f$ with \f$ j <= t_c\f$ for particular time step \f$ t_0 \f$
* > to a given matrix class.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Time index
* @param M
* > Matrix to which the lesser component is given.
*/
  template<typename T> template <class Matrix> 
  void herm_matrix_timestep_moving<T>::get_ret(int j,Matrix &M){
    assert(0<=j && j<=tc_);
    cplx *x;
    x=retptr(j);
    herm_matrix_timestep_moving_READ_ELEMENT
      }
/** \brief <b> Returns the greater component at a given time index. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the greater component \f$ C^\mathrm{>}(t_0,t_0-j) \f$
* > at a given time \f$ t_0-j\f$ with \f$ j <= t_c\f$ for particular time step \f$ t_0 \f$
* > to a given matrix class.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Time index
* @param M
* > Matrix to which the lesser component is given.
*/
  template<typename T> template <class Matrix> 
  void herm_matrix_timestep_moving<T>::get_gtr(int j,Matrix &M){
    Matrix M1;
    get_ret(j,M);
    get_les(j,M1);
    M += M1;
  }
/** \brief <b> Returns the retarded component at a given time index. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the scalar-valued retarded component \f$ C^\mathrm{R}(t_0,t_0-j) \f$
* > at a given time \f$ t_0-j\f$ with \f$ j <= t_c\f$ for particular time step \f$ t_0 \f$
* > to a given complex number.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Time index
* @param x
* > complex number to which the retarded component is given.
*/
  template<typename T> 
  inline void herm_matrix_timestep_moving<T>::get_ret(int j,cplx &x){
    assert(0<=j && j<=tc_);
    x=*retptr(j);
  }
/** \brief <b> Returns the lesser component at a given time index. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the scalar-valued lesser component \f$ C^\mathrm{<}(t_0,t_0-j) \f$
* > at a given time \f$ t_0-j\f$ with \f$ j <= t_c\f$ for particular time step \f$ t_0 \f$
* > to a given complex number.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Time index
* @param x
* > complex number to which the retarded component is given.
*/
  template<typename T> 
  inline void herm_matrix_timestep_moving<T>::get_les(int j,cplx &x){
    assert(0<=j && j<=tc_);
    x=*lesptr(j);
  }
/** \brief <b> Returns the greater component at a given time index. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the scalar-valued greater component \f$ C^\mathrm{>}(t_0,t_0-j) \f$
* > at a given time \f$ t_0-j\f$ with \f$ j <= t_c\f$ for particular time step \f$ t_0 \f$
* > to a given complex number.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Time index
* @param x
* > complex number to which the retarded component is given.
*/
  template<typename T> 
  inline void herm_matrix_timestep_moving<T>::get_gtr(int j,cplx &x){
    cplx x1;
    get_ret(j,x);
    get_les(j,x1);
    x+=x1;
  }
/** \brief <b> Returns the density matrix at given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the scalar-valued density matrix (occupation, that is) at
* > given time step `t0`. Returns \f$ \rho(t) = i \eta C^<(t,t) \f$.
* > The return value is formally complex.
*
* <!-- ARGUMENTS
*      ========= -->
*/
  template<typename T> 
  std::complex<T> herm_matrix_timestep_moving<T>::density_matrix(){
    cplx x1;
    get_les(0,x1);
    return std::complex<T>(0.0,sig_)*x1;
  }
/** \brief <b> Returns the density matrix at given time step. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Returns the  matrix-valued density matrix (occupation, that is) at
* > given time step `t0`. Returns \f$ \rho(t) = i \eta C^<(t,t) \f$.
* > The return value is formally complex.
*
* <!-- ARGUMENTS
*      ========= -->
* @param M
* > Matrix to which the density matrix is given.
*/
  template<typename T> template<class Matrix> 
  void herm_matrix_timestep_moving<T>::density_matrix(Matrix &M){
    get_les(0,M);
    M *= std::complex<T>(0.0,1.0*sig_);
  }
/** \brief <b> Set the bosonic/fermionic symmetry. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* >Manipulates the internal bosonic/fermionic symmetry to match the given 'sig'.
* <!-- ARGUMENTS
*      ========= -->
* @param s
* > Given symmetry +1=bosonic, -1=fermionic.
*/
  template<typename T> 
  inline void herm_matrix_timestep_moving<T>::set_sig(int s){
    assert(s*s==1);
    sig_=s;
  }

/** \brief <b> Increase the value of the  `herm_matrix_timestep_moving` by `weight` \f$ * g(t) \f$ </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Increase the value of `herm_matrix_timestep_moving` by a value of `weight`*\f$ g(t)\f$,
* > where \f$ g(t)\f$ is a `herm_matrix_timestep_moving` and `weight` is a constant 
* > scalar number.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g1
* > The `herm_matrix_timestep_moving` which is added
* @param alpha
* > Scalar constant weight factor
*
*/
  template <typename T>
  void herm_matrix_timestep_moving<T>::incr_timestep(herm_matrix_timestep_moving<T> &g,std::complex<T> alpha){
    assert(size1_==g.size1());
    assert(size2_==g.size2());
    assert(tc_==g.tc());
    int ndata1=(tc_+1)*element_size_;
    for(int l=0;l<ndata1;l++){
      retptr(0)[l] +=  alpha*g.retptr(0)[l];
      lesptr(0)[l] +=  alpha*g.lesptr(0)[l];
    }
  }
/** \brief <b> Left-multiplication of the `herm_matrix_timestep_moving` with a moving contour function </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Left-multiplication of the `herm_matrix_timestep_moving` with a time dependent 
* > moving contour function \f$ F(t) \f$ 
* > i.e. it performs operation \f$G(t,t') \rightarrow w F(t)G(t,t')\f$
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param ft
* > the contour function \f$ F(t) \f$
* @param weight
* > some number (weight)
*/
  template <typename T>
  void herm_matrix_timestep_moving<T>::left_multiply(function_moving<T> &ft,T weight){
    assert(size1_==ft.size1());
    assert(size2_==ft.size2());
    assert(tc_==ft.tc());
    cplx *xtemp,*ftemp;
    xtemp=new cplx [element_size_];
    ftemp=ft.ptr(0);
    for(int j=0;j<=tc_;j++){
      element_mult<T,LARGESIZE>(size1_,xtemp,ftemp,retptr(j));
      element_smul<T,LARGESIZE>(size1_,xtemp,weight);
      element_set<T,LARGESIZE>(size1_,retptr(j),xtemp);
      element_mult<T,LARGESIZE>(size1_,xtemp,ftemp,lesptr(j));
      element_smul<T,LARGESIZE>(size1_,xtemp,weight);
      element_set<T,LARGESIZE>(size1_,lesptr(j),xtemp);
    }
    delete [] xtemp;
  }
/** \brief <b> Right-multiplication of the `herm_matrix_timestep_moving` with a moving contour function </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Right-multiplication of the `herm_matrix_timestep_moving` with a time dependent 
* > moving contour function \f$ F(t) \f$
* > i.e. it performs operation \f$G(t,t') \rightarrow w G(t,t')F(t')\f$
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param ft
* > the contour function \f$ F(t) \f$
* @param weight
* > some number (weight)
*/
  template <typename T>
  void herm_matrix_timestep_moving<T>::right_multiply(function_moving<T> &ft,T weight){
    assert(size1_==ft.size1());
    assert(size2_==ft.size2());
    assert(tc_==ft.tc());
    cplx *xtemp;
    xtemp=new cplx [element_size_];
    for(int j=0;j<=tc_;j++){
      element_mult<T,LARGESIZE>(size1_,xtemp,retptr(j),ft.ptr(j));
      element_smul<T,LARGESIZE>(size1_,xtemp,weight);
      element_set<T,LARGESIZE>(size1_,retptr(j),xtemp);
      element_mult<T,LARGESIZE>(size1_,xtemp,lesptr(j),ft.ptr(j));
      element_smul<T,LARGESIZE>(size1_,xtemp,weight);
      element_set<T,LARGESIZE>(size1_,lesptr(j),xtemp);
    }
    delete [] xtemp;
  }
/** \brief <b> Multiply  `herm_matrix_timestep_moving` with scalar `weight`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Multiply `herm_matrix_timestep_moving` with a scalar.
* > Works for scalar or matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param weight
* > The `template argument` multiplication factor
*/
  template <typename T>
  void herm_matrix_timestep_moving<T>::smul(T weight)
  {
    int ndata = element_size_ * (tc_ + 1) * 2;
    for (int i = 0; i < ndata; ++i)
      {
	data_[i] *= weight;
      }
  }
/** \brief <b> Multiply  `herm_matrix_timestep_moving` with a complex scalar `weight`.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Multiply `herm_matrix_timestep_moving` with a scalar.
* > Works for scalar or matrix contour objects.
*
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param weight
* > The complex multiplication factor
*/
  template <typename T>
  void herm_matrix_timestep_moving<T>::smul(cplx weight)
  {
    int ndata = element_size_ * (tc_ + 1) * 2;
    for (int i = 0; i < ndata; ++i)
      {
	data_[i] *= weight;
      }
  }
/** \brief <b> Increase the value of the  `herm_matrix_timestep_moving` by `weight` \f$ * f(t) \f$ </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Increase the value of `herm_matrix_timestep_moving` by a value of `weight`*\f$ f(t)\f$,
* > where \f$ f(t)\f$ is a `function_moving` and `weight` is a 
* > constant scalar number
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g1
* > The `function_moving` which is added
* @param weight
* > Scalar constant multiplication factor
*
*/  
  template <typename T>
  void herm_matrix_timestep_moving<T>::incr(function_moving<T> &g1, T weight){
    int m;
    int total_size_=2*(tc_ + 1)*element_size_;
    assert(g1.size1_ == size1_ && g1.tc_ == tc_);
    for (m = 0; m < total_size_; m++)
      data_[m] += weight * g1.data_[m];
  }
/** \brief <b> Increase the value of the  `herm_matrix_timestep_moving` by `weight` \f$ * g(t) \f$ </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Increase the value of `herm_matrix_timestep_moving` by a value of `weight`*\f$ g(t)\f$,
* > where \f$ g(t)\f$ is a `herm_matrix_timestep_moving` and `weight` is a 
* > constant scalar number
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g1
* > The `herm_matrix_timestep_moving` which is added
* @param weight
* > Scalar constant multiplication factor
*
*/
  template <typename T>
  void herm_matrix_timestep_moving<T>::incr(herm_matrix_timestep_moving<T> &g1, T weight){
    int m;
    int total_size_=2*(tc_ + 1)*element_size_;
    assert(g1.size1_ == size1_ && g1.tc_ == tc_);
    for (m = 0; m < total_size_; m++)
      data_[m] += weight * g1.data_[m];
  }

#if CNTR_USE_MPI == 1
/** \brief <b> MPI reduce for the `herm_matrix_timestep_moving` </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > MPI reduce for the `herm_matrix_timestep_moving` to the `root`
* > Works for scalar or square-matrix contour objects.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param root
* > Index of root
*/
  template <typename T>
  void herm_matrix_timestep_moving<T>::MPI_Reduce(int root) {
    int taskid;
    //int len1 = 2 * ((tstp_ + 1) * 2 + (ntau_ + 1)) * size1_ * size1_; This is the original version from timestep check if the reduction to the data window is correct!!
    int len1 = ((tc_ + 1) * 2 ) * element_size_;
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
      std::cerr << "herm_matrix_timestep_moving<T>::MPI_Reduce only for double "
		<< std::endl;
      exit(0);
    }
  }
#endif
} // namespace cntr

#endif  // CNTR_HERM_MATRIX_TIMESTEP_IMPL_H
