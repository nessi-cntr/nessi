#ifndef CNTR_HERM_MATRIX_MOVING_IMPL_H
#define CNTR_HERM_MATRIX_MOVING_IMPL_H

#include "cntr_herm_matrix_decl.hpp"
#include "cntr_herm_matrix_moving_decl.hpp"
#include "cntr_herm_pseudo_decl.hpp"
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
  herm_matrix_moving<T>::herm_matrix_moving(){
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
  herm_matrix_moving<T>::~herm_matrix_moving(){
    if (data_!=0) delete [] data_;
    if (ret_!=0)  delete [] ret_;
    if (les_!=0)  delete [] les_;
  }
  /** \brief <b> Initializes the `herm_matrix_moving` class for a square-matrix two-time contour function.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Initializes the `herm_matrix_moving` class for a square-matrix two-time contour function.
* 
* <!-- ARGUMENTS
*      ========= -->
*
* @param tc
* > Number of time steps
* @param t0
* > Real time step of the object
* @param size1
* > Matrix rank of the contour function
* @param sig
* > Set `sig = -1` for fermions or `sig = +1` for bosons.
*/
  template <typename T>
  herm_matrix_moving<T>::herm_matrix_moving(int tc,int t0,int size1,int sig){
    assert(-1<=tc);
    assert(1==sig*sig);
    // CNTR_ASSERT(herm_matrix_moving_ASSERT_LEVEL,!(tc==-1 && size1>0),__PRETTY_FUNCTION__)
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
      long ndata2=(tc_+1)*(tc_+1)*element_size_;
      long ndata1=(tc_+1)*element_size_;
      data_ = new cplx [2*ndata2];
      ret_ = new cplx* [tc_+1];
      les_ = new cplx* [tc_+1];
      memset(data_, 0, 2*sizeof(cplx)*ndata2);
      for(int t=0;t<=tc_;t++){
	ret_[t]=data_+t*ndata1;
	les_[t]=data_+(t+tc_+1)*ndata1;
      }
    }
  }
  
/** \brief <b> Initializes the `herm_matrix_moving` class with the same layout as a given `herm_matrix_moving`.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Initializes the `herm_matrix` class with the same number of time steps `tc`,
* real time step `t0`, matrix rank `size1` and
* bosonic/fermionic symmetry `sig`. Works for scalar or square-matrix contour objects
* only.
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `herm_matrix_moving` according to which the class should be initialized
*/
  template <typename T>
  herm_matrix_moving<T>::herm_matrix_moving(const herm_matrix_moving &g){
    tc_=g.tc_;
    t0_=g.t0_;
    sig_=g.sig_;
    size1_=g.size1_;
    size2_=g.size2_;
    element_size_=g.element_size_;
    if(tc_==-1){
      data_=0;
      les_=0;
      ret_=0;
    }else{
      // here tc>0 AND size>0
      long ndata2=(tc_+1)*(tc_+1)*element_size_;
      long ndata1=(tc_+1)*element_size_;
      data_ = new cplx [2*ndata2];
      ret_ = new cplx* [tc_+1];
      les_ = new cplx* [tc_+1];
      memcpy(data_, g.data_, 2*sizeof(cplx)*ndata2);
      // correctly redirect the pointers
      for(int t=0;t<=tc_;t++){
	ret_[t]=data_+(g.ret_[t]-g.data_);
	les_[t]=data_+(g.les_[t]-g.data_);
      }
    }
  }
  template <typename T>
  herm_matrix_moving<T> &  herm_matrix_moving<T>::operator=(const  herm_matrix_moving &g){
    if(this==&g) return *this;
    sig_=g.sig_;
    t0_=g.t0_;
    if( tc_!=g.tc_ || size1_!=g.size1_ || size2_!=g.size2_){
      // reallocate
      if (data_!=0) delete [] data_;
      if (ret_!=0)  delete [] ret_;
      if (les_!=0)  delete [] les_;
      tc_=g.tc_;
      size1_=g.size1_;
      size2_=g.size2_;
      element_size_=g.element_size_;
      if(tc_>=0){
	// here tc>0 AND size>0
	long ndata2=(tc_+1)*(tc_+1)*element_size_;
	data_ = new cplx [2*ndata2];
	ret_ = new cplx* [tc_+1];
	les_ = new cplx* [tc_+1];
      }else{
	data_=0;
	les_=0;
	ret_=0;
      }
    }
    if(tc_>=0){
      memcpy(data_, g.data_, 2*sizeof(cplx)*(tc_+1)*(tc_+1)*element_size_);
      for(int t=0;t<=tc_;t++){
	ret_[t]=data_+(g.ret_[t]-g.data_);
	les_[t]=data_+(g.les_[t]-g.data_);
      }
    }
    return *this;
  }
/** \brief <b> Set all data to zero for the `herm_matrix_moving` class</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set all data to zero for the `herm_matrix_moving` class.
*
* <!-- ARGUMENTS
*      ========= -->
*
*/  
  template <typename T>
  void herm_matrix_moving<T>::clear(void){
    if(tc_==-1) return;
    memset(data_, 0, 2*sizeof(cplx)*(tc_+1)*(tc_+1)*element_size_);
    long ndata1=(tc_+1)*element_size_;
    for(int t=0;t<=tc_;t++){
      ret_[t]=data_+t*ndata1;
      les_[t]=data_+(t+tc_+1)*ndata1;
    }
  }
/** \brief <b> Set the current timestep of the `herm_matrix_moving` class object</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set the current timestep `t0` of the `herm_matrix_moving` class to t0.
*
* <!-- ARGUMENTS
*      ========= -->
* @param t0
* > the new internal timestep t0.
*/
  template <typename T>
  void herm_matrix_moving<T>::set_t0(int t0){
    assert(tc_<=t0);
    t0_=t0;
  }
/** \brief <b> Resizes `herm_matrix_moving` object with respect to the number of
 * time points `tc` or the matrix size `size1`.  </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 *
 * > Resizes `herm_matrix_moving` class with respect to number of time steps `tc`. 
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tc
 * > New number of time steps.
 * @param size1
 * > size of the matrix
 */
  template <typename T>
  void  herm_matrix_moving<T>::resize(int tc,int size1){
    if( tc!=tc_ || size1!=size1_){
      // reallocate
      if (data_!=0) delete [] data_;
      if (ret_!=0)  delete [] ret_;
      if (les_!=0)  delete [] les_;
      tc_=tc;
      size1_=size1;
      size2_=size1;
      element_size_=size1*size1;
      if(tc_>0){
	long ndata2=(tc_+1)*(tc_+1)*element_size_;
	long ndata1=(tc_+1)*element_size_;
	data_ = new cplx [2*ndata2];
	ret_ = new cplx* [tc_+1];
	les_ = new cplx* [tc_+1];
	memset(data_, 0, 2*sizeof(cplx)*ndata2);
	for(int t=0;t<=tc_;t++){
	  ret_[t]=data_+t*ndata1;
	  les_[t]=data_+(t+tc_+1)*ndata1;
	}
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
#define herm_matrix_moving_READ_ELEMENT {int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++) M(r,s)=x[r*dim+s];}
#define herm_matrix_moving_READ_ELEMENT_MINUS_CONJ {cplx w;int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++){ w=x[s*dim+r];M(r,s)=std::complex<T>(-w.real(),w.imag());}}
/** \brief <b> Returns the lesser component at given time indices.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the lesser component \f$ C^<(t_0-i,t_0-i-j) \f$ for the internal time \f$ t_0 \f$
* > and given indices \f$ i \f$ and \f$ j \f$ to a given matrix class. 
* > Note that the access is only valid on the upper triangular matrix 
* > for \f$ i>=0 \f$ and \f$j>=0 \f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > time index \f$ i \f$ .
* @param j
* > time index \f$ j \f$ .
* @param M
* > Matrix to which the lesser component is given.
*/
  template<typename T> template <class Matrix> 
  void herm_matrix_moving<T>::get_les(int i,int j,Matrix &M) const{
    assert(0<=i && i<=tc_);
    assert(0<=j && j<=tc_);
    cplx *x;
    x=lesptr(i,j);
    herm_matrix_moving_READ_ELEMENT
      }
  /** \brief <b> Returns the retarded component at given time indices.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the lesser component \f$ C^\mathrm{R}(t_0-i,t_0-i-j) \f$ for the internal 
* > time \f$ t_0 \f$ and given indices \f$ i \f$ and \f$ j \f$ to a given matrix class. 
* > Note that the access is only valid on the upper triangular matrix 
* > for \f$ i>=0 \f$ and \f$j>=0 \f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > time index \f$ i \f$ .
* @param j
* > time index \f$ j \f$ .
* @param M
* > Matrix to which the retarded component is given.
*/
  template<typename T> template <class Matrix> 
  void herm_matrix_moving<T>::get_ret(int i,int j,Matrix &M) const{
    assert(0<=i && i<=tc_);
    assert(0<=j && j<=tc_);
    cplx *x;
    x=retptr(i,j);
    herm_matrix_moving_READ_ELEMENT
      }
/** \brief <b> Returns the greater component at given time indices.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the greater component \f$ C^>(t_0-i,t_0-i-j) \f$ for the internal 
* > time \f$ t_0 \f$ and given indices \f$ i \f$ and \f$ j \f$ to a given matrix class. 
* > Note that the access is only valid on the upper triangular matrix 
* > for \f$ i>=0 \f$ and \f$j>=0 \f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > time index \f$ i \f$ .
* @param j
* > time index \f$ j \f$ .
* @param M
* > Matrix to which the greater component is given.
*/
  template<typename T> template <class Matrix> 
  void herm_matrix_moving<T>::get_gtr(int i,int j,Matrix &M) const{
    Matrix M1;
    get_ret(i,j,M);
    get_les(i,j,M1);
    M += M1;
  }
/** \brief <b> Returns the retarded component at given time indices.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the lesser component \f$ C^\mathrm{R}(t_0-i,t_0-i-j) \f$ for the internal 
* > time \f$ t_0 \f$ and given indices \f$ i \f$ and \f$ j \f$ to a given matrix class. 
* > Note that the access is only valid on the upper triangular matrix 
* > for \f$ i>=0 \f$ and \f$j>=0 \f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > time index \f$ i \f$ .
* @param j
* > time index \f$ j \f$ .
* @param x
* > The value of the retarded component.
*/
  template<typename T> 
  inline void herm_matrix_moving<T>::get_ret(int i,int j,cplx &x) const{ 
    assert(0<=i && i<=tc_);
    assert(0<=j && j<=tc_);
    x=*retptr(i,j);
  }
  /** \brief <b> Returns the lesser component at given time indices.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the lesser component \f$ C^<(t_0-i,t_0-i-j) \f$ for the internal time \f$ t_0 \f$
* > and given indices \f$ i \f$ and \f$ j \f$ to a given matrix class. 
* > Note that the access is only valid on the upper triangular matrix 
* > for \f$ i>=0 \f$ and \f$j>=0 \f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > time index \f$ i \f$ .
* @param j
* > time index \f$ j \f$ .
* @param x
* > The value of the retarded component.
*/
  template<typename T> 
  inline void herm_matrix_moving<T>::get_les(int i,int j,cplx &x) const{ 
    assert(0<=i && i<=tc_);
    assert(0<=j && j<=tc_);
    x=*lesptr(i,j);
  }
  /** \brief <b> Returns the greater component at given time indices.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the greater component \f$ C^>(t_0-i,t_0-i-j) \f$ for the internal 
* > time \f$ t_0 \f$ and given indices \f$ i \f$ and \f$ j \f$ to a given matrix class. 
* > Note that the access is only valid on the upper triangular matrix 
* > for \f$ i>=0 \f$ and \f$j>=0 \f$.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > time index \f$ i \f$ .
* @param j
* > time index \f$ j \f$ .
* @param x
* > The value of the greater component.
*/ 
  template<typename T> 
  inline void herm_matrix_moving<T>::get_gtr(int i,int j,cplx &x) const{
    cplx x1;
    get_ret(i,j,x);
    get_les(i,j,x1);
    x+=x1;
  }
  /** \brief <b> Returns the density matrix at given time index. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Returns the scalar-valued density matrix (occupation, that is) at
* given time index `i`. Since `i` is always on the real time branch of the Keldysh contour 
* the function returns \f$ \rho(t) = i \eta C^<(t,t) \f$.
* The return value is formally complex.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > The time index at which the density matrix is returned.
*/ 
  template<typename T> 
  std::complex<T> herm_matrix_moving<T>::density_matrix(int i){
    cplx x1;
    get_les(i,0,x1);
    return std::complex<T>(0.0,sig_)*x1;
  }
  /** \brief <b> Returns the density matrix at given time index. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Returns the matrix-valued density matrix (occupation, that is) at
* given time index `i`. Since `i` is always on the real time branch of the Keldysh contour 
* the function returns \f$ \rho(t) = i \eta C^<(t,t) \f$.
* The return value is formally complex.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > The time index at which the density matrix is returned.
* @param M
* > The density matrix at time index `i`. 
*/
  template<typename T> 
  template<class Matrix> void herm_matrix_moving<T>::density_matrix(int i,Matrix &M){
    get_les(i,0,M);
    M *= std::complex<T>(0.0,1.0*sig_);
  }
/* #######################################################################################
#
#   INPUT/OUTPUT FROM/TO FILES
#
########################################################################################*/

/** \brief <b> Saves the `herm_matrix_moving` to file in text format. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Saves the `herm_matrix_moving` to file in plain text format, which can
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
  template<typename T> 
  void herm_matrix_moving<T>::print_to_file(const char *file,int precision){
    int i,j,l,sg=element_size_;
    std::ofstream out;
    out.open(file,std::ios::out);
    out.precision(precision);
    out << "# " << t0_ << " " << tc_ << " " << size1_ << " " << " " << sig_ << std::endl;
    if(tc_>=0){
      for(i=0;i<=tc_;i++){
	for(j=0;j<=tc_;j++){
	  out << "ret: " << i << " " << j;
	  for(l=0;l<sg;l++) out << " " << retptr(i,j)[l].real() << " " << retptr(i,j)[l].imag();
	  out << std::endl;
	}
	out << std::endl;
      }
      out << std::endl;
      for(i=0;i<=tc_;i++){
	for(j=0;j<=tc_;j++){
	  out << "les: " << i << " " << j;
	  for(l=0;l<sg;l++) out << " " << lesptr(i,j)[l].real() << " " << lesptr(i,j)[l].imag();
	  out << std::endl;
	}
	out << std::endl;
      }
      out << std::endl;
    }
    out.close();
  }
 /** \brief <b> Reads the `herm_matrix_moving` from file in text format. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Reads the `herm_matrix_moving` to file in plain text format, which happens
 * > been created with the corresponding routine `herm_matrix::print_to_file`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param file
 * > The input file name.
 */
  template<typename T> 
  void herm_matrix_moving<T>::read_from_file(const char *file){
    int i,t0,tc,j,l,size1,sg,sig;
    double real, imag;
    std::string s;
    std::ifstream out;
    out.open(file,std::ios::in);
    if(!(out >> s >> t0 >> tc >> size1 >> sig)){
      std::cerr << "read G from file " << file << " error in file" << std::endl; 
      abort();
    }
    if(tc!=tc_ || size1!=size1_) resize(tc,size1);
    set_t0(t0);
    sig_=sig;
    sg=element_size_;
    if(tc_>=0){
      for(i=0;i<=tc_;i++){
	for(j=0;j<=tc_;j++){
	  out >> s >> s >> s ;
	  for(l=0;l<sg;l++){
	    if(!( out >> real >> imag )){
	      std::cerr << "read G from file " << file << " error at ret (" << i<< "," << j << ")"<< std::endl; 
	      abort();
	    }
	    retptr(i,j)[l] = std::complex<T>(real, imag);
	  }
	}
      }
      for(i=0;i<=tc_;i++){
	for(j=0;j<=tc_;j++){
	  out >> s >> s >> s ;
	  for(l=0;l<sg;l++){
	    if(!( out >> real >> imag )){
	      std::cerr << "read G from file " << file << " error at ret (" << i<< "," << j << ")"<< std::endl; 
	      abort();
	    }
	    lesptr(i,j)[l] = std::complex<T>(real, imag);
	  }
	}
      }
    }
    out.close();
  }
/** \brief <b> Propagates the `herm_matrix_moving` timeslices by one index.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Shifts the timeslices in the `herm_matrix_moving` truncated time window by one index.
* > The final timeslice at `tc` is shifted cyclically to the index `0`.
*
* <!-- ARGUMENTS
*      ========= -->
*
*/
  template <typename T>
  void herm_matrix_moving<T>::forward(void){
    if(tc_>0){
      cplx* tmp1=les_[tc_];
      cplx* tmp2=ret_[tc_];
      for(int t=tc_;t>0;t--){
	les_[t]=les_[t-1];
	ret_[t]=ret_[t-1];
      }
      les_[0]=tmp1;
      ret_[0]=tmp2;
      t0_=t0_+1;
    }
    
  }
/* #######################################################################################
#
#   SIMPLE OPERATIONS ON TIMESTEPS
#   NOTE: i IS A TIME INDEX RELATED TO THE PHYSICAL TIME t0-i
#
########################################################################################*/
/** \brief <b> Sets all components at time index `i` to zero. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets all components of the `herm_matrix_moving` at time index `i` to zero.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param i
 * > The time index at which the components are set to zero.
 *
 */
  template <typename T>
  void herm_matrix_moving<T>::clear_timestep(int i){
    for(int t1=0;t1<=tc_;t1++){
      element_set_zero<T,LARGESIZE>(size1_,retptr(i,t1));
      element_set_zero<T,LARGESIZE>(size1_,lesptr(i,t1));
    }
  }
  /** \brief <b> Sets all components at timestep with index `i` to the components of
 *  a given `herm_matrix_timestep_view`. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets all components of the `herm_matrix_moving` at time index `i` to
 * > the components of given `herm_matrix_timestep_view`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param i
 * > The time index at which the components are set.
 *
 * @param g
 * > The `herm_matrix_timestep_view` from which the time step is copied.
 *
 * @param gcc
 * > The hermetian conjugate of the `herm_matrix_timestep_view` from 
 * > which the time step is copied.
 *
 */


  template <typename T>
  void herm_matrix_moving<T>::set_timestep(int i,herm_matrix_timestep_view<T> &g,herm_matrix_timestep_view<T> &gcc){
    // this.ret(i,j)=g.ret(tstp,i-j) for j=0...min(tc,i)
    // this.les(i,j)=g.les(tstp,i-j) for j=0...min(tc,i)
    // remaining entries are filled with zeros, if any
    assert(size1_==g.size1());
    assert(0<=i && i <=tc_);
    int t1;
    int tstp=g.tstp();
    clear_timestep(i);
    int smax=(tstp > tc_ ? tc_ : tstp); // note that at this point t>=0
    for(t1=0;t1<=smax;t1++){
      element_set<T,LARGESIZE>(size1_,retptr(i,t1),g.retptr(tstp-t1));
      element_minusconj<T,LARGESIZE>(size1_,lesptr(i,t1),gcc.lesptr(tstp-t1));
    }
  }
    /** \brief <b> Sets all components at timestep with index `i` to the components of
 *  a given `herm_matrix`. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets all components of the `herm_matrix_moving` at time index `i` to
 * > the components of given `herm_matrix`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param i
 * > The time index at which the components are set.
 *
 * @param tstp
 * > The timestep of the `herm_matrix` from which the components are copied from.
 *
 * @param g
 * > The `herm_matrix` from which the time step is copied.
 *
 * @param gcc
 * > The hermetian conjugate of the `herm_matrix` from 
 * > which the time step is copied.
 *
 */
  template <typename T>
  void herm_matrix_moving<T>::set_timestep(int i,int tstp,herm_matrix<T> &g,herm_matrix<T> &gcc){
    herm_matrix_timestep_view<T> tmp(tstp,g);
    herm_matrix_timestep_view<T> tmp1(tstp,gcc);
    set_timestep(i,tmp,tmp1);
  }
/** \brief <b>Set all values of `herm_matrix_moving` class from a `herm_matrix` class object. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Copies the tc timesteps of a `herm_matrix` object into this `herm_matrix_moving` object. 
* > Starting at \f$ tstp\rightarrow 0\f$ to \f$tstp-tc\rightarrow tc\f$
*
* <!-- ARGUMENTS
*      ========= -->
* @param tstp
* > Latest timestep of the `herm_matrix` that is transferred to the `herm_matrix_moving`.
* @param g
* > The `herm_matrix` from which the time step is copied.
* @param gcc
* > The hermetian conjugate of the `herm_matrix` from 
* > which the time step is copied.
*/
  template <typename T>
  void herm_matrix_moving<T>::set_from_G_backward(int tstp,herm_matrix<T> &g,herm_matrix<T> &gcc){
    assert(tc_<=tstp && tstp <= g.nt());
    for(int i=0;i<=tc_;i++) set_timestep(i,tstp-i,g,gcc);
  }
/** \brief <b>Set all values of `herm_matrix_moving` class from a `herm_matrix` class object. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Copies the tc timesteps of a `herm_matrix` object into this `herm_matrix_moving` object. 
* > Starting at \f$ tstp\rightarrow 0\f$ to \f$tstp-tc\rightarrow tc\f$
*
* <!-- ARGUMENTS
*      ========= -->
* @param tstp
* > Latest timestep of the `herm_matrix` that is transferred to the `herm_matrix_moving`.
* @param g
* > The `herm_matrix` from which the time step is copied.
*/
  template <typename T>
  void herm_matrix_moving<T>::set_from_G_backward(int tstp,herm_matrix<T> &g){
    set_from_G_backward(tstp,g,g);
  }


/** \brief <b>Set all values of `herm_matrix_moving` class from a `herm_pseudo` class object. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Copies the tc timesteps of a `herm_pseudo` object into this `herm_matrix_moving` object. 
* > Starting at \f$ tstp\rightarrow 0\f$ to \f$tstp-tc\rightarrow tc\f$
*
* <!-- ARGUMENTS
*      ========= -->
* @param tstp
* > Latest timestep of the `herm_matrix` that is transferred to the `herm_matrix_moving`.
* @param g
* > The `herm_pseudo` from which the time step is copied.
*/
  template <typename T>
  void herm_matrix_moving<T>::set_from_G_backward(int tstp, const herm_pseudo<T> &g){
    for(int i = 0; i <= tc_; ++i) {
	assert(size1_==g.size1());
	int t1;
	clear_timestep(i);
	int smax=(tstp > tc_ ? tc_ : tstp); // note that at this point t>=0
	for(t1=0;t1<=smax;t1++){
	  element_set<T,LARGESIZE>(size1_,retptr(i,t1),g.retptr(tstp - i, tstp - i - t1));
	  element_minusconj<T,LARGESIZE>(size1_,lesptr(i,t1),g.lesptr(tstp - i - t1, tstp - i)); 
	}
    }
  }
  /** \brief <b> Sets all components at timestep with index `i` to the components of
 *  a given `herm_matrix_timestep_moving`. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets all components of the `herm_matrix_moving` at time index `i` to
 * > the components of given `herm_matrix_timestep_moving`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param i
 * > The time index at which the components are set.
 *
 * @param g
 * > The `herm_matrix_timestep_moving` from which the timestep is copied.
 *
 */

  template <typename T>
  void herm_matrix_moving<T>::set_timestep(int i,herm_matrix_timestep_moving<T> &g){
    assert(size1_=g.size1());
    assert(size2_=g.size2());
    assert(tc_=g.tc());
    assert(0<=i && i<=tc_);

    int ndata1=(tc_+1)*element_size_;
    memcpy(retptr(i,0),g.retptr(0), sizeof(cplx)*ndata1);
    memcpy(lesptr(i,0),g.lesptr(0), sizeof(cplx)*ndata1);
  }
/** \brief <b> Copies all components from timestep with index `i` to a given 
 * `herm_matrix_timestep_moving`. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Copies all components form the `herm_matrix_moving` at time index `i` to
 * > the components of a given `herm_matrix_timestep_moving`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param i
 * > The time index at which the components are set.
 *
 * @param g
 * > The `herm_matrix_timestep_moving` to which the timestep is copied.
 *
 */
  template <typename T>
  void herm_matrix_moving<T>::get_timestep(int i,herm_matrix_timestep_moving<T> &g){
    assert(0<=i && i<=tc_);
    g.resize(tc_,size1_);
    g.set_sig(sig_);
    int ndata1=(tc_+1)*element_size_;
    memcpy(g.retptr(0),retptr(i,0), sizeof(cplx)*ndata1);
    memcpy(g.lesptr(0),lesptr(i,0), sizeof(cplx)*ndata1);
  }
/** \brief <b> Sets all components at timestep with index `i` to the components 
 *  of a timestep with index `j` of a given `herm_matrix_moving`. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Sets all components of this `herm_matrix_moving` at time index `i` to
 * > the components of given `herm_matrix_moving` at time index `j`.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param i
 * > The time index at which the components are set.
 *
 * @param g
 * > The `herm_matrix_timestep_moving` from which the timestep is copied.
 *
 */
  template <typename T>
  void herm_matrix_moving<T>::set_timestep(int i,herm_matrix_moving<T> &g,int j){
    assert(size1_=g.size1());
    assert(size2_=g.size2());
    assert(tc_=g.tc());
    assert(0<=i && i<=tc_);
    assert(0<=j && j<=tc_);

    int ndata1=(tc_+1)*element_size_;
    memcpy(retptr(i,0),g.retptr(j,0), sizeof(cplx)*ndata1);
    memcpy(lesptr(i,0),g.lesptr(j,0), sizeof(cplx)*ndata1);
  }

/** \brief <b> Adds a timestep of a given `herm_matrix_moving` with given weight to the `herm_matrix_moving`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Performs the operation \f$C \rightarrow C + \alpha A\f$, where \f$C\f$ is the
* > `herm_matrix_moving`, \f$A\f$ is a timestep obtained from a given `herm_matrix_moving` 
* > at time index `j` and \f$\alpha\f$ is a complex weight. 
* > The operation is performed at given time index `i`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] The time index where the `herm_matrix_moving` is increased.
* @param g
* > [herm_matrix_moving] The `herm_matrix_moving` which is added to the `herm_matrix`.
* @param j
* > [int] The time index of the timestep from which the data are added.
* @param alpha
* > [complex<T>] The weight in front of `g`.
*/
  template <typename T>
  void herm_matrix_moving<T>::incr_timestep(int i,herm_matrix_moving<T> &g,int j,cplx alpha){
    assert(size1_=g.size1());
    assert(size2_=g.size2());
    assert(tc_=g.tc());
    assert(0<=i && i<=tc_);
    assert(0<=j && j<=tc_);

    int ndata1=(tc_+1)*element_size_;
    for(int l=0;l<ndata1;l++){
      retptr(i,0)[l] +=  alpha*g.retptr(j,0)[l];
      lesptr(i,0)[l] +=  alpha*g.lesptr(j,0)[l];
    }
  }
/** \brief <b> Adds a timestep of a given `herm_matrix_timestep_moving` with given weight to the `herm_matrix_moving`. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Performs the operation \f$C \rightarrow C + \alpha A\f$, where \f$C\f$ is the
* > `herm_matrix_moving`, \f$A\f$ is  a given `herm_matrix_timestep_moving` and \f$\alpha\f$
* > is a complex weight. 
* > The operation is performed at given time index `i`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > [int] The time index where the `herm_matrix_moving` is increased.
* @param g
* > [herm_matrix_moving] The `herm_matrix_moving` which is added to the `herm_matrix`.
* @param alpha
* > [complex<T>] The weight in front of `g`.
*/
  template <typename T>
  void herm_matrix_moving<T>::incr_timestep(int i,herm_matrix_timestep_moving<T> &g,cplx alpha){
    assert(size1_=g.size1());
    assert(size2_=g.size2());
    assert(tc_=g.tc());
    assert(0<=i && i<=tc_);

    int ndata1=(tc_+1)*element_size_;
    for(int l=0;l<ndata1;l++){
      retptr(i,0)[l] +=  alpha*g.retptr(0)[l];
      lesptr(i,0)[l] +=  alpha*g.lesptr(0)[l];
    }
  }


/** \brief <b> Left-multiplies the `herm_matrix_moving` with contour function. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Performs the operation \f$C(t,t^\prime) \rightarrow w F(t)C(t,t^\prime)\f$, where \f$C(t,t^\prime)\f$ is the
* `herm_matrix`, \f$F(t)\f$ is a `function` given in pointer format and \f$w\f$ is a real weight. The operation is performed at given time index `i`.
*
* <!-- ARGUMENTS
*      ========= -->
* @param ft
* > [function] The truncated contour function \f$F(t)\f$. 
* @param weight
* > [T] The weight as above.
* @param i
* > [int] The time index at which \f$F(t)\f$ and \f$C(t,t^\prime)\f$ are multiplied.
*/
  template <typename T>
  void herm_matrix_moving<T>::left_multiply(function_moving<T> &ft,T weight,int i){
    assert(size1_==ft.size1());
    assert(size2_==ft.size2());
    assert(size1_==size2_);
    assert(tc_==ft.tc());

    cplx *xtemp,*ftemp;
    xtemp=new cplx [element_size_];
    ftemp=ft.ptr(0);
    for(int j=0;j<=tc_;j++){
      element_mult<T,LARGESIZE>(size1_,xtemp,ftemp,retptr(i,j));
      element_smul<T,LARGESIZE>(size1_,xtemp,weight);
      element_set<T,LARGESIZE>(size1_,retptr(i,j),xtemp);
      element_mult<T,LARGESIZE>(size1_,xtemp,ftemp,lesptr(i,j));
      element_smul<T,LARGESIZE>(size1_,xtemp,weight);
      element_set<T,LARGESIZE>(size1_,lesptr(i,j),xtemp);
    }
    delete [] xtemp;
  }
  /** \brief <b> Right-multiplies the `herm_matrix_moving` with contour function. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Performs the operation \f$C(t,t^\prime) \rightarrow w C(t,t^\prime)F(t^\prime)\f$, 
* where \f$C(t,t^\prime)\f$ is the `herm_matrix_moving`, \f$F(t)\f$ is a contour 
* `function_moving` and \f$w\f$ is a real weight. The operation is performed
* at given time index `i`.
*
* <!-- ARGUMENTS
*      ========= -->
* @param ft
* > [function] The truncated contour function \f$F(t)\f$. 
* @param weight
* > [T] The weight as above.
* @param i
* > [int] The time index at which \f$F(t)\f$ and \f$C(t,t^\prime)\f$ are multiplied.
*/  

  template <typename T>
  void herm_matrix_moving<T>::right_multiply(function_moving<T> &ft,T weight,int i){
    assert(size1_==ft.size1());
    assert(size2_==ft.size2());
    assert(size1_==size2_);
    assert(tc_==ft.tc());
    
    cplx *xtemp;
    xtemp=new cplx [element_size_];
    for(int j=0;j<=tc_;j++){
      element_mult<T,LARGESIZE>(size1_,xtemp,retptr(i,j),ft.ptr(j));
      element_smul<T,LARGESIZE>(size1_,xtemp,weight);
      element_set<T,LARGESIZE>(size1_,retptr(i,j),xtemp);
      element_mult<T,LARGESIZE>(size1_,xtemp,lesptr(i,j),ft.ptr(j));
      element_smul<T,LARGESIZE>(size1_,xtemp,weight);
      element_set<T,LARGESIZE>(size1_,lesptr(i,j),xtemp);
    }
    delete [] xtemp;
  }
/** \brief <b> Multiplies all component of the `herm_matrix_moving` with a real scalar. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * Multiplies all component of the `herm_matrix_moving` at a given time step
 * `tstp` with a real scalar `weight`.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > [int] The time step at which the contour object is weighted with a factor.
 * @param weight
 * > [T] The weight.
 */
  template <typename T>
  void herm_matrix_moving<T>::smul(int tstp, T weight)
  {
    int ndata1 = element_size_ * (tc_ + 1);
    for(int i = 0; i < ndata1; ++i)
      {
	les_[tstp][i] *= weight;
	ret_[tstp][i] *= weight;
      }
  }
/** \brief <b> Multiplies all component of the `herm_matrix_moving` with a complex scalar. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * Multiplies all component of the `herm_matrix_moving` at a given time step
 * `tstp` with a real scalar `weight`.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param tstp
 * > [int] The time step at which the contour object is weighted with a factor.
 * @param weight
 * > [T] The weight.
 */

  template <typename T>
  void herm_matrix_moving<T>::smul(int tstp, cplx weight)
  {
    int ndata1 = element_size_ * (tc_ + 1);
    for(int i = 0; i < ndata1; ++i)
      {
	les_[tstp][i] *= weight;
	ret_[tstp][i] *= weight;
      }
  }

#if CNTR_USE_HDF5 == 1
/** \brief <b> Stores `herm_matrix_moving` to a given group in HDF5 format. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Stores the `herm_matrix_moving`, including the cutoff time `tc`, 
 * > the current timestep `t0`, matrix size, elementsize and fermionic/bosonic character,
 * > to a given HDF5 group.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param group_id
 * > The HDF5 group handle under which the `herm_matrix_moving` is stored.
 */
  template <typename T>
  void herm_matrix_moving<T>::write_to_hdf5(hid_t group_id) {
    store_int_attribute_to_hid(group_id, std::string("tc"), tc_);
    store_int_attribute_to_hid(group_id, std::string("t0"), t0_);
    store_int_attribute_to_hid(group_id, std::string("size1"), size1_);
    store_int_attribute_to_hid(group_id, std::string("size2"), size2_);
    store_int_attribute_to_hid(group_id, std::string("element_size"),element_size_);
    store_int_attribute_to_hid(group_id, std::string("sig"),sig_);
    hsize_t len_shape = 3, shape[3],slice[3],start[3];
    shape[1] = size1_;
    shape[2] = size2_;
    slice[1] = size1_;
    slice[2] = size2_;
    start[1] = 0;
    start[2] = 0;
    if (tc_==-1) {
      // continue;
    }else{
      shape[0] = (tc_+1)*(tc_+1)*element_size_;
      slice[0] = (tc_+1)*element_size_;
      //store_cplx_array just to init datastructure
      store_cplx_array_to_hid(group_id, std::string("ret"), retptr(0, 0), shape, len_shape);
      store_cplx_array_to_hid(group_id, std::string("les"), lesptr(0, 0), shape, len_shape);
      for(int t=0;t<=tc_;t++){
	//actual storage routine, check interface of function for simplifications!!
	start[0] = t*(tc_+1)*element_size_;
	store_cplx_slice_to_hid(group_id, std::string("ret"), retptr(t, 0), start, slice, shape, len_shape);
	store_cplx_slice_to_hid(group_id, std::string("les"), lesptr(t, 0), start, slice, shape, len_shape);
	//std::cout<<"lesptr("<<t<<",0)"<<*lesptr(t,0)<<std::endl;
      }
    }
  }

/** \brief <b> Stores `herm_matrix_moving` to a given group in HDF5 format. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Stores the `herm_matrix_moving`, including the cutoff time `tc`, 
 * > the current timestep `t0`, matrix size, elementsize and fermionic/bosonic character,
 * > to a given HDF5 group with given groupname.

 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param group_id
 * > [hid_t] The HDF5 group handle under which the `herm_matrix_moving` is stored.
 * @param groupname
 * > [char*] The name of the HDF5 group.
 */
  template <typename T>
  void herm_matrix_moving<T>::write_to_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = create_group(group_id, groupname);
    this->write_to_hdf5(sub_group_id);
    close_group(sub_group_id);
  }
/** \brief <b> Stores `herm_matrix_moving` to a given file in HDF5 format. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Stores the `herm_matrix_moving`, including the cutoff time `tc`, 
 * > the current timestep `t0`, matrix size, elementsize and fermionic/bosonic character,
 * > to a given HDF5 group with given groupname.

 * <!-- ARGUMENTS
 *      ========= -->
 * @param filename
 * > [char*] The name of the file to which the `herm_matrix_moving` is stored.
 * @param groupname
 * > [char*] The name of the HDF5 group.
 */

  template <typename T>
  void herm_matrix_moving<T>::write_to_hdf5(const char *filename,
					    const char *groupname) {
    hid_t file_id = open_hdf5_file(filename);
    this->write_to_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
  }

/** \brief <b> Reads `herm_matrix_moving` from a given HDF5 group handle. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Reads the `herm_matrix_moving`, including the cutoff time `tc`, 
 * > the current timestep `t0`, matrix size, elementsize and fermionic/bosonic character,
 * > from a given HDF5 group handle.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param group_id
 * > [hid_t] The HDF5 group handle from which the `herm_matrix_moving` is read.
 */
  template <typename T>
  void herm_matrix_moving<T>::read_from_hdf5(hid_t group_id) {
    // -- Read dimensions
    int tc_ = read_primitive_type<int>(group_id, "tc");
    int t0 = read_primitive_type<int>(group_id, "t0");
    int size1_ = read_primitive_type<int>(group_id, "size1");
    int size2_ = read_primitive_type<int>(group_id, "size2");
    int element_size_ = read_primitive_type<int>(group_id, "element_size");
    int sig = read_primitive_type<int>(group_id, "sig");
    // RESIZE G
    this->resize(tc_,size1_);
    sig_ = sig;
    t0_=t0;
    if (tc_!=-1) {
      hsize_t ret_size = ((tc_ + 1) * (tc_ + 1))  * element_size_;
      hsize_t les_size = ((tc_ + 1) * (tc_ + 1))  * element_size_;
      read_primitive_type_array(group_id, "ret", ret_size, retptr(0, 0));
      read_primitive_type_array(group_id, "les", les_size, lesptr(0, 0));
    }
  }
  /** \brief <b> Reads `herm_matrix_moving` from a given HDF5 group handle. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > Reads the `herm_matrix_moving`, including the cutoff time `tc`, 
 * > the current timestep `t0`, matrix size, elementsize and fermionic/bosonic character,
 * > from a given HDF5 group handle.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param group_id
 * > The HDF5 group handle from which the `herm_matrix_moving` is read.
 * @param groupname
 * > The group name from which the `herm_matrix_moving` is read,
 */
  template <typename T>
  void herm_matrix_moving<T>::read_from_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(sub_group_id);
    close_group(sub_group_id);
  }

#endif

}

#endif  // CNTR_HERM_MATRIX_IMPL_H
