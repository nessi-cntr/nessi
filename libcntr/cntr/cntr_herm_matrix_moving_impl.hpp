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
  template <typename T>
  herm_matrix_moving<T>::herm_matrix_moving(){
    data_=0;
    les_=0;
    ret_=0;
    tc_=-1;
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
  
/**
 * @brief Constructor for herm_matrix_moving.
 *
 * This constructor creates a new instance of `herm_matrix_moving`
 * it allocates new data, which are initialized with zero
 *
 * @tparam T Must be either `double` or `float`.
 * @param tc memory truncation
 * @param size1 matrix dimension
 * @param sig statistical sign (bose/fermi)
 */
  template <typename T>
  herm_matrix_moving<T>::herm_matrix_moving(int tc,int size1,int sig){
    assert(-1<=tc);
    assert(1==sig*sig);
    tc_=tc;
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
  
/**
 * @brief constructor for herm_matrix_moving.
 *
 * This constructor creates a new instance of `herm_matrix_moving`
 * it allocates new data, and copies  the internal data pointers and attributes
 * from another instance.
 *
 * @tparam T Must be either `double` or `float`.
 * @param g The `herm_matrix_moving` object to copy from.
 */
  template <typename T>
  herm_matrix_moving<T>::herm_matrix_moving(const herm_matrix_moving &g){
    tc_=g.tc_;
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

/**
 * @brief Assignment operator for herm_matrix_moving.
 *
 * creates a new instance (allocates memory) and copies the data
 */

  template <typename T>
  herm_matrix_moving<T> &  herm_matrix_moving<T>::operator=(const  herm_matrix_moving &g){
    if(this==&g) return *this;
    sig_=g.sig_;
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
  
/**
 * @brief Retrieves the lesser Green's function component for a given time slice and relative time step.
 *
 * This function extracts the lesser component at a specific time slice `delti` and  difference `deltj`
 * and stores it in the provided matrix `M`. The function resizes `M` to match
 * the matrix dimensions (`size1_ × size2_`) before copying the values:
 * `M` set to lesser component at times (`t0-delti,t0-delti-deltj`)
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the output matrix, which must support `resize()` and element access via `operator()`.
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Output matrix where the extracted values will be stored.
 *
 * @pre `delti`, `deltj` must be within the allowed range (`0 <= delti,deltj <= tc_`).
 * @post `M` is resized to `size1_ × size2_` and filled with the extracted values.
 */  
  template<typename T> template <class Matrix>
  void herm_matrix_moving<T>::get_les(int delti,int deltj,Matrix &M){
    assert(0<=delti && delti<=tc_);
    assert(0<=deltj && deltj<=tc_);
    cplx *x;
    x=lesptr(delti,deltj);
    herm_matrix_moving_READ_ELEMENT
      }

/**
 * @brief Retrieves the retarded Green's function component for a given time slice and relative time step.
 *
 * This function extracts the retarded component at a specific time slice `delti` and  difference `deltj`
 * and stores it in the provided matrix `M`. The function resizes `M` to match
 * the matrix dimensions (`size1_ × size2_`) before copying the values:
 * `M` set to retarded component at times (`t0-delti,t0-delti-deltj`)
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the output matrix, which must support `resize()` and element access via `operator()`.
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Output matrix where the extracted values will be stored.
 *
 * @pre `delti`, `deltj` must be within the allowed range (`0 <= delti,deltj <= tc_`).
 * @post `M` is resized to `size1_ × size2_` and filled with the extracted values.
 */
  template<typename T> template <class Matrix>
  void herm_matrix_moving<T>::get_ret(int delti,int deltj,Matrix &M){
    assert(0<=delti && delti<=tc_);
    assert(0<=deltj && deltj<=tc_);
    cplx *x;
    x=retptr(delti,deltj);
    herm_matrix_moving_READ_ELEMENT
      }

/**
 * @brief Retrieves the greater Green's function component for a given time slice and relative time step.
 *
 * This function extracts the greater component at a specific time slice `delti` and  difference `deltj`
 * and stores it in the provided matrix `M`. The function resizes `M` to match
 * the matrix dimensions (`size1_ × size2_`) before copying the values:
 * `M` set to greater component at times (`t0-delti,t0-delti-deltj`)
 * We use the lelation  `greater = retarded + lesser`.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the output matrix, which must support `resize()` and element access via `operator()`.
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Output matrix where the extracted values will be stored.
 *
 * @pre `delti`, `deltj` must be within the allowed range (`0 <= delti,deltj <= tc_`).
 * @post `M` is resized to `size1_ × size2_` and filled with the extracted values.
 */
  template<typename T> template <class Matrix>
  void herm_matrix_moving<T>::get_gtr(int delti,int deltj,Matrix &M){
    Matrix M1;
    get_ret(delti,deltj,M);
    get_les(delti,deltj,M1);
    M += M1;
  }

/**
 * @brief Retrieves the retarded Green's function component for a given time slice and relative time step.
 *
 * This function extracts the greater component at a specific time slice `delti` and  difference `deltj`
 * and stores it in the provided variable `x`. The value of `x` corresponds to the
 * retarded component at the time (`t0-delti,t0-delti-deltj`)
 * This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Output variable where the extracted value will be stored.
 *
 * @pre `delti`, `deltj` must be within the allowed range (`0 <= delti,deltj <= tc_`).
 * @post `x` is updated with the extracted value.
 */
  template<typename T>
  inline void herm_matrix_moving<T>::get_ret(int delti,int deltj,cplx &x) { 
    assert(0<=delti && delti<=tc_);
    assert(0<=deltj && deltj<=tc_);
    x=*retptr(delti,deltj);
  }
/**
 * @brief Retrieves the lesser Green's function component for a given time slice and relative time step.
 *
 * This function extracts the lesser component at a specific time slice `delti` and  difference `deltj`
 * and stores it in the provided variable `x`. The value of `x` corresponds to the
 * lesser component at the time (`t0-delti,t0-delti-deltj`)
 * This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Output variable where the extracted value will be stored.
 *
 * @pre `delti`, `deltj` must be within the allowed range (`0 <= delti,deltj <= tc_`).
 * @post `x` is updated with the extracted value.
 */
  template<typename T>
  inline void herm_matrix_moving<T>::get_les(int delti,int deltj,cplx &x){ 
    assert(0<=delti && delti<=tc_);
    assert(0<=deltj && deltj<=tc_);
    x=*lesptr(delti,deltj);
  }
/**
 * @brief Retrieves the greater Green's function component for a given time slice and relative time step.
 *
 * This function extracts the greater component at a specific time slice `delti` and  difference `deltj`
 * and stores it in the provided variable `x`. The value of `x` corresponds to the
 * greater component at the time (`t0-delti,t0-delti-deltj`)
 * We use the lelation  `greater = retarded + lesser`.
 * This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Output variable where the extracted value will be stored.
 *
 * @pre `delti`, `deltj` must be within the allowed range (`0 <= delti,deltj <= tc_`).
 * @post `x` is updated with the extracted value.
 */
  template<typename T>
  inline void herm_matrix_moving<T>::get_gtr(int delti,int deltj,cplx &x){
    cplx x1;
    get_ret(delti,deltj,x);
    get_les(delti,deltj,x1);
    x+=x1;
  }
/**
 * @brief Computes the density matrix from the lesser Green's function component.
 *
 * This function retrieves the lesser Green's function component at the time step specified by
 * `t0` and `delti` (i.e., with `deltj = 0`),
 * The density matrix is computed as:
 * `density_matrix = imag_i * (Bose/Fermi) * lesser_component`.
 * This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @return The computed density matrix as a complex number.
 *
 * @post The density matrix is calculated based on the retrieved lesser Green's function component.
 */
  template<typename T>
  std::complex<T> herm_matrix_moving<T>::density_matrix(int delti){
    cplx x1;
    get_les(delti,0,x1);
    return std::complex<T>(0.0,sig_)*x1;
  }
/**
 * @brief Computes the density matrix from the lesser Green's function component.
 *
 * This function retrieves the lesser Green's function component at the time step specified by
 * `t0` and `delti` (i.e., with `deltj = 0`), and computes the density matrix by multiplying it
 * with `imag_i * (Bose/Fermi)`. The result is stored in the provided matrix `M`.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the output matrix, which must support element-wise multiplication.
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param M Output matrix where the computed density matrix will be stored.
 *
 * @post The matrix `M` is updated with the computed density matrix.
 */
  template<typename T>
  template<class Matrix> void herm_matrix_moving<T>::density_matrix(int delti,Matrix &M){
    get_les(delti,0,M);
    M *= std::complex<T>(0.0,1.0*sig_);
  }
  
#define herm_matrix_SET_ELEMENT_MATRIX					\
  {									\
    int r, s;								\
    for (r = 0; r < size1_; r++)					\
      for (s = 0; s < size2_; s++)					\
	x[r * size2_ + s] = M(r, s);                                    \
  }
/**
 * @brief Sets the lesser Green's function component for a given time slice and relative time step.
 *
 * This function assigns values from the provided matrix `M` to the lesser component
 * at a specific time slice `delti` and  difference `deltj`. The matrix elements are copied into the
 * corresponding positions of the internal storage: lesser component at `(t0-delti,t0-delti-deltj)` set to `M`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the input matrix, which must support element access via `operator()`.
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Input matrix containing the new values.
 *
 * @pre `delti`,`deltj` must be within the allowed range (`0 <= delti,deltj <= tc_`).
 * @pre The dimensions of `M` must match `size1_ × size2_`.
 */
  template<typename T>
  template<class Matrix>
  void herm_matrix_moving<T>::set_les(int delti,int deltj,Matrix &M){
    assert(delti<=tc_ && deltj<=tc_);
    cplx *x = lesptr(delti,deltj);
    herm_matrix_SET_ELEMENT_MATRIX
      }
  
/**
 * @brief Sets the retarded Green's function component for a given time slice and relative time step.
 *
 * This function assigns values from the provided matrix `M` to the retarded component
 * at a specific time slice `delti` and  difference `deltj`. The matrix elements are copied into the
 * corresponding positions of the internal storage: retarded component at `(t0-delti,t0-delti-deltj)` set to `M`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the input matrix, which must support element access via `operator()`.
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Input matrix containing the new values.
 *
 * @pre `delti`,`deltj` must be within the allowed range (`0 <= delti,deltj <= tc_`).
 * @pre The dimensions of `M` must match `size1_ × size2_`.
 */
  template<typename T>
  template<class Matrix>
  void herm_matrix_moving<T>::set_ret(int delti,int deltj,Matrix &M){
    assert(delti<=tc_ && deltj<=tc_);
    cplx *x = retptr(delti,deltj);
    herm_matrix_SET_ELEMENT_MATRIX
      }

/**
 * @brief Sets the lesser Green's function component for a given time slice and relative time step.
 *
 * This function assigns the scalar value `x` to the lesser component at a specific time slice `delti` and  difference `deltj`. Specifically, it sets the lesser component at `(t0-delti,t0-delti-deltj)`
 * to `x`. This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Scalar value to assign.
 *
 * @pre `delti`,`deltj` must be within the allowed range (`0 <= delti,deltj <= tc_`).
 */
  template<typename T>
  inline void herm_matrix_moving<T>::set_les(int delti,int deltj,cplx x){
    assert(delti<=tc_ && deltj<=tc_);
    *lesptr(delti,deltj) = x;
  }

/**
 * @brief Sets the retarded Green's function component for a given time slice and relative time step.
 *
 * This function assigns the scalar value `x` to the retarded component at a specific time slice `delti` and  difference `deltj`. Specifically, it sets the retarded component at `(t0-delti,t0-delti-deltj)`
 * to `x`. This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Scalar value to assign.
 *
 * @pre `delti`,`deltj` must be within the allowed range (`0 <= delti,deltj <= tc_`).
 */
  template<typename T>
  inline void herm_matrix_moving<T>::set_ret(int delti,int deltj,cplx x){
    assert(delti<=tc_ && deltj<=tc_);
    *retptr(delti,deltj) = x;  
  }
  
  template<typename T> 
  void herm_matrix_moving<T>::print_to_file(const char *file,int precision){
    int i,j,l,sg=element_size_;
    std::ofstream out;
    out.open(file,std::ios::out);
    out.precision(precision);
    out << "# " << " " << tc_ << " " << size1_ << " " << " " << sig_ << std::endl;
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
  template<typename T> 
  void herm_matrix_moving<T>::read_from_file(const char *file){
    int i,tc,j,l,size1,sg,sig;
    double real, imag;
    std::string s;
    std::ifstream out;
    out.open(file,std::ios::in);
    if(!(out >> s >> tc >> size1 >> sig)){
      std::cerr << "read G from file " << file << " error in file" << std::endl; 
      abort();
    }
    if(tc!=tc_ || size1!=size1_) resize(tc,size1);

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
	      std::cerr << "read G from file " << file << " error at les (" << i<< "," << j << ")"<< std::endl; 
	      abort();
	    }
	    lesptr(i,j)[l] = std::complex<T>(real, imag);
	  }
	}
      }
    }
    out.close();
  }

/**
 * @brief Move truncated time window forward.
 *
 * Move truncated time window forward by one timestep.
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
      
    }    
  }

/**
 * @brief Clears the timestep data by setting all elements to zero.
 *
 * This function sets all elements of the retarded and lesser Green’s function
 * components to zero for the entire range of the timestep `delti`.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 *
 */
  template <typename T>
  void herm_matrix_moving<T>::clear_timestep(int delti){
    for(int t1=0;t1<=tc_;t1++){
      element_set_zero<T,LARGESIZE>(size1_,retptr(delti,t1));
      element_set_zero<T,LARGESIZE>(size1_,lesptr(delti,t1));
    }
  }
/**
 * @brief Clears the timestep data by setting all elements to zero.
 *
 * This function sets all elements of the retarded and lesser Green’s function
 * components to zero for the entire range of the timestep `delti`.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 *
 */
  template <typename T>
  void herm_matrix_moving<T>::set_timestep_zero(int delti){
    for(int t1=0;t1<=tc_;t1++){
      element_set_zero<T,LARGESIZE>(size1_,retptr(delti,t1));
      element_set_zero<T,LARGESIZE>(size1_,lesptr(delti,t1));
    }
  }

/**
 * @brief Copies the data from a `herm_matrix_timestep_view` object.
 *
 * This function copies the retarded and lesser Green’s function components
 * from the given `herm_matrix_timestep_view` object into the current instance at a given time step `delti`.
 * \f$ this \f$(`t0-delti,t0-delti-j`) set to \f$ G \f$(`tstp_G,tstp_G-j`) for j on the timestep.
 * `tstp_G` is the time step of the `herm_matrix_timestep_view` object, not a parameter.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param g Source `herm_matrix_timestep_view` object.
 * @param gcc Hermitian conjugate of g.
 *
 * @pre The time cut-off (`tc_`) of both objects must be the same.
 * @pre The matrix sizes (`size1_`) must match.
 * @pre `delti` must be within the allowed range (`0 <= delti <= tc_`).
 *
 * @note The function uses `memcpy` for efficient memory copying.
 */
  template <typename T>
  void herm_matrix_moving<T>::set_timestep(int delti,herm_matrix_timestep_view<T> &g,herm_matrix_timestep_view<T> &gcc){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.set_timestep(g,gcc);
  }

/**
 * @brief Copies the data from a `herm_matrix_timestep` object.
 *
 * This function copies the retarded and lesser Green’s function components
 * from the given `herm_matrix_timestep` object into the current instance at a given time step `delti`.
 * \f$ this \f$(`t0-delti,t0-delti-j`) set to \f$ G \f$(`tstp_G,tstp_G-j`) for j on the timestep.
 * `tstp_G` is the time step of the `herm_matrix_timestep` object, not a parameter.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param g Source `herm_matrix_timestep` object.
 * @param gcc Hermitian conjugate of g.
 *
 * @pre The time cut-off (`tc_`) of both objects must be the same.
 * @pre The matrix sizes (`size1_`) must match.
 * @pre `delti` must be within the allowed range (`0 <= delti <= tc_`).
 *
 * @note The function uses `memcpy` for efficient memory copying.
 */
  template <typename T>
  void herm_matrix_moving<T>::set_timestep(int delti,herm_matrix_timestep<T> &g,herm_matrix_timestep<T> &gcc){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.set_timestep(g,gcc);
  }

/**
 * @brief Copies the data from a `herm_matrix` object.
 *
 * This function copies the retarded and lesser Green’s function components
 * from the given `herm_matrix` object at time step `tstp_G` into the current instance at a given time difference `delti`:
 * \f$ this \f$(`t0-delti,t0-delti-j`) set to \f$ G \f$(`tstp_G,tstp_G-j`) for j on the timestep
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param g Source `herm_matrix` object.
 * @param gcc Hermitian conjugate of g.
 * @param tstp_G Time step of g.
 *
 * @pre The time cut-off (`tc_`) of both objects must be the same.
 * @pre The matrix sizes (`size1_`) must match.
 * @pre `delti` must be within the allowed range (`0 <= delti <= tc_`).
 *
 * @note The function uses `memcpy` for efficient memory copying.
 */
  template <typename T>
  void herm_matrix_moving<T>::set_timestep(int delti,herm_matrix<T> &g,herm_matrix<T> &gcc,int tstp_G){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.set_timestep(g,gcc,tstp_G);
  }

/**
 * @brief Initialize a moving window to start  a truncated time evolution.
 *
 * Equivalent to `C.set_timestep(i,g,gcc,tstp-i)` for `i=0,...,tc`;
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g Source `herm_matrix` object.
 * @param gcc Hermitian conjugate of g.
 * @param tstp Time step of g.
 *
 * @pre The matrix sizes (`size1_`) must match.
 * @pre `tc <= tstp <= g.nt()`
 */
  template <typename T>
  void herm_matrix_moving<T>::set_from_G_backward(herm_matrix<T> &g,herm_matrix<T> &gcc,int tstp){
    assert(tc_<=tstp && tstp <= g.nt());
    for(int i=0;i<=tc_;i++) set_timestep(i,g,gcc,tstp-i);
  }

/**
 * @brief Initialize a moving window to start  a truncated time evolution assuming Hermiticity.
 *
 * Equivalent to `C.set_timestep(i,g,g,tstp-i)` for `i=0,...,tc`;
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g Hermitian source `herm_matrix` object.
 * @param tstp Time step of g.
 *
 * @pre The matrix sizes (`size1_`) must match.
 * @pre `tc <= tstp <= g.nt()`
 */
  template <typename T>
  void herm_matrix_moving<T>::set_from_G_backward(herm_matrix<T> &g,int tstp){
    set_from_G_backward(g,g,tstp);
  }

  template <typename T>
  void herm_matrix_moving<T>::set_from_G_backward(const herm_pseudo<T> &g,int tstp){
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
  template <typename T>
  void herm_matrix_moving<T>::set_from_G_backward(int tstp, const herm_pseudo<T> &g){
    set_from_G_backward(g,tstp);
  }
  
/**
 * @brief Sets a specific matrix element from another `herm_matrix_moving` object at a given relative timestep `delti_g`.
 *
 * This function copies the matrix element located at `(j1, j2)` in the source object `g`
 * to the position `(i1, i2)` in the given time step `delti` of the current object for all times on the timestep `delti_g` of `g`:
 * \f$ this_{i1,i2}\f$ (`t0-delti,t0-delti-j`)set to \f$ G_{j1,j2}\f$ (`t0-delti_g,t0-delti_g-j`) for all j
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step of this object (must be within the time cutoff `tc_`).
 * @param i1 Row index of the destination matrix element.
 * @param i2 Column index of the destination matrix element.
 * @param g Source `herm_matrix_timestep_moving_view` object.
 * @param delti_g The relative time step at which to create the view. Must be in the range `[0, g.tc()]`.
 * @param j1 Row index of the source matrix element.
 * @param j2 Column index of the source matrix element.
 *
 * @note The function asserts that the time cutoffs (`tc_`) of both objects match.
 *       It also checks that indices are within valid bounds.
 */
  template <typename T>
  void herm_matrix_moving<T>::set_matrixelement(int delti, int i1, int i2,herm_matrix_moving<T> &g,int delti_g, int j1, int j2) {
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.set_matrixelement(i1,i2,g,delti_g,j1,j2);
  }
  
/**
 * @brief Sets a specific matrix element from a `herm_matrix_timestep_moving` object.
 *
 * This function copies the matrix element located at `(j1, j2)` in the source object `g`
 * to the position `(i1, i2)` in the current object  on the timestep `delti`:
 * \f$ this_{i1,i2} \f$ (`t0-delti,t0-delti-j`) set to \f$ G_{j1,j2} \f$ (`t0,t0-j`)for all `j`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param i1 Row index of the destination matrix element.
 * @param i2 Column index of the destination matrix element.
 * @param g Source `herm_matrix_timestep_moving` object.
 * @param j1 Row index of the source matrix element.
 * @param j2 Column index of the source matrix element.
 *
 * @note The function asserts that the time cutoffs (`tc_`) of both objects match.
 *       It also checks that indices are within valid bounds.
 */
  template <typename T>
  void herm_matrix_moving<T>::set_matrixelement(int delti, int i1, int i2,herm_matrix_timestep_moving<T> &g, int j1,int j2) {
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.set_matrixelement(i1,i2,g,j1,j2);
  }

/**
 * @brief Sets a specific matrix element from a `herm_matrix_timestep_moving_view` object.
 *
 * This function copies the matrix element located at `(j1, j2)` in the source object `g`
 * to the position `(i1, i2)` in the current object  on the timestep `delti`:
 * \f$ this_{i1,i2} \f$ (`t0-delti,t0-delti-j`) set to \f$ G_{j1,j2} \f$ (`t0,t0-j`)for all `j`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param i1 Row index of the destination matrix element.
 * @param i2 Column index of the destination matrix element.
 * @param g Source `herm_matrix_timestep_moving_view` object.
 * @param j1 Row index of the source matrix element.
 * @param j2 Column index of the source matrix element.
 *
 * @note The function asserts that the time cutoffs (`tc_`) of both objects match.
 *       It also checks that indices are within valid bounds.
 */
  template <typename T>
  void herm_matrix_moving<T>::set_matrixelement(int delti, int i1, int i2,herm_matrix_timestep_moving_view<T> &g, int j1,int j2) {
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.set_matrixelement(i1,i2,g,j1,j2);
  }
  
  
/**
 * @brief Copies the data from a `herm_matrix_timestep_moving_view` object.
 *
 * This function copies the retarded and lesser Green’s function components
 * from the given `herm_matrix_timestep_moving_view` object into the current instance at a given time step `delti`.
 * \f$ this \f$(`t0-delti,t0-delti-j`) set to \f$ G \f$(`t0,t0-j`) for j on the timestep
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param g Source `herm_matrix_timestep_moving_view` object.
 * @param gcc Hermitian conjugate of g.
 *
 * @pre The time cut-off (`tc_`) of both objects must be the same.
 * @pre The matrix sizes (`size1_`) must match.
 * @pre `delti` must be within the allowed range (`0 <= delti <= tc_`).
 *
 * @note The function uses `memcpy` for efficient memory copying.
 */
  template <typename T>
  void herm_matrix_moving<T>::set_timestep(int delti,herm_matrix_timestep_moving_view<T> &g){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.set_timestep(g);
  }

/**
 * @brief Copies the data from a `herm_matrix_timestep_moving` object.
 *
 * This function copies the retarded and lesser Green’s function components
 * from the given `herm_matrix_timestep_moving` object into the current instance at a given time step `delti`.
 * \f$ this \f$(`t0-delti,t0-delti-j`) set to \f$ G \f$(`t0,t0-j`) for j on the timestep
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param g Source `herm_matrix_timestep_moving` object.
 * @param gcc Hermitian conjugate of g.
 *
 * @pre The time cut-off (`tc_`) of both objects must be the same.
 * @pre The matrix sizes (`size1_`) must match.
 * @pre `delti` must be within the allowed range (`0 <= delti <= tc_`).
 *
 * @note The function uses `memcpy` for efficient memory copying.
 */
  template <typename T>
  void herm_matrix_moving<T>::set_timestep(int delti,herm_matrix_timestep_moving<T> &g){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.set_timestep(g);
  }

/**
 * @brief Copies the data from a `herm_matrix_moving` object.
 *
 * This function copies the retarded and lesser Green’s function components
 * from the given `herm_matrix_moving` object at time step `delti_g` into the current instance at a given time step `delti`:
 * \f$ this \f$(`t0-delti,t0-delti-j`) set to \f$ G \f$(`t0-delti_g,t0-delti_g-j`) for j on the timestep
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time step (must be within the time cutoff `tc_`).
 * @param g Source `herm_matrix_moving` object.
 * @param gcc Hermitian conjugate of g.
 * @param delti_g Time step of g (must be within the time cutoff `tc_`).
 *
 * @pre The time cut-off (`tc_`) of both objects must be the same.
 * @pre The matrix sizes (`size1_`) must match.
 * @pre `delti_g`,`delti` must be within the allowed range (`0 <= delti_g,delti <= tc_`).
 *
 * @note The function uses `memcpy` for efficient memory copying.
 */
  template <typename T>
  void herm_matrix_moving<T>::set_timestep(int delti,herm_matrix_moving<T> &g,int delti_g){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.set_timestep(g,delti_g);
  }
  
  
/**
 * @brief Increments the Green's function at a specified timestep by another Green's function at the given timestep, scaled by a factor.
 *
 * This function performs the operation:
 * \f$ this\f$ (`t0-delti,t0-delti-j`) -->   \f$ this\f$(`t0-delti,t0-delti-j`) + `\alpha` * \f$ G\f$(`t0,t0-j`) for all j on the timestep
 * where g is another Green’s function of the same dimensions,
 * and `\alpha` is a complex scaling factor. The increment is applied to both the
 * retarded and lesser Green’s function components.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time slice of the current object to add to.
 * @param g The Green’s function G(t,t') at given time slice to be added.
 * @param alpha the complex scaling factor
 *
 * @pre `g` must have the same `tc_`, `size1_`, and `size2_` as this object.
 * @post The retarded and lesser components are incremented by `\alpha * G(t,t')`.
 */
  template <typename T>
  void herm_matrix_moving<T>::incr_timestep(int delti,herm_matrix_timestep_moving_view<T> &g,cplx alpha){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.incr_timestep(g,alpha);
  }
  
/**
 * @brief Increments the Green's function at a specified timestep by another Green's function at the given timestep, scaled by a factor.
 *
 * This function performs the operation:
 * \f$ this\f$ (`t0-delti,t0-delti-j`) -->   \f$ this\f$(`t0-delti,t0-delti-j`) + `\alpha` * \f$ G\f$(`t0,t0-j`) for all j on the timestep
 * where g is another Green’s function of the same dimensions,
 * and `\alpha` is a complex scaling factor. The increment is applied to both the
 * retarded and lesser Green’s function components.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time slice of the current object to add to.
 * @param g The Green’s function G(t,t') at given time slice to be added.
 * @param alpha the complex scaling factor
 *
 * @pre `g` must have the same `tc_`, `size1_`, and `size2_` as this object.
 * @post The retarded and lesser components are incremented by `\alpha * G(t,t')`.
 */
  template <typename T>
  void herm_matrix_moving<T>::incr_timestep(int delti,herm_matrix_timestep_moving<T> &g,cplx alpha){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.incr_timestep(g,alpha);
  }

/**
 * @brief Increments the Green's function at a specified timestep by another Green's function at a specified timestep, scaled by a factor.
 *
 * This function performs the operation:
 * \f$ this\f$ (`t0-delti,t0-delti-j`) -->   \f$ this\f$(`t0-delti,t0-delti-j`) + `\alpha` * \f$ G\f$(`t0-delti_g,t0-delti_g-j`) for all j on the timestep
 * where g is another Green’s function of the same dimensions,
 * and `\alpha` is a complex scaling factor. The increment is applied to both the
 * retarded and lesser Green’s function components.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param delti Time slice of the current object to add to.
 * @param g The Green’s function G(t,t') at specified time slice to be added.
 * @param delti_g Time slice of the object to be added..
 * @param alpha the complex scaling factor
 *
 * @pre `g` must have the same `tc_`, `size1_`, and `size2_` as this object.
 * @post The retarded and lesser components are incremented by `\alpha * G(t,t')`.
 */
  template <typename T>
  void herm_matrix_moving<T>::incr_timestep(int delti,herm_matrix_moving<T> &g,int delti_g,cplx alpha){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.incr_timestep(g,delti_g,alpha);
  }
  
/**
 * @brief Left-multiplies the Green's function by a time-dependent contour function.
 *
 * This function performs the operation \f$ G(t,t') \rightarrow F(t) G(t,t') \f$
 * at the leading time step, where `F(t)` is a time-dependent contour function.
 * The multiplication is applied element-wise to both the retarded and lesser
 * Green's function components.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param ft The time-dependent contour function `F(t)`, which must have the same time cutoff (`tc_`).
 * @param alpha  a scalar
 *
 * @pre The time cutoff `tc_` of `ft` must match `tc_` of the Green’s function.
 * @post The retarded and lesser Green’s function components are left-multiplied by `F(t)`.
 */
  template <typename T>
  void herm_matrix_moving<T>::left_multiply(function_moving<T> &ft){
    herm_matrix_timestep_moving_view<T> selfview(*this,0);
    selfview.left_multiply(ft);
  }

/**
 * @brief Left-multiplies the Green's function by the hermitian conjugate of a  time-dependent contour function.
 *
 * This function performs the operation \f$ G(t,t') \rightarrow F(t)^\dagger G(t,t') \f$
 * at the leading time step, where `F(t)` is a time-dependent contour function.
 * The multiplication is applied element-wise to both the retarded and lesser
 * Green's function components.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param ft The time-dependent contour function `F(t)`, which must have the same time cutoff (`tc_`).
 *
 * @pre The time cutoff `tc_` of `ft` must match `tc_` of the Green’s function.
 * @post The retarded and lesser Green’s function components are left-multiplied by `F(t)`.
*/
  template <typename T>
  void herm_matrix_moving<T>::left_multiply_hermconj(function_moving<T> &ft){
    herm_matrix_timestep_moving_view<T> selfview(*this,0);
    selfview.left_multiply_hermconj(ft);
  }
  
/**
 * @brief Right-multiplies the Green's function by a time-dependent contour function.
 *
 * This function performs the operation \f$ G(t,t') \rightarrow G(t,t') F(t')  \f$
 * at the leading time step, where `F(t')` is a time-dependent contour function
 * The multiplication is applied element-wise to both the retarded and
 * lesser Green's function components.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param ft The time-dependent contour function `F(t')`, which must have a time cutoff (`tc_`)
 *           greater than or equal to that of the Green's function.
 *
 * @pre The time cutoff `tc_` of `ft` must be at least as large as `tc_` of the Green’s function.
 * @post The retarded and lesser Green’s function components are right-multiplied by `F(t')`,
 *       and then scaled by `\alpha` if `\alpha ≠ 1.0`.
 */
  template <typename T>
  void herm_matrix_moving<T>::right_multiply(function_moving<T> &ft){
    herm_matrix_timestep_moving_view<T> selfview(*this,0);
    selfview.right_multiply(ft);
  }
  
/**
 * @brief Right-multiplies the Green's function by a time-dependent contour function.
 *
 * This function performs the operation \f$ G(t,t') \rightarrow G(t,t') F(t')^\dagger  \f$
 * at the leading time step, where `F(t')` is a time-dependent contour function. The
 * multiplication is applied element-wise to both the retarded and
 * lesser Green's function components.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param ft The time-dependent contour function `F(t')`, which must have a time cutoff (`tc_`)
 *           greater than or equal to that of the Green's function.
 *
 * @pre The time cutoff `tc_` of `ft` must be at least as large as `tc_` of the Green’s function.
 * @post The retarded and lesser Green’s function components are right-multiplied by `F(t')`,
 *       and then scaled by `\alpha` if `\alpha ≠ 1.0`.
 */
  
  template <typename T>
  void herm_matrix_moving<T>::right_multiply_hermconj(function_moving<T> &ft){
    herm_matrix_timestep_moving_view<T> selfview(*this,0);
    selfview.right_multiply_hermconj(ft);
  }
  
/**
 * @brief Scales the Green's function at a given time step components by a given weight.
 *
 * This function multiplies all elements of the retarded and lesser Green's function
 * components by the provided scalar `weight` at a specified time step `delti`.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param weight The scalar factor by which all elements are multiplied.
 *
 * @post Both the retarded and lesser Green's function components are scaled by `weight`.
 */
  template <typename T>
  void herm_matrix_moving<T>::smul(int delti, T weight){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.smul(weight);
  }

/**
 * @brief Scales the Green's function at a given time step components by a given weight.
 *
 * This function multiplies all elements of the retarded and lesser Green's function
 * components by the provided scalar `weight` at a specified time step `delti`.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param weight The scalar complex type factor by which all elements are multiplied.
 *
 * @post Both the retarded and lesser Green's function components are scaled by `weight`.
 */
  template <typename T>
  void herm_matrix_moving<T>::smul(int delti, cplx weight){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.smul(weight);
  }
  
#if CNTR_USE_HDF5 == 1
  
  template <typename T>
  void herm_matrix_moving<T>::write_to_hdf5(hid_t group_id) {
    store_int_attribute_to_hid(group_id, std::string("herm_matrix_moving"), 1);
    store_int_attribute_to_hid(group_id, std::string("tc"), tc_);
    store_int_attribute_to_hid(group_id, std::string("size1"), size1_);
    store_int_attribute_to_hid(group_id, std::string("size2"), size2_);
    store_int_attribute_to_hid(group_id, std::string("element_size"),element_size_);
    store_int_attribute_to_hid(group_id, std::string("sig"),sig_);
    hsize_t len_shape = 3, shape[3],start[3];
    shape[1] = size1_;
    shape[2] = size2_;
    if (tc_==-1) {
      // continue;
    }else{
      shape[0] = (tc_+1)*(tc_+1);
      //store_cplx_array just to init datastructure
      store_cplx_array_to_hid(group_id, std::string("ret"), retptr(0, 0), shape, len_shape);
      store_cplx_array_to_hid(group_id, std::string("les"), lesptr(0, 0), shape, len_shape);
    }
  }
  template <typename T>
  void herm_matrix_moving<T>::write_to_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = create_group(group_id, groupname);
    this->write_to_hdf5(sub_group_id);
    close_group(sub_group_id);
  }
  template <typename T>
  void herm_matrix_moving<T>::write_to_hdf5(const char *filename,
					    const char *groupname, h5_mode mode) {
    hid_t file_id = open_hdf5_file(filename,mode);
    this->write_to_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
  }
  template <typename T>
  void herm_matrix_moving<T>::write_timestep_to_hdf5(int delti,hid_t group_id){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.write_to_hdf5(group_id);    
  }
  template <typename T>
  void herm_matrix_moving<T>::write_timestep_to_hdf5(int delti,hid_t group_id, const char *groupname){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.write_to_hdf5(group_id,groupname);
  }
  template <typename T>
  void herm_matrix_moving<T>::write_timestep_to_hdf5(int delti,const char *filename, const char *groupname, h5_mode mode){
    herm_matrix_timestep_moving_view<T> selfview(*this,delti);
    selfview.write_to_hdf5(filename,groupname,mode);
  }
  
  template <typename T>
  void herm_matrix_moving<T>::read_from_hdf5(hid_t group_id) {
    // -- Read dimensions
    int tc_ = read_primitive_type<int>(group_id, "tc");
    int size1_ = read_primitive_type<int>(group_id, "size1");
    int size2_ = read_primitive_type<int>(group_id, "size2");
    int element_size_ = read_primitive_type<int>(group_id, "element_size");
    int sig = read_primitive_type<int>(group_id, "sig");
    // RESIZE G
    this->resize(tc_,size1_);
    sig_ = sig;

    if (tc_!=-1) {
      hsize_t ret_size = ((tc_ + 1) * (tc_ + 1))  * element_size_;
      hsize_t les_size = ((tc_ + 1) * (tc_ + 1))  * element_size_;
      read_primitive_type_array(group_id, "ret", ret_size, retptr(0, 0));
      read_primitive_type_array(group_id, "les", les_size, lesptr(0, 0));
    }
  }
  template <typename T>
  void herm_matrix_moving<T>::read_from_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(sub_group_id);
    close_group(sub_group_id);
  }
  template <typename T>
  void herm_matrix_moving<T>::read_from_hdf5(const char *filename,const char *groupname) {
    hid_t file_id = read_hdf5_file(filename);
    this->read_from_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
  }
  
  
#endif
  
}

#endif  // CNTR_HERM_MATRIX_IMPL_H
