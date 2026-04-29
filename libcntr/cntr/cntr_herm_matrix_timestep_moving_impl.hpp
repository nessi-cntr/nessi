#ifndef CNTR_HERM_MATRIX_TIMESTEP_MOVING_IMPL_H
#define CNTR_HERM_MATRIX_TIMESTEP_MOVING_IMPL_H

#include "cntr_herm_matrix_timestep_decl.hpp"
#include "cntr_herm_matrix_timestep_moving_decl.hpp"
#include "cntr_elements.hpp"
#include "cntr_function_decl.hpp"
#include "cntr_herm_matrix_decl.hpp"
#include "cntr_herm_matrix_timestep_view_impl.hpp"

namespace cntr {
  template <typename T> 
  herm_matrix_timestep_moving<T>::herm_matrix_timestep_moving(){
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
  herm_matrix_timestep_moving<T>::~herm_matrix_timestep_moving(){
    if (data_!=0) delete [] data_;
  }

/**
 * @brief Constructor for herm_matrix_timestep_moving.
 *
 * This constructor creates a new instance of `herm_matrix_timestep_moving`
 * it allocates new data, which are initialized with zero
 *
 * @tparam T Must be either `double` or `float`.
 * @param tc memory truncation
 * @param size1 matrix dimension
 * @param sig statistical sign (bose/fermi)
 */
  template <typename T>
  herm_matrix_timestep_moving<T>::herm_matrix_timestep_moving(int tc,int size1,int sig){
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
      long ndata1=(tc_+1)*element_size_;
      data_ = new cplx [2*ndata1];
      ret_ = data_;
      les_ = data_+ndata1;
      memset(data_, 0, 2*sizeof(cplx)*ndata1);
    }
  }

/**
 * @brief constructor for herm_matrix_timestep_moving.
 *
 * This constructor creates a new instance of `herm_matrix_timestep_moving`
 * it allocates new data, and copies  the internal data pointers and attributes 
 * from another instance.
 *
 * @tparam T Must be either `double` or `float`.
 * @param g The `herm_matrix_timestep_moving` object to copy from.
 */ 
  template <typename T>
  herm_matrix_timestep_moving<T>::herm_matrix_timestep_moving(const herm_matrix_timestep_moving &g){
    tc_=g.tc();
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

/**
 * @brief constructor for herm_matrix_timestep_moving.
 *
 * This constructor creates a new instance of `herm_matrix_timestep_moving`
 * it allocates new data, and copies  the internal data pointers and attributes 
 * from timestep `delti_g` of `g`
 *
 * @tparam T Must be either `double` or `float`.
 * @param g The `herm_matrix_timestep_moving` object to copy from.
 * @param `delti_g` the timestep to be copied from
 */
  template <typename T>
  herm_matrix_timestep_moving<T>::herm_matrix_timestep_moving(herm_matrix_moving<T> &g,int delti_g){
    int tc=g.tc();
    assert(0<=delti_g && delti_g<=tc);
    tc_=g.tc();
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
      memcpy(ret_, g.retptr(delti_g,0), sizeof(cplx)*ndata1);
      memcpy(les_, g.lesptr(delti_g,0), sizeof(cplx)*ndata1);
    }
  }
/**
 * @brief Assignment operator for herm_matrix_timestep_moving.
 *
 * creates a new instance (allocates memory) and copies the data
 */

  template <typename T>
  herm_matrix_timestep_moving<T> &  herm_matrix_timestep_moving<T>::operator=(const  herm_matrix_timestep_moving &g){
    if(this==&g) return *this;
    sig_=g.sig_;
    if( tc_!=g.tc_ || size1_!=g.size1_ || size2_!=g.size2_){
      // reallocate
      if (data_!=0) delete [] data_;
      tc_=g.tc_;
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
    clear_timestep();
  }
  
  
  // READING ELEMENTS TO ANY MATRIX TYPE OR TO COMPLEX NUMBERS
  // (then only the (0,0) element is addressed for dim>0)
#define herm_matrix_timestep_moving_READ_ELEMENT {int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++) M(r,s)=x[r*dim+s];}
#define herm_matrix_timestep_moving_READ_ELEMENT_MINUS_CONJ {cplx w;int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++){ w=x[s*dim+r];M(r,s)=std::complex<T>(-w.real(),w.imag());}}
/**
 * @brief Retrieves the lesser Green's function component for a given relative time step.
 *
 * This function extracts the lesser component at a specific time difference `deltj` 
 * and stores it in the provided matrix `M`. The function resizes `M` to match 
 * the matrix dimensions (`size1_ × size2_`) before copying the values:
 * `M` set to lesser component at times (`t0,t0-deltj`)
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the output matrix, which must support `resize()` and element access via `operator()`.
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Output matrix where the extracted values will be stored.
 *
 * @pre `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 * @post `M` is resized to `size1_ × size2_` and filled with the extracted values.
 */  
  template<typename T> template <class Matrix> 
  void herm_matrix_timestep_moving<T>::get_les(int j,Matrix &M){
    assert(0<=j && j<=tc_);
    cplx *x;
    x=lesptr(j);
    herm_matrix_timestep_moving_READ_ELEMENT
  }
/**
 * @brief Retrieves the retarded Green's function component for a given relative time step.
 *
 * This function extracts the retarded component at a specific time difference `deltj` 
 * and stores it in the provided matrix `M`. The function resizes `M` to match 
 * the matrix dimensions (`size1_ × size2_`) before copying the values:
 * `M` set to retarded component at times (`t0,t0-deltj`)
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the output matrix, which must support `resize()` and element access via `operator()`.
 * @param `deltj` Relative time index (must be within the time cutoff `tc_`).
 * @param `M` Output matrix where the extracted values will be stored.
 *
 * @pre `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 * @post `M` is resized to `size1_ × size2_` and filled with the extracted values.
 */
  template<typename T> template <class Matrix> 
  void herm_matrix_timestep_moving<T>::get_ret(int j,Matrix &M){
    assert(0<=j && j<=tc_);
    cplx *x;
    x=retptr(j);
    herm_matrix_timestep_moving_READ_ELEMENT
  }
/**
 * @brief Retrieves the greater Green's function component for a given relative time step.
 *
 * This function extracts the greater component at a specific time difference `deltj` 
 * and stores it in the provided matrix `M`. The function resizes `M` to match 
 * the matrix dimensions (`size1_ × size2_`) before copying the values:
 * `M` set to greater component at times (`t0,t0-deltj`)
 * We use the lelation  `greater = retarded + lesser`.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the output matrix, which must support `resize()` and element access via `operator()`.
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param `M` Output matrix where the extracted values will be stored.
 *
 * @pre `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 * @post `M` is resized to `size1_ × size2_` and filled with the extracted values.
 */  
  template<typename T> template <class Matrix> 
  void herm_matrix_timestep_moving<T>::get_gtr(int j,Matrix &M){
    Matrix M1;
    get_ret(j,M);
    get_les(j,M1);
    M += M1;
  }
/**
 * @brief Retrieves the retarded Green's function component for a given relative time step.
 *
 * This function extracts the retarded component at a specific time difference `deltj` 
 * and stores it in the provided variable `x`. The value of `x` corresponds to the 
 * retarded component at the time `(t0,t0-deltj)`.
 * This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Output variable where the extracted value will be stored.
 *
 * @pre `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 * @post `x` is updated with the extracted value.
 */  
  
  template<typename T> 
  inline void herm_matrix_timestep_moving<T>::get_ret(int j,cplx &x){
    assert(0<=j && j<=tc_);
    x=*retptr(j);
  }
/**
 * @brief Retrieves the lesser Green's function component for a given relative time step.
 *
 * This function extracts the lesser component at a specific time difference `deltj` 
 * and stores it in the provided variable `x`. The value of `x` corresponds to the 
 * lesser component at the time `(t0,t0-deltj)`.
 * This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Output variable where the extracted value will be stored.
 *
 * @pre `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 * @post `x` is updated with the extracted value.
 */

  template<typename T> 
  inline void herm_matrix_timestep_moving<T>::get_les(int j,cplx &x){
    assert(0<=j && j<=tc_);
    x=*lesptr(j);
  }
/**
 * @brief Retrieves the greater Green's function component for a given relative time step.
 *
 * This function extracts the greater component at a specific time difference `deltj` 
 * and stores it in the provided variable `x`. The value of `x` corresponds to the 
 * greater component at the time `(t0,t0-deltj)`.
 * We use the lelation  `greater = retarded + lesser`.
 * This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Output variable where the extracted value will be stored.
 *
 * @pre `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 * @post `x` is updated with the extracted value.
 */    
  template<typename T> 
  inline void herm_matrix_timestep_moving<T>::get_gtr(int j,cplx &x){
    cplx x1;
    get_ret(j,x);
    get_les(j,x1);
    x+=x1;
  }
/**
 * @brief Computes the density matrix from the lesser Green's function component.
 *
 * This function retrieves the lesser Green's function component at the time step 
 * `t0` (i.e., with `deltj = 0`),
 * The density matrix is computed as: 
 * `density_matrix = imag_i * (Bose/Fermi) * lesser_component`.
 * This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @return The computed density matrix as a complex number.
 *
 * @post The density matrix is calculated based on the retrieved lesser Green's function component.
 */  
  template<typename T> 
  std::complex<T> herm_matrix_timestep_moving<T>::density_matrix(void){
    cplx x1;
    get_les(0,x1);
    return std::complex<T>(0.0,sig_)*x1;
  }

/**
 * @brief Computes the density matrix from the lesser Green's function component.
 *
 * This function retrieves the lesser Green's function component at the time step 
 * `t0` (i.e., with `deltj = 0`) and computes the density matrix by multiplying it 
 * with `imag_i * (Bose/Fermi)`. The result is stored in the provided matrix `M`.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the output matrix, which must support element-wise multiplication.
 * @param M Output matrix where the computed density matrix will be stored.
 *
 * @post The matrix `M` is updated with the computed density matrix.
 */  
  template<typename T> template<class Matrix> 
  void herm_matrix_timestep_moving<T>::density_matrix(Matrix &M){
    get_les(0,M);
    M *= std::complex<T>(0.0,1.0*sig_);
  }

#define herm_matrix_SET_ELEMENT_MATRIX		\
  {						\
    int r, s;					\
    for (r = 0; r < size1_; r++)		\
      for (s = 0; s < size2_; s++)		\
	x[r * size2_ + s] = M(r, s);		\
  }
  
/**
 * @brief Sets the lesser Green's function component for a given relative time step.
 *
 * This function assigns values from the provided matrix `M` to the lesser component 
 * at a specific time difference `deltj`. The matrix elements are copied into the 
 * corresponding positions of the internal storage: lesser component at `(t0,t0-deltj)` set to `M
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the input matrix, which must support element access via `operator()`.
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Input matrix containing the new values.
 *
 * @pre `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 * @pre The dimensions of `M` must match `size1_ × size2_`.
 */  
  template<typename T>
  template<class Matrix>
  void herm_matrix_timestep_moving<T>::set_les(int j,Matrix &M){
    assert(j<=tc_);
    cplx *x=lesptr(j);
    herm_matrix_SET_ELEMENT_MATRIX
      }
/**
 * @brief Sets the retarded Green's function component for a given relative time step.
 *
 * This function assigns values from the provided matrix `M` to the retarded component 
 * at a specific time difference `deltj`. The matrix elements are copied into the 
 * corresponding positions of the internal storage: retarded component at `(t0,t0-deltj)` set to `M`
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @tparam Matrix Type of the input matrix, which must support element access via `operator()`.
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param M Input matrix containing the new values.
 *
 * @pre `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 * @pre The dimensions of `M` must match `size1_ × size2_`.
 */
  template<typename T>
  template<class Matrix>
  void herm_matrix_timestep_moving<T>::set_ret(int j,Matrix &M){
    assert(j<=tc_);
    cplx *x=retptr(j);
    herm_matrix_SET_ELEMENT_MATRIX
      }
  
/**
 * @brief Sets the lesser Green's function component for a given relative time step.
 *
 * This function assigns the scalar value `x` to the lesser component at a specific 
 * time difference `deltj`. Specifically, it sets the lesser component at `(t0, t0 - deltj)
 * to `x`. This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Scalar value to assign.
 *
 * @pre `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 */  
  template<typename T>
  inline void herm_matrix_timestep_moving<T>::set_les(int j,cplx x){
    assert(j<=tc_);
    x=*lesptr(j);
  }
/**
 * @brief Sets the retarded Green's function component for a given relative time step.
 *
 * This function assigns the scalar value `x` to the retarded component at a specific 
 * time difference `deltj`. Specifically, it sets the retarded component at (t0, t0 - deltj) 
 * to `x`. This function is intended for cases where the matrix size is 1×1.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param deltj Relative time index (must be within the time cutoff `tc_`).
 * @param x Scalar value to assign.
 *
 * @pre `deltj` must be within the allowed range (`0 <= deltj <= tc_`).
 */  
  template<typename T> 
  inline void herm_matrix_timestep_moving<T>::set_ret(int j,cplx x){
    assert(j<=tc_);
    x=*retptr(j);
  }

/**
 * @brief Clears the timestep data by setting all elements to zero.
 *
 * This function sets all elements of the retarded and lesser Green’s function 
 * components to zero for the entire range of the timestep
 *
 * @tparam T Floating-point type (`double` or `float`).
 *
 */

  template <typename T> 
  void herm_matrix_timestep_moving<T>::clear_timestep(void){
    if(tc_==-1) return;
    memset(data_, 0, 2*sizeof(cplx)*(tc_+1)*element_size_);
  }
/**
 * @brief Clears the timestep data by setting all elements to zero.
 *
 * This function sets all elements of the retarded and lesser Green’s function 
 * components to zero for the entire range of the timestep
 *
 * @tparam T Floating-point type (`double` or `float`).
 *
 */
  template <typename T> 
  void herm_matrix_timestep_moving<T>::set_timestep_zero(void){
    if(tc_==-1) return;
    memset(data_, 0, 2*sizeof(cplx)*(tc_+1)*element_size_);
  }


  
/**
 * @brief Copies the data from another `herm_matrix_timestep_moving_view` object.
 *
 * This function copies the retarded and lesser Green’s function components 
 * from the given `herm_matrix_timestep_moving_view` object into the current instance.
 * \f$ this \f$(`t0,t0-j`) set to \f$ G \f$(`t0,t0-j`) for j on the timestep
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g1 Source `herm_matrix_timestep_moving_view` object.
 *
 * @pre The time cut-off (`tc_`) of both objects must be the same.
 * @pre The matrix sizes (`size1_`) must match.
 *
 * @note The function uses `memcpy` for efficient memory copying.
 */
  template <typename T>
  void herm_matrix_timestep_moving<T>::set_timestep(herm_matrix_timestep_moving_view<T> &g1){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.set_timestep(g1);      
  }
  
/**
 * @brief Copies the data from another `herm_matrix_timestep_moving` object.
 *
 * This function copies the retarded and lesser Green’s function components 
 * from the given `herm_matrix_timestep_moving` object into the current instance.
 * \f$ this \f$(`t0,t0-j`) set to \f$ G \f$(`t0,t0-j`) for j on the timestep
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g1 Source `herm_matrix_timestep_moving` object.
 *
 * @pre The time cut-off (`tc_`) of both objects must be the same.
 * @pre The matrix sizes (`size1_`) must match.
 *
 * @note The function uses `memcpy` for efficient memory copying.
 */
  template <typename T>
  void herm_matrix_timestep_moving<T>::set_timestep(herm_matrix_timestep_moving<T> &g1){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.set_timestep(g1);      
  }
/**
 * @brief Copies the data from another `herm_matrix_moving` object at its relative timestep `delti_g`
 *
 * This function copies the retarded and lesser Green’s function components 
 * from the given `herm_matrix_moving` object into the current instance:
 * \f$ this \f$(`t0,t0-j`) set to \f$ G \f$(`t0-delti_g,t0-delti_g-j`) for j on the timestep
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g1 Source `herm_matrix_moving` object.
 * @param delti_g The relative time step at which to create the view. Must be in the range `[0, g.tc()]`.
 *
 * @pre The time cut-off (`tc_`) of both objects must be the same.
 * @pre The matrix sizes (`size1_`) must match.
 *
 * @note The function uses `memcpy` for efficient memory copying.
 */  
  template <typename T>
  void herm_matrix_timestep_moving<T>::set_timestep(herm_matrix_moving<T> &g1,int delti_g){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.set_timestep(g1,delti_g);      
  }
  
  
/**
 * @brief Copies the data from a `herm_matrix_timestep_view` object.
 *
 * This function copies the retarded and lesser Green’s function components
 * from the given `herm_matrix_timestep_view` object into the current instance.
 * \f$ this \f$(`t0,t0-j`) set to \f$ G \f$(`tstp_G,tstp_G-j`) for j on the timestep
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g Source `herm_matrix_timestep_view` object.
 * @param gcc Hermitian conjugate of g.
 *
 * @pre The matrix sizes (`size1_`) must match.
 *
 * @note The function uses `memcpy` for efficient memory copying.
 */
  template <typename T>
  void herm_matrix_timestep_moving<T>::set_timestep(herm_matrix_timestep_view<T> &G,herm_matrix_timestep_view<T> &Gcc){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.set_timestep(G,Gcc);
  }

/**
 * @brief Copies the data from a `herm_matrix_timestep` object.
 *
 * This function copies the retarded and lesser Green’s function components
 * from the given `herm_matrix_timestep` object into the current instance.
 * \f$ this \f$(`t0,t0-j`) set to \f$ G \f$(`tstp_G,tstp_G-j`) for j on the timestep
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g Source `herm_matrix_timestep` object.
 * @param gcc Hermitian conjugate of g.
 *
 * @pre The matrix sizes (`size1_`) must match.
 *
 * @note The function uses `memcpy` for efficient memory copying.
 */
  template <typename T>
  void herm_matrix_timestep_moving<T>::set_timestep(herm_matrix_timestep<T> &G,herm_matrix_timestep<T> &Gcc){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.set_timestep(G,Gcc);
  }

/**
 * @brief Copies the data from a `herm_matrix` object at its timestep `delti_g`
 *
 * This function copies the retarded and lesser Green’s function components
 * from the given `herm_matrix` object at time step `tstp_G` into the current instance:
 * \f$ this \f$(`t0,t0-j`) set to \f$ G \f$(`tstp_G,tstp_G-j`) for j on the timestep
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param g Source `herm_matrix_timestep` object.
 * @param gcc Hermitian conjugate of g.
 * @param tstp_G Time step of g.
 *
 * @pre The matrix sizes (`size1_`) must match.
 *
 * @note The function uses `memcpy` for efficient memory copying.
 */
  template <typename T>
  void herm_matrix_timestep_moving<T>::set_timestep(herm_matrix<T> &G,herm_matrix<T> &Gcc,int tstp_g){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.set_timestep(G,Gcc,tstp_g);
  } 
  
/**
 * @brief Sets a specific matrix element from a `herm_matrix_timestep_moving_view` object.
 *
 * This function copies the matrix element located at `(j1, j2)` in the source object `g`
 * to the position `(i1, i2)` in the current object  for all times on the timestep:
 * \f$ this_{i1,i2} \f$ (`t0,t0-j`) set to \f$ G_{j1,j2} \f$ (`t0,t0-j`)for all `j`
 *
 * @tparam T Floating-point type (`double` or `float`).
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
  void herm_matrix_timestep_moving<T>::set_matrixelement(int i1, int i2, herm_matrix_timestep_moving_view<T> &G,int j1, int j2){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.set_matrixelement(i1,i2,G,j1,j2);      
  }
/**
 * @brief Sets a specific matrix element from another `herm_matrix_timestep_moving` object.
 *
 * This function copies the matrix element located at `(j1, j2)` in the source object `g`
 * to the position `(i1, i2)` in the current object for all times on the timestep:
 * \f$ this_{i1,i2}\f$ (`t0,t0-j`)set to \f$ G_{j1,j2}\f$ (`t0,t0-j`) for all `j`
 *
 * @tparam T Floating-point type (`double` or `float`).
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
  void herm_matrix_timestep_moving<T>::set_matrixelement(int i1, int i2, herm_matrix_timestep_moving<T> &G,int j1, int j2){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.set_matrixelement(i1,i2,G,j1,j2);    }
/**
 * @brief Sets a specific matrix element from another `herm_matrix_moving` object at a given relative timestep `delti_g`
 *
 * This function copies the matrix element located at `(j1, j2)` in the source object `g`
 * to the position `(i1, i2)` in the current object for all times on the timestep:
 * \f$ this_{i1,i2}\f$ (`t0,t0-j`)set to \f$ G_{j1,j2}\f$ (`t0-delti_g,t0-delti_g-j`) for all j
 *
 * @tparam T Floating-point type (`double` or `float`).
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
  void herm_matrix_timestep_moving<T>::set_matrixelement(int i1, int i2, herm_matrix_moving<T> &G,int delti_g,int j1, int j2){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.set_matrixelement(i1,i2,G,delti_g,j1,j2);
  }
/**
 * @brief Increments the Green's function by another Green's function at the given timestep, scaled by a factor.
 *
 * This function performs the operation:
 * \f$ this\f$ (`t0,t0-j`) -->   \f$ this\f$ (`t0,t0-j`) + `\alpha` * \f$ G\f$(`t0,t0-j`) for all j on the timestep
 * where g is another Green’s function of the same dimensions,
 * and `\alpha` is a complex scaling factor. The increment is applied to both the 
 * retarded and lesser Green’s function components.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param G The Green’s function G(t,t') at given time slice to be added.
 * @param alpha the complex scaling factor
 *
 * @pre `G` must have the same `tc_`, `size1_`, and `size2_` as this object.
 * @post The retarded and lesser components are incremented by `\alpha * G(t,t')`.
 */
  template <typename T>
  void herm_matrix_timestep_moving<T>::incr_timestep(herm_matrix_timestep_moving_view<T> &G,cplx alpha){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.incr_timestep(G,alpha);
  }
/**
 * @brief Increments the Green's function by another Green's function at the given timestep, scaled by a factor.
 *
 * This function performs the operation:
 * \f$ this\f$ (`t0,t0-j`) -->   \f$ this\f$ (`t0,t0-j`) + \alpha * \f$ G\f$ (`t0,t0-j`) for all j on the tiemstep
 * where g is another Green’s function of the same dimensions,
 * and `\alpha` is a complex scaling factor. The increment is applied to both the 
 * retarded and lesser Green’s function components.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param G The Green’s function G(t,t') at given time slice to be added.
 * @param alpha the complex scaling factor
 *
 * @pre `G` must have the same `tc_`, `size1_`, and `size2_` as this object.
 * @post The retarded and lesser components are incremented by `\alpha \cdot g(t,t')`.
 */
  template <typename T>
  void herm_matrix_timestep_moving<T>::incr_timestep(herm_matrix_timestep_moving<T> &G,cplx alpha){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.incr_timestep(G,alpha);
  }
/**
 * @brief Increments the current timeslive by the timeslice `delti_g` of another Green's function
 *
 * This function performs the operation:
 * \f$ this\f$ (`t0,t0-j`) -->   \f$ this\f$ (`t0,t0-j`) + `\alpha` * \f$ G \f$(`t0-delti_g,t0-delti_g-j`) for all `j` on the timestep
 * where g is another Green’s function of the same dimensions,
 * and `\alpha` is a complex scaling factor. The increment is applied to both the 
 * retarded and lesser Green’s function components.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param G The Green’s function G(`t,t'`) to be added.
 * @param delti_g timeslice of G
 * @param alpha the complex scaling factor
 *
 * @pre `g` must have the same `tc_`, `size1_`, and `size2_` as this object.
 * @post The retarded and lesser components are incremented by `\alpha * g(t,t')`.
 */
  
  template <typename T>
  void herm_matrix_timestep_moving<T>::incr_timestep(herm_matrix_moving<T> &G,int delti_g,cplx alpha){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.incr_timestep(G,delti_g,alpha);
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
  void herm_matrix_timestep_moving<T>::left_multiply(function_moving<T> &ft){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.left_multiply(ft);
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
  void herm_matrix_timestep_moving<T>::right_multiply(function_moving<T> &ft){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.right_multiply(ft);
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
  void herm_matrix_timestep_moving<T>::left_multiply_hermconj(function_moving<T> &ft){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.left_multiply_hermconj(ft);
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
  void herm_matrix_timestep_moving<T>::right_multiply_hermconj(function_moving<T> &ft){
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.right_multiply_hermconj(ft);
  }
/**
 * @brief Scales the Green's function components by a given weight.
 *
 * This function multiplies all elements of the retarded and lesser Green's function 
 * components by the provided scalar `weight`.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param weight The scalar factor by which all elements are multiplied.
 *
 * @post Both the retarded and lesser Green's function components are scaled by `weight`.
 */  
  template <typename T>
  void herm_matrix_timestep_moving<T>::smul(cplx weight)
  {
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.smul(weight);
  }
/**
 * @brief Scales the Green's function components by a given weight.
 *
 * This function multiplies all elements of the retarded and lesser Green's function 
 * components by the provided scalar `weight`.
 *
 * @tparam T Floating-point type (`double` or `float`).
 * @param weight The scalar factor by which all elements are multiplied.
 *
 * @post Both the retarded and lesser Green's function components are scaled by `weight`.
 */  
  template <typename T>
  void herm_matrix_timestep_moving<T>::smul(T weight)
  {
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.smul(weight);
  }

#if CNTR_USE_MPI == 1
  template <typename T>
  void herm_matrix_timestep_moving<T>::MPI_Reduce(int root) {
    herm_matrix_timestep_moving_view<T> selfview(*this);
    selfview.MPI_Reduce(root);
  }
#endif



  template<typename T> 
  void herm_matrix_timestep_moving<T>::print_to_file(const char *file,int precision){
    int i,l,sg=element_size_;
    std::ofstream out;
    out.open(file,std::ios::out);
    out.precision(precision);
    out << "# " << " " << tc_ << " " << size1_ << " " << " " << sig_ << std::endl;
    if(tc_>=0){
      for(i=0;i<=tc_;i++){
	out << "ret: " << i ;
	for(l=0;l<sg;l++) out << " " << retptr(i)[l].real() << " " << retptr(i)[l].imag();
	out << std::endl;
      }
      out << std::endl;
      for(i=0;i<=tc_;i++){
	out << "les: " << i ;
	for(l=0;l<sg;l++) out << " " << lesptr(i)[l].real() << " " << lesptr(i)[l].imag();
	out << std::endl;
      }
      out << std::endl;
    }
    out.close();
  }
  template<typename T> 
  void herm_matrix_timestep_moving<T>::read_from_file(const char *file){
    int i,tc,l,size1,sg,sig;
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
	out >> s >> s ;
	for(l=0;l<sg;l++){
	  if(!( out >> real >> imag )){
	    std::cerr << "read G from file " << file << " error at ret (" << i  << ")"<< std::endl; 
	    abort();
	  }
	  retptr(i)[l] = std::complex<T>(real, imag);
	}
      }
      for(i=0;i<=tc_;i++){
	out >> s >> s ;
	for(l=0;l<sg;l++){
	  if(!( out >> real >> imag )){
	    std::cerr << "read G from file " << file << " error at les (" << i << ")"<< std::endl; 
	    abort();
	  }
	  lesptr(i)[l] = std::complex<T>(real, imag);
	}
      }
    }
    out.close();
  }
  
  
#if CNTR_USE_HDF5 == 1
  
  template <typename T>
  void herm_matrix_timestep_moving<T>::write_to_hdf5(hid_t group_id) {
    store_int_attribute_to_hid(group_id, std::string("herm_matrix_timestep_moving"), 1);
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
      shape[0] = (tc_+1)*element_size_;
      //store_cplx_array just to init datastructure
      store_cplx_array_to_hid(group_id, std::string("ret"), retptr(0), shape, len_shape);
      store_cplx_array_to_hid(group_id, std::string("les"), lesptr(0), shape, len_shape);
    }
  }
  template <typename T>
  void herm_matrix_timestep_moving<T>::write_to_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = create_group(group_id, groupname);
    this->write_to_hdf5(sub_group_id);
    close_group(sub_group_id);
  }
  template <typename T>
  void herm_matrix_timestep_moving<T>::write_to_hdf5(const char *filename,
					    const char *groupname, h5_mode mode) {
    hid_t file_id = open_hdf5_file(filename,mode);
    this->write_to_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
  }
  template <typename T>
  void herm_matrix_timestep_moving<T>::read_from_hdf5(hid_t group_id) {
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
      hsize_t ret_size = (tc_ + 1)   * element_size_;
      hsize_t les_size = (tc_ + 1)   * element_size_;
      read_primitive_type_array(group_id, "ret", ret_size, retptr(0));
      read_primitive_type_array(group_id, "les", les_size, lesptr(0));
    }
  }
  template <typename T>
  void herm_matrix_timestep_moving<T>::read_from_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(sub_group_id);
    close_group(sub_group_id);
  }
  template <typename T>
  void herm_matrix_timestep_moving<T>::read_from_hdf5(const char *filename,const char *groupname) {
    hid_t file_id = read_hdf5_file(filename);
    this->read_from_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
  }
  
  
#endif


  
} // namespace cntr

#endif  // CNTR_HERM_MATRIX_TIMESTEP_IMPL_H


