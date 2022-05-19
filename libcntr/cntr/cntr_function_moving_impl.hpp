#ifndef CNTR_FUNCTION_MOVING_IMPL_H
#define CNTR_FUNCTION_MOVING_IMPL_H

#include "cntr_function_decl.hpp"
#include "linalg.hpp"
#include "cntr_elements.hpp"


namespace cntr {

/* #######################################################################################
#
#   CONSTRUCTION/DESTRUCTION
#
########################################################################################*/
  template <typename T> 
  function_moving<T>::function_moving(){
   data_=0;
   value_=0;
   tc_=-1;
   t0_=0;
   size1_=0;
   size2_=0;
   element_size_=0;
  }

  template <typename T> 
  function_moving<T>::~function_moving(){
   if (data_!=0) delete [] data_;
   if (value_!=0)  delete [] value_;
  }
/** \brief <b> Initializes the `function_moving` class for a function \f$ f(t)\f$ with one real time argument.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Initializes the `function_moving` class for a function \f$ f(t)\f$ with one real time argument.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tc
* > Number of timesteps that are stored in the function.
* @param size1
* > Matrix rank of the function.
*/
  template <typename T> 
  function_moving<T>::function_moving(int tc,int size1){
    assert(-1<=tc);
    // CNTR_ASSERT(MOVING_HERM_ASSERT_LEVEL,!(tc==-1 && size1>0),__PRETTY_FUNCTION__)
    tc_=tc;
    t0_=0;
    size1_=size1;
    size2_=size1;
    element_size_=size1*size1;
    if(tc_==-1){
        data_=0;
        value_=0;
    }else{
        // here tc>=0 AND size>0
        long ndata1=(tc_+1)*element_size_;
        data_ = new cplx [ndata1];
        value_ = new cplx* [tc_+1];
        memset(data_, 0, sizeof(cplx)*ndata1);
        for(int t=0;t<=tc_;t++) value_[t]=data_+t;
    }
  }
/** \brief <b> Initializes the `function_moving` class with the same layout as a given `function_moving`.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Initializes the `function_moving` class with the same cutoff time `tc`, timestep `t0`, number of colums `size1`, number of rows `size2` and `element_size`.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param g
* > The `function_moving` according to which the class should be initialized.
*/
  template <typename T> 
  function_moving<T>::function_moving(const function_moving &g){
    tc_=g.tc_;
    t0_=g.t0_;
    size1_=g.size1_;
    size2_=g.size2_;
    element_size_=g.element_size_;
    if(tc_==-1){
        data_=0;
        value_=0;
    }else{
        // here tc>0 AND size>0
        long ndata1=(tc_+1)*element_size_;
        data_ = new cplx [ndata1];
        value_ = new cplx* [tc_+1];
        memcpy(data_, g.data_, sizeof(cplx)*ndata1);
        // correctly redirect the pointers
        for(int t=0;t<=tc_;t++){
            value_[t]=data_+(g.value_[t]-g.data_);
        }
    }
  }

  template <typename T>  
  function_moving<T> &  function_moving<T>::operator=(const  function_moving &g){
    if(this==&g) return *this;
    t0_=g.t0_;
    if( tc_!=g.tc_ || size1_!=g.size1_ || size2_!=g.size2_){
        // reallocate
        if (data_!=0) delete [] data_;
        if (value_!=0)  delete [] value_;
        tc_=g.tc_;
        size1_=g.size1_;
        size2_=g.size2_;
        element_size_=g.element_size_;
        if(tc_>=0){
            // here tc>0 AND size>0
            long ndata1=(tc_+1)*element_size_;
            data_ = new cplx [ndata1];
            value_ = new cplx* [tc_+1];
        }else{
            data_=0;
            value_=0;
        }
    }
    if(tc_>=0){
        memcpy(data_, g.data_, sizeof(cplx)*(tc_+1)*element_size_);
        for(int t=0;t<=tc_;t++){
            value_[t]=data_+(g.value_[t]-g.data_);
        }
    }
    return *this;
  }
  
/** \brief <b> Set all data to zero for the `function_moving` class</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set all data to zero for the `function` class.
*
* <!-- ARGUMENTS
*      ========= -->
*
*/

  template <typename T> void function_moving<T>::clear(void){
  if(tc_==-1) return;
  memset(data_, 0, sizeof(cplx)*(tc_+1)*element_size_);
    for(int t=0;t<=tc_;t++) value_[t]=data_+t;
}

/** \brief <b> Set the current timestep of the `function_moving` class object</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set the current timestep `t0` of the `function_moving` class to `t0`.
*
* <!-- ARGUMENTS
*      ========= -->
* @param t0
* > the new internal timestep `t0`.
*/
template <typename T> void function_moving<T>::set_t0(int t0){
  assert(tc_<=t0);
    t0_=t0;
}

/** \brief <b> Resize the time length and/or matrix dimension of a square-matrix function</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Resize the time length and/or matrix dimension of a square-matrix function
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tc
* > the new number of time-steps.
* @param size1
* > the new leading dimension of the matrix
*/
 
template <typename T>  void  function_moving<T>::resize(int tc,int size1){
    if( tc!=tc_ || size1!=size1_){
        // reallocate
        if (data_!=0) delete [] data_;
        if (value_!=0)  delete [] value_;
        tc_=tc;
        size1_=size1;
        size2_=size1;
        element_size_=size1*size1;
        if(tc_>0){
            long ndata1=(tc_+1)*element_size_;
            data_ = new cplx [ndata1];
            value_ = new cplx* [tc_+1];
            memset(data_, 0, sizeof(cplx)*ndata1);
            for(int t=0;t<=tc_;t++) value_[t]=data_+t;
        }else{
            data_=0;
            value_=0;
        }
    }
    clear();
}

// READING ELEMENTS TO ANY MATRIX TYPE OR TO COMPLEX NUMBERS
// (then only the (0,0) element is addressed for dim>0)
  #define function_moving_READ_ELEMENT {int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++) M(r,s)=x[r*dim+s];}
  #define function_moving_READ_ELEMENT_MINUS_CONJ {cplx w;int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++){ w=x[s*dim+r];M(r,s)=std::complex<T>(-w.real(),w.imag());}}

/** \brief <b> Get matrix value of this `function_moving` at a specific time step inde`</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Get matrix value of this `function_moving` object at a specific time index.* > 
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > the index of the timestep 
* @param M
* > the Eigen matrix in which the value will be stored
*/
  template<typename T> template <class Matrix>
  void function_moving<T>::get_value(int i,Matrix &M) const{
    assert(0<=i && i<=tc());
    cplx *x;
    x=ptr(i);
    function_moving_READ_ELEMENT
  }
/** \brief <b> Set matrix value at a specific time point</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Set scalar value at a specific time point.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param i
* > the index of the timestep to be set
* @param M
* > the new matrix value of function at timestep index i as an EigenMatrix.
*/
  template<typename T> template <class Matrix> 
  void function_moving<T>::set_value(int i,Matrix &M){
    assert(0<=i && i<=tc());
    assert(M.rows()==size1_);
    assert(M.cols()==size2_);
    int i1,i2;
    cplx *ft;
    ft=ptr(i);
    for(i1=0;i1<size1_;i1++) for(i2=0;i2<size2_;i2++) ft[i1*size2_+i2]=M(i1,i2);
  }

/** \brief <b> Store the `function_moving` data to a txt file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Store the `function_moving` data to a file. Data are stored in txt form.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param file
* > the file name
* @param precision
* > the precision of stored `function` data
*/
  template<typename T> 
  void function_moving<T>::print_to_file(const char *file,int precision){
    int i,l,sg=element_size_;
    std::ofstream out;
    out.open(file,std::ios::out);
    out.precision(precision);
    out << "# " << t0_ << " " << tc_ << " " << size1_ << " " << std::endl;
    if(tc_>=0){
      for(i=0;i<=tc_;i++){
                for(l=0;l<sg;l++) out << ptr(i)[l].real() << " " << ptr(i)[l].imag();
                out << std::endl;
      }
        }
    out.close();
  }

/** \brief <b> Read the `function_moving` data from a txt file</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Read the `function_moving` from a txt file. The data file is probably written by print_to_file method,
* > or is assumed to have the same format.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param file
* > the file name
*/
 
  template<typename T> 
  void function_moving<T>::read_from_file(const char *file){
    int i,t0,tc,l,size1,sg;
    double real, imag;
    std::string s;
    std::ifstream out;
    out.open(file,std::ios::in);
    if(!(out >> s >> t0 >> tc >> size1)){
      std::cerr << "read MOVING_FUNC from file " << file << " error in file" << std::endl;
      abort();
    }
    if(tc!=tc_ || size1!=size1_) resize(tc,size1);
    set_t0(t0);
    sg=element_size_;
    if(tc_>=0){
        for(i=0;i<=tc_;i++){
           for(l=0;l<sg;l++){
             if(!( out >> real >> imag )){
               std::cerr << "read MOVING_FUNC from file " << file << " error at ret (" << i << ")"<< std::endl;
               abort();
             }
             ptr(i)[l] = std::complex<T>(real, imag);
           }
        }
    }
    out.close();
  }
/** \brief <b> Propagates the `function_moving` data by one index.</b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Shifts the data values in the `function_moving` value vector by one index.
* > The final value `value[tc]` is shifted cyclically to `value[0]`
*
* <!-- ARGUMENTS
*      ========= -->
*
*/
  template <typename T>
  void function_moving<T>::forward(void){
    if(tc_>0){
        cplx* tmp1=value_[tc_];
        for(int t=tc_;t>0;t--){
            value_[t]=value_[t-1];
        }
        value_[0]=tmp1;
    }
  }
/** \brief <b>Set the value of `funcion_moving` class from a `function` class object. </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > Copies the first tc timesteps of a `function` object into this `function_moving` object. 
*
* <!-- ARGUMENTS
*      ========= -->
* @param tstp
* > Highest tstp read from the 'function' object.
* @param f
* > `function` object from which the data are copied from.
*/

  template <typename T> 
  void function_moving<T>::set_from_function_backward(int tstp, function<T>& f){
    assert(tc_<=tstp && tstp <= f.nt());
    for(int i=0;i<=tc_;i++){
      MatrixXcd tmp;
      f.get_value(i,tmp);
      set_value(tstp-i,tmp);
    }
  }

} // namespace cntr

#endif  // CNTR_FUNCTION_MOVING_IMPL_H
