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

  template <typename T> 
  void herm_matrix_timestep_moving<T>::clear(void){
    if(tc_==-1) return;
    memset(data_, 0, 2*sizeof(cplx)*(tc_+1)*element_size_);
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
    clear();
  }
  // READING ELEMENTS TO ANY MATRIX TYPE OR TO COMPLEX NUMBERS
  // (then only the (0,0) element is addressed for dim>0)
#define herm_matrix_timestep_moving_READ_ELEMENT {int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++) M(r,s)=x[r*dim+s];}
#define herm_matrix_timestep_moving_READ_ELEMENT_MINUS_CONJ {cplx w;int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++){ w=x[s*dim+r];M(r,s)=std::complex<T>(-w.real(),w.imag());}}
  template<typename T> template <class Matrix> 
  void herm_matrix_timestep_moving<T>::get_les(int j,Matrix &M){
    assert(0<=j && j<=tc_);
    cplx *x;
    x=lesptr(j);
    herm_matrix_timestep_moving_READ_ELEMENT
      }
  template<typename T> template <class Matrix> 
  void herm_matrix_timestep_moving<T>::get_ret(int j,Matrix &M){
    assert(0<=j && j<=tc_);
    cplx *x;
    x=retptr(j);
    herm_matrix_timestep_moving_READ_ELEMENT
      }

  template<typename T> template <class Matrix> 
  void herm_matrix_timestep_moving<T>::get_gtr(int j,Matrix &M){
    Matrix M1;
    get_ret(j,M);
    get_les(j,M1);
    M += M1;
  }

  template<typename T> 
  inline void herm_matrix_timestep_moving<T>::get_ret(int j,cplx &x){
    assert(0<=j && j<=tc_);
    x=*retptr(j);
  }

  template<typename T> 
  inline void herm_matrix_timestep_moving<T>::get_les(int j,cplx &x){
    assert(0<=j && j<=tc_);
    x=*lesptr(j);
  }

  template<typename T> 
  inline void herm_matrix_timestep_moving<T>::get_gtr(int j,cplx &x){
    cplx x1;
    get_ret(j,x);
    get_les(j,x1);
    x+=x1;
  }

  template<typename T> 
  std::complex<T> herm_matrix_timestep_moving<T>::density_matrix(){
    cplx x1;
    get_les(0,x1);
    return std::complex<T>(0.0,sig_)*x1;
  }

  template<typename T> template<class Matrix> 
  void herm_matrix_timestep_moving<T>::density_matrix(Matrix &M){
    get_les(0,M);
    M *= std::complex<T>(0.0,1.0*sig_);
  }

  template<typename T> 
  inline void herm_matrix_timestep_moving<T>::set_sig(int s){
    assert(s*s==1);
    sig_=s;
  }


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

  template <typename T>
  void herm_matrix_timestep_moving<T>::smul(T weight)
  {
    int ndata = element_size_ * (tc_ + 1) * 2;
    for (int i = 0; i < ndata; ++i)
      {
	data_[i] *= weight;
      }
  }

  template <typename T>
  void herm_matrix_timestep_moving<T>::smul(cplx weight)
  {
    int ndata = element_size_ * (tc_ + 1) * 2;
    for (int i = 0; i < ndata; ++i)
      {
	data_[i] *= weight;
      }
  }
  
  template <typename T>
  void herm_matrix_timestep_moving<T>::incr(function_moving<T> &g1, T weight){
    int m;
    int total_size_=2*(tc_ + 1)*element_size_;
    assert(g1.size1_ == size1_ && g1.tc_ == tc_);
    for (m = 0; m < total_size_; m++)
      data_[m] += weight * g1.data_[m];
  }

  template <typename T>
  void herm_matrix_timestep_moving<T>::incr(herm_matrix_timestep_moving<T> &g1, T weight){
    int m;
    int total_size_=2*(tc_ + 1)*element_size_;
    assert(g1.size1_ == size1_ && g1.tc_ == tc_);
    for (m = 0; m < total_size_; m++)
      data_[m] += weight * g1.data_[m];
  }

#if CNTR_USE_MPI == 1

  template <typename T>
  void herm_matrix_timestep_moving<T>::MPI_Reduce(int root) {
    int taskid;
    //int len1 = 2 * ((tstp_ + 1) * 2 + (ntau_ + 1)) * size1_ * size1_; This is the original version from timestep check if the reduction to the data window is correct!!
    int len1 = ((tc_ + 1) * 2 ) * element_size;
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
