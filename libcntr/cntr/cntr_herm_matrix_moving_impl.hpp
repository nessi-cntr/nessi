#ifndef CNTR_HERM_MATRIX_MOVING_IMPL_H
#define CNTR_HERM_MATRIX_MOVING_IMPL_H

#include "cntr_herm_matrix_decl.hpp"
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
  void herm_matrix_moving<T>::set_t0(int t0){
    assert(tc_<=t0);
    t0_=t0;
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
  template<typename T> template <class Matrix> 
  void herm_matrix_moving<T>::get_les(int i,int j,Matrix &M) const{
    assert(0<=i && i<=tc_);
    assert(0<=j && j<=tc_);
    cplx *x;
    x=lesptr(i,j);
    herm_matrix_moving_READ_ELEMENT
      }
  template<typename T> template <class Matrix> 
  void herm_matrix_moving<T>::get_ret(int i,int j,Matrix &M) const{
    assert(0<=i && i<=tc_);
    assert(0<=j && j<=tc_);
    cplx *x;
    x=retptr(i,j);
    herm_matrix_moving_READ_ELEMENT
      }

  template<typename T> template <class Matrix> 
  void herm_matrix_moving<T>::get_gtr(int i,int j,Matrix &M) const{
    Matrix M1;
    get_ret(i,j,M);
    get_les(i,j,M1);
    M += M1;
  }

  template<typename T> 
  inline void herm_matrix_moving<T>::get_ret(int i,int j,cplx &x) const{ 
    assert(0<=i && i<=tc_);
    assert(0<=j && j<=tc_);
    x=*retptr(i,j);
  }
  template<typename T> 
  inline void herm_matrix_moving<T>::get_les(int i,int j,cplx &x) const{ 
    assert(0<=i && i<=tc_);
    assert(0<=j && j<=tc_);
    x=*lesptr(i,j);
  }
  template<typename T> 
  inline void herm_matrix_moving<T>::get_gtr(int i,int j,cplx &x) const{
    cplx x1;
    get_ret(i,j,x);
    get_les(i,j,x1);
    x+=x1;
  }
  template<typename T> 
  std::complex<T> herm_matrix_moving<T>::density_matrix(int i){
    cplx x1;
    get_les(i,0,x1);
    return std::complex<T>(0.0,sig_)*x1;
  }
  template<typename T> 
  template<class Matrix> void herm_matrix_moving<T>::density_matrix(int i,Matrix &M){
    get_les(i,0,M);
    M *= std::complex<T>(0.0,1.0*sig_);
  }
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

  // DATA exchange with HERM_MATRIX
  // read data to slice i (relative to t0)
  template <typename T>
  void herm_matrix_moving<T>::clear_timestep(int i){
    for(int t1=0;t1<=tc_;t1++){
      element_set_zero<T,LARGESIZE>(size1_,retptr(i,t1));
      element_set_zero<T,LARGESIZE>(size1_,lesptr(i,t1));
    }
  }

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
  template <typename T>
  void herm_matrix_moving<T>::set_timestep(int i,int tstp,herm_matrix<T> &g,herm_matrix<T> &gcc){
    herm_matrix_timestep_view<T> tmp(tstp,g);
    herm_matrix_timestep_view<T> tmp1(tstp,gcc);
    set_timestep(i,tmp,tmp1);
  }

  template <typename T>
  void herm_matrix_moving<T>::set_from_G_backward(int tstp,herm_matrix<T> &g,herm_matrix<T> &gcc){
    assert(tc_<=tstp && tstp <= g.nt());
    for(int i=0;i<=tc_;i++) set_timestep(i,tstp-i,g,gcc);
  }

  template <typename T>
  void herm_matrix_moving<T>::set_from_G_backward(int tstp,herm_matrix<T> &g){
    set_from_G_backward(tstp,g,g);
  }


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

  template <typename T>
  void herm_matrix_moving<T>::get_timestep(int i,herm_matrix_timestep_moving<T> &g){
    assert(0<=i && i<=tc_);
    g.resize(tc_,size1_);
    g.set_sig(sig_);
    int ndata1=(tc_+1)*element_size_;
    memcpy(g.retptr(0),retptr(i,0), sizeof(cplx)*ndata1);
    memcpy(g.lesptr(0),lesptr(i,0), sizeof(cplx)*ndata1);
  }

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


  template <typename T>
  void herm_matrix_moving<T>::write_to_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = create_group(group_id, groupname);
    this->write_to_hdf5(sub_group_id);
    close_group(sub_group_id);
  }

  template <typename T>
  void herm_matrix_moving<T>::write_to_hdf5(const char *filename,
					    const char *groupname) {
    hid_t file_id = open_hdf5_file(filename);
    this->write_to_hdf5(file_id, groupname);
    close_hdf5_file(file_id);
  }


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
  template <typename T>
  void herm_matrix_moving<T>::read_from_hdf5(hid_t group_id, const char *groupname) {
    hid_t sub_group_id = open_group(group_id, groupname);
    this->read_from_hdf5(sub_group_id);
    close_group(sub_group_id);
  }

#endif


}



#endif  // CNTR_HERM_MATRIX_IMPL_H
