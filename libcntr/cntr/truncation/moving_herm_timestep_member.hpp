namespace cntr{
template <typename T> moving_herm_timestep<T>::moving_herm_timestep(){
   data_=0;
   les_=0;
   ret_=0;
   tc_=-1;
   size1_=0;
   size2_=0;
   element_size_=0;
   sig_=-1;
}
template <typename T> moving_herm_timestep<T>::~moving_herm_timestep(){
   if (data_!=0) delete [] data_;
}
template <typename T> moving_herm_timestep<T>::moving_herm_timestep(int tc,int size1,int sig){
    CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,-1,tc,__PRETTY_FUNCTION__)
    CNTR_ASSERT(MOVING_HERM_ASSERT_LEVEL,!(tc==-1 && size1>0),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,1,sig*sig,__PRETTY_FUNCTION__)
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
template <typename T> moving_herm_timestep<T>::moving_herm_timestep(const moving_herm_timestep &g){
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
        long ndata1=(tc_+1)*element_size_;
        data_ = new cplx [2*ndata1];
        ret_ = data_;
        les_ = data_+ndata1;
        memcpy(data_, g.data_, 2*sizeof(cplx)*ndata1);
    }
}
template <typename T> moving_herm_timestep<T>::moving_herm_timestep(int n,moving_herm<T> &g){
    int tc=g.tc();
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,n,tc,__PRETTY_FUNCTION__)
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
        memcpy(ret_, g.retptr(n,0), sizeof(cplx)*ndata1);
        memcpy(les_, g.lesptr(n,0), sizeof(cplx)*ndata1);
    }
}
template <typename T>  moving_herm_timestep<T> &  moving_herm_timestep<T>::operator=(const  moving_herm_timestep &g){
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
template <typename T> void moving_herm_timestep<T>::clear(void){
	if(tc_==-1) return;
	memset(data_, 0, 2*sizeof(cplx)*(tc_+1)*element_size_);
}
template <typename T>  void  moving_herm_timestep<T>::resize(int tc,int size1){
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
#define moving_herm_timestep_READ_ELEMENT {int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++) M(r,s)=x[r*dim+s];}
#define moving_herm_timestep_READ_ELEMENT_MINUS_CONJ {cplx w;int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++){ w=x[s*dim+r];M(r,s)=std::complex<T>(-w.real(),w.imag());}}
template<typename T> template <class Matrix> void moving_herm_timestep<T>::get_les(int j,Matrix &M){
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,j,tc_,__PRETTY_FUNCTION__)
    cplx *x;
    x=lesptr(j);
    moving_herm_timestep_READ_ELEMENT
}
template<typename T> template <class Matrix> void moving_herm_timestep<T>::get_ret(int j,Matrix &M){
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,j,tc_,__PRETTY_FUNCTION__)
    cplx *x;
    x=retptr(j);
    moving_herm_timestep_READ_ELEMENT
}
template<typename T> template <class Matrix> void moving_herm_timestep<T>::get_gtr(int j,Matrix &M){
  Matrix M1;
  get_ret(j,M);
  get_les(j,M1);
  M += M1;
}
template<typename T> inline void moving_herm_timestep<T>::get_ret(int j,cplx &x){
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,j,tc_,__PRETTY_FUNCTION__)
    x=*retptr(j);
}
template<typename T> inline void moving_herm_timestep<T>::get_les(int j,cplx &x){
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,j,tc_,__PRETTY_FUNCTION__)
    x=*lesptr(j);
}
template<typename T> inline void moving_herm_timestep<T>::get_gtr(int j,cplx &x){
   cplx x1;
   get_ret(j,x);
   get_les(j,x1);
   x+=x1;
}
template<typename T> std::complex<T> moving_herm_timestep<T>::density_matrix(){
   cplx x1;
   get_les(0,x1);
   return std::complex<T>(0.0,sig_)*x1;
}
template<typename T> template<class Matrix> void moving_herm_timestep<T>::density_matrix(Matrix &M){
    get_les(0,M);
    M *= std::complex<T>(0.0,1.0*sig_);
}
template<typename T> inline void moving_herm_timestep<T>::set_sig(int s){
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,s*s,1,__PRETTY_FUNCTION__)
    sig_=s;
}


template <typename T>
void moving_herm_timestep<T>::incr_timestep(moving_herm_timestep<T> &g,std::complex<T> alpha){
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size1_,g.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size2_,g.size2(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc_,g.tc(),__PRETTY_FUNCTION__)
    int ndata1=(tc_+1)*element_size_;
    for(int l=0;l<ndata1;l++){
        retptr(0)[l] +=  alpha*g.retptr(0)[l];
        lesptr(0)[l] +=  alpha*g.lesptr(0)[l];
    }
}

template <typename T>
void moving_herm_timestep<T>::left_multiply(moving_function<T> &ft,T weight){
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size1_,ft.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size2_,ft.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc_,ft.tc(),__PRETTY_FUNCTION__)
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
void moving_herm_timestep<T>::right_multiply(moving_function<T> &ft,T weight){
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size1_,ft.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size2_,ft.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc_,ft.tc(),__PRETTY_FUNCTION__)
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



}

