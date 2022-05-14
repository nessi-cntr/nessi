namespace cntr{

template <typename T> moving_herm<T>::moving_herm(){
   data_=0;
   les_=0;
   ret_=0;
   tc_=-1;
   nt_=-1;
   size1_=0;
   size2_=0;
   element_size_=0;
   sig_=-1;
}
template <typename T> moving_herm<T>::~moving_herm(){
   if (data_!=0) delete [] data_;
   if (ret_!=0)  delete [] ret_;
   if (les_!=0)  delete [] les_;
}
template <typename T> moving_herm<T>::moving_herm(int tc,int nt,int size1,int sig){
    CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,-1,tc,__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,-1,nt,__PRETTY_FUNCTION__)
    CNTR_ASSERT(MOVING_HERM_ASSERT_LEVEL,!(tc==-1 && size1>0),__PRETTY_FUNCTION__)
    CNTR_ASSERT(MOVING_HERM_ASSERT_LEVEL,!(tc==-1 && nt>=0),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,1,sig*sig,__PRETTY_FUNCTION__)
    tc_=tc;
    nt_=nt;
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
        long ndata2=(tc_+1)*(nt_+1)*element_size_;
        long ndata1=(tc_+1)*element_size_;
        data_ = new cplx [2*ndata2];
        ret_ = new cplx* [nt_+1];
        les_ = new cplx* [nt_+1];
        memset(data_, 0, 2*sizeof(cplx)*ndata2);
        for(int t=0;t<=nt_;t++){
            ret_[t]=data_+t*ndata1;
            les_[t]=data_+(t+nt_+1)*ndata1;
        }
    }
}
template <typename T> moving_herm<T>::moving_herm(const moving_herm &g){
    tc_=g.tc_;
    nt_=g.nt_;
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
        long ndata2=(tc_+1)*(nt_+1)*element_size_;
        long ndata1=(tc_+1)*element_size_;
        data_ = new cplx [2*ndata2];
        ret_ = new cplx* [nt_+1];
        les_ = new cplx* [nt_+1];
        memcpy(data_, g.data_, 2*sizeof(cplx)*ndata2);
        // correctly redirect the pointers
        for(int t=0;t<=nt_;t++){
            ret_[t]=data_+(g.ret_[t]-g.data_);
            les_[t]=data_+(g.les_[t]-g.data_);
        }
    }
}
template <typename T>  moving_herm<T> &  moving_herm<T>::operator=(const  moving_herm &g){
    if(this==&g) return *this;
    sig_=g.sig_;
    if( tc_!=g.tc_ || nt_!=g.nt_ || size1_!=g.size1_ || size2_!=g.size2_){
        // reallocate
        if (data_!=0) delete [] data_;
        if (ret_!=0)  delete [] ret_;
        if (les_!=0)  delete [] les_;
        tc_=g.tc_;
        nt_=g.nt_;
        size1_=g.size1_;
        size2_=g.size2_;
        element_size_=g.element_size_;
        if(tc_>=0){
            // here tc>0 AND size>0
            long ndata2=(tc_+1)*(nt_+1)*element_size_;
            data_ = new cplx [2*ndata2];
            ret_ = new cplx* [nt_+1];
            les_ = new cplx* [nt_+1];
        }else{
            data_=0;
            les_=0;
            ret_=0;
        }
    }
    if(tc_>=0){
        memcpy(data_, g.data_, 2*sizeof(cplx)*(tc_+1)*(nt_+1)*element_size_);
        for(int t=0;t<=nt_;t++){
            ret_[t]=data_+(g.ret_[t]-g.data_);
            les_[t]=data_+(g.les_[t]-g.data_);
        }
    }
    return *this;
}
template <typename T> void moving_herm<T>::clear(void){
	if(tc_==-1) return;
	memset(data_, 0, 2*sizeof(cplx)*(tc_+1)*(nt_+1)*element_size_);
    long ndata1=(tc_+1)*element_size_;
    for(int t=0;t<=nt_;t++){
        ret_[t]=data_+t*ndata1;
        les_[t]=data_+(t+nt_+1)*ndata1;
    }
}
template <typename T>  void  moving_herm<T>::resize(int tc,int nt,int size1){
    CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,-1,tc,__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,-1,nt,__PRETTY_FUNCTION__)
    CNTR_ASSERT(MOVING_HERM_ASSERT_LEVEL,!(tc==-1 && size1>0),__PRETTY_FUNCTION__)
    CNTR_ASSERT(MOVING_HERM_ASSERT_LEVEL,!(tc==-1 && nt>=0),__PRETTY_FUNCTION__)
    if( tc!=tc_ ||nt!=nt_ || size1!=size1_){
        // reallocate
        if (data_!=0) delete [] data_;
        if (ret_!=0)  delete [] ret_;
        if (les_!=0)  delete [] les_;
        tc_=tc;
        nt_=nt;
        size1_=size1;
        size2_=size1;
        element_size_=size1*size1;
        if(tc_>=0){
            long ndata2=(tc_+1)*(nt_+1)*element_size_;
            long ndata1=(tc_+1)*element_size_;
            data_ = new cplx [2*ndata2];
            ret_ = new cplx* [nt_+1];
            les_ = new cplx* [nt_+1];
            memset(data_, 0, 2*sizeof(cplx)*ndata2);
            for(int t=0;t<=nt_;t++){
                ret_[t]=data_+t*ndata1;
                les_[t]=data_+(t+nt_+1)*ndata1;
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
#define moving_herm_READ_ELEMENT {int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++) M(r,s)=x[r*dim+s];}
#define moving_herm_READ_ELEMENT_MINUS_CONJ {cplx w;int r,s,dim=size1_;M.resize(dim,dim);for(r=0;r<dim;r++) for(s=0;s<dim;s++){ w=x[s*dim+r];M(r,s)=std::complex<T>(-w.real(),w.imag());}}
template<typename T> template <class Matrix> void moving_herm<T>::get_les(int i,int j,Matrix &M){
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,i,tc_,__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,j,tc_,__PRETTY_FUNCTION__)
    cplx *x;
    x=lesptr(i,j);
    moving_herm_READ_ELEMENT
}
template<typename T> template <class Matrix> void moving_herm<T>::get_ret(int i,int j,Matrix &M){
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,i,tc_,__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,j,tc_,__PRETTY_FUNCTION__)
    cplx *x;
    x=retptr(i,j);
    moving_herm_READ_ELEMENT
}
template<typename T> template <class Matrix> void moving_herm<T>::get_gtr(int i,int j,Matrix &M){
  Matrix M1;
  get_ret(i,j,M);
  get_les(i,j,M1);
  M += M1;
}
template<typename T> inline void moving_herm<T>::get_ret(int i,int j,cplx &x){ 
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,i,tc_,__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,j,tc_,__PRETTY_FUNCTION__)
    x=*retptr(i,j);
}
template<typename T> inline void moving_herm<T>::get_les(int i,int j,cplx &x){ 
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,i,tc_,__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,j,tc_,__PRETTY_FUNCTION__)
    x=*lesptr(i,j);
}
template<typename T> inline void moving_herm<T>::get_gtr(int i,int j,cplx &x){
   cplx x1;
   get_ret(i,j,x);
   get_les(i,j,x1);
   x+=x1;
}
template<typename T> std::complex<T> moving_herm<T>::density_matrix(int i){
   cplx x1;
   get_les(i,0,x1);
   return std::complex<T>(0.0,sig_)*x1;
}
template<typename T> template<class Matrix> void moving_herm<T>::density_matrix(int i,Matrix &M){
    get_les(i,0,M);
    M *= std::complex<T>(0.0,1.0*sig_);
}
////////////////////////////////
template<typename T> void moving_herm<T>::print_to_file(const char *file,int precision){
		int i,j,l,sg=element_size_;
		std::ofstream out;
		out.open(file,std::ios::out);
		out.precision(precision);
		out << "# " << tc_ << " " << nt_ << " " << size1_ << " " << " " << sig_ << std::endl;
		if(tc_>=0){
			for(i=0;i<=nt_;i++){
			   for(j=0;j<=tc_;j++){
				  out << "ret: " << i << " " << j;
				  for(l=0;l<sg;l++) out << " " << retptr(i,j)[l].real() << " " << retptr(i,j)[l].imag();
				  out << std::endl;
			   }
			   out << std::endl;
			}
			out << std::endl;
			for(i=0;i<=nt_;i++){
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
template<typename T> void moving_herm<T>::read_from_file(const char *file){
    int i,tc,j,l,size1,sg,sig,nt;
    double real, imag;
    std::string s;
    std::ifstream out;
    out.open(file,std::ios::in);
    if(!(out >> s >> tc >> nt >> size1 >> sig)){
      std::cerr << "read G from file " << file << " error in file" << std::endl; 
      abort();
    }
    if(tc!=tc_ || nt!=nt_ || size1!=size1_) resize(tc,size1);
    sig_=sig;
    sg=element_size_;
    if(tc_>=0){
        for(i=0;i<=nt_;i++){
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
        for(i=0;i<=nt_;i++){
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
void moving_herm<T>::forward(void){
    if(tc_>0){
        cplx* tmp1=les_[nt_];
        cplx* tmp2=ret_[nt_];
        for(int t=nt_;t>0;t--){
            les_[t]=les_[t-1];
            ret_[t]=ret_[t-1];
        }
        les_[0]=tmp1;
        ret_[0]=tmp2;
    }
}

// DATA exchange with HERM_MATRIX
    // read data to slice i (relative to t0)
template <typename T>
void moving_herm<T>::clear_timestep(int i){
    for(int t1=0;t1<=tc_;t1++){
        element_set_zero<T,LARGESIZE>(size1_,retptr(i,t1));
        element_set_zero<T,LARGESIZE>(size1_,lesptr(i,t1));
    }
}
template <typename T>
void moving_herm<T>::set_timestep(int i,herm_matrix_timestep_view<T> &g,herm_matrix_timestep_view<T> &gcc){
    // this.ret(i,j)=g.ret(tstp,i-j) for j=0...min(tc,i)
    // this.les(i,j)=g.les(tstp,i-j) for j=0...min(tc,i)
    // remaining entries are filled with zeros, if any
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size1_,g.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,i,tc_,__PRETTY_FUNCTION__)
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
void moving_herm<T>::set_timestep(int i,int tstp,herm_matrix<T> &g,herm_matrix<T> &gcc){
    herm_matrix_timestep_view<T> tmp(tstp,g);
    herm_matrix_timestep_view<T> tmp1(tstp,gcc);
	set_timestep(i,tmp,tmp1);
}

template <typename T>
void moving_herm<T>::set_from_G_backward(int tstp,herm_matrix<T> &g,herm_matrix<T> &gcc){
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,tc_,tstp,g.nt(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,nt_,tstp,g.nt(),__PRETTY_FUNCTION__)
    for(int i=0;i<=nt_;i++) set_timestep(i,tstp-i,g,gcc);
}

template <typename T>
void moving_herm<T>::set_from_G_backward(int tstp,herm_matrix<T> &g){
    set_from_G_backward(tstp,g,g);
}


template <typename T>
void moving_herm<T>::set_timestep(int i,moving_herm<T> &g,int j){
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size1_,g.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size2_,g.size2(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc_,g.tc(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,i,nt_,__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,j,g.nt(),__PRETTY_FUNCTION__)
    int ndata1=(tc_+1)*element_size_;
    memcpy(retptr(i,0),g.retptr(j,0), sizeof(cplx)*ndata1);
    memcpy(lesptr(i,0),g.lesptr(j,0), sizeof(cplx)*ndata1);
}

template <typename T>
void moving_herm<T>::get_timestep(int i,moving_herm_timestep<T> &g){
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,i,nt_,__PRETTY_FUNCTION__)
    g.resize(tc_,size1_);
    int ndata1=(tc_+1)*element_size_;
    memcpy(g.retptr(0),retptr(i,0), sizeof(cplx)*ndata1);
    memcpy(g.lesptr(0),lesptr(i,0), sizeof(cplx)*ndata1);
}



template <typename T>
void moving_herm<T>::incr_timestep(int i,moving_herm<T> &g,int j,cplx alpha){
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size1_,g.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size2_,g.size2(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc_,g.tc(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,i,nt_,__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,j,g.nt(),__PRETTY_FUNCTION__)
    int ndata1=(tc_+1)*element_size_;
    for(int l=0;l<ndata1;l++){
        retptr(i,0)[l] +=  alpha*g.retptr(j,0)[l];
        lesptr(i,0)[l] +=  alpha*g.lesptr(j,0)[l];
    }
}

template <typename T>
void moving_herm<T>::left_multiply(moving_function<T> &ft,T weight){
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size1_,ft.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size2_,ft.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc_,ft.tc(),__PRETTY_FUNCTION__)
    cplx *xtemp,*ftemp;
    xtemp=new cplx [element_size_];
    ftemp=ft.ptr(0);
    for(int j=0;j<=tc_;j++){
        element_mult<T,LARGESIZE>(size1_,xtemp,ftemp,retptr(0,j));
        element_smul<T,LARGESIZE>(size1_,xtemp,weight);
        element_set<T,LARGESIZE>(size1_,retptr(0,j),xtemp);
        element_mult<T,LARGESIZE>(size1_,xtemp,ftemp,lesptr(0,j));
        element_smul<T,LARGESIZE>(size1_,xtemp,weight);
        element_set<T,LARGESIZE>(size1_,lesptr(0,j),xtemp);
    }
   delete [] xtemp;
}
template <typename T>
void moving_herm<T>::right_multiply(moving_function<T> &ft,T weight){
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size1_,ft.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,size2_,ft.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc_,ft.tc(),__PRETTY_FUNCTION__)
    cplx *xtemp;
    xtemp=new cplx [element_size_];
    for(int j=0;j<=tc_;j++){
        element_mult<T,LARGESIZE>(size1_,xtemp,retptr(0,j),ft.ptr(j));
        element_smul<T,LARGESIZE>(size1_,xtemp,weight);
        element_set<T,LARGESIZE>(size1_,retptr(0,j),xtemp);
        element_mult<T,LARGESIZE>(size1_,xtemp,lesptr(0,j),ft.ptr(j));
        element_smul<T,LARGESIZE>(size1_,xtemp,weight);
        element_set<T,LARGESIZE>(size1_,lesptr(0,j),xtemp);
    }
   delete [] xtemp;
}

}

