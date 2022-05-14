namespace cntr{

//  k-th order polynomila extrapolate of realtime functions (G^ret, G^vt, G^les)
//  to time 0, using information at times t=-j  [j=0...k]
template <typename T,int SIZE1>
void extrapolate_timestep_dispatch(moving_herm<T> &G,integration::Integrator<T> &I){
    typedef std::complex<T> cplx;
    int tc=G.tc();
    int k=I.k();
    int sg=G.element_size();
    int size1=G.size1();
    cplx *gtemp=new cplx [sg];
    T *p1=new T [k+1];
    if(k==0){
        for(int j=0;j<=tc;j++){
            element_set<T,SIZE1>(size1,G.retptr(0,j),G.retptr(1,j));
            element_set<T,SIZE1>(size1,G.lesptr(0,j),G.lesptr(1,j));
        }
    }else{
        for(int m=0;m<=k;m++){
            p1[m]=0.0;
            for(int l=0;l<=k;l++) p1[m]+=(1-2*(l%2))*I.poly_interpolation(l,m);
        }
        for(int j=0;j<=tc;j++){
            element_set_zero<T,SIZE1>(size1,gtemp);
            for(int l=0;l<=k;l++) element_incr<T,SIZE1>(size1,gtemp,p1[l],G.retptr(l,j));
            element_set<T,SIZE1>(size1,G.retptr(0,j),gtemp);
            element_set_zero<T,SIZE1>(size1,gtemp);
            for(int l=0;l<=k;l++) element_incr<T,SIZE1>(size1,gtemp,p1[l],G.lesptr(l,j));
            element_set<T,SIZE1>(size1,G.lesptr(0,j),gtemp);
        }
    }
    delete [] gtemp;
    delete [] p1;
}
template <typename T> 
void extrapolate_timestep(moving_herm<T> &G,integration::Integrator<T> &I){
    int k=I.k();
    int tc=G.tc();
    int size1=G.size1();
    CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,k+1,tc,__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,k+1,G.nt(),__PRETTY_FUNCTION__)
    if(size1) extrapolate_timestep_dispatch<T,1>(G,I);
    else extrapolate_timestep_dispatch<T,LARGESIZE>(G,I);
}

template <typename T>
T distance_norm2_moving_herm_dispatch(int n,int size,std::complex<T> *ret1,std::complex<T> *ret2){
	T err=0.0;
    int sg=size*size;
	std::complex<T> *temp= new std::complex<T> [sg];
	for(int i=0;i<=n;i++){
		element_set<T,LARGESIZE>(size,temp,ret1+i*sg);
		element_incr<T,LARGESIZE>(size,temp,-1.0,ret2+i*sg);
		err+=element_norm2<T,LARGESIZE>(size,temp);
	}
	delete [] temp;
	return err;
}

template <typename T>
T distance_norm2(int j,moving_herm<T> &g1,int j2,moving_herm<T> &g2){
	int size1=g1.size1();
    int tc=g1.tc();
	T err=0.0;
	CNTR_ASSERT_LESEQ_3(CNTR_ASSERT_LEVEL_0,0,j,g1.nt(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(CNTR_ASSERT_LEVEL_0,0,j2,g2.nt(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(CNTR_ASSERT_LEVEL_0,tc,g2.tc(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(CNTR_ASSERT_LEVEL_0,size1,g2.size1(),__PRETTY_FUNCTION__)
    err+=distance_norm2_moving_herm_dispatch(tc,size1,g1.retptr(j,0),g2.retptr(j2,0));
    err+=distance_norm2_moving_herm_dispatch(tc,size1,g1.lesptr(j,0),g2.lesptr(j2,0));
	return err;
}
}
