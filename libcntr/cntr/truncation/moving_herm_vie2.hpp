namespace cntr{



template < typename T, int SIZE1 >
void vie2_timestep_dispatch_ret(moving_herm<T> &G,moving_herm<T> &F,moving_herm<T> &Fcc,moving_herm<T> &Q, integration::Integrator<T> &I, T h){
    typedef std::complex<T> cplx;
    int k=I.get_k();
    // various asserts are one level higher
    cplx *mm,*qq,*one,*gtemp,*stemp,weight;
    int i,j,p,l,q;
    int size1=G.size1();
    int tc=G.tc();
    int sg=G.element_size();
    mm = new cplx [k*k*sg];
    qq = new cplx [k*sg];
    one = new cplx [sg];
    gtemp = new cplx [k*sg];
    stemp = new cplx [sg];
    element_set<T,SIZE1>(size1,one,1);
    // SET ENTRIES IN TIMESTEP TO 0
    for(i=0;i<=tc;i++) element_set_zero<T,SIZE1>(size1,G.retptr(0,i));
    // INITIAL VALUE t' = t
    element_set<T,SIZE1>(size1,G.retptr(0,0),Q.retptr(0,0));
    // START VALUES  t' = n-j, j = 1...k: solve a kxk problem
    for(i=0;i<k*k*sg;i++) mm[i]=0;
    for(i=0;i<k*sg;i++) qq[i]=0;
    for(j=1;j<=k;j++){
        p=j-1;
        element_incr<T,SIZE1>(size1,qq+p*sg,Q.retptr(0,j));
        element_incr<T,SIZE1>(size1,mm+sg*(p+k*p),one);
        // integral
        for(l=0;l<=k;l++){
            weight=-h*I.gregory_weights(j,l);
            if(j>=l){
                // ... F(t-l,t-j)=F(t-l,t-l-(j-l))
                element_set<T,SIZE1>(size1,stemp,F.retptr(l,j-l));
            }else{
                // ... Fcc(t-j,t-l)=F(t-j,t-j-(l-j))
                element_set<T,SIZE1>(size1,stemp,Fcc.retptr(j,l-j));
                element_conj<T,SIZE1>(size1,stemp);
                weight *= -1;
            }
            if(l==0){
                element_incr<T,SIZE1>(size1,qq+p*sg,weight,G.retptr(0,0),stemp);
            }else{
                q=l-1;
                element_incr<T,SIZE1>(size1,mm+sg*(p*k+q),-weight,stemp);
            }
        }
    }
    element_linsolve_left<T,SIZE1>(size1,k,gtemp,mm,qq); // gtemp * mm = qq 
    for(j=1;j<=k;j++) element_set<T,SIZE1>(size1,G.retptr(0,j),gtemp + (j-1)*sg);
    // REMAINING VALUES  t' = n-j, j = k+1,...,tc: solve a 1x1 problem
    for(j=k+1;j<=tc;j++){
        // qq <- Q(t,t-j)
        element_set<T,SIZE1>(size1,qq,Q.retptr(0,j));
        // qq += contribution from -\int_{t'}^t dt1 Fret(t,t1)Gret(t1,t')
        //       which does not contain Gret(t,t')
        element_set_zero<T,SIZE1>(size1,stemp);
        if(j<2*k+2){
            for(l=1;l<=j;l++){
                element_incr<T,SIZE1>(size1,stemp,I.gregory_weights(j,l),F.retptr(0,l),G.retptr(l,j-l));
            }
        }else{
            for(l=1;l<=k;l++){
                element_incr<T,SIZE1>(size1,stemp,I.gregory_omega(l),F.retptr(0,l),G.retptr(l,j-l));
            }
            for(l=k+1;l<j-k;l++){
                element_incr<T,SIZE1>(size1,stemp,F.retptr(0,l),G.retptr(l,j-l));
            }
            for(l=j-k;l<=j;l++){
                element_incr<T,SIZE1>(size1,stemp,I.gregory_omega(j-l),F.retptr(0,l),G.retptr(l,j-l));
            }
        }
        element_incr<T,SIZE1>(size1,qq,-h,stemp);
        // get prefactor of G on the other side
        for(i=0;i<sg;i++) mm[i] = one[i]  + h*I.gregory_omega(0)*F.retptr(0,0)[i];
        element_linsolve_right<T,SIZE1>(size1,G.retptr(0,j),mm,qq);
    }
    delete [] stemp;
    delete [] qq;
    delete [] mm;
    delete [] gtemp;
    delete [] one;
    return;
}

template < typename T, int SIZE1 >
void vie2_timestep_dispatch_les(moving_herm<T> &G,moving_herm<T> &F,moving_herm<T> &Fcc,moving_herm<T> &Q, integration::Integrator<T> &I, T h){
    typedef std::complex<T> cplx;
    int k=I.get_k();
    // various asserts are one level higher
    cplx *mm,*qq,*one,*gtemp,*stemp,weight;
    int i,j,l,dl,dl1;
    int size1=G.size1();
    int tc=G.tc();
    int sg=G.element_size();
    mm = new cplx [sg];
    qq = new cplx [sg];
    one = new cplx [sg];
    gtemp = new cplx [sg];
    stemp = new cplx [sg];
    element_set<T,SIZE1>(size1,one,1);
    for(i=0;i<=tc;i++) element_set_zero<T,SIZE1>(size1,G.lesptr(0,i));
    // DO Gles(t,t-j), j=tc,...,1
    element_set_zero<T,SIZE1>(size1,mm);
    element_set_zero<T,SIZE1>(size1,qq);
    // NOTE: THE SEQUENCE IS VERY IMPORTANT: THE LAST STEP (j=0) NEEDS j=tc,...,1 AS INPUT
    for(j=tc;j>=0;j--){
        // qq <- Q(t,t-j)
        element_set<T,SIZE1>(size1,qq,Q.lesptr(0,j));
        // qq += contribution from -\int_{t'}^t dt1 Fret(t,t1)Gles(t1,t') which does not contain Gles(t,t')
        element_set_zero<T,SIZE1>(size1,stemp);
        for(l=1;l<=tc;l++){
            if(j>l){
                 element_set<T,SIZE1>(size1,gtemp,G.lesptr(l,j-l));
            }else{
                 element_minusconj<T,SIZE1>(size1,gtemp,G.lesptr(j,l-j));
            }
            element_incr<T,SIZE1>(size1,stemp,I.gregory_weights(tc,l),F.retptr(0,l),gtemp);
        }
        element_incr<T,SIZE1>(size1,qq,-h,stemp);
        // qq += contribution from -\int_{t'}^t dt1 Fles(t,t1)Gadv(t1,t')
        element_set_zero<T,SIZE1>(size1,stemp);
        dl=tc-j;
        dl1=(dl>k ? dl : k);
        for(l=0;l<=dl1;l++){
            // Gadv(t-tc+l,t-j) = Gret(t-j,t-tc+l)^*
            if(tc-l>j){
                 element_conj<T,SIZE1>(size1,gtemp,G.retptr(j,tc-l-j));
            }else{
                 element_set<T,SIZE1>(size1,gtemp,G.retptr(tc-l,j-tc+l));
            }
            element_incr<T,SIZE1>(size1,stemp,I.gregory_weights(tc-j,l),F.lesptr(0,tc-l),gtemp);
        }
        element_incr<T,SIZE1>(size1,qq,-h,stemp);
        // get prefactor of G on the other side
        for(i=0;i<sg;i++) mm[i] = one[i]  + h*I.gregory_omega(0)*F.retptr(0,0)[i];
        element_linsolve_right<T,SIZE1>(size1,G.lesptr(0,j),mm,qq);
    }
    delete [] stemp;
    delete [] qq;
    delete [] mm;
    delete [] gtemp;
    delete [] one;
    return;
}

template < typename T>
void vie2_timestep(moving_herm<T> &G,moving_herm<T> &F,moving_herm<T> &Fcc,moving_herm<T> &Q, integration::Integrator<T> &I, T h){
    int kt=I.get_k();
    int size1=G.size1();
    CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,kt*2+2,G.tc(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,G.size1(),F.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,G.size1(),Fcc.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,G.size1(),Q.size1(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,G.tc(),F.tc(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,G.tc(),Fcc.tc(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,G.tc(),Q.tc(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,kt,F.nt(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,kt,Fcc.nt(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,0,Q.nt(),__PRETTY_FUNCTION__)
    CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,G.tc(),G.nt(),__PRETTY_FUNCTION__)
    if(size1==1){
        vie2_timestep_dispatch_ret<T,1>(G,F,Fcc,Q,I,h);
        vie2_timestep_dispatch_les<T,1>(G,F,Fcc,Q,I,h);
    }else{
        vie2_timestep_dispatch_ret<T,LARGESIZE>(G,F,Fcc,Q,I,h);
        vie2_timestep_dispatch_les<T,LARGESIZE>(G,F,Fcc,Q,I,h);
    }
}


template < typename T>
void vie2_timestep(moving_herm<T> &G,moving_herm<T> &F,moving_herm<T> &Fcc,moving_herm<T> &Q, int kt, T h){
    vie2_timestep(G,F,Fcc,Q,integration::I<T>(kt),h);
}



}