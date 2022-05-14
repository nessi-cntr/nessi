namespace cntr{



template < typename T, int SIZE1 >
void dyson_timestep_dispatch_ret(herm_matrix_moving<T> &G,herm_matrix_moving<T> &Sigma,function_moving<T> eps,T mu,integration::Integrator<T> &I, T h){
    typedef std::complex<T> cplx;
    int k=I.get_k(),k1=k+1;
    // various asserts are one level higher
    cplx *mm,*qq,*one,*gtemp,*stemp,weight,cplx_i,cweight,*diffw;
    int i,j,p,l,q;
    int size1=G.size1();
    int tc=G.tc();
    int sg=G.element_size();
    cplx_i=cplx(0,1);
    mm = new cplx [k*k*sg];
    qq = new cplx [k*sg];
    one = new cplx [sg];
    gtemp = new cplx [k*sg];
    stemp = new cplx [sg];
    diffw=new cplx [k1+1];
  	element_set<T,SIZE1>(size1,one,1);
    // SET ENTRIES IN TIMESTEP TO 0
    for(i=0;i<=tc;i++) element_set_zero<T,SIZE1>(size1,G.retptr(0,i));
    // INITIAL VALUE t' = t
    element_set<T,SIZE1>(size1,G.retptr(0,0),-cplx_i);
    // START VALUES  t' = n-j, j = 1...k: solve a kxk problem
    for(i=0;i<k*k*sg;i++) mm[i]=0;
    for(i=0;i<k*sg;i++) qq[i]=0;
    for(j=1;j<=k;j++){
        p=j-1;
        //derivatives:
        for(l=0;l<=k;l++){
            cweight=cplx_i/h*I.poly_differentiation(j,l);
            if(l==0){
                element_incr<T,SIZE1>(size1,qq + p*sg,-cweight,G.retptr(0,0));
            }else{
                q=l-1;
                element_incr<T,SIZE1>(size1,mm + sg*(p*k+q),cweight);
            }
        }
        // H
        element_set<T,SIZE1>(size1,gtemp,eps.ptr(j));
        element_smul<T,SIZE1>(size1,gtemp,-1);
        for(i=0;i<sg;i++) gtemp[i] += mu*one[i];
        element_incr<T,SIZE1>(size1,mm+sg*(p+k*p),gtemp);
        // integral
        for(l=0;l<=k;l++){
            weight=h*I.gregory_weights(j,l);
            if(j>=l){
                // ... F(t-l,t-j)=F(t-l,t-l-(j-l))
                element_set<T,SIZE1>(size1,stemp,Sigma.retptr(l,j-l));
            }else{
                // ... Fcc(t-j,t-l)=F(t-j,t-j-(l-j))
                element_set<T,SIZE1>(size1,stemp,Sigma.retptr(j,l-j));
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
    for(p=0;p<=k+1;p++) diffw[p]=I.bd_weights(p)*cplx_i/h; // use BD(k+1!!)
    for(j=k+1;j<=tc;j++){
        // qq <- Q(t,t-j)
        element_set_zero<T,SIZE1>(size1,qq);
        // qq += contribution from -\int_{t'}^t dt1 Fret(t,t1)Gret(t1,t')
        //       which does not contain Gret(t,t')
        element_set_zero<T,SIZE1>(size1,stemp);
        if(j<2*k+2){
            for(l=1;l<=j;l++){
                element_incr<T,SIZE1>(size1,stemp,I.gregory_weights(j,l),Sigma.retptr(0,l),G.retptr(l,j-l));
            }
        }else{
            for(l=1;l<=k;l++){
                element_incr<T,SIZE1>(size1,stemp,I.gregory_omega(l),Sigma.retptr(0,l),G.retptr(l,j-l));
            }
            for(l=k+1;l<j-k;l++){
                element_incr<T,SIZE1>(size1,stemp,Sigma.retptr(0,l),G.retptr(l,j-l));
            }
            for(l=j-k;l<=j;l++){
                element_incr<T,SIZE1>(size1,stemp,I.gregory_omega(j-l),Sigma.retptr(0,l),G.retptr(l,j-l));
            }
        }
        element_incr<T,SIZE1>(size1,qq,h,stemp);
        // contribution from derivative
        for(p=1;p<=k+1;p++) element_incr<T,SIZE1>(size1,qq,-diffw[p],G.retptr(p,j-p)); // G(n-p,n-j)=(n-p,n-p-(j-p))
        element_set<T,SIZE1>(size1,mm,diffw[0]+mu);
        element_incr<T,SIZE1>(size1,mm,-h*I.gregory_omega(0),Sigma.retptr(0,0));
        element_incr<T,SIZE1>(size1,mm,cplx(-1.0,0.0),eps.ptr(0));
        element_linsolve_right<T,SIZE1>(size1,G.retptr(0,j),mm,qq);
    }
    delete [] stemp;
    delete [] qq;
    delete [] mm;
    delete [] gtemp;
    delete [] one;
    delete [] diffw;
    return;
}

template < typename T, int SIZE1 >
void dyson_timestep_dispatch_les(herm_matrix_moving<T> &G,herm_matrix_moving<T> &Sigma,function_moving<T> &eps,T mu, integration::Integrator<T> &I, T h){
    typedef std::complex<T> cplx;
    int k=I.get_k(),k1=k+1;;
    // various asserts are one level higher
    cplx *mm,*qq,*one,*gtemp,*stemp,weight,*diffw;
    cplx cplx_i=cplx(0,1);
    int i,j,l,dl,dl1,p;
    int size1=G.size1();
    int tc=G.tc();
    int sg=G.element_size();
    mm = new cplx [sg];
    qq = new cplx [sg];
    one = new cplx [sg];
    gtemp = new cplx [sg];
    stemp = new cplx [sg];
    diffw=new cplx [k1+1];
    element_set<T,SIZE1>(size1,one,1);
    for(p=0;p<=k+1;p++) diffw[p]=I.bd_weights(p)*cplx_i/h; // use BD(k+1!!)
    for(i=0;i<=tc;i++) element_set_zero<T,SIZE1>(size1,G.lesptr(0,i));
    // DO Gles(t,t-j), j=tc,...,1
    element_set_zero<T,SIZE1>(size1,mm);
    element_set_zero<T,SIZE1>(size1,qq);
    // NOTE: THE SEQUENCE IS VERY IMPORTANT: THE LAST STEP (j=0) NEEDS j=tc,...,1 AS INPUT
    for(j=tc;j>k;j--){
        element_set_zero<T,SIZE1>(size1,qq);
        // qq += contribution from -\int_{t'}^t dt1 Sigma^ret(t,t-l)Gles(t-l,t-j) which does not contain Gles(t,t')
        //if(j==20){
          //  for(int i=0;i<=tc;i++){
            //    for(int m=0;m<=tc;m++){
              //      //Sigma.retptr(i,m)[0]=0.0;
                //    G.retptr(i,m)[0]=0.0;
                  //  //G.lesptr(i,m)[0]=0.0;
                    //Sigma.lesptr(i,m)[0]=0.0;
                //}
            //}
        //}
        
        element_set_zero<T,SIZE1>(size1,stemp);
        for(l=1;l<=tc;l++){
            if(j>=l){
                 element_set<T,SIZE1>(size1,gtemp,G.lesptr(l,j-l));
            }else{
                 element_minusconj<T,SIZE1>(size1,gtemp,G.lesptr(j,l-j));
            }
            element_incr<T,SIZE1>(size1,stemp,I.gregory_weights(tc,l),Sigma.retptr(0,l),gtemp);
        }
        element_incr<T,SIZE1>(size1,qq,h,stemp);
        // qq += contribution from -\int_{t'}^t dt1 Fles(t,t1)Gadv(t1,t')
        element_set_zero<T,SIZE1>(size1,stemp);
        dl=tc-j;
        dl1=(dl>k ? dl : k);
        for(l=0;l<=dl1;l++){
            // Gadv(t-tc+l,t-j) = Gret(t-j,t-tc+l)^*
            if(tc-l>=j){
                 element_conj<T,SIZE1>(size1,gtemp,G.retptr(j,tc-l-j));
            }else{
                 element_set<T,SIZE1>(size1,gtemp,G.retptr(tc-l,j-tc+l));
		 element_smul<T,SIZE1>(size1,gtemp,std::complex<double>(-1.0,0));//added by cstahl
            }
            element_incr<T,SIZE1>(size1,stemp,I.gregory_weights(tc-j,l),Sigma.lesptr(0,tc-l),gtemp);
        }
        element_incr<T,SIZE1>(size1,qq,h,stemp);
        /*if(j==20) std::cout << "TTT XXX " <<  Sigma.lesptr(0,j)[0].real() << " " << Sigma.lesptr(0,j)[0].imag() <<  std::endl;
        if(j==20) std::cout << "TTT XXX " <<  qq[0].real() << " " << qq[0].imag() <<  std::endl;*/
       
         // contribution from derivative:
        // contribution from derivative
        for(p=1;p<=k+1;p++){
            // G(n-p,n-j)=(n-p,n-p-(j-p))
            if(j>=p){
                element_set<T,SIZE1>(size1,gtemp,G.lesptr(p,j-p));
            }else{
                element_minusconj<T,SIZE1>(size1,gtemp,G.lesptr(j,p-j));
            }
            element_incr<T,SIZE1>(size1,qq,-diffw[p],gtemp); // G(n-p,n-j)=(n-p,n-p-(j-p))
        }
        element_set<T,SIZE1>(size1,mm,diffw[0]+mu);
        element_incr<T,SIZE1>(size1,mm,-h*I.gregory_omega(0),Sigma.retptr(0,0));
        element_incr<T,SIZE1>(size1,mm,cplx(-1.0,0.0),eps.ptr(0));
        element_linsolve_right<T,SIZE1>(size1,G.lesptr(0,j),mm,qq);
        /*if(1){
            std::cout.precision(10);
            std::cout << "j " << 30-j;
            std::cout << " G " << G.lesptr(0,j)[0].real() << " " << G.lesptr(0,j)[0].imag();
            std::cout << " mm " << mm[0].real() << " " << mm[0].imag();
            std::cout << " qq " << qq[0].real() << " " << qq[0].imag();
            std::cout<< std::endl;
        }*/
    }
    //// steps kt ... 0 from d/dt' equation:
    for(j=k;j>=0;j--){
        element_set_zero<T,SIZE1>(size1,qq);
        // qq += contribution from \int_{t'}^t dt1 Gret(t,t-l)Sigmales(t-l,t-j) which does not contain Gles(t,t')
        element_set_zero<T,SIZE1>(size1,stemp);
        for(l=0;l<=tc;l++){
            if(j>l){
                 element_set<T,SIZE1>(size1,gtemp,Sigma.lesptr(l,j-l));
            }else{
                 element_minusconj<T,SIZE1>(size1,gtemp,Sigma.lesptr(j,l-j));
            }
            element_incr<T,SIZE1>(size1,stemp,I.gregory_weights(tc,l),G.retptr(0,l),gtemp);
        }
        element_incr<T,SIZE1>(size1,qq,h,stemp);
        // qq += contribution from \int^{t-j} dt1 Gles(t,t-l)Sigma^adv(t-l,t-j)
        element_set_zero<T,SIZE1>(size1,stemp);
        for(l=tc;l>j;l--){
            // Gadv(t-tc+l,t-j) = Gret(t-j,t-tc+l)^*
            element_conj<T,SIZE1>(size1,gtemp,Sigma.retptr(j,l-j));
            element_incr<T,SIZE1>(size1,stemp,I.gregory_weights(tc-j,l-j),G.lesptr(0,l),gtemp);
        }
        element_incr<T,SIZE1>(size1,qq,h,stemp);
        // contribution from derivative -i d/dt' = -i/h sum_p a_p G(t,t-j-p)
        for(p=1;p<=k+1;p++){
            // G(n-p,n-j)=(n-p,n-p-(j-p))
            element_set<T,SIZE1>(size1,gtemp,G.lesptr(0,j+p));
            element_incr<T,SIZE1>(size1,qq,diffw[p],gtemp); // G(n-p,n-j)=(n-p,n-p-(j-p))
        }
        element_set<T,SIZE1>(size1,mm,-diffw[0]+mu);
        element_conj<T,SIZE1>(size1,gtemp,Sigma.retptr(j,0));
        element_incr<T,SIZE1>(size1,mm,-h*I.gregory_omega(0),gtemp);
        element_conj<T,SIZE1>(size1,gtemp,eps.ptr(j));
        element_incr<T,SIZE1>(size1,mm,cplx(-1.0,0.0),gtemp);
        element_linsolve_left<T,SIZE1>(size1,G.lesptr(0,j),mm,qq);
    }
    //
    
    delete [] stemp;
    delete [] qq;
    delete [] mm;
    delete [] gtemp;
    delete [] one;
    delete [] diffw;
    return;
}

#define MOVING_HERM_ASSERT_LEVEL 1
template < typename T>
void dyson_timestep(herm_matrix_moving<T> &G,herm_matrix_moving<T> &Sigma,function_moving<T> &eps,T mu, integration::Integrator<T> &I, T h){
    int kt=I.get_k();
    int size1=G.size1();
    // CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,kt*2+2,G.tc(),__PRETTY_FUNCTION__);
    // CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,G.size1(),Sigma.size1(),__PRETTY_FUNCTION__);
    // CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,G.size1(),eps.size1(),__PRETTY_FUNCTION__);
    // CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,G.tc(),Sigma.tc(),__PRETTY_FUNCTION__);
    // CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,kt,Sigma.tc(),__PRETTY_FUNCTION__);
    //CNTR_ASSERT_LESEQ(MOVING_HERM_ASSERT_LEVEL,G.tc(),G.nt(),__PRETTY_FUNCTION__);
    // CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,G.tc(),eps.tc(),__PRETTY_FUNCTION__);
    if(size1==1){
        dyson_timestep_dispatch_ret<T,1>(G,Sigma,eps,mu,I,h);
        dyson_timestep_dispatch_les<T,1>(G,Sigma,eps,mu,I,h);
    }else{
        dyson_timestep_dispatch_ret<T,LARGESIZE>(G,Sigma,eps,mu,I,h);
        dyson_timestep_dispatch_les<T,LARGESIZE>(G,Sigma,eps,mu,I,h);
    }
}

template < typename T>
void dyson_timestep(herm_matrix_moving<T> &G,herm_matrix_moving<T> &Sigma,function_moving<T> &eps,T mu, int kt, T h){
    dyson_timestep(G,Sigma,eps,mu,integration::I<T>(kt),h);
}




}
