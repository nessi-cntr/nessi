namespace cntr{



/*####################################################################################
#
#   Utilities for computing diagrams, such as bubbles, etc.
#
######################################################################################*/


/*/////////////////////////////////////////////////////////////////////////////////

#   BUBBLE1:  C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b2,b1}(t2,t1) 
#
#   BUBBLE2:  C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b1,b2}(t1,t2) 

//////////////////////////////////////////////////////////////////////////////////*/



// -------------  Auxiliary routines: -------------

// C(t-tstp,t-tstp-m)= ii * A(t-tstp,t-tstp-m) * B(t-tstp-m,t-tstp), l=0...tc
// where aret[m] -> A(t-tstp,t-tstp-m) etc.
template<typename T>
void get_bubble_1_timestep_moving_herm(int tc,std::complex<T> *cret,std::complex<T> *cles,int sc,int c1,int c2,std::complex<T> *aret,std::complex<T> *ales,std::complex<T> *accret,std::complex<T> *accles,int sa,int a1,int a2,std::complex<T> *bret,std::complex<T> *bles,std::complex<T> *bccret,std::complex<T> *bccles,int sb,int b1,int b2,int sigb){
	int m;
	int a12=a1*sa+a2,a21=a2*sa+a1,sa2=sa*sa;
	int b12=b1*sb+b2,b21=b2*sb+b1,sb2=sb*sb;
	int c12=c1*sc+c2,sc2=sc*sc;
	std::complex<T> msigb=-1.0*sigb,ii=std::complex<T>(0,1.0);
	std::complex<T> bgtr21_tt1,cgtr_tt1,cles_tt1,ales12_tt1,bgtr21_t1t,agtr12_tt1,bles21_t1t;
	for(m=0;m<=tc;m++){
        //Ales{12}(t,t-m)
        ales12_tt1 = ales[m*sa2+a12];
        //Agtr{12}(t,t-m) = Ales{12}(t,t-m) + Aret{12}(t,t-m)
        agtr12_tt1 = ales12_tt1 + aret[m*sa2+a12];
        //Bles{21}(t-m,t)
        bles21_t1t = -conj(bccles[m*sb2+b12]);
        //Bgtr{21}(t-m,t) =
        bgtr21_t1t = bles21_t1t-conj(bccret[m*sb2+b12]);
	cgtr_tt1 = ii*agtr12_tt1*bles21_t1t;
        cles_tt1 = ii*ales12_tt1*bgtr21_t1t;
	cret[m*sc2+c12] = cgtr_tt1-cles_tt1;
	cles[m*sc2+c12] = ii*ales12_tt1*bgtr21_t1t;
	}
}

//  BUBBLE 2 : 
//  C(t,t') = ii * A(t,t') * B(t,t') 

template<typename T>
void get_bubble_2_timestep_moving_herm(int tc,std::complex<T> *cret,std::complex<T> *cles,int sc,int c1,int c2,std::complex<T> *aret,std::complex<T> *ales,std::complex<T> *accret,std::complex<T> *accles,int sa,int a1,int a2,std::complex<T> *bret,std::complex<T> *bles,std::complex<T> *bccret,std::complex<T> *bccles,int sb,int b1,int b2){
	int m;
	int a12=a1*sa+a2,a21=a2*sa+a1,sa2=sa*sa;
	int b12=b1*sb+b2,b21=b2*sb+b1,sb2=sb*sb;
	int c12=c1*sc+c2,sc2=sc*sc;
	std::complex<T> ii=std::complex<T>(0,1.0);
	std::complex<T> bgtr12_tt1,agtr12_tt1,cgtr_tt1,cles_tt1;
	for(m=0;m<=tc;m++){
		bgtr12_tt1 = bret[m*sb2+b12] + bles[m*sb2+b12];
		agtr12_tt1 = aret[m*sa2+a12] + ales[m*sa2+a21];
		// Cgtr_{12}(tstp,m) = ii * Agtr_{12}(tstp,m)*Bles_{21}(m,tstp)
		cgtr_tt1 = ii * agtr12_tt1 * bgtr12_tt1;
		// Cles_{12}(tstp,m) = ii * Ales_{12}(tstp,m)*Bles_{12}(tstp,m)
		cles_tt1 = ii * ales[m*sa2+a21]*bles[m*sb2+b21];
		// Cret_{12}(tstp,m) = Cgtr_{12}(tstp,m) - Cles_{12}(tstp,m)
		cret[m*sc2+c12] = cgtr_tt1-cles_tt1;
		// Cles_{12}(m,tstp) = ii * Ales_{12}(m,tstp)*Bgtr_{12}(m,tstp)
		cles[m*sc2+c12] = ii* ales[m*sa2+a12] * bles[m*sb2+b12]; 
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
////// MAIN FUNCTIONS

// A,B not necesarily hermitian, orbital dimension > 1:
// C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b2,b1}(t2,t1) 


template<typename T>
void MovBubble1(int tstp, moving_herm<T>  &C,int c1,int c2, int tstpA,moving_herm<T>  &A,moving_herm<T>  &Acc,int a1,int a2, int tstpB,moving_herm<T>  &B,moving_herm<T>  &Bcc,int b1,int b2){
	int tc=C.tc();
    int nt=C.nt();
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc,C.tc(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc,B.tc(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc,A.tc(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc,Acc.tc(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc,Bcc.tc(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,tstp,nt,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,tstpA,A.nt(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,tstpB,B.nt(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,a1,A.size1()-1,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,a2,A.size2()-1,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,b1,B.size1()-1,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,b2,B.size2()-1,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,c1,C.size1()-1,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,c2,C.size2()-1,__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,Acc.size1(),A.size1(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,Bcc.size1(),B.size1(),__PRETTY_FUNCTION__)
    get_bubble_1_timestep_moving_herm(tc,C.retptr(tstp,0),C.lesptr(tstp,0),C.size1(),c1,c2,
    A.retptr(tstpA,0),A.lesptr(tstpA,0),Acc.retptr(tstpA,0),Acc.lesptr(tstpA,0),A.size1(),a1,a2,
    B.retptr(tstpB,0),B.lesptr(tstpB,0),Bcc.retptr(tstpB,0),Bcc.lesptr(tstpB,0),B.size1(),b1,b2,
    B.sig());
}


template<typename T>
void MovBubble1(moving_herm<T>  &C,int c1,int c2,moving_herm<T>  &A, moving_herm<T>  &Acc,int a1,int a2,moving_herm<T>  &B, moving_herm<T>  &Bcc,int b1,int b2){
    return MovBubble1(0,C,c1,c2,0,A,Acc,a1,a2,0,B,Bcc,b1,b2);
}
template<typename T>
void MovBubble1(moving_herm<T>  &C,int c1,int c2,moving_herm<T>  &A,int a1,int a2,moving_herm<T>  &B,int b1,int b2){
    return MovBubble1(0,C,c1,c2,0,A,A,a1,a2,0,B,B,b1,b2);
}
template<typename T>
void MovBubble1(moving_herm<T>  &C,moving_herm<T>  &A,moving_herm<T>  &B){
    return MovBubble1(0,C,0,0,0,A,A,0,0,0,B,B,0,0);
}
template<typename T>
void MovBubble1(moving_herm<T> &C,moving_herm<T>  &A, moving_herm<T>  &Acc,moving_herm<T>  &B, moving_herm<T>  &Bcc){
    return MovBubble1(0,C,0,0,0,A,Acc,0,0,0,B,Bcc,0,0);
}

// A,B hermitian, orbital dimension > 1:
// C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b1,b2}(t1,t2) 

template<typename T>
void MovBubble2(int tstp, moving_herm<T>  &C,int c1,int c2, int tstpA,moving_herm<T>  &A, moving_herm<T>  &Acc,int a1,int a2,
int tstpB,moving_herm<T>  &B, moving_herm<T>  &Bcc,int b1,int b2){
	int tc=C.tc();
    CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc,C.tc(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc,B.tc(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc,A.tc(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc,Acc.tc(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,tc,Bcc.tc(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,tstp,C.nt(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,tstpA,A.nt(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,tstpB,B.nt(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,a1,A.size1()-1,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,a2,A.size2()-1,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,b1,B.size1()-1,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,b2,B.size2()-1,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,c1,C.size1()-1,__PRETTY_FUNCTION__)
	CNTR_ASSERT_LESEQ_3(MOVING_HERM_ASSERT_LEVEL,0,c2,C.size2()-1,__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,Acc.size1(),A.size1(),__PRETTY_FUNCTION__)
	CNTR_ASSERT_EQ(MOVING_HERM_ASSERT_LEVEL,Bcc.size1(),B.size1(),__PRETTY_FUNCTION__)

    get_bubble_2_timestep_moving_herm(tc,C.retptr(tstp,0),C.lesptr(tstp,0),C.size1(),c1,c2,
    A.retptr(tstpA,0),A.lesptr(tstpA,0),Acc.retptr(tstpA,0),Acc.lesptr(tstpA,0),A.size1(),a1,a2,
    B.retptr(tstpB,0),B.lesptr(tstpB,0),Bcc.retptr(tstpB,0),Bcc.lesptr(tstpB,0),B.size1(),b1,b2);
}

template<typename T>
void MovBubble2(moving_herm<T>  &C,int c1,int c2,moving_herm<T>  &A, moving_herm<T>  &Acc,int a1,int a2,moving_herm<T>  &B, moving_herm<T>  &Bcc,int b1,int b2){
    return MovBubble2(0,C,c1,c2,0,A,Acc,a1,a2,0,B,Bcc,b1,b2);
}
template<typename T>
void MovBubble2(moving_herm<T>  &C,int c1,int c2,moving_herm<T>  &A,int a1,int a2,moving_herm<T>  &B,int b1,int b2){
    return MovBubble2(0,C,c1,c2,0,A,A,a1,a2,0,B,B,b1,b2);
}
template<typename T>
void MovBubble2(moving_herm<T>  &C,moving_herm<T>  &A,moving_herm<T>  &B){
    return MovBubble2(0,C,0,0,0,A,A,0,0,0,B,B,0,0);
}
template<typename T>
void MovBubble2(moving_herm<T> &C,moving_herm<T>  &A, moving_herm<T>  &Acc,moving_herm<T>  &B, moving_herm<T>  &Bcc){
    return MovBubble2(0,C,0,0,0,A,Acc,0,0,0,B,Bcc,0,0);
}

}

