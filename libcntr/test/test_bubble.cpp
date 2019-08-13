#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include "cntr.hpp"

#define CPLX std::complex<double>  
using namespace std;
 
 
#define GREEN cntr::herm_matrix<double> 

/*///////////////////////////////////////////////////////////////////////////////////////

Test of bubble, comparison with Routines written by Sharareh

///////////////////////////////////////////////////////////////////////////////////////*/ 
void get_bubble_FF(int tstp, GREEN &Cbose, GREEN &Afermi, GREEN &Bfermi){
	int m;
	int nt,ntau;
	CPLX tempX,tempf,temploc,tempn;
	CPLX Icplx;
	Icplx.real(0.0);
	Icplx.imag(1.0);
	
	ntau=Cbose.ntau();
	nt=Cbose.nt();
	assert( tstp<=nt && tstp>=-1 );
	assert( Afermi.nt()==nt );
	assert( Bfermi.nt()==nt );
	assert( Afermi.ntau()==ntau );
	assert( Bfermi.ntau()==ntau );
	assert( Afermi.sig()==Bfermi.sig());
	assert( Afermi.sig()==-Cbose.sig());
  
	//CF =  Af* Bf 
	if(tstp==-1){
		//Matsubara part: sigma(-i\tau)= mat* mat
		for(m=0;m<=ntau;m++){
		  
			//Bfermi.get_matminus(m,tempf);//CHECK
			//Afermi.get_mat(m,temploc);
			Bfermi.get_mat(m,tempf);//CHECK
			Afermi.get_matminus(m,temploc);
			
			tempn= - temploc * tempf ;
			Cbose.set_mat(m,tempn);
		}
	}else{
		// retarded part: sigma(tstp,m)
		
		for(m=0;m<=tstp;m++){
		
			Bfermi.get_ret(tstp,m,tempf);
			Afermi.get_les(m,tstp,temploc);
			
			tempn= tempf * temploc  ;
			
			// Bgtr(t,t1)*Ales(t1,t) - Bles(t,t1)*Ales(t1,t) 
			
			Bfermi.get_les(tstp,m,tempf);
			Afermi.get_ret(m,tstp,temploc);
			temploc *= -1;
			
			// Bgtr(t,t1)*Ales(t1,t) - Bles(t,t1)*Ales(t1,t) - Bles(t,t1)*Agtr(t1,t) +  Bles(t,t1)*Ales(t1,t) 
			// Bgtr(t,t1)*Ales(t1,t) - Bles(t,t1)*Agtr(t1,t) - Bles(t,t1)*Ales(t1,t)  +  Bles(t,t1)*Ales(t1,t) 
			
			tempn += tempf * temploc  ;
			tempn *= Icplx;
			Cbose.set_ret(tstp,m,tempn);
		}
		// lesser part: sigma(m,tstp)
		for(m=0;m<=tstp;m++){
		  
			Bfermi.get_les(m,tstp,tempf);
			Afermi.get_gtr(tstp,m,temploc);
			
			// tempn= tempf * temploc;
			tempn= Icplx * tempf * temploc;
			Cbose.set_les(m,tstp,tempn);
		}
		// tv part: 
		for(m=0;m<=ntau;m++){
		  
			Bfermi.get_tv(tstp,m,tempf);
			Afermi.get_vt(m,tstp,temploc);
			
			// tempn= tempf * temploc;
			tempn= Icplx * tempf * temploc;
			Cbose.set_tv(tstp,m,tempn);
		}
	}
}
void get_bubble_BF(int tstp,GREEN &Cfermi,GREEN &Afermi,GREEN &Bbose){
	int m;
	int nt,ntau;
	CPLX tempX,tempf,temploc,tempn;
	CPLX Icplx;
	Icplx.real(0.0);
	Icplx.imag(1.0);
	
	ntau=Cfermi.ntau();
	nt=Cfermi.nt();
	//assert( tstp<=nt && tstp>=-1 );
	//assert( Afermi.nt()==nt );
	//assert( Bbose.nt()==nt );
	//assert( Afermi.ntau()==ntau );
	//assert( Bbose.ntau()==ntau );
	//assert( Afermi.sig()==Cfermi.sig());
	//assert( Bbose.sig()==-Cfermi.sig());
	
	
	//CF =  Af* Bb 
	if(tstp==-1){
		//Matsubara part: C(-i\tau)
		for(m=0;m<=ntau;m++){
			
			Bbose.get_mat(m,tempX);
			Afermi.get_mat(m,temploc);
			
			tempn= tempX * temploc  ;
			Cfermi.set_mat(m,tempn);
			
		}
	}else{
		// retarded part: C(tstp,m)
		
		for(m=0;m<=tstp;m++){

			
			Bbose.get_gtr(tstp,m,tempX);
			Afermi.get_gtr(tstp,m,temploc);
			
			tempn= tempX * temploc  ;
			
			Bbose.get_les(tstp,m,tempX);
			Afermi.get_les(tstp,m,temploc);
			
			tempn -= tempX * temploc  ;
			Cfermi.set_ret(tstp,m,Icplx*tempn);
			
			

		}
		// lesser part: C(m,tstp)
		for(m=0;m<=tstp;m++){
		  
			Bbose.get_les(m,tstp,tempX);
			Afermi.get_les(m,tstp,temploc);
			
			tempn= tempX * temploc;
			Cfermi.set_les(m,tstp,Icplx*tempn);
			

		}
		// tv part: 
		for(m=0;m<=ntau;m++){
		  
			Bbose.get_tv(tstp,m,tempX);
			Afermi.get_tv(tstp,m,temploc);
			
			tempn= tempX * temploc;
			Cfermi.set_tv(tstp,m,Icplx*tempn);
			
		}
	}
}
void get_bubble_BB(int tstp,GREEN &Cbose,GREEN &Abose,GREEN &Bbose){
	int m;
	int nt,ntau;
	CPLX tempX,tempf,temploc,tempn;
	CPLX Icplx;
	Icplx.real(0.0);
	Icplx.imag(1.0);
	
	ntau=Cbose.ntau();
	nt=Cbose.nt();
	//assert( tstp<=nt && tstp>=-1 );
	//assert( Abose.nt()==nt );
	//assert( Bbose.nt()==nt );
	//assert( Abose.ntau()==ntau );
	//assert( Bbose.ntau()==ntau );
	///assert( Abose.sig()==Bbose.sig());
	//assert( Bbose.sig()==Cbose.sig());
	
	
		//Cb =  Ab* Bb 
	if(tstp==-1){
		//Matsubara part: sigma(-i\tau)= mat* mat
		for(m=0;m<=ntau;m++){
		  
			Bbose.get_mat(m,tempf);//CHECK
			Abose.get_matminus(m,temploc);
			
			tempn= - temploc * tempf ;
			Cbose.set_mat(m,tempn);
		}
	}else{
		// retarded part: sigma(tstp,m)
		
		for(m=0;m<=tstp;m++){
		
			Bbose.get_ret(tstp,m,tempf);
			Abose.get_les(m,tstp,temploc);
			
			tempn= tempf * temploc  ;
			
			Bbose.get_les(tstp,m,tempf);
			Abose.get_ret(m,tstp,temploc);
			temploc *= -1;
			
			tempn += tempf * temploc  ;
			Cbose.set_ret(tstp,m,Icplx*tempn);
		}
		// lesser part: sigma(m,tstp)
		for(m=0;m<=tstp;m++){
		  
			Bbose.get_les(m,tstp,tempf);
			Abose.get_gtr(tstp,m,temploc);
			
			tempn= tempf * temploc;
			Cbose.set_les(m,tstp,Icplx*tempn);
		}
		// tv part: 
		for(m=0;m<=ntau;m++){
		  
			Bbose.get_tv(tstp,m,tempf);
			Abose.get_vt(m,tstp,temploc);
			
			tempn= tempf * temploc;
			Cbose.set_tv(tstp,m,Icplx*tempn);
		}
	}
}






int main(int argc,char *argv[]){
	int nt,ntau,kt,tstp;
	double wa,wb,h,beta,fac;
	double err;
	
	GREEN A,B,C,C1,A2,B2,A2cc,B2cc,A1,A1cc,B1,B1cc;
	
	nt=10;
	ntau=400;
	beta=0.1;
	h=0.02;
	wa=1.123;
	wb=0.345;    
	kt=5;

	A=GREEN(nt,ntau,1,-1);
	B=GREEN(nt,ntau,1,-1);
	C=GREEN(nt,ntau,1,1);
	C1=GREEN(nt,ntau,1,1);
	
	cntr::green_from_eps(A,0.0,&wa,beta,h);
	cntr::green_from_eps(B,0.0,&wb,beta,h);
	
	for(tstp=-1;tstp<=nt;tstp++){
		get_bubble_BB(tstp,C1,B,A);
		cntr::Bubble1(tstp,C,A,B);
		cout << "BUBBLE 1 tstp: " << tstp;
		err=cntr::distance_norm2(tstp,C,C1);		
		cout << " |C1-C| " << err;
		err=cntr::distance_norm2_ret(tstp,C,C1);		
		cout << " |A1-A| ret " << err;
		err=cntr::distance_norm2_tv(tstp,C,C1);		
		cout << " |A1-A| tv " << err;
		err=cntr::distance_norm2_les(tstp,C,C1);		
		cout << " |A1-A| les " << err;
		cout << endl;
	}
	for(tstp=-1;tstp<=nt;tstp++){
		get_bubble_BF(tstp,C1,B,A);
		cntr::Bubble2(tstp,C,A,B);
		cout << "BUBBLE 2 tstp: " << tstp;
		err=cntr::distance_norm2(tstp,C,C1);		
		cout << " |C1-C| " << err;
		err=cntr::distance_norm2_ret(tstp,C,C1);		
		cout << " |A1-A| ret " << err;
		err=cntr::distance_norm2_tv(tstp,C,C1);		
		cout << " |A1-A| tv " << err;
		err=cntr::distance_norm2_les(tstp,C,C1);		
		cout << " |A1-A| les " << err;
		cout << endl;
	}
	
	////////////////////////////////////////////////////////////////////
	// BOSONIC GREEN FUNCTIONS
	A=GREEN(nt,ntau,1,1);
	B=GREEN(nt,ntau,1,1);
	C=GREEN(nt,ntau,1,1);
	C1=GREEN(nt,ntau,1,1);
	cntr::green_single_pole_bose(A,&wa,beta,h);
	cntr::green_single_pole_bose(B,&wb,beta,h);
	
	for(tstp=-1;tstp<=nt;tstp++){
		get_bubble_BB(tstp,C1,B,A);
		cntr::Bubble1(tstp,C,A,B);
		cout << "BUBBLE BOSE 1 tstp: " << tstp;
		err=cntr::distance_norm2(tstp,C,C1);		
		cout << " |C1-C| " << err;
		err=cntr::distance_norm2_ret(tstp,C,C1);		
		cout << " |A1-A| ret " << err;
		err=cntr::distance_norm2_tv(tstp,C,C1);		
		cout << " |A1-A| tv " << err;
		err=cntr::distance_norm2_les(tstp,C,C1);		
		cout << " |A1-A| les " << err;
		cout << endl;
	}
	for(tstp=-1;tstp<=nt;tstp++){
		get_bubble_BF(tstp,C1,B,A);
		cntr::Bubble2(tstp,C,A,B);
		cout << "BUBBLE BOSE 2 tstp: " << tstp;
		err=cntr::distance_norm2(tstp,C,C1);		
		cout << " |C1-C| " << err;
		err=cntr::distance_norm2_ret(tstp,C,C1);		
		cout << " |A1-A| ret " << err;
		err=cntr::distance_norm2_tv(tstp,C,C1);		
		cout << " |A1-A| tv " << err;
		err=cntr::distance_norm2_les(tstp,C,C1);		
		cout << " |A1-A| les " << err;
		cout << endl;
	}

	////////////////////////////////////////////////////////////////////
	// SIMPLE TEST OF MATRIX VALUED BUBBLES
	A2=GREEN(nt,ntau,2,-1);
	A2cc=GREEN(nt,ntau,2,-1);
	B2=GREEN(nt,ntau,2,-1);
	B2cc=GREEN(nt,ntau,2,-1);
	A1=GREEN(nt,ntau,1,-1);
	A1cc=GREEN(nt,ntau,1,-1);
	B1=GREEN(nt,ntau,1,-1);
	B1cc=GREEN(nt,ntau,1,-1);
	
	cntr::green_from_eps(A1,0.0,&wa,beta,h);
	for(tstp=-1;tstp<=nt;tstp++){
	    A2.set_matrixelement(tstp,0,1,A1,0,0);
	    A2cc.set_matrixelement(tstp,1,0,A1,0,0);
	}
	wa=0.24679;
	cntr::green_from_eps(B1,0.0,&wa,beta,h);
	for(tstp=-1;tstp<=nt;tstp++){
	    B2.set_matrixelement(tstp,1,0,B1,0,0);
	    B2cc.set_matrixelement(tstp,0,1,B1,0,0);
	}
	for(tstp=-1;tstp<=nt;tstp++){
		cntr::Bubble1(tstp,C,0,0,A2,A2cc,0,1,B2,B2cc,0,1);
		cntr::Bubble1(tstp,C1,A1,B1);
		cout << "BUBBLE 1 tstp: " << tstp;
		err=cntr::distance_norm2(tstp,C,C1);		
		cout << " |C1-C| " << err;
		err=cntr::distance_norm2_ret(tstp,C,C1);		
		cout << " |A1-A| ret " << err;
		err=cntr::distance_norm2_tv(tstp,C,C1);		
		cout << " |A1-A| tv " << err;
		err=cntr::distance_norm2_les(tstp,C,C1);		
		cout << " |A1-A| les " << err;
		cout << endl;
	}

	for(tstp=-1;tstp<=nt;tstp++){
		cntr::Bubble2(tstp,C,0,0,A2,A2cc,0,1,B2,B2cc,1,0);
		cntr::Bubble2(tstp,C1,A1,B1);
		cout << "BUBBLE 2 tstp: " << tstp;
		err=cntr::distance_norm2(tstp,C,C1);		
		cout << " |C1-C| " << err;
		err=cntr::distance_norm2_ret(tstp,C,C1);		
		cout << " |A1-A| ret " << err;
		err=cntr::distance_norm2_tv(tstp,C,C1);		
		cout << " |A1-A| tv " << err;
		err=cntr::distance_norm2_les(tstp,C,C1);		
		cout << " |A1-A| les " << err;
		cout << endl;
	}


	



	return 0;
}




