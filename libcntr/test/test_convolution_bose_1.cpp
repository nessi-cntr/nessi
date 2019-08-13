#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include "cntr.hpp"

#define CPLX std::complex<double>  
using namespace std;
 

/*///////////////////////////////////////////////////////////////////////////////////////


This is a minimal test program to solve the convolution

A*B=C 

for bosonic Greens functions:

A= 1/(iwn-wa)
B= 1/(iwn-wb)
C= (A-B)/(wa-wb)

///////////////////////////////////////////////////////////////////////////////////////*/ 


int main(int argc,char *argv[]){
	int nt,ntau,kt,tstp;
	double wa,wb,h,beta,fac;
	cntr::herm_matrix<double> A,B,C,AB;

	nt=10;
	ntau=400;
	beta=0.1;
	h=0.02;
	wa=1.123;
	wb=0.345;    
	kt=5;

	A=cntr::herm_matrix<double>(nt,ntau,1,1);
	B=cntr::herm_matrix<double>(nt,ntau,1,1);
	C=cntr::herm_matrix<double>(nt,ntau,1,1);
	AB=cntr::herm_matrix<double>(nt,ntau,1,1);
	
	cntr::green_single_pole_bose(A,wa,beta,h);
	cntr::green_single_pole_bose(B,wb,beta,h);
	// C=(A-B)/(wa-wb)
	fac=1.0/(wa-wb);
	C=A;
	for(tstp=-1;tstp<=nt;tstp++){
		C.incr_timestep(tstp,B,CPLX(-1.0,0.0));
		C.smul(tstp,fac);
	}
	cout << "SIMPLE TEST WITH GIVEN h=0.02, dtau=0.005" << endl;
	cntr::convolution(AB,A,A,B,B,integration::I<double>(kt),beta,h);
	//cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
	for(tstp=-1;tstp<=nt;tstp++){
		double err;
		cout << "tstp: " << tstp;
		err=cntr::distance_norm2(tstp,C,AB);		
		cout << " |AB-C| " << err;
		err=cntr::distance_norm2_ret(tstp,C,AB);		
		cout << " |AB-C| ret " << err;
		err=cntr::distance_norm2_tv(tstp,C,AB);		
		cout << " |AB-C| tv " << err;
		err=cntr::distance_norm2_les(tstp,C,AB);		
		cout << " |AB-C| les " << err;
		cout << endl;
	}
	
	C.print_to_file("C.out");
	AB.print_to_file("AB.out");
	A.print_to_file("A.out");
	B.print_to_file("B.out");
	
	cout << "TESTING CONVERGENCE (SEE WHETHER BOSONIC SIGNS ARE WRONG IN BOUNDARY TERMS)" << endl;
	// TEST REAL-TIME ... SET MATSUBARA TO EXACT 
	{
		double tmax=2.0;
		int nt0=10,fac2;
		beta=0.1;
		for(fac2=1;fac2<=128;fac2*=2){
			nt=nt0*fac2;
			
			A=cntr::herm_matrix<double>(nt,ntau,1,1);
			B=cntr::herm_matrix<double>(nt,ntau,1,1);
			C=cntr::herm_matrix<double>(nt,ntau,1,1);
			AB=cntr::herm_matrix<double>(nt,ntau,1,1);
			
			h=tmax/nt;
			
			cntr::green_single_pole_bose(A,wa,beta,h);
			cntr::green_single_pole_bose(B,wb,beta,h);
			// C=(A-B)/(wa-wb)
			fac=1.0/(wa-wb);
			C=A;
			for(tstp=-1;tstp<=nt;tstp++){
				C.incr_timestep(tstp,B,CPLX(-1.0,0.0));
				C.smul(tstp,fac);
			}
			cntr::convolution(AB,A,A,B,B,integration::I<double>(kt),beta,h);
			//cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
			for(tstp=nt;tstp<=nt;tstp++){
				double err;
				cout << "tstp: " << tstp;
				err=cntr::distance_norm2(tstp,C,AB);		
				cout << " |AB-C| " << err;
				err=cntr::distance_norm2_ret(tstp,C,AB);		
				cout << " |AB-C| ret " << err;
				err=cntr::distance_norm2_tv(tstp,C,AB);		
				cout << " |AB-C| tv " << err;
				err=cntr::distance_norm2_les(tstp,C,AB);		
				cout << " |AB-C| les " << err;
				cout << endl;
			}
		}
	
	}
	return 0;
}




