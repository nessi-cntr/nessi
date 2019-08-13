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

solve A + A*B*A = C

for A, using some single-mode bosonic Greens functions:

A = phiA(t) * [1/(iwn-wa)] * phiA(t')^dagger
B = phiB(t) * [1/(iwn-wa)] * phiB(t')^dagger




where C is computed from the convolution above, and then
vie2 is used to solve

(1+F)*A=C
F=A*B
Fcc=B*A

///////////////////////////////////////////////////////////////////////////////////////*/ 


int main(int argc,char *argv[]){
	int nt,ntau,kt,tstp,nt0=10,fac2;
	double wa,wb,h,beta,tmax;
	nt=10;
	ntau=400;
	beta=0.1;
	h=0.02;
	wa=1.123;
	wb=0.345;    
	kt=5;
	tmax=2.0;
	beta=1.0;
	for(fac2=1;fac2<=128;fac2*=2){
		nt=nt0*fac2;
		cntr::herm_matrix<double> A(nt,ntau,1,1);
		cntr::herm_matrix<double> B(nt,ntau,1,1);
		cntr::herm_matrix<double> C(nt,ntau,1,1);
		cntr::herm_matrix<double> F(nt,ntau,1,1);
		cntr::herm_matrix<double> Fcc(nt,ntau,1,1);
		cntr::herm_matrix<double> A1(nt,ntau,1,1);
		cntr::function<double> fA(nt,1),fB(nt,1);
		h=tmax/nt;
		fA[-1]=CPLX(1.0,0.0);
		fB[-1]=CPLX(1.0,0.0);
		for(tstp=0;tstp<=nt;tstp++){
			double time=h*tstp;
			fA[tstp]=CPLX(cos(cos(0.77*time)),sin(cos(0.77*time)));
			fB[tstp]=CPLX(cos(cos(0.532*time)),sin(cos(0.532*time)));
		}
		cntr::green_single_pole_bose(A,wa,beta,h);
		cntr::green_single_pole_bose(B,wb,beta,h);
		for(tstp=-1;tstp<=nt;tstp++){
			A.left_multiply(tstp,fA,1.0);
			A.right_multiply_hermconj(tstp,fA,1.0);
			B.left_multiply_hermconj(tstp,fB,1.0);
			B.right_multiply(tstp,fB,1.0);
		}
		// C =  A + A*B*A
		// F   = A*B
		// Fcc = B*A
		cntr::convolution(F,A,A,B,B,integration::I<double>(kt),beta,h);
		cntr::convolution(Fcc,B,B,A,A,integration::I<double>(kt),beta,h);
		cntr::convolution(C,A,A,Fcc,F,integration::I<double>(kt),beta,h);
		for(tstp=-1;tstp<=nt;tstp++) C.incr_timestep(tstp,A,CPLX(1.0,0.0));
		// solve for (1+F)*A1=C
		cntr::vie2(A1,F,Fcc,C,integration::I<double>(kt),beta,h);
		for(tstp=-1;tstp<=nt;tstp++){
			double err;
			cout << "NT: " << nt << " tstp: " << tstp;
			err=cntr::distance_norm2(tstp,A,A1);		
			cout << " |A1-A| " << err;
			err=cntr::distance_norm2_ret(tstp,A,A1);		
			cout << " |A1-A| ret " << err;
			err=cntr::distance_norm2_tv(tstp,A,A1);		
			cout << " |A1-A| tv " << err;
			err=cntr::distance_norm2_les(tstp,A,A1);		
			cout << " |A1-A| les " << err;
			cout << endl;
		}
		
		//C.print_to_file("C.out");
		//AB.print_to_file("AB.out");
		//A.print_to_file("A.out");
		//B.print_to_file("B.out");
	}
	return 0;
}




