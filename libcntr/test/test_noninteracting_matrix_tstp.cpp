#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include "cntr.hpp"

#define CPLX std::complex<double>  
using namespace std;
#define GREEN cntr::herm_matrix<double> 
#define GREEN_TSTP cntr::herm_matrix_timestep<double> 

// For a moment  only diagonal matrices are implemented

int main(int argc,char *argv[]){
	int nt,ntau,kt,tstp,size=5;
	double wa[size];
	std::complex<double> waa[size*size];
	double err,beta,h;
	
	GREEN A,B,C,D;
	
	nt=100;
	ntau=400;
	beta=10;
	h=0.02;

	// Test diagonal
	A=GREEN(nt,ntau,size,-1);
	cdmatrix a(size,size);
	a.setZero();
	for (int i=0;i<size;i++){
		a(i,i)=i+1.0;
		wa[i]=i+1.0;
	}
	cntr::green_from_H(A,0.0,a,beta,h);

	B=GREEN(nt,ntau,size,-1);
	for(int tstp=-1;tstp<nt;tstp++){
		GREEN_TSTP tmp=GREEN_TSTP(tstp,ntau,size);
		green_from_eps(tmp,0.0,wa,beta,h);
		B.set_timestep(tstp,tmp);
	}

	std::cout << "Difference of greens function: " << std::endl;
	for(int tstp=-1;tstp<nt;tstp++){
		std::cout << " " << tstp  << " " << cntr::distance_norm2(tstp,A,B)   << std::endl;	 
	}


	// Non-interacting with integration -- full green's function
	C=GREEN(nt,ntau,size,-1);
	cntr::function<double> ene(nt,size);
	ene.set_constant(a);
	cntr::green_from_eps(C,0.0,ene,beta,h,5);
	for(int tstp=-1;tstp<nt;tstp++){		
		std::cout << " " << tstp  << " " << cntr::distance_norm2(tstp,A,C) << std::endl;	 
	}
	// Non-interacting with integration -- time-slice
	D=GREEN(nt,ntau,size,-1);
	for(int tstp=-1;tstp<nt;tstp++){
		GREEN_TSTP tmp=GREEN_TSTP(tstp,ntau,size);
		green_from_eps(tmp,0.0,ene,beta,h,5);
		D.set_timestep(tstp,tmp);
	}
	for(int tstp=-1;tstp<nt;tstp++){
		std::cout << " " << tstp  << " " << cntr::distance_norm2(tstp,A,D) << std::endl;	 
	}

	return 0;
}




