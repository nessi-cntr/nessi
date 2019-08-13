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


int main(int argc,char *argv[]){
	int nt,ntau,kt,tstp,size=5;
	double wa[size];
	std::complex<double> waa[size*size];
	double err,beta,h;
	
	GREEN A,Atest,A2test,A3test,B,C;
	
	nt=100;
	ntau=400;
	beta=10;
	h=0.02;

	// Test diagonal
	A=GREEN(nt,ntau,size,-1);
	Atest=GREEN(nt,ntau,size,-1);
	A2test=GREEN(nt,ntau,size,-1);
	A3test=GREEN(nt,ntau,size,-1);
	cdmatrix a(size,size);
	a.setZero();
	for (int i=0;i<size;i++){
		a(i,i)=i+1.0;
		wa[i]=i+1.0;
		waa[i+size*i]=i+1.0;
	}
	cntr::green_from_eps(Atest,0.0,wa,beta,h);
	cntr::green_from_H(A2test,0.0,waa,beta,h);
	cntr::green_from_H(A3test,0.0,a,beta,h);

	std::cout << "Difference of greens function: " << std::endl;
	for(int tstp=-1;tstp<nt;tstp++){
		std::cout << " " << tstp  << " " << cntr::distance_norm2(tstp,Atest,A2test) << " " << cntr::distance_norm2(tstp,A2test,A3test)  << std::endl;	 
	}
	nt=1000;
	size=2;
	// Test off-diagonal
	B=GREEN(nt,ntau,size,-1);
	cdmatrix b(size,size);
	b.setZero();
	b(0,0)=0.0;
	b(0,1)=0.5;
	b(1,0)=0.5;
	b(1,1)=1.0;

	cntr::green_from_H(B,0.0,b,beta,h);
	B.write_to_hdf5("B.h5","B");

	// Test green free from herm-matrix timestep
	// C=GREEN(nt,ntau,size,-1);
	// for(int tstp=-1;tstp<nt;tstp++){
	// 	GREEN_TSTP tmp=GREEN_TSTP(tstp,ntau,size);
	// 	green_from_eps(tmp,0.0,b,beta,h);
	// 	C.set_timestep(tstp,tmp);
	// }
	// for(int tstp=-1;tstp<nt;tstp++){
	// 	std::cout << " " << tstp  << " " << cntr::distance_norm2(tstp,B,C) << std::endl;	 
	// }
	return 0;
}




