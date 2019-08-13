#include "catch.hpp"
#include "cntr.hpp"
// #include <sys/stat.h>
// #include <iostream>
// #include <complex>
// #include <cmath>
// #include <cstring>



#define CPLX std::complex<double>  
#define GREEN cntr::herm_matrix<double>
#define GREENTSTP cntr::herm_matrix_timestep<double> 
#define FUNCTION cntr::function<double> 
// Check access to the submatrix of the bigger herm_matrix Green's function
// Check on the block diagonal matrix that the submatrix of green's functionis Green's function of submatrix


TEST_CASE("Herm matrix submatrix","[Herm_submatrix]"){
	int nt,ntau,kt,tstp,size=5,subsize=2;
	double err=0,beta,h,eps=1e-7;
	GREEN A,B,Asub,B1;
	nt=100;
	ntau=400;
	beta=10;
	h=0.02;

	// Test diagonal
	A=GREEN(nt,ntau,size,-1);
	Asub=GREEN(nt,ntau,subsize,-1);
	B=GREEN(nt,ntau,subsize,-1);
	B1=GREEN(nt,ntau,3,-1);

	cdmatrix a(size,size);
	a.setZero();
	cdmatrix b(2,2);
	
	for (int i=0;i<subsize;i++){
		for (int j=0;j<subsize;j++){
			a(i,j)=i+j;
		}
	}
	for (int i=size-3;i<size;i++){
		for (int j=size-3;j<size;j++){
			a(i,j)=i+j+1;
		}
	}
	cntr::green_from_H(A,0.0,a,beta,h);
	std::vector<int> i1{0,0,1,1};
	std::vector<int> i2{0,1,0,1};
	std::vector<int> j1{0,0,1,1};
	std::vector<int> j2{0,1,0,1};
	Asub.set_submatrix(j1,j2,A,i1,i2);


	cdmatrix b1(3,3);
	b1.setZero();
	std::vector<int> i11{2,2,2,3,3,3,4,4,4};
	std::vector<int> i12{2,3,4,2,3,4,2,3,4};
	std::vector<int> j11{0,0,0,1,1,1,2,2,2};
	std::vector<int> j12{0,1,2,0,1,2,0,1,2};
	for (int i=0;i<i11.size();i++){
		b1(j11[i],j12[i])=a(i11[i],i12[i]);
	}
	
	cntr::green_from_H(B1,0.0,b1,beta,h);

	SECTION ("Green of submatrix = Submatrix of Green"){
		// Take elements of first submatrix
		b.setZero();
		for (int i=0;i<i1.size();i++){
			b(j1[i],j2[i])=a(i1[i],i2[i]);
		}
		cntr::green_from_H(B,0.0,b,beta,h);
		for(int tstp=-1;tstp<nt;tstp++){
			err+=cntr::distance_norm2(tstp,B,Asub);
		}
		REQUIRE(err<eps);
	}
	
	SECTION ("Green of submatrix = Submatrix of Green, Test 2"){
		GREEN Asub1=GREEN(nt,ntau,3,-1);
		Asub1.set_submatrix(j11,j12,A,i11,i12);
		
		for(int tstp=-1;tstp<nt;tstp++){
			err+=cntr::distance_norm2(tstp,B1,Asub1);
		}
		REQUIRE(err<eps);
	}

}





