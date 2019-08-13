#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include "cntr.hpp"


#define CPLX std::complex<double>  
using namespace std;
#define GREEN cntr::herm_matrix<double>
#define GREENTSTP cntr::herm_matrix_timestep<double> 
#define FUNCTION cntr::function<double> 
// Access to the submatrix of the bigger herm_matrix Green's function
// Check on the block diagonal matrix where the submatrix of green's function
// is Green's function of submatrix

int main(int argc,char *argv[]){
	int nt,ntau,kt,tstp,size=5,subsize=2;
	double err,beta,h;
	GREEN A,B,B1,Asub,Asub1,tmp;
	nt=100;
	ntau=400;
	beta=10;
	h=0.02;
	// Test diagonal
	A=GREEN(nt,ntau,size,-1);
	tmp=GREEN(nt,ntau,size,-1);
	Asub=GREEN(nt,ntau,subsize,-1);
	B=GREEN(nt,ntau,subsize,-1);
	cdmatrix a(size,size);
	a.setZero();
	cdmatrix b(2,2);
	b.setZero();
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
	// Take elements of first submatrix
	std::vector<int> i1{0,0,1,1};
	std::vector<int> i2{0,1,0,1};
	std::vector<int> j1{0,0,1,1};
	std::vector<int> j2{0,1,0,1};
	Asub.set_submatrix(j1,j2,A,i1,i2);
	for (int i=0;i<i1.size();i++){
		b(j1[i],j2[i])=a(i1[i],i2[i]);
	}
	cntr::green_from_H(B,0.0,b,beta,h);
	
	std::cout << "Difference of greens function: " << std::endl;
	for(int tstp=-1;tstp<nt;tstp++){
		std::cout << tstp << " " << cntr::distance_norm2(tstp,B,Asub) << std::endl;
	}

	// Take elements of second submatrix
	Asub1=GREEN(nt,ntau,3,-1);
	B1=GREEN(nt,ntau,3,-1);
	cdmatrix b1(3,3);
	b1.setZero();
	std::vector<int> i11{2,2,2,3,3,3,4,4,4};
	std::vector<int> i12{2,3,4,2,3,4,2,3,4};
	std::vector<int> j11{0,0,0,1,1,1,2,2,2};
	std::vector<int> j12{0,1,2,0,1,2,0,1,2};
	Asub1.set_submatrix(j11,j12,A,i11,i12);
	for (int i=0;i<i11.size();i++){
		b1(j11[i],j12[i])=a(i11[i],i12[i]);
	}
	std::cout << "A: " << a << std::endl;
	std::cout << "B: " << b1 << std::endl;
	cntr::green_from_H(B1,0.0,b1,beta,h);
	

	std::cout << "Difference of greens function: " << std::endl;
	for(int tstp=-1;tstp<nt;tstp++){
		std::cout << tstp << " " << cntr::distance_norm2(tstp,B1,Asub1) << std::endl;
	}

	// Cfunction matrix element
	FUNCTION E,Esub;
	E=FUNCTION(nt,3);
	Esub=FUNCTION(nt,1);

	E.set_constant(b1);
	Esub.set_matrixelement(0,0,E,2,2);
	for(int tstp=-1;tstp<nt;tstp++){
		cdmatrix tmp,tmpsub;
		E.get_value(tstp,tmp);
		Esub.get_value(tstp,tmpsub);
		std::cout << tstp <<" " << tmp << std::endl;
		std::cout << std::endl;
		std::cout  << tmpsub  << std::endl;
	}

	// Multiplication of matrices and cfunctions with different sizes
	FUNCTION F=FUNCTION(nt,size,1);
	a.setZero();
	for (int i=0;i<size;i++){
		for (int j=0;j<size;j++){
			a(i,j)=i+j;
		}
	}
	cdmatrix offdiag=a.block(0,0,size,1).eval();
	F.set_constant(offdiag);
	{
		GREENTSTP G=GREENTSTP(-1,ntau,size,1,-1);
		GREENTSTP Atstp=GREENTSTP(-1,ntau,size,-1);
		// std::cout << "tu smo 2"<< std::endl;
		cdmatrix tmp,tmp1;
		A.get_timestep(-1,Atstp);
		A.get_mat(5,tmp);
		std::cout << "tu smo 2"<< std::endl;
		std::cout << "A green " << tmp << std::endl;
		std::cout << "tu smo 3 " << Atstp.size1() << " " << Atstp.size2() << " " << F.size1() << " " << F.size2()  << std::endl;
		G=Atstp.right_multiply_eigen(F,1.0);
		std::cout << "tu smo 4 "<< G.size1() << " " << G.size2() <<  std::endl;
		G.get_mat(5,tmp1);
		std::cout << "A*F green after " << tmp1 << std::endl;
	}
	return 0;
}




