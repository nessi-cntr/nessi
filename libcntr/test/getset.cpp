#include "catch.hpp"
#include "cntr.hpp"

#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>

TEST_CASE("Global get","[Global get]"){
	int nt=50;
	int ntau=500;
	int size=9;
	double eps=1e-7;
	double beta=5.0;
	double h=0.01;

	// scalar Green's function
	GREEN A_scalar=GREEN(nt,ntau,1,FERMION);
	cdmatrix a1x1(1,1);
	a1x1(0,0) = 1.0;
	cntr::green_from_H(A_scalar,0.0,a1x1,beta,h);

	// Square Green's function
	GREEN A_matrix=GREEN(nt,ntau,size,FERMION);

	cdmatrix a(size,size);
	a.setZero();
	// a(0,0)=sqrt(2.0);
	// a(0,1)=sqrt(2.0)*std::complex<double>(0.0,1.0);
	// a(1,0)=sqrt(2.0)*std::complex<double>(0.0,-1.0);
	// a(1,1)=-sqrt(2.0);
	for(int i=0; i<size; i++) a(i,i) = 0.01*(i-4);

	cntr::green_from_H(A_matrix,0.0,a,beta,h);

	SECTION ("Get for herm_matrix"){
		std::complex<double> les;
		cdmatrix les_matrix(size,size);
		double err=0.0;
		// for(int i = 0; i <= nt; i++){
		// 	for (int j = 0; j <= i; j++)
		// 	{
		// 		cntr::get_les(i, j, les, A_scalar);
		// 		std::cout << i << "  " << j << "  " << les << std::endl;
		// 	}
		// }

		for(int i = 0; i <= nt; i++){
			for (int j = 0; j <= nt; j++)
			{
				cntr::get_les(i, j, les_matrix, A_matrix);
				std::cout << i << "  " << j << "  " << les_matrix << std::endl;
				A_matrix.get_les(i, j, les_matrix);
				std::cout << i << "  " << j << "  " << les_matrix << std::endl;
			}
		}


	}

}
