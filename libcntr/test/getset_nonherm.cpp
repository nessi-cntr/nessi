#include "catch.hpp"
#include "cntr.hpp"

#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define GREEN_TSTP_VIEW cntr::herm_matrix_timestep_view<double>
#define CINTEG integration::I<double>

TEST_CASE("cntr::get_XXX for nonhermitian functions","[cntr::get_XXX for nonhermitian functions]"){
	int nt=50;
	int ntau=500;
	int size=2;
	double eps=1e-7;
	double beta=5.0;
	double h=0.01;

	// scalar Green's function
	GREEN A_scalar=GREEN(nt,ntau,1,FERMION);
	GREEN B_scalar=GREEN(nt,ntau,1,FERMION);
	GREEN C_scalar=GREEN(nt,ntau,1,FERMION);
	GREEN Ccc_scalar=GREEN(nt,ntau,1,FERMION);
	cdmatrix a1x1(1,1);
	cdmatrix b1x1(1,1);
	a1x1(0,0) = 1.0;
	b1x1(0,0) = -0.3;
	cntr::green_from_H(A_scalar,0.0,a1x1,beta,h);
	cntr::green_from_H(B_scalar,0.0,b1x1,beta,h);

	cntr::convolution(C_scalar, A_scalar, A_scalar, B_scalar, B_scalar, CINTEG(MAX_SOLVE_ORDER), beta, h);
	cntr::convolution(Ccc_scalar, B_scalar, B_scalar, A_scalar, A_scalar, CINTEG(MAX_SOLVE_ORDER), beta, h);

	// Square Green's function
	GREEN A_matrix=GREEN(nt,ntau,size,FERMION);
	GREEN B_matrix=GREEN(nt,ntau,size,FERMION);
	GREEN C_matrix=GREEN(nt,ntau,size,FERMION);
	GREEN Ccc_matrix=GREEN(nt,ntau,size,FERMION);

	cdmatrix a(size,size);
	a.setZero();
	a(0,0)=sqrt(2.0);
	a(0,1)=sqrt(2.0)*std::complex<double>(0.0,1.0);
	a(1,0)=sqrt(2.0)*std::complex<double>(0.0,-1.0);
	a(1,1)=-sqrt(2.0);

	cdmatrix b(size,size);
	b.setZero();
	b(0,0)=sqrt(2.0);
	b(0,1)=sqrt(2.0)*std::complex<double>(0.2,0.3);
	b(1,0)=sqrt(2.0)*std::complex<double>(0.2,-0.3);
	b(1,1)=-sqrt(2.0);

	cntr::green_from_H(A_matrix,0.0,a,beta,h);
	cntr::green_from_H(B_matrix,0.0,b,beta,h);

	cntr::convolution(C_matrix, A_matrix, A_matrix, B_matrix, B_matrix, CINTEG(MAX_SOLVE_ORDER), beta, h);
	cntr::convolution(Ccc_matrix, B_matrix, B_matrix, A_matrix, A_matrix, CINTEG(MAX_SOLVE_ORDER), beta, h);

	SECTION ("Get for herm_matrix - scalar"){
		std::complex<double> les1,les2,ret1,ret2,gtr1,gtr2;
		std::complex<double> mat1,mat2,tv1,tv2,vt1,vt2;
		double err=0.0;
		for (int m=0; m <= ntau; m++){
			cntr::get_mat(m, mat1, C_scalar, Ccc_scalar);
			C_scalar.get_mat(m, mat2);
			err += std::abs(mat1-mat2);
			cntr::get_mat(m, mat1, C_scalar);
			C_scalar.get_mat(m, mat2);
			err += std::abs(mat1-mat2);
		}	

		for(int i = 0; i <= nt; i++){
			for (int j = 0; j <= nt; j++)
			{
				cntr::get_les(i, j, les1, C_scalar, Ccc_scalar);
				if(i <= j) {
					C_scalar.get_les(i, j, les2);
				}
				else {
					Ccc_scalar.get_les(j, i, les2);
					les2 = -std::conj(les2);
				}
				err += std::abs(les1-les2);

				cntr::get_ret(i, j, ret1, C_scalar, Ccc_scalar);
				if(j <= i) {
					C_scalar.get_ret(i, j, ret2);
				}
				else {
					Ccc_scalar.get_ret(j, i, ret2);
					ret2 = -std::conj(ret2);
				}
				err += std::abs(ret1-ret2);
				cntr::get_gtr(i, j, gtr1, C_scalar, Ccc_scalar);
				if(j <= i) {
					C_scalar.get_gtr(i, j, gtr2);
				}
				else {
					Ccc_scalar.get_gtr(j, i, gtr2);
					gtr2 = -std::conj(gtr2);
				}
				err += std::abs(gtr1-gtr2);
			}
		}
		for(int i = 0; i <= nt; i++){
			for (int m=0; m <= ntau; m++){
				cntr::get_tv(i, m, tv1, C_scalar, Ccc_scalar);
				C_scalar.get_tv(i, m, tv2);
				err += std::abs(tv1-tv2);
				cntr::get_vt(m, i, vt1, C_scalar, Ccc_scalar);
				Ccc_scalar.get_vt(m, i, vt2);
				err += std::abs(vt1-vt2);
			}
		}

		REQUIRE(err < eps);
	}

		SECTION ("Get for herm_matrix - matrix"){
		cdmatrix les1,les2,ret1,ret2,gtr1,gtr2;
		cdmatrix mat1,mat2,tv1,tv2,vt1,vt2;
		double err=0.0;
		for (int m=0; m <= ntau; m++){
			cntr::get_mat(m, mat1, C_matrix, Ccc_matrix);
			C_matrix.get_mat(m, mat2);
			err += (mat1-mat2).norm();
			cntr::get_mat(m, mat1, C_matrix);
			C_matrix.get_mat(m, mat2);
			err += (mat1-mat2).norm();
		}	

		for(int i = 0; i <= nt; i++){
			for (int j = 0; j <= nt; j++)
			{
				cntr::get_les(i, j, les1, C_matrix, Ccc_matrix);
				if(i <= j) {
					C_matrix.get_les(i, j, les2);
				}
				else {
					Ccc_matrix.get_les(j, i, les2);
					les2.adjointInPlace();
					les2 = -les2;
				}
				// std::cout << i << " " << j << " " << les1 << " " << les2 << std::endl;
				err += (les1-les2).norm();

				cntr::get_ret(i, j, ret1, C_matrix, Ccc_matrix);
				if(j <= i) {
					C_matrix.get_ret(i, j, ret2);
				}
				else {
					Ccc_matrix.get_ret(j, i, ret2);
					ret2.adjointInPlace();
					ret2 = -ret2;
				}
				err += (ret1-ret2).norm();
				// std::cout << i << " " << j << " " << (ret1 - ret2).norm() << std::endl;

				cntr::get_gtr(i, j, gtr1, C_matrix, Ccc_matrix);
				gtr2 = ret2 + les2;
				err += (gtr1-gtr2).norm();
			}
		}
		for(int i = 0; i <= nt; i++){
			for (int m=0; m <= ntau; m++){
				cntr::get_tv(i, m, tv1, C_matrix, Ccc_matrix);
				C_matrix.get_tv(i, m, tv2);
				err += (tv1-tv2).norm();
				cntr::get_vt(m, i, vt1, C_matrix, Ccc_matrix);
				Ccc_matrix.get_vt(m, i, vt2);
				err += (vt1-vt2).norm();
			}
		}

		REQUIRE(err < eps);
	}


	SECTION ("Get for herm_matrix_timestep - scalar"){
		std::complex<double> les1,les2,ret1,ret2,gtr1,gtr2;
		std::complex<double> mat1,mat2,tv1,tv2,vt1,vt2;

		GREEN_TSTP C_step(nt,ntau,1,FERMION);
		GREEN_TSTP Ccc_step(nt,ntau,1,FERMION);
		C_scalar.get_timestep(-1, C_step);
		Ccc_scalar.get_timestep(-1, Ccc_step);
		double err=0.0;
		for (int m=0; m <= ntau; m++){
			cntr::get_mat(m, mat1, C_step, Ccc_step);
			C_scalar.get_mat(m, mat2);
			err += std::abs(mat1-mat2);
			cntr::get_mat(m, mat1, C_scalar);
			C_scalar.get_mat(m, mat2);
			err += std::abs(mat1-mat2);
		}	

		for(int i = 0; i <= nt; i++){
			for (int j = 0; j <= nt; j++)
			{
				if(j <= i) {
					C_scalar.get_timestep(i, C_step);
					Ccc_scalar.get_timestep(i, Ccc_step);
				}
				else {
					C_scalar.get_timestep(j, C_step);
					Ccc_scalar.get_timestep(j, Ccc_step);
				}

				cntr::get_les(i, j, les1, C_step, Ccc_step);
				if(i <= j) {					
					C_scalar.get_les(i, j, les2);
				}
				else {
					Ccc_scalar.get_les(j, i, les2);
					les2 = -std::conj(les2);
				}
				err += std::abs(les1-les2);

				cntr::get_ret(i, j, ret1, C_step, Ccc_step);
				if(j <= i) {
					C_scalar.get_ret(i, j, ret2);
				}
				else {
					Ccc_scalar.get_ret(j, i, ret2);
					ret2 = -std::conj(ret2);
				}
				err += std::abs(ret1-ret2);

				cntr::get_gtr(i, j, gtr1, C_step, Ccc_step);
				if(j <= i) {
					C_scalar.get_gtr(i, j, gtr2);
				}
				else {
					Ccc_scalar.get_gtr(j, i, gtr2);
					gtr2 = -std::conj(gtr2);
				}
				err += std::abs(gtr1-gtr2);
			}
		}
		for(int i = 0; i <= nt; i++){
			C_scalar.get_timestep(i, C_step);
			Ccc_scalar.get_timestep(i, Ccc_step);
			for (int m=0; m <= ntau; m++){
				cntr::get_tv(i, m, tv1, C_step, Ccc_step);
				C_scalar.get_tv(i, m, tv2);
				err += std::abs(tv1-tv2);
				cntr::get_vt(m, i, vt1, C_step, Ccc_step);
				Ccc_scalar.get_vt(m, i, vt2);
				err += std::abs(vt1-vt2);
			}
		}

		REQUIRE(err < eps);
	}


	SECTION ("Get for herm_matrix_timestep - matrix"){
		cdmatrix les1,les2,ret1,ret2,gtr1,gtr2;
		cdmatrix mat1,mat2,tv1,tv2,vt1,vt2;

		GREEN_TSTP C_step(nt,ntau,size,FERMION);
		GREEN_TSTP Ccc_step(nt,ntau,size,FERMION);
		C_matrix.get_timestep(-1, C_step);
		Ccc_matrix.get_timestep(-1, Ccc_step);
		double err=0.0;
		for (int m=0; m <= ntau; m++){
			cntr::get_mat(m, mat1, C_step, Ccc_step);
			C_matrix.get_mat(m, mat2);
			err += (mat1-mat2).norm();
			// std::cout << m << " " << mat1 << " " << mat2 << std::endl;
			cntr::get_mat(m, mat1, C_matrix);
			C_matrix.get_mat(m, mat2);
			err += (mat1-mat2).norm();
		}	

		for(int i = 0; i <= nt; i++){
			for (int j = 0; j <= nt; j++)
			{
				if(j <= i) {
					C_matrix.get_timestep(i, C_step);
					Ccc_matrix.get_timestep(i, Ccc_step);
				}
				else {
					C_matrix.get_timestep(j, C_step);
					Ccc_matrix.get_timestep(j, Ccc_step);
				}

				cntr::get_les(i, j, les1, C_step, Ccc_step);
				if(i <= j) {					
					C_matrix.get_les(i, j, les2);
				}
				else {
					Ccc_matrix.get_les(j, i, les2);
					les2.adjointInPlace();
					les2 = -les2;
				}
				err += (les1-les2).norm();

				cntr::get_ret(i, j, ret1, C_step, Ccc_step);
				if(j <= i) {
					C_matrix.get_ret(i, j, ret2);
				}
				else {
					Ccc_matrix.get_ret(j, i, ret2);
					ret2.adjointInPlace();
					ret2 = -ret2;
				}
				err += (ret1-ret2).norm();

				cntr::get_gtr(i, j, gtr1, C_step, Ccc_step);
				gtr2 = ret2 + les2;
				err += (gtr1-gtr2).norm();
			}
		}
		for(int i = 0; i <= nt; i++){
			C_matrix.get_timestep(i, C_step);
			Ccc_matrix.get_timestep(i, Ccc_step);
			for (int m=0; m <= ntau; m++){
				cntr::get_tv(i, m, tv1, C_step, Ccc_step);
				C_matrix.get_tv(i, m, tv2);
				err += (tv1-tv2).norm();
				cntr::get_vt(m, i, vt1, C_step, Ccc_step);
				Ccc_matrix.get_vt(m, i, vt2);
				err += (vt1-vt2).norm();
			}
		}

		REQUIRE(err < eps);
	}




	SECTION ("Get for herm_matrix_timestep_view - scalar"){
		std::complex<double> les1,les2,ret1,ret2,gtr1,gtr2;
		std::complex<double> mat1,mat2,tv1,tv2,vt1,vt2;

		GREEN_TSTP_VIEW C_step = GREEN_TSTP_VIEW(-1, C_scalar);
		GREEN_TSTP_VIEW Ccc_step = GREEN_TSTP_VIEW(-1, Ccc_scalar);
		double err=0.0;
		for (int m=0; m <= ntau; m++){
			cntr::get_mat(m, mat1, C_step, Ccc_step);
			C_scalar.get_mat(m, mat2);
			err += std::abs(mat1-mat2);
			cntr::get_mat(m, mat1, C_scalar);
			C_scalar.get_mat(m, mat2);
			err += std::abs(mat1-mat2);
		}	

		for(int i = 0; i <= nt; i++){
			GREEN_TSTP_VIEW C_step = GREEN_TSTP_VIEW(i, C_scalar);
			GREEN_TSTP_VIEW Ccc_step = GREEN_TSTP_VIEW(i, Ccc_scalar);
			for (int j = 0; j <= i; j++)
			{

				cntr::get_les(j, i, les1, C_step, Ccc_step);
				C_scalar.get_les(j, i, les2);
				err += std::abs(les1-les2);

				cntr::get_ret(i, j, ret1, C_step, Ccc_step);
				C_scalar.get_ret(i, j, ret2);
				err += std::abs(ret1-ret2);
			}
		}
		for(int i = 0; i <= nt; i++){
			for (int j = i+1; j <= nt; j++)
			{
				GREEN_TSTP_VIEW C_step = GREEN_TSTP_VIEW(j, C_scalar);
				GREEN_TSTP_VIEW Ccc_step = GREEN_TSTP_VIEW(j, Ccc_scalar);
				cntr::get_les(i, j, les1, C_step, Ccc_step);
				C_scalar.get_les(i, j, les2);
				err += std::abs(les1-les2);

				cntr::get_ret(j, i, ret1, C_step, Ccc_step);
				C_scalar.get_ret(j, i, ret2);
				err += std::abs(ret1-ret2);
			}
		}
		for(int i = 0; i <= nt; i++){
			GREEN_TSTP_VIEW C_step = GREEN_TSTP_VIEW(i, C_scalar);
			GREEN_TSTP_VIEW Ccc_step = GREEN_TSTP_VIEW(i, Ccc_scalar);
			for (int m=0; m <= ntau; m++){
				cntr::get_tv(i, m, tv1, C_step, Ccc_step);
				C_scalar.get_tv(i, m, tv2);
				err += std::abs(tv1-tv2);
				cntr::get_vt(m, i, vt1, C_step, Ccc_step);
				Ccc_scalar.get_vt(m, i, vt2);
				err += std::abs(vt1-vt2);
			}
		}

		REQUIRE(err < eps);
	}


	SECTION ("Get for herm_matrix_timestep_view - matrix"){
		cdmatrix les1,les2,ret1,ret2,gtr1,gtr2;
		cdmatrix mat1,mat2,tv1,tv2,vt1,vt2;

		GREEN_TSTP_VIEW C_step = GREEN_TSTP_VIEW(-1, C_matrix);
		GREEN_TSTP_VIEW Ccc_step = GREEN_TSTP_VIEW(-1, Ccc_matrix);
		double err=0.0;
		for (int m=0; m <= ntau; m++){
			cntr::get_mat(m, mat1, C_step, Ccc_step);
			C_matrix.get_mat(m, mat2);
			err += (mat1-mat2).norm();
			cntr::get_mat(m, mat1, C_matrix);
			C_matrix.get_mat(m, mat2);
			err += (mat1-mat2).norm();
		}	

		for(int i = 0; i <= nt; i++){
			GREEN_TSTP_VIEW C_step = GREEN_TSTP_VIEW(i, C_matrix);
			GREEN_TSTP_VIEW Ccc_step = GREEN_TSTP_VIEW(i, Ccc_matrix);
			for (int j = 0; j <= i; j++)
			{

				cntr::get_les(j, i, les1, C_step, Ccc_step);
				C_matrix.get_les(j, i, les2);
				err += (les1-les2).norm();

				cntr::get_ret(i, j, ret1, C_step, Ccc_step);
				C_matrix.get_ret(i, j, ret2);
				err += (ret1-ret2).norm();
			}
		}
		for(int i = 0; i <= nt; i++){
			for (int j = i+1; j <= nt; j++)
			{
				GREEN_TSTP_VIEW C_step = GREEN_TSTP_VIEW(j, C_matrix);
				GREEN_TSTP_VIEW Ccc_step = GREEN_TSTP_VIEW(j, Ccc_matrix);
				cntr::get_les(i, j, les1, C_step, Ccc_step);
				C_matrix.get_les(i, j, les2);
				err += (les1-les2).norm();

				cntr::get_ret(j, i, ret1, C_step, Ccc_step);
				C_matrix.get_ret(j, i, ret2);
				err += (ret1-ret2).norm();
			}
		}
		for(int i = 0; i <= nt; i++){
			GREEN_TSTP_VIEW C_step = GREEN_TSTP_VIEW(i, C_matrix);
			GREEN_TSTP_VIEW Ccc_step = GREEN_TSTP_VIEW(i, Ccc_matrix);
			for (int m=0; m <= ntau; m++){
				cntr::get_tv(i, m, tv1, C_step, Ccc_step);
				C_matrix.get_tv(i, m, tv2);
				err += (tv1-tv2).norm();
				cntr::get_vt(m, i, vt1, C_step, Ccc_step);
				Ccc_matrix.get_vt(m, i, vt2);
				err += (vt1-vt2).norm();
			}
		}

		REQUIRE(err < eps);
	}

}
