#include "catch.hpp"
#include "cntr.hpp"

#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define GREEN_TSTP_VIEW cntr::herm_matrix_timestep_view<double>

TEST_CASE("Global get","[Global get]"){
	int nt=50;
	int ntau=500;
	int size=2;
	double eps=1e-7;
	double beta=5.0;
	double h=0.01;

	// scalar Green's function
	GREEN A_scalar=GREEN(nt,ntau,1,BOSON);
	cdmatrix a1x1(1,1);
	a1x1(0,0) = 1.0;
	cntr::green_from_H(A_scalar,0.0,a1x1,beta,h);

	// Square Green's function
	GREEN A_matrix=GREEN(nt,ntau,size,BOSON);

	cdmatrix a(size,size);
	a.setZero();
	a(0,0)=sqrt(2.0);
	a(0,1)=sqrt(2.0)*std::complex<double>(0.0,1.0);
	a(1,0)=sqrt(2.0)*std::complex<double>(0.0,-1.0);
	a(1,1)=-sqrt(2.0);
	// for(int i=0; i<size; i++) a(i,i) = 0.1*(i-4);

	cntr::green_from_H(A_matrix,0.0,a,beta,h);

	SECTION ("Get for herm_matrix - scalar"){
		std::complex<double> les1,les2,ret1,ret2,gtr1,gtr2;
		std::complex<double> mat1,mat2,tv1,tv2,vt1,vt2;
		double err=0.0;
		for (int m=0; m <= ntau; m++){
			cntr::get_mat(m, mat1, A_scalar, A_scalar);
			A_scalar.get_mat(m, mat2);
			err += std::abs(mat1-mat2);
			cntr::get_mat(m, mat1, A_scalar);
			A_scalar.get_mat(m, mat2);
			err += std::abs(mat1-mat2);
		}	

		for(int i = 0; i <= nt; i++){
			for (int j = 0; j <= nt; j++)
			{
				cntr::get_les(i, j, les1, A_scalar, A_scalar);
				A_scalar.get_les(i, j, les2);
				err += std::abs(les1-les2);
				cntr::get_ret(i, j, ret1, A_scalar, A_scalar);
				A_scalar.get_ret(i, j, ret2);
				err += std::abs(ret1-ret2);
				cntr::get_gtr(i, j, gtr1, A_scalar, A_scalar);
				A_scalar.get_gtr(i, j, gtr2);
				err += std::abs(gtr1-gtr2);
			}
		}
		for(int i = 0; i <= nt; i++){
			for (int m=0; m <= ntau; m++){
				cntr::get_tv(i, m, tv1, A_scalar, A_scalar);
				A_scalar.get_tv(i, m, tv2);
				err += std::abs(tv1-tv2);
				cntr::get_vt(m, i, vt1, A_scalar, A_scalar);
				A_scalar.get_vt(m, i, vt2);
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
			cntr::get_mat(m, mat1, A_matrix, A_matrix);
			A_matrix.get_mat(m, mat2);
			err += (mat1-mat2).norm();
			cntr::get_mat(m, mat1, A_matrix);
			A_matrix.get_mat(m, mat2);
			err += (mat1-mat2).norm();
		}	

		for(int i = 0; i <= nt; i++){
			for (int j = 0; j <= nt; j++)
			{
				cntr::get_les(i, j, les1, A_matrix, A_matrix);
				A_matrix.get_les(i, j, les2);
				err += (les1-les2).norm();
				cntr::get_ret(i, j, ret1, A_matrix, A_matrix);
				A_matrix.get_ret(i, j, ret2);
				err += (ret1-ret2).norm();
				cntr::get_gtr(i, j, gtr1, A_matrix, A_matrix);
				A_matrix.get_gtr(i, j, gtr2);
				err += (gtr1-gtr2).norm();
			}
		}
		for(int i = 0; i <= nt; i++){
			for (int m=0; m <= ntau; m++){
				cntr::get_tv(i, m, tv1, A_matrix, A_matrix);
				A_matrix.get_tv(i, m, tv2);
				err += (tv1-tv2).norm();
				cntr::get_vt(m, i, vt1, A_matrix, A_matrix);
				A_matrix.get_vt(m, i, vt2);
				err += (vt1-vt2).norm();
			}
		}

		REQUIRE(err < eps);
	}

	SECTION ("herm_matrix_timestep.get_XXX - scalar"){
		std::complex<double> les1,les2,ret1,ret2;
		std::complex<double> mat1,mat2,tv1,tv2,vt1,vt2;
		GREEN_TSTP A_step(-1,ntau,1,BOSON);
		A_scalar.get_timestep(-1, A_step);
		double err=0.0;
		for (int m=0; m <= ntau; m++){
			A_scalar.get_mat(m, mat1);
			A_step.get_mat(m, mat2);
			err += std::abs(mat1-mat2);
		}	

		for(int i = 0; i <= nt; i++){
			GREEN_TSTP A_step(i,ntau,1,BOSON);
			A_scalar.get_timestep(i, A_step);
			for (int j = 0; j <= i; j++)
			{
				A_scalar.get_les(i, j, les1);
				A_step.get_les(i, j, les2);
				err += std::abs(les1-les2);
				A_scalar.get_ret(i, j, ret1);
				A_step.get_ret(i, j, ret2);
				err += std::abs(ret1-ret2);
			}
		}
		for(int i = 0; i <= nt; i++){
			GREEN_TSTP A_step(i,ntau,1,BOSON);
			A_scalar.get_timestep(i, A_step);
			for (int m=0; m <= ntau; m++){
				A_scalar.get_tv(i, m, tv1);
				A_step.get_tv(i, m, tv2);
				err += std::abs(tv1-tv2);
				A_scalar.get_vt(m, i, vt1);
				A_step.get_vt(m, i, vt2);
				err += std::abs(vt1-vt2);
			}
		}

		REQUIRE(err < eps);
	}

	SECTION ("herm_matrix_timestep.get_XXX - matrix"){
		cdmatrix les1,les2,ret1,ret2;
		cdmatrix mat1,mat2,tv1,tv2;
		GREEN_TSTP A_step(-1,ntau,size);
		A_matrix.get_timestep(-1, A_step);
		double err=0.0;
		for (int m=0; m <= ntau; m++){
			A_matrix.get_mat(m, mat1);
			A_step.get_mat(m, mat2);
			err += (mat1-mat2).norm();
		}	

		for(int i = 0; i <= nt; i++){
			GREEN_TSTP A_step(i,ntau,size);
			A_matrix.get_timestep(i, A_step);
			for (int j = 0; j <= i; j++)
			{
				A_matrix.get_les(j, i, les1);
				A_step.get_les(j, i, les2);
				err += (les1-les2).norm();
				A_matrix.get_ret(i, j, ret1);
				A_step.get_ret(i, j, ret2);
				err += (ret1-ret2).norm();
			}
		}
		for(int i = 0; i <= nt; i++){
			GREEN_TSTP A_step(i,ntau,size);
			A_matrix.get_timestep(i, A_step);
			for (int m=0; m <= ntau; m++){
				A_matrix.get_tv(i, m, tv1);
				A_step.get_tv(i, m, tv2);
				err += (tv1-tv2).norm();
			}
		}

		REQUIRE(err < eps);
	}

	SECTION ("Get herm_matrix_timestep_view - scalar"){
		std::complex<double> les1,les2,ret1,ret2,gtr1,gtr2;
		std::complex<double> mat1,mat2,tv1,tv2,vt1,vt2;
		GREEN_TSTP_VIEW A_step = GREEN_TSTP_VIEW(-1, A_scalar);
		double err=0.0;
		for (int m=0; m <= ntau; m++){
			A_scalar.get_mat(m, mat1);
			cntr::get_mat(m, mat2, A_step);
			err += std::abs(mat1-mat2);
			cntr::get_mat(m, mat2, A_step, A_step);
			err += std::abs(mat1-mat2);
		}	

		for(int i = 0; i <= nt; i++){
			GREEN_TSTP_VIEW A_step = GREEN_TSTP_VIEW(i, A_scalar);
			for (int j = 0; j <= i; j++)
			{
				A_scalar.get_les(i, j, les1);
				cntr::get_les(i, j, les2, A_step, A_step);
				err += std::abs(les1-les2);
				A_scalar.get_ret(i, j, ret1);
				cntr::get_ret(i, j, ret2, A_step, A_step);
				err += std::abs(ret1-ret2);
				A_scalar.get_gtr(i, j, gtr1);
				cntr::get_gtr(i, j, gtr2, A_step, A_step);
				err += std::abs(gtr1-gtr2);
			}
		}

		for(int i = 0; i <= nt; i++){
			for (int j = i+1; j <= nt; j++)
			{
				GREEN_TSTP_VIEW A_step = GREEN_TSTP_VIEW(j, A_scalar);
				A_scalar.get_les(i, j, les1);
				cntr::get_les(i, j, les2, A_step, A_step);
				err += std::abs(les1-les2);
				A_scalar.get_ret(i, j, ret1);
				cntr::get_ret(i, j, ret2, A_step, A_step);
				err += std::abs(ret1-ret2);
				A_scalar.get_gtr(i, j, gtr1);
				cntr::get_gtr(i, j, gtr2, A_step, A_step);
				err += std::abs(gtr1-gtr2);
			}
		}

		for(int i = 0; i <= nt; i++){
			GREEN_TSTP_VIEW A_step = GREEN_TSTP_VIEW(i, A_scalar);
			for (int m=0; m <= ntau; m++){
				A_scalar.get_tv(i, m, tv1);
				cntr::get_tv(i, m, tv2, A_step, A_step);
				err += std::abs(tv1-tv2);
				A_scalar.get_vt(m, i, vt1);
				cntr::get_vt(m, i, vt2, A_step, A_step);
				err += std::abs(vt1-vt2);
			}
		}

		REQUIRE(err < eps);
	}


	SECTION ("Get herm_matrix_timestep_view - matrix"){
		cdmatrix les1,les2,ret1,ret2,gtr1,gtr2;
		cdmatrix mat1,mat2,tv1,tv2,vt1,vt2;
		GREEN_TSTP_VIEW A_step = GREEN_TSTP_VIEW(-1, A_matrix);
		double err=0.0;
		for (int m=0; m <= ntau; m++){
			A_matrix.get_mat(m, mat1);
			cntr::get_mat(m, mat2, A_step);
			err += (mat1-mat2).norm();
			cntr::get_mat(m, mat2, A_step, A_step);
			err += (mat1-mat2).norm();
		}	

		for(int i = 0; i <= nt; i++){
			GREEN_TSTP_VIEW A_step = GREEN_TSTP_VIEW(i, A_matrix);
			for (int j = 0; j <= i; j++)
			{
				A_matrix.get_les(i, j, les1);
				cntr::get_les(i, j, les2, A_step, A_step);
				err += (les1-les2).norm();
				A_matrix.get_ret(i, j, ret1);
				cntr::get_ret(i, j, ret2, A_step, A_step);
				err += (ret1-ret2).norm();
				A_matrix.get_gtr(i, j, gtr1);
				cntr::get_gtr(i, j, gtr2, A_step, A_step);
				err += (gtr1-gtr2).norm();
			}
		}

		for(int i = 0; i <= nt; i++){
			for (int j = i+1; j <= nt; j++)
			{
				GREEN_TSTP_VIEW A_step = GREEN_TSTP_VIEW(j, A_matrix);
				A_matrix.get_les(i, j, les1);
				cntr::get_les(i, j, les2, A_step, A_step);
				err += (les1-les2).norm();
				A_matrix.get_ret(i, j, ret1);
				cntr::get_ret(i, j, ret2, A_step, A_step);
				err += (ret1-ret2).norm();
				A_matrix.get_gtr(i, j, gtr1);
				cntr::get_gtr(i, j, gtr2, A_step, A_step);
				err += (gtr1-gtr2).norm();
			}
		}

		for(int i = 0; i <= nt; i++){
			GREEN_TSTP_VIEW A_step = GREEN_TSTP_VIEW(i, A_matrix);
			for (int m=0; m <= ntau; m++){
				A_matrix.get_tv(i, m, tv1);
				cntr::get_tv(i, m, tv2, A_step, A_step);
				err += (tv1-tv2).norm();
				A_matrix.get_vt(m, i, vt1);
				cntr::get_vt(m, i, vt2, A_step, A_step);
				err += (vt1-vt2).norm();
			}
		}

		REQUIRE(err < eps);
	}


}
