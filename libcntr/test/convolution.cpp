#include "catch.hpp"
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

#include "cntr.hpp"

#define GREEN cntr::herm_matrix<double>
#define CPLX std::complex<double>
#define cfunction cntr::function<double>
using namespace std;

/*///////////////////////////////////////////////////////////////////////////////////////


  This is a minimal test program to solve the convolution

  A*B=C

  for bosonic Greens functions:

  A= 1/(iwn-wa)
  B= 1/(iwn-wb)
  C= (A-B)/(wa-wb)

///////////////////////////////////////////////////////////////////////////////////////*/

TEST_CASE("convolution: size1x1","[convolution: size1x1]"){

  int nt,ntau,kt,tstp;
  double wa,wb,h,beta,fac,mu;
  GREEN D;
  GREEN A_fer,B_fer,C_fer,AB_fer;
  GREEN A_bos,B_bos,C_bos,AB_bos,ABtemp_bos;

  nt=7;
  ntau=500;
  beta=05.0;
  h=0.01;
  wa=1.123;
  wb=0.345;
  kt=5;
	mu=0.0;
  double eps=1e-6;

  SECTION("convolution (bose)"){

    A_bos=GREEN(nt,ntau,1,1);
    B_bos=GREEN(nt,ntau,1,1);
    C_bos=GREEN(nt,ntau,1,1);
    AB_bos=GREEN(nt,ntau,1,1);

		//ABtemp_bos=GREEN(nt,ntau,1,1);

		cdmatrix eps_a(1,1);
		cdmatrix eps_b(1,1);

		eps_a.setZero();
		eps_b.setZero();

		eps_a(0,0)=wa;
		eps_b(0,0)=wb;

		cntr::green_from_H(A_bos,mu,eps_a,beta,h);
		cntr::green_from_H(B_bos,mu,eps_b,beta,h);

    //cntr::green_single_pole_bose(ABtemp_bos,&wa,beta,h);

		// C=(A-B)/(wa-wb)
    fac=1.0/(wa-wb);
    C_bos=A_bos;
    for(tstp=-1;tstp<=nt;tstp++){
      C_bos.incr_timestep(tstp,B_bos,CPLX(-1.0,0.0));
      C_bos.smul(tstp,fac);
    }
    //cout << "SIMPLE TEST WITH GIVEN h=0.02, dtau=0.005" << endl;

		cntr::convolution(AB_bos,A_bos,A_bos,B_bos,B_bos,integration::I<double>(kt),beta,h);
    //cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);

    double err=0.0;
		for(tstp=-1;tstp<=nt;tstp++){
      err+=cntr::distance_norm2(tstp,C_bos,AB_bos);
      err+=cntr::distance_norm2_ret(tstp,C_bos,AB_bos);
      err+=cntr::distance_norm2_tv(tstp,C_bos,AB_bos);
      err+=cntr::distance_norm2_les(tstp,C_bos,AB_bos);

			//err+=cntr::distance_norm2(tstp,A_bos,ABtemp_bos);
    }
    //err+=-A_bos.sig();

    REQUIRE(err<eps);
  }

  // TEST REAL-TIME ... SET MATSUBARA TO EXACT
  SECTION("convolution (bosonic signs in boundary terms)"){
    double tmax=2.0;
    int nt0=10,fac2;
    beta=0.1;
    double err=0.0;

    for(fac2=1;fac2<=64;fac2*=2){
      nt=nt0*fac2;

      A_bos=GREEN(nt,ntau,1,1);
      B_bos=GREEN(nt,ntau,1,1);
      C_bos=GREEN(nt,ntau,1,1);
      AB_bos=GREEN(nt,ntau,1,1);

      h=tmax/nt;

			cdmatrix eps_a(1,1);
			cdmatrix eps_b(1,1);

			eps_a.setZero();
			eps_b.setZero();

			eps_a(0,0)=wa;
			eps_b(0,0)=wb;

			cntr::green_from_H(A_bos,0.0,eps_a,beta,h);
			cntr::green_from_H(B_bos,0.0,eps_b,beta,h);

      // C=(A-B)/(wa-wb)
      fac=1.0/(wa-wb);
      C_bos=A_bos;
      for(tstp=-1;tstp<=nt;tstp++){
					C_bos.incr_timestep(tstp,B_bos,CPLX(-1.0,0.0));
					C_bos.smul(tstp,fac);
      }
      cntr::convolution(AB_bos,A_bos,A_bos,B_bos,B_bos,integration::I<double>(kt),beta,h);
      //cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
      for(tstp=nt;tstp<=nt;tstp++){
					err=cntr::distance_norm2(tstp,C_bos,AB_bos);
					err=cntr::distance_norm2_ret(tstp,C_bos,AB_bos);
					err=cntr::distance_norm2_tv(tstp,C_bos,AB_bos);
					err=cntr::distance_norm2_les(tstp,C_bos,AB_bos);
      }
    }
    REQUIRE(err<eps);
  }

  SECTION("convolution (fermi)"){

    A_fer=GREEN(nt,ntau,1,-1);
    B_fer=GREEN(nt,ntau,1,-1);
    C_fer=GREEN(nt,ntau,1,-1);
    AB_fer=GREEN(nt,ntau,1,-1);

    cdmatrix eps_a(1,1);
    cdmatrix eps_b(1,1);

    eps_a.setZero();
    eps_b.setZero();

    eps_a(0,0)=wa;
    eps_b(0,0)=wb;

    cntr::green_from_H(A_fer,mu,eps_a,beta,h);
    cntr::green_from_H(B_fer,mu,eps_b,beta,h);

    // C=(A-B)/(wa-wb)
    fac=1.0/(wa-wb);
    C_fer=A_fer;
    for(tstp=-1;tstp<=nt;tstp++){
      C_fer.incr_timestep(tstp,B_fer,CPLX(-1.0,0.0));
      C_fer.smul(tstp,fac);
    }
    //cout << "SIMPLE TEST WITH GIVEN h=0.02, dtau=0.005" << endl;
    cntr::convolution(AB_fer,A_fer,A_fer,B_fer,B_fer,integration::I<double>(kt),beta,h);
    //cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
    double err=0.0;
    for(tstp=-1;tstp<=nt;tstp++){
      err+=cntr::distance_norm2(tstp,C_fer,AB_fer);
      err+=cntr::distance_norm2_ret(tstp,C_fer,AB_fer);
      err+=cntr::distance_norm2_tv(tstp,C_fer,AB_fer);
      err+=cntr::distance_norm2_les(tstp,C_fer,AB_fer);
    }

    REQUIRE(err<eps);
  }

	SECTION("convolution_new (bose)"){

    A_bos=GREEN(nt,ntau,1,1);
    B_bos=GREEN(nt,ntau,1,1);
    C_bos=GREEN(nt,ntau,1,1);
    AB_bos=GREEN(nt,ntau,1,1);

		cdmatrix eps_a(1,1);
		cdmatrix eps_b(1,1);

		eps_a.setZero();
		eps_b.setZero();

		eps_a(0,0)=wa;
		eps_b(0,0)=wb;

		cntr::green_from_H(A_bos,0.0,eps_a,beta,h);
		cntr::green_from_H(B_bos,0.0,eps_b,beta,h);

    // C=(A-B)/(wa-wb)
    fac=1.0/(wa-wb);
    C_bos=A_bos;
    for(tstp=-1;tstp<=nt;tstp++){
      C_bos.incr_timestep(tstp,B_bos,CPLX(-1.0,0.0));
      C_bos.smul(tstp,fac);
    }
    //cout << "SIMPLE TEST WITH GIVEN h=0.02, dtau=0.005" << endl;
    cntr::convolution_new(AB_bos,A_bos,A_bos,B_bos,B_bos,integration::I<double>(kt),beta,h);
    //cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
    double err=0.0;
    for(tstp=-1;tstp<=nt;tstp++){
      err+=cntr::distance_norm2(tstp,C_bos,AB_bos);
      err+=cntr::distance_norm2_ret(tstp,C_bos,AB_bos);
      err+=cntr::distance_norm2_tv(tstp,C_bos,AB_bos);
      err+=cntr::distance_norm2_les(tstp,C_bos,AB_bos);
    }

    REQUIRE(err<eps);
  }

  // TEST REAL-TIME ... SET MATSUBARA TO EXACT
  SECTION("convolution_new (bosonic signs in boundary terms)"){
    double tmax=2.0;
    int nt0=10,fac2;
    beta=0.1;
    double err=0.0;

    for(fac2=1;fac2<=64;fac2*=2){
      nt=nt0*fac2;

      A_bos=GREEN(nt,ntau,1,1);
      B_bos=GREEN(nt,ntau,1,1);
      C_bos=GREEN(nt,ntau,1,1);
      AB_bos=GREEN(nt,ntau,1,1);

      h=tmax/nt;

			cdmatrix eps_a(1,1);
			cdmatrix eps_b(1,1);

			eps_a.setZero();
			eps_b.setZero();

			eps_a(0,0)=wa;
			eps_b(0,0)=wb;

			cntr::green_from_H(A_bos,0.0,eps_a,beta,h);
			cntr::green_from_H(B_bos,0.0,eps_b,beta,h);

      // C=(A-B)/(wa-wb)
      fac=1.0/(wa-wb);
      C_bos=A_bos;
      for(tstp=-1;tstp<=nt;tstp++){
				C_bos.incr_timestep(tstp,B_bos,CPLX(-1.0,0.0));
				C_bos.smul(tstp,fac);
      }
      cntr::convolution_new(AB_bos,A_bos,A_bos,B_bos,B_bos,integration::I<double>(kt),beta,h);
      //cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
      for(tstp=nt;tstp<=nt;tstp++){
				err=cntr::distance_norm2(tstp,C_bos,AB_bos);
				err=cntr::distance_norm2_ret(tstp,C_bos,AB_bos);
				err=cntr::distance_norm2_tv(tstp,C_bos,AB_bos);
				err=cntr::distance_norm2_les(tstp,C_bos,AB_bos);
      }
    }
    REQUIRE(err<eps);
  }

	SECTION("convolution_new (fermi)"){

    A_fer=GREEN(nt,ntau,1,-1);
    B_fer=GREEN(nt,ntau,1,-1);
    C_fer=GREEN(nt,ntau,1,-1);
    AB_fer=GREEN(nt,ntau,1,-1);

    cdmatrix eps_a(1,1);
    cdmatrix eps_b(1,1);

    eps_a.setZero();
    eps_b.setZero();

    eps_a(0,0)=wa;
    eps_b(0,0)=wb;

    cntr::green_from_H(A_fer,mu,eps_a,beta,h);
    cntr::green_from_H(B_fer,mu,eps_b,beta,h);

		// C=(A-B)/(wa-wb)
    fac=1.0/(wa-wb);
    C_fer=A_fer;
    for(tstp=-1;tstp<=nt;tstp++){
      C_fer.incr_timestep(tstp,B_fer,CPLX(-1.0,0.0));
      C_fer.smul(tstp,fac);
    }

    //cout << "SIMPLE TEST WITH GIVEN h=0.02, dtau=0.005" << endl;
    cntr::convolution_new(AB_fer,A_fer,A_fer,B_fer,B_fer,integration::I<double>(kt),beta,h);
		//cntr::convolution(AB_fer,A_fer,A_fer,B_fer,B_fer,integration::I<double>(kt),beta,h);
    //cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
    double err=0.0;
    for(tstp=-1;tstp<=nt;tstp++){
      err+=cntr::distance_norm2(tstp,C_fer,AB_fer);
      err+=cntr::distance_norm2_ret(tstp,C_fer,AB_fer);
      err+=cntr::distance_norm2_tv(tstp,C_fer,AB_fer);
      err+=cntr::distance_norm2_les(tstp,C_fer,AB_fer);
    }
    //err+=A_fer.sig();

    REQUIRE(err<eps);
  }

#if CNTR_USE_OMP==1
	SECTION("convolution_omp (bose)"){

		A_bos=GREEN(nt,ntau,1,1);
		B_bos=GREEN(nt,ntau,1,1);
		C_bos=GREEN(nt,ntau,1,1);
		AB_bos=GREEN(nt,ntau,1,1);

		int nthreads =4;

		cdmatrix eps_a(1,1);
		cdmatrix eps_b(1,1);

		eps_a.setZero();
		eps_b.setZero();

		eps_a(0,0)=wa;
		eps_b(0,0)=wb;

		cntr::green_from_H(A_bos,0.0,eps_a,beta,h);
		cntr::green_from_H(B_bos,0.0,eps_b,beta,h);

		// C=(A-B)/(wa-wb)
		fac=1.0/(wa-wb);
		C_bos=A_bos;
		for(tstp=-1;tstp<=nt;tstp++){
			C_bos.incr_timestep(tstp,B_bos,CPLX(-1.0,0.0));
			C_bos.smul(tstp,fac);
		}
		//cout << "SIMPLE TEST WITH GIVEN h=0.02, dtau=0.005" << endl;
		cntr::convolution_omp(nthreads,AB_bos,A_bos,A_bos,B_bos,B_bos,integration::I<double>(kt),beta,h);
		//cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
		double err=0.0;
		for(tstp=-1;tstp<=nt;tstp++){
			err+=cntr::distance_norm2(tstp,C_bos,AB_bos);
			err+=cntr::distance_norm2_ret(tstp,C_bos,AB_bos);
			err+=cntr::distance_norm2_tv(tstp,C_bos,AB_bos);
			err+=cntr::distance_norm2_les(tstp,C_bos,AB_bos);
		}

		REQUIRE(err<eps);
	}

	// TEST REAL-TIME ... SET MATSUBARA TO EXACT
	SECTION("convolution_omp (bosonic signs in boundary terms)"){
		double tmax=2.0;
		int nt0=10,fac2;
		beta=0.1;
		double err=0.0;

		int nthreads =4;

		for(fac2=1;fac2<=64;fac2*=2){
			nt=nt0*fac2;

			A_bos=GREEN(nt,ntau,1,1);
			B_bos=GREEN(nt,ntau,1,1);
			C_bos=GREEN(nt,ntau,1,1);
			AB_bos=GREEN(nt,ntau,1,1);

			h=tmax/nt;

			cdmatrix eps_a(1,1);
			cdmatrix eps_b(1,1);

			eps_a.setZero();
			eps_b.setZero();

			eps_a(0,0)=wa;
			eps_b(0,0)=wb;

			cntr::green_from_H(A_bos,0.0,eps_a,beta,h);
			cntr::green_from_H(B_bos,0.0,eps_b,beta,h);

			// C=(A-B)/(wa-wb)
			fac=1.0/(wa-wb);
			C_bos=A_bos;
			for(tstp=-1;tstp<=nt;tstp++){
				C_bos.incr_timestep(tstp,B_bos,CPLX(-1.0,0.0));
				C_bos.smul(tstp,fac);
			}
			cntr::convolution_omp(nthreads,AB_bos,A_bos,A_bos,B_bos,B_bos,integration::I<double>(kt),beta,h);
			//cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
			for(tstp=nt;tstp<=nt;tstp++){
				err=cntr::distance_norm2(tstp,C_bos,AB_bos);
				err=cntr::distance_norm2_ret(tstp,C_bos,AB_bos);
				err=cntr::distance_norm2_tv(tstp,C_bos,AB_bos);
				err=cntr::distance_norm2_les(tstp,C_bos,AB_bos);
			}
		}
		REQUIRE(err<eps);
	}

	SECTION("convolution_omp (fermi)"){

		A_fer=GREEN(nt,ntau,1,-1);
		B_fer=GREEN(nt,ntau,1,-1);
		C_fer=GREEN(nt,ntau,1,-1);
		AB_fer=GREEN(nt,ntau,1,-1);

		cdmatrix eps_a(1,1);
		cdmatrix eps_b(1,1);

		int nthreads =4;

/*#pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }*/

		eps_a.setZero();
		eps_b.setZero();

		eps_a(0,0)=wa;
		eps_b(0,0)=wb;

		cntr::green_from_H(A_fer,mu,eps_a,beta,h);
		cntr::green_from_H(B_fer,mu,eps_b,beta,h);

		// C=(A-B)/(wa-wb)
		fac=1.0/(wa-wb);
		C_fer=A_fer;
		for(tstp=-1;tstp<=nt;tstp++){
			C_fer.incr_timestep(tstp,B_fer,CPLX(-1.0,0.0));
			C_fer.smul(tstp,fac);
		}

		//cout << "SIMPLE TEST WITH GIVEN h=0.02, dtau=0.005" << endl;
		cntr::convolution_omp(nthreads, AB_fer,A_fer,A_fer,B_fer,B_fer,integration::I<double>(kt),beta,h);
		//cntr::convolution(AB_fer,A_fer,A_fer,B_fer,B_fer,integration::I<double>(kt),beta,h);
		//cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
		double err=0.0;
		for(tstp=-1;tstp<=nt;tstp++){
			err+=cntr::distance_norm2(tstp,C_fer,AB_fer);
			err+=cntr::distance_norm2_ret(tstp,C_fer,AB_fer);
			err+=cntr::distance_norm2_tv(tstp,C_fer,AB_fer);
			err+=cntr::distance_norm2_les(tstp,C_fer,AB_fer);
		}
		//err+=A_fer.sig();

		REQUIRE(err<eps);
	}
#endif

}


TEST_CASE("convolution: size2x2","[convolution: size2x2]"){

  int nt,ntau,kt,tstp;
  double wa,wb,h,beta,fac,mu;
  GREEN D;
  GREEN A_fer,B_fer,C_fer,AB_fer;
	GREEN A00_fer,A01_fer,A10_fer,A11_fer;
	GREEN B00_fer,B01_fer,B10_fer,B11_fer;
	GREEN AB00_fer,AB01_fer,AB10_fer,AB11_fer;
	GREEN C00_fer,C01_fer,C10_fer,C11_fer,Ctemp_fer;

  GREEN A_bos,B_bos,C_bos,AB_bos;
	GREEN A00_bos,A01_bos,A10_bos,A11_bos;
	GREEN B00_bos,B01_bos,B10_bos,B11_bos;
	GREEN AB00_bos,AB01_bos,AB10_bos,AB11_bos;
	GREEN C00_bos,C01_bos,C10_bos,C11_bos,Ctemp_bos;

  int size_=2;

  nt=10;
  ntau=400;
  beta=0.1;
  h=0.02;
  wa=1.123;
  wb=0.345;
  kt=5;
	mu=0.0;
  double eps=1e-6;

  SECTION("convolution (bose)"){

    A_bos=GREEN(nt,ntau,size_,1);
    B_bos=GREEN(nt,ntau,size_,1);
    AB_bos=GREEN(nt,ntau,size_,1);

		cdmatrix eps_a(size_,size_);
		cdmatrix eps_b(size_,size_);

		C00_bos=GREEN(nt,ntau,1,1);
		C01_bos=GREEN(nt,ntau,1,1);
		C10_bos=GREEN(nt,ntau,1,1);
		C11_bos=GREEN(nt,ntau,1,1);
		Ctemp_bos=GREEN(nt,ntau,1,1);

		A00_bos=GREEN(nt,ntau,1,1);
		A01_bos=GREEN(nt,ntau,1,1);
		A10_bos=GREEN(nt,ntau,1,1);
		A11_bos=GREEN(nt,ntau,1,1);

		B00_bos=GREEN(nt,ntau,1,1);
		B01_bos=GREEN(nt,ntau,1,1);
		B10_bos=GREEN(nt,ntau,1,1);
		B11_bos=GREEN(nt,ntau,1,1);

		AB00_bos=GREEN(nt,ntau,1,1);
		AB01_bos=GREEN(nt,ntau,1,1);
		AB10_bos=GREEN(nt,ntau,1,1);
		AB11_bos=GREEN(nt,ntau,1,1);

		eps_a.setZero();
		eps_b.setZero();

		eps_a(0,0)=1.123;
    eps_a(0,1)=0.1;
    eps_a(1,0)=0.1;
    eps_a(1,1)=0.567;

    eps_b(0,0)=0.345;
    eps_b(0,1)=0.2;
    eps_b(1,0)=0.2;
    eps_b(1,1)=0.876;

		cntr::green_from_H(A_bos,mu,eps_a,beta,h);
		cntr::green_from_H(B_bos,mu,eps_b,beta,h);

		A00_bos.set_matrixelement(0,0,A_bos,0,0);
		A01_bos.set_matrixelement(0,0,A_bos,0,1);
		A10_bos.set_matrixelement(0,0,A_bos,1,0);
		A11_bos.set_matrixelement(0,0,A_bos,1,1);

		B00_bos.set_matrixelement(0,0,B_bos,0,0);
		B01_bos.set_matrixelement(0,0,B_bos,0,1);
		B10_bos.set_matrixelement(0,0,B_bos,1,0);
		B11_bos.set_matrixelement(0,0,B_bos,1,1);

		//cout << "SIMPLE TEST WITH GIVEN h=0.02, dtau=0.005" << endl;

		cntr::convolution(C00_bos,A00_bos,A00_bos,B00_bos,B00_bos,integration::I<double>(kt),beta,h);
		cntr::convolution(Ctemp_bos,A01_bos,A01_bos,B10_bos,B10_bos,integration::I<double>(kt),beta,h);
		for(tstp=-1;tstp<=nt;tstp++){
				C00_bos.incr_timestep(tstp,Ctemp_bos,CPLX(1.0,0.0));
		}
		cntr::convolution(C01_bos,A00_bos,A00_bos,B01_bos,B01_bos,integration::I<double>(kt),beta,h);
		cntr::convolution(Ctemp_bos,A01_bos,A01_bos,B11_bos,B11_bos,integration::I<double>(kt),beta,h);
		for(tstp=-1;tstp<=nt;tstp++){
				C01_bos.incr_timestep(tstp,Ctemp_bos,CPLX(1.0,0.0));
		}
		cntr::convolution(C10_bos,A10_bos,A10_bos,B00_bos,B00_bos,integration::I<double>(kt),beta,h);
		cntr::convolution(Ctemp_bos,A11_bos,A11_bos,B10_bos,B10_bos,integration::I<double>(kt),beta,h);
		for(tstp=-1;tstp<=nt;tstp++){
				C10_bos.incr_timestep(tstp,Ctemp_bos,CPLX(1.0,0.0));
		}
		cntr::convolution(C11_bos,A10_bos,A10_bos,B01_bos,B01_bos,integration::I<double>(kt),beta,h);
		cntr::convolution(Ctemp_bos,A11_bos,A11_bos,B11_bos,B11_bos,integration::I<double>(kt),beta,h);
		for(tstp=-1;tstp<=nt;tstp++){
				C11_bos.incr_timestep(tstp,Ctemp_bos,CPLX(1.0,0.0));
		}

		cntr::convolution(AB_bos,A_bos,A_bos,B_bos,B_bos,integration::I<double>(kt),beta,h);

		AB00_bos.set_matrixelement(0,0,AB_bos,0,0);
		AB01_bos.set_matrixelement(0,0,AB_bos,0,1);
		AB10_bos.set_matrixelement(0,0,AB_bos,1,0);
		AB11_bos.set_matrixelement(0,0,AB_bos,1,1);

		//cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
		double err=0.0;
		for(tstp=-1;tstp<=nt;tstp++){
			err+=cntr::distance_norm2(tstp,C00_bos,AB00_bos);
			err+=cntr::distance_norm2(tstp,C01_bos,AB01_bos);
			err+=cntr::distance_norm2(tstp,C10_bos,AB10_bos);
			err+=cntr::distance_norm2(tstp,C11_bos,AB11_bos);

			err+=cntr::distance_norm2_ret(tstp,C00_bos,AB00_bos);
			err+=cntr::distance_norm2_ret(tstp,C01_bos,AB01_bos);
			err+=cntr::distance_norm2_ret(tstp,C10_bos,AB10_bos);
			err+=cntr::distance_norm2_ret(tstp,C11_bos,AB11_bos);

			err+=cntr::distance_norm2_tv(tstp,C00_bos,AB00_bos);
			err+=cntr::distance_norm2_tv(tstp,C01_bos,AB01_bos);
			err+=cntr::distance_norm2_tv(tstp,C10_bos,AB10_bos);
			err+=cntr::distance_norm2_tv(tstp,C11_bos,AB11_bos);

			err+=cntr::distance_norm2_les(tstp,C00_bos,AB00_bos);
			err+=cntr::distance_norm2_les(tstp,C01_bos,AB01_bos);
			err+=cntr::distance_norm2_les(tstp,C10_bos,AB10_bos);
			err+=cntr::distance_norm2_les(tstp,C11_bos,AB11_bos);

		}

    REQUIRE(err<eps);
  }

  // TEST REAL-TIME ... SET MATSUBARA TO EXACT
  SECTION("convolution (bosonic signs in boundary terms)"){
    double tmax=2.0;
    int nt0=10,fac2;
    beta=0.1;
    double err=0.0;

    for(fac2=1;fac2<=8;fac2*=2){
      nt=nt0*fac2;

      A_bos=GREEN(nt,ntau,size_,1);
      B_bos=GREEN(nt,ntau,size_,1);
      C_bos=GREEN(nt,ntau,size_,1);
      AB_bos=GREEN(nt,ntau,size_,1);

      h=tmax/nt;

			cdmatrix eps_a(size_,size_);
			cdmatrix eps_b(size_,size_);

			C00_bos=GREEN(nt,ntau,1,1);
			C01_bos=GREEN(nt,ntau,1,1);
			C10_bos=GREEN(nt,ntau,1,1);
			C11_bos=GREEN(nt,ntau,1,1);
			Ctemp_bos=GREEN(nt,ntau,1,1);

			A00_bos=GREEN(nt,ntau,1,1);
			A01_bos=GREEN(nt,ntau,1,1);
			A10_bos=GREEN(nt,ntau,1,1);
			A11_bos=GREEN(nt,ntau,1,1);

			B00_bos=GREEN(nt,ntau,1,1);
			B01_bos=GREEN(nt,ntau,1,1);
			B10_bos=GREEN(nt,ntau,1,1);
			B11_bos=GREEN(nt,ntau,1,1);

			AB00_bos=GREEN(nt,ntau,1,1);
			AB01_bos=GREEN(nt,ntau,1,1);
			AB10_bos=GREEN(nt,ntau,1,1);
			AB11_bos=GREEN(nt,ntau,1,1);

			eps_a.setZero();
			eps_b.setZero();

			eps_a(0,0)=1.123;
			eps_a(0,1)=0.1;
			eps_a(1,0)=0.1;
			eps_a(1,1)=0.567;

			eps_b(0,0)=0.345;
			eps_b(0,1)=0.2;
			eps_b(1,0)=0.2;
			eps_b(1,1)=0.876;

			cntr::green_from_H(A_bos,mu,eps_a,beta,h);
			cntr::green_from_H(B_bos,mu,eps_b,beta,h);

			A00_bos.set_matrixelement(0,0,A_bos,0,0);
			A01_bos.set_matrixelement(0,0,A_bos,0,1);
			A10_bos.set_matrixelement(0,0,A_bos,1,0);
			A11_bos.set_matrixelement(0,0,A_bos,1,1);

			B00_bos.set_matrixelement(0,0,B_bos,0,0);
			B01_bos.set_matrixelement(0,0,B_bos,0,1);
			B10_bos.set_matrixelement(0,0,B_bos,1,0);
			B11_bos.set_matrixelement(0,0,B_bos,1,1);

			//cout << "SIMPLE TEST WITH GIVEN h=0.02, dtau=0.005" << endl;

			cntr::convolution(C00_bos,A00_bos,A00_bos,B00_bos,B00_bos,integration::I<double>(kt),beta,h);
			cntr::convolution(Ctemp_bos,A01_bos,A01_bos,B10_bos,B10_bos,integration::I<double>(kt),beta,h);
			for(tstp=-1;tstp<=nt;tstp++){
					C00_bos.incr_timestep(tstp,Ctemp_bos,CPLX(1.0,0.0));
			}
			cntr::convolution(C01_bos,A00_bos,A00_bos,B01_bos,B01_bos,integration::I<double>(kt),beta,h);
			cntr::convolution(Ctemp_bos,A01_bos,A01_bos,B11_bos,B11_bos,integration::I<double>(kt),beta,h);
			for(tstp=-1;tstp<=nt;tstp++){
					C01_bos.incr_timestep(tstp,Ctemp_bos,CPLX(1.0,0.0));
			}
			cntr::convolution(C10_bos,A10_bos,A10_bos,B00_bos,B00_bos,integration::I<double>(kt),beta,h);
			cntr::convolution(Ctemp_bos,A11_bos,A11_bos,B10_bos,B10_bos,integration::I<double>(kt),beta,h);
			for(tstp=-1;tstp<=nt;tstp++){
					C10_bos.incr_timestep(tstp,Ctemp_bos,CPLX(1.0,0.0));
			}
			cntr::convolution(C11_bos,A10_bos,A10_bos,B01_bos,B01_bos,integration::I<double>(kt),beta,h);
			cntr::convolution(Ctemp_bos,A11_bos,A11_bos,B11_bos,B11_bos,integration::I<double>(kt),beta,h);
			for(tstp=-1;tstp<=nt;tstp++){
					C11_bos.incr_timestep(tstp,Ctemp_bos,CPLX(1.0,0.0));
			}

			cntr::convolution(AB_bos,A_bos,A_bos,B_bos,B_bos,integration::I<double>(kt),beta,h);

			AB00_bos.set_matrixelement(0,0,AB_bos,0,0);
			AB01_bos.set_matrixelement(0,0,AB_bos,0,1);
			AB10_bos.set_matrixelement(0,0,AB_bos,1,0);
			AB11_bos.set_matrixelement(0,0,AB_bos,1,1);

			//cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
			double err=0.0;
			for(tstp=-1;tstp<=nt;tstp++){
				err+=cntr::distance_norm2(tstp,C00_bos,AB00_bos);
				err+=cntr::distance_norm2(tstp,C01_bos,AB01_bos);
				err+=cntr::distance_norm2(tstp,C10_bos,AB10_bos);
				err+=cntr::distance_norm2(tstp,C11_bos,AB11_bos);

				err+=cntr::distance_norm2_ret(tstp,C00_bos,AB00_bos);
				err+=cntr::distance_norm2_ret(tstp,C01_bos,AB01_bos);
				err+=cntr::distance_norm2_ret(tstp,C10_bos,AB10_bos);
				err+=cntr::distance_norm2_ret(tstp,C11_bos,AB11_bos);

				err+=cntr::distance_norm2_tv(tstp,C00_bos,AB00_bos);
				err+=cntr::distance_norm2_tv(tstp,C01_bos,AB01_bos);
				err+=cntr::distance_norm2_tv(tstp,C10_bos,AB10_bos);
				err+=cntr::distance_norm2_tv(tstp,C11_bos,AB11_bos);

				err+=cntr::distance_norm2_les(tstp,C00_bos,AB00_bos);
				err+=cntr::distance_norm2_les(tstp,C01_bos,AB01_bos);
				err+=cntr::distance_norm2_les(tstp,C10_bos,AB10_bos);
				err+=cntr::distance_norm2_les(tstp,C11_bos,AB11_bos);

			}

    }
    REQUIRE(err<eps);
  }

	SECTION("convolution (fermi)"){

		A_fer=GREEN(nt,ntau,size_,-1);
		B_fer=GREEN(nt,ntau,size_,-1);
		C_fer=GREEN(nt,ntau,size_,-1);
		AB_fer=GREEN(nt,ntau,size_,-1);

		cdmatrix eps_a(size_,size_);
		cdmatrix eps_b(size_,size_);

		C00_fer=GREEN(nt,ntau,1,-1);
		C01_fer=GREEN(nt,ntau,1,-1);
		C10_fer=GREEN(nt,ntau,1,-1);
		C11_fer=GREEN(nt,ntau,1,-1);
		Ctemp_fer=GREEN(nt,ntau,1,-1);

		A00_fer=GREEN(nt,ntau,1,-1);
		A01_fer=GREEN(nt,ntau,1,-1);
		A10_fer=GREEN(nt,ntau,1,-1);
		A11_fer=GREEN(nt,ntau,1,-1);

		B00_fer=GREEN(nt,ntau,1,-1);
		B01_fer=GREEN(nt,ntau,1,-1);
		B10_fer=GREEN(nt,ntau,1,-1);
		B11_fer=GREEN(nt,ntau,1,-1);

		AB00_fer=GREEN(nt,ntau,1,-1);
		AB01_fer=GREEN(nt,ntau,1,-1);
		AB10_fer=GREEN(nt,ntau,1,-1);
		AB11_fer=GREEN(nt,ntau,1,-1);

		eps_a.setZero();
		eps_b.setZero();

		eps_a(0,0)=1.123;
    eps_a(0,1)=0.1;
    eps_a(1,0)=0.1;
    eps_a(1,1)=0.567;

    eps_b(0,0)=0.345;
    eps_b(0,1)=0.2;//0.2;
    eps_b(1,0)=0.2;
    eps_b(1,1)=0.876;

		cntr::green_from_H(A_fer,mu,eps_a,beta,h);
		cntr::green_from_H(B_fer,mu,eps_b,beta,h);

		A00_fer.set_matrixelement(0,0,A_fer,0,0);
		A01_fer.set_matrixelement(0,0,A_fer,0,1);
		A10_fer.set_matrixelement(0,0,A_fer,1,0);
		A11_fer.set_matrixelement(0,0,A_fer,1,1);

		B00_fer.set_matrixelement(0,0,B_fer,0,0);
		B01_fer.set_matrixelement(0,0,B_fer,0,1);
		B10_fer.set_matrixelement(0,0,B_fer,1,0);
		B11_fer.set_matrixelement(0,0,B_fer,1,1);

		//cout << "SIMPLE TEST WITH GIVEN h=0.02, dtau=0.005" << endl;

		cntr::convolution(C00_fer,A00_fer,A00_fer,B00_fer,B00_fer,integration::I<double>(kt),beta,h);
		cntr::convolution(Ctemp_fer,A01_fer,A01_fer,B10_fer,B10_fer,integration::I<double>(kt),beta,h);
		for(tstp=-1;tstp<=nt;tstp++){
				C00_fer.incr_timestep(tstp,Ctemp_fer,CPLX(1.0,0.0));
		}
		cntr::convolution(C01_fer,A00_fer,A00_fer,B01_fer,B01_fer,integration::I<double>(kt),beta,h);
		cntr::convolution(Ctemp_fer,A01_fer,A01_fer,B11_fer,B11_fer,integration::I<double>(kt),beta,h);
		for(tstp=-1;tstp<=nt;tstp++){
				C01_fer.incr_timestep(tstp,Ctemp_fer,CPLX(1.0,0.0));
		}
		cntr::convolution(C10_fer,A10_fer,A10_fer,B00_fer,B00_fer,integration::I<double>(kt),beta,h);
		cntr::convolution(Ctemp_fer,A11_fer,A11_fer,B10_fer,B10_fer,integration::I<double>(kt),beta,h);
		for(tstp=-1;tstp<=nt;tstp++){
				C10_fer.incr_timestep(tstp,Ctemp_fer,CPLX(1.0,0.0));
		}
		cntr::convolution(C11_fer,A10_fer,A10_fer,B01_fer,B01_fer,integration::I<double>(kt),beta,h);
		cntr::convolution(Ctemp_fer,A11_fer,A11_fer,B11_fer,B11_fer,integration::I<double>(kt),beta,h);
		for(tstp=-1;tstp<=nt;tstp++){
				C11_fer.incr_timestep(tstp,Ctemp_fer,CPLX(1.0,0.0));
		}

		cntr::convolution(AB_fer,A_fer,A_fer,B_fer,B_fer,integration::I<double>(kt),beta,h);

		AB00_fer.set_matrixelement(0,0,AB_fer,0,0);
		AB01_fer.set_matrixelement(0,0,AB_fer,0,1);
		AB10_fer.set_matrixelement(0,0,AB_fer,1,0);
		AB11_fer.set_matrixelement(0,0,AB_fer,1,1);

		//cntr::convolution(C,B,B,A,A,integration::I<double>(kt),beta,h);
		double err=0.0;
		for(tstp=-1;tstp<=nt;tstp++){
			err+=cntr::distance_norm2(tstp,C00_fer,AB00_fer);
			err+=cntr::distance_norm2(tstp,C01_fer,AB01_fer);
			err+=cntr::distance_norm2(tstp,C10_fer,AB10_fer);
			err+=cntr::distance_norm2(tstp,C11_fer,AB11_fer);

			err+=cntr::distance_norm2_ret(tstp,C00_fer,AB00_fer);
			err+=cntr::distance_norm2_ret(tstp,C01_fer,AB01_fer);
			err+=cntr::distance_norm2_ret(tstp,C10_fer,AB10_fer);
			err+=cntr::distance_norm2_ret(tstp,C11_fer,AB11_fer);

			err+=cntr::distance_norm2_tv(tstp,C00_fer,AB00_fer);
			err+=cntr::distance_norm2_tv(tstp,C01_fer,AB01_fer);
			err+=cntr::distance_norm2_tv(tstp,C10_fer,AB10_fer);
			err+=cntr::distance_norm2_tv(tstp,C11_fer,AB11_fer);

			err+=cntr::distance_norm2_les(tstp,C00_fer,AB00_fer);
			err+=cntr::distance_norm2_les(tstp,C01_fer,AB01_fer);
			err+=cntr::distance_norm2_les(tstp,C10_fer,AB10_fer);
			err+=cntr::distance_norm2_les(tstp,C11_fer,AB11_fer);

		}
		//err+=A_fer.sig();

		REQUIRE(err<eps);
	}

  SECTION("time-diagonal convolution") {
    CPLX I(0.0,1.0);
    int size1=1;
    double err=0.0,err_tot=0.0;
    A_fer=GREEN(nt,ntau,1,-1);
    B_fer=GREEN(nt,ntau,1,-1);
    AB_fer=GREEN(nt,ntau,1,-1);

    cdmatrix eps_a(1,1);
    cdmatrix eps_b(1,1);

    eps_a.setZero();
    eps_b.setZero();

    eps_a(0,0)=wa;
    eps_b(0,0)=wb;

    cntr::green_from_H(A_fer,mu,eps_a,beta,h);
    cntr::green_from_H(B_fer,mu,eps_b,beta,h);

    cntr::convolution(AB_fer,A_fer,A_fer,B_fer,B_fer,integration::I<double>(kt),beta,h);

    complex<double> *AB_diag = new complex<double>[size1*size1];

    for(tstp=-1; tstp<= nt; tstp++){
      cntr::convolution_density_matrix(tstp,AB_diag,A_fer,B_fer,integration::I<double>(kt),beta,h);
      cdmatrix C_diag(size1,size1);
      if(tstp == -1){
				AB_fer.get_mat(0,C_diag);
				err = fabs(C_diag(0,0) - AB_diag[0]);
      } else {
				AB_fer.get_les(tstp,tstp,C_diag);
				err = fabs(C_diag(0,0) - I*AB_diag[0]);
      }
      err_tot += err;
    }
    REQUIRE(err<eps);
  }

  SECTION("correlation energy"){
    const double lam=0.4,eps1=-1.0,eps2=0.2,mu=0.0;
    std::complex<double> I(0.0,1.0);
    cdmatrix h2x2(2,2),h2x2_diag(2,2),h11(1,1),h22(1,1);
    GREEN G2x2(nt,ntau,2,-1);
    GREEN G1x1(nt,ntau,1,-1);
    GREEN SGM(nt,ntau,1,-1);
    cntr::function<double> h11_func(nt,1);
    double err=0.0;

    h2x2(0,0) = eps1;
    h2x2(0,1) = I*lam;
    h2x2(1,0) = -I*lam;
    h2x2(1,1) = eps2;
    h2x2_diag(0,0) = eps1;
    h2x2_diag(1,1) = eps2;
    h2x2_diag(0,1) = 0.0;
    h2x2_diag(1,0) = 0.0;

    h11(0,0) = eps1;
    h22(0,0) = eps2;

    h11_func.set_constant(h11);

    cntr::green_from_H(SGM,mu,h22,beta,h);
    for(tstp=-1;tstp<=nt;tstp++)
      SGM.smul(tstp,lam*lam);

    cntr::dyson_mat(G1x1, SGM, mu, h11_func, integration::I<double>(kt),
			     beta, CNTR_MAT_FIXPOINT);
    cntr::dyson_start(G1x1,mu,h11_func,SGM,integration::I<double>(kt),beta,h);
    for(tstp=kt+1;tstp<=nt;tstp++){
      cntr::dyson_timestep(tstp, G1x1, mu, h11_func, SGM, integration::I<double>(kt), beta, h);
    }

    cntr::green_from_H(G2x2,mu,h2x2,beta,h);
    cdmatrix rho2x2(2,2);
    double E2x2,E2x2_diag,Ecorr1,Ecorr2;
    for(tstp=-1;tstp<=nt;tstp++){
      G2x2.density_matrix(tstp,rho2x2);
      E2x2 = ((h2x2*rho2x2).trace()).real();
      E2x2_diag = ((h2x2_diag*rho2x2).trace()).real();
      Ecorr1 = cntr::correlation_energy(tstp, G1x1, SGM, integration::I<double>(kt), beta, h);
      Ecorr2 = 0.25*(E2x2 - E2x2_diag);
      err += fabs(Ecorr1-Ecorr2);
      //cout << tstp << "  " << Ecorr1 << "  " <<  Ecorr2 <<   endl;
    }
    REQUIRE(err<eps);
  }

#if CNTR_USE_OMP==1
  SECTION("convolution (openmp)"){

    A_fer=GREEN(nt,ntau,size_,-1);
    B_fer=GREEN(nt,ntau,size_,-1);
    C_fer=GREEN(nt,ntau,size_,-1);
    GREEN C_fer_omp(nt,ntau,size_,-1);

    cdmatrix eps_a(size_,size_);
    cdmatrix eps_b(size_,size_);

    std::complex<double> I(0.0,1.0);
    int nthreads;

#pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }

    eps_a.setZero();
    eps_b.setZero();

    eps_a(0,0)=1.123;
    eps_a(0,1)=0.1;
    eps_a(1,0)=0.1;
    eps_a(1,1)=0.567;

    eps_b(0,0)=0.345;
    eps_b(0,1)=0.2;
    eps_b(1,0)=0.2;
    eps_b(1,1)=0.876;

    cntr::green_from_H(A_fer,mu,eps_a,beta,h);
    cntr::green_from_H(B_fer,mu,eps_b,beta,h);

    cntr::convolution_matsubara(C_fer,A_fer,B_fer,integration::I<double>(kt),beta);


    cntr::convolution_matsubara_nomp(nthreads,C_fer_omp,A_fer,B_fer,integration::I<double>(kt),beta);

    double err = cntr::distance_norm2(-1,C_fer,C_fer_omp);

    REQUIRE(err<eps);
  }
#endif
}
