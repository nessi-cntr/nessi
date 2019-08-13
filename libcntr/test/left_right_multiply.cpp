#include "catch.hpp"
#include "cntr.hpp"
#include <cmath>

using namespace std;
#define GREEN cntr::herm_matrix<double>
#define CFUNC cntr::function<double>

TEST_CASE("left and right multiply","[Herm_matrix_left_right_multiply]"){
  int size=2;
  int nt=100, ntau=50;
  double eps=1e-6;
  double h=0.01, mu=0.0, beta=10.0;
  double tmax=h*nt;
  double eps1=-0.4,eps2=0.6;
  double pi=4.0*atan(1.0);
  std::complex<double> I(0.0,1.0);
  cdmatrix h0(2,2),vmat(size,size);
  CFUNC Vfunc(nt,size);
  GREEN G0(nt,ntau,size,-1);
  GREEN G1(nt,ntau,size,-1);
  
  h0(0,0) = eps1;
  h0(1,1) = eps2;
  h0(0,1) = 0.0;
  h0(1,0) = 0.0;

  vmat(0,0) = 1.0;
  vmat(0,1) = 0.0;
  vmat(1,0) = 0.0;
  vmat(1,1) = 1.0;

  Vfunc.set_value(-1,vmat);
  
  for(int i=0; i<=nt; i++){
    vmat(0,0) = cos(pi*h*i/tmax);
    vmat(0,1) = -sin(pi*h*i/tmax);
    vmat(1,0) = sin(pi*h*i/tmax);
    vmat(1,1) = cos(pi*h*i/tmax);
    Vfunc.set_value(i,vmat);
  }
  
  cntr::green_from_H(G0,mu,h0,beta,h,5,4,true);
  cntr::green_from_H(G1,mu,h0,beta,h,5,4,true);

  SECTION ("left multiply"){
    double err;
    cdmatrix vmat(2,2);
    cdmatrix mat1(2,2),ret1(2,2),les1(2,2),tv1(2,2);
    cdmatrix mat2(2,2),ret2(2,2),les2(2,2),tv2(2,2);
  
    for(int tstp=-1; tstp<=nt; tstp++){
      G1.left_multiply(tstp, Vfunc, 1.0);
    }

    err = 0.0;

    Vfunc.get_value(-1,vmat);
    for(int itau=0; itau<=ntau; itau++){
      G0.get_mat(itau,mat1);
      G1.get_mat(itau,mat2);
      mat1 = vmat * mat1 - mat2;
      err += mat1.norm();
    }

    for(int i=0; i<=nt; i++){

      for(int j=0; j<=i; j++){
	Vfunc.get_value(j,vmat);
      	G0.get_les(j,i,les1);
      	G1.get_les(j,i,les2);
      	les1 = vmat*les1 - les2;
      	err += les1.norm();

	Vfunc.get_value(i,vmat);
      	G0.get_ret(i,j,ret1);
      	G1.get_ret(i,j,ret2);
      	ret1 = vmat*ret1 - ret2;
      	err += ret1.norm();
	
      }
      for(int itau=0; itau<=ntau; itau++){
    	G0.get_tv(i,itau,tv1);
    	G1.get_tv(i,itau,tv2);
    	tv1 = vmat*tv1 - tv2;
    	err += tv1.norm();
      }
    }

    REQUIRE(err<eps);
    
  }

  SECTION ("left hermconjg multiply"){
    double err;
    cdmatrix vmat(2,2);
    cdmatrix mat1(2,2),ret1(2,2),les1(2,2),tv1(2,2);
    cdmatrix mat2(2,2),ret2(2,2),les2(2,2),tv2(2,2);

    cntr::green_from_H(G1,mu,h0,beta,h,5,4,true);
    
    for(int tstp=-1; tstp<=nt; tstp++){
      G1.left_multiply_hermconj(tstp, Vfunc, 1.0);
    }

    err = 0.0;

    Vfunc.get_value(-1,vmat);
    for(int itau=0; itau<=ntau; itau++){
      G0.get_mat(itau,mat1);
      G1.get_mat(itau,mat2);
      mat1 = vmat.adjoint() * mat1 - mat2;
      err += mat1.norm();
    }

    for(int i=0; i<=nt; i++){
      for(int j=0; j<=i; j++){
	Vfunc.get_value(j,vmat);
	G0.get_les(j,i,les1);
	G1.get_les(j,i,les2);
	les1 = vmat.adjoint() * les1 - les2;
	err += les1.norm();

	Vfunc.get_value(i,vmat);
	G0.get_ret(i,j,ret1);
	G1.get_ret(i,j,ret2);
	ret1 = vmat.adjoint() *ret1 - ret2;
	err += ret1.norm();	
      }
      Vfunc.get_value(i,vmat);
      for(int itau=0; itau<=ntau; itau++){
	G0.get_tv(i,itau,tv1);
	G1.get_tv(i,itau,tv2);
	tv1 = vmat.adjoint() *tv1 - tv2;
	err += tv1.norm();
      }
    }

    REQUIRE(err<eps);
    
  }
  
  SECTION ("right multiply"){
    double err;
    cdmatrix vmat(2,2);
    cdmatrix mat1(2,2),ret1(2,2),les1(2,2),tv1(2,2);
    cdmatrix mat2(2,2),ret2(2,2),les2(2,2),tv2(2,2);

    cntr::green_from_H(G1,mu,h0,beta,h,5,4,true);
    
    for(int tstp=-1; tstp<=nt; tstp++){
      G1.right_multiply(tstp, Vfunc, 1.0);
    }

    err = 0.0;

    Vfunc.get_value(-1,vmat);
    for(int itau=0; itau<=ntau; itau++){
      G0.get_mat(itau,mat1);
      G1.get_mat(itau,mat2);
      mat1 = mat1 * vmat - mat2;
      err += mat1.norm();
    }

    for(int i=0; i<=nt; i++){
      for(int j=0; j<=i; j++){
      	Vfunc.get_value(i,vmat);
      	G0.get_les(j,i,les1);
      	G1.get_les(j,i,les2);
      	les1 = les1 * vmat - les2;
      	err += les1.norm();

      	Vfunc.get_value(j,vmat);
      	G0.get_ret(i,j,ret1);
      	G1.get_ret(i,j,ret2);
      	ret1 = ret1 * vmat - ret2;
      	err += ret1.norm();	
      }
      
      Vfunc.get_value(-1,vmat);
      for(int itau=0; itau<=ntau; itau++){
      	G0.get_tv(i,itau,tv1);
      	G1.get_tv(i,itau,tv2);
      	tv1 = tv1 * vmat - tv2;
      	err += tv1.norm();
      }
    }

    REQUIRE(err<eps);
    
  }
  
  SECTION ("right hermconj multiply"){
    double err;
    cdmatrix vmat(2,2);
    cdmatrix mat1(2,2),ret1(2,2),les1(2,2),tv1(2,2);
    cdmatrix mat2(2,2),ret2(2,2),les2(2,2),tv2(2,2);

    cntr::green_from_H(G1,mu,h0,beta,h,5,4,true);
    
    for(int tstp=-1; tstp<=nt; tstp++){
      G1.right_multiply_hermconj(tstp, Vfunc, 1.0);
    }

    err = 0.0;

    Vfunc.get_value(-1,vmat);
    for(int itau=0; itau<=ntau; itau++){
      G0.get_mat(itau,mat1);
      G1.get_mat(itau,mat2);
      mat1 = mat1 * vmat.adjoint() - mat2;
      err += mat1.norm();
    }

    for(int i=0; i<=nt; i++){
      for(int j=0; j<=i; j++){
      	Vfunc.get_value(i,vmat);
      	G0.get_les(j,i,les1);
      	G1.get_les(j,i,les2);
      	les1 = les1 * vmat.adjoint() - les2;
      	err += les1.norm();

      	Vfunc.get_value(j,vmat);
      	G0.get_ret(i,j,ret1);
      	G1.get_ret(i,j,ret2);
      	ret1 = ret1 * vmat.adjoint() - ret2;
      	err += ret1.norm();	
      }
      
      Vfunc.get_value(-1,vmat);
      for(int itau=0; itau<=ntau; itau++){
      	G0.get_tv(i,itau,tv1);
      	G1.get_tv(i,itau,tv2);
      	tv1 = tv1 * vmat.adjoint() - tv2;
      	err += tv1.norm();
      }
    }

    REQUIRE(err<eps);
    
  }
  
}
