#include "catch.hpp"
#include "cntr.hpp"


using namespace std;
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>


TEST_CASE("herm_matrix algebra","[herm_matrix_algebra]"){
  int size=2;
  int nt=100, ntau=50;
  double eps=1e-6;
  double dt=0.01, mu=0.0, beta=10.0;
  double tmax=dt*nt;
  double eps1=-0.4,eps2=0.6,lam1=0.1;
  double eps3=0.435,eps4=0.5676,lam2=0.1566;
  double wr=0.3;
  std::complex<double> I(0.0,1.0);
  std::complex<double> wz=(1.0,-0.3);
  cdmatrix h1(2,2);
  cdmatrix h2(2,2);
  GREEN G1(nt,ntau,size,-1);
  GREEN G2(nt,ntau,size,-1);
  GREEN G3(nt,ntau,size,-1);
  GREEN G4(nt,ntau,size,-1);
  
  h1(0,0) = eps1;
  h1(1,1) = eps2;
  h1(0,1) = I*lam1;
  h1(1,0) = -I*lam1;

  h2(0,0) = eps3;
  h2(1,1) = eps4;
  h2(0,1) = I*lam2;
  h2(1,0) = -I*lam2;

  cntr::green_from_H(G1,mu,h1,beta,dt);
  cntr::green_from_H(G2,mu,h2,beta,dt);

  SECTION("incr_timestep"){
    double err=0.0;
    cdmatrix mat1(2,2),tv1(2,2),les1(2,2),ret1(2,2);
    cdmatrix mat2(2,2),tv2(2,2),les2(2,2),ret2(2,2);
    cdmatrix mat3(2,2),tv3(2,2),les3(2,2),ret3(2,2);
    
    for(int q; q<=ntau; q++){
      G1.get_mat(q,mat1);
      G2.get_mat(q,mat2);
      mat3 = mat1 + wz*mat2;
      G3.set_mat(q,mat3);
    }
  
    for(int i=0; i<=nt; i++){
      for(int j=0; j<=i; j++){
	G1.get_les(j,i,les1);
	G2.get_les(j,i,les2);
	les3 = les1 + wz*les2;
	G3.set_les(j,i,les3);

	G1.get_ret(i,j,ret1);
	G2.get_ret(i,j,ret2);
	ret3 = ret1 + wz*ret2;
	G3.set_ret(i,j,ret3);	
      }
      for(int q=0; q<=ntau; q++){
	G1.get_tv(i,q,tv1);
	G2.get_tv(i,q,tv2);
	tv3 = tv1 + wz*tv2;
	G3.set_tv(i,q,tv3);
      }
    }

    cntr::green_from_H(G4,mu,h1,beta,dt);
    err = 0.0;
    G4.incr_timestep(G2,wz);
    for(int tstp=-1; tstp<=nt; tstp++){
      err += cntr::distance_norm2(tstp,G3,G4);
    }
    REQUIRE(err<eps);   
    
    cntr::green_from_H(G4,mu,h1,beta,dt);  
    err = 0.0;  
    for(int tstp=-1; tstp<=nt; tstp++){
      GREEN_TSTP A(tstp,ntau,size);
      G2.get_timestep(tstp,A);
      G4.incr_timestep(tstp,A,wz);
      err += cntr::distance_norm2(tstp,G3,G4);
    }
    REQUIRE(err<eps);

    cntr::green_from_H(G4,mu,h1,beta,dt);  
    err = 0.0;  
    for(int tstp=-1; tstp<=nt; tstp++){
      G4.incr_timestep(tstp,G2,wz);
      err += cntr::distance_norm2(tstp,G3,G4);
    }
    REQUIRE(err<eps);
 
  }

  SECTION("smul (complex weight)"){
    double err=0.0;
    cdmatrix mat1(2,2),tv1(2,2),les1(2,2),ret1(2,2);
    cdmatrix mat3(2,2),tv3(2,2),les3(2,2),ret3(2,2);

    for(int q; q<=ntau; q++){
      G1.get_mat(q,mat1);
      mat3 = wz*mat1;
      G3.set_mat(q,mat3);
    }
  
    for(int i=0; i<=nt; i++){
      for(int j=0; j<=i; j++){
	G1.get_les(j,i,les1);
	les3 = wz*les1;
	G3.set_les(j,i,les3);

	G1.get_ret(i,j,ret1);
	ret3 = wz*ret1;
	G3.set_ret(i,j,ret3);	
      }
      for(int q=0; q<=ntau; q++){
	G1.get_tv(i,q,tv1);
	tv3 = wz*tv1;
	G3.set_tv(i,q,tv3);
      }
    }

    err = 0.0;
    cntr::green_from_H(G4,mu,h1,beta,dt);
    for(int tstp=-1; tstp<=nt; tstp++){
      G4.smul(tstp,wz);
      err += cntr::distance_norm2(tstp,G3,G4);
    }
    REQUIRE(err<eps);
  }

  SECTION("smul (real weight)"){
    double err=0.0;
    cdmatrix mat1(2,2),tv1(2,2),les1(2,2),ret1(2,2);
    cdmatrix mat3(2,2),tv3(2,2),les3(2,2),ret3(2,2);
    
    for(int q; q<=ntau; q++){
      G1.get_mat(q,mat1);
      mat3 = wr*mat1;
      G3.set_mat(q,mat3);
    }
  
    for(int i=0; i<=nt; i++){
      for(int j=0; j<=i; j++){
	G1.get_les(j,i,les1);
	les3 = wr*les1;
	G3.set_les(j,i,les3);

	G1.get_ret(i,j,ret1);
	ret3 = wr*ret1;
	G3.set_ret(i,j,ret3);	
      }
      for(int q=0; q<=ntau; q++){
	G1.get_tv(i,q,tv1);
	tv3 = wr*tv1;
	G3.set_tv(i,q,tv3);
      }
    }

    err = 0.0;
    cntr::green_from_H(G4,mu,h1,beta,dt);
    for(int tstp=-1; tstp<=nt; tstp++){
      G4.smul(tstp,wr);
      err += cntr::distance_norm2(tstp,G3,G4);
    }
    REQUIRE(err<eps);
  }
}

