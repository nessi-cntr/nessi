#include "catch.hpp"
#include "cntr.hpp"

// TODO:
// Set matrixelement
// Incr
// Right/left multiply
// Interaction with herm_matrix_timestep 

#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define cfunction cntr::function<double>

// Test get_set for herm_matrix,herm_matrix_timestep,cfunction
double setget(GREEN &A,cdmatrix &a){
	double toterr=0.0;
	// Matsubara
	for(int m=0;m<A.ntau();m++){
		cdmatrix tmp;
		A.set_mat(m,a);
		A.get_mat(m,tmp);
		toterr+=(tmp-a).norm();
		// tv
		for(int n=0;n<A.nt();n++){
			A.set_tv(n,m,a);
			A.get_tv(n,m,a);
			toterr+=(tmp-a).norm();
		}
	}

	for(int m=0;m<A.nt();m++){
		for(int n=0;n<m;n++){
			cdmatrix ret,les;
			// Ret
			A.set_ret(m,n,a);
			A.get_ret(m,n,ret);
			toterr+=(ret-a).norm();
			// Les
			A.set_les(n,m,a);
			A.get_les(n,m,les);
			toterr+=(les-a).norm();
		}
	}
	return toterr;
}


void exact_rightmultiply(double beta,double dt,GREEN &G){
	int ntau=G.ntau();
	int nt=G.nt();
	int size=G.size1();
	assert(size==2);

	double dtau=beta/ntau,tau;
	std::complex<double> I(0.0,1.0);

	// Matsubara
	cdmatrix mat(size,size),tv(size,size);
	for(int m=0;m<=ntau;m++){
		tau=m*dtau;
		mat(0,0)=std::complex<double>(-1.7071067811865475,-0.17677669529663675)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-0.2928932188134524,+0.17677669529663687)*exp(2.0*tau);
		mat(0,1)=std::complex<double>(-0.4267766952966371,-1.0606601717798207)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-0.0732233047033631,+1.0606601717798212)*exp(2.0*tau);
		mat(1,0)=std::complex<double>(-0.07322330470336319,0.7071067811865475)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-0.4267766952966368,-0.7071067811865475)*exp(2.0*tau);
		mat(1,1)=std::complex<double>(-0.43933982822017864,0.17677669529663687)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-2.560660171779821,-0.17677669529663687)*exp(2.0*tau);
		mat=mat/(1.0+exp(2.0*beta));
		// std::cout << "Mat " << tau << " " << mat << std::endl;

		G.set_mat(m,mat);
		for(int n=0;n<=nt;n++){
			double t1=n*dt;
			tv(0,0)=std::complex<double>(0.17677669529663687,0.2928932188134524)*exp(2.0*I*t1 + 2.0*(beta-tau)) - std::complex<double>(0.1767766952966371,-1.707106781186548)*exp(-2.0*I*t1 + tau*2.0);
			tv(0,1)=std::complex<double>(1.0606601717798212,0.07322330470336319)*exp(2.0*I*t1 + 2.0*(beta-tau)) - std::complex<double>(1.0606601717798216,-0.426776695296637)*exp(-2.0*I*t1 + tau*2.0);
			tv(1,0)=-1.0*std::complex<double>(0.7071067811865475,-0.4267766952966369)*exp(2.0*I*t1 + 2.0*(beta-tau)) + std::complex<double>(0.7071067811865475,0.07322330470336313)*exp(-2.0*I*t1 + tau*2.0);
			tv(1,1)=-1.0*std::complex<double>(0.17677669529663675,-2.5606601717798214)*exp(2.0*I*t1 + 2.0*(beta-tau)) + std::complex<double>(0.1767766952966369,0.4393398282201787)*exp(-2.0*I*t1 + tau*2.0);
			tv=tv/(1.0+exp(2.0*beta));
			G.set_tv(n,m,tv);
		}
	}
	// Les + ret
	cdmatrix les(size,size),ret(size,size);
	for(int m=0;m<=nt;m++){
		for(int n=0;n<=m;n++){
			double t1=m*dt;
			double t2=n*dt;

			// Ret
			ret(0,0)=exp(-2.0*I*(t2+t1))*cos(t2)*(exp(4.0*I*t1)*std::complex<double>(-0.17677669529663687,-0.2928932188134524)+exp(4.0*I*t2)*std::complex<double>(0.1767766952966371,-1.707106781186548));
			ret(0,1)=exp(-2.0*I*(t2+t1))*cos(t2)*(exp(4.0*I*t1)*std::complex<double>(-1.0606601717798212,-0.07322330470336319)+exp(4.0*I*t2)*std::complex<double>(1.0606601717798216,-0.426776695296637));
			ret(1,0)=exp(-2.0*I*(t2+t1))*cos(t2)*(exp(4.0*I*t1)*std::complex<double>(0.7071067811865475,-0.4267766952966369)-exp(4.0*I*t2)*std::complex<double>(0.7071067811865475,0.07322330470336313));
			ret(1,1)=exp(-2.0*I*(t2+t1))*cos(t2)*(exp(4.0*I*t1)*std::complex<double>(0.17677669529663675,-2.5606601717798214)-exp(4.0*I*t2)*std::complex<double>(0.1767766952966369,0.4393398282201787));
			G.set_ret(m,n,ret);

			// Les
			t1=n*dt;
			t2=m*dt;
			les(0,0)=exp(-2.0*I*(t2+t1))*cos(t2)*(std::complex<double>(-0.1767766952966371,1.707106781186548)*exp(4.0*I*t2)+std::complex<double>(0.17677669529663687,0.2928932188134524)*exp(4.0*I*t1+2.0*beta));
			les(0,1)=exp(-2.0*I*(t2+t1))*cos(t2)*(std::complex<double>(-1.0606601717798216,0.426776695296637)*exp(4.0*I*t2)+std::complex<double>(1.0606601717798212,0.07322330470336319)*exp(4.0*I*t1+2.0*beta));
			les(1,0)=cos(t2)*(std::complex<double>(0.7071067811865475,0.07322330470336313)*exp(2.0*I*(t2-t1))+std::complex<double>(-0.7071067811865475,0.4267766952966369)*exp(2.0*I*(t1-t2)+2.0*beta));
			les(1,1)=cos(t2)*(std::complex<double>(0.1767766952966369,0.4393398282201787)*exp(2.0*I*(t2-t1))+std::complex<double>(-0.17677669529663675,2.5606601717798214)*exp(2.0*I*(t1-t2)+2.0*beta));
			les=les/(1.0+exp(2.0*beta));
			G.set_les(n,m,les);
		}
	}
}


void exact_leftmultiply(double beta,double dt,GREEN &G){
	int ntau=G.ntau();
	int nt=G.nt();
	int size=G.size1();
	assert(size==2);

	double dtau=beta/ntau,tau;
	std::complex<double> I(0.0,1.0);

	// Matsubara
	cdmatrix mat(size,size),tv(size,size);
	for(int m=0;m<=ntau;m++){
		tau=m*dtau;
		mat(0,0)=std::complex<double>(-1.7071067811865475,0.17677669529663675)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-0.2928932188134524,-0.17677669529663687)*exp(2.0*tau);
		mat(0,1)=std::complex<double>(-0.07322330470336319,-0.7071067811865475)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-0.4267766952966368,0.7071067811865475)*exp(2.0*tau);
		mat(1,0)=std::complex<double>(-0.4267766952966371,1.0606601717798207)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-0.0732233047033631,-1.0606601717798212)*exp(2.0*tau);
		mat(1,1)=std::complex<double>(-0.43933982822017864,-0.17677669529663687)*exp(-2*tau)*exp(beta*2.0)+std::complex<double>(-2.560660171779821,0.17677669529663687)*exp(2.0*tau);
		mat=mat/(1.0+exp(2.0*beta));
		// std::cout << "Mat " << tau << " " << mat << std::endl;

		G.set_mat(m,mat);
		for(int n=0;n<=nt;n++){
			double t1=n*dt;
			tv(0,0)=cos(t1)*(std::complex<double>(-0.17677669529663687,0.2928932188134524)*exp(2.0*I*t1 + 2.0*(beta-tau)) + std::complex<double>(0.1767766952966371,1.707106781186548)*exp(-2.0*I*t1 + tau*2.0));
			tv(0,1)=cos(t1)*(std::complex<double>(0.7071067811865475,0.4267766952966369)*exp(2.0*I*t1 + 2.0*(beta-tau)) + std::complex<double>(-0.7071067811865475,0.07322330470336313)*exp(-2.0*I*t1 + tau*2.0));
			tv(1,0)=cos(t1)*(std::complex<double>(-1.0606601717798212,0.07322330470336319)*exp(2.0*I*t1 + 2.0*(beta-tau)) + std::complex<double>(1.0606601717798216,0.426776695296637)*exp(-2.0*I*t1 + tau*2.0));
			tv(1,1)=cos(t1)*(std::complex<double>(0.17677669529663675,2.5606601717798214)*exp(2.0*I*t1 + 2.0*(beta-tau)) + std::complex<double>(-0.1767766952966369,0.4393398282201787)*exp(-2.0*I*t1 + tau*2.0));
			tv=tv/(1.0+exp(2.0*beta));
			G.set_tv(n,m,tv);
		}
	}
	// Les + ret
	cdmatrix les(size,size),ret(size,size);
	for(int m=0;m<=nt;m++){
		for(int n=0;n<=m;n++){
			double t1=m*dt;
			double t2=n*dt;

			// Ret
			ret(0,0)=exp(-2.0*I*(t2+t1))*cos(t1)*(exp(4.0*I*t1)*std::complex<double>(0.17677669529663687,-0.2928932188134524)-exp(4.0*I*t2)*std::complex<double>(0.1767766952966371,1.707106781186548));
			ret(0,1)=exp(-2.0*I*(t2+t1))*cos(t1)*(exp(4.0*I*t1)*std::complex<double>(-0.7071067811865475,-0.4267766952966369)+exp(4.0*I*t2)*std::complex<double>(0.7071067811865475,-0.07322330470336313));
			ret(1,0)=exp(-2.0*I*(t2+t1))*cos(t1)*(exp(4.0*I*t1)*std::complex<double>(1.0606601717798212,-0.07322330470336319)-exp(4.0*I*t2)*std::complex<double>(1.0606601717798216,0.426776695296637));
			ret(1,1)=exp(-2.0*I*(t2+t1))*cos(t1)*(exp(4.0*I*t1)*std::complex<double>(-0.17677669529663675,-2.5606601717798214)+exp(4.0*I*t2)*std::complex<double>(0.1767766952966369,-0.4393398282201787));
			G.set_ret(m,n,ret);

			// Les
			t1=n*dt;
			t2=m*dt;
			les(0,0)=cos(t1)*(std::complex<double>(0.1767766952966371,1.707106781186548)*exp(2.0*I*(t2-t1))-std::complex<double>(0.17677669529663687,-0.2928932188134524)*exp(2.0*I*(t1-t2)+2.0*beta));
			les(0,1)=cos(t1)*exp(-2.0*I*(t2+t1))*(std::complex<double>(-0.7071067811865475,0.07322330470336325)*exp(4.0*I*t2)+std::complex<double>(0.7071067811865476,0.4267766952966368)*exp(4.0*I*t1+2.0*beta));
			les(1,0)=cos(t1)*(std::complex<double>(1.0606601717798219,0.4267766952966371)*exp(2.0*I*(t2-t1))+std::complex<double>(-1.0606601717798212,0.07322330470336325)*exp(2.0*I*(t1-t2)+2.0*beta));
			les(1,1)=cos(t1)*exp(-2.0*I*(t2+t1))*(std::complex<double>(-0.1767766952966369,0.4393398282201787)*exp(4.0*I*t2)+std::complex<double>(0.17677669529663675,2.5606601717798214)*exp(4.0*I*t1+2.0*beta));
			les=les/(1.0+exp(2.0*beta));
			G.set_les(n,m,les);
		}
	}
}



TEST_CASE("Herm member","[Herm_member]"){
	int nt=50;
	int ntau=500;
	int size=2;
	int size2=5;
	double eps=1e-7;
	double beta=5.0;
	double h=0.01;

	GREEN g=GREEN(nt,ntau,size,-1);
	REQUIRE(g.nt()==nt);

	// Square Green's function
	GREEN A=GREEN(nt,ntau,size,-1);
	cdmatrix a(size,size);
	a.setZero();
	a(0,0)=sqrt(2.0);
	a(0,1)=sqrt(2.0)*std::complex<double>(0.0,1.0);
	a(1,0)=sqrt(2.0)*std::complex<double>(0.0,-1.0);
	a(1,1)=-sqrt(2.0);
	cntr::green_from_H(A,0.0,a,beta,h,5,4,true);

	// Rectangular Green's function
	GREEN B=GREEN(nt,ntau,size,size2,-1);
	cdmatrix b(size,size2);
	b.setZero();
	
	for (int i=0;i<size;i++){
		for (int j=0;j<size2;j++){
			b(i,j)=i+j*2;
		}
	}

	// Function
	cfunction funcC(nt,size);
	cdmatrix c(size,size);
	c.setZero();
	double t;
	for (int tstp=-1;tstp<nt;tstp++){
				
		if(tstp==-1){
			t=0.0;
		}else{
			t=tstp*h;
		}
		c(0,0)=2.0*cos(t);
		c(0,1)=0.5*cos(t);
		c(1,0)=0.5*cos(t);
		c(1,1)=3.0*cos(t);
		funcC.set_value(tstp,c);
	}

	SECTION ("Set/get for herm_matrix"){
		// Writing/reading from normal green's function
		REQUIRE(setget(A,a)<eps);
		// Writing/reading from rectangular Green's function
		REQUIRE(setget(B,b)<eps);
	}

	SECTION ("Left/right multiply"){
		GREEN Ar,Al;
		//Left/Right multiply
		Ar=A;
		Al=A;
		for(int tstp=-1;tstp<=nt;tstp++){
			Ar.right_multiply(tstp,funcC);
			Al.left_multiply(tstp,funcC);
		}

		double err=0.0;
		GREEN exactR(nt,ntau,size,-1);
		exact_rightmultiply(beta,h,exactR);
		for(int tstp=-1;tstp<nt;tstp++){
			err+=cntr::distance_norm2(tstp,Ar,exactR);
		}
		REQUIRE(err<eps);

		err=0.0;
		GREEN exactL(nt,ntau,size,-1);
		exact_leftmultiply(beta,h,exactL);
		for(int tstp=-1;tstp<nt;tstp++){
			err+=cntr::distance_norm2(tstp,Al,exactL);
		}
		REQUIRE(err<eps);
		
		// REQUIRE(simple_readwrite(A,a)<eps);
		// // Writing/reading from rectangular Green's function
		// REQUIRE(simple_readwrite(B,b)<eps);
	}

}
