#include "catch.hpp"
#include "cntr.hpp"


#define CPLX std::complex<double>  
using namespace std;
#define GREEN cntr::herm_matrix<double> 
#define GREEN_TSTP cntr::herm_matrix_timestep<double> 

void timedependent_from_vie(GREEN &G,double mu,cntr::function<double> &epst,double beta,double h){
  int nt=G.nt();
  int ntau=G.ntau();
  int size=G.size1();

  GREEN G0(nt,ntau,size,-1);
  GREEN G0eps(nt,ntau,size,-1);
  GREEN epsG0(nt,ntau,size,-1);
  cdmatrix eps0(size,size);
  eps0.setZero();

  cntr::green_from_H(G0,mu,eps0,beta,h);
  cdmatrix tmp;
  G0.get_mat(0,tmp);
  G0eps=G0;
  epsG0=G0;

  for(int tstp=-1;tstp<=nt;tstp++){
    G0eps.right_multiply(tstp,epst,-1.0);
    epsG0.left_multiply(tstp,epst,-1.0);	
  }

  cntr::vie2(G,G0eps,epsG0,G0,integration::I<double>(5),beta,h);
}


void exact_timeindependent(double beta,double dt,GREEN &G){
  int ntau=G.ntau();
  int nt=G.nt();
  int size=G.size1();

  double dtau=beta/ntau,tau;
  std::complex<double> I(0.0,1.0);

  // Matsubara
  cdmatrix mat(size,size),tv(size,size);
  for(int m=0;m<=ntau;m++){
    tau=m*dtau;
    mat(0,0)=exp(-2.0*tau)*(-0.8535533905932737*exp(beta*2.0)-0.1464466094067262*exp(4.0*tau));
    mat(1,1)=exp(-2.0*tau)*(-0.1464466094067262*exp(beta*2.0)-0.8535533905932737*exp(4.0*tau));
    mat(1,0)=-1.0*I*0.35355339059327373*(exp(2.0*(beta-tau)) - exp(2.0*tau));
    mat(0,1)=1.0*I*0.35355339059327373*(exp(2.0*(beta-tau)) - exp(2.0*tau));
    mat=mat/(1.0+exp(2.0*beta));
    // std::cout << "Mat " << tau << " " << mat << std::endl;

    G.set_mat(m,mat);
    for(int n=0;n<=nt;n++){
      double t1=n*dt;
      tv(0,0)=I*0.8535533905932737*exp(-2.0*I*t1 + 2.0*tau) + I*0.1464466094067262*exp(2.0*I*t1 - tau*2.0 + 2.0*beta);
      tv(0,1)=-0.35355339059327373*(exp(2.0*I*t1+2.0*beta-2.0*tau)-exp(-2.0*I*t1+2.0*tau));
      tv(1,0)=0.35355339059327373*(exp(2.0*I*t1+2.0*beta-2.0*tau)-exp(-2.0*I*t1+2.0*tau));
      tv(1,1)=I*0.8535533905932737*(exp(2.0*I*t1+2.0*beta-2.0*tau))+0.1464466094067262*I*exp(-2.0*I*t1+2.0*tau);
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
      ret(0,0)=-1.0*I*cos(2.0*(t1 -t2)) - 0.7071067811865475*sin(2.0*(t1 - t2));
      ret(1,0)=-0.7071067811865475*I*sin(2.0*(t1 - t2));
      ret(0,1)=0.7071067811865475*I*sin(2.0*(t1 - t2));
      ret(1,1)=-1.0*I*cos(2.0*(t1 - t2)) + 0.7071067811865475*sin(2.0*(t1 - t2));
      G.set_ret(m,n,ret);

      // Les
      t1=n*dt;
      t2=m*dt;
      les(0,0)=0.8535533905932737*I*std::exp(2.0*I*(t2-t1)) + 0.1464466094067262*I*std::exp(2.0*I*(t1-t2) + 2.0*beta);
      les(1,0)=-0.35355339059327373*std::exp(2.0*I*(t2-t1)) + 0.35355339059327373*std::exp(2.0*I*(t1-t2)+2.0*beta);
      les(0,1)=-0.35355339059327373*std::exp(2.0*I*(t1-t2)+2.0*beta) + 0.35355339059327373*std::exp(-2.0*I*(t1-t2));
      les(1,1)=0.8535533905932737*I*std::exp(2.0*I*(t1-t2)+ 2.0*beta) + 0.1464466094067262*I*std::exp(2.0*I*(t2-t1));
      les=les/(1.0+exp(2.0*beta));
      G.set_les(n,m,les);
    }
  }
}

// -------------------------------
// Driven harmonic oscillator
// H(t)= 
// [Delta,V exp(-2 I om t)]
// [V exp(2 I om t), -Delta]
// and the exact propagator  is given by
// see Eq. 50 in A.Alvermann in Journal of Computational Physics 230, 5930-5956, doi:10.1016/j.jcp.2011.04.006
// U(t,0)=
// [exp(-I omega t)(cos(Om*t)-I*(Delta-omega)/Om*Sin(Om*t)) ; -I V/Om * exp(-I*omega*t)*sin(Om*t)]
// [-I V/Om * exp(I*omega*t)*sin(Om*t)] ;  exp(I omega t)(cos(Om*t)+I*(Delta-omega)/Om*Sin(Om*t))
// -------------------------------

cdmatrix U(double t,int size,double Delta,double V,double omega){
  cdmatrix U(size,size);
  std::complex<double> I(0.0,1.0);
  double Om=sqrt((Delta-omega)*(Delta-omega)+V*V);

  U(0,0)=exp(-I*omega*t)*(cos(Om*t)-I*(Delta-omega)*sin(Om*t)/Om);
  U(0,1)=-I*V*exp(-I*omega*t)*sin(Om*t)/Om;
  U(1,0)=-I*V*exp(I*omega*t)*sin(Om*t)/Om;
  U(1,1)=exp(I*omega*t)*(cos(Om*t)+I*(Delta-omega)*sin(Om*t)/Om);

  return U;
}

void exact_time(double beta,double dt,GREEN &G,double Delta,double V,double omega){
  int ntau=G.ntau();
  int nt=G.nt();
  int size=G.size1();

  double dtau=beta/ntau,tau;
  std::complex<double> I(0.0,1.0);

  // Eigenvectors
  cdmatrix vec(2,2),fermiDiag(2,2);
  vec(0,0)=0.22975292054736116;
  vec(0,1)=-0.9732489894677301;
  vec(1,0)=-0.9732489894677301;
  vec(1,1)=-0.22975292054736116;

  fermiDiag(0,0)=1.0/(1.0+exp(-1.118033988749895*beta));
  fermiDiag(0,1)=0.0;
  fermiDiag(1,0)=0.0;
  fermiDiag(1,1)=1.0/(1.0+exp(1.118033988749895*beta));

  // Matsubara
  cdmatrix mat(size,size),tv(size,size);
  for(int m=0;m<=ntau;m++){
    tau=m*dtau;
    mat(0,0)=exp(-1.118033988749895*tau)*(-0.9472135954999577*exp(1.118033988749895*beta)- 0.05278640450004204*exp(2.23606797749979*tau));
    mat(0,1)=exp(-1.118033988749895*tau)*(-0.22360679774997894*exp(1.118033988749895*beta)+ 0.22360679774997894*exp(2.23606797749979*tau));
    mat(1,0)=exp(-1.118033988749895*tau)*(-0.22360679774997894*exp(1.118033988749895*beta)+ 0.22360679774997894*exp(2.23606797749979*tau));
    mat(1,1)=exp(-1.118033988749895*tau)*(-0.052786404500042045*exp(1.118033988749895*beta)- 0.9472135954999577*exp(2.23606797749979*tau));
    mat=mat/(1.0+exp(1.118033988749895*beta));
    // std::cout << "Mat " << tau << " " << mat << std::endl;

    G.set_mat(m,mat);
    for(int n=0;n<=nt;n++){
      double t1=n*dt;
      cdmatrix expon(2,2);
      expon(0,0)=exp(-1.118033988749895*tau);expon(0,1)=0.0;expon(1,0)=0.0;expon(1,1)=exp(1.118033988749895*tau);
      tv=U(t1,size,Delta,V,omega)*vec.transpose()*fermiDiag*expon*vec.conjugate();
      tv=tv*I;
      G.set_tv(n,m,tv);
    }
  }
  // Les + ret
  cdmatrix les(size,size),ret(size,size),Ut1(size,size),Ut2(size,size);
  for(int m=0;m<=nt;m++){
    for(int n=0;n<=m;n++){
      double t1=m*dt;
      double t2=n*dt;
      Ut1=U(t1,size,Delta,V,omega);
      Ut2=U(t2,size,Delta,V,omega);

      // Ret
      ret=-I*Ut1*Ut2.adjoint();
      G.set_ret(m,n,ret);

      // Les (mark that arguments are switched)
      les=I*Ut2*vec.transpose()*fermiDiag*vec.conjugate()*Ut1.adjoint();
      G.set_les(n,m,les);
    }
  }
}


void exact_time_boson(double beta,double dt,GREEN &G,double Delta,double V,double omega){
  int ntau=G.ntau();
  int nt=G.nt();
  int size=G.size1();

  double dtau=beta/ntau,tau;
  std::complex<double> I(0.0,1.0);

  // Eigenvectors
  cdmatrix vec(2,2),boseDiag(2,2);
  vec(0,0)=0.22975292054736116;
  vec(0,1)=-0.9732489894677301;
  vec(1,0)=-0.9732489894677301;
  vec(1,1)=-0.22975292054736116;

  boseDiag(0,0)=1.0/(-1.0+exp(-1.118033988749895*beta));
  boseDiag(0,1)=0.0;
  boseDiag(1,0)=0.0;
  boseDiag(1,1)=1.0/(-1.0+exp(1.118033988749895*beta));

  // Matsubara
  cdmatrix mat(size,size),tv(size,size);
  for(int m=0;m<=ntau;m++){
    tau=m*dtau;
    mat(0,0)=exp(-1.118033988749895*tau)*(-0.9472135954999577*exp(1.118033988749895*beta)+0.05278640450004204*exp(2.23606797749979*tau));
    mat(0,1)=exp(-1.118033988749895*tau)*(-0.22360679774997894*exp(1.118033988749895*beta)- 0.22360679774997894*exp(2.23606797749979*tau));
    mat(1,0)=exp(-1.118033988749895*tau)*(-0.22360679774997894*exp(1.118033988749895*beta)- 0.22360679774997894*exp(2.23606797749979*tau));
    mat(1,1)=exp(-1.118033988749895*tau)*(-0.052786404500042045*exp(1.118033988749895*beta)+ 0.9472135954999577*exp(2.23606797749979*tau));
    mat=mat/(-1.0+exp(1.118033988749895*beta));

    G.set_mat(m,mat);
    for(int n=0;n<=nt;n++){
      double t1=n*dt;
      cdmatrix expon(2,2);
      expon(0,0)=exp(-1.118033988749895*tau);expon(0,1)=0.0;expon(1,0)=0.0;expon(1,1)=exp(1.118033988749895*tau);
      tv=U(t1,size,Delta,V,omega)*vec.transpose()*boseDiag*expon*vec.conjugate();
      tv=tv*I*(-1.0);
      G.set_tv(n,m,tv);
    }
  }
  // Les + ret
  cdmatrix les(size,size),ret(size,size),Ut1(size,size),Ut2(size,size);
  for(int m=0;m<=nt;m++){
    for(int n=0;n<=m;n++){
      double t1=m*dt;
      double t2=n*dt;
      Ut1=U(t1,size,Delta,V,omega);
      Ut2=U(t2,size,Delta,V,omega);

      // Ret
      ret=-I*Ut1*Ut2.adjoint();
      G.set_ret(m,n,ret);

      // Les (mark that arguments are switched)
      les=(-1.0)*I*Ut2*vec.transpose()*boseDiag*vec.conjugate()*Ut1.adjoint();
      G.set_les(n,m,les);
    }
  }
}

TEST_CASE("Non-interacting","[Herm_matrix_noniteracting]"){
  int nt,ntau,kt=5,tstp,size=2;
  double wa[size];
  std::complex<double> waa[size*size];
  double err,err2,err4,err4B,beta,h,eps=1e-6;
  GREEN A,ACF2,ACF4,B,C;
  cntr::function<double> epsilon;
  double mu=0.0;
  nt=100;
  ntau=500;
  beta=5.0;
  h=0.01;

  // Test diagonal
  epsilon=cntr::function<double>(nt,size,size);

  SECTION ("Time independent"){
    GREEN Atest(nt,ntau,size,-1),A2test(nt,ntau,size,-1),A3test(nt,ntau,size,-1),A4test(nt,ntau,size,-1),A5test(nt,ntau,size,-1);
    GREEN exact(nt,ntau,size,-1);

    cdmatrix a(size,size);
    a.setZero();
    a(0,0)=sqrt(2.0);
    a(1,0)=sqrt(2.0)*std::complex<double>(0.0,1.0);
    a(0,1)=sqrt(2.0)*std::complex<double>(0.0,-1.0);
    a(1,1)=-sqrt(2.0);

    epsilon.set_constant(a);
    exact_timeindependent(beta,h,exact);

    // cntr::green_from_eps(Atest,mu,wa,beta,h);
    // cntr::green_from_H(A2test,mu,waa,beta,h);
    cntr::green_from_H(Atest,mu,a,beta,h);
    cntr::green_from_H(A2test,mu,epsilon,beta,h,true);
    // timedependent_from_vie(A4test,mu,epsilon,beta,h);

    double err=0.0;
    for(int tstp=-1;tstp<nt;tstp++){
      err+=cntr::distance_norm2(tstp,Atest,exact);
    }
    REQUIRE(err<eps);


    // std::cout << "Eps " <<err <<std::endl;
    err=0.0;
    for(int tstp=-1;tstp<nt;tstp++){
      err+=cntr::distance_norm2(tstp,A2test,exact);
    }
    REQUIRE(err<eps);
    // std::cout << "Eps " <<err <<std::endl;
    // err=0.0;
    // for(int tstp=-1;tstp<nt;tstp++){
    // 	err+=cntr::distance_norm2(tstp,A4test,exact);
    // }
    // // std::cout << "Eps " <<err <<std::endl;
    // REQUIRE(err<eps);
  }

  SECTION ("Time-dependent-fermions"){
    A=GREEN(nt,ntau,size,-1);
    ACF2=GREEN(nt,ntau,size,-1);
    ACF4=GREEN(nt,ntau,size,-1);
    B=GREEN(nt,ntau,size,-1);
    GREEN exact(nt,ntau,size,-1);
    // These parameters should be fixed, since exact comparison point assumes it
    double delta=1.0,V=0.5,omega=0.7;
    std::complex<double> I(0.0,1.0);

		

    double t;
    cntr::function<double> epsT(nt,size),epsTrun(nt,size),UCN(nt,size),Uexp2(nt,size),Uexp4(nt,size);
    for(int tstp=-1;tstp<=nt;tstp++){
      cdmatrix tmp(size,size);
      t=tstp*h;
      if(tstp==-1){
	tmp(0,0)=delta;
	tmp(0,1)=V;
	tmp(1,0)=V;
	tmp(1,1)=-delta;
      }else{
	tmp(0,0)=delta;
	tmp(0,1)=V*exp(-2.0*I*omega*t);
	tmp(1,0)=V*exp(2.0*I*omega*t);
	tmp(1,1)=-delta;
      }
      epsT.set_value(tstp,tmp);
    }

    // -------------------------------------
    // Test propagators - fix Hamiltonian
    // -------------------------------------
    err=0.0;
    // cntr::Schroedinger_propagator_CrayN<double,2>(size,UCN.ptr(0),epsT.ptr(0),nt,h,kt,1);
    for(int tstp=-1;tstp<nt;tstp++){
      cntr::propagator_exp<double>(tstp,Uexp2,epsT,h,2,kt,true);
      cntr::propagator_exp<double>(tstp,Uexp4,epsT,h,4,kt,true);
      cdmatrix tmp,tmp2,tmp4;
      Uexp2.get_value(tstp,tmp2);
      Uexp4.get_value(tstp,tmp4);
      if(tstp==-1){
	err+=(U(0,size,delta,V,omega)-tmp4).norm();	
      }else{
	err+=(U(tstp*h,size,delta,V,omega)-tmp4).norm();
      }
    }
    REQUIRE(err<eps);

    // -------------------------------------
    // Test propagators - running Hamiltonian
    // Idea is to have hamiltonian until tstp-1, extrapolate and obtain the propagator for the new time-step
    // This is useful for predictor corrector like methods. 
    // TODO: write test which does predictor-corrector for Hartree-Fock or something like this
    // -------------------------------------
    err=0.0;
    for(int tstp=-1;tstp<nt;tstp++){
      int n1=(tstp<=kt ? kt : tstp-1);

      for(int trun=-1;trun<=n1;trun++){
	cdmatrix tmp;
	epsT.get_value(trun,tmp);
	epsTrun.set_value(trun,tmp);
      }
      cntr::propagator_exp<double>(tstp,Uexp4,epsTrun,h,4,kt,false);
      cdmatrix tmp;
      Uexp4.get_value(tstp,tmp);
      if(tstp==-1){
	err+=(U(0,size,delta,V,omega)-tmp).norm();	
      }else{
	err+=(U(tstp*h,size,delta,V,omega)-tmp).norm();
      }
    }
    REQUIRE(err<eps);


    // -------------------------------------
    // Test Green - fix Hamiltonian
    // -------------------------------------
    green_from_H(ACF4,mu,epsT,beta,h,true);
    exact_time(beta,h,exact,delta,V,omega);
    // Interesting observation is that for these parameters vie2 solution is more accurate, but much slower
    // Please check if this is a general case
    // timedependent_from_vie(B,mu,epsT,beta,h);

    err=0.0;
    for(int tstp=-1;tstp<=nt;tstp++){
      err+=cntr::distance_norm2(tstp,ACF4,exact);
    }
    REQUIRE(err<eps);



    // -------------------------------------
    // Test Green_tstp - fix Hamiltonian
    // -------------------------------------
    err=0.0;
    for(int tstp=-1;tstp<=nt;tstp++){
      GREEN_TSTP CF4tstp(tstp,ntau,size);
      green_from_H(CF4tstp,mu,epsT,beta,h,true);
      err+=cntr::distance_norm2(tstp,CF4tstp,exact);
    }
    REQUIRE(err<eps);
  }

  SECTION ("Time-dependent-bosons"){
    A=GREEN(nt,ntau,size,1);
    ACF2=GREEN(nt,ntau,size,1);
    ACF4=GREEN(nt,ntau,size,1);
    B=GREEN(nt,ntau,size,1);
    GREEN exact(nt,ntau,size,1);
    // These parameters should be fixed, since exact comparison point assumes it
    double delta=1.0,V=0.5,omega=0.7;
    std::complex<double> I(0.0,1.0);

    double t;
    cntr::function<double> epsT(nt,size),epsTrun(nt,size),UCN(nt,size),Uexp2(nt,size),Uexp4(nt,size);
    for(int tstp=-1;tstp<=nt;tstp++){
      cdmatrix tmp(size,size);
      t=tstp*h;
      if(tstp==-1){
        tmp(0,0)=delta;
        tmp(0,1)=V;
        tmp(1,0)=V;
        tmp(1,1)=-delta;
      }else{
        tmp(0,0)=delta;
        tmp(0,1)=V*exp(-2.0*I*omega*t);
        tmp(1,0)=V*exp(2.0*I*omega*t);
        tmp(1,1)=-delta;
      }
      epsT.set_value(tstp,tmp);
    }

    // -------------------------------------
    // Test Green - fix Hamiltonian
    // -------------------------------------
    green_from_H(ACF4,mu,epsT,beta,h,true);
    exact_time_boson(beta,h,exact,delta,V,omega);
    // Interesting observation is that for these parameters vie2 solution is more accurate, but much slower
    // Please check if this is a general case
    // timedependent_from_vie(B,mu,epsT,beta,h);

    err=0.0;
    for(int tstp=-1;tstp<=nt;tstp++){
      err+=cntr::distance_norm2(tstp,ACF4,exact);
    }
    REQUIRE(err<eps);
  }

}




