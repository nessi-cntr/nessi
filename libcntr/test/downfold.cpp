#include "catch.hpp"
#include "cntr.hpp"


using namespace std;
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double> 
#define CFUNCTION cntr::function<double>

/*///////////////////////////////////////////////////////////////////////////////////////

Test of the downfolding on the three level system

///////////////////////////////////////////////////////////////////////////////////////*/ 

// Embeded self energy for subspace S elements (i,j)
// from the total space T, where 
// the rest of the space R is integrated out 
// Assumption that the off-diagonal terms are not retarded
void downfold(int tstp,int i,int j,GREEN &R,GREEN &S,CFUNCTION &epsilon){
	int nt=R.nt();
	int ntau=R.ntau();
	int ssize=S.size1();
	int Rsize=R.size1();

	GREEN_TSTP etatmp(tstp,ntau,ssize);
	cntr::function<double> epsR(nt,ssize),epsL(nt,ssize);
	etatmp.clear();
	for(int k=0;k<Rsize;k++){
		for(int l=0;l<Rsize;l++){
			GREEN_TSTP tmp(tstp,ntau,1);
			tmp.set_matrixelement(0,0,R,k,l);
			epsR.set_matrixelement(0,0,epsilon,ssize+l,j);
			epsL.set_matrixelement(0,0,epsilon,i,ssize+k);
			tmp.right_multiply(epsR,1.0);
			tmp.left_multiply(epsL,1.0);
			etatmp.incr(tmp,1.0);
		}
	}
	S.set_matrixelement(tstp,i,j,etatmp);
}

// Embeded self energy for subspace S elements (i,j)
// from the total space T, where  the rest of the space R is integrated out 
// and the self energy of total system is Sigma

void downfold3(int tstp,int i,int j,GREEN &R,GREEN &S,CFUNCTION &epsilon,GREEN &Sigma,GREEN &GSigma,GREEN &SigmaG,int kt,double beta,double h){
	int nt=R.nt();
	int ntau=R.ntau();
	int ssize=S.size1();
	int Rsize=R.size1();
	
	int n1=(tstp<=kt && tstp>=0 ? 0 : tstp);
	int n2=(tstp<=kt && tstp>=0 ? kt : tstp);
	GREEN_TSTP etatmp(tstp,ntau,ssize);
	cntr::function<double> epsR(nt,ssize),epsL(nt,ssize);
	etatmp.clear();
	for(int k=0;k<Rsize;k++){
		for(int l=0;l<Rsize;l++){
			GREEN_TSTP tmp(tstp,ntau,1);
			tmp.set_matrixelement(0,0,R,k,l);
			epsR.set_matrixelement(0,0,epsilon,ssize+l,j);
			epsL.set_matrixelement(0,0,epsilon,i,ssize+k);
			tmp.right_multiply(epsR,1.0);
			tmp.left_multiply(epsL,1.0);
			etatmp.incr(tmp,1.0);
		}
	}
	 
	// Auxiliary quantity
	// [G * Sigma ]_{k0}=\sum_l G_{kl} * \sigma_{l0}
	for(int n=n1;n<=n2;n++){
	  for(int k=0;k<Rsize;k++){
	    GREEN_TSTP sum(n,ntau,1);
	    sum.clear();
	    for(int l=0;l<Rsize;l++){
	      GREEN gr(nt,ntau,1);
	      gr.set_matrixelement(0,0,R,k,l);
	      GREEN Sigmak_nl_rest(nt,ntau,1);
	      Sigmak_nl_rest.set_matrixelement(0,0,Sigma,ssize+l,j);
	      GREEN Sigmak_nl_rest_cc(nt,ntau,1);
	      Sigmak_nl_rest_cc.set_matrixelement(0,0,Sigma,j,ssize+l);
	      GREEN tmp(nt,ntau,1);
	      cntr::convolution_timestep(n,tmp,gr,gr,Sigmak_nl_rest,Sigmak_nl_rest_cc,integration::I<double>(kt),beta,h);
	      sum.incr(tmp,1.0);
	    }
	    GSigma.set_matrixelement(n,k,0,sum);
	  }
	}
	
	// cc conjugate of [G * Sigma ]_{k0} = [Sigma * G ]_{0k}=\sum_l \Sigma_{0l} * G_{l k}
	for(int n=n1;n<=n2;n++){
	  for(int k=0;k<Rsize;k++){
	    GREEN_TSTP sum(n,ntau,1);
	    sum.clear();
	    for(int l=0;l<Rsize;l++){
	      GREEN gr(nt,ntau,1);
	      gr.set_matrixelement(0,0,R,l,k);
	      GREEN Sigmak_nl_rest(nt,ntau,1);
	      Sigmak_nl_rest.set_matrixelement(0,0,Sigma,i,ssize+l);
	      GREEN Sigmak_nl_rest_cc(nt,ntau,1);
	      Sigmak_nl_rest_cc.set_matrixelement(0,0,Sigma,ssize+l,i);
	      GREEN tmp(nt,ntau,1);
	      cntr::convolution_timestep(n,tmp,Sigmak_nl_rest,Sigmak_nl_rest_cc,gr,gr,integration::I<double>(kt),beta,h);
	      sum.incr(tmp,1.0);
	    }
	    SigmaG.set_matrixelement(n,0,k,sum);
	  }
	}
	// Mixed terms 
	// [eta * G * Sigma]_{00}= [eta]_ik [[G]_kl [Sigma]_lj]
	for(int k=0;k<Rsize;k++){
	  GREEN_TSTP tmp(tstp,ntau,1);
	  tmp.set_matrixelement(0,0,GSigma,k,j);
	  epsL.set_matrixelement(0,0,epsilon,i,ssize+k);
	  tmp.left_multiply(epsL,1.0);
	  etatmp.incr(tmp,1.0);
	}
	// Mixed terms 
	// [sigma * G * eta]_{00}= [Sigma]_ik [G]_kl [eta]_lj
	for(int l=0;l<Rsize;l++){
	  GREEN_TSTP tmp(tstp,ntau,1);
	  tmp.set_matrixelement(0,0,SigmaG,i,l);
	  epsR.set_matrixelement(0,0,epsilon,ssize+l,j);
	  tmp.left_multiply(epsR,1.0);
	  etatmp.incr(tmp,1.0);
	}

	// [Sigma * G * Sigma ]_{00}=\sum_k  \Sigma_0k * G_kl * \Sigma_l0
	// In practice : Sigma * [G * Sigma], where [G*Sigma] is stored auxiliary quantity
	// [Sigma * G * Sigma ]_{00}=\sum_k  \Sigma_{0k} * [G * \Sigma]_{k 0}
	
	GREEN_TSTP SkGkSk(tstp,ntau,1);
	SkGkSk.clear();
	for(int k=0;k<Rsize;k++){
	  GREEN GStmp(nt,ntau,1);
	  GStmp.set_matrixelement(0,0,GSigma,k,j);

	  GREEN SGtmp(nt,ntau,1);
	  SGtmp.set_matrixelement(0,0,SigmaG,j,k);
	    
	  GREEN Sigmak_nl_rest(nt,ntau,1);
	  Sigmak_nl_rest.set_matrixelement(0,0,Sigma,i,ssize+k);

	  GREEN Sigmak_nl_rest_cc(nt,ntau,1);
	  Sigmak_nl_rest_cc.set_matrixelement(0,0,Sigma,ssize+k,i);
	    
	  GREEN tmp(nt,ntau,1);
	  cntr::convolution_timestep(tstp,tmp,Sigmak_nl_rest,Sigmak_nl_rest_cc,GStmp,SGtmp,integration::I<double>(kt),beta,h);
	  etatmp.incr(tmp,1.0);
	    //SkGkSk.incr(tmp,1.0);
	}
	//etatmp.incr(SkGkSk,1.0);
	//Add Sigmak_nl_[q]_{dd}
	GREEN_TSTP Sigmak_nl_tmp(tstp,ntau,1);
	Sigmak_nl_tmp.set_matrixelement(0,0,Sigma,i,j);
	etatmp.incr(Sigmak_nl_tmp,1.0);
	S.set_matrixelement(tstp,i,j,etatmp);
	
}



// Older routine with a bit of doubling of convolutions
// Embeded self energy for subspace S elements (i,j)
// from the total space T, where  the rest of the space R is integrated out 
// and the self energy of total system is Sigma
// void downfold2(int tstp,int i,int j,GREEN &R,GREEN &S,CFUNCTION &epsilon,GREEN &Sigma,GREEN &GSigma,GREEN &SigmaG,int kt,double beta,double h){
// 	int nt=R.nt();
// 	int ntau=R.ntau();
// 	int ssize=S.size1();
// 	int Rsize=R.size1();
	
// 	int n1=(tstp<=kt && tstp>=0 ? 0 : tstp);
// 	int n2=(tstp<=kt && tstp>=0 ? kt : tstp);
// 	GREEN_TSTP etatmp(tstp,ntau,ssize);
// 	cntr::function<double> epsR(nt,ssize),epsL(nt,ssize);
// 	etatmp.clear();
// 	for(int k=0;k<Rsize;k++){
// 		for(int l=0;l<Rsize;l++){
// 			GREEN_TSTP tmp(tstp,ntau,1);
// 			tmp.set_matrixelement(0,0,R,k,l);
// 			epsR.set_matrixelement(0,0,epsilon,ssize+l,j);
// 			epsL.set_matrixelement(0,0,epsilon,i,ssize+k);
// 			tmp.right_multiply(epsR,1.0);
// 			tmp.left_multiply(epsL,1.0);
// 			etatmp.incr(tmp,1.0);
// 		}
// 	}
// 	//S.set_matrixelement(tstp,i,j,etatmp);
// 	// Mixed terms 
// 	// [eta * G * Sigma]_{00}= [eta]_ik [G]_kl [Sigma]_lj
// 	// [eta * G]_il= [eta]_ik [G]_kl and cc [G * eta]_li= G_lk [eta]_ki
// 	GREEN_TSTP etaGS(tstp,ntau,1);
// 	etaGS.clear();
// 	for(int l=0;l<Rsize;l++){
// 	  GREEN sumEtaG(nt,ntau,1),sumGEta(nt,ntau,1);
// 	  sumEtaG.clear();sumGEta.clear();
// 	  for(int k=0;k<Rsize;k++){
// 	    GREEN etaG(nt,ntau,1),Geta(nt,ntau,1);
// 	    etaG.set_matrixelement(0,0,R,k,l);
// 	    Geta.set_matrixelement(0,0,R,l,k);
	      
// 	    epsR.set_matrixelement(0,0,epsilon,ssize+k,i);
// 	    epsL.set_matrixelement(0,0,epsilon,i,ssize+k);
// 	    for(int n=-1;n<=n2;n++){
// 	      etaG.left_multiply(n,epsL);
// 	      Geta.right_multiply(n,epsR);
// 	    }
// 	    //Incr for all timesteps
	     
// 	    sumEtaG.incr_timestep(etaG,1.0);
// 	    sumGEta.incr_timestep(Geta,1.0);
	      
// 	  }
// 	  GREEN Sigmak_nl_rest(nt,ntau,1);
// 	  Sigmak_nl_rest.set_matrixelement(0,0,Sigma,ssize+l,j);
// 	  GREEN Sigmak_nl_rest_cc(nt,ntau,1);
// 	  Sigmak_nl_rest_cc.set_matrixelement(0,0,Sigma,j,ssize+l);
	  
// 	  GREEN tmp(nt,ntau,1);
// 	  cntr::convolution_timestep(tstp,tmp,sumEtaG,sumGEta,Sigmak_nl_rest,Sigmak_nl_rest_cc,integration::I<double>(kt),beta,h);
	  
// 	  etaGS.incr(tmp,1.0);
// 	}
// 	etatmp.incr(etaGS,1.0);
// 	// Mixed terms 
// 	// [sigma * G * eta]_{00}= [Sigma]_ik [G]_kl [eta]_lj
// 	// [G * eta]_kj= [G]_kl [eta]_lj and cc [eta * G]_jk= [eta]_jl [G]_lk 
// 	GREEN_TSTP SGeta(tstp,ntau,1);
// 	SGeta.clear();
// 	for(int k=0;k<Rsize;k++){
// 	  GREEN sumEtaG(nt,ntau,1),sumGEta(nt,ntau,1);
// 	  sumEtaG.clear();sumGEta.clear();
// 	  for(int l=0;l<Rsize;l++){
// 	    GREEN etaG(nt,ntau,1),Geta(nt,ntau,1);
// 	    etaG.set_matrixelement(0,0,R,l,k);
// 	    Geta.set_matrixelement(0,0,R,k,l);
// 	    epsR.set_matrixelement(0,0,epsilon,ssize+l,j);
// 	    epsL.set_matrixelement(0,0,epsilon,j,ssize+l);
// 	    for(int n=-1;n<=n2;n++){
// 	      etaG.left_multiply(n,epsL);
// 	      Geta.right_multiply(n,epsR);
// 	    }
// 	    //Incr for all timesteps
// 	    sumEtaG.incr_timestep(etaG,1.0);
// 	    sumGEta.incr_timestep(Geta,1.0);
// 	  }
// 	  GREEN Sigmak_nl_rest(nt,ntau,1);
// 	  Sigmak_nl_rest.set_matrixelement(0,0,Sigma,i,ssize+k);
// 	  GREEN Sigmak_nl_rest_cc(nt,ntau,1);
// 	  Sigmak_nl_rest_cc.set_matrixelement(0,0,Sigma,ssize+k,i);
// 	  GREEN tmp(nt,ntau,1);
// 	  cntr::convolution_timestep(tstp,tmp,sumEtaG,sumGEta,Sigmak_nl_rest,Sigmak_nl_rest_cc,integration::I<double>(kt),beta,h);
// 	  SGeta.incr(tmp,1.0);
// 	}
// 	etatmp.incr(SGeta,1.0);
// 	// [Sigma * G * Sigma ]_{00}=\sum_k  \Sigma_0k * G_kl * \Sigma_l0
// 	// In practice : Sigma * [G * Sigma], where [G*Sigma] is stored auxiliary quantity
	  
// 	// Auxiliary quantity
// 	// [G * Sigma ]_{k0}=\sum_l G_{kl} * \sigma_{l0}
// 	for(int n=n1;n<=n2;n++){
// 	  for(int k=0;k<Rsize;k++){
// 	    GREEN_TSTP sum(n,ntau,1);
// 	    sum.clear();
// 	    for(int l=0;l<Rsize;l++){
// 	      GREEN gr(nt,ntau,1);
// 	      gr.set_matrixelement(0,0,R,k,l);
// 	      GREEN Sigmak_nl_rest(nt,ntau,1);
// 	      Sigmak_nl_rest.set_matrixelement(0,0,Sigma,ssize+l,j);
// 	      GREEN Sigmak_nl_rest_cc(nt,ntau,1);
// 	      Sigmak_nl_rest_cc.set_matrixelement(0,0,Sigma,j,ssize+l);
// 	      GREEN tmp(nt,ntau,1);
// 	      cntr::convolution_timestep(n,tmp,gr,gr,Sigmak_nl_rest,Sigmak_nl_rest_cc,integration::I<double>(kt),beta,h);
// 	      sum.incr(tmp,1.0);
// 	    }
// 	    GSigma.set_matrixelement(n,k,0,sum);
// 	  }
// 	}
// 	// cc conjugate of [G * Sigma ]_{k0} = [Sigma * G ]_{0k}=\sum_l \Sigma_{0l} * G_{l k}
// 	for(int n=n1;n<=n2;n++){
// 	  for(int k=0;k<Rsize;k++){
// 	    GREEN_TSTP sum(n,ntau,1);
// 	    sum.clear();
// 	    for(int l=0;l<Rsize;l++){
// 	      GREEN gr(nt,ntau,1);
// 	      gr.set_matrixelement(0,0,R,l,k);
// 	      GREEN Sigmak_nl_rest(nt,ntau,1);
// 	      Sigmak_nl_rest.set_matrixelement(0,0,Sigma,i,ssize+l);
// 	      GREEN Sigmak_nl_rest_cc(nt,ntau,1);
// 	      Sigmak_nl_rest_cc.set_matrixelement(0,0,Sigma,ssize+l,i);
// 	      GREEN tmp(nt,ntau,1);
// 	      cntr::convolution_timestep(n,tmp,Sigmak_nl_rest,Sigmak_nl_rest_cc,gr,gr,integration::I<double>(kt),beta,h);
// 	      sum.incr(tmp,1.0);
// 	    }
// 	    SigmaG.set_matrixelement(n,0,k,sum);
// 	  }
// 	}
// 	// [Sigma * G * Sigma ]_{00}=\sum_k  \Sigma_{0k} * [G * \Sigma]_{k 0}
	
// 	GREEN_TSTP SkGkSk(tstp,ntau,1);
// 	SkGkSk.clear();
// 	for(int k=0;k<Rsize;k++){
// 	  GREEN GStmp(nt,ntau,1);
// 	  GStmp.set_matrixelement(0,0,GSigma,k,j);

// 	  GREEN SGtmp(nt,ntau,1);
// 	  SGtmp.set_matrixelement(0,0,SigmaG,j,k);
	    
// 	  GREEN Sigmak_nl_rest(nt,ntau,1);
// 	  Sigmak_nl_rest.set_matrixelement(0,0,Sigma,i,ssize+k);

// 	  GREEN Sigmak_nl_rest_cc(nt,ntau,1);
// 	  Sigmak_nl_rest_cc.set_matrixelement(0,0,Sigma,ssize+k,i);
	    
// 	  GREEN tmp(nt,ntau,1);
// 	  cntr::convolution_timestep(tstp,tmp,Sigmak_nl_rest,Sigmak_nl_rest_cc,GStmp,SGtmp,integration::I<double>(kt),beta,h);
	    
// 	  SkGkSk.incr(tmp,1.0);
// 	}
// 	etatmp.incr(SkGkSk,1.0);
// 	//Add Sigmak_nl_[q]_{dd}
// 	GREEN_TSTP Sigmak_nl_tmp(tstp,ntau,1);
// 	Sigmak_nl_tmp.set_matrixelement(0,0,Sigma,i,j);
// 	etatmp.incr(Sigmak_nl_tmp,1.0);
// 	S.set_matrixelement(tstp,i,j,etatmp);
	
// }



void dyson(GREEN &G,GREEN &G0,GREEN &S,CFUNCTION &eps,int kt,double beta,double h){
	// Temp functions
	int nt=S.nt();
	int size=S.size1();
	int ntau=S.ntau();
	CFUNCTION unity=CFUNCTION(nt,size,size);
	cdmatrix one=cdmatrix::Identity(size,size);
	unity.set_constant(one);

	GREEN SG0=GREEN(nt,ntau,size,-1);
	GREEN G0S=GREEN(nt,ntau,size,-1); 
	cntr::convolution(G0S,G0,G0,S,S,integration::I<double>(kt),beta,h);
	cntr::convolution(SG0,S,S,G0,G0,integration::I<double>(kt),beta,h);
	
	for(int tstp=-1;tstp<=nt;tstp++){
		G0S.right_multiply(tstp,unity,-1.0);
		SG0.right_multiply(tstp,unity,-1.0);
	}

	for(int tstp=-1;tstp<=nt;tstp++){
		GREEN_TSTP tmpR(tstp,ntau,size),tmpL(tstp,ntau,size);
		G0.get_timestep(tstp,tmpR);
		tmpL=tmpR;
		tmpR.right_multiply(eps);
		tmpL.left_multiply(eps);
		G0S.incr_timestep(tstp,tmpR,std::complex<double>(-1.0,0.0));
		SG0.incr_timestep(tstp,tmpL,std::complex<double>(-1.0,0.0));
	}

	cntr::vie2(G,G0S,SG0,G0,integration::I<double>(kt),beta,h);

}


TEST_CASE("Comparison of downfolded and full solution","[Downfold]"){
	int size=3,ssize=1;
	int nt=20;
	int ntau=300;
	double beta=1.0,eps=1e-6;
	double h=0.01;
	int kt=5;

	const std::complex<double> i(0, 1);
	const double pi = std::acos(-1);

	GREEN A=GREEN(nt,ntau,size,-1);
	GREEN A11=GREEN(nt,ntau,ssize,-1);
	GREEN A011=GREEN(nt,ntau,ssize,-1);
	GREEN A022=GREEN(nt,ntau,size-ssize,-1);
	GREEN eta=GREEN(nt,ntau,ssize,-1);
	GREEN Gdownfold(nt,ntau,ssize,-1);

	CFUNCTION epsilon(nt,size,size),epsilon11(nt,size-ssize,size-ssize),epsilon00(nt,ssize,ssize);
	CFUNCTION zero(nt,ssize,ssize);
	cdmatrix nic=cdmatrix::Zero(ssize,ssize);
	zero.set_constant(nic);

	SECTION ("Time-Independent"){
		cdmatrix a(size,size);
		a(0,0)=1.0;
		a(1,1)=2.0;
		a(2,2)=3.0;
		a(0,1)=std::exp(i*pi);
		a(0,2)=std::exp(i*pi);
		a(1,2)=std::exp(i*pi);
		a(1,0)=std::exp(-i*pi);
		a(2,0)=std::exp(-i*pi);
		a(2,1)=std::exp(-i*pi);

		epsilon.set_constant(a);
		// To the full one we put zero, so that the Dyson is properly checked
		epsilon00.set_matrixelement(0,0,epsilon,0,0);
		epsilon11.set_matrixelement(0,0,epsilon,1,1);
		epsilon11.set_matrixelement(0,1,epsilon,1,2);
		epsilon11.set_matrixelement(1,0,epsilon,2,1);
		epsilon11.set_matrixelement(1,1,epsilon,2,2);

		cntr::green_from_H(A,0.0,epsilon,beta,h,"true");
		cntr::green_from_H(A011,0.0,zero,beta,h,"true");
		cntr::green_from_H(A022,0.0,epsilon11,beta,h,"true");

		for(int tstp=-1;tstp<nt;tstp++){
			downfold(tstp,0,0,A022,eta,epsilon);
		}

		dyson(Gdownfold,A011,eta,epsilon00,kt,beta,h);

		A11.set_matrixelement(0,0,A,0,0);
		double err=0.0;
		for(int tstp=-1;tstp<nt;tstp++){
			err+=distance_norm2(tstp,Gdownfold,A11);
		}
		REQUIRE(err<eps);

	}

	SECTION ("Time-Dependent"){
		cdmatrix a(size,size);
		for(int tstp=-1;tstp<=nt;tstp++){
		  a(0,0)=1.0+sin(1.0*h*double(tstp));
		  a(1,1)=2.0+sin(1.0*h*double(tstp));
		  a(2,2)=3.0+sin(1.0*h*double(tstp));
		  a(0,1)=std::exp(i*pi*double(tstp)*h);
		  a(0,2)=std::exp(i*pi*double(tstp)*h);
		  a(1,2)=std::exp(i*pi*double(tstp)*h);
		  a(1,0)=std::exp(-i*pi*double(tstp)*h);
		  a(2,0)=std::exp(-i*pi*double(tstp)*h);
		  a(2,1)=std::exp(-i*pi*double(tstp)*h);
		  epsilon.set_value(tstp,a);
		}
		
		// To the full one we put zero, so that the Dyson is properly checked
		epsilon00.set_matrixelement(0,0,epsilon,0,0);
		epsilon11.set_matrixelement(0,0,epsilon,1,1);
		epsilon11.set_matrixelement(0,1,epsilon,1,2);
		epsilon11.set_matrixelement(1,0,epsilon,2,1);
		epsilon11.set_matrixelement(1,1,epsilon,2,2);

		
		cntr::green_from_H(A,0.0,epsilon,beta,h,"true");
		cntr::green_from_H(A011,0.0,zero,beta,h,"true");
		cntr::green_from_H(A022,0.0,epsilon11,beta,h,"true");
		
		for(int tstp=-1;tstp<nt;tstp++){
			downfold(tstp,0,0,A022,eta,epsilon);
		}
		
		dyson(Gdownfold,A011,eta,epsilon00,kt,beta,h);

		A11.set_matrixelement(0,0,A,0,0);
		double err=0.0;
		for(int tstp=-1;tstp<nt;tstp++){
			err+=distance_norm2(tstp,Gdownfold,A11);
		}
		REQUIRE(err<eps);
	}

	SECTION ("Interacting"){
	  
	  cdmatrix a(size,size);
	  for(int tstp=-1;tstp<=nt;tstp++){
	    a(0,0)=1.0+0.0*sin(1.0*h*double(tstp));
	    a(1,1)=2.0+0.0*sin(1.0*h*double(tstp));
	    a(2,2)=3.0+0.0*sin(1.0*h*double(tstp));
	    a(0,1)=std::exp(i*pi*0.0*double(tstp)*h);
	    a(0,2)=std::exp(i*pi*0.0*double(tstp)*h);
	    a(1,2)=std::exp(i*pi*0.0*double(tstp)*h);
	    a(1,0)=std::exp(-i*pi*0.0*double(tstp)*h);
	    a(2,0)=std::exp(-i*pi*0.0*double(tstp)*h);
	    a(2,1)=std::exp(-i*pi*0.0*double(tstp)*h);
	    epsilon.set_value(tstp,a);
	  }
	  
	  // To the full one we put zero, so that the Dyson is properly checked
	  epsilon00.set_matrixelement(0,0,epsilon,0,0);
	  epsilon11.set_matrixelement(0,0,epsilon,1,1);
	  epsilon11.set_matrixelement(0,1,epsilon,1,2);
	  epsilon11.set_matrixelement(1,0,epsilon,2,1);
	  epsilon11.set_matrixelement(1,1,epsilon,2,2);
	  
	  GREEN Sigma=GREEN(nt,ntau,size,-1);
	  GREEN GSigma=GREEN(nt,ntau,size-ssize,ssize,-1);
	  GREEN SigmaG=GREEN(nt,ntau,ssize,size-ssize,-1);
	  // Define self-energy
	  cntr::green_from_H(Sigma,0.0,epsilon,beta,h,"true");
	  GREEN SigmaRest=GREEN(nt,ntau,size-ssize,-1);
	  std::vector<int> rest11{0,0,1,1};
	  std::vector<int> rest12{0,1,0,1};
	  std::vector<int> full11{1,1,2,2};
	  std::vector<int> full12{1,2,1,2};
	  SigmaRest.set_submatrix(rest11,rest12,Sigma,full11,full12);
	  
	  // Solve with Vie since dyson solvers are not accurate enought for these kind of problems
	  // Full Green's function
	  CFUNCTION zero3(nt,size,size);
	  zero3.set_zero();
	  GREEN G0eps=GREEN(nt,ntau,size,-1);
	  cntr::green_from_H(G0eps,0.0,epsilon,beta,h,"true");
	  dyson(A,G0eps,Sigma,zero3,kt,beta,h);
	  cntr::green_from_H(A011,0.0,zero,beta,h,"true");
	  
	  GREEN GRest=GREEN(nt,ntau,size-ssize,-1);
	  CFUNCTION zero2(nt,size-ssize,size-ssize);
	  zero2.set_zero();
	  GREEN G0epsRest=GREEN(nt,ntau,size-ssize,-1);
	  cntr::green_from_H(G0epsRest,0.0,epsilon11,beta,h,"true");
	  dyson(GRest,G0epsRest,SigmaRest,zero2,kt,beta,h);
	  
	  for(int tstp=-1;tstp<=nt;tstp++){
	    downfold3(tstp,0,0,GRest,eta,epsilon,Sigma,GSigma,SigmaG,kt,beta,h);
	  }
	  
	  dyson(Gdownfold,A011,eta,epsilon00,kt,beta,h);
	  
	  A11.set_matrixelement(0,0,A,0,0);
	  double err=0.0;
	  for(int tstp=-1;tstp<nt;tstp++){
	    err+=distance_norm2(tstp,Gdownfold,A11);
	  }
	  REQUIRE(err<eps);

	}
}
