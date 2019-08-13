#include "catch.hpp"
#include "cntr.hpp"

#define CNTR_USE_OMP
#define CNTR_USE_MPI
#include "nca_impurity_model.hpp"
#include "read_inputfile.hpp"

#define CPLX std::complex<double>  
using namespace std;
 


// ME NOTE : Following routines are copied from nca_Impurity_Problem_L1_Hirsch.hpp, line 208ff.
cdouble phi(cdouble t, cdouble tprime, double beta, double w0, double g) {
	cdouble ii;
	ii.real(0.0);
	ii.imag(1.0);
	return (g/w0)*(g/w0)*( cosh((beta/2 - (ii*t- ii*tprime))*w0) - cosh(beta*w0/2) )/sinh(beta*w0/2);
	// t must be later on contout than tprime
}
void get_phonon_correction(cntr::herm_matrix<double> &Phonon_factor,double g,double w0,double beta,double h){
	// g(t)=g
	// -ii * < TC e^{bdag-b}   e^{b-bdag} >
	// Matsubara - < Ttau e^{bdag-b}   e^{bdag-b} >
	// dU_ph[tstp=-1...1]	
	int ntau=Phonon_factor.ntau();
	int nt=Phonon_factor.nt();
	cdouble phi1,zt,zt1;
	cdouble ii;
	ii.real(0.0);
	ii.imag(1.0);
	
	for (int tstp=-1; tstp<=nt; tstp++) {
		
		if(tstp==-1){
			//Matsubara part
			for(int m=0;m<=ntau;m++){
				zt1=0.0;
				zt.real(0.0);
				zt.imag(-m*beta/ntau);
				phi1=phi(zt,zt1,beta,w0,g);
				Phonon_factor.set_mat(m,-exp(phi1));
			}
		}else{
			// retarded part: C(tstp,m)
			zt1.real(h*tstp);
			zt1.imag(0.0);
			for(int m=0;m<=tstp;m++){
				zt.real(m*h);
				zt.imag(0.);

				// ME NOTE !! : The choce of arguments seems to be inconsonsistent with
				// the previous definition in ncanew::update_g_phonon_dispatch<double,1>
				// which was Bles(t,t') = exp(sig*phi(t',t)), Bgtr(t,t') = exp(sig*phi(t,t'))
				// plus ii factor
				// (cf. nca_phonon_correction.hpp, line 39ff. )
				// the old version is also consistent with the symmetry of the function, s=+/- 1
				// B(t,t') = < TC exp(-s*ii*P(t)) exp(s*ii*P(t')) >
				// using that the HAmiltonian is symmetric under X->-X and P->-P
				
				// FROM CURRENT NCA_IMPURITY_PROBLEM_L1_HIRSCH
				//cdouble tmp_gtr=phi(zt,zt1,beta,w0,g); // t later;
				//cdouble tmp_les=phi(zt1,zt,beta,w0,g);
				//Phonon_factor.set_ret(tstp,m,-ii*(exp(tmp_gtr)-exp(tmp_les))); // HOW DO THIS PROPERLY?
				
				// MY SUGGESTION 
				cdouble tmp_gtr=phi(zt1,zt,beta,w0,g); // t later;
				cdouble tmp_les=phi(zt,zt1,beta,w0,g);
				Phonon_factor.set_ret(tstp,m,-ii*(exp(tmp_gtr)-exp(tmp_les))); // HOW DO THIS PROPERLY?


			}                
			// lesser part: C(m,tstp)
			//zt1.real()=h*tstp;
			//zt1.imag()=0.;
			for(int m=0;m<=tstp;m++){
				zt.real(m*h);
				zt.imag(0.0);
				cdouble tmp_les=phi(zt1,zt,beta,w0,g);
				//Phonon_factor_.lesptr(m,tstp,-Icplx*phi(zt1,zt,beta,w0,g));
				Phonon_factor.set_les(m,tstp,-ii*exp(tmp_les));
			}
			// tv part: 
			//zt1=CPLX(h*tstp,0);
			for(int m=0;m<=ntau;m++){
				//zt=CPLX(t*h,0);
				zt.real(0.0);
				zt.imag(-m*beta/ntau);
				phi1=phi(zt,zt1,beta,w0,g);
				//Phonon_factor_.set_tvptr(tstp,m,-Icplx*phi1);
				Phonon_factor.set_tv(tstp,m,-ii*exp(phi1));
			}
			
		}
		
	}
	
}


TEST_CASE("Lang-firsov phonon factors","[LF]"){
	int nt=50,tstp;
	int ntau=50;
	double beta=2.0;
	double wph=1.0;
	double gph=0.7;
	double h=0.05;
	cntr::herm_matrix<double> Phonon_factor(nt,ntau,1,1);
	cntr::herm_matrix<double> Delta(nt,ntau,1,-1);
	cntr::herm_matrix<double> DeltaLF(nt,ntau,1,-1);
	cntr::herm_matrix<double> DeltaLF_2(nt,ntau,1,-1);
	cntr::full_matrix<double> Dfull(nt,ntau,1,-1);
	cntr::full_matrix<double> DfullLF(nt,ntau,1,-1);
	cntr::full_matrix<double> DfullLF_2(nt,ntau,1,-1);
	
	cntr::green_equilibrium_bethe(Delta,beta,h);
	
	for(tstp=-1;tstp<=nt;tstp++) Dfull.set_timestep(tstp,Delta,Delta);
	DfullLF=Dfull;
	for(tstp=-1;tstp<=nt;tstp++) ncanew::update_g_phonon_dispatch<double,1>(tstp,-1,gph,wph,beta,h,DfullLF);
	for(tstp=-1;tstp<=nt;tstp++) DeltaLF_2.set_timestep(tstp,DfullLF);


	get_phonon_correction(Phonon_factor,gph,wph,beta,h);
	for(tstp=-1;tstp<=nt;tstp++) Bubble2(tstp,DeltaLF,Delta,Phonon_factor);
	for(tstp=-1;tstp<=nt;tstp++) DfullLF_2.set_timestep(tstp,DeltaLF,DeltaLF);
		
	double err=0.0;
	for(tstp=-1;tstp<=nt;tstp++){
		
		// cout << "NT: " << nt << " tstp: " << tstp;
		err+=cntr::distance_norm2(tstp,DeltaLF,DeltaLF_2);		
		// cout << " |A1-A| " << err;
		err+=cntr::distance_norm2_ret(tstp,DeltaLF,DeltaLF_2);		
		// cout << " |A1-A| ret " << err;
		err+=cntr::distance_norm2_tv(tstp,DeltaLF,DeltaLF_2);		
		// cout << " |A1-A| tv " << err;
		err+=cntr::distance_norm2_les(tstp,DeltaLF,DeltaLF_2);		
		// cout << " |A1-A| les " << err;
		// cout << endl;
	}
	REQUIRE(err<1e-7);
}




