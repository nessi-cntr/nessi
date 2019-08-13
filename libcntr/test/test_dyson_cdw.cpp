#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>


#include "cntr.hpp"
#include "read_inputfile.hpp"

#define CPLX std::complex<double>  
using namespace std;
 

/*///////////////////////////////////////////////////////////////////////////////////////


This is a minimal test program to solve the self-consistent equation 

[ii*d/dt + mu - f(t) ] G(t,t') - \int_C dt1 \Delta(t,t1) G(t1,t') = \delta_C(t,t')

with self-consistent hybridization function 

\Delta(t,t') = G(t,t')

and arbitrary function h(t)

///////////////////////////////////////////////////////////////////////////////////////*/ 


int main(int argc,char *argv[]){
	int itermax,iter_rtime;
	int nt,ntau,kt,iter,size,sig,tstp;
	double mu,beta,h,errmax,err,v1,v2,aa;
	cntr::herm_matrix<double> G,Delta,Ga,Gb,Da,Db;
	cntr::herm_matrix<double> fDa,Dfa,fDb,Dfb,Ga0,Gb0,Ga1,Gb1;
	cntr::herm_matrix<double> Gbethe;
	cntr::function<double> f,fa,fb,fa0,fb0,dfa,dfb;
	
	try{
	////////////////////////////////////////////////////////////////////////////////////	
	// (II) READ INPUT
	{
		if(argc<2) throw("COMMAND LINE ARGUMENT MISSING");
		// This is some hand-made way of reading in input from a text-file: argv[1] 
		// (the first command-line argument)
		//
		// void find_param(char *file,const char *flag,T &x)
		// defined in code/utils/read_inputfile.hpp:
		// scans the file for a line beginning with flag and x to the second entry in that line
		// meaning of parameters explained below, at their first occurence
		find_param(argv[1],"__nt=",nt);
		find_param(argv[1],"__ntau=",ntau);
		find_param(argv[1],"__beta=",beta);
		find_param(argv[1],"__h=",h);    
		//find_param(argv[1],"__mu=",mu);
		find_param(argv[1],"__itermax=",itermax);
		find_param(argv[1],"__errmax=",errmax);
		find_param(argv[1],"__iter_rtime=",iter_rtime);
        find_param(argv[1],"__kt=",kt);
        find_param(argv[1],"__v1=",v1);
        find_param(argv[1],"__v2=",v2);
        find_param(argv[1],"__aa=",aa);
	}	
	////////////////////////////////////////////////////////////////////////////////////	
	// (IV) INITIALIZE GREEN'S FUNCTIONS
	{
		// class herm_matrix<double>: defined in code/cntr/cntr.hpp
		// Green's function object that stores retarded, matsubara, lesser, and mixed 
		// (tau,t) component
		size=1;
		sig=-1;
		G=cntr::herm_matrix<double>(nt,ntau,size,sig);
		// parameters of constructor:
		// nt: real-time branch has discrete time-points 0....nt
		// ntau: imag-time branch has discrete time-points 0....ntau
		// size: orbital dimension, i.e., G(t,t') is a size*size matrix
		// sig: -1 for fermionic G, 1 for bosonic G 
		// By call f the constructor, G is initialized with 0
		Delta=cntr::herm_matrix<double>(nt,ntau,size,sig);
		Gbethe=cntr::herm_matrix<double>(nt,ntau,size,sig); 
		
		////////////////////////////////////
		// Gbethe will be the exact solution (which know for this test, what a surprise!):
		// The exact GF (for f=0) is given by Eq.(27) in review, with A(w) ~ sqrt(4-w^2)
		// The call to compute this is definined in code/cntr/cntr.hpp
		// there is a formalism to set up GF for any density-of states, but this is not important
		cntr::green_equilibrium_bethe(Gbethe,beta,h);

		Da=Gbethe;
		Db=Gbethe;
		Ga=Gbethe;
		Gb=Gbethe;
		Ga1=Gbethe;
		Gb1=Gbethe;
		Gb0=Gbethe;
		Ga0=Gbethe;
		Dfa=Gbethe;
		fDa=Gbethe;
		Dfb=Gbethe;
		fDb=Gbethe;

		// beta: inverse temperature
		// h: timestep
		////////////////////////////////////
		f=cntr::function<double>(nt,size);
		fa=cntr::function<double>(nt,size);
		fb=cntr::function<double>(nt,size);
		dfa=cntr::function<double>(nt,size);
		dfb=cntr::function<double>(nt,size);
		fa0=cntr::function<double>(nt,size);
		fb0=cntr::function<double>(nt,size);
		// function<double> is defined in defined in code/cntr/cntr.hpp
		// function on the contour, i.e., it stores one value for the 
		// matsubara branch (timespep=-1), and values on timesteps 0...nt
		// each entry is a size*size matrix
		mu=0.0;
		for(tstp=-1;tstp<=nt;tstp++) f[tstp]=1.0;
		for(tstp=-1;tstp<=nt;tstp++) fa[tstp]=v2*sqrt(2.0)*(1.0+aa*cos(0.2*h*tstp))/(1.0+aa);
		for(tstp=-1;tstp<=nt;tstp++) fb[tstp]=-v2*sqrt(2.0)*(1.0+aa*cos(0.2*h*tstp))/(1.0+aa);
		fa[-1]=v1*sqrt(2.0);
		fb[-1]=-v1*sqrt(2.0);
		for(tstp=-1;tstp<=nt;tstp++) fa0[tstp]=fa[-1];
		for(tstp=-1;tstp<=nt;tstp++) fb0[tstp]=fb[-1];
		for(tstp=-1;tstp<=nt;tstp++) dfa[tstp]=-(fa[tstp]-fa0[tstp]);
		for(tstp=-1;tstp<=nt;tstp++) dfb[tstp]=-(fb[tstp]-fb0[tstp]);
		//fa0=fa;
		//fb0=fb;
		
	}
	////////////////////////////////////////////////////////////////////////////////////
	// FIRST TEST: NO SELF-CONSISTENCY
	// SET Delta=Gbethe, f=0, and solve [id/dt + mu  - f - Delta ] G = 1
	// => result should be Gbethe
	/*
	{
		Delta=Gbethe;
		// STEP 1) solve equation on the Matsubara branch:
		cntr::dyson_mat(G,Delta,mu,f,beta);
		// STEP 2) solve equation on the first kt timesteps, where kt is the integration order:
		// use kt=5 unless there is a good reason not to do so
		cntr::dyson_start(G,mu,f,Delta,integration::I<double>(kt),beta,h);
		// integration::I<double>(kt) is a strange way of setting the integration order
		// (to be replaced by simpler call to kt in later revision of code)
		// STEP 3)
		// solve equation successively on timesteps kt+1,...,nt:
		for(tstp=kt+1;tstp<=nt;tstp++){
			cntr::dyson_timestep(tstp,G,mu,f,Delta,integration::I<double>(kt),beta,h);
		}
		// THAT WAS IT!
		
		// COMPARE Gbethe and the solution G
		// (the difference norm is defined in )
		for(tstp=-1;tstp<=nt;tstp++){
			err=cntr::distance_norm2(tstp,G,Gbethe);
			cout << "timestep: " << tstp << " |G-Gbethe|= " << err << endl;
		}
	}
	*/
	////////////////////////////////////////////////////////////////////////////////////
	// SECOND TEST: SELF-CONSISTET SOLUTION OF SAME THING
	// SET DELTA=0 IN THE BEGINNING, FIND DELTA
	{
		bool matsubara_converged=false;
		Delta.clear();  // set zero
		// STEP 1) solve equation on the Matsubara branch iteratively:
		tstp=-1;
		cntr::herm_matrix<double> gtemp; // to store last timestep
		gtemp=cntr::herm_matrix<double>(kt,ntau,size,-1); // contains only kt timestep (nt-parameter=-1)
		for(int iter=0;iter<=itermax;iter++){
			cntr::dyson_mat(Ga,Da,mu,fa,beta);	
			cntr::dyson_mat(Gb,Db,mu,fb,beta);	
			err=cntr::distance_norm2(tstp,Ga,gtemp);
			cout << "iteration : " << iter << " |G-Gtemp|= " << err << endl;
			if(err<errmax){
				matsubara_converged=true;
				break;
			}
			gtemp.set_timestep(-1,Ga); // save timestep -1 to gtemp
			// self-consistency:
			Da.set_timestep(-1,Gb); // save timestep -1 to gtemp
			Db.set_timestep(-1,Ga); // save timestep -1 to gtemp
		}
		if(!matsubara_converged){
			cout << "cannot cnverge matsubara; end here " << endl;
			// should end here ....
		}
		// STEP 2): first kt timesteps: do iteration like above:
		for(int iter=0;iter<=itermax;iter++){
			cntr::dyson_start(Ga,mu,fa,Da,integration::I<double>(kt),beta,h);
			cntr::dyson_start(Gb,mu,fb,Db,integration::I<double>(kt),beta,h);
			err=0.0;
			for(tstp=0;tstp<=kt;tstp++) err += cntr::distance_norm2(tstp,Ga,gtemp);
			cout << "START iteration : " << iter << " |G-Gtemp|= " << err << endl;
			if(err<errmax){
				matsubara_converged=true;
				break;
			}
			for(tstp=0;tstp<=kt;tstp++) gtemp.set_timestep(tstp,Ga); // save timestep 0...kt to gtemp
			// self-consistency:
			for(tstp=0;tstp<=kt;tstp++) Da.set_timestep(tstp,Gb); // save timestep 0...kt to gtemp
			for(tstp=0;tstp<=kt;tstp++) Db.set_timestep(tstp,Ga); // save timestep 0...kt to gtemp
		}
		// STEP 3): time-stepping kt+1 ... nt
		for(tstp=kt+1;tstp<=nt;tstp++){
			// extrapoaltion of Delta by one timestep:
			cntr::extrapolate_timestep(tstp-1,Da,integration::I<double>(kt));
			cntr::extrapolate_timestep(tstp-1,Db,integration::I<double>(kt));
			// iteration on the timestep, fixed iteration number, no error check!
			cout << " ... at timestep " << tstp << endl;
			for(int iter=1;iter<=iter_rtime;iter++){
				cntr::dyson_timestep(tstp,Gb,mu,fb,Db,integration::I<double>(kt),beta,h);
				cntr::dyson_timestep(tstp,Ga,mu,fa,Da,integration::I<double>(kt),beta,h);
				Da.set_timestep(tstp,Gb); 
				Db.set_timestep(tstp,Ga); 
			}
		}
		Gb.print_to_file("gb.out");
		Ga.print_to_file("ga.out");
		
	/////////////////////////////////////////////////////////////////////////////////////////////////////
		Da=Gbethe;
		Db=Gbethe;
		matsubara_converged=false;
		// STEP 1) solve equation on the Matsubara branch iteratively:
		tstp=-1;
		//cntr::herm_matrix<double> gtemp; // to store last timestep
		gtemp=cntr::herm_matrix<double>(kt,ntau,size,-1); // contains only kt timestep (nt-parameter=-1)
		for(int iter=0;iter<=itermax;iter++){
			cntr::dyson_mat(Ga0,Da,mu,fa0,beta);	
			cntr::dyson_mat(Gb0,Db,mu,fb0,beta);	
			err=cntr::distance_norm2(tstp,Ga0,gtemp);
			cout << "iteration : " << iter << " |G-Gtemp|= " << err << endl;
			if(err<errmax){
				matsubara_converged=true;
				break;
			}
			gtemp.set_timestep(-1,Ga0); // save timestep -1 to gtemp
			// self-consistency:
			Da.set_timestep(-1,Gb0); // save timestep -1 to gtemp
			Db.set_timestep(-1,Ga0); // save timestep -1 to gtemp
		}
		if(!matsubara_converged){
			cout << "cannot cnverge matsubara; end here " << endl;
			// should end here ....
		}
		Ga1.set_timestep(-1,Ga0);
		Gb1.set_timestep(-1,Gb0);
		Dfa.set_timestep_zero(-1);
		Dfb.set_timestep_zero(-1);
		fDa.set_timestep_zero(-1);
		fDb.set_timestep_zero(-1);
		
		
		// STEP 2): first kt timesteps: do iteration like above:
		for(int iter=0;iter<=itermax;iter++){
			cntr::dyson_start(Ga0,mu,fa0,Da,integration::I<double>(kt),beta,h);
			cntr::dyson_start(Gb0,mu,fb0,Db,integration::I<double>(kt),beta,h);
			for(tstp=0;tstp<=kt;tstp++){
				Dfa.set_timestep(tstp,Ga0);
				fDa.set_timestep(tstp,Ga0);
				Dfb.set_timestep(tstp,Gb0);
				fDb.set_timestep(tstp,Gb0);
			
				Dfa.right_multiply(tstp,dfa);
				fDa.left_multiply(tstp,dfa);
				Dfb.right_multiply(tstp,dfb);
				fDb.left_multiply(tstp,dfb);
			
			}
			cntr::vie2_start(Ga1,Dfa,fDa,Ga0,integration::I<double>(kt),beta,h);
			cntr::vie2_start(Gb1,Dfb,fDb,Gb0,integration::I<double>(kt),beta,h);
			
			
			err=0.0;
			for(tstp=0;tstp<=kt;tstp++) err += cntr::distance_norm2(tstp,Ga0,gtemp);
			cout << "START iteration : " << iter << " |G-Gtemp|= " << err << endl;
			if(err<errmax){
				matsubara_converged=true;
				break;
			}
			for(tstp=0;tstp<=kt;tstp++) gtemp.set_timestep(tstp,Ga0); // save timestep 0...kt to gtemp
			// self-consistency:
			for(tstp=0;tstp<=kt;tstp++) Da.set_timestep(tstp,Gb1); // save timestep 0...kt to gtemp
			for(tstp=0;tstp<=kt;tstp++) Db.set_timestep(tstp,Ga1); // save timestep 0...kt to gtemp
		}
		// STEP 3): time-stepping kt+1 ... nt
		for(tstp=kt+1;tstp<=nt;tstp++){
			// extrapoaltion of Delta by one timestep:
			cntr::extrapolate_timestep(tstp-1,Da,integration::I<double>(kt));
			cntr::extrapolate_timestep(tstp-1,Db,integration::I<double>(kt));
			// iteration on the timestep, fixed iteration number, no error check!
			cout << " ... at timestep " << tstp << endl;
			for(int iter=1;iter<=iter_rtime;iter++){
				cntr::dyson_timestep(tstp,Gb0,mu,fb0,Db,integration::I<double>(kt),beta,h);
				cntr::dyson_timestep(tstp,Ga0,mu,fa0,Da,integration::I<double>(kt),beta,h);
				Dfa.set_timestep(tstp,Ga0);
				fDa.set_timestep(tstp,Ga0);
				Dfb.set_timestep(tstp,Gb0);
				fDb.set_timestep(tstp,Gb0);
				Dfa.right_multiply(tstp,dfa);
				fDa.left_multiply(tstp,dfa);
				Dfb.right_multiply(tstp,dfb);
				fDb.left_multiply(tstp,dfb);
				cntr::vie2_timestep(tstp,Ga1,Dfa,fDa,Ga0,integration::I<double>(kt),beta,h);
				cntr::vie2_timestep(tstp,Gb1,Dfb,fDb,Gb0,integration::I<double>(kt),beta,h);
				
				Da.set_timestep(tstp,Gb1); 
				Db.set_timestep(tstp,Ga1); 
			}
		}
		Gb0.print_to_file("gb0.out");
		Ga0.print_to_file("ga0.out");
		
		/*
		////////////////////////////////////////////////////////
		// recompute real time part using vie2
		Ga1=Ga;
		Gb1=Gb;
		for(tstp=0;tstp<=nt;tstp++) Da.set_timestep(tstp,Gbethe);
		for(tstp=0;tstp<=nt;tstp++) Db.set_timestep(tstp,Gbethe);
		Ga0.set_timestep(-1,Ga);
		Gb0.set_timestep(-1,Gb);
		Dfa.set_timestep_zero(-1);
		Dfb.set_timestep_zero(-1);
		fDa.set_timestep_zero(-1);
		fDb.set_timestep_zero(-1);
		// STEP 2): first kt timesteps: do iteration like above:
		for(int iter=0;iter<=itermax;iter++){
			cntr::dyson_start(Ga0,mu,fa0,Da,integration::I<double>(kt),beta,h);
			cntr::dyson_start(Gb0,mu,fb0,Db,integration::I<double>(kt),beta,h);
			
			for(tstp=0;tstp<=kt;tstp++){
				Dfa.set_timestep(tstp,Ga0);
				fDa.set_timestep(tstp,Ga0);
				Dfb.set_timestep(tstp,Gb0);
				fDb.set_timestep(tstp,Gb0);
			
				Dfa.right_multiply(tstp,dfa);
				fDa.left_multiply(tstp,dfa);
				Dfb.right_multiply(tstp,dfb);
				fDb.left_multiply(tstp,dfb);
			
			}
			
			cntr::vie2_start(Ga1,Dfa,fDa,Ga0,integration::I<double>(kt),beta,h);
			cntr::vie2_start(Gb1,Dfb,fDb,Gb0,integration::I<double>(kt),beta,h);
			
			err=0.0;
			for(tstp=0;tstp<=kt;tstp++) err += cntr::distance_norm2(tstp,Ga,gtemp);
			cout << "START iteration : " << iter << " |G-Gtemp|= " << err << endl;
			if(err<errmax){
				matsubara_converged=true;
				break;
			}
			for(tstp=0;tstp<=kt;tstp++) gtemp.set_timestep(tstp,Ga0); // save timestep 0...kt to gtemp
			// self-consistency:
			for(tstp=0;tstp<=kt;tstp++) Da.set_timestep(tstp,Gb0); // save timestep 0...kt to gtemp
			for(tstp=0;tstp<=kt;tstp++) Db.set_timestep(tstp,Ga0); // save timestep 0...kt to gtemp
		}
		// STEP 3): time-stepping kt+1 ... nt
		for(tstp=kt+1;tstp<=nt;tstp++){
			// extrapoaltion of Delta by one timestep:
			cntr::extrapolate_timestep(tstp-1,Da,integration::I<double>(kt));
			cntr::extrapolate_timestep(tstp-1,Db,integration::I<double>(kt));
			// iteration on the timestep, fixed iteration number, no error check!
			cout << " ... at timestep " << tstp << endl;
			for(int iter=1;iter<=iter_rtime;iter++){
				cntr::dyson_timestep(tstp,Gb0,mu,fb0,Db,integration::I<double>(kt),beta,h);
				cntr::dyson_timestep(tstp,Ga0,mu,fa0,Da,integration::I<double>(kt),beta,h);
				Dfa.set_timestep(tstp,Ga0);
				fDa.set_timestep(tstp,Ga0);
				Dfb.set_timestep(tstp,Gb0);
				fDb.set_timestep(tstp,Gb0);
				Dfa.right_multiply(tstp,dfa);
				fDa.left_multiply(tstp,dfa);
				Dfb.right_multiply(tstp,dfb);
				fDb.left_multiply(tstp,dfb);
				cntr::vie2_timestep(tstp,Ga1,Dfa,fDa,Ga0,integration::I<double>(kt),beta,h);
				cntr::vie2_timestep(tstp,Gb1,Dfb,fDb,Gb0,integration::I<double>(kt),beta,h);
				
				
				
				Da.set_timestep(tstp,Gb0); 
				Db.set_timestep(tstp,Ga0); 
			}
		}
		*/
		

		Gb1.print_to_file("gb1.out");
		Ga1.print_to_file("ga1.out");
		//Ga1=Ga0;
		//Gb1=Gb0;
		
		for(tstp=-1;tstp<=nt;tstp++){
			double err;
			cout << "NT: " << nt << " tstp: " << tstp;
			err=cntr::distance_norm2(tstp,Ga1,Ga);		
			cout << " |A1-A| " << err;
			err=cntr::distance_norm2_ret(tstp,Ga1,Ga);		
			cout << " |A1-A| ret " << err;
			err=cntr::distance_norm2_tv(tstp,Ga1,Ga);		
			cout << " |A1-A| tv " << err;
			err=cntr::distance_norm2_les(tstp,Ga1,Ga);		
			cout << " |A1-A| les " << err;
			cout << endl;
			
			cout << "NT: " << nt << " tstp: " << tstp;
			err=cntr::distance_norm2(tstp,Gb1,Gb);		
			cout << " |B1-B| " << err;
			err=cntr::distance_norm2_ret(tstp,Gb1,Gb);		
			cout << " |B1-B| ret " << err;
			err=cntr::distance_norm2_tv(tstp,Gb1,Gb);		
			cout << " |B1-B| tv " << err;
			err=cntr::distance_norm2_les(tstp,Gb1,Gb);		
			cout << " |B1-B| les " << err;
			cout << endl;
		}

		
		
		//for(tstp=-1;tstp<=nt;tstp++){
	// err=cntr::distance_norm2(tstp,G,Gbethe);
		//	cout << "timestep: " << tstp << " |G-Gbethe|= " << err << endl;
		//}

	}
	} // try
	catch(char *message){
	    cerr << "exception\n**** " << message << " ****" << endl;
		cerr << "CDMFT input_file [ --test ]\n" << endl;
	}
	catch(...){
	    cerr << "unspecified exception " << endl;
		cerr << "\nCDMFT input_file [ --test ]\n" << endl;
	}
   return 0;
}




