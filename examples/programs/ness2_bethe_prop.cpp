#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <cmath>
#include <sstream>


// INCLUDES ------------------------------------------------------------

#include "cntr/cntr.hpp"
#include "cntr/utils/read_inputfile.hpp"


// NAMESPACES ------------------------------------------------------------

using namespace std;
using namespace cntr;


// UTILITIES --------------------------------------------------------------

// compute Kernel K = G + Sigma, with a second order self-energy
// Sigma(t,t')=U(t)Gweiss(t,t')Gweiss(t,t')Gweiss(t',t)U(t')
// for full Green's functions
// use Gweiss=G for self-consitent perturbation theory
void get_Sigma_timestep(int tstp, GREEN &Sigma, GREEN &Gweiss, CFUNC &U){
  GREEN_TSTP W(tstp,Sigma.ntau(),Sigma.size1(),BOSON); // temporary variable;
  // Note: a bubble of two Fermion GF is bosonic
  cntr::Bubble1(tstp,W,Gweiss,Gweiss); // W(t,t')  set to chi=ii*G(t,t')G(t',t);
  W.left_multiply(tstp,U);
  W.right_multiply(tstp,U); // now W  set to W(t,t')=U(t)chi(t,t')U(t')
  cntr::Bubble2(tstp,Sigma,Gweiss,W); // Sigma(t,t') is set to Sigma=ii*G(t,t')W(t',t);
  Sigma.smul(tstp,-1.0); // the final -1 sign.
}

void get_K_timestep(int tstp, GREEN &K, GREEN &Sigma, GREEN &G){
  K.set_timestep(tstp,Sigma); // K=Sigma
  K.incr_timestep(tstp,G,1.0);  // K += G
}


//
// MAIN #######################################################################
//

int main(int argc, char** argv){
    
    // PARAMETERS ---------------------------------------------------------------
    // General params:
    int size = 1;   // orbital size, here 1
    int nt;    // number of timesteps
    double h;  // timestep h
    double beta; //inverse temperature
    double U, mu, err;
    // NESSi params:
    int MatsMaxIter,CorrectorSteps,BootstrapMaxIter;
    double MatsMaxErr,BootstrapMaxErr;
    int SolveOrder = MAX_SOLVE_ORDER;
    int ntau;
    int tmax,iter,tstp,ipt_flag;
    char *flin,*flout;
    //  read params from input file
    {
        cout << endl;
        cout << " reading input file ..." << endl;
        cout << endl;
        if (argc < 3) {
            std::cerr << "Usuage: ness2_bethe_prop.x param-file out-file" << std::endl;
            throw("COMMAND LINE ARGUMENT MISSING");
        }
        flin=argv[1];
        flout=argv[2];
        find_param(flin, "__h=", h);
        find_param(flin, "__nt=", nt);
        find_param(flin, "__beta=", beta);
        find_param(flin, "__U=", U);
        find_param(flin, "__mu=", mu);
        find_param(flin,"__ntau=",ntau);
        find_param(flin,"__ipt_flag=",ipt_flag);
        find_param(flin,"__BootstrapMaxIter=",BootstrapMaxIter);
        find_param(flin,"__BootstrapMaxErr=",BootstrapMaxErr);
        find_param(flin,"__CorrectorSteps=",CorrectorSteps);
    }
    
    // INITIALIZE GREENS FUNCTIONS ------------------------------------
    
    GREEN G(nt,ntau,1,FERMION);
    GREEN Gweiss(nt,ntau,1,FERMION);
    GREEN Sigma(nt,ntau,1,FERMION);
    GREEN K(nt,ntau,1,FERMION);
    CFUNC Uvec(nt,1);
    CFUNC eloc(nt,1); // =0
    
    //start output (will write matsubara iterations )
    hid_t file_id = open_hdf5_file(flout);
    int itermax=BootstrapMaxIter;
    // INITIALIZE U AND BATH ------------------------------------
    for(tstp=-1;tstp<=nt;tstp++){
      cdmatrix ut = cdmatrix::Constant(1, 1, U);
      Uvec.set_value(tstp,ut);
    }
    
    // init G on matsubara branch with noninteracting GF:
    cntr::green_equilibrium_mat_bethe(G,beta);
    Gweiss.set_timestep(-1,G);
    Sigma.set_timestep_zero(-1);
    get_K_timestep(-1,K,Sigma,G);
    
    // MATSUBARA:
    {
      tstp=-1;
      GREEN_TSTP gtemp(tstp,ntau,1); // to store last timestep
      for(iter=0;iter<=BootstrapMaxIter;iter++){
	gtemp.set_timestep(tstp,G); //store last iteration
	if(ipt_flag){
	  cntr::dyson_mat(Gweiss,mu,eloc,G,beta,SolveOrder, CNTR_MAT_FIXPOINT);
	}else{
	  Gweiss.set_timestep(-1,G);
	}	    
	get_Sigma_timestep(tstp,Sigma,Gweiss,Uvec);
	get_K_timestep(tstp,K,Sigma,G);
	cntr::dyson_mat(G,mu,eloc,K,beta,SolveOrder, CNTR_MAT_FIXPOINT);
	//Self-concistency check
	double err=cntr::distance_norm2(tstp,gtemp,G);
	cout << "t= " << tstp << " iteration : " << iter << "  |  Error = " << err << endl;
	// write err to output:
	{
	  hid_t sub_group_id = create_group(file_id, "it"+std::to_string(iter));
	  store_double_attribute_to_hid(sub_group_id, "err", err);
	  close_group(sub_group_id);
	}
	if(err<BootstrapMaxErr && iter>3){
	  itermax=iter;
	  break;
	}
      }
    }
    store_int_attribute_to_hid(file_id, "itermax", itermax);
    
    // BOOTSTRAPPING ------------------------------------    
    // "real time bootstrapping": solution on steps 0...SolveOrder
    {
      tstp=SolveOrder;
      GREEN_TSTP gtemp(tstp,ntau,1); // to store last timestep
      cntr::set_t0_from_mat(G);
      for(iter=0;iter<=BootstrapMaxIter;iter++){
	gtemp.set_timestep(tstp,G); //store last iteration
	if(ipt_flag){
	  cntr::dyson_start(Gweiss,mu,eloc,G,beta,h,SolveOrder);
	}else{
	  for(int n=0;n<=SolveOrder;n++) Gweiss.set_timestep(n,G);
	}	    	
	for(int n=0;n<=SolveOrder;n++){
	  get_Sigma_timestep(n,Sigma,Gweiss,Uvec);
	  get_K_timestep(n,K,Sigma,G);
	}
	cntr::dyson_start(G,mu,eloc,K,beta,h,SolveOrder);
	//Self-concistency check
	double err=cntr::distance_norm2(tstp,gtemp,G);
	cout << "t= " << tstp << " iteration : " << iter << "  |  Error = " << err << endl;
	if(err<BootstrapMaxErr && iter>3) break;
      }
    }
    
    // TIMESTEPPING ------------------------------------
    
    // timestepping using full KBE
    for(tstp=SolveOrder+1;tstp<=nt;tstp++){
      GREEN_TSTP gtemp(tstp,ntau,1); // to store last timestep
      cntr::extrapolate_timestep(tstp-1,G,SolveOrder); // extrapolate to timestep tstp
      for(iter=0;iter<=CorrectorSteps;iter++){
        G.get_timestep(tstp,gtemp); //store last iteration
	if(ipt_flag){
	  cntr::dyson_timestep(tstp,Gweiss,mu,eloc,G,beta,h,SolveOrder);
	}else{
	  Gweiss.set_timestep(tstp,G);
	}	    		
        get_Sigma_timestep(tstp,Sigma,Gweiss,Uvec);
        get_K_timestep(tstp,K,Sigma,G);
        cntr::dyson_timestep(tstp,G,mu,eloc,K,beta,h,SolveOrder);
        //Self-concistency check
        double err=cntr::distance_norm2(tstp,gtemp,G);
        cout << "t= " << tstp << " iteration : " << iter << "  |  Error = " << err << endl;
      }
    }
        
    // OUTPUT ------------------------------------

    store_int_attribute_to_hid(file_id, std::string("nt"), nt);
    store_double_attribute_to_hid(file_id, std::string("U"), U);
    store_double_attribute_to_hid(file_id, std::string("h"), h);
    store_double_attribute_to_hid(file_id, std::string("beta"), beta);
    G.write_to_hdf5(file_id, "G");
    Sigma.write_to_hdf5(file_id, "S");
    close_hdf5_file(file_id);    
    return 0;
}
