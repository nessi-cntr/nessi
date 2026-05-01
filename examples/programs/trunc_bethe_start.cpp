#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

// contour library headers
#include "cntr/cntr.hpp"
#include "cntr/utils/read_inputfile.hpp"

// local headers to include
#include "formats.hpp"

using namespace std;


//------------------------------------------------------------------------------
// compute Kernel K = G + Sigma, with a second order self-energy
// Sigma(t,t')=U(t)G(t,t')G(t,t')G(t',t)U(t')

// for full Green's functions
void get_Sigma_timestep(int tstp, GREEN &Sigma, GREEN &G, CFUNC &U){
  GREEN_TSTP W(tstp,Sigma.ntau(),Sigma.size1(),BOSON); // temporary variable;
  // Note: a bubble of two Fermion GF is bosonic
  cntr::Bubble1(tstp,W,G,G); // W(t,t')  set to chi=ii*G(t,t')G(t',t);
  W.left_multiply(tstp,U);
  W.right_multiply(tstp,U); // now W  set to W(t,t')=U(t)chi(t,t')U(t')
  cntr::Bubble2(tstp,Sigma,G,W); // Sigma(t,t') is set to Sigma=ii*G(t,t')W(t',t);
  Sigma.smul(tstp,-1.0); // the final -1 sign.
}

void get_K_timestep(int tstp, GREEN &K, GREEN &Sigma, GREEN &G){
  K.set_timestep(tstp,Sigma); // K=Sigma
  K.incr_timestep(tstp,G,1.0);  // K += G
}

//==============================================================================
//         main program
//==============================================================================
int main(int argc,char *argv[]){    
  int MatsMaxIter,CorrectorSteps,BootstrapMaxIter;
  double MatsMaxErr,BootstrapMaxErr;
  int SolveOrder = MAX_SOLVE_ORDER;
  int nt,ntau;
  double beta,h,U1,mu=0;
  char *flin,*flout;
  int tmax,iter,tstp;
  int nk; // number of kpoints
  std::vector<double> ek;
  
  try{
    {
      cout << endl;
      cout << " reading input file ..." << endl;
      cout << endl;
      if (argc < 3) {
	std::cerr << "Usuage: test_truncation_bethe_quench.x param-file out-file  epsk_1 epsk_2 .... epsk_n " << std::endl;
	throw("COMMAND LINE ARGUMENT MISSING");	  
      }
      flin=argv[1];
      flout=argv[2];
      // solver parameters
      find_param(flin,"__beta=",beta);
      find_param(flin,"__h=",h);	
      find_param(flin,"__U1=",U1);
      find_param(flin,"__nt=",nt);
      find_param(flin,"__ntau=",ntau);
      find_param(flin,"__BootstrapMaxIter=",BootstrapMaxIter);
      find_param(flin,"__BootstrapMaxErr=",BootstrapMaxErr);      
      find_param(flin,"__CorrectorSteps=",CorrectorSteps);	
      if(nt<=SolveOrder){
	cerr << "nt = " << nt << " <= SolveOrder (=5) not allowed" << endl;
	throw("Invalid Input");
      }
      // read in epsk_values from command line
      nk= argc-3;
      ek.resize(nk);
      for(int k=0;k<nk;k++) stringstream(argv[3+k]) >> ek[k];
    }
    // start from a full KBE simulations for steps [-1,...,tc]:
    GREEN G(nt,ntau,1,FERMION);
    GREEN Sigma(nt,ntau,1,FERMION);
    GREEN K(nt,ntau,1,FERMION);
    CFUNC U(nt,1);
    CFUNC eloc(nt,1); // =0
    
    std::vector<GREEN> Gk(nk);
    std::vector<CFUNC> epsk_func(nk);
    for(int k=0;k<nk;k++){
      epsk_func[k]=CFUNC(nt,1);
      Gk[k]=GREEN(nt,ntau,1,FERMION);
    }
    // init function:
    // U=0 on matsubara by initialization, set real time values
    for(tstp=0;tstp<=nt;tstp++){
      cdmatrix ut = cdmatrix::Constant(1, 1, U1);
      U.set_value(tstp,ut);
    }	  
    // init G on matsubara branch with noninteracting GF:      
    cntr::green_equilibrium_mat_bethe(G,beta);
    Sigma.set_timestep_zero(-1);
    get_K_timestep(-1,K,Sigma,G);
    // "real time bootstrapping": solution on steps 0...SolveOrder
    {
      tstp=SolveOrder;
      GREEN_TSTP gtemp(tstp,ntau,1); // to store last timestep
      cntr::set_t0_from_mat(G);
      for(iter=0;iter<=BootstrapMaxIter;iter++){	  
	gtemp.set_timestep(tstp,G); //store last iteration
	for(int n=0;n<=SolveOrder;n++){
	  get_Sigma_timestep(n,Sigma,G,U);
	  get_K_timestep(n,K,Sigma,G);
	}
	cntr::dyson_start(G,mu,eloc,K,beta,h,SolveOrder);
	//Self-concistency check
	double err=cntr::distance_norm2(tstp,gtemp,G);
	cout << "t= " << tstp << " iteration : " << iter << "  |  Error = " << err << endl;	
	if(err<BootstrapMaxErr && iter>3) break;
      }
      if(iter==BootstrapMaxIter){
	throw("Bootstrapping not converged.");
      }
    }    
    // timestepping using full KBE
    for(tstp=SolveOrder+1;tstp<=nt;tstp++){
      GREEN_TSTP gtemp(tstp,ntau,1); // to store last timestep
      cntr::extrapolate_timestep(tstp-1,G,SolveOrder); // extrapolate to timestep tstp
      for(iter=0;iter<=CorrectorSteps;iter++){
	G.get_timestep(tstp,gtemp); //store last iteration
	get_Sigma_timestep(tstp,Sigma,G,U);
	get_K_timestep(tstp,K,Sigma,G);
	cntr::dyson_timestep(tstp,G,mu,eloc,K,beta,h,SolveOrder);
	//Self-concistency check
	double err=cntr::distance_norm2(tstp,gtemp,G);
	cout << "t= " << tstp << " iteration : " << iter << "  |  Error = " << err << endl;	
      }
    }
    // solve dyson equastion Gk = [idt + mu - epsk - Sigma ]^{-1} for some k values:	    
    for (int k=0;k<nk;k++){	  
      cdmatrix ekval = cdmatrix::Constant(1, 1, ek[k]);
      for (int tstp=-1;tstp<=nt;tstp++){
	epsk_func[k].set_value(tstp,ekval);
      }
      cntr::dyson(Gk[k],mu,epsk_func[k],Sigma,beta,h,SolveOrder,CNTR_MAT_FIXPOINT);
    }
    {
      // output:
      hid_t file_id = open_hdf5_file(flout);
      store_int_attribute_to_hid(file_id, std::string("nt"), nt);
      store_double_attribute_to_hid(file_id, std::string("U"), U1);
      store_double_attribute_to_hid(file_id, std::string("h"), h);      
      store_double_attribute_to_hid(file_id, std::string("beta"), beta);      
      G.write_to_hdf5(file_id, "G"); 
      Sigma.write_to_hdf5(file_id, "Sigma");      
      store_int_attribute_to_hid(file_id, std::string("nk"), nk);
      hid_t group_id = create_group(file_id, "lattice_functions");
      for(int k=0;k<nk;k++){
	std::string groupname=to_string(k);
        hid_t sub_group_id = create_group(group_id,groupname.c_str());
	store_double_attribute_to_hid(sub_group_id, std::string("epsk"), ek[k]);
	Gk[k].write_to_hdf5(sub_group_id,"Gk");
	close_group(sub_group_id);
      }
      close_group(group_id);
      close_hdf5_file(file_id);
    }
  } // try
  catch(char *message){
    cerr << "exception\n**** " << message << " ****" << endl;
    cerr << " No input file found. Exiting ... " << endl;
  }
  catch(...){
    cerr << " No input file found. Exiting ... " << endl;
  }
  return 0;
}
