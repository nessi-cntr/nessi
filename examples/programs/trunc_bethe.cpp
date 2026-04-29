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
using namespace cntr;


//------------------------------------------------------------------------------
// compute Kernel K = G + Sigma, with a second order self-energy
// Sigma(t,t')=U(t)G(t,t')G(t,t')G(t',t)U(t')

void get_Sigma_timestep(GTRUNC &Sigma, GTRUNC &G, CTRUNC &U){
  GTRUNC_TSTP W(Sigma.tc(),Sigma.size1(),BOSON); // temporary variable; 
  cntr::Bubble1(W,G,G); // W(t,t') is set to chi=ii*G(t,t')G(t',t);
  W.left_multiply(U);
  W.right_multiply(U); // now W is set to W(t,t')=U(t)chi(t,t')U(t') at leading timestep
  cntr::Bubble2(Sigma,G,W); // Sigma(t,t') is set to Sigma=ii*G(t,t')W(t',t);
  Sigma.smul(0,-1.0); // the final -1 sign at the leading timestep 0
}

void get_K_timestep(GTRUNC &K, GTRUNC &Sigma, GTRUNC &G){
  K.set_timestep(0,Sigma,0); // K=Sigma at leading timestep 0
  K.incr_timestep(0,G,0,1.0);  // K += G  at leading timestep 0
}

void write_slice(int tstp, hid_t file_id, GTRUNC &G_t,GTRUNC &Sigma_t){
  hid_t sub_group = create_group(file_id,"t"+std::to_string(tstp));
  G_t.write_timestep_to_hdf5(0,sub_group,"G"); // write leading timestep "0"
  Sigma_t.write_timestep_to_hdf5(0,sub_group,"S");
  close_group(sub_group); 
}


//==============================================================================
//         main program
//==============================================================================
int main(int argc,char *argv[]){    
  int CorrectorSteps;
  int SolveOrder = MAX_SOLVE_ORDER;
  int tc;
  double mu=0;
  char *flin,*fldata,*flout,*flout1;
  int tmax,iter,tstp,out_every;
  
  try{
    {
      cout << endl;
      cout << " reading input file ..." << endl;
      cout << endl;
      if (argc < 5) {
	std::cerr << "Usuage: trunc_bethe.x param-file GF_input nk_out-file hd5_out ..." << std::endl;
	throw("COMMAND LINE ARGUMENT MISSING");	  
      }
      flin=argv[1];
      fldata=argv[2];
      flout=argv[3];
      flout1=argv[4];      
      // solver parameters
      find_param(flin,"__CorrectorSteps=",CorrectorSteps);	
      find_param(flin,"__tmax=",tmax);		
      find_param(flin,"__tc=",tc);
      find_param(flin,"__out_every=",out_every);
      if(tc<=SolveOrder){
	cerr << "tc = " << tc << " <= SolveOrder (=5) not allowed" << endl;
	throw("Invalid Input");
      }	
    }
    hid_t file_id = read_hdf5_file(fldata);
    int nt=read_primitive_type<int>(file_id, "nt");
    int nk=read_primitive_type<int>(file_id, "nk");
    double U1=read_primitive_type<double>(file_id, "U");
    double h=read_primitive_type<double>(file_id, "h");
    if(tc>nt){
      cerr << "Stored GF have nt=" << nt<< "; must be larger than tc = " << tc << endl;
      throw("Invalid Input");
    }
    cout << "simulation with:\n\ttc=" << tc << endl;
    cout << "\tnt=" << nt << endl;
    cout << "\tnk=" << nk << endl;
    cout << "\tU=" << U1 << endl;
    cout << "\th=" << h << endl;    
    GTRUNC G_t(tc,1,FERMION),Sigma_t(tc,1,FERMION),K_t(tc,1,FERMION);
    CTRUNC U_t(tc,1),eloc_t(tc,1);
    // read G's from file
    GREEN Gtmp;      
    Gtmp.read_from_hdf5(file_id, "G"); //read the harm_matrix
    G_t.set_from_G_backward(Gtmp,Gtmp,tc); 
    Gtmp.read_from_hdf5(file_id, "Sigma");
    Sigma_t.set_from_G_backward(Gtmp,Gtmp,tc);    
    for(int t1=0;t1<=tc;t1++){
      cdmatrix mtmp = cdmatrix::Constant(1, 1, U1);
      U_t.set_value(t1,mtmp);
    }
    vector<GTRUNC> Gk_t(nk);
    vector<CTRUNC> epsk_t(nk);
    // init all functions (which are constant here)
    std::vector<double> ek(nk);
    hid_t group_id = open_group(file_id,"lattice_functions");
    for(int k=0;k<nk;k++){
      std::string groupname=to_string(k);
      hid_t sub_group_id = open_group(group_id,groupname.c_str());
      ek[k]=read_primitive_type<double>(sub_group_id, "epsk");
      cdmatrix mtmp = cdmatrix::Constant(1, 1, ek[k]);      
      Gtmp.read_from_hdf5(sub_group_id, "Gk"); 
      Gk_t[k]=GTRUNC(tc,1,FERMION);
      Gk_t[k].set_from_G_backward(Gtmp,Gtmp,tc);
      epsk_t[k]=CTRUNC(tc,1);
      for(int t1=0;t1<=tc;t1++){
	epsk_t[k].set_value(t1,mtmp);
      }
    }
    // init K=G+Sigma on all steps
    K_t=G_t;
    for(int t1=0;t1<=tc;t1++){
      K_t.incr_timestep(t1,Sigma_t,t1,1.0);
    }
    close_hdf5_file(file_id);
    // open file to store slices
    file_id=open_hdf5_file(flout1,h5_mode::create_truncate);
    // write first slice (time tc)
    write_slice(tc, file_id, G_t,Sigma_t);
    
    // 
    vector<vector<double> > densk(nk); // to store results
    vector<double>  dens(tmax+1); // to store results
    for(int t1=0;t1<=tc;t1++){
      cdmatrix mtmp;
      G_t.density_matrix(t1,mtmp);
      dens[tc-t1]=mtmp(0,0).real();	
    }
    for(int k=0;k<nk;k++){
      densk[k].resize(tmax+1);
      for(int t1=0;t1<=tc;t1++){
	cdmatrix mtmp;
	Gk_t[k].density_matrix(t1,mtmp);
	densk[k][tc-t1]=mtmp(0,0).real();	
      }
    }
    // truncatedtime-stepping
    for(tstp=tc+1;tstp<=tmax;tstp++){
      cdmatrix mtmp;
      double err;
      GTRUNC_TSTP gtmp(tc,1,FERMION);      
      Sigma_t.forward();
      K_t.forward();
      G_t.forward();
      // Note: not pushing functions, because they are constant anyway
      // otherwhise: U_t.forward() + set U_t.set_value(0, ....) etc.
      //Predictor
      cntr::extrapolate_timestep(G_t,SolveOrder);      
      // propagation of local quantities
      for(iter=0;iter<=CorrectorSteps;iter++){
	//Store the last iteration
	gtmp.set_timestep(G_t,0);
	//Update self-energy
	get_Sigma_timestep(Sigma_t,G_t,U_t);
	get_K_timestep(K_t,Sigma_t,G_t);
	// Solve Dyson
	cntr::dyson_timestep(G_t,K_t,eloc_t,mu,SolveOrder,h);
	// Self-consistency check
	err=cntr::distance_norm2(G_t,0,gtmp);	
	cout << "t= " << tstp << " iteration : " << iter << "  |  Error = " << err << endl;	
      }
      // solve for Geps
      for(int k=0;k<nk;k++){
	Gk_t[k].forward();
	cntr::dyson_timestep(Gk_t[k],Sigma_t,epsk_t[k],mu,SolveOrder,h);
	Gk_t[k].density_matrix(0,mtmp);
	densk[k][tstp]=mtmp(0,0).real();	
      }
      G_t.density_matrix(0,mtmp);
      dens[tstp]=mtmp(0,0).real();

      // write slices of G and Sigma to hd5 output:
      if((tstp-tc)%out_every==0){
	write_slice(tstp, file_id, G_t,Sigma_t);
      }
    }
    // write last slice (if not written alread)
    if((tmax-tc)%out_every != 0){
      write_slice(tmax, file_id, G_t,Sigma_t);
    }
    close_hdf5_file(file_id);
    
    //==============================
    //  write densities to a  text-file
    //==============================
    {
      FILE *outfile = fopen(flout,"w");
      fprintf(outfile,"\n");
      for(int k=0;k<nk;k++){
	for(int tstp=0;tstp<=tmax;tstp++){
	  fprintf(outfile,"%.10g ", ek[k]);
	  fprintf(outfile," %d ", tstp);
	  fprintf(outfile,"%.10g ", densk[k][tstp]); 
	  fprintf(outfile,"\n");
	}
	fprintf(outfile,"\n");
      }
      fclose(outfile);
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
