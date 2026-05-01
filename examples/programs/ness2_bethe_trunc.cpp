#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <cmath>
#include <sstream>


// INCLUDES ------------------------------------------------------------

#include "cntr/cntr.hpp"
#include "cntr/ness2.hpp"
#include "cntr/utils/read_inputfile.hpp"


// NAMESPACES ------------------------------------------------------------

using namespace std;
using namespace cntr;
using namespace ness2;



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


void write_slice(hid_t file_id,int tstp,GTRUNC &G){
  int tc=G.tc();
  int Nft=4;
  while(Nft/2-1<tc) Nft*=2;
  herm_matrix_ness Gness(Nft,1);
  cntr2ness(Gness,G);
  hid_t sub_group_id = create_group(file_id, "t"+std::to_string(tstp));
  Gness.write_to_hdf5(sub_group_id, "G");
  close_group(sub_group_id);
}

int main(int argc, char** argv){
    
  int size = 1;  
  int SolveOrder = MAX_SOLVE_ORDER;
  int CorrectorSteps;
  int tmax,iter,tstp,tc,downsampling,out_every;
  char *ness_filename,*out_filename;
  if (argc < 4) {
    std::cerr << "Usage: ness2_bethe_trunc.x param-file ness_input out-file" << std::endl;
    throw("COMMAND LINE ARGUMENT MISSING");
  }
  find_param(argv[1], "__downsampling=", downsampling);
  find_param(argv[1], "__tc=", tc);
  find_param(argv[1], "__tmax=", tmax);
  find_param(argv[1], "__out_every=", out_every);	
  find_param(argv[1], "__CorrectorSteps=",CorrectorSteps);
  if(tc<=SolveOrder){
    cerr << "tc = " << tc << " <= SolveOrder (=5) not allowed" << endl;
    throw("Invalid Input");
  }
  ness_filename=argv[2];
  out_filename=argv[3];
  // read further commmon params from ness file:
  std::cout << "trying to read " << ness_filename << std::endl;
  hid_t file_id = read_hdf5_file(ness_filename);
  double h_ness=read_primitive_type<double>(file_id, "h");
  double U=read_primitive_type<double>(file_id, "U");
  double mu=read_primitive_type<double>(file_id, "mu");
  double beta=read_primitive_type<double>(file_id, "beta");
  int ipt_flag=read_primitive_type<int>(file_id, "ipt_flag");
  int nft_ness=read_primitive_type<int>(file_id, "Nft");
  double h_cntr=h_ness*downsampling;
  
  GTRUNC Gcntr(tc,1,FERMION),Scntr(tc,1,FERMION),Kcntr(tc,1,FERMION),Gweisscntr(tc,1,FERMION);
  CTRUNC U_t(tc,1),eloc_t(tc,1);
  // read G's from file
  herm_matrix_ness Gness,Gness1;
  Gness.read_from_hdf5(ness_filename,"G");
  Gness1=downsample(Gness,downsampling);
  ness2cntr(Gcntr, Gness1);
  Gness.read_from_hdf5(ness_filename,"S");
  Gness1=downsample(Gness,downsampling);
  ness2cntr(Scntr, Gness1);
  Gness.read_from_hdf5(ness_filename,"Gweiss");
  Gness1=downsample(Gness,downsampling);
  ness2cntr(Gweisscntr, Gness1);
  // init U, eloc, K=G+Sigma on all steps
  Kcntr=Gcntr;
  for(int t1=0;t1<=tc;t1++){
    cdmatrix mtmp;
    mtmp = cdmatrix::Constant(1, 1, U);
    U_t.set_value(t1,mtmp);
    mtmp = cdmatrix::Zero(1, 1);
    eloc_t.set_value(t1,mtmp);
    Kcntr.incr_timestep(t1,Scntr,t1,1.0);
  }
  // start output file:
  std::cout << "trying to write " << out_filename << std::endl;					   
  hid_t out_file_id = open_hdf5_file(out_filename);
  store_int_attribute_to_hid(out_file_id, "tc", tc);
  store_double_attribute_to_hid(out_file_id, "U", U);
  store_double_attribute_to_hid(out_file_id, "mu", mu);
  store_double_attribute_to_hid(out_file_id, "beta", beta);
  write_slice(out_file_id,0,Gcntr);
  // truncated time-stepping (tstp is relative to tc)
  for(tstp=1;tstp<=tmax;tstp++){
    cdmatrix mtmp;
    double err;
    GTRUNC_TSTP gtmp(tc,1,FERMION);      
    Scntr.forward();
    Kcntr.forward();
    Gcntr.forward();
    Gweisscntr.forward();
    cntr::extrapolate_timestep(Gcntr,SolveOrder);      
    for(iter=0;iter<=CorrectorSteps;iter++){
      gtmp.set_timestep(Gcntr,0);
      if(ipt_flag){
	cntr::dyson_timestep(Gweisscntr,Gcntr,eloc_t,mu,SolveOrder,h_cntr);
      }else{
	Gweisscntr.set_timestep(0,Gcntr,0);
      }	
      //Update self-energy
      get_Sigma_timestep(Scntr,Gweisscntr,U_t);
      get_K_timestep(Kcntr,Scntr,Gcntr);
      cntr::dyson_timestep(Gcntr,Kcntr,eloc_t,mu,SolveOrder,h_cntr);
      // Self-consistency check
      err=cntr::distance_norm2(Gcntr,0,gtmp);	
      cout << "tstp= " << tstp << " iteration : " << iter << "  |  Error = " << err << endl;	
    }
    if(tstp%out_every==0) write_slice(out_file_id,tstp,Gcntr);
  }
  close_hdf5_file(file_id);
  return 0;
}
