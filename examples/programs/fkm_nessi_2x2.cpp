#include <sys/stat.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cmath>
#include <cstring>
#include <chrono>

// contour library headers
#include <Eigen/Eigen>
#include "cntr/cntr.hpp"
#include "cntr/hdf5/hdf5_interface.hpp"
#include "cntr/hdf5/hdf5_interface_cntr.hpp"

#include "cntr/utils/read_inputfile.hpp"

using namespace std::chrono;
using ZMatrixMap = Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;

/*///////////////////////////////////////////////////////////////////////////////////////


This is a minimal test program to solve the self-consistent equation 

G0 = [ii*d/dt + mu - Delta] ^{-1}
G1 = [ii*d/dt + mu - U(t) - Delta] ^{-1}
G = (G0+G1)/2
Delta = G

U(t) = U0 for t<=(kt+1)*dt, U1 for t>0

mu=U/2

///////////////////////////////////////////////////////////////////////////////////////*/ 

void force_nessi_tti(cntr::herm_matrix<double> &G, int k) {
  int size = G.size1();

  for( int t = 0; t <= k; t++) {
    for( int i = 1; i <= k-t; i++) {
      ZMatrixMap(G.retptr(t+i,i),size,size) = ZMatrixMap(G.retptr(t,0), size, size);
      ZMatrixMap(G.lesptr(i,t+i),size,size) = ZMatrixMap(G.lesptr(0,t), size, size);
    }
  }
}


//==============================================================================
//         main program
//==============================================================================
int main(int argc,char *argv[]){
  double beta,mu,h,svdtol;
  double lambda,epsdlr; //dlr convergence parameters
  int nt,ntau,ntauN,nlvl,xi,r;
  int size=2;
  int SolverOrder=5;
  int matsMaxIter,BootstrapMaxIter,correctorSteps;
  double matsMaxErr,timeMaxErr,errC,err,err1R,err2,err1TV,err1L,errM;
  char* flin;
  std::vector<double> U,V;
  cntr::function<double> h0,h1;
  
  flin=argv[1];
  // Parameters
  
  find_param(flin,"__beta=",beta);
  find_param(flin,"__mu=",mu);
  
  // Convergence parameters
  find_param(flin,"__Ntau=",ntau);
  find_param(flin,"__Nt=",nt);
  find_param(flin,"__h=",h);

  find_param(flin,"__MatsMaxIter=",matsMaxIter);
  find_param(flin,"__BootstrapMaxIter=",BootstrapMaxIter);
  find_param(flin,"__CorrectorSteps=",correctorSteps);
  find_param(flin,"__MatsMaxErr=",matsMaxErr);
  find_param(flin,"__TimeMaxErr=",timeMaxErr);
  find_param_tvector(flin, "__U=", U, nt);
  find_param_tvector(flin, "__V=", V, nt);
  
  r=ntau;

  // Set Hamiltonian - currently fixed mu,U
  h0=cntr::function<double>(nt,size); // -> - mu - U/2.0 
  h1=cntr::function<double>(nt,size); // -> - mu + U/2.0
  cdmatrix tmp(size,size),tmp1(size,size);

  for(int tstp=-1;tstp<=nt;tstp++){
    tmp.setZero();
    tmp1.setZero();
    // h0
    tmp(0,0)=-mu-U[tstp+1]/2.0;
    tmp(1,1)=-mu-U[tstp+1]/2.0;
    h0.set_value(tstp,tmp);
    // h1
    tmp1(0,0)=-mu+U[tstp+1]/2.0;
    tmp1(1,1)=-mu+U[tstp+1]/2.0;
    tmp1(0,1)=V[tstp+1];
    tmp1(1,0)=V[tstp+1];
    h1.set_value(tstp,tmp);
  }

  // Set Nessi propagators
  cntr::herm_matrix<double> G0(nt,ntau,size,-1);
  cntr::herm_matrix<double> G1(nt,ntau,size,-1);
  cntr::herm_matrix<double> G(nt,ntau,size,-1);
  cntr::herm_matrix<double> Delta(nt,ntau,size,-1);

  

  // Matsubara
  errM=0.0;
  for(int iter=0;iter<=matsMaxIter;iter++){
    // HODLR
    errC=0.0;err=0;

    // NESSI
    cntr::dyson_mat(G0,Delta,mu,h0,integration::I<double>(SolverOrder), beta); 
    cntr::dyson_mat(G1,Delta,mu,h1,integration::I<double>(SolverOrder), beta); 
    G.set_timestep(-1,G0); 
    G.smul(-1,0.5);
    G.incr_timestep(-1,G1,0.5);
    err += cntr::distance_norm2(-1,G,Delta);
    Delta.set_timestep(-1,G);

    // Diff for hodlr
    std::cout << "Mat iter: " << iter << " errC: " << errC << " errN: " << err << std::endl;

    if(err<matsMaxErr && iter>2){
      break;
    }
  }

  err1R=0.0;err1TV=0.0;err1L=0.0; // Measure difference between Nessi and Hodlr results
  // Bootstrap
  for(int iter=0;iter<=BootstrapMaxIter;iter++){
    err=0.0;errC=0.0;

    // NESSI
    cntr::dyson_start(G0,mu,h0,Delta,integration::I<double>(SolverOrder),beta,h);
    cntr::dyson_start(G1,mu,h1,Delta,integration::I<double>(SolverOrder),beta,h);
    force_nessi_tti(G0,SolverOrder);
    force_nessi_tti(G1,SolverOrder);
    for(int tstp=0;tstp<=SolverOrder;tstp++){
      G.set_timestep(tstp,G0);
      G.smul(tstp,0.5);
      G.incr_timestep(tstp,G1,0.5);
    }

    // Diff for NESSI
    for(int tstp=0;tstp<=SolverOrder;tstp++) err += cntr::distance_norm2(tstp,G,Delta);
    
    // self-consistency 
    for(int tstp=0;tstp<=SolverOrder;tstp++) Delta.set_timestep(tstp,G);
    G.get_les(0,0,tmp);

    std::cout << "Bootstrap iter: " << iter << " errC: " << errC << " errN: " << err  << std::endl;
    if(err<timeMaxErr){
      break;
    }

  }

  auto startT = high_resolution_clock::now();
  // Timestepping
  for(int tstp = SolverOrder+1; tstp < nt; tstp++) {

    for(int iter=0;iter<=correctorSteps;iter++){
      errC=0.0;err=0.0;
      // Solve dyson NESSI
      cntr::dyson_timestep(tstp,G0,mu,h0,Delta,integration::I<double>(SolverOrder),beta,h);
      cntr::dyson_timestep(tstp,G1,mu,h1,Delta,integration::I<double>(SolverOrder),beta,h);
      G.set_timestep(tstp,G0);
      G.smul(tstp,0.5);
      G.incr_timestep(tstp,G1,0.5);
      
      // Diff for NESSI
      err += cntr::distance_norm2(tstp,G,Delta);
      Delta.set_timestep(tstp,G);

      std::cout << "Step " <<tstp <<" iter: " << iter << " errC: " << errC << " err: " <<err   <<  std::endl;
      if(err<timeMaxErr){
        std::cout << "ERRL " << tstp << " " << err1R << " " << err1L << " " << err1TV << std::endl;
        break;
      }
    }
    std::cout << std::endl;
  }

  // svdtol = 1e-5;

  // hodlr::herm_matrix_hodlr G_hodlr(G,nlvl,svdtol,1,1);
  // G_hodlr.print_memory_usage();

  hid_t file_id = open_hdf5_file("G.h5");
  hid_t group_id = create_group(file_id, "G");
  store_herm_greens_function(group_id, G);
  close_group(group_id);
  close_hdf5_file(file_id);

  auto endT = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(endT - startT);
  std::cout << "Whole timesteping time: " <<  duration.count()/1000000.0 << std::endl;
  std::cout << std::endl;

  return 0;
}
