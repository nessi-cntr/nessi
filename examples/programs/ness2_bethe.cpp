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


void save_iteration(int out_every,int iter,hid_t file_id, herm_matrix_ness &G, herm_matrix_ness &S,double err){
  // write iteration to an open h5 file
  hid_t sub_group_id = create_group(file_id, "it"+std::to_string(iter));
  if(out_every){
    G.write_to_hdf5(sub_group_id, "G");
    S.write_to_hdf5(sub_group_id, "S");
  }
  store_double_attribute_to_hid(sub_group_id, "err", err);
  close_group(sub_group_id);
}

// MAIN #######################################################################

int main(int argc, char** argv){
    
    // PARAMETERS ---------------------------------------------------------------
    // General params:
    int size = 1;   // orbital size, here 1
    int Nft;        // Fourier domain size
    double h;       // timestep 
    double beta;    //inverse temperature
    double U, mu;
    double mix, err,errmax;
    int N_it, out_every,ipt_flag;
    char *flin,*flout;

    //  read params from input file
    {
        cout << endl;
        cout << " reading input file ..." << endl;
        cout << endl;
        if (argc < 3) {
            std::cerr << "Usuage: ness_bethe.x input-file out-file" << std::endl;
            throw("COMMAND LINE ARGUMENT MISSING");
        }        
        flin=argv[1];
        flout=argv[2];
        find_param(flin, "__h=", h);
        find_param(flin, "__Nft=", Nft);
        find_param(flin, "__beta=", beta);
        find_param(flin, "__U=", U);
        find_param(flin, "__mu=", mu);
        find_param(flin, "__Niter=", N_it);
        find_param(flin, "__errmax=", errmax);
        find_param(flin, "__mixing=", mix);
        find_param(flin, "__ipt_flag=", ipt_flag);
        find_param(flin, "__out_every=", out_every);
    }
    
    // INITIALIZE GREENS FUNCTIONS ------------------------------------

    herm_matrix_ness  G_ness(Nft,size);
    herm_matrix_ness  G_old(Nft,size);
    herm_matrix_ness  SU_ness(Nft,size);
    herm_matrix_ness  W_ness(Nft,size);
    herm_matrix_ness  K_ness(Nft,size);
    herm_matrix_ness  Gweiss_ness(Nft,size);
        
    
    // INITIALIZE HAMILTONIAN ------------------------------------
    cdmatrix ham=cdmatrix::Zero(size, size);
    
    
    // INITIALIZE INITIAL GUESS (noninteracting Green's function) --
    ness2::bethedos dos(-2,2); // create bethe dos object
    green_equilibrium_ness(FERMION,G_ness, dos, beta,mu, h, FFT_TRAPEZ);    

    // start the output file:
    hid_t file_id = open_hdf5_file(flout);
    store_int_attribute_to_hid(file_id, "Nft", Nft);
    store_int_attribute_to_hid(file_id, "ipt_flag", ipt_flag);
    store_double_attribute_to_hid(file_id, "U", U);
    store_double_attribute_to_hid(file_id, "mu", mu);
    store_double_attribute_to_hid(file_id, "beta", beta);
    store_double_attribute_to_hid(file_id, "h", h);
    save_iteration(out_every,0,file_id,G_ness,SU_ness,0.0);      
           
    // DMFT LOOP ------------------------------------
    int itermax=N_it-1;
    for(int iter = 1; iter <= N_it; iter++){
      std::cout << "Iteration " + std::to_string(iter) + ":" << std::endl;
      // SAVE OLD ITERATION ----------
      G_old=G_ness;
      
      // GET SIGMA -------------
      // Self-consistent PT:
      if(ipt_flag){
	dyson(Gweiss_ness, mu, ham, G_ness, h, FFT_TRAPEZ);
      }else{
	Gweiss_ness=G_ness;
      }
      // SOLVE DYSON EQ WITH KERNEL K (no eta regularization)
      

      Bubble1_ness(W_ness, Gweiss_ness, Gweiss_ness);  // W = ii * G * G
      W_ness.smul(U*U,fft_domain::time);  // W ==> W*U**2
      Bubble2_ness(SU_ness, Gweiss_ness, W_ness);  // SU = -ii * G * W
      SU_ness.smul(-1.0, fft_domain::time);

      
      // CALCULATE KERNEL K = Sigma + G
      K_ness.set_zero(fft_domain::time);
      K_ness.incr(SU_ness,1.0,fft_domain::time);
      K_ness.incr(G_ness,1.0,fft_domain::time);
            
      dyson(G_ness, mu, ham, K_ness, h, FFT_TRAPEZ);
      err = distance_norm2(G_ness,G_old,fft_domain::time);
      // linear mixing for convergence
      G_ness.smul(mix,fft_domain::time);
      G_ness.incr(G_old,1.0-mix,fft_domain::time);
      
      // print out convergence error to previous iteration
      std::cout << "iter: " << iter << ",  GF error = " + std::to_string(err) << std::endl;      
      save_iteration(out_every,iter,file_id,G_ness,SU_ness,err);
      if(err<errmax && iter>3){
	itermax=iter;
	std::cout << "converged after " << itermax << "iterations " << std::endl;      
	break;
      }
    }
    
    // write output, postprocessing:

    // === Observables ===
    cdmatrix SUG,rho;
    // ii*[SU*G]^<(t=0):
    ness2::convolution_density_matrix(SUG,FERMION,SU_ness,G_ness,h); 
    ness2::density_matrix(rho,FERMION,G_ness); //   = ii*G^<(t=0)
    double Eint = (SUG).trace().real(); // interaction energy
    double dens = 2*rho.trace().real(); // <n_up + n_do >
    store_double_attribute_to_hid(file_id, "dens", dens);
    store_double_attribute_to_hid(file_id, "Eint", Eint);
    store_int_attribute_to_hid(file_id, "itermax", itermax);
    
    // === Fourier transform for spectra ===
    G_ness.integral_transform_to_freq(h, FFT_TRAPEZ);
    SU_ness.integral_transform_to_freq(h, FFT_TRAPEZ);    
    G_ness.write_to_hdf5(file_id, "G");
    SU_ness.write_to_hdf5(file_id, "S");
    Gweiss_ness.write_to_hdf5(file_id, "Gweiss");
    close_hdf5_file(file_id);
    
    std::cout << "Simulation complete. Data written to " << flout << "\n";
    return 0;
    
}

