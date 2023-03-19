#include <sys/stat.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cmath>
#include <cstring>
#include <chrono>

// contour library headers
#include "cntr/cntr.hpp"
#include "cntr/utils/read_inputfile.hpp"

#include "cntr/hdf5/hdf5_interface.hpp"
#include "cntr/hdf5/hdf5_interface_cntr.hpp"

#define CPLX std::complex<double>  
using namespace std;
 

/*///////////////////////////////////////////////////////////////////////////////////////


This is a minimal test program to solve the self-consistent equation 

G0 = [ii*d/dt + mu - Delta] ^{-1}
G1 = [ii*d/dt + mu - U(t) - Delta] ^{-1}
G = (G0+G1)/2
Delta = G

U(t) = U0 for t<=0, U1 for t>0

mu=U/2

///////////////////////////////////////////////////////////////////////////////////////*/ 


int main(int argc,char *argv[]){

	int nt=50,ntau=200,tstp,itermax=100,iter_rtime=3,kt=5,size,sig;
	double U0=1.0,U1=2.0,beta=5.0,h=0.01,mu,err,errmax=1e-7;
	cntr::herm_matrix<double> G0,G1,G,Delta;
	cntr::function<double> f1,f0;
	std::vector<double> densityOLD0(nt+2),densityOLD1(nt+2),density0(nt+2),density1(nt+2);

	
	size=1;
	sig=-1;
	G0=cntr::herm_matrix<double>(nt,ntau,size,sig);
	G1=cntr::herm_matrix<double>(nt,ntau,size,sig);
	Delta=cntr::herm_matrix<double>(nt,ntau,size,sig);
	G=cntr::herm_matrix<double>(nt,ntau,size,sig);
	// init with bethe, only for matsubara branch
	cntr::green_equilibrium_mat_bethe(Delta,beta);
	f0=cntr::function<double>(nt,size); // -> - mu 
	f1=cntr::function<double>(nt,size); // -> - mu + U 
	mu=0.0;
	f0[-1]=-U0*0.5;
	for(tstp=0;tstp<=nt;tstp++) f0[tstp]=-U1*0.5;
	f1[-1]=U0*0.5;
	for(tstp=0;tstp<=nt;tstp++) f1[tstp]=U1*0.5;		

	bool matsubara_converged=false;
	tstp=-1;
	cntr::herm_matrix<double> gtemp; // to store last timestep
	gtemp=cntr::herm_matrix<double>(kt,ntau,size,-1); // contains only kt timesteps
	for(int iter=0;iter<=itermax;iter++){
		cntr::dyson_mat(G0,Delta,mu,f0,integration::I<double>(kt), beta, CNTR_MAT_FOURIER);	
		cntr::dyson_mat(G1,Delta,mu,f1,integration::I<double>(kt), beta, CNTR_MAT_FOURIER);	
		G.set_timestep(tstp,G0); 
		G.smul(tstp,0.5);
		G.incr_timestep(tstp,G1,0.5);
		//
		err=cntr::distance_norm2(tstp,G,gtemp);
		// cout << "iteration : " << iter << " |G-Gtemp|= " << err << endl;
		if(err<errmax && iter>2){
			matsubara_converged=true;
			break;
		}
		gtemp.set_timestep(tstp,G); // save timestep -1 to gtemp
		// self-consistency 
		Delta.set_timestep(tstp,G);
	}
	if(!matsubara_converged){
		cout << "cannot cnverge matsubara; end here " << endl;
		// should end here ....
	}
	// STEP 2): first kt timesteps: do iteration like above:
	for(int iter=0;iter<=itermax;iter++){
		cntr::dyson_start(G0,mu,f0,Delta,integration::I<double>(kt),beta,h);
		cntr::dyson_start(G1,mu,f1,Delta,integration::I<double>(kt),beta,h);
		for(tstp=0;tstp<=kt;tstp++){
			G.set_timestep(tstp,G0); 
			G.smul(tstp,0.5);
			G.incr_timestep(tstp,G1,0.5);
		}
		//
		err=0.0;
		for(tstp=0;tstp<=kt;tstp++) err += cntr::distance_norm2(tstp,G,gtemp);
				// cout << "START iteration : " << iter << " |G-Gtemp|= " << err << endl;
		if(err<errmax) break;
		for(tstp=0;tstp<=kt;tstp++) gtemp.set_timestep(tstp,G); // save timestep -1 to gtemp			
				// self-consistency 
		for(tstp=0;tstp<=kt;tstp++) Delta.set_timestep(tstp,G);
	}
	// STEP 3): time-stepping kt+1 ... nt
	for(tstp=kt+1;tstp<=nt;tstp++){
		// extrapoaltion of Delta by one timestep:
		cntr::extrapolate_timestep(tstp-1,Delta,integration::I<double>(kt));
		// cout << " ... at timestep " << tstp << endl;
		for(int iter=1;iter<=iter_rtime;iter++){
			cntr::dyson_timestep(tstp,G0,mu,f0,Delta,integration::I<double>(kt),beta,h);
			cntr::dyson_timestep(tstp,G1,mu,f1,Delta,integration::I<double>(kt),beta,h);
			G.set_timestep(tstp,G0); 
			G.smul(tstp,0.5);
			G.incr_timestep(tstp,G1,0.5);
			Delta.set_timestep(tstp,G); 
		}
	}

	hid_t file_id = open_hdf5_file("G1.h5");
	hid_t group_id = create_group(file_id, "G");
    store_herm_greens_function(group_id, G);
    close_group(group_id);
    close_hdf5_file(file_id);

}
