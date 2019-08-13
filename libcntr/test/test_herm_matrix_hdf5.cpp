#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>


#include "cntr.hpp"
#include "read_inputfile.hpp"
#include "hdf5_interface.hpp"
#include "hdf5_interface_cntr.hpp"


#define CPLX std::complex<double>  
using namespace std;
 

int main(int argc,char *argv[]){
	int nt,ntau,size,sig;
	double err;
	cntr::herm_matrix<double> G,G1;
	
	try{
	////////////////////////////////////////////////////////////////////////////////////	
	// (II) READ INPUT
	{
		if(argc<2) throw("COMMAND LINE ARGUMENT MISSING");
		find_param(argv[1],"__nt=",nt);
		find_param(argv[1],"__ntau=",ntau);
		find_param(argv[1],"__size=",size);
	}		
	// TEST 1
	{
		G=cntr::herm_matrix<double>(10,20,1,-1);
		G1=cntr::herm_matrix<double>(10,20,1,-1);
		cntr::green_equilibrium_bethe(G,1.0,0.1);
		write_gf_to_hdf5_file<cntr::herm_matrix<double> >("gbethe_hdf5.out","G",G);
		//G1=read_gf_from_hdf5_file<cntr::herm_matrix<double> >("gbethe_hdf5.out","G");
		G1=read_gf_from_hdf5_file<cntr::herm_matrix<double> >("gbethe_hdf5.out","G");
		for(int tstp=-1;tstp<=10;tstp++){
			err=cntr::distance_norm2(tstp,G,G1);
			cout << "timestep: " << tstp << " |G-G1|= " << err << endl;
		}
	}
	// TEST Gtav-trel
	{
		write_gf_tavtrel_to_hdf5_file("gtavtrel.out","G",G,1);
	}
	// TEST 2
	{
		double eps=2.5648;
		G=cntr::herm_matrix<double>(nt,ntau,size,-1);
		cntr::green_from_eps(G,0.0,&eps,10.0,0.1);

		write_gf_to_hdf5_file<cntr::herm_matrix<double> >("ghdf5.out","G",G);
		G.print_to_file("g.out");
	}
	// TEST3
	
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




