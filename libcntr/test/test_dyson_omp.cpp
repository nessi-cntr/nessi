#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#define CNTR_USE_OMP 1
#define CNTR_USE_MPI 1
#include "cntr.hpp"
#include "read_inputfile.hpp"

#define CPLX std::complex<double>  
using namespace std;
 
#include <time.h>
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

/*///////////////////////////////////////////////////////////////////////////////////////


This solves the self-consistent equation 

ii*d/dt  G(t,t') - \int_C dt1 \Delta(t,t1) G(t1,t') = \delta_C(t,t')

with self-consistent hybridization function  \Delta(t,t') = G(t,t')


using OMP-parallel or serial version;

OUTPUT: cpu time per timestep, comparison to serial version and exact result


///////////////////////////////////////////////////////////////////////////////////////*/ 


int main(int argc,char *argv[]){
	int nomp=omp_get_max_threads();
	int iter_rtime;
	int nt,ntau,kt,size,sig,tstp;
	double mu,beta,h;
	cntr::herm_matrix<double> G,Delta,Gomp;
	cntr::herm_matrix<double> Gbethe;
	cntr::function<double> f;
	std::vector<double> time_serial_wall,time_omp_wall;
	std::vector<double> time_serial_cpu,time_omp_cpu;
	
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
		find_param(argv[1],"__iter_rtime=",iter_rtime);
        find_param(argv[1],"__kt=",kt);
	}	
	////////////////////////////////////////////////////////////////////////////////////	
	// (IV) INITIALIZE GREEN'S FUNCTIONS; get exact result
	{
		size=1;
		sig=-1;
		G=cntr::herm_matrix<double>(nt,ntau,size,sig);
		Gomp=cntr::herm_matrix<double>(nt,ntau,size,sig);
		Delta=cntr::herm_matrix<double>(nt,ntau,size,sig);
		Gbethe=cntr::herm_matrix<double>(nt,ntau,size,sig); 
		cntr::green_equilibrium_bethe(Gbethe,beta,h);
		f=cntr::function<double>(nt,size);
		mu=0.0;
		for(tstp=-1;tstp<=nt;tstp++) f[tstp]=0.0;
		time_serial_cpu.resize(nt+1,0);
		time_omp_cpu.resize(nt+1,0);
		time_serial_wall.resize(nt+1,0);
		time_omp_wall.resize(nt+1,0);
		Delta=Gbethe; // for startig steps
		G=Gbethe;     // for startig steps
		Gomp=Gbethe;     // for startig steps
	}
	////////////////////////////////////////////////////////////////////////////////////
	// SERIAL SOLUTION FOR TIME 0... KT
	{
		for(tstp=kt+1;tstp<=nt;tstp++){
			// extrapoaltion of Delta by one timestep:
			cntr::extrapolate_timestep(tstp-1,Delta,integration::I<double>(kt));
			cout << " serial solution ... at timestep " << tstp << endl;
			double t1w=get_wall_time();
			double t1c=get_cpu_time();
			for(int iter=1;iter<=iter_rtime;iter++){
				cntr::dyson_timestep(tstp,G,mu,f,Delta,integration::I<double>(kt),beta,h);
				Delta.set_timestep(tstp,G);
			}
			double t2w=get_wall_time();
			double t2c=get_cpu_time();
			time_serial_cpu[tstp]=t2c-t1c;
			time_serial_wall[tstp]=t2w-t1w;
		}
	}
	////////////////////////////////////////////////////////////////////////////////////
	// OMP SOLUTION FOR TIME 0... KT
	{
		for(tstp=kt+1;tstp<=nt;tstp++){
			// extrapoaltion of Delta by one timestep:
			cntr::extrapolate_timestep(tstp-1,Delta,integration::I<double>(kt));
			cout << " omp solution ... at timestep " << tstp << endl;
			double t1w=get_wall_time();
			double t1c=get_cpu_time();
			for(int iter=1;iter<=iter_rtime;iter++){
				cntr::dyson_timestep_omp(nomp,tstp,Gomp,mu,f,Delta,integration::I<double>(kt),beta,h);
				Delta.set_timestep(tstp,Gomp);
			}
			double t2w=get_wall_time();
			double t2c=get_cpu_time();
			time_omp_cpu[tstp]=t2c-t1c;
			time_omp_wall[tstp]=t2w-t1w;
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////
	for(tstp=0;tstp<=nt;tstp++){
		cout << "t: " << tstp;
		cout << " |serial-exact|= " << cntr::distance_norm2(tstp,G,Gbethe);
		cout << " |omp-exact|= " << cntr::distance_norm2(tstp,Gomp,Gbethe);
		cout << " |serial-omp|= " << cntr::distance_norm2(tstp,G,Gomp);
		cout << " walltime(serial)= " << time_serial_wall[tstp];
		cout << " cputime(serial)= " << time_serial_cpu[tstp];
		cout << " walltime(serial)= " << time_omp_wall[tstp];
		cout << " cputime(omp)= " << time_omp_cpu[tstp];
		cout << endl;
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




