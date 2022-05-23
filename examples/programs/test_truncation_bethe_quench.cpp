#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <chrono>

// contour library headers
#include "cntr/cntr.hpp"
#include "cntr/utils/read_inputfile.hpp"

// local headers to include
#include "formats.hpp"
#include "test_truncation_utils_decl.hpp"
#include "test_truncation_utils_impl.hpp"

using namespace std;

#define CPLX complex<double>
#define GREEN cntr::herm_matrix<double>
#define CFUNC cntr::function<double>
#define GTRUNC cntr::herm_matrix_moving<double>

//==============================================================================
//         main program
//==============================================================================
int main(int argc,char *argv[]){
  //..................................................
  //                internal
  //..................................................
  GREEN G,Sigma,Delta;
  CFUNC eps0,ufunc;
  cntr::function_moving<double> eps0_t,ufunc_t;
  GTRUNC G_t,Delta_t,Sigma_t,chi_t;
  std::vector<GTRUNC> Geps_t;
  std::vector<cntr::function_moving<double> > eps_t;
  std::vector<std::vector<double> > nk_val;
  std::vector<std::vector<double> > nk_val_t;
  //std::ofstream output_nk;
  std::vector<double> ek,wk,n_val,n_val_t, elapsed_time, elapsed_time_t;
  int tc=200,nt,size,tstp,ntau,kt,tmax,neps,itermax=100,uf_precision=2;
  double beta,h,u0,u1,mu,errmax=1e-6,max_energy,min_energy;
  //..................................................
  //                output
  //..................................................
  char* flout;
  //..................................................
  //                input
  //..................................................
  char* flin;
  
  bool prima_volta_energy=false;//??
  //cout.precision(12);
  //..................................................
  //                timer
  //..................................................
  chrono::time_point<std::chrono::steady_clock> start, end, start_tot, end_tot;
  print_line_star(60);
  cout << "   Test program: Truncation on the Bethe lattice" << endl;

  start_tot = chrono::steady_clock::now();

  cout << endl;
  cout << " reading input file ..." << endl;
  cout << endl;

  try{
    //============================================================================
    //                          (II) READ INPUT
    //============================================================================
    {
      if(argc<2) throw("COMMAND LINE ARGUMENT MISSING");

      if (argc < 3) {
	// Tell the user how to run the program
	std::cerr << " Please provide a prefix for the output files. Exiting ..." << std::endl;
	/* "Usage messages" are a conventional way of telling the user
	 * how to run a program if they enter the command incorrectly.
	 */
	return 1;
      }
      flin=argv[1];
      
      //system parameters
      find_param(flin,"__u0=",u0);
      find_param(flin,"__u1=",u1);
      find_param(flin,"__neps=",neps);
      find_param(flin,"__max_energy=",max_energy);
      find_param(flin,"__min_energy=",min_energy);
      
      // solver parameters
      find_param(flin,"__nt=",nt);
      find_param(flin,"__beta=",beta);
      find_param(flin,"__ntau=",ntau);
      find_param(flin,"__h=",h);
      find_param(flin,"__tmax=",tmax);//??
      find_param(flin,"__tc=",tc);

      //output file prefix
      flout=argv[2];
		 
    }
    int time_interval=tc;
    std::vector <double> Ekin;
    //============================================================================
    //               INITIALIZE GREEN'S FUNCTIONS
    //============================================================================
    {
    size=1;
    kt=MAX_SOLVE_ORDER;
    mu=0.0;
    n_val.resize(nt+1);
    n_val_t.resize(tmax+1);
    elapsed_time.resize(nt+1);
    elapsed_time_t.resize(tmax+1);
    Ekin.resize(tmax+1);

    G=GREEN(nt,ntau,size,-1);
    Sigma=GREEN(nt,ntau,size,-1);
    Delta=GREEN(nt,ntau,size,-1);
    eps0=CFUNC(nt,1);
    ufunc=CFUNC(nt,1);
    }
    //============================================================================
    //            SOLUTION FOR EQUILIBRIUM WITH SIGMA = 0
    //============================================================================
    {
    for (tstp=-1;tstp<=nt;tstp++) eps0[tstp]=0.0;
    for (tstp=0;tstp<=nt;tstp++) ufunc[tstp]=u1;
    ufunc[-1]=u0;
    cntr::green_equilibrium_mat_bethe(G,beta);
    }
    //============================================================================
    //             SELF-CONSISTENT FULL TIME PROPAGATION UNTIL NT
    //============================================================================
    print_line_minus(50);
    cout << "        Full Time propagation" << endl;
    print_line_minus(50);

    for(tstp=-1;tstp<=nt;tstp++){

      start = std::chrono::steady_clock::now(); //start measuring the time
      
      int nt1=(tstp>=0 && tstp<=kt ? 0 : tstp);
      int nt2=(tstp>=0 && tstp<=kt ? kt : tstp);

      //Predictor
      if(tstp>kt) cntr::extrapolate_timestep(tstp-1,G,integration::I<double>(kt));
      
      cntr::herm_matrix_timestep<double> gtemp(tstp,ntau,size); // to store last timestep

      //Corrector
      for(int iter=0;iter<=itermax;iter++){
	G.get_timestep(tstp,gtemp); //store last iteration
	for(int n=nt1;n<=tstp;n++){
	  //update self-energy
	  get_sigma(n,Sigma,ufunc,G);
	  //update hybridization
	  Delta.set_timestep(n,G);
	  Delta.incr_timestep(n,Sigma);
	}
	
	if(tstp==-1){
	  // Solve the equilibrium problem
	  cntr::dyson_mat(G,Delta,mu,eps0,integration::I<double>(kt),beta,3,6);
	}else if(tstp==0){
	  cntr::set_t0_from_mat(G);
	}else if(tstp<=kt){
	  // Bootstrapping phase
	  cntr::dyson_start(G,mu,eps0,Delta,integration::I<double>(tstp),beta,h);
	}else{
	  // Solve Dyson
	  cntr::dyson_timestep(tstp,G,mu,eps0,Delta, integration::I<double>(kt),beta,h);
	}
	//Self-concistency check
	double err=cntr::distance_norm2(tstp,gtemp,G);
	cout << "t= " << tstp << " iteration : " << iter << "  |  Error = " << err << endl;

	if(err<errmax && iter>3) break;
      }
      //store the current density
      if (tstp==-1)  n_val[tstp+1]=-1.0*G.matptr(ntau)[0].real();	 
      if (tstp>=1){
	  n_val[tstp]=G.density_matrix(tstp).real();
	  //cout <<" n= "<<n_val[tstp];
	  if (tstp<=tc) n_val_t[tstp]=n_val[tstp];            	
      }

      end = std::chrono::steady_clock::now(); //end time here
      std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "[µs]" << std::endl;
      if (tstp>=1){
	elapsed_time[tstp]=  std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() ;
	if (tstp<=tc) elapsed_time_t[tstp]=elapsed_time[tstp];            	
      }


      print_line_minus(50);
    }// end loop over tstp, where tstp \in [-1,nt]
    //============================================================================
    //             OUTPUT OF NONTRUNCATED CALCULATION
    //============================================================================

    string folder_out=std::string(flout);

    std::stringstream stream_u0;
    stream_u0 << std::fixed << std::setprecision(uf_precision) << u0;
    std::string s_u0 = stream_u0.str();
    std::stringstream stream_u1;
    stream_u1 << std::fixed << std::setprecision(uf_precision) << u1;
    std::string s_u1 = stream_u1.str();


    std::stringstream stream_h;
    stream_h << std::fixed << std::setprecision(uf_precision) << h;
    std::string s_h = stream_h.str();

    string body="_tc_"+ std::to_string(tc)+"_nt_"+ std::to_string(nt)+"_h_"+s_h+ "_ui_"+s_u0+ "_uf_" +s_u1+"_tmax_"+std::to_string(tmax);
    string G_name=folder_out+"temp_G"+body+".out";
    string Sigma_name=folder_out+"temp_Sigma"+body+".out";

    G.print_to_file(G_name.c_str());
    Sigma.print_to_file(Sigma_name.c_str());


    //============================================================================
    //               INITIALIZE TRUNCATED GREEN'S FUNCTIONS
    //============================================================================
    print_line_minus(50);
    cout << "        Initializing truncated objects with tc="<<tc << endl;
    print_line_minus(50);

    G_t=GTRUNC(tc,tc,1,-1);
    Delta_t=GTRUNC(tc,tc,1,-1);
    Sigma_t=GTRUNC(tc,tc,1,-1);
    chi_t=GTRUNC(tc,tc,1,+1);
    eps0_t=cntr::function_moving<double>(tc,1);
    ufunc_t=cntr::function_moving<double>(tc,1);
    for(tstp=0;tstp<=tc;tstp++){
      eps0_t[tstp]=0.0;
      ufunc_t[tstp]=u1;
    }
    G_t.set_from_G_backward(tc,G);
    Sigma_t.set_from_G_backward(tc,Sigma);
    Delta_t.set_from_G_backward(tc,Delta);

    //============================================================================
    //             INITIALIZE THE K-DEPENDENT PROBLEM
    //============================================================================

    Geps_t.resize(neps);
    eps_t.resize(neps);

    ek.resize(neps);
    wk.resize(neps);

    calculate_weights_linear_eps(neps,wk,ek,max_energy,min_energy);
    nk_val.resize(nt+1);	
    for(tstp=0;tstp<=nt;tstp++) nk_val[tstp].resize(neps);
    nk_val_t.resize(tmax+1);
    for(tstp=0;tstp<=tmax;tstp++) nk_val_t[tstp].resize(neps);
    string G_eps_name;
    //we are here at t=nt
    for(int k=0;k<neps;k++){
      // Initaliaze temporary storage
      CPLX rho_k;
      GREEN Geps_tmp(nt,ntau,1,-1);
      CFUNC eps_tmp(nt,1);
      for(tstp=-1;tstp<=nt;tstp++) eps_tmp[tstp]=ek[k];

      // Initialize truncated objects 
      Geps_t[k]=GTRUNC(tc,tc,1,-1);
      eps_t[k]=cntr::function_moving<double>(tc,1);
      for(tstp=0;tstp<=tc;tstp++) eps_t[k][tstp]=ek[k];

      // Solve the Dyson equation
      cntr::dyson_mat(Geps_tmp,Sigma,mu,eps_tmp,integration::I<double>(kt),beta,3,6);
      cntr::dyson_start(Geps_tmp,mu,eps_tmp,Sigma,integration::I<double>(kt),beta,h);
      for(tstp=kt+1;tstp<=nt;tstp++)
	cntr::dyson_timestep(tstp,Geps_tmp,mu,eps_tmp,Sigma,
			     integration::I<double>(kt),beta,h);
      // print a selected number of k-points
      if (k%10 ==0){
	G_eps_name=folder_out+"temp_G_k_"+std::to_string(k)+body+".out";//??
	Geps_tmp.print_to_file(G_eps_name.c_str());//??
      }
      
      // Transfer the solution to the truncated object
      Geps_t[k].set_from_G_backward(tc,Geps_tmp);

      // Calculate observables and store them 
      for(int tstp=1;tstp<=nt;tstp++){
	rho_k=Geps_tmp.density_matrix(tstp);
	nk_val[tstp][k]=rho_k.real();//calculate nk at given tstp and given k
	Ekin[tstp]+=wk[k]*ek[k]*nk_val[tstp][k]; //calculate Ekin at given tstp

	if(tstp<=tc){
	  nk_val_t[tstp][k]=Geps_t[k].lesptr(tc-tstp,0)[0].imag();
	}else{
	  nk_val_t[tstp][k]=-0.1;//??
	}
      }
    } //end of loop over k 
    
    //============================================================================
    //             TRUNCATED TIME PROPAGATION UNTIL TMAX
    //============================================================================
    print_line_minus(50);
    cout << "        Truncated time propagation"<< endl;
    print_line_minus(50);

    for (tstp=tc+1;tstp<=tmax;tstp++){

      start = std::chrono::steady_clock::now(); //start measuring the time 
      double err;
      cntr::herm_matrix_timestep_moving<double> gtmp(tc,1,-1);

      //Push the truncated windows one timestep
      Sigma_t.forward();
      Delta_t.forward();
      G_t.forward();
      chi_t.forward();
      ufunc_t.forward();
      
      ufunc_t[0]=u1;//??
      // eps0=0 anyway, no need to update
      
      // propagation of local quantities
      for(int iter=0;iter<=itermax;iter++){
	//Store the last iteration
	G_t.get_timestep(0,gtmp);
	//Predictor
	if(iter==0) cntr::extrapolate_timestep(G_t,integration::I<double>(kt));
	//Update self-energy
	get_sigma(tstp,Sigma_t,ufunc_t,G_t,chi_t);
	//update hybridization
	Delta_t.set_timestep(0,G_t,0);
	Delta_t.incr_timestep(0,Sigma_t,0,1.0);
	// Solve Dyson
	cntr::dyson_timestep(G_t,Delta_t,eps0_t,mu,integration::I<double>(kt),h);
	// Self-consistency check
	err=cntr::distance_norm2(0,G_t,gtmp);
	cout << "t= " << tstp << " iteration : " << iter << "  |  Error = " << err << endl;

	if(err<errmax && iter>3) break;
      }
      //store the current density
      n_val_t[tstp]=G.lesptr(0,0)[0].imag();

      // propagation of G_eps:
      for(int k=0;k<neps;k++){
	// eps_t costant, no need to update
	Geps_t[k].forward();
	// solve k-dependent Dyson equation
	cntr::dyson_timestep(Geps_t[k],Sigma_t,eps_t[k],mu,integration::I<double>(kt),h);
	// store k-dependent occupation
	nk_val_t[tstp][k]=Geps_t[k].lesptr(0,0)[0].imag();
	// calculate the kinetic energy
	Ekin[tstp]+=wk[k]*ek[k]*nk_val_t[tstp][k]; //calculate Ekin at given tstp
      }  //loop over k

      end = std::chrono::steady_clock::now(); //end time here
      std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "[µs]" << std::endl;
      
      elapsed_time_t[tstp]=  std::chrono::duration_cast<std::chrono::microseconds>(end - start).count(); 

      //Print the stored quantities to a file
      if (tstp%time_interval==0)
	print_energy(G, Sigma, nt, kt, beta, h, folder_out, Ekin,
		     body, prima_volta_energy, tstp, time_interval,
		     tc, G_t, n_val, n_val_t);
 
      print_line_minus(50);
    } //loop over tstp from tc+1 to tmax
    
    //============================================================================
    //             OUTPUT OF TRUNCATED SOLUTION
    //============================================================================
    {
    print_n(n_val, n_val_t, elapsed_time, elapsed_time_t, u0, u1, tmax, h, tc,
	    uf_precision, nt, folder_out, body);
    print_nk(nk_val, nk_val_t, neps, ek, u0, u1, tmax, h, tc, uf_precision, nt, folder_out);
    } // end output
    
    end_tot = std::chrono::steady_clock::now();
    std::chrono::duration<double> runtime_seconds = end_tot-start_tot;
    cout << endl;
    cout << "Program execution completed"<<endl;
    cout << "Time [total] = " << runtime_seconds.count() << "s\n";
    print_line_star(60);
   
    // now calculate Geps and update.
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
//==============================================================================



