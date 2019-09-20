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

// local headers to include
#include "formats.hpp"
// #include "read_inputfile.hpp"
#include "gw_latt_decl.hpp"
#include "gw_kpoints_decl.hpp"
#include "gw_selfene_decl.hpp"

using namespace std;

#define CFUNC cntr::function<double>
#define GREEN cntr::herm_matrix<double>
#define DIST_TIMESTEP cntr::distributed_timestep_array<double>


// -----------------------------------------------------------------------
// double KineticEnergy(int tstp, cdmatrix &h0, GREEN &G){
  // int nst = G.size1();
  // cdmatrix rho(nst,nst);
  // G.density_matrix(tstp, rho);
  // return ((h0*rho).trace()).real();
// }


// -----------------------------------------------------------------------
// double CorrelationEnergy(int tstp,cdmatrix &h0, CFUNC &hmf, GREEN &G, GREEN &Sigma,
// 			 integration::Integrator<double> &I, double beta, double dt){
  // int nst = G.size1();
  // double Emf,Ecorr;
  // cdmatrix rho(nst,nst),vmf(nst,nst);

  // hmf.get_value(tstp,vmf);
  // vmf = vmf - h0;

  // G.density_matrix(tstp, rho);
  // Emf = 0.5*((vmf*rho).trace()).real();

  // Ecorr = cntr::correlation_energy(tstp, G, Sigma, I, beta, dt);
//   return Emf + Ecorr;
// }



//==============================================================================
//         main program
//==============================================================================
int main(int argc,char *argv[]){
  //..................................................
  //                input
  //..................................................
  int SaveGreen,SaveMomentum,output;
  int Nt,Ntau,MatsMaxIter,CorrectorSteps,SolverOrder,Nk;
  int BootstrapMaxIter;
  double HoppingT,HubbardU,V,beta,dt,MuChem,MatsMaxErr,BootstrapMaxErr,TimeMaxErr;
  double RampV0;
  char* flin;
  char* flsave;
  char* flout;
  std::vector<double> Epulse;
  //..................................................
  //                internal
  //..................................................
  int tstp;
  double err;
  int Norb=1; // Single-band case
  cdmatrix h0;
  GREEN gtemp,Wtemp;
  CFUNC hmf;
  int tid,ntasks,tid_root;
  std::vector<int> tid_map;

  //..................................................
  //                MPI
  //..................................................

  {
    MPI::Init(argc,argv);
    ntasks=MPI::COMM_WORLD.Get_size();
    tid=MPI::COMM_WORLD.Get_rank();
    tid_root=0;
  }
  //..................................................
  //                timer
  //..................................................
  chrono::time_point<std::chrono::system_clock> start, end, start_tot, end_tot;
  if(tid==tid_root){
    print_line_star(60);
    cout << "   Test program: Hubbard 1d in GW approximation" << endl;
    print_line_star(60);
    start_tot = chrono::system_clock::now();

    cout << endl;
    cout << " reading input file ..." << endl;
    cout << endl;
  }

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

      // system parameters
      find_param(flin,"__Nk=",Nk);
      find_param(flin,"__HoppingT=",HoppingT);
      find_param(flin,"__HubbardU=",HubbardU);
      find_param(flin,"__V=",V); //Nearest neighbor interaction
      find_param(flin,"__MuChem=",MuChem);
      find_param(flin,"__beta=",beta);
      

      // solver parameters
      find_param(flin,"__SaveGreen=",SaveGreen);
      find_param(flin,"__SaveMomentum=",SaveMomentum);
      find_param(flin,"__Nt=",Nt);
      find_param(flin,"__Ntau=",Ntau);
      find_param(flin,"__dt=",dt);
      find_param(flin,"__MatsMaxIter=",MatsMaxIter);
      find_param(flin,"__MatsMaxErr=",MatsMaxErr);
      find_param(flin,"__BootstrapMaxIter=",BootstrapMaxIter);
      find_param(flin,"__BootstrapMaxErr=",BootstrapMaxErr);
      find_param(flin,"__TimeMaxErr=",TimeMaxErr);
      find_param(flin,"__CorrectorSteps=",CorrectorSteps);
      find_param(flin,"__SolverOrder=",SolverOrder);

      if(SaveMomentum==1){
        find_param(flin,"__output=",output);
      }
      // Pulse 
      find_param_tvector(argv[1],"__Epulse=",Epulse,Nt);

      // output file prefix
      flout=argv[2];

    }
    //============================================================================
    //                   MEMORY REQUIREMENTS
    //============================================================================
    if(tid==tid_root){
      print_line_minus(50);
      cout << "     Memory requirements" << endl;
      print_line_minus(50);

      const size_t size_MB=1024*1024;
      size_t mem_1time=0, mem_2time=0;

      mem_1time += cntr::mem_function<double>(Nt,Norb)*Nk; // Ut
      mem_1time += cntr::mem_function<double>(Nt,Norb)*Nk; // hmf
      mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Norb)*Nk; // G
      mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Norb)*Nk; // G0
      mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Norb)*Nk; // G0Sigma
      mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Norb)*Nk; // SigmaG0
      mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Norb)*Nk; // Sigma
      mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Norb)*Nk; // PV
      mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Norb)*Nk; // VP
      mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Norb)*Nk; // W
      // convert to MB
      mem_1time = ceil(mem_1time/(double)size_MB);
      mem_2time = ceil(mem_2time/(double)size_MB);

      cout << "Total" << endl;
      cout << "Hamiltonian : " << mem_1time << " MB" << endl;
      cout << "Propagators : " << mem_2time << " MB" << endl;

      cout << "Per rank" << endl;
      cout << "Hamiltonian : " << mem_1time/double(ntasks) << " MB" << endl;
      cout << "Propagators : " << mem_2time/double(ntasks) << " MB" << endl;

      print_line_minus(50);
      cout << "\n\n";
    }
    //============================================================================
    //               (IV) INITIALIZE lattice_1d_1b AND CORRESPONDING PROPAGATORS
    //============================================================================

    CFUNC Ut(Nt,Norb),Vt(Nt,Norb),tt(Nt,Norb);
    Ut.set_constant(HubbardU*MatrixXcd::Identity(Norb,Norb));
    Vt.set_constant(V*MatrixXcd::Identity(Norb,Norb));
    tt.set_constant(HoppingT*MatrixXcd::Identity(Norb,Norb));
    lattice_1d_1b lattice(Nk,Nt,tt,Ut,Vt,Epulse,MuChem,Norb,SolverOrder,dt);
    // Remember k dependent density and vertex on every rank
    std::vector<CFUNC> density_k(Nk);
    std::vector<CFUNC> vertex(Nk);
    CFUNC rho_loc(Nt,Norb);
    GREEN Gloc,Wloc;

    //Timestep of all k dependent Green's functions in array
    DIST_TIMESTEP gk_all_timesteps(Nk,Nt,Ntau,Norb,-1,true);
    DIST_TIMESTEP wk_all_timesteps(Nk,Nt,Ntau,Norb,+1,true);

    tid_map=gk_all_timesteps.data().tid_map();
    assert(tid==gk_all_timesteps.tid());
    int Nk_rank=0;
    for(int k=0;k<Nk;k++) if(tid_map[k]==tid) Nk_rank++;
    // Information on one rank
    // Map to global lattice_1d_1b
    std::vector<int> kindex_rank;kindex_rank.resize(Nk_rank);
    // Momentum depedent correlators on one rank
    std::vector<gw::kpoint> corrK_rank;corrK_rank.resize(Nk_rank);
    // Convergence error
    cntr::distributed_array<double> convergence_error_ele(Nk,1,true),convergence_error_bos(Nk,1,true);

    if(tid==tid_root){
      Gloc=GREEN(Nt,Ntau,Norb,-1);
      Wloc=GREEN(Nt,Ntau,Norb,1);
    }
    int iq=0;
    for(int k=0;k<Nk;k++){
      vertex[k]=CFUNC(Nt,Norb);
      density_k[k]=CFUNC(Nt,Norb);
      for(int tstp=-1;tstp<=Nt;tstp++){
        cdmatrix vtmp;
        lattice.V(tstp,lattice.kpoints_[k],vtmp);
        vertex[k].set_value(tstp,vtmp);
      }
      if(tid_map[k]==tid){
        kindex_rank[iq]=k;
        corrK_rank[iq]=gw::kpoint(Nt,Ntau,Norb,beta,dt,lattice.kpoints_[k],lattice,MuChem);
        iq++;
      }

    }
    //============================================================================
    //            SELF-CONSISTET SOLUTION, Sigma=0 IN THE BEGINNING
    //============================================================================
    // Solve Sigma=0
    {
      diag::init_G_mat_nointeraction(Nk_rank,corrK_rank,lattice,SolverOrder);
      diag::gather_gk_timestep(-1,Nk_rank,gk_all_timesteps,corrK_rank,kindex_rank);
      // diag::gather_wk_timestep(-1,Nk_rank,wk_all_timesteps,corrK_rank,kindex_rank);
      diag::set_density_k(-1,Norb,gk_all_timesteps,lattice,density_k,kindex_rank,rho_loc);
      diag::get_loc(-1,Ntau,Norb,lattice,Gloc,gk_all_timesteps);
      // diag::get_loc(-1,Ntau,Norb,Nk,lattice,Wloc,wk_all_timesteps);
    }
    // for(int k=0;k<Nk_rank;k++){
    //     {
    //       std::ostringstream name;
    //       name << "test/Gk_0_0" << "_" << kindex_rank[k] <<  ".h5";
    //       corrK_rank[k].G_.write_to_hdf5(name.str().c_str(),"G");
    //     }
    // }

    // {
    //   for(int k=0;k<Nk;k++){
    //     std::ostringstream name;
    //     name << "test/denk_0_" << k << "_" << tid <<  ".h5";
    //     density_k[k].write_to_hdf5(name.str().c_str(),"G");
    //   }
    // }

    { // begin Matsubara Dyson iteration
      if(tid==tid_root){
        print_line_minus(50);
        cout << "     Solution of equilibrium problem " << endl;
        print_line_minus(50);
        start = std::chrono::system_clock::now();
      }
      bool matsubara_converged=false;
      tstp=-1;
      for(int iter=0;iter<=MatsMaxIter;iter++){
        // update propagators via MPI
        diag::gather_gk_timestep(tstp,Nk_rank,gk_all_timesteps,corrK_rank,kindex_rank);
        diag::gather_wk_timestep(tstp,Nk_rank,wk_all_timesteps,corrK_rank,kindex_rank);
        // for(int k=0;k<Nk_rank;k++){
        //   std::ostringstream name;
        //   name << "test/Gk1_-1" << "_" << iter << "_"  << kindex_rank[k] <<  ".h5";
        //   corrK_rank[k].G_.write_to_hdf5(name.str().c_str(),"G");

        //   std::ostringstream name1;
        //   name1 << "test/Wk1_-1" << "_" << iter << "_" << kindex_rank[k] <<  ".h5";
        //   corrK_rank[k].W_.write_to_hdf5(name1.str().c_str(),"G");
        // }
        diag::set_density_k(-1,Norb,gk_all_timesteps,lattice,density_k,kindex_rank,rho_loc);
        diag::get_loc(-1,Ntau,Norb,lattice,Gloc,gk_all_timesteps);
        diag::get_loc(-1,Ntau,Norb,lattice,Wloc,wk_all_timesteps);
        // update mean field and self-energy
	      for(int k=0;k<Nk_rank;k++){
            diag::sigma_Hartree(-1,Norb,corrK_rank[k].SHartree_,lattice,density_k,vertex);
            diag::sigma_Fock(-1,Norb,kindex_rank[k],corrK_rank[k].SFock_,lattice,density_k,vertex);
            diag::sigma_GW(-1,kindex_rank[k],corrK_rank[k].Sigma_,gk_all_timesteps,wk_all_timesteps,lattice,Ntau,Norb);

          // {
          //   std::ostringstream nameH,nameF,nameE,nameS,nameHk;
          //   nameHk << "test/hk_" << iter << "_" << kindex_rank[k]  <<  ".h5";
          //   corrK_rank[k].hk_.write_to_hdf5(nameHk.str().c_str(),"G");

          //   nameH << "test/Hartree_" << iter << "_" << kindex_rank[k]  <<  ".h5";
          //   corrK_rank[k].SHartree_.write_to_hdf5(nameH.str().c_str(),"G");

          //   nameF << "test/Fock_" << iter << "_" << kindex_rank[k] <<  ".h5";
          //   corrK_rank[k].SFock_.write_to_hdf5(nameF.str().c_str(),"G");

          //   nameE << "test/heff_" << iter << "_" << kindex_rank[k] <<  ".h5";
          //   corrK_rank[k].hkeff_.write_to_hdf5(nameE.str().c_str(),"G");

          //   nameS << "test/Sigmak_-1_" << iter << "_" << kindex_rank[k] <<  ".h5";
          //   corrK_rank[k].Sigma_.write_to_hdf5(nameS.str().c_str(),"G");
          // }
        }
	      // solve Dyson equation
        for(int k=0;k<Nk_rank;k++){
          convergence_error_ele.block(kindex_rank[k])[0]=corrK_rank[k].step_dyson_with_error(tstp,iter,SolverOrder,lattice);
          diag::get_Polarization_Bubble(tstp,Norb,Ntau,kindex_rank[k],corrK_rank[k].P_,gk_all_timesteps,lattice);
          convergence_error_bos.block(kindex_rank[k])[0]=corrK_rank[k].step_W_with_error(tstp,SolverOrder,lattice);
        }


        convergence_error_ele.mpi_bcast_all();
        convergence_error_bos.mpi_bcast_all();
        double err_ele=0.0,err_bos=0.0;
        for(int k=0;k<Nk;k++){
          err_ele+=convergence_error_ele.block(k)[0];
          err_bos+=convergence_error_bos.block(k)[0];
        }

        for(int k=0;k<Nk_rank;k++){
          cntr::force_matsubara_hermitian(corrK_rank[k].G_);
          cntr::force_matsubara_hermitian(corrK_rank[k].Sigma_);
          cntr::force_matsubara_hermitian(corrK_rank[k].W_);
          cntr::force_matsubara_hermitian(corrK_rank[k].P_);
        }

        diag::set_density_k(-1,Norb,gk_all_timesteps,lattice,density_k,kindex_rank,rho_loc);
        diag::get_loc(-1,Ntau,Norb,lattice,Gloc,gk_all_timesteps);
        diag::get_loc(-1,Ntau,Norb,lattice,Wloc,wk_all_timesteps);

        if(tid==tid_root){
          cdmatrix tmp;
          rho_loc.get_value(tstp,tmp);
          cout << "tstp= " << tstp << " iter:  " << iter << " err ele: " << err_ele << " err bos: " << err_bos << " - Local density matrix: " << tmp << " and trace: " << tmp.trace() << endl;
        }
        if(err_ele<MatsMaxErr && err_bos<MatsMaxErr){
            matsubara_converged=true;
            break;
        }


      }

      if(!matsubara_converged && tid==tid_root){
        cout << endl;
        cout << " Matsubara iteration not converged! Exiting ... " << endl;
        // should end here ....
        return 0;
      }

      if(tid==tid_root){
        print_line_dot(50);
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        cout << "Time [equilibrium calculation] = " << elapsed_seconds.count() << "s\n\n";
      }
    } // end Matsubara Dyson iteration
    
    // TODO: Set guess propagator to free + HartreeFock 

    //============================================================================
    //           BOOTSTRAPPING PHASE
    //============================================================================
    { // begin bootstrapping
      if(tid==tid_root){
        print_line_minus(50);
        cout << "     Time propagation: bootstrapping phase " << endl;
        print_line_minus(50);
      }
      bool bootstrap_converged;
      if(tid==tid_root){
        start = std::chrono::system_clock::now();
        bootstrap_converged=false;
      }

      // For T=0 just set up all propagators
      for(tstp=0; tstp<=SolverOrder; tstp++){
        int n;
        int n1= 0;
        int n2= tstp;

        diag::extrapolate_rho(tstp,Nk,density_k);
        diag::extrapolate_timestep_G(tstp-1,Nk_rank,(SolverOrder<tstp-1 ? SolverOrder : tstp-1),Nt,corrK_rank);
        diag::extrapolate_timestep_W(tstp-1,Nk_rank,(SolverOrder<tstp-1 ? SolverOrder : tstp-1),Nt,corrK_rank);

        for (int iter = 0; iter <= BootstrapMaxIter; iter++){
          double err_ele=0.0,err_bos=0.0;
          // for(int k=0;k<Nk_rank;k++){
          //     // Tmp print 
          //     std::ostringstream name;
          //     name << "test/Gk1_" << tstp << "_" << iter << "_"  << kindex_rank[k] <<  ".h5";
          //     corrK_rank[k].G_.write_to_hdf5(name.str().c_str(),"G");

          //     std::ostringstream name1;
          //     name1 << "test/Wk1_" << tstp << "_" << iter << "_" << kindex_rank[k] <<  ".h5";
          //     corrK_rank[k].W_.write_to_hdf5(name1.str().c_str(),"G");
          // }
          // Internal loop to update previous times
          for(int n=n1;n<=n2;n++){
            diag::gather_gk_timestep(n,Nk_rank,gk_all_timesteps,corrK_rank,kindex_rank);
            diag::gather_wk_timestep(n,Nk_rank,wk_all_timesteps,corrK_rank,kindex_rank);

            
            diag::set_density_k(n,Norb,gk_all_timesteps,lattice,density_k,kindex_rank,rho_loc);          
            diag::get_loc(n,Ntau,Norb,lattice,Gloc,gk_all_timesteps);

            for(int k=0;k<Nk_rank;k++){
              diag::sigma_Hartree(n,Norb,corrK_rank[k].SHartree_,lattice,density_k,vertex);
              diag::sigma_Fock(n,Norb,kindex_rank[k],corrK_rank[k].SFock_,lattice,density_k,vertex);
              diag::sigma_GW(n,kindex_rank[k],corrK_rank[k].Sigma_,gk_all_timesteps,wk_all_timesteps,lattice,Ntau,Norb);
              diag::get_Polarization_Bubble(n,Norb,Ntau,kindex_rank[k],corrK_rank[k].P_,gk_all_timesteps,lattice);

              // {
              //   std::ostringstream nameH,nameF,nameE,nameS,nameHk;
              //   nameHk << "test/hk_" << tstp << "_" << iter << "_" << kindex_rank[k]  <<  ".h5";
              //   corrK_rank[k].hk_.write_to_hdf5(nameHk.str().c_str(),"G");

              //   nameH << "test/Hartree_" << tstp << "_" << iter << "_" << kindex_rank[k]  <<  ".h5";
              //   corrK_rank[k].SHartree_.write_to_hdf5(nameH.str().c_str(),"G");

              //   nameF << "test/Fock_" << tstp << "_" << iter << "_" << kindex_rank[k] <<  ".h5";
              //   corrK_rank[k].SFock_.write_to_hdf5(nameF.str().c_str(),"G");

              //   nameE << "test/heff_" << tstp << "_" << iter << "_" << kindex_rank[k] <<  ".h5";
              //   corrK_rank[k].hkeff_.write_to_hdf5(nameE.str().c_str(),"G");

              //   nameS << "test/Sigmak_" << tstp << "_" << iter << "_" << kindex_rank[k] <<  ".h5";
              //   corrK_rank[k].Sigma_.write_to_hdf5(nameS.str().c_str(),"G");
                
              //   cdmatrix tmp1;
              //   corrK_rank[k].P_.density_matrix(tstp,tmp1);
              //   std::cout << "Pol outside2 " << tstp << " " << iter << " " << kindex_rank[k]   << tmp1 << std::endl;
                
              //   std::ostringstream nameP;
              //   nameP << "test/P_" << tstp  <<"_" << iter << "_" << kindex_rank[k] <<  ".h5";
              //   corrK_rank[k].P_.write_to_hdf5(nameP.str().c_str(),"G");

              // }

              convergence_error_ele.block(kindex_rank[k])[0]=corrK_rank[k].step_dyson_with_error(n,iter,SolverOrder,lattice);
              convergence_error_bos.block(kindex_rank[k])[0]=corrK_rank[k].step_W_with_error(n,SolverOrder,lattice);
              
              // {
              //   std::ostringstream name0;
              //   name0 << "test/Gk2_" << tstp << "_" << iter << "_"  << kindex_rank[k] <<  ".h5";
              //   corrK_rank[k].G_.write_to_hdf5(name0.str().c_str(),"G");
                  
              //   std::ostringstream name1;
              //   name1 << "test/Wk2_" << tstp   <<"_" << iter << "_" << kindex_rank[k] <<  ".h5";
              //   corrK_rank[k].W_.write_to_hdf5(name1.str().c_str(),"G");  
              // }
            }

            convergence_error_ele.mpi_bcast_all();
            convergence_error_bos.mpi_bcast_all();
            
            for(int k=0;k<Nk;k++){
              err_ele+=convergence_error_ele.block(k)[0];
              err_bos+=convergence_error_bos.block(k)[0];
            }


            diag::gather_gk_timestep(n,Nk_rank,gk_all_timesteps,corrK_rank,kindex_rank);
            diag::gather_wk_timestep(n,Nk_rank,wk_all_timesteps,corrK_rank,kindex_rank);
            diag::set_density_k(n,Norb,gk_all_timesteps,lattice,density_k,kindex_rank,rho_loc);
            diag::get_loc(n,Ntau,Norb,lattice,Gloc,gk_all_timesteps);
            diag::get_loc(n,Ntau,Norb,lattice,Wloc,wk_all_timesteps);
            

            if(tid==tid_root){
              cdmatrix tmp;
              rho_loc.get_value(tstp,tmp);
              cout << "bootstrap iteration tstp: " << tstp   << " iter:  " << iter << " err ele: " << err_ele << " err bos: " << err_bos << " - " << tmp << " " << tmp.trace() << endl;
            }
          }

          if(err_ele<BootstrapMaxErr && err_bos<BootstrapMaxErr && iter>2){
            bootstrap_converged=true;
             break;
          }
        }

        if(!bootstrap_converged && tid==tid_root){
           cout << endl;
            cout << " Bootstrap iteration not converged! Exiting ... " << endl;
            // should end here ....
           return 0;
        }
        if(tid==tid_root){
          print_line_dot(50);
          end = std::chrono::system_clock::now();
          std::chrono::duration<double> elapsed_seconds= end -start;
          cout << "Time [bootstrapping] = " << elapsed_seconds.count() << "s\n";
        }
      }
    }
      // } // end iter
    // }// end bootstrapping
    //============================================================================
    //             TIME PROPAGATION
    //============================================================================
    { // begin propagation loop
      if(tid==tid_root){
        print_line_minus(50);
        cout << "               Time propagation" << endl;
        print_line_minus(50);
        start = std::chrono::system_clock::now();
      }

      for(tstp = SolverOrder+1; tstp <= Nt; tstp++){
	      // Predictor: extrapolation
        diag::extrapolate_rho(tstp,Nk,density_k);
        diag::extrapolate_timestep_G(tstp-1,Nk_rank,SolverOrder,Nt,corrK_rank);
        diag::extrapolate_timestep_W(tstp-1,Nk_rank,SolverOrder,Nt,corrK_rank);
	      // Corrector
	      for (int iter=0; iter < CorrectorSteps; iter++){
          // for(int k=0;k<Nk_rank;k++){
          //     // Tmp print 
          //     std::ostringstream name;
          //     name << "test/Gk1_" << tstp << "_" << iter << "_"  << kindex_rank[k] <<  ".h5";
          //     corrK_rank[k].G_.write_to_hdf5(name.str().c_str(),"G");

          //     std::ostringstream name1;
          //     name1 << "test/Wk1_" << tstp << "_" << iter << "_" << kindex_rank[k] <<  ".h5";
          //     corrK_rank[k].W_.write_to_hdf5(name1.str().c_str(),"G");
          // }
          // Gather propagators
          diag::gather_gk_timestep(tstp,Nk_rank,gk_all_timesteps,corrK_rank,kindex_rank);
          diag::gather_wk_timestep(tstp,Nk_rank,wk_all_timesteps,corrK_rank,kindex_rank);

          diag::set_density_k(tstp,Norb,gk_all_timesteps,lattice,density_k,kindex_rank,rho_loc);          
          diag::get_loc(tstp,Ntau,Norb,lattice,Gloc,gk_all_timesteps);
          diag::get_loc(tstp,Ntau,Norb,lattice,Wloc,wk_all_timesteps);

	  // update mean field and SigmaGW
          for(int k=0;k<Nk_rank;k++){
            diag::sigma_Hartree(tstp,Norb,corrK_rank[k].SHartree_,lattice,density_k,vertex);
            diag::sigma_Fock(tstp,Norb,kindex_rank[k],corrK_rank[k].SFock_,lattice,density_k,vertex);
            diag::sigma_GW(tstp,kindex_rank[k],corrK_rank[k].Sigma_,gk_all_timesteps,wk_all_timesteps,lattice,Ntau,Norb);
          }
          // Solve dyson, update polarization and solve two particle self-consistency
          for(int k=0;k<Nk_rank;k++){
            // {
            //     std::ostringstream nameH,nameF,nameE,nameS,nameHk;
            //     nameHk << "test/hk_" << tstp << "_" << iter << "_" << kindex_rank[k]  <<  ".h5";
            //     corrK_rank[k].hk_.write_to_hdf5(nameHk.str().c_str(),"G");

            //     nameH << "test/Hartree_" << tstp << "_" << iter << "_" << kindex_rank[k]  <<  ".h5";
            //     corrK_rank[k].SHartree_.write_to_hdf5(nameH.str().c_str(),"G");

            //     nameF << "test/Fock_" << tstp << "_" << iter << "_" << kindex_rank[k] <<  ".h5";
            //     corrK_rank[k].SFock_.write_to_hdf5(nameF.str().c_str(),"G");

            //     nameE << "test/heff_" << tstp << "_" << iter << "_" << kindex_rank[k] <<  ".h5";
            //     corrK_rank[k].hkeff_.write_to_hdf5(nameE.str().c_str(),"G");

            //     nameS << "test/Sigmak_" << tstp << "_" << iter << "_" << kindex_rank[k] <<  ".h5";
            //     corrK_rank[k].Sigma_.write_to_hdf5(nameS.str().c_str(),"G");
                
            //     cdmatrix tmp1;
            //     corrK_rank[k].P_.density_matrix(tstp,tmp1);
            //     std::cout << "Pol outside2 " << tstp << " " << iter << " " << kindex_rank[k]   << tmp1 << std::endl;
                
            //     std::ostringstream nameP;
            //     nameP << "test/P_" << tstp  <<"_" << iter << "_" << kindex_rank[k] <<  ".h5";
            //     corrK_rank[k].P_.write_to_hdf5(nameP.str().c_str(),"G");

            // }

            convergence_error_ele.block(kindex_rank[k])[0]=corrK_rank[k].step_dyson_with_error(tstp,iter,SolverOrder,lattice);
            diag::get_Polarization_Bubble(tstp,Norb,Ntau,kindex_rank[k],corrK_rank[k].P_,gk_all_timesteps,lattice);
            convergence_error_bos.block(kindex_rank[k])[0]=corrK_rank[k].step_W_with_error(tstp,SolverOrder,lattice);

            // {
            //     std::ostringstream name0;
            //     name0 << "test/Gk2_" << tstp << "_" << iter << "_"  << kindex_rank[k] <<  ".h5";
            //     corrK_rank[k].G_.write_to_hdf5(name0.str().c_str(),"G");
                  
            //     std::ostringstream name1;
            //     name1 << "test/Wk2_" << tstp   <<"_" << iter << "_" << kindex_rank[k] <<  ".h5";
            //     corrK_rank[k].W_.write_to_hdf5(name1.str().c_str(),"G");  
            // }

          }
          convergence_error_ele.mpi_bcast_all();
          convergence_error_bos.mpi_bcast_all();
          double err_ele=0.0,err_bos=0.0;
          for(int k=0;k<Nk;k++){
            err_ele+=convergence_error_ele.block(k)[0];
            err_bos+=convergence_error_bos.block(k)[0];
          }
          diag::set_density_k(tstp,Norb,gk_all_timesteps,lattice,density_k,kindex_rank,rho_loc);
          diag::get_loc(tstp,Ntau,Norb,lattice,Gloc,gk_all_timesteps);
          diag::get_loc(tstp,Ntau,Norb,lattice,Wloc,wk_all_timesteps);
          if(tid==tid_root){
            cdmatrix tmp;
            rho_loc.get_value(tstp,tmp);
            cout << "Iteration tstp: " << tstp << " iter:  " << iter << " err ele: " << err_ele << " err bos: " << err_bos << " - " << tmp << " " << tmp.trace() << endl;
          }
          if(err_ele<TimeMaxErr && err_bos<TimeMaxErr && iter>1){
            if(tid==tid_root){
              print_line_dot(50);
            }
            break;
          }
        } // end propagation loop
      } // Iteration over time
    }
    {

      CFUNC Ekin(Nt,1),Epot(Nt,1),Epulsetmp(Nt,1);
      cdmatrix tmp(Norb,Norb);
      for(tstp=-1; tstp <= Nt; tstp++){
          if(tid==tid_root){
            tmp(0,0)=diag::KineticEnergy(tstp,lattice,density_k);
            Ekin.set_value(tstp,tmp);
            tmp(0,0)=Epulse[tstp+1];
            Epulsetmp.set_value(tstp,tmp);
          }
          tmp(0,0) = diag::CorrelationEnergy(tstp,Nk_rank,SolverOrder,beta,dt,corrK_rank,lattice);
          if(tid==tid_root){
            Epot.set_value(tstp,tmp);
          }
      }
      // output
      if(tid==tid_root){
        hid_t file_id = open_hdf5_file("out/data_gw.h5");
        hid_t group_id;
        // Param
        group_id = create_group(file_id, "parm");
        store_int_attribute_to_hid(group_id, "Nk", Nk); 
        store_int_attribute_to_hid(group_id, "Nt", Nt);
        store_int_attribute_to_hid(group_id, "Ntau", Ntau);
        store_double_attribute_to_hid(group_id, "beta", beta);
        store_double_attribute_to_hid(group_id, "dt", dt);
        store_double_attribute_to_hid(group_id, "HubbardU", HubbardU);
        store_double_attribute_to_hid(group_id, "HoppingT", HoppingT);
        store_double_attribute_to_hid(group_id, "MuChem", MuChem);
        lattice.Apulse_.write_to_hdf5(group_id,"Apulse");
        Epulsetmp.write_to_hdf5(group_id,"Epulse");
        close_group(group_id);
        group_id = create_group(file_id, "obs");
        Ekin.write_to_hdf5(group_id,"Ekin");
        Epot.write_to_hdf5(group_id,"Epot");
        close_group(group_id);
        // Print local green's function
        if(SaveGreen==1){
	        
	        group_id = create_group(file_id, "Gloc");
	        store_herm_greens_function(group_id, Gloc);
	        close_group(group_id);
	  
	        group_id = create_group(file_id, "Wloc");
	        store_herm_greens_function(group_id, Wloc);
	        close_group(group_id);     
	     }
       close_hdf5_file(file_id);
      }
      if(SaveMomentum==1){
        for(int k=0;k<Nk_rank;k++){
          if(tid_map[k]==tid){
            char fnametmp[1000];
            std::sprintf(fnametmp,"out/k%d.h5",k);
            std::cout << "writing hdf5 data to " << fnametmp << std::endl;
            hid_t file_id = open_hdf5_file(std::string(fnametmp));
            corrK_rank[k].G_.write_to_hdf5_slices(file_id,"G",output);
            corrK_rank[k].W_.write_to_hdf5_slices(file_id,"W",output);
            close_hdf5_file(file_id);
          }
        }
      }
    } // end output

    if(tid=tid_root){
      end_tot = std::chrono::system_clock::now();
      std::chrono::duration<double> runtime_seconds = end_tot-start_tot;

      cout << endl;
      cout << endl;
      cout << "Time [total] = " << runtime_seconds.count() << "s\n";
      print_line_star(60);
    }
    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
  }
