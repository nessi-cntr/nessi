#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <cstring>
#include <chrono>

// contour library headers
#include "cntr/cntr.hpp"
#include "cntr/hdf5/hdf5_interface.hpp"
#include "cntr/hdf5/hdf5_interface_cntr.hpp"

// local headers to include
#include "formats.hpp"
#include "cntr/utils/read_inputfile.hpp"
#include "Holstein_impurity_decl.hpp"
#include "Holstein_utils_decl.hpp"


using namespace std;

#define CFUNC cntr::function<double>
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define CPLX complex<double>


//==============================================================================
//         main program
//==============================================================================
int main(int argc,char *argv[]){
	// use maximum order
  	int SolverOrder=MAX_SOLVE_ORDER;
	//..................................................
	//              Constant
	//..................................................
	cdmatrix sig3=cdmatrix::Zero(2,2);
	sig3(0,0)=1.0;
	sig3(1,1)=-1.0;
	cdmatrix sig1=cdmatrix::Zero(2,2);
	sig1(0,1)=1.0;
	sig1(1,0)=1.0;
	//..................................................
	//                input
	//..................................................
	int SaveGreen;
	int Nt,Ntau,MatsMaxIter,CorrectorSteps;
	int BootstrapMaxIter;
	double Hopping,El_Ph_g,Phfreq_w0,beta,dt,MuChem_MF,MatsMaxErr,BootstrapMaxErr;
	double MuCHem;
	std::vector<double> dHopping,sc_field;
	char* flin;

	//..................................................
	//                output
	//..................................................
	char* flout;
	int OutEvery;
	//..................................................
	//                internal
	//..................................................
	int tstp;
	double err;
	int Norb=2;
	GREEN G,Sigma,Hyb,gtemp,Hyb_Sig; //For Fermion in the Nambu form
	GREEN D,Pi,D0,D0_Pi,Pi_D0,dtemp; //For the phonon propagator
	CFUNC g_elph_t,J_hop_t,h0_imp_t; // note g_elph_t=El_Ph_g\sig_3
	CFUNC Ekin_t,Enx_MF_t,Enx_corr_t,Eph_t,X_t;
	cdmatrix h0_imp;
	//..................................................
	//                timer
	//..................................................
	chrono::time_point<std::chrono::system_clock> start, end, start_tot, end_tot;

	print_line_star(60);
	cout << "   Test program: Holstein model on Bethe lattice in Midal approximation" << endl;
	print_line_star(60);

	start_tot = chrono::system_clock::now();

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

			// system parameters
			find_param(flin,"__Hopping=",Hopping);
			find_param(flin,"__El_Ph_g=",El_Ph_g);
			find_param(flin,"__Phfreq_w0=",Phfreq_w0);
			find_param(flin,"__MuChem_MF=",MuChem_MF);
			find_param(flin,"__beta=",beta);
			// ramp parameters

			// solver parameters
			find_param(flin,"__SaveGreen=",SaveGreen);
			find_param(flin,"__Nt=",Nt);
			find_param(flin,"__Ntau=",Ntau);
			find_param(flin,"__dt=",dt);
			find_param(flin,"__MatsMaxIter=",MatsMaxIter);
			find_param(flin,"__MatsMaxErr=",MatsMaxErr);
			find_param(flin,"__BootstrapMaxIter=",BootstrapMaxIter);
			find_param(flin,"__BootstrapMaxErr=",BootstrapMaxErr);
			find_param(flin,"__CorrectorSteps=",CorrectorSteps);

			// excitation parameters
			find_param_tvector(flin, "__dHopping=",dHopping,Nt);// The vector size is set to Nt+2,i.e. [-t,Nt]
			find_param_tvector(flin, "__sc_field=",sc_field,Nt);

			// output parameters
			find_param(flin,"__OutEvery=",OutEvery);
			// output file prefix
			flout=argv[2];

		}
		//============================================================================
		//               (IV) INITIALIZE GREEN'S FUNCTIONS
		//============================================================================
		{
			//allocation
			//system parameters
			g_elph_t = CFUNC(Nt,Norb);
			J_hop_t = CFUNC(Nt,Norb);
			h0_imp=cdmatrix::Zero(Norb,Norb);
			h0_imp_t = CFUNC(Nt,Norb);


			//electron
			G = GREEN(Nt,Ntau,Norb,FERMION);
			Hyb = GREEN(Nt,Ntau,Norb,FERMION);
			Sigma = GREEN(Nt,Ntau,Norb,FERMION);
			Hyb_Sig = GREEN(Nt,Ntau,Norb,FERMION); //Delta+Sig_imp

			//phonon
			D = GREEN(Nt,Ntau,1,BOSON);
			Pi = GREEN(Nt,Ntau,1,BOSON);
			D0 = GREEN(Nt,Ntau,1,BOSON);
			D0_Pi = GREEN(Nt,Ntau,1,BOSON);
			Pi_D0 = GREEN(Nt,Ntau,1,BOSON);

			//observables:energies
			Ekin_t = CFUNC(Nt,1);
			Enx_MF_t = CFUNC(Nt,1);
			Enx_corr_t = CFUNC(Nt,1);
			Eph_t = CFUNC(Nt,1);
			X_t = CFUNC(Nt,1);


			//initialization
			//system parameters in equilibrium
			g_elph_t.set_constant(El_Ph_g*sig3);
			J_hop_t.set_constant(Hopping*sig3);
			h0_imp=-MuChem_MF*sig3;//[note]: This includes chemical potential.
			h0_imp_t.set_constant(h0_imp);

			//system parameters with excitations
			for(int it=-1 ; it<=Nt ; it++) J_hop_t.set_value(it,(Hopping+dHopping[it+1])*sig3);

			for(int it=-1 ; it<=Nt ; it++) h0_imp_t.set_value(it,h0_imp+sc_field[it+1]*sig1);

			//electron[note]:this part can be improved for better guess.
			//cntr::green_from_H(G, 0.0, h0_imp, beta, dt);
			cntr::green_equilibrium_mat_bethe(G, beta);

			//Hybridization function.
			tstp=-1;
			Hyb.set_timestep(tstp,G);
			Hyb.right_multiply(tstp,J_hop_t);
			Hyb.left_multiply(tstp,J_hop_t);

			//initialize D0
			cntr::green_single_pole_XX(D0,Phfreq_w0,beta,dt);
			for(int it=-1 ; it<=Nt ; it++) D0.smul(it,2.0);

		}
		//============================================================================
		//            SELF-CONSISTET SOLUTION, Sigma=0 IN THE BEGINNING
		//============================================================================
		{ // begin Matsubara Dyson iteration
			print_line_minus(50);
			cout << "     Solution of equilibrium problem " << endl;
			print_line_minus(50);

			start = std::chrono::system_clock::now();

			bool matsubara_converged=false;
			tstp=-1;

			gtemp = GREEN(SolverOrder,Ntau,Norb,-1);
			gtemp.set_timestep(tstp,G);

			bool FIXPOINT = false;

			for(int iter=0;iter<=MatsMaxIter;iter++){

				//////////////////
				//Solving impurity
				//////////////////

				//update self-energy
				Hols::Sigma_uMig_sc(tstp,G, D0, g_elph_t,Sigma);

				//solve Dyson for impurity G_imp=[(i\partial_t I + mu*sig3)delta_c-Hyb-Sigma_imp]^{-1}
				Hyb_Sig.set_timestep(tstp,Hyb);
				Hyb_Sig.incr_timestep(tstp,Sigma,1.0);

				double SC_seed=(iter<5 ? 0.1 : 0.0);

				cdmatrix h0_imp_sc=h0_imp+SC_seed*sig1;
				h0_imp_t.set_value(tstp,h0_imp_sc);

				if(FIXPOINT) cntr::dyson_mat(G, 0.0, h0_imp_t, Hyb_Sig, beta, SolverOrder, CNTR_MAT_FIXPOINT);
				else cntr::dyson_mat(G, 0.0, h0_imp_t, Hyb_Sig, beta, SolverOrder, CNTR_MAT_FOURIER);

				//mixing
				G.smul(tstp,0.5);
				G.incr_timestep(tstp,gtemp,0.5);

				////////////////////////////////////////////
				//Lattice self-consistency for Bethe lattice
				////////////////////////////////////////////
				//Update hybridization
				Hyb.set_timestep(tstp,G);
				Hyb.right_multiply(tstp,J_hop_t);
				Hyb.left_multiply(tstp,J_hop_t);

				// self-consistency check
				err = cntr::distance_norm2(tstp,G,gtemp);
				if(FIXPOINT) cout << "Fixpoint method: iteration : " << iter << "  |  Error = " << err << endl;
				else cout << "Fourier method: iteration : " << iter << "  |  Error = " << err << endl;
				// Output temporal values
				cdmatrix rho_M=cdmatrix::Zero(2,2);
				G.get_mat(Ntau,rho_M);
				cout << "rho_00:"<<-rho_M(0,0)<<" rho_01"<<-rho_M(0,1)<<endl;

				if(err<MatsMaxErr && FIXPOINT){
					matsubara_converged=true;
					break;
				}
				else if(err<MatsMaxErr && !FIXPOINT){
					FIXPOINT=true;
				}
				gtemp.set_timestep(tstp,G);
			}

			if(!matsubara_converged){
				cout << endl;
				cout << " Matsubara iteration not converged! Exiting ... " << endl;
				// should end here ....
				return 0;
			}

			print_line_dot(50);
			end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end-start;
			cout << "Time [equilibrium calculation] = " << elapsed_seconds.count() << "s\n\n";

		} // end Matsubara Dyson iteration

		//============================================================================
		//           SET TSTP=kt
		//============================================================================

		{
			cntr::set_tk_from_mat(G,SolverOrder);
			cntr::set_tk_from_mat(Hyb,SolverOrder);
		}

		//============================================================================
		//           BOOTSTRAPPING PHASE
		//============================================================================
		{ // begin bootstrapping

			print_line_minus(50);
			cout << "     Time propagation: bootstrapping phase " << endl;
			print_line_minus(50);

			start = std::chrono::system_clock::now();

			bool bootstrap_converged=false;

			for(tstp=0; tstp<=SolverOrder; tstp++)
				gtemp.set_timestep(tstp,G);

			for (int iter = 0; iter <= BootstrapMaxIter; iter++) {

				//////////////////
				//Solving impurity
				//////////////////

				// update self-energy
                Hols::Sigma_uMig_sc(G, D0, g_elph_t, Sigma, SolverOrder);

				//solve Dyson for impurity
				for(tstp=0; tstp<=SolverOrder; tstp++){
					Hyb_Sig.set_timestep(tstp,Hyb);
					Hyb_Sig.incr_timestep(tstp,Sigma,1.0);
				}

				cntr::dyson_start(G, 0.0, h0_imp_t, Hyb_Sig, beta, dt, SolverOrder);

				////////////////////////////////////////////
				//Lattice self-consistency for Bethe lattice
				////////////////////////////////////////////
				//Update hybridization
				for(tstp=0; tstp<=SolverOrder; tstp++){
					Hyb.set_timestep(tstp,G);
					Hyb.right_multiply(tstp,J_hop_t);
					Hyb.left_multiply(tstp,J_hop_t);
				}

				// self-consistency check
				err=0.0;
				for(tstp=0; tstp<=SolverOrder; tstp++) {
					err += cntr::distance_norm2(tstp,G,gtemp);
				}
				cout << "bootstrap iteration : " << iter << "  |  Error = " << err << endl;
				if(err<BootstrapMaxErr && iter>2){
					bootstrap_converged=true;
					break;
				}
				for(tstp=0; tstp<=SolverOrder; tstp++) {
					gtemp.set_timestep(tstp,G);
				}
			}

			print_line_dot(50);
			end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds= end -start;
			cout << "Time [bootstrapping] = " << elapsed_seconds.count() << "s\n";

			if(!bootstrap_converged){
				cout << endl;
				cout << " Bootstrap iteration not converged! Exiting ... " << endl;
				// should end here ....
				return 0;
			}

		} // end bootstrapping

		//============================================================================
		//             TIME PROPAGATION
		//============================================================================
		{ // begin propagation loop
			print_line_minus(50);
			cout << "               Time propagation" << endl;
			print_line_minus(50);

			start = std::chrono::system_clock::now();

			for(tstp = SolverOrder+1; tstp <= Nt; tstp++){

				if (tstp%50 ==0) cout<<"tstp:"<<tstp<<endl;

				// Predictor: extrapolation
				cntr::extrapolate_timestep(tstp-1,G,SolverOrder);
				cntr::extrapolate_timestep(tstp-1,Hyb,SolverOrder);

				// Corrector
				for (int iter=0; iter < CorrectorSteps; iter++){


					//////////////////
					//Solving impurity
					//////////////////

					//update self-energy
					Hols::Sigma_uMig_sc(tstp,G, D0, g_elph_t,Sigma);

					//solve Dyson for impurity
					Hyb_Sig.set_timestep(tstp,Hyb);
					Hyb_Sig.incr_timestep(tstp,Sigma,1.0);
					cntr::dyson_timestep(tstp, G, 0.0, h0_imp_t, Hyb_Sig, beta, dt, SolverOrder);

					////////////////////////////////////////////
					//Lattice self-consistency for Bethe lattice
					////////////////////////////////////////////

					//Update hybridization
					Hyb.set_timestep(tstp,G);
					Hyb.right_multiply(tstp,J_hop_t);
					Hyb.left_multiply(tstp,J_hop_t);
				}
			}

			print_line_dot(50);
			end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end-start;
			cout << "Time [KBEs] = " << elapsed_seconds.count() << "s\n";

		} // end propagation loop
        
        //============================================================================
		//             Post Process
		//============================================================================
        {
            for(tstp = -1 ; tstp<=Nt ; tstp++) D.set_timestep(tstp,D0);
        }
        
		//============================================================================
		//             Evaluateion of energies
		//============================================================================
		{
			print_line_minus(50);
			cout << "               Evaluation of Energies" << endl;
			print_line_minus(50);

			cdmatrix Mat_tmp(1,1);
			cdmatrix XX(1,1),PP(1,1);

			// Displacement X and Enk_MF
			cdmatrix rho_M=cdmatrix::Zero(2,2);
			G.get_mat(Ntau,rho_M);
			rho_M *= -1.0;

			Mat_tmp(0,0) = -2.0*El_Ph_g/Phfreq_w0*(rho_M(0,0)-rho_M(1,1));
			X_t.set_constant(Mat_tmp);

			Mat_tmp(0,0) *= El_Ph_g * (rho_M(0,0)-rho_M(1,1));


			// Ekin for two spin
			for(tstp = -1; tstp <= Nt ; tstp++){
				Mat_tmp(0,0) = 2.0*cntr::correlation_energy(tstp, G, Hyb, beta, dt, SolverOrder);
				Ekin_t.set_value(tstp,Mat_tmp);
			}
			// Enx_corr
			for(tstp = -1; tstp <= Nt ; tstp++){
				Mat_tmp(0,0) = 2.0*cntr::correlation_energy(tstp, G, Sigma, beta, dt, SolverOrder);
				Enx_corr_t.set_value(tstp,Mat_tmp);
			}

			// Eph
            for(tstp = -1; tstp <= Nt ; tstp++) Pi.set_timestep_zero(tstp);
            Hols::evaluate_phonon_energy_qu(Eph_t, D, Pi, SolverOrder, beta, dt, Phfreq_w0);
        }
		//============================================================================
		//             OUT PUT
		//============================================================================

		{   char file_occ[255];
			strcpy(file_occ,flout);
			strcat(file_occ,"_obs.dat");
			ofstream f_out;
			f_out.open (file_occ);

			//observables
			cdmatrix rho(Norb,Norb);
			cdmatrix XX(1,1);
			cdmatrix ekin(1,1);
			cdmatrix enx_mf(1,1);
			cdmatrix enx_corr(1,1);
			cdmatrix eph(1,1);
			cdmatrix etot(1,1);

			int nobs=0;
			f_out<<"#0:time ";

			for(int a=0; a<Norb; a++){
				for(int b=0; b<Norb; b++){
					nobs++;
					f_out << nobs << ":rho(" << a << "," << b << ").re ";
					nobs++;
					f_out << nobs << ":rho(" << a << "," << b << ").im ";
				}
			}

			nobs++;
			f_out<<nobs<<":XX ";
			nobs++;
			f_out<<nobs<<":Ekin ";
			nobs++;
			f_out<<nobs<<":Exn_MF ";
			nobs++;
			f_out<<nobs<<":Exn_corr ";
			nobs++;
			f_out<<nobs<<":Eph ";
			nobs++;
			f_out<<nobs<<":Etot ";

			f_out<<endl;
			for(tstp=0; tstp <= Nt; tstp++){

				G.density_matrix(tstp, rho);
				D.get_les(tstp,tstp,XX);
				XX(0,0)*=CPLX(0.0,1.0);
				Ekin_t.get_value(tstp,ekin);
				Enx_MF_t.get_value(tstp,enx_mf);
				Enx_corr_t.get_value(tstp,enx_corr);
				Eph_t.get_value(tstp,eph);
				etot = ekin + enx_mf + enx_corr + eph;

				f_out << tstp*dt << "  " ;

				for(int a=0; a<Norb; a++){
					for(int b=0; b<Norb; b++){
						f_out << rho(a,b).real() << " " << rho(a,b).imag()<<" ";
					}
				}

				f_out<< XX(0,0).real()<<" ";
				f_out<< ekin(0,0).real()<<" ";
				f_out<< enx_mf(0,0).real()<<" ";
				f_out<< enx_corr(0,0).real()<<" ";
				f_out<< eph(0,0).real()<<" ";
				f_out<< etot(0,0).real()<<" ";

				f_out << endl;
			}
			f_out.close();
		}
			////////////
			/// HDF5////
			////////////
		{
			char fnametmp[1000];
			std::sprintf(fnametmp,"%s_green.h5",flout);
			hid_t file_id = open_hdf5_file(std::string(fnametmp));
			hid_t group_id = create_group(file_id,"parm");

			//parameters
			store_int_attribute_to_hid(group_id,"Nt",Nt);
			store_int_attribute_to_hid(group_id,"Ntau",Ntau);
			store_double_attribute_to_hid(group_id,"dt",dt);
			store_double_attribute_to_hid(group_id,"beta",beta);
			close_group(group_id);

			//electron
			//G at each time step t_n: G^R(t_n,t), G^<(t,t_n), G^tv(t_n,\tau) for t<t_n
			group_id = create_group(file_id, "Gloc_slice");
            G.write_to_hdf5_slices(group_id,OutEvery);
            close_group(group_id);
            // G(t_rel,t_av) for t_av = [0,OutEvery*dt,2*OutEvery*dt...].
            group_id = create_group(file_id, "Gloc_tavrel");
            G.write_to_hdf5_tavtrel(group_id,OutEvery);
            close_group(group_id);

			//phonon
			group_id = create_group(file_id, "Dloc_slice");
            D.write_to_hdf5_slices(group_id,OutEvery);
            close_group(group_id);

            group_id = create_group(file_id, "Dloc_tavrel");
            D.write_to_hdf5_tavtrel(group_id,OutEvery);
            close_group(group_id);

            close_hdf5_file(file_id);


		} // end output

		end_tot = std::chrono::system_clock::now();
		std::chrono::duration<double> runtime_seconds = end_tot-start_tot;

		cout << endl;
		cout << endl;
		cout << "Time [total] = " << runtime_seconds.count() << "s\n";

		print_line_star(60);


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
