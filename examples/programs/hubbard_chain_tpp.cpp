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
#include "hubbard_chain_selfen_decl.hpp"

using namespace std;
//------------------------------------------------------------------------------

// -----------------------------------------------------------------------
double KineticEnergy(int tstp, cdmatrix &eps0, GREEN &G){
  int nst = G.size1();
  cdmatrix rho(nst,nst);
  G.density_matrix(tstp, rho);
  return ((eps0*rho).trace()).real();
}
// -----------------------------------------------------------------------
double InteractionEnergy(int tstp,cdmatrix &eps0, CFUNC &eps_mf, GREEN &G, GREEN &Sigma,
  double beta, double h, int SolveOrder){
    int nst = G.size1();
    double Emf,Ecorr;
    cdmatrix rho(nst,nst),vmf(nst,nst);

    eps_mf.get_value(tstp,vmf);
    vmf = vmf - eps0;

    G.density_matrix(tstp, rho);
    Emf = 0.5*((vmf*rho).trace()).real();

    Ecorr = cntr::correlation_energy(tstp, G, Sigma, beta, h, SolveOrder);
    return Emf + Ecorr;
  }
  // -----------------------------------------------------------------------

  //==============================================================================
  //         main program
  //==============================================================================
  int main(int argc,char *argv[]){
    int SolveOrder = MAX_SOLVE_ORDER;
    //..................................................
    //                input
    //..................................................
    int Nt,Ntau,MatsMaxIter,CorrectorSteps,Nsites;
    int BootstrapMaxIter;
    double HoppingT,HubbardU,beta,h,MuChem,MatsMaxErr,BootstrapMaxErr;
    int RampSite;
    double RampW0;
    char* flin;
    char* flsave;
    char* flout;
    //..................................................
    //                internal
    //..................................................
    int tstp;
    double npart,err;
    cdmatrix eps0,densm;
    GREEN G,Sigma,Phi,PHIxU,UxPHI,TPP,gtemp;
    CFUNC Ut,eps_mf;
    //..................................................
    //                timer
    //..................................................
    chrono::time_point<std::chrono::system_clock> start, end, start_tot, end_tot;

    print_line_star(60);
    cout << "   Test program: Hubbard chain in T-matrix approximation" << endl;
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
        find_param(flin,"__Nsites=",Nsites);
        find_param(flin,"__HoppingT=",HoppingT);
        find_param(flin,"__HubbardU=",HubbardU);
        find_param(flin,"__MuChem=",MuChem);
        find_param(flin,"__beta=",beta);

        // ramp parameters
        find_param(flin,"__RampSite=",RampSite);
        find_param(flin,"__RampW0=",RampW0);

        // solver parameters
        find_param(flin,"__Nt=",Nt);
        find_param(flin,"__Ntau=",Ntau);
        find_param(flin,"__h=",h);
        find_param(flin,"__MatsMaxIter=",MatsMaxIter);
        find_param(flin,"__MatsMaxErr=",MatsMaxErr);
        find_param(flin,"__BootstrapMaxIter=",BootstrapMaxIter);
        find_param(flin,"__BootstrapMaxErr=",BootstrapMaxErr);
        find_param(flin,"__CorrectorSteps=",CorrectorSteps);

        // output file prefix
        flout=argv[2];

      }
      //============================================================================
      //                   MEMORY REQUIREMENTS
      //============================================================================
      {
        print_line_minus(50);
        cout << "     Memory requirements" << endl;
        print_line_minus(50);

        const size_t size_MB=1024*1024;
        size_t mem_1time=0, mem_2time=0;

        mem_1time += cntr::mem_function<double>(Nt,Nsites); // Ut
        mem_1time += cntr::mem_function<double>(Nt,Nsites); // eps_mf

        mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Nsites); // G
        mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Nsites); // Sigma
        mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Nsites); // Phi
        mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Nsites); // TPP
        mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Nsites); // UxPHI
        mem_2time += cntr::mem_herm_matrix<double>(Nt,Ntau,Nsites); // PHIxU

        // convert to MB
        mem_1time = ceil(mem_1time/(double)size_MB);
        mem_2time = ceil(mem_2time/(double)size_MB);

        cout << "Hamiltonian : " << mem_1time << " MB" << endl;
        cout << "two-time quantities : " << mem_2time << " MB" << endl;

        print_line_minus(50);
        cout << "\n\n";
      }
      //============================================================================
      //               (IV) INITIALIZE GREEN'S FUNCTIONS
      //============================================================================
      {
        G = GREEN(Nt,Ntau,Nsites,FERMION);
        Phi = GREEN(Nt,Ntau,Nsites,BOSON);
        UxPHI = GREEN(Nt,Ntau,Nsites,BOSON);
        PHIxU = GREEN(Nt,Ntau,Nsites,BOSON);
        TPP = GREEN(Nt,Ntau,Nsites,BOSON);
        Sigma = GREEN(Nt,Ntau,Nsites,FERMION);

        // --- mean field ---
        eps0.resize(Nsites,Nsites);
        eps0.setZero();
        eps_mf = CFUNC(Nt,Nsites);
        for (int i = 0; i < Nsites-1; i++){
          eps0(i,i+1) = -HoppingT;
          eps0(i+1,i) = -HoppingT;
        }

        Ut = CFUNC(Nt,Nsites);
        Ut.set_constant(HubbardU*MatrixXcd::Identity(Nsites,Nsites));

        cntr::green_from_H(G, MuChem, eps0, beta, h);

        // shift the chemical potential according to the expected average
        // occupation
        G.density_matrix(-1,densm);
        npart = (densm.trace()).real();
        MuChem = MuChem + HubbardU * npart / Nsites;

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

        gtemp = GREEN(SolveOrder,Ntau,Nsites,FERMION);
        gtemp.set_timestep(tstp,G);

        for(int iter=0;iter<=MatsMaxIter;iter++){
          // update mean field
          hubb::Ham_MF(tstp, G, Ut, eps0, eps_mf);

          // update T-matrix
          hubb::GenTPP(tstp, h, beta, G, Phi, Ut, UxPHI, PHIxU, TPP, SolveOrder);

          // update self-energy
          hubb::Sigma_TPP(tstp, G, Ut, TPP, Sigma);

          // solve Dyson equation
          cntr::dyson_mat(G, MuChem, eps_mf, Sigma, beta, SolveOrder);

          // compute number of particles
          G.density_matrix(-1,densm);
          npart = (densm.trace()).real();

          // self-consistency check
          err = cntr::distance_norm2(tstp,G,gtemp);
          cout << "iteration : " << iter << "  | N =  " << npart << " |  Error = " << err << endl;

          if(err<MatsMaxErr){
            matsubara_converged=true;
            break;
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
        //           BOOTSTRAPPING PHASE
        //============================================================================
        { // begin bootstrapping

          print_line_minus(50);
          cout << "     Time propagation: bootstrapping phase " << endl;
          print_line_minus(50);

          start = std::chrono::system_clock::now();

          bool bootstrap_converged=false;

          // to represent the quench, the free Hamiltonian is updated
          eps0(RampSite-1,RampSite-1) = RampW0;

          for(tstp=0; tstp<=SolveOrder; tstp++)
          gtemp.set_timestep(tstp,G);

          for (int iter = 0; iter <= BootstrapMaxIter; iter++) {
            // update mean field
            for(tstp=0; tstp<=SolveOrder; tstp++){
              hubb::Ham_MF(tstp, G, Ut, eps0, eps_mf);
            }

            // update T-matrix
            hubb::GenTPP(h, beta, G, Phi, Ut, UxPHI, PHIxU, TPP, SolveOrder);

            // update self-energy
            for(tstp=0; tstp<=SolveOrder; tstp++){
              hubb::Sigma_TPP(tstp, G, Ut, TPP, Sigma);
            }

            // solve Dyson equation
            cntr::dyson_start(G, MuChem, eps_mf, Sigma, beta, h, SolveOrder);

            // self-consistency check
            err=0.0;
            for(tstp=0; tstp<=SolveOrder; tstp++) {
              err += cntr::distance_norm2(tstp,G,gtemp);
            }
            cout << "bootstrap iteration : " << iter << "  |  Error = " << err << endl;
            if(err<BootstrapMaxErr && iter>2){
              bootstrap_converged=true;
              break;
            }
            for(tstp=0; tstp<=SolveOrder; tstp++) {
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

          for(tstp = SolveOrder+1; tstp <= Nt; tstp++){
            // Predictor: extrapolation
            cntr::extrapolate_timestep(tstp-1,G,SolveOrder);
            // Corrector
            for (int iter=0; iter < CorrectorSteps; iter++){
              // update mean field
              hubb::Ham_MF(tstp, G, Ut, eps0, eps_mf);

              // update T-matrix
              hubb::GenTPP(tstp, h, beta, G, Phi, Ut, UxPHI, PHIxU, TPP, SolveOrder);

              // update self-energy
              hubb::Sigma_TPP(tstp, G, Ut, TPP, Sigma);

              // solve Dyson equation
              cntr::dyson_timestep(tstp,G,MuChem,eps_mf,Sigma,beta,h,SolveOrder);

            }
          }

          print_line_dot(50);
          end = std::chrono::system_clock::now();
          std::chrono::duration<double> elapsed_seconds = end-start;
          cout << "Time [KBEs] = " << elapsed_seconds.count() << "s\n";

        } // end propagation loop

        { // output

          char file_occ[255];
          strcpy(file_occ,flout);
          strcat(file_occ,"_occupation.dat");
          ofstream f_occ;
          f_occ.open (file_occ);

          // compute density matrix
          cdmatrix rho(Nsites,Nsites);
          vector<double> occ(Nsites);
          for(tstp=0; tstp <= Nt; tstp++){
            G.density_matrix(tstp, rho);
            f_occ << tstp*h << "  " ;
            for(int i=0; i<Nsites; i++){
              occ[i] = (rho(i,i)).real();
              f_occ << occ[i] << "  ";
            }
            f_occ << endl;
          }
          f_occ.close();

          char file_en[255];
          strcpy(file_en,flout);
          strcat(file_en,"_energy.dat");
          ofstream f_en;
          f_en.open (file_en);

          // compute energy
          double Ekin,Epot,Etot;
          f_en << setprecision(10);
          for(tstp=0; tstp <= Nt; tstp++){
            Ekin = KineticEnergy(tstp, eps0, G);
            Epot = InteractionEnergy(tstp, eps0, eps_mf, G, Sigma, beta, h, SolveOrder);
              Etot = Ekin + Epot;
              f_en << tstp*h << "  " << Ekin << "  " << Epot << "  " << Etot << endl;
            }
            f_en.close();


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
