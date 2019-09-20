// -----------------------------------------------------------------------
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

#include "cntr/cntr.hpp"

#include "hubbard_chain_selfen_decl.hpp"
// -----------------------------------------------------------------------
namespace hubb{
  // -----------------------------------------------------------------------
  void Polarization(int tstp, GREEN &G, GREEN &Pol){
    int nsites=G.size1();

    for(int i=0; i<nsites; i++){
      for(int j=0; j<nsites; j++){
        cntr::Bubble1(tstp,Pol,i,j,G,i,j,G,i,j);
      }
    }
    Pol.smul(tstp,-1.0);
  }
  // -----------------------------------------------------------------------
  void Polarization(int tstp, GREEN &G, GREEN_TSTP &Pol){
    int nsites=G.size1();

    for(int i=0; i<nsites; i++){
      for(int j=0; j<nsites; j++){
        cntr::Bubble1(tstp,Pol,i,j,G,i,j,G,i,j);
      }
    }
    Pol.smul(-1.0);
  }
  // -----------------------------------------------------------------------
  void Sigma_MF(int tstp, GREEN &G, CFUNC &U, cdmatrix &Sigma_mf){
    int nst=G.size1();
    cdmatrix rho(nst,nst), Ut(nst,nst);

    G.density_matrix(tstp, rho);
    U.get_value(tstp,Ut);

    Sigma_mf.setZero();
    for(int n1=0; n1<nst; n1++){
      Sigma_mf(n1,n1) = Ut(n1,n1)*rho(n1,n1);
    }

  }
  // -----------------------------------------------------------------------
  void Ham_MF(int tstp, GREEN &G, CFUNC &U, cdmatrix &h0, CFUNC &hmf){
    int nsites=G.size1();
    cdmatrix rho(nsites,nsites), Ut(nsites,nsites);
    cdmatrix hmf_val(nsites,nsites);

    G.density_matrix(tstp, rho);
    U.get_value(tstp,Ut);
    hmf_val = h0;
    for(int i=0; i<nsites; i++) hmf_val(i,i) += Ut(i,i)*rho(i,i);
    hmf.set_value(tstp, hmf_val);

  }
  // -----------------------------------------------------------------------
  void Sigma_2B(int tstp, GREEN &G, CFUNC &U, GREEN &Sigma){
    int nst=G.size1();
    int ntau=G.ntau();
    GREEN_TSTP Pol(tstp,ntau,nst,BOSON);

    Polarization(tstp, G, Pol);

    Pol.right_multiply(tstp, U);
    Pol.left_multiply(tstp, U);

    for(int i=0; i<nst; i++){
      for(int j=0; j<nst; j++){
        cntr::Bubble2(tstp,Sigma,i,j,G,i,j,Pol,i,j);
      }
    }

  }
  // -----------------------------------------------------------------------
  void Sigma_GW(int tstp, GREEN &G, CFUNC &U, GREEN &Chi, GREEN &Sigma){
    int nsites=G.size1();
    int ntau=G.ntau();
    GREEN_TSTP deltaW(tstp,ntau,nsites,BOSON);

    Chi.get_timestep(tstp,deltaW);
    deltaW.left_multiply(tstp, U);
    deltaW.right_multiply(tstp, U);

    for(int i=0; i<nsites; i++){
      for(int j=0; j<nsites; j++){
        cntr::Bubble2(tstp,Sigma,i,j,G,i,j,deltaW,i,j);
      }
    }
  }
  // -----------------------------------------------------------------------
  void Sigma_TPP(int tstp, GREEN &G, CFUNC &U, GREEN &TPP, GREEN &Sigma){
    int nsites=G.size1();
    int ntau=G.ntau();
    GREEN_TSTP UTU(tstp,ntau,nsites,BOSON);

    TPP.get_timestep(tstp,UTU);
    UTU.left_multiply(U);
    UTU.right_multiply(U);

    for(int i=0; i<nsites; i++){
      for(int j=0; j<nsites; j++){
        cntr::Bubble1(tstp,Sigma,i,j,UTU,i,j,G,i,j);
      }
    }
  }
  // -----------------------------------------------------------------------
  void GenChi(double dt, double beta, GREEN &Pol, CFUNC &U,
    GREEN &PxU, GREEN &UxP, GREEN &Chi, int SolveOrder){

      for(int n = 0; n <= SolveOrder; n++){
        PxU.set_timestep(n, Pol);
        UxP.set_timestep(n, Pol);
        PxU.right_multiply(n, U);
        UxP.left_multiply(n, U);
        PxU.smul(n,-1.0);
        UxP.smul(n,-1.0);
      }

      cntr::vie2_start(Chi,PxU,UxP,Pol,beta,dt,SolveOrder);

  }
  // -----------------------------------------------------------------------
  void GenChi(int tstp, double dt, double beta, GREEN &Pol, CFUNC &U, GREEN &PxU,
      GREEN &UxP, GREEN &Chi, int SolveOrder){

      PxU.set_timestep(tstp, Pol);
      UxP.set_timestep(tstp, Pol);
      PxU.right_multiply(tstp, U);
      UxP.left_multiply(tstp, U);
      PxU.smul(tstp,-1.0);
      UxP.smul(tstp,-1.0);

      if(tstp==-1){
        cntr::vie2_mat(Chi,PxU,UxP,Pol,beta,SolveOrder);
      } else{
        cntr::vie2_timestep(tstp,Chi,PxU,UxP,Pol,beta,dt,SolveOrder);
      }
  }
  // -----------------------------------------------------------------------
  void GenTPP(double dt, double beta, GREEN &G, GREEN &Phi, CFUNC &U, GREEN &UxPhi,
      GREEN &PhixU, GREEN &TPP, int SolveOrder){
      int nst=G.size1();

      for(int tstp=0; tstp<=SolveOrder; tstp++){
        for(int n1=0; n1<nst; n1++){
          for(int n2=0; n2<nst; n2++){
            cntr::Bubble2(tstp,Phi,n1,n2,G,n1,n2,G,n1,n2);
          }
        }
        Phi.smul(tstp, -1.0);

        PhixU.set_timestep(tstp, Phi);
        UxPhi.set_timestep(tstp, Phi);

        PhixU.right_multiply(tstp, U);
        UxPhi.left_multiply(tstp, U);
      }

      cntr::vie2_start(TPP,PhixU,UxPhi,Phi,beta,dt,SolveOrder);
  }
    // -----------------------------------------------------------------------
    void GenTPP(int tstp, double dt, double beta, GREEN &G, GREEN &Phi, CFUNC &U,
      GREEN &UxPhi, GREEN &PhixU, GREEN &TPP, int SolveOrder){
      int nst=G.size1();

      for(int n1=0; n1<nst; n1++){
        for(int n2=0; n2<nst; n2++){
          cntr::Bubble2(tstp,Phi,n1,n2,G,n1,n2,G,n1,n2);
        }
      }
      Phi.smul(tstp, -1.0);

      PhixU.set_timestep(tstp, Phi);
      UxPhi.set_timestep(tstp, Phi);

      PhixU.right_multiply(tstp, U);
      UxPhi.left_multiply(tstp, U);

      if(tstp == -1) {
        cntr::vie2_mat(TPP,PhixU,UxPhi,Phi,beta,SolveOrder);
      } else {
        cntr::vie2_timestep(tstp,TPP,PhixU,UxPhi,Phi,beta,dt,SolveOrder);
      }
    }
    // -----------------------------------------------------------------------
  }
  // -----------------------------------------------------------------------
