// -----------------------------------------------------------------------
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <omp.h>

#include "cntr/cntr.hpp"

#include "Holstein_impurity_decl.hpp"
// -----------------------------------------------------------------------
#define CFUNC cntr::function<double>
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define CINTEG integration::I<double>
// -----------------------------------------------------------------------
namespace Hols{
    
    ////////////////////////////////////////////// 
    // self-consistent Migdal for normal states //
    //////////////////////////////////////////////
    
    // << Matsubara and Timestep >> ---------------------------------------
    void Sigma_Mig(int tstp, GREEN &G, GREEN &Sigma, GREEN &D0, GREEN &D,
                      GREEN &Pi, GREEN &D0_Pi, GREEN &Pi_D0, CFUNC &g_el_ph, double beta,
                      double h, int SolverOrder, int MAT_METHOD){
        assert(G.size1()==1);
        int Norb=G.size1();//Nambu index=1
        int Ntau=G.ntau();
        
        GREEN_TSTP gGg(tstp,Ntau,Norb,FERMION);
        G.get_timestep(tstp,gGg);//copy time step from G
        gGg.right_multiply(tstp,g_el_ph);
        gGg.left_multiply(tstp,g_el_ph);
        
        //Get Pi(t,t') = -i2g^2 G(t,t')G(t',t)
        
        GREEN_TSTP Pol(tstp,Ntau,1,BOSON);//=i gGg(t,t') G(t',t)
        
        Pi.set_timestep_zero(tstp);
        
        Bubble1(tstp,Pol,0,0,gGg,0,0,G,0,0);
        Pi.incr_timestep(tstp,Pol,-2.0);
        
        //Solve the dyson equiation for D=D_0+D_0*Pi*D for the time step
        step_D(tstp, D0, D, Pi, D0_Pi, Pi_D0, beta, h, SolverOrder, MAT_METHOD);
        
        //Get Sig(t,t')=ig^2 D(t,t') G(t,t') 
        
        Bubble2(tstp,Sigma,0,0,D,0,0,gGg,0,0);
    }
    
    // << Bootstrapping >> ----------------------------------------------
    void Sigma_Mig(GREEN &G, GREEN &Sigma, GREEN &D0, GREEN &D, GREEN &Pi,
                      GREEN &D0_Pi, GREEN &Pi_D0, CFUNC &g_el_ph, double beta, double h,
                      int SolverOrder){
        // this is the starting version of the routine
        assert(G.size1()==1);
        int Norb=G.size1();//Nambu index=1
        int Ntau=G.ntau();
        
        GREEN gGg(SolverOrder,Ntau,Norb,FERMION);
        GREEN Pol(SolverOrder,Ntau,1,BOSON); //=i gGg(t,t') G(t',t)
        
        for(int n=0; n<=SolverOrder; n++){
            gGg.set_timestep(n,G);//copy time step from G
            gGg.right_multiply(n,g_el_ph);
            gGg.left_multiply(n,g_el_ph);
            
            //Get Pi(t,t') = -i2g^2 G(t,t')G(t',t)
            Pi.set_timestep_zero(n);
            
            Bubble1(n,Pol,0,0,gGg,0,0,G,0,0);
            Pi.incr_timestep(n,Pol,-2.0);
        }
        
        //Solve the dyson equiation for D=D_0+D_0*Pi*D for first time steps
        start_D(D0, D, Pi, D0_Pi, Pi_D0, beta, h, SolverOrder);
        
        //Get Sig(t,t')=ig^2 D(t,t')G(t,t')
        for(int n=0; n<=SolverOrder; n++) Bubble2(n,Sigma,0,0,D,0,0,gGg,0,0);
    }
    
    ///////////////////////////////////////////////////////
    // self-consistent Migdal for superconducting states //
    ///////////////////////////////////////////////////////
    
    // << Matsubara and Timestep >> ------------------------------------
    void Sigma_Mig_sc(int tstp, GREEN &G, GREEN &Sigma, GREEN &D0, GREEN &D,
      GREEN &Pi, GREEN &D0_Pi, GREEN &Pi_D0, CFUNC &g_el_ph, double beta,
      double h, int SolverOrder, int MAT_METHOD){
        assert(G.size1()==2);
        int Norb=G.size1();//Nambu index=2
        int Ntau=G.ntau();

        GREEN_TSTP gGg(tstp,Ntau,Norb,FERMION);
        G.get_timestep(tstp,gGg);//copy time step from G
        gGg.right_multiply(g_el_ph);
        gGg.left_multiply(g_el_ph);

        //Get Pi(t,t') = -ig^2 tr[sig_3 G(t,t') sig_3 G(t',t)]

        GREEN_TSTP Pol(tstp,Ntau,1,BOSON);//=i gGg_{ab}(t,t') G_{ba}(t',t)

        Pi.set_timestep_zero(tstp);

        for(int a_=0;a_<Norb;a_++){
            for(int b_=0;b_<Norb;b_++){
                Bubble1(tstp,Pol,0,0,gGg,a_,b_,G,a_,b_);
                Pi.incr_timestep(tstp,Pol,-1.0);
            }
        }

        //Solve the dyson equiation for D=D_0+D_0*Pi*D for the time step
        step_D(tstp, D0, D, Pi, D0_Pi, Pi_D0, beta, h, SolverOrder, MAT_METHOD);

        //Get Sig(t,t')=ig^2 D(t,t') sig_3 G(t,t') sig_3
        for(int a_=0;a_<Norb;a_++){
            for(int b_=0;b_<Norb;b_++){
                Bubble2(tstp,Sigma,a_,b_,D,0,0,gGg,a_,b_);
            }
        }

    }
    
    // << Bootstrapping >> ----------------------------------------------
    void Sigma_Mig_sc(GREEN &G, GREEN &Sigma, GREEN &D0, GREEN &D, GREEN &Pi,
      GREEN &D0_Pi, GREEN &Pi_D0, CFUNC &g_el_ph, double beta, double h,
      int SolverOrder){
        // this is the starting version of the routine
        assert(G.size1()==2);
        int Norb=G.size1();//Nambu index=2
        int Ntau=G.ntau();

        GREEN gGg(SolverOrder,Ntau,Norb,FERMION);
        GREEN Pol(SolverOrder,Ntau,1,BOSON); //=i gGg_{ab}(t,t') G_{ba}(t',t)

        for(int n=0; n<=SolverOrder; n++){
          gGg.set_timestep(n,G);//copy time step from G
          gGg.right_multiply(n,g_el_ph);
          gGg.left_multiply(n,g_el_ph);

          //Get Pi(t,t') = -ig^2 tr[sig_3 G(t,t') sig_3 G(t',t)]
          Pi.set_timestep_zero(n);

          for(int a_=0;a_<Norb;a_++){
              for(int b_=0;b_<Norb;b_++){
                  Bubble1(n,Pol,0,0,gGg,a_,b_,G,a_,b_);
                  Pi.incr_timestep(n,Pol,-1.0);
              }
          }
        }

        //Solve the dyson equiation for D=D_0+D_0*Pi*D for first time steps
        start_D(D0, D, Pi, D0_Pi, Pi_D0, beta, h, SolverOrder);

        //Get Sig(t,t')=ig^2 D(t,t') sig_3 G(t,t') sig_3

        for(int n=0; n<=SolverOrder; n++){
          for(int a_=0;a_<Norb;a_++){
              for(int b_=0;b_<Norb;b_++){
                  Bubble2(n,Sigma,a_,b_,D,0,0,gGg,a_,b_);
              }
          }
        }
    }
    
    
    ///////////////////////////////////////////////////////
    // Unrenormalized Migdal for normal states //
    ///////////////////////////////////////////////////////
    
    // <<Matsubara,timestep>>-----------------------------------------
    void Sigma_uMig(int tstp, GREEN &G, GREEN &D0, CFUNC &g_el_ph, GREEN &Sigma){
        
        assert(G.size1()==1);
        int Norb=G.size1();
        int Ntau=G.ntau();
        
        GREEN_TSTP gGg(tstp,Ntau,Norb,-1);
        G.get_timestep(tstp,gGg);//copy time step from G
        gGg.right_multiply(tstp,g_el_ph);
        gGg.left_multiply(tstp,g_el_ph);
        
        //Get Sig(t,t')=ig^2 D_0(t,t') G(t,t') 

        Bubble2(tstp,Sigma,0,0,D0,0,0,gGg,0,0);

    }
    // <<bootstrap>>-------------------------------------------------
    void Sigma_uMig(GREEN &G, GREEN &D0, CFUNC &g_el_ph, GREEN &Sigma,int SolverOrder){
        
        assert(G.size1()==1);
        int Norb=G.size1();
        int Ntau=G.ntau();
        
        for(int n=0 ; n<=SolverOrder ;n++){
            GREEN_TSTP gGg(n,Ntau,Norb,-1);
            G.get_timestep(n,gGg);//copy time step from G
            gGg.right_multiply(n,g_el_ph);
            gGg.left_multiply(n,g_el_ph);
        
            //Get Sig(t,t')=ig^2 D_0(t,t') G(t,t') 
            Bubble2(n,Sigma,0,0,D0,0,0,gGg,0,0);
        }
        
    }
    ///////////////////////////////////////////////////////
    // Unrenormalized Migdal for superconducting states //
    ///////////////////////////////////////////////////////
    
    // <<Matsubara,timestep>>-----------------------------------------
    void Sigma_uMig_sc(int tstp, GREEN &G, GREEN &D0, CFUNC &g_el_ph, GREEN &Sigma){
        
        assert(G.size1()==2);
        int Norb=G.size1();//Nambu index=2
        int Ntau=G.ntau();
        
        GREEN_TSTP gGg(tstp,Ntau,Norb,-1);
        G.get_timestep(tstp,gGg);//copy time step from G
        gGg.right_multiply(tstp,g_el_ph);
        gGg.left_multiply(tstp,g_el_ph);
        
        //Get Sig(t,t')=ig^2 D_0(t,t') sig_3 G(t,t') sig_3
        for(int a_=0;a_<Norb;a_++){
            for(int b_=0;b_<Norb;b_++){
                Bubble2(tstp,Sigma,a_,b_,D0,0,0,gGg,a_,b_);
            }
        }
    }
    
    // <<bootstrap>>-------------------------------------------------
    void Sigma_uMig_sc(GREEN &G, GREEN &D0, CFUNC &g_el_ph, GREEN &Sigma,int SolverOrder){
        
        assert(G.size1()==2);
        int Norb=G.size1();//Nambu index=2
        int Ntau=G.ntau();
        
        for(int n=0 ; n<=SolverOrder ; n++){
            GREEN_TSTP gGg(n,Ntau,Norb,-1);
            G.get_timestep(n,gGg);//copy time step from G
            gGg.right_multiply(n,g_el_ph);
            gGg.left_multiply(n,g_el_ph);
        
            //Get Sig(t,t')=ig^2 D_0(t,t') sig_3 G(t,t') sig_3
            for(int a_=0;a_<Norb;a_++){
                for(int b_=0;b_<Norb;b_++){
                    Bubble2(n,Sigma,a_,b_,D0,0,0,gGg,a_,b_);
                }
            }
        }
    }
        
    ////////////////////////////////////////////////
    // Dyson equation for phonon Green's function //
    ////////////////////////////////////////////////
    // <<Matsubara,timestep>>-------------------------------------
    void step_D(int tstp,GREEN &D0, GREEN &D, GREEN &Pi, GREEN &D0_Pi, GREEN &Pi_D0,
      double beta, double h, int SolverOrder, int MAT_METHOD){
        //set D0_Pi=-D0*Pi
        cntr::convolution_timestep(tstp,D0_Pi,D0,Pi,beta,h,SolverOrder);
        D0_Pi.smul(tstp,-1.0);

        cntr::convolution_timestep(tstp,Pi_D0,Pi,D0,beta,h,SolverOrder);
        Pi_D0.smul(tstp,-1.0);


        //solve [1-D0*Pi]*D=[1+D0_Pi]*D=D0
        if(tstp==-1){
            cntr::vie2_mat(D,D0_Pi,Pi_D0,D0,beta,SolverOrder,MAT_METHOD);
        }else {
            cntr::vie2_timestep(tstp,D,D0_Pi,Pi_D0,D0,beta,h,SolverOrder);
        }

    }
    // <<bootstrap>>----------------------------------------------
    void start_D(GREEN &D0, GREEN &D, GREEN &Pi, GREEN &D0_Pi, GREEN &Pi_D0,
      double beta, double h, int SolverOrder){
        //set D0_Pi=-D0*Pi
        for(int n=0; n<=SolverOrder; n++){
          cntr::convolution_timestep(n,D0_Pi,D0,Pi,beta,h,SolverOrder);
          D0_Pi.smul(n,-1.0);

          cntr::convolution_timestep(n,Pi_D0,Pi,D0,beta,h,SolverOrder);
          Pi_D0.smul(n,-1.0);
        }

        //solve [1-D0*Pi]*D=[1+D0_Pi]*D=D0
        cntr::vie2_start(D,D0_Pi,Pi_D0,D0,beta,h,SolverOrder);
    }
    
    ////////////////////////////
    // Phonon displacement    //
    ////////////////////////////
    // <<Matsubara,timestep>>---------------------------------
    void get_phonon_displace(int tstp,CFUNC &X_ph, CFUNC &n_tot,CFUNC &g_el_ph,GREEN &D0, double w0, int SolverOrder, double dt){
        
        assert(X_ph.size1()==1);
        assert(X_ph.size2()==1);
        assert(n_tot.size1()==1);
        assert(n_tot.size2()==1);
        assert(g_el_ph.size1()==1);
        assert(g_el_ph.size2()==1);
        
        cdmatrix n_g_0(1,1),X_ph_0(1,1);
        cdmatrix n_g_tmp(1,1),X_ph_tmp(1,1);
        cdmatrix D0_ret_tmp(1,1);
        
        CFUNC n_g_el_ph(n_tot);
        n_g_el_ph.left_multiply(g_el_ph,1.0);
        
        n_g_el_ph.get_value(-1,n_g_0);
        
        X_ph_0 =-2.0*n_g_0/w0;
        
        if(tstp==-1) X_ph_tmp = X_ph_0;
        
        else{
            vector<complex<double> > D0_gn;
            D0_gn.assign(tstp+1,0.0);
            integration::Integrator<double> Integ(SolverOrder); //[Check <<= should this be complex double?]
            
            for(int it =0 ; it<=tstp ; it++){
                n_g_el_ph.get_value(it,n_g_tmp);
                D0.get_ret(tstp,it,D0_ret_tmp);
                
                D0_gn[it] = D0_ret_tmp(0,0)*(n_g_tmp(0,0)-n_g_0(0,0));
            }
            
            X_ph_tmp(0,0) = Integ.integrate(D0_gn,tstp)*dt; //[check]
            
            X_ph_tmp += X_ph_0;
        }
        
        X_ph.set_value(tstp,X_ph_tmp);
    }
    
    // <<bootstrap>>---------------------------------
    void get_phonon_displace(CFUNC &X_ph, CFUNC &n_tot,CFUNC &g_el_ph,GREEN &D0, double w0, int SolverOrder, double dt){
        
        assert(X_ph.size1()==1);
        assert(X_ph.size2()==1);
        assert(n_tot.size1()==1);
        assert(n_tot.size2()==1);
        assert(g_el_ph.size1()==1);
        assert(g_el_ph.size2()==1);
        
        cdmatrix n_g_0(1,1),X_ph_0(1,1);
        cdmatrix n_g_tmp(1,1),X_ph_tmp(1,1);
        cdmatrix D0_ret_tmp(1,1);
        
        CFUNC n_g_el_ph(n_tot);
        n_g_el_ph.left_multiply(g_el_ph,1.0);
        
        n_g_el_ph.get_value(-1,n_g_0);
        
        X_ph_0 =-2.0*n_g_0/w0;

        integration::Integrator<double> Integ(SolverOrder); 
    
        for(int n=0 ; n<=SolverOrder ; n++){
            vector<complex<double> > D0_gn;
            D0_gn.assign(SolverOrder+1,0.0);

            //\int^tk_0 dt' D^R_0(t_k-t') (gn(t'+tn-tk) -gn(0)) 
            // with gn(t') = gn(0) for t'<0
            for(int it =0 ; it<=SolverOrder ; it++){
                
                if(it+n-SolverOrder>=0) n_g_el_ph.get_value(it+n-SolverOrder,n_g_tmp);
                else n_g_tmp = n_g_0;
                
                D0.get_ret(SolverOrder,it,D0_ret_tmp);
                
                D0_gn[it] = D0_ret_tmp(0,0)*(n_g_tmp(0,0)-n_g_0(0,0));
            }
            
            X_ph_tmp(0,0) = Integ.integrate(D0_gn,SolverOrder)*dt;//[check]
            
            X_ph_tmp += X_ph_0;
            
            X_ph.set_value(n,X_ph_tmp);
        }
    }
}
// -----------------------------------------------------------------------
