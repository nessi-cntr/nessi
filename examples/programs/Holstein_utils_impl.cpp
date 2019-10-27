// -----------------------------------------------------------------------
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include <omp.h>

#include "cntr/cntr.hpp"
#include "Holstein_utils_decl.hpp"
// -----------------------------------------------------------------------
#define CFUNC cntr::function<double>
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define CPLX complex<double>
// -----------------------------------------------------------------------
namespace Hols{
    
    void evaluate_phonon_energy_qu(CFUNC & Eph_t, GREEN & D, GREEN & Pi,int SolverOrder, double beta, double h, double w0){
        
        CPLX ii = CPLX(0,1.0);
        
        int nt = D.nt();
        int ntau = D.ntau();
        double dtau = beta/(double)ntau;
        
        cdmatrix Mat_tmp(1,1);
        cdmatrix XX(1,1),PP(1,1);
        GREEN D0_d1, D0_d2, D_d1, D_d2, D0_d1_Pi, Pi_D0_d2;
        
        D0_d1 = GREEN(nt,ntau,1,BOSON);
        D0_d2 = GREEN(nt,ntau,1,BOSON);
        D_d1 = GREEN(nt,ntau,1,BOSON);
        D_d2 = GREEN(nt,ntau,1,BOSON);
        D0_d1_Pi = GREEN(nt,ntau,1,BOSON);
        Pi_D0_d2 = GREEN(nt,ntau,1,BOSON);
        
        //--------------------
        //set D0_d1 and D0_d2
        //--------------------
        double fb;
        vector<CPLX> expp;
        expp.resize(nt+1);
        for(int it=0 ; it<=nt ; it++) expp[it] = exp(-ii*(double)it*h*w0);
        fb = bose(beta,w0);
        
        //matsubara and left-mixing
        for(int itau=0 ; itau<=ntau ; itau++){
            double tau = double(itau)*dtau;
            CPLX coeff1 = exp(-tau*w0);
            CPLX coeff2 = exp(tau*w0);
            CPLX d0_d1_mat = -(1.0 + fb)*coeff1 + fb*coeff2;

            D0_d1.set_mat(itau, d0_d1_mat);
            D0_d2.set_mat(itau, -d0_d1_mat);
            
            for(int it=0 ; it<=nt ; it++){
                CPLX d0_d1_tv = -fb*expp[it]*coeff2 + (1.0+fb)*conj(expp[it])*coeff1;
                D0_d1.set_tv(it,itau,d0_d1_tv);
                D0_d2.set_tv(it,itau,-d0_d1_tv);
            }
        }
        
        //retarded and lesser
        for(int it1=0 ; it1<=nt ; it1++){
            for(int it2=0 ; it2<=it1 ; it2++){
                CPLX d0_d1_ret = -(expp[it1-it2]+conj(expp[it1-it2]));
                CPLX d0_d1_les = -fb*conj(expp[it1-it2]) + (1.0+fb)*expp[it1-it2];
                D0_d1.set_ret(it1,it2,d0_d1_ret);
                D0_d2.set_ret(it1,it2,-d0_d1_ret);
                D0_d1.set_les(it2,it1,d0_d1_les);
                D0_d2.set_les(it2,it1,-d0_d1_les);
            }
        }
        
        //--------------------
        //evaluate D0_d1_Pi
        //--------------------
        cntr::convolution(D0_d1_Pi, D0_d1, D0_d2, Pi, Pi, beta, h, SolverOrder);
        
        //--------------------
        //evaluate Pi_D0_d2
        //--------------------
        cntr::convolution(Pi_D0_d2, Pi, Pi, D0_d2, D0_d1, beta, h, SolverOrder);
        
        //--------------------
        // evalaute Dd1
        //--------------------
        cntr::convolution(D_d1, D0_d1_Pi, Pi_D0_d2, D, D, beta, h, SolverOrder);
        for(int tstp=-1; tstp <= nt; tstp++) D_d1.incr_timestep(tstp,D0_d1,1.0);
        
        //--------------------
        // evaluate Dd2
        //--------------------
        cntr::convolution(D_d2, D, D, Pi_D0_d2, D0_d1_Pi, beta, h, SolverOrder);
        for(int tstp=-1; tstp <= nt; tstp++) D_d2.incr_timestep(tstp,D0_d2,1.0);
        
        for(int tstp=-1; tstp <= nt; tstp++){
            
            if(tstp==-1){
                D.get_mat(0,XX);
                XX(0,0)*=-1.0;
            }
            else{
                D.get_les(tstp,tstp,XX);
                XX(0,0)*=CPLX(0.0,1.0);
            }
            //PP_int = i*[D0_d1 * Pi * D_d2]^<(t,t)
            cdmatrix PP_(1,1);
            cntr::convolution_density_matrix(tstp,PP,D0_d1_Pi,Pi_D0_d2,D_d2,D_d1, beta, h, SolverOrder);
            PP *= -1.0;
            //PP_corr = PP_int + i*[D0_d1d2]^<(t,t)
            PP(0,0) += 1.0 + 2.0*cntr::bose(beta,w0);
            
            
            Mat_tmp = (XX+PP)*w0/4.0;
            
            Eph_t.set_value(tstp,Mat_tmp);
            
        }
        
    }
    
    void evaluate_phonon_energy_cl(CFUNC & Eph_t, CFUNC &X_ph, CFUNC &P_ph ,double w0){
        
        cdmatrix Eph_tmp(1,1),Xph_tmp(1,1),Pph_tmp(1,1);
        int Nt = Eph_t.nt();
        for(int tstp=-1 ; tstp<=Nt ; tstp++){
            X_ph.get_value(tstp,Xph_tmp);
            P_ph.get_value(tstp,Pph_tmp);
            
            Eph_tmp(0,0) = 0.25*w0*(Xph_tmp(0,0)*Xph_tmp(0,0)+Pph_tmp(0,0)*Pph_tmp(0,0));
            
            Eph_t.set_value(tstp,Eph_tmp);
        }
    }
    
    void get_phonon_momentum(CFUNC &P_ph, CFUNC &n_tot,CFUNC &g_el_ph, GREEN &D0, double w0, int SolverOrder, double h){
        
        CPLX ii = CPLX(0,1.0);
        int tstp;
        cdmatrix Mat_tmp(1,1);
        int Nt = D0.nt();
                
        CFUNC n_g_el_ph(n_tot);
        n_g_el_ph.left_multiply(g_el_ph,1.0);
        
        cdmatrix n_g_0(1,1),n_g_tmp(1,1);
        n_g_el_ph.get_value(-1,n_g_0);
        
        //Rtarded part of D0_d1 = \partial_t D_0^R/\omega_0
        vector<complex<double> > D0_d1_ret;
        D0_d1_ret.assign(Nt+1,0.0);
        for(int it=0 ; it<=Nt ; it++) D0_d1_ret[it] = -(exp(-ii*(double)it*h*w0)+exp(ii*(double)it*h*w0));

        //--------------------
        // Matsubara
        //--------------------
        tstp = -1;
        Mat_tmp(0,0) = 0.0;
        P_ph.set_value(tstp,Mat_tmp);

        //--------------------
        // Bootstrap
        //--------------------
        for(tstp=0 ; tstp<=SolverOrder ; tstp++){
            vector<complex<double> > D0_d1_gn;
            D0_d1_gn.assign(SolverOrder+1,0.0);
            integration::Integrator<double> Integ(SolverOrder); 
            
            //\int^tk_0 h' D^R_0_d1(t_k-t') (gn(t'+tn-tk) -gn(0)) 
            // with gn(t') = gn(0) for t'<0
            for(int it =0 ; it<=SolverOrder ; it++){
                
                if(it+tstp-SolverOrder>=0) n_g_el_ph.get_value(it+tstp-SolverOrder,n_g_tmp);
                else n_g_tmp = n_g_0;
                                
                D0_d1_gn[it] = D0_d1_ret[SolverOrder-it]*(n_g_tmp(0,0)-n_g_0(0,0));
            }
            
            Mat_tmp(0,0) = Integ.integrate(D0_d1_gn,SolverOrder)*h;//[check]
                        
            P_ph.set_value(tstp,Mat_tmp);
        }
        
        //--------------------
        // timestep
        //--------------------
        for(tstp=SolverOrder+1 ; tstp<=Nt ; tstp++){
            
            vector<complex<double> > D0_d1_gn;
            D0_d1_gn.assign(tstp+1,0.0);
            integration::Integrator<double> Integ(SolverOrder); //[Check <<= should this be complex double?]
            
            for(int it =0 ; it<=tstp ; it++){
                n_g_el_ph.get_value(it,n_g_tmp);
                D0_d1_gn[it] = D0_d1_ret[tstp-it]*(n_g_tmp(0,0)-n_g_0(0,0));
            }
            
            Mat_tmp(0,0) = Integ.integrate(D0_d1_gn,tstp)*h; //[check]
            
            P_ph.set_value(tstp,Mat_tmp);
        }
                
    }
}
// -----------------------------------------------------------------------
