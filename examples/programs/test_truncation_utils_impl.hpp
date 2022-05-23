#ifndef TRUNCATION_TEST_UTILS_IMPL
#define TRUNCATION_TEST_UTILS_IMPL

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
#include "test_truncation_utils_decl.hpp"

void calculate_weights_linear_eps(int ntheta,vector<double>& wk_,vector<double>& eps_k,double max_energy,double min_energy){
  assert(ntheta>1 && (ntheta%2)==0);
  double norm=0.0;
  //double max_energy=2.0, min_energy=-2.0;
  double interval=max_energy-min_energy;
  double Delta=interval/ntheta; 
  for(int i=0;i<ntheta;i++){
    
    double simpson;    
    if (i==0 or i== ntheta) simpson=0.0;
    else if (i==1 or i==ntheta-1) simpson=55.0/24.0;
    else if (i==2 or i==ntheta-2) simpson=-1.0/6.0;
    else if (i==3 or i==ntheta-3) simpson=11.0/8.0;
    else simpson=1.0;

    eps_k[i]=min_energy+i*Delta;
    wk_[i]= Delta *simpson* ( 1.0/(2.0*M_PI) * sqrt(4.0 - eps_k[i]*eps_k[i]) ) ;  // integral(f) = h * Somma_su_tutte_le_i di  (simpson[i]*f[i])

    norm+=wk_[i];
  } 
  // renormalize weights to sum_k wk =1
  for(int i=0;i<ntheta;i++) wk_[i]=wk_[i]/norm; 
  //cout<<"wk "<<i<<" = "<<wk_[i]<<endl; 
        
}



void calculate_weights_sin_square(int ntheta,vector<double>& wk_,vector<double>& eps_k ){
  assert(ntheta>1 && (ntheta%2)==0);
  double norm=0.0;
  double Delta_theta=M_PI/ntheta; //this is because eps in the PM goes from -2 to +2
  for(int i=0;i<ntheta;i++){
    double theta=i*Delta_theta;
    
    double simpson;    
    if (i==0 or i== ntheta) simpson=0.0;
    else if (i==1 or i==ntheta-1) simpson=55.0/24.0;
    else if (i==2 or i==ntheta-2) simpson=-1.0/6.0;
    else if (i==3 or i==ntheta-3) simpson=11.0/8.0;
    else simpson=1.0;

    eps_k[i]=8.0/M_PI*(i* (M_PI/2.0)/ntheta*0.5 - 0.125* sin(4.0*  i* (M_PI/2.0)/ntheta )) -2.0;    //1) case one: max of derivative is 8/M_PI
    wk_[i]= Delta_theta *simpson* (    4.0/(M_PI*M_PI) * sin(2*theta)*sin(2*theta) *   sqrt ( 4.0 -  eps_k[i]*eps_k[i])         );  // integral(f) = h * Somma_su_tutte_le_i di  (simpson[i]*f[i])

    norm+=wk_[i];
  } 
  // renormalize weights to sum_k wk =1
  for(int i=0;i<ntheta;i++) wk_[i]=wk_[i]/norm; 
  //cout<<"wk "<<i<<" = "<<wk_[i]<<endl; 
        
}

CPLX timeder(GREEN & G, int up_down,int tstp, int order,double h){
  CPLX timeder_out(0.0,0.0);
  CPLX plusi(0,1), minusi(0,-1);
  if (tstp-order <= 0) { return timeder_out;}
  integration::Integrator<double> obj(order); //Integrator <double>
  for (int l=0;l<= (order+1);++l) timeder_out += obj.bd_weights(l)*G.lesptr(tstp-l,tstp)[up_down];
  timeder_out= minusi*1./h*plusi*timeder_out;
  return timeder_out;
}

CPLX timeder(GTRUNC & G, int tstp_max,int up_down,int tstp, int order,double h){
  CPLX timeder_out(0.0,0.0);
  CPLX plusi(0,1),minusi(0,-1);
  integration::Integrator<double> obj(order);
  for (int l=0;l<= (order+1);++l) timeder_out += obj.bd_weights(l)*G.lesptr(tstp_max-tstp,l)[up_down];
  return (minusi*1./h*plusi*timeder_out);
}

//void print_energy(double u0,double u1,double beta, int nt,int tc, int tmax,int tstp_max, vector <double>& Ekin_up, vector <double> & Ekin_down,vector <double>& rho_up, vector <double>& rho_down, CFUNC & ufunc, int time_interval,bool & prima_volta_energy, GREEN & G, GTRUNC & G_t,double h,int ntheta,int uf_prec){
void print_energy(GREEN & G,GREEN & Sigma,int nt,int kt,double beta,double h,string folder_out,vector <double>& Ekin, string corpus,bool & prima_volta_energy,int tstp_max,int time_interval,int tc,GTRUNC & G_t,vector <double> & n_val,vector <double> & n_val_t){

  std::ofstream output_n;
  string os_n=folder_out+"temp_energy_"+corpus+".out";
  output_n.open(os_n,std::ios_base::app);
  output_n.precision(10);
  int order=5,up=0;//,down=3;
  CPLX timeder_up(0.0,0.0),timeder_do(0.0,0.0);
  CPLX rho1,rho2,rho3;

  int n1= prima_volta_energy == false ? kt : tstp_max-time_interval+1;
  for(int tstp=n1;tstp<=tstp_max;++tstp){
    if(tstp<=tc){
      timeder_up= timeder(G,up,tstp,order,h);
      //timeder_do= timeder(G,down,n,order,h);
    }
    else{
      timeder_up= -timeder(G_t,tstp_max,up,tstp,order,h);
      //timeder_do= -timeder(G_t,tstp_max,down,n,order,h);
    }


    if (tstp<=nt){
      cntr::convolution_density_matrix(tstp,&rho1,G,Sigma,integration::I<double>(kt),beta,h);
      cntr::convolution_density_matrix(tstp,&rho2,G,G,integration::I<double>(kt),beta,h);
      rho3=G.density_matrix(tstp);
      output_n << " t= " <<tstp;
      output_n << " dens= " << rho3.real();
      output_n << " ekin/2= " << rho2.real();
      output_n << " eint= " << rho1.real();
      output_n << " Dens= " << n_val[tstp];
      output_n << " Ekin/2= " << Ekin[tstp];
      output_n << " dt_G= " <<timeder_up.real();
      output_n<< endl;
    }

    else{
      output_n << " t= " <<tstp;
      output_n << " dens= " << 0;
      output_n << " ekin/2= " << 0;
      output_n << " eint= " << 0;
      //truncated code
      output_n << " Dens= " << n_val_t[tstp];
      output_n << " Ekin/2= " << Ekin[tstp];
      output_n << " dt_G= " <<timeder_up.real();
      output_n<< endl;
    }

  }
  output_n.close();
  prima_volta_energy=true;
}


void print_n( vector<double> & n_val, vector<double> n_val_t,  vector<double> & elapsed_time, vector<double> elapsed_time_t,  double u0,double u1, int tmax,double h,int tc,int uf_prec,int nt,string folder,string corpus){
  std::ofstream output_n;

  string os_n = folder+"temp_n_"+corpus+".out";
  //output_nk.open("nk.out");
  output_n.open(os_n);
  output_n.precision(10);

  for(int tstp=1;tstp<=tmax;tstp++){
    output_n << " tstp " << tstp;
    if(tstp<=nt){
      output_n << " n " << n_val[tstp];
      output_n << " elapsed_time " << elapsed_time[tstp];

    }else{
      output_n << " n " << -0.1;
      output_n << " t " << -0.1;

    }
    output_n << " n_t " << n_val_t[tstp];
    output_n << " elapsed_time_trunc " << elapsed_time_t[tstp];

    output_n << endl;
  }
  output_n << endl;
  output_n.close();

  //output_nk.close();
}


void print_nk( vector <vector<double>> & nk_val, vector <vector<double>> &nk_val_t, int neps,vector <double> & ek,double u0,double u1, int tmax,double h,int tc,int uf_prec,int nt,string folder){
  std::ofstream output_nk;
  std::stringstream stream_u0;
  stream_u0 << std::fixed << std::setprecision(1) << u0;
  std::string s_u0 = stream_u0.str();
  std::stringstream stream_u1;
  stream_u1 << std::fixed << std::setprecision(uf_prec) << u1;
  std::string s_u1 = stream_u1.str();

  std::stringstream stream_tmax;
  stream_tmax << std::fixed << std::setprecision(1) << tmax;
  std::string s_tmax = stream_tmax.str();

  std::stringstream stream_h;
  stream_h << std::fixed << std::setprecision(2) << h;
  std::string s_h = stream_h.str();

  //string folder="/nfs/tfkp00/picanoan/Desktop/C++/codes/libcntr/workdir/pm/";
  for(int k=0;k<neps;k++){
    string os_n = folder+"temp_nk_tc_"+ std::to_string(tc)+"_h_"+s_h+ "_ui_"+ s_u0+ "_uf_" +s_u1+ "_k_"+std::to_string(k)+"_tmax_"+s_tmax +".out";
    //output_nk.open("nk.out");
    output_nk.open(os_n);
    output_nk.precision(10);

    for(int tstp=1;tstp<=tmax;tstp++){
      output_nk << "k " << k;
      output_nk << " ek " << ek[k];
      output_nk << " t " << tstp;
      if(tstp<=nt){
	output_nk << " nk " << nk_val[tstp][k];
      }else{
	output_nk << " nk " << -0.1;
      }
      output_nk << " nk_t " << nk_val_t[tstp][k];
      output_nk << endl;
    }
    output_nk << endl;
    output_nk.close();
  }
  //output_nk.close();
}

void get_sigma(int n,GREEN &Sigma,CFUNC &ufunc,GREEN &G){
  cntr::herm_matrix_timestep<double> chi(n,G.ntau(),1,1);
  cntr::Bubble1(n,chi,G,G);
  chi.left_multiply(ufunc,1.0);
  chi.right_multiply(ufunc,-1.0);
  cntr::Bubble2(n,Sigma,G,chi);
}

void get_sigma(int tstp, GTRUNC &Sigma,cntr::function_moving<double> &ufunc,GTRUNC &G,GTRUNC &chi){
  cntr::Bubble1_moving(tstp,chi,G,G);
  chi.left_multiply(ufunc,1.0);
  chi.right_multiply(ufunc,-1.0);
  cntr::Bubble2_moving(tstp,Sigma,G,chi);
}

#endif //TRUNCATION_TEST_UTILS_IMPL
