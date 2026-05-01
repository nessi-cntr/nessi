#ifndef NESS2_DYSON_IMPL_H
#define NESS2_DYSON_IMPL_H

#include "ness2_global_settings.hpp"
//#include "ness2_matrix_decl.hpp"
#include "ness2_dyson_decl.hpp"


namespace ness2{

  void dyson(herm_matrix_ness &G, double mu, cdmatrix &epsilon, herm_matrix_ness &Sigma, double h,
	     fft_integral_method method,dyson_bath_type bath_type,double bath_eta,double bath_beta,double bath_mu){
    int size1=G.size1_;
    int Nft=G.Nft_;
    if(Sigma.size1_!=size1 || Sigma.Nft_!=Nft || epsilon.rows()!=size1 ||  epsilon.cols()!=size1 ){
      throw std::runtime_error("dyson_fft:input param mismatch");
    }    
    if(method==FFT_TRAPEZ){
      dyson_fft(G,Sigma,epsilon,mu,h,bath_type,bath_eta,bath_beta,bath_mu);
    }else{
      throw std::runtime_error("ness::dyson: Only FFT_TRAPEZ implemented");
    }
  }
  //@private
  void sigma_bath(cdouble &ret,cdouble &les, double omega,double h,dyson_bath_type bath_type,double bath_eta,double bath_beta,double bath_mu){
    double fm,arg;
    cdouble ii(0,1);
    if(bath_type==BATH_CONST){
      ret=-ii*bath_eta;      
      arg = bath_beta * (omega - bath_mu);
      if (arg > 100) fm=0.0;
      else if (arg < -100) fm=1.0;
      else fm=1.0/(1.0 + std::exp(arg));      
      les=ii*2.0*fm*bath_eta;      
    }else if(bath_type==BATH_GAUSS){
      double omegac=M_PI/(h*10.0);
      double expfac=exp(-omega*omega/omegac*omegac);
      ret=-ii*bath_eta*expfac;      
      arg = bath_beta * (omega - bath_mu);
      if (arg > 100) fm=0.0;
      else if (arg < -100) fm=1.0;
      else fm=1.0/(1.0 + std::exp(arg));      
      les=ii*2.0*fm*bath_eta*expfac;      
    }else if(bath_type==BATH_OHMIC){
      double omegac=M_PI/(h*10.0);
      double expfac=exp(-omega*omega/omegac*omegac);
      ret=-ii*bath_eta*omega*expfac;
      arg = bath_beta * omega;
      if (arg > 100) les=0.0;
      else if (arg < -100) les=ii*2.0*bath_eta*omega*expfac;
      else if (std::abs(arg)<1e-6)  les=-ii*2.0*bath_eta*expfac/bath_beta;
      else les=-ii*2.0*bath_eta*omega*expfac/(std::exp(arg)-1.0);
    }else{
      ret=0;
      les=0;    
    }
  }
  
  
  //@private
  void dyson_fft(herm_matrix_ness  & G, herm_matrix_ness  & Sigma, cdmatrix & H,double mu,double h,
		 dyson_bath_type bath_type,double bath_eta,double bath_beta,double bath_mu){
    int size1=G.size1_;
    int Nft=G.Nft_;
    int nt=Nft/2-1;
    int element_size=size1*size1;
    grid_info grid(Nft,h);
    cdmatrix Id = cdmatrix::Identity(size1, size1);
    // the following integral transform used only a factor 1/2 at Sigma(0) as "boundary correction"
    Sigma.integral_transform_to_freq(h,FFT_TRAPEZ);
    // only points  -(nfreq-1) ...  0...(nfreq+1), are used wheer nfreq=nt
    int nfreq=(Nft/2-1)/3;
    //int nfreq=(Nft/2-1);
    G.ret_.set_zero(fft_domain::freq);
    G.les_.set_zero(fft_domain::freq);        
    for(int w = -nfreq+1; w < nfreq; ++w){
      cdmatrix gret,sret,sles;
      cdouble bath_ret=0.0,bath_les=0;
      double omega=grid.freq_at(w);
      Sigma.get_ret(w,sret,fft_domain::freq);
      if(bath_type){
	sigma_bath(bath_ret,bath_les, omega,h,bath_type,bath_eta,bath_beta,bath_mu);
      }
      cdmatrix ginv=Id*(omega+mu)-H-sret-bath_ret*Id;
      gret=ginv.inverse();
      G.set_ret(w,gret,fft_domain::freq);
      Sigma.get_les(w,sles,fft_domain::freq);
      G.set_les(w,gret*(sles+bath_les*Id)*gret.adjoint(),fft_domain::freq);
    }
    // this transforms back the imaginary part of G(w) only
    G.integral_transform_to_time(h,FFT_TRAPEZ);
  }

} //namespace ness
#endif  // NESS_EQUILIBRIUM_IMPL_H
