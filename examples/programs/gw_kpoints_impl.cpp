#include "gw_kpoints_decl.hpp"

// -----------------------------------------------------------------------
#define CFUNC cntr::function<double>
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
// -----------------------------------------------------------------------
namespace gw{

/* #######################################################################################
#
#   CONSTRUCTION/DESTRUCTION
#
########################################################################################*/
    
    kpoint::kpoint(void){}
    kpoint::~kpoint(void){}

/** \brief <b> Initializes the `kpoint` class   </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Initializes the `kpoint` class for propagators at a given k point.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param nt
* > Number of time steps
* @param ntau
* > Number of points on Matsubara axis
* @param size
* > Matrix rank of the contour function
* @param beta
* > Inverse temperature
* @param h
* > timestep
* @param kk
* > momentum
* @param latt
* > The corresponding lattice_1d_1b with momentum points and energies 
* @param mu
* > chemical potential
*/
  kpoint::kpoint(int nt,int ntau,int size,double beta,double h,double kk,lattice_1d_1b &latt,double mu){
  
    beta_=beta;
    h_=h;
    nt_=nt;
    ntau_=ntau;
    nrpa_=size;
    mu_=mu;
    //Set up instantaneous quantities
    SHartree_=CFUNC(nt_,nrpa_);
    SFock_=CFUNC(nt_,nrpa_);
    rho_=CFUNC(nt_,nrpa_);
    hk_=CFUNC(nt_,nrpa_);
    hkeff_=CFUNC(nt_,nrpa_);
    vertex_=CFUNC(nt_,nrpa_);

    kk_=kk;
    for(int tstp=-1;tstp<=nt_;tstp++) set_hk(tstp,1,latt);
    for(int tstp=-1;tstp<=nt_;tstp++) set_vertex(tstp,latt);

    //Set up propagators
    Sigma_=GREEN(nt_,ntau_,nrpa_,-1);
    G_=GREEN(nt_,ntau_,nrpa_,-1);
    G0_=GREEN(nt_,ntau_,nrpa_,-1);

    G0Sigma_=GREEN(nt_,ntau_,nrpa_,-1);
    SigmaG0_=GREEN(nt_,ntau_,nrpa_,-1);
    P_=GREEN(nt_,ntau_,nrpa_,+1);
    PV_=GREEN(nt_,ntau_,nrpa_,+1);
    VP_=GREEN(nt_,ntau_,nrpa_,+1);
    W_=GREEN(nt_,ntau_,nrpa_,+1);

    cdmatrix tmp(nrpa_,nrpa_);
    tmp.setZero();
    cntr::green_from_H(G0_,mu_,tmp,beta_,h_);
  }

  
  void kpoint::set_hk(int tstp,int iter,lattice_1d_1b &latt){
    assert(-1<=tstp && tstp<=nt_);
    cdmatrix hktmp(nrpa_,nrpa_);
    latt.hk(hktmp,tstp,kk_,iter);
    hk_.set_value(tstp,hktmp);
  }

  
  void kpoint::set_vertex(int tstp,lattice_1d_1b &latt){
    assert(-1<=tstp && tstp<=nt_);
    cdmatrix vktmp(nrpa_,nrpa_);
    latt.V(tstp,kk_,vktmp);
    vertex_.set_value(tstp,vktmp);
  }

  void kpoint::get_Density_matrix(int tstp){
    cdmatrix tmp(nrpa_,nrpa_);
    G_.density_matrix(tstp,tmp);
    rho_.set_value(tstp,tmp);
  }

  void kpoint::init_G_mat_nointeraction(lattice_1d_1b &latt,int kt){
    // set G = (idt + mu - hk(U=V=0=A=V01) )^{-1}
    int tstp=-1;
    cdmatrix hktmp;
    // latt.hkfree(hktmp,tstp,kk_);
    hkeff_=hk_;
    cdmatrix tmp;
    cntr::green_from_H(G_,mu_,hkeff_,beta_,h_,5,4,false);
    
  }

  void kpoint::step_W(int tstp,int kt,lattice_1d_1b &latt){
    // solve W = V + V*Pi*W, assuming P is set on all relevant timesteps
    int n;
    int n1=(tstp==-1 || tstp>kt ? tstp : 0);
    int n2=(tstp==-1 || tstp>kt ? tstp : kt);   
    //for(n=n1;n<=n2;n++) chi_.set_timestep(n,P_);
    //return;
    // set PV=-P*V etc.
    for(n=n1;n<=n2;n++){
      PV_.set_timestep(n,P_);
      PV_.right_multiply(n,vertex_);
      PV_.smul(n,-1.0);
      VP_.set_timestep(n,P_);
      VP_.left_multiply(n,vertex_);
      VP_.smul(n,-1.0);
    }
    // solve [1-UP]*W=U
    CFUNC tmpFsin(nt_,nrpa_);
    tmpFsin.set_zero();
    GREEN Q(nt_,ntau_,nrpa_,1);

    cntr::vie2_timestep_sin(tstp,W_,vertex_,VP_,PV_,tmpFsin,Q,vertex_,beta_,h_,kt);
  }

  void kpoint::step_dyson(int tstp,int iter,int kt,lattice_1d_1b &latt){
    // solve G=(-idt+mu-hk-Hartree-Fock-Sigma)^{-1}, assuming Sigma is set
    int n;
    int n1=(tstp==-1 || tstp>kt ? tstp : 0);
    int n2=(tstp==-1 || tstp>kt ? tstp : kt);   
    // set hkeff=hk+Hartree+Fock+Symmetry breaking
    for(n=n1;n<=n2;n++){
      cdmatrix tmp,tmp1,tmp2,rtmp(nrpa_,nrpa_);
      set_hk(n,iter,latt);
      SHartree_.get_value(n,tmp);
      SFock_.get_value(n,tmp1);
      // std::cout <<  "Hartree " <<tstp << " " << iter << " " << kk_ << " " << tmp <<std::endl;
      // std::cout <<  "Fock " <<tstp << " " << iter << " " << kk_ << " " << tmp1 <<std::endl;
      tmp+=tmp1;
      hk_.get_value(n,tmp2);
      
      tmp+=tmp2;
      hkeff_.set_value(n,tmp);
      // std::cout <<  "hkeff " << tmp <<std::endl;
    }
    
    // solve Dyson
    if(tstp==-1 && iter==1){
      GREEN_TSTP tmp(-1,ntau_,nrpa_,G_.sig());
      cntr::green_from_H(tmp,mu_,hkeff_,beta_,h_,kt,4,false);
      G_.set_timestep(-1,tmp);
    }else if(tstp==-1){
      cntr::dyson_mat(G_,Sigma_,mu_,hkeff_,integration::I<double>(kt),beta_,0,true);
    }else if (tstp==0){
      // TODO: Add epsilon to the initial condition to check the linear response
        cntr::set_t0_from_mat(G_);
    }else if(tstp<=kt){
      cntr::dyson_start(G_,mu_,hkeff_,Sigma_,integration::I<double>(kt),beta_,h_);
    }else{
        cdmatrix rtmp(nrpa_,nrpa_);
        cdmatrix Xtmp_prev(1,1),Ptmp_prev(1,1);
        cntr::dyson_timestep(tstp,G_,mu_,hkeff_,Sigma_,integration::I<double>(kt),beta_,h_);      
    }
    for(n=n1;n<=n2;n++) get_Density_matrix(n);
  }
  
  double kpoint::step_dyson_with_error(int tstp,int iter,int kt,lattice_1d_1b &latt){
    GREEN_TSTP gtmp;
    G_.get_timestep(tstp,gtmp);
    step_dyson(tstp,iter,kt,latt);
    return cntr::distance_norm2(tstp,gtmp,G_);
  }

  double kpoint::step_W_with_error(int tstp,int kt,lattice_1d_1b &latt){
    GREEN_TSTP gtmp;
    W_.get_timestep(tstp,gtmp);
    step_W(tstp,kt,latt);
    return cntr::distance_norm2(tstp,gtmp,W_);
  }
}