#include "gw_latt_decl.hpp"
lattice_1d_1b::lattice_1d_1b(void){
}

lattice_1d_1b::lattice_1d_1b(int nk,int nt,CFUNC &tt,CFUNC &Ut,CFUNC &Vt,std::vector<double> &Epulse,double mu,int Norb,int kt,double h):
	nt_(nt),Norb_(Norb),mu_(0.0)
{
  assert(-1<=nt);
  assert(-1<=nt);
  assert(nt==(int)Ut.nt());
  tt_=tt;
  U_=Ut;
  V_=Vt;
  mu_=mu;
  cdmatrix tmp(Norb,Norb);
  Apulse_=CFUNC(nt,1);
  efield_to_afield(nt,h,Epulse,kt);  
  
  init_kk(nk);
}

void lattice_1d_1b::efield_to_afield(int nt,double h,std::vector<double> &efield,int kt){
  int kt1=(nt>=kt ? kt : nt),n,n1,tstp;
  cdmatrix At(1,1);
  At.setZero();

  assert((int)efield.size()>=nt+2);
  At(0,0)=0.0;
  Apulse_.set_value(-1,At); // tstp=-1
  for(tstp=0;tstp<=nt;tstp++){
    At(0,0)=0.0;
    n1=(tstp<kt1 ? kt1 : tstp);
    for(n=0;n<=n1;n++){
      At(0,0) += integration::I<double>(kt1).gregory_weights(tstp,n)*efield[n+1];
    }
    At(0,0)*=(-h);
    //std::cout << "n: " << efield[tstp+1] << " " << At << std::endl; 
    Apulse_.set_value(tstp,At);
  }
}


int lattice_1d_1b::add_kpoints(int k1,int s1,int k2,int s2){
  int k12=s1*(k1-G_)+s2*(k2-G_)+G_;
  while (k12<0){k12 += nk_;}
  while (k12>=nk_){k12 -= nk_;}
  return k12;
}

void lattice_1d_1b::init_kk(int nk){
  double dk,kw;
  assert(nk>1);
  nk_=nk;
  G_=nk_/2;
  dk=2*PI/nk_; // kk= (index-G_)*dk 
  kw=1.0/nk_;
  kpoints_.resize(nk_);
  kweight_.resize(nk_);	
  for(int i1=0;i1<nk_;i1++){
    kpoints_[i1]=(i1-G_)*dk;
  }
  for(int i1=0;i1<nk_;i1++) kweight_[i1]=kw;
}


void lattice_1d_1b::hk(cdmatrix &hkmatrix,int tstp,double kk,int iter){
  cdmatrix tmpTT,tmpA;
  tt_.get_value(tstp,tmpTT);
  Apulse_.get_value(tstp,tmpA);
  double epsk=-2.0*std::real(tmpTT(0,0))*cos(kk-std::real(tmpA(0,0)));
  hkmatrix.resize(Norb_,Norb_);
  hkmatrix.setZero();
  hkmatrix(0,0)=epsk;
}

void lattice_1d_1b::vk(cdmatrix &vkmatrix,int tstp,double kk){
	cdmatrix tmpTT;
  tt_.get_value(tstp,tmpTT);
  double vk=2.0*std::real(tmpTT(0,0))*sin(kk);
	vkmatrix.resize(Norb_,Norb_);
	vkmatrix.setZero();
	vkmatrix(0,0)=vk;
}

// interaction vertex V_{a1;a2}(q) = V_{a1,a2}(q) 
void lattice_1d_1b::V(int tstp,double k,cdmatrix &V){
  cdmatrix Utmp(Norb_,Norb_),Vtmp(Norb_,Norb_);
  U_.get_value(tstp,Utmp);
  V_.get_value(tstp,Vtmp);
  V=Utmp+2.0*Vtmp*cos(k);
}
