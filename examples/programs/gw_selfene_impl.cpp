#include <sys/stat.h>
#include <complex>
#include "cntr/cntr.hpp"
#include "gw_selfene_decl.hpp"

using namespace cntr;


namespace diag {
   // Initialize noninteracting green's functions
   void init_G_mat_nointeraction(int Nk_rank,std::vector<gw::kpoint> &corrK_rank,lattice_1d_1b &lattice,int SolverOrder){
      for(int k=0;k<Nk_rank;k++){
         corrK_rank[k].init_G_mat_nointeraction(lattice,SolverOrder);
      }
   }

   // update fermionic propagators via MPI
   void gather_gk_timestep(int tstp,int Nk_rank,DIST_TIMESTEP &gk_all_timesteps,std::vector<gw::kpoint> &corrK_rank,std::vector<int> &kindex_rank){
      gk_all_timesteps.reset_tstp(tstp);
      for(int k=0;k<Nk_rank;k++){
         gk_all_timesteps.G()[kindex_rank[k]].get_data(corrK_rank[k].G_);
      }
      cdmatrix tmp;
      gk_all_timesteps.G()[0].density_matrix(tstp,tmp);
   // distribute to all nodes
      gk_all_timesteps.mpi_bcast_all();
   }

   // update bosonic propagators via MPI
   void gather_wk_timestep(int tstp,int Nk_rank,DIST_TIMESTEP &wk_all_timesteps,std::vector<gw::kpoint> &corrK_rank,std::vector<int> &kindex_rank){
      wk_all_timesteps.reset_tstp(tstp);
      for(int k=0;k<Nk_rank;k++){
         wk_all_timesteps.G()[kindex_rank[k]].get_data(corrK_rank[k].W_);
      }
   // distribute to all nodes
      wk_all_timesteps.mpi_bcast_all();
   }

   // set densities
   void set_density_k(int tstp, int Norb,DIST_TIMESTEP &gk_all_timesteps,lattice_1d_1b &lattice,std::vector<CFUNC> & density_k,std::vector<int> &kindex_rank,CFUNC &rho_loc){
      cdmatrix tmp(Norb,Norb);
      cdmatrix local(Norb,Norb);
      local.setZero();
      for(int k=0;k<lattice.nk_;k++){
         gk_all_timesteps.G()[k].density_matrix(tstp,tmp);
         double wt=lattice.kweight_[k];
         local+=wt*tmp;
         density_k[k].set_value(tstp,tmp);
      }
      rho_loc.set_value(tstp,local);
   }

   // get local propagator
   void get_loc(int tstp,int Ntau,int Norb,lattice_1d_1b &lattice,GREEN &Gloc,DIST_TIMESTEP &gk_all_timesteps){
      cntr::herm_matrix_timestep<double> gtmp(tstp,Ntau,Norb,-1);
      for(int k=0;k<lattice.nk_;k++){
         double wt=lattice.kweight_[k];
         herm_matrix_timestep_view<double> tview(tstp,gtmp);
         tview.incr_timestep(gk_all_timesteps.G()[k],std::complex<double>(wt,0.0));
      }
      Gloc.set_timestep(tstp,gtmp);
   }

   // Evaluate hartree
   void sigma_Hartree(int tstp,int Norb,CFUNC &S,lattice_1d_1b &lattice,std::vector<CFUNC> &density_k,std::vector<CFUNC> &vertex,CFUNC &Ut){
      cdmatrix stmp(Norb,Norb),vloc,vnl;
      stmp.setZero();
      vertex[lattice.G_].get_value(tstp,vnl);
      Ut.get_value(tstp,vloc); //Local component of interaction
      vnl-=vloc; //Get only non-local part of the interaction

      for(int k=0;k<lattice.nk_;k++){
         double wk=lattice.kweight_[k];
         cdmatrix rtmp;
         density_k[k].get_value(tstp,rtmp);
         // Separate local and non-local part
         for(int i1=0;i1<Norb;i1++){
            for(int itmp=0;itmp<Norb;itmp++){
               stmp(i1,i1)+=vloc(i1,itmp)*rtmp(itmp,itmp)*wk; // Local component with Paulli exclusion
               stmp(i1,i1)+=vnl(i1,itmp)*rtmp(itmp,itmp)*wk*2.0; // Non-local component				}
            }
         }
         S.set_value(tstp,stmp);
      }
   }

   // Evaluate fock
   void sigma_Fock(int tstp,int Norb,int kk,CFUNC &S,lattice_1d_1b &lattice,std::vector<CFUNC> &density_k,
      std::vector<CFUNC> &vertex,CFUNC &Ut){
      cdmatrix stmp(Norb,Norb),Vloc(Norb,Norb);
      cdmatrix Xtmp(1,1);
      stmp.setZero();
      Vloc.setZero();
      Ut.get_value(tstp,Vloc); //Local component of interaction
      // Electrons
      for(int q=0;q<lattice.nk_;q++){
         double wk=lattice.kweight_[q];
         cdmatrix rtmp,Vnl;
         int kq;
         kq=lattice.add_kpoints(kk,1,q,-1);
         density_k[kq].get_value(tstp,rtmp);
         vertex[q].get_value(tstp,Vnl);
         Vnl-=Vloc; //Get only non-local part of the interaction
         for(int i1=0;i1<Norb;i1++){
            for(int i2=0;i2<Norb;i2++){
               stmp(i1,i2)-= Vnl(i1,i2)*rtmp(i1,i2)*wk; //Only non-local interaction contribute to the Fock due to exclusion principle
            }
         }
      }
      S.set_value(tstp,stmp);
   }

   //Evaluate GW self-energy
   void sigma_GW(int tstp,int kk,GREEN &S,DIST_TIMESTEP &gk_all_timesteps,DIST_TIMESTEP &wk_all_timesteps,
      lattice_1d_1b &lattice,int Ntau,int Norb){
      assert(tstp==gk_all_timesteps.tstp());
      assert(tstp==wk_all_timesteps.tstp());
      GREEN_TSTP stmp(tstp,Ntau,Norb,FERMION);
      S.set_timestep_zero(tstp);
      for(int q=0;q<lattice.nk_;q++){
         double wk=lattice.kweight_[q];
         int kq=lattice.add_kpoints(kk,1,q,-1);
         stmp.clear();
         for(int i1=0;i1<Norb;i1++){
            for(int i2=0;i2<Norb;i2++){
               cntr::Bubble2(tstp,stmp,i1,i2,gk_all_timesteps.G()[kq],gk_all_timesteps.G()[kq],i1,i2,
                  wk_all_timesteps.G()[q],wk_all_timesteps.G()[q],i1,i2);
            }
         }
         S.incr_timestep(tstp,stmp,wk);

      }
   }

   //Evaluate polarization
   void get_Polarization_Bubble(int tstp,int Norb,int Ntau,int qq,GREEN &P,DIST_TIMESTEP &gk_all_timesteps,lattice_1d_1b &lattice){
      assert(tstp==gk_all_timesteps.tstp());	
      // compute P_{a,a'}(q;t,t') on a timestep
      int kq;
      GREEN_TSTP ptmp(tstp,Ntau,Norb,BOSON);
      P.set_timestep_zero(tstp);
      for(int kk=0;kk<lattice.nk_;kk++){
      double wk=lattice.kweight_[kk]; // factor 2: ksum normalized for RBZ
      // get bubble 
      // -ii*sum_{j1,j2} G_{(a1),(a2)}(k-q;t,t') G_{(a2),(a1)}(k;t',t)
      kq=lattice.add_kpoints(kk,1,qq,-1);
      ptmp.clear();
      for(int i1=0;i1<Norb;i1++){
         for(int i2=0;i2<Norb;i2++){
            cntr::Bubble1(tstp,ptmp,i1,i2,gk_all_timesteps.G()[kq],gk_all_timesteps.G()[kq],i1,i2,gk_all_timesteps.G()[kk],
            gk_all_timesteps.G()[kk],i1,i2);										  
         }
      }
      // increment P += -wt*Ptmp
      P.incr_timestep(tstp,ptmp,-wk);
      }
   }

   //Extrapolate single-particle density matrix
   void extrapolate_rho(int tstp,int Nk,std::vector<CFUNC> &density_k){
      if(tstp<-1){
         return;	
      } 
      else if(tstp==0){
         for(int k=0;k<Nk;k++){
            cdmatrix tmp;
            density_k[k].get_value(-1,tmp);
            density_k[k].set_value(0,tmp);
         }
      }
      else if(tstp>0){
         for(int k=0;k<Nk;k++){
            cdmatrix tmp;
            density_k[k].get_value(tstp-1,tmp);
            density_k[k].set_value(tstp,tmp);
         }
      }
   }

   //Extrapolate Green's function
   void extrapolate_timestep_G(int tstp,int Nk_rank,int SolverOrder,int Nt,std::vector<gw::kpoint> &corrK_rank){
      assert(tstp<=Nt-1);
   // extrapolate from tstp to tstp+1
      if(tstp<-1){
         return;	
      }else if(tstp==-1){
         for(int k=0;k<Nk_rank;k++){
            cntr::set_t0_from_mat(corrK_rank[k].G_);
         }
      }else if(tstp>=0){
         for(int k=0;k<Nk_rank;k++){
            cntr::extrapolate_timestep(tstp,corrK_rank[k].G_,SolverOrder);
         }
      }
   }

   //Extrapolate bosonic propagator
   void extrapolate_timestep_W(int tstp,int Nk_rank,int SolverOrder,int Nt,std::vector<gw::kpoint> &corrK_rank){
      assert(tstp<=Nt-1);
   // extrapolate from tstp to tstp+1
      if(tstp<-1){
         return;	
      }else if(tstp==-1){
         for(int k=0;k<Nk_rank;k++){
            cntr::set_t0_from_mat(corrK_rank[k].W_);
         }
      }else if(tstp>=0){
         for(int k=0;k<Nk_rank;k++){
            cntr::extrapolate_timestep(tstp,corrK_rank[k].W_,SolverOrder);
         }
      }
   }

   //Evaluate the kinetic energy
   double KineticEnergy(int tstp,lattice_1d_1b &lattice,std::vector<CFUNC> & density_k){
      cdmatrix rtmp,hktmp;
      double ekin=0.0;
      for(int k=0;k<lattice.nk_;k++){
         double wt=lattice.kweight_[k];
         density_k[k].get_value(tstp,rtmp);
         lattice.hk(hktmp,tstp,lattice.kpoints_[k],0);
         ekin+=std::real((hktmp*rtmp).trace())*wt;
      }
      return ekin;
   }

   //Evaluate the potential energy
   double CorrelationEnergy(int tstp,int Nk_rank,int SolverOrder,double beta, double h,std::vector<gw::kpoint> &corrK_rank,
      lattice_1d_1b &lattice){
      cdmatrix tmpH,tmpF,tmpR;
      double etmp=0.0,etot=0.0;
      for(int k=0;k<Nk_rank;k++){
         double wk=lattice.kweight_[k]; 
         corrK_rank[k].rho_.get_value(tstp,tmpR);
         corrK_rank[k].SHartree_.get_value(tstp,tmpH);
         corrK_rank[k].SFock_.get_value(tstp,tmpF);
         etmp+=std::real(((tmpH+tmpF)*tmpR).trace())*wk; //Actually *2 due to spin and /2 from the Migdal formula
         etmp+=cntr::correlation_energy(tstp, corrK_rank[k].G_,corrK_rank[k].Sigma_,
            beta,h,SolverOrder)*wk; //Actually *2 due to spin and /2 from the Migdal formula

      }
      MPI_Allreduce(&etmp,&etot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD);
      // MPI::COMM_WORLD.Allreduce(&etmp,&etot,1,MPI::DOUBLE,MPI_SUM);
      return etot;
      }
   }
