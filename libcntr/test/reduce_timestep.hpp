#include "cntr.hpp"

using namespace std;
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>


TEST_CASE("Reduce_timestep","[Reduce_timestep]"){
  int ntasks,taskid,ierr;
  int master=0;
  
  int size=2;
  int nt=50, ntau=50;
  double eps=1e-6;
  double dt=0.01, mu=0.0, beta=10.0;
  double tmax=dt*nt;
  double eps1=-0.4,eps2=0.6,lam1=0.1;
  std::complex<double> I(0.0,1.0);
  cdmatrix h1(size,size);
  cdmatrix iden(size,size);
  iden = MatrixXd::Identity(size, size);
  
  ntasks = MPI::COMM_WORLD.Get_size();
  taskid = MPI::COMM_WORLD.Get_rank();

  h1(0,0) = eps1;
  h1(1,1) = eps2;
  h1(0,1) = I*lam1;
  h1(1,0) = -I*lam1;

  SECTION("herm_matrix: member function"){
    
    GREEN Gk;
    GREEN Gloc;

    Gk = GREEN(nt,ntau,size,-1);

    if(taskid == master) {
      Gloc = GREEN(nt,ntau,size,-1);
      for(int tstp=-1; tstp<=nt; tstp++) Gloc.set_timestep_zero(tstp);
      GREEN Gk_master(nt,ntau,size,-1);
      cdmatrix hk(size, size);
      for(int tid=0; tid<ntasks; tid++){
        hk = h1 + (double)tid * iden;
        cntr::green_from_H(Gk_master,mu,hk,beta,dt);
        for(int tstp=-1; tstp<=nt; tstp++) Gloc.incr_timestep(tstp, Gk_master);
      }      
    } 

    cdmatrix hk(size, size);
    hk = h1 + (double)taskid * iden;
    cntr::green_from_H(Gk,mu,hk,beta,dt);

    for(int tstp=-1; tstp<=nt; tstp++){
      Gk.Reduce_timestep(tstp, master);
    }

    if(taskid == master){
      double err=0.0;
      for(int tstp=-1; tstp<=nt; tstp++){
        err += cntr::distance_norm2(tstp, Gk, Gloc);
      }
      REQUIRE(err < eps);
    }

    
  }

  SECTION("herm_matrix_timestep: member function"){
    
    double err = 0.0;
    GREEN Gk;
    GREEN Gloc;

    Gk = GREEN(nt,ntau,size,-1);

    if(taskid == master) {
      Gloc = GREEN(nt,ntau,size,-1);
      for(int tstp=-1; tstp<=nt; tstp++) Gloc.set_timestep_zero(tstp);
      GREEN Gk_master(nt,ntau,size,-1);
      cdmatrix hk;
      for(int tid=0; tid<ntasks; tid++){
        hk = h1 + tid * iden;
        cntr::green_from_H(Gk_master,mu,hk,beta,dt);
        for(int tstp=-1; tstp<=nt; tstp++) Gloc.incr_timestep(tstp, Gk_master);
      }      
    } 

    cdmatrix hk;
    hk = h1 + taskid * iden;
    cntr::green_from_H(Gk,mu,hk,beta,dt);

    for(int tstp=-1; tstp<=nt; tstp++){
      GREEN_TSTP G_slice(tstp,ntau,size);
      Gk.get_timestep(tstp, G_slice);
      G_slice.Reduce_timestep(tstp, master);

      if(taskid == master){
        err += cntr::distance_norm2(tstp, G_slice, Gloc);
      }
    }

    REQUIRE(err < eps);
    
  }

  SECTION("Reduce_timestep: herm_matrix_timestep -> herm_matrix_timestep"){
    
    double err = 0.0;
    GREEN Gk;
    GREEN Gloc;
    GREEN Gloc_ref;

    Gk = GREEN(nt,ntau,size,-1);

    if(taskid == master) {
      Gloc = GREEN(nt,ntau,size,-1);
      Gloc_ref = GREEN(nt,ntau,size,-1);
      for(int tstp=-1; tstp<=nt; tstp++) Gloc_ref.set_timestep_zero(tstp);
      GREEN Gk_master(nt,ntau,size,-1);
      cdmatrix hk;
      for(int tid=0; tid<ntasks; tid++){
        hk = h1 + tid * iden;
        cntr::green_from_H(Gk_master,mu,hk,beta,dt);
        for(int tstp=-1; tstp<=nt; tstp++) Gloc_ref.incr_timestep(tstp, Gk_master);
      }      
    } 

    cdmatrix hk;
    hk = h1 + taskid * iden;
    cntr::green_from_H(Gk,mu,hk,beta,dt);

    for(int tstp=-1; tstp<=nt; tstp++){
      GREEN_TSTP Gk_step(tstp, ntau, size);
      GREEN_TSTP Gloc_step;
      if(taskid == master) Gloc_step = GREEN_TSTP(tstp, ntau, size);
      Gk.get_timestep(tstp, Gk_step);
      cntr::Reduce_timestep(tstp, master, Gloc_step, Gk_step);
      if(taskid == master){
        Gloc.set_timestep(tstp, Gloc_step);
        err += cntr::distance_norm2(tstp, Gloc, Gloc_ref);
      }
    }

    REQUIRE(err < eps);
    
  }

  SECTION("Reduce_timestep: herm_matrix -> herm_matrix_timestep"){
    
    double err = 0.0;
    GREEN Gk;
    GREEN Gloc;
    GREEN Gloc_ref;

    Gk = GREEN(nt,ntau,size,-1);

    if(taskid == master) {
      Gloc = GREEN(nt,ntau,size,-1);
      Gloc_ref = GREEN(nt,ntau,size,-1);
      for(int tstp=-1; tstp<=nt; tstp++) Gloc_ref.set_timestep_zero(tstp);
      GREEN Gk_master(nt,ntau,size,-1);
      cdmatrix hk;
      for(int tid=0; tid<ntasks; tid++){
        hk = h1 + tid * iden;
        cntr::green_from_H(Gk_master,mu,hk,beta,dt);
        for(int tstp=-1; tstp<=nt; tstp++) Gloc_ref.incr_timestep(tstp, Gk_master);
      }      
    } 

    cdmatrix hk;
    hk = h1 + taskid * iden;
    cntr::green_from_H(Gk,mu,hk,beta,dt);

    for(int tstp=-1; tstp<=nt; tstp++){
      GREEN_TSTP Gloc_step;
      if(taskid == master) Gloc_step = GREEN_TSTP(tstp, ntau, size);
      cntr::Reduce_timestep(tstp, master, Gloc_step, Gk);
      if(taskid == master){
        Gloc.set_timestep(tstp, Gloc_step);
        err += cntr::distance_norm2(tstp, Gloc, Gloc_ref);
      }
    }

    REQUIRE(err < eps);
    
  }

  SECTION("Reduce_timestep: herm_matrix_timestep -> herm_matrix"){
    
    double err = 0.0;
    GREEN Gk;
    GREEN Gloc;
    GREEN Gloc_ref;

    Gk = GREEN(nt,ntau,size,-1);

    if(taskid == master) {
      Gloc = GREEN(nt,ntau,size,-1);
      Gloc_ref = GREEN(nt,ntau,size,-1);
      for(int tstp=-1; tstp<=nt; tstp++) Gloc_ref.set_timestep_zero(tstp);
      GREEN Gk_master(nt,ntau,size,-1);
      cdmatrix hk;
      for(int tid=0; tid<ntasks; tid++){
        hk = h1 + tid * iden;
        cntr::green_from_H(Gk_master,mu,hk,beta,dt);
        for(int tstp=-1; tstp<=nt; tstp++) Gloc_ref.incr_timestep(tstp, Gk_master);
      }      
    } 

    cdmatrix hk;
    hk = h1 + taskid * iden;
    cntr::green_from_H(Gk,mu,hk,beta,dt);

    for(int tstp=-1; tstp<=nt; tstp++){
      GREEN_TSTP Gk_step(tstp, ntau, size);
      Gk.get_timestep(tstp, Gk_step);
      cntr::Reduce_timestep(tstp, master, Gloc, Gk_step);
      if(taskid == master){
        err += cntr::distance_norm2(tstp, Gloc, Gloc_ref);
      }
    }

    REQUIRE(err < eps);
    
  }

  SECTION("Reduce_timestep: herm_matrix -> herm_matrix"){
    
    double err = 0.0;
    GREEN Gk;
    GREEN Gloc;
    GREEN Gloc_ref;

    Gk = GREEN(nt,ntau,size,-1);

    if(taskid == master) {
      Gloc = GREEN(nt,ntau,size,-1);
      Gloc_ref = GREEN(nt,ntau,size,-1);
      for(int tstp=-1; tstp<=nt; tstp++) Gloc_ref.set_timestep_zero(tstp);
      GREEN Gk_master(nt,ntau,size,-1);
      cdmatrix hk;
      for(int tid=0; tid<ntasks; tid++){
        hk = h1 + tid * iden;
        cntr::green_from_H(Gk_master,mu,hk,beta,dt);
        for(int tstp=-1; tstp<=nt; tstp++) Gloc_ref.incr_timestep(tstp, Gk_master);
      }      
    } 

    cdmatrix hk;
    hk = h1 + taskid * iden;
    cntr::green_from_H(Gk,mu,hk,beta,dt);

    for(int tstp=-1; tstp<=nt; tstp++){
      cntr::Reduce_timestep(tstp, master, Gloc, Gk);
      if(taskid == master){
        err += cntr::distance_norm2(tstp, Gloc, Gloc_ref);
      }
    }

    REQUIRE(err < eps);
    
  }

     
}

