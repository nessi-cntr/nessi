#include "cntr.hpp"

using namespace std;
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>


TEST_CASE("herm_matrix MPI","[herm_matrix_mpi]"){
  int ntasks,taskid,ierr;
  int master=0;
  
  int size=2;
  int nt=100, ntau=50;
  double eps=1e-6;
  double dt=0.01, mu=0.0, beta=10.0;
  double tmax=dt*nt;
  double eps1=-0.4,eps2=0.6,lam1=0.1;
  std::complex<double> I(0.0,1.0);
  cdmatrix h1(2,2);
  
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

  h1(0,0) = eps1;
  h1(1,1) = eps2;
  h1(0,1) = I*lam1;
  h1(1,0) = -I*lam1;

  SECTION("broadcast + send/receive"){
    double err_loc=0.0;
    double err_glob=0.0;
    
    GREEN G;
    GREEN G_ref;

    G = GREEN(nt,ntau,size,-1);
    G_ref = GREEN(nt,ntau,size,-1);
    
    if(taskid == master) {
      cntr::green_from_H(G,mu,h1,beta,dt);
    } else {
      cntr::green_from_H(G_ref,mu,h1,beta,dt);
    }

    for(int tstp=-1; tstp<=nt; tstp++){
      G.Bcast_timestep(tstp, master);
      if(taskid > master) err_loc += cntr::distance_norm2(tstp,G,G_ref);
    }
  
    MPI_Reduce(&err_loc,&err_glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD);
    if(taskid == master){
      err_glob = err_glob / ntasks;
      REQUIRE(err_glob<eps);
    }

    if(taskid == master) {
      cntr::green_from_H(G_ref,mu,h1,beta,dt);
    } else {
      cntr::green_from_H(G,mu,h1,beta,dt);
    }

    err_loc=0.0;
    for(int tstp=-1; tstp<=nt; tstp++){
      if(taskid > master) {
  	G.Send_timestep(tstp, master, ierr);
      }
      if(taskid == master) {
  	for(int iwork=0; iwork<ntasks; iwork++){
  	  G.Recv_timestep(tstp, iwork, ierr);
  	  err_loc += cntr::distance_norm2(tstp,G,G_ref);
  	}
      }

    }

    if(taskid == master){
      err_loc = err_loc / ntasks;
      REQUIRE(err_loc<eps);
    }
    
  }

     
}

