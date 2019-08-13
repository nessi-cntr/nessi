#include "cntr.hpp"


using namespace std;
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>


TEST_CASE("distributed_timestep_array MPI","[distributed_timestep_array_mpi]"){
  int ntasks,taskid,ierr;
  int master=0;
  double eps=1e-6;
  int nblock=2;
  double err=0.0;

  int size=2;
  int nt=10, ntau=50;
  double dt=0.01, mu=0.0, beta=10.0;
  double lam=0.1;
  int npoints=6;
  std::complex<double> I(0.0,1.0);
  cdmatrix h1(2,2);

  ntasks = MPI::COMM_WORLD.Get_size();
  taskid = MPI::COMM_WORLD.Get_rank();

  std::vector<GREEN> Gvec(npoints);
  for(int i=0;i<npoints;i++){
      h1(0,0) = i;
      h1(1,1) = i;
      h1(0,1) = I*lam;
      h1(1,0) = -I*lam;
      Gvec[i]=GREEN(nt,ntau,size,-1);
      cntr::green_from_H(Gvec[i],mu,h1,beta,dt);
  } 
  /*
  Idea of test: The are globally npoints of Greens function and in
  distributed_timestep_array each mpi thread has information about some subset of those.
  Via mpi_bcast_all [in usual MPI  nomenclature this is Allgather] we distribute the information
  to all threads and test it. 
  */
  SECTION("BcastAll"){
    cntr::distributed_timestep_array<double> Gall(npoints,nt,ntau,size,-1,true);
    std::vector<int> mpi_pid=Gall.data().tid_map();

    for(int tstp=-1;tstp<=nt;tstp++){
      Gall.reset_tstp(tstp);
      // Set 
      for(int i=0;i<npoints;i++){
        if(taskid==mpi_pid[i]){
          cdmatrix mat;
          Gall.G()[i].get_data(Gvec[i]);
        }
      }

      Gall.mpi_bcast_all();

      for(int i=0;i<npoints;i++){
          err+=distance_norm2(tstp,Gall.G()[i],Gvec[i]);

      }
    }
    REQUIRE(err<eps);
  }
}