#include "catch.hpp"
#include "cntr.hpp"
#include <mpi.h>

using namespace std;
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>


TEST_CASE("distributed_array MPI","[distributed_array_mpi]"){
  int ntasks,taskid,ierr;
  int master=0;
  int n=10;
  
  MPI::Init();
  ntasks = MPI::COMM_WORLD.Get_size();
  taskid = MPI::COMM_WORLD.Get_rank();

  cntr::distributed_array<double> A(n,1,true);
  
  SECTION("broadcast + send/receive"){
    double err_loc=0.0;
    double err_glob=0.0;
    
    if(taskid == master) {
      for(int i=0;i<n;i++){
	double *x = A.block(i);
	for(int j=0;j<maxlen;j++){
	  std::cout << " Inside 1 M" << j << " " << j << " " << j+i*nblock <<std::endl;
	  x[j]=j+i*nblock;
	  std::cout << " Inside 1 MA" << x[j] <<std::endl;
	}
      }
    } else {
      A.clear();
    }

    for(int j=0; j<=nblock; j++){
      A.mpi_bcast_block(j);
      if(taskid > master){
	for(int i=0;i<nblock;i++){
	  double *x=A.block(i);
	  for(int j=0;j<size;j++){
	    std::cout << " Inside 2 M" << taskid << j << " " << j << " " << j+i*nblock <<std::endl;
	    err_loc += (j+i*nblock-*x);
	    std::cout << " Inside 2 MA" << *x << " " << j+i*nblock  <<std::endl;
	  }
	}
      }
    }
  
    MPI_Reduce(&err_loc,&err_glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD);
    if(taskid == master){
      err_glob = err_glob / ntasks;
      REQUIRE(err_glob<eps);
    }
  }
  MPI::Finalize();
}

