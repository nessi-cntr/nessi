#include "cntr.hpp"


using namespace std;
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>


TEST_CASE("distributed_array MPI","[distributed_array_mpi]"){
  int ntasks,taskid,ierr;
  int master=0;
  int nblock=10;
  int blocksize=5;
  double eps=1e-6;

  ntasks = MPI::COMM_WORLD.Get_size();
  taskid = MPI::COMM_WORLD.Get_rank();

  /*
    The goal of the test is to gather data on all blocks to certain rank
    For example on 2 ranks at the initiation the whole array should be
    rank0: [abcd,0000]
    rank1: [0000,efgh]
    and when we gather on rank 0 we get
    rank0: [abcd,efgh]
    rank1: [0000,efgh]
  */

  SECTION("Gather"){
    cntr::distributed_array<double> A(nblock,blocksize,true);
    double err_loc=0.0;
    double err_glob=0.0;
    

    for(int j=0;j<nblock;j++){
      if(taskid==A.tid_map()[j]){
        for(int i=0;i<blocksize;i++){
          A.block(j)[i]=j+i*nblock;
        }
      }
    }

    A.mpi_gather(0);
    
    if(taskid==0){
      for(int i=0;i<nblock;i++){
        for(int j=0;j<blocksize;j++){
          err_loc += (i+j*nblock-A.block(i)[j]);
        }
      }
      REQUIRE(err_loc<eps);
    }

    A.mpi_gather(1);
    
    if(taskid==1){
      for(int i=0;i<nblock;i++){
        for(int j=0;j<blocksize;j++){
          err_loc += (i+j*nblock-A.block(i)[j]);
        }
      }
      REQUIRE(err_loc<eps);
    }
    
    
  }


  /*
    The goal of the test is to set each block with some data
    and after broadcast the all blocks should be filled. 
    For example on 2 ranks at the initiation the whole array should be
    rank0: [abcd,0000]
    rank0: [0000,efgh]
    and after the broadcast
    rank0: [abcd,efgh]
    rank0: [abcd,efgh]
  */

  SECTION("broadcast"){
    cntr::distributed_array<double> A(nblock,blocksize,true);
    double err_loc=0.0;
    double err_glob=0.0;
    

    for(int j=0;j<nblock;j++){
      if(taskid==A.tid_map()[j]){
        for(int i=0;i<blocksize;i++){
          A.block(j)[i]=j+i*nblock;
        }
      }
    }
    
    for(int j=0;j<nblock;j++){
      A.mpi_bcast_block(j);
    }
     
    
	  for(int i=0;i<nblock;i++){
	    for(int j=0;j<blocksize;j++){
	      err_loc += (i+j*nblock-A.block(i)[j]);
	    }
    }
    
    MPI_Reduce(&err_loc,&err_glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD);
    REQUIRE(err_glob<eps);
  }

  /*
    The same test just using bcast all
  */
  SECTION("All-gather"){
    cntr::distributed_array<double> A(nblock,blocksize,true);
    double err_loc=0.0;
    double err_glob=0.0;
    

    for(int j=0;j<nblock;j++){
      if(taskid==A.tid_map()[j]){
        for(int i=0;i<blocksize;i++){
          A.block(j)[i]=j+i*nblock;
        }
      }
    }

    A.mpi_bcast_all();
    
    for(int i=0;i<nblock;i++){
      for(int j=0;j<blocksize;j++){
        err_loc += (i+j*nblock-A.block(i)[j]);
      }
    }
    
    MPI_Reduce(&err_loc,&err_glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD);
    REQUIRE(err_glob<eps);
  }

}