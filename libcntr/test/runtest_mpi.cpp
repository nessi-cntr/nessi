#include <mpi.h>
// Next line tells CATCH we will use our own main function
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "herm_matrix_mpi.hpp"
#include "distributed_array_mpi.hpp"
#include "distributed_timestep_array_mpi.hpp"


int main(int argc, char * argv[]) {
    int ierr;
    MPI_Init(&argc,&argv);
    int result = Catch::Session().run(argc, argv);
    MPI_Finalize();
    return result;
}