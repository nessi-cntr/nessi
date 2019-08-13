

#include "catch.hpp"

#include <cmath>
#include <complex>
#include <cstring>
#include <iostream>
#include <sys/stat.h>

#include "cntr.hpp"

#define CPLX std::complex<double>
using namespace std;

TEST_CASE("Herm matrix hdf5 file read write", "[Herm_matrix_hdf5]") {

  int nt = 10;
  int ntau = 40;
  int size = 1;
  int sig = -1;
  double err;
  double eps = 1e-16;

  cntr::herm_matrix<double> G, G1;
  cntr::herm_matrix_timestep<double> Gtstp;
  cntr::herm_matrix_timestep_view<double> Gview;
  cntr::function<double> f, f1;

  G = cntr::herm_matrix<double>(nt, ntau, size, sig);
  cntr::green_equilibrium_bethe(G, 1.0, 0.1);

  {
    // Read entire Green's function from file

    G.write_to_hdf5("gbethe.h5", "G");
    G1.read_from_hdf5("gbethe.h5", "G");

    {
      double err = 0.;
      for (int tstp = -1; tstp <= nt; tstp++)
        err += cntr::distance_norm2(tstp, G, G1);
      REQUIRE(err < eps);
    }
  }

  {
    // Read from t=0 to read_nt from file

    G1.clear();
    int read_nt = 2;
    G1.read_from_hdf5(read_nt, "gbethe.h5", "G");

    {
      double err = 0.;
      for (int tstp = -1; tstp <= read_nt; tstp++)
        err += cntr::distance_norm2(tstp, G, G1);
      REQUIRE(err < eps);
    }
  }

  {
    // Write and read a single time step

    int select_timestep = 4;
    Gview = cntr::herm_matrix_timestep_view<double>(select_timestep, G);
    Gview.write_to_hdf5("gbethe_tstp4.h5", "G");
    Gtstp.read_from_hdf5("gbethe_tstp4.h5", "G");
    G1.set_timestep(Gtstp.tstp_, Gtstp);

    {
      double err = cntr::distance_norm2(select_timestep, G, G1);
      REQUIRE(err < eps);
    }
  }

  {
    // Write individual time slices out to disk
    // and read only time slice no 8

    G.write_to_hdf5_slices("g_slices.h5", "G", 2);

    hid_t file_id = read_hdf5_file("g_slices.h5");
    hid_t group_id = open_group(file_id, "G");
    hid_t sub_group_id = open_group(group_id, "t8");
    Gview = cntr::herm_matrix_timestep_view<double>(8, G1);
    Gview.read_from_hdf5(sub_group_id);
    close_group(sub_group_id);
    close_group(group_id);
    close_hdf5_file(file_id);

    {
      double err = cntr::distance_norm2(8, G, G1);
      REQUIRE(err < eps);
    }
  }

  /*
  {
  // The method append_to_hdf5_slices does not exist

    hid_t file_id = open_hdf5_file("test.out");
    G.append_to_hdf5_slices(file_id, -1);
    close_hdf5_file(file_id);
    file_id = H5Fopen("test.out", H5F_ACC_RDWR, H5P_DEFAULT);
    G.append_to_hdf5_slices(file_id, 8);
    G.append_to_hdf5_slices(file_id, 5);
    close_hdf5_file(file_id);
    file_id = H5Fopen("test.out", H5F_ACC_RDWR, H5P_DEFAULT);
    G.append_to_hdf5_slices(file_id, 8);
    close_hdf5_file(file_id);
  }
  */

  {
    // Write and read some contour functions to hdf5 file
    
    size = 2;
    f = cntr::function<double>(nt, size);
    for (int tstp = -1; tstp <= nt; tstp++) {
      f.ptr(tstp)[0] = 100 * tstp + 0;
      f.ptr(tstp)[1] = 100 * tstp + 1;
      f.ptr(tstp)[2] = 100 * tstp + 2;
      f.ptr(tstp)[3] = 100 * tstp + 3;
    }
    
    f.write_to_hdf5("func.h5", "f");
    f1.read_from_hdf5("func.h5", "f");

    {
      double err = 0.;
      for (int tstp = -1; tstp <= nt; tstp++)
        err += std::abs(f.ptr(tstp)[1] - f1.ptr(tstp)[1]);
      REQUIRE(err < eps);
    }
    
    f1 = cntr::function<double>(nt, size);

    int max_timestep = 4;
    f1.read_from_hdf5(max_timestep, "func.h5", "f");

    {
      double err = 0.;
      for (int tstp = -1; tstp <= max_timestep; tstp++)
        err += std::abs(f.ptr(tstp)[1] - f1.ptr(tstp)[1]);
      REQUIRE(err < eps);
    }
  }
  
}
