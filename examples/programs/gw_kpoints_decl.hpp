#pragma once

#include <sys/stat.h>
#include <complex>
#include "cntr/cntr.hpp"
#include "gw_latt_decl.hpp"

using namespace cntr;

#define CFUNC cntr::function<double>
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>


 /** \brief <b> Class `kpoint` includes propagators and correspoding routines for each momentum point.</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 *  Due to translational invariance all momentum points are independent and this 
 *  class includes information per k-point [instantaneous and retarded propagators] and
 *  corresponding solvers. 
 *  The communication between different k point takes place on higher level in the 
 *  class step.
 *
 */

namespace gw{

 class kpoint{
 public:
    kpoint();
    ~kpoint();
    kpoint(int nt,int ntau,int size,double beta,double h,double kk,lattice_1d_1b &latt,double mu);
    void init(int nt,int ntau,int size,double beta,double h,double kk,lattice_1d_1b &latt,double mu);
    void init_G_mat_nointeraction(lattice_1d_1b &latt,int SolverOrder);
    void set_hk(int tstp,int iter,lattice_1d_1b &latt);
    void set_vertex(int tstp,lattice_1d_1b &latt);    
    void write_to_hdf5(hid_t group_id);
    void write_to_hdf5(const char *filename);
    void write_to_hdf5_slices(hid_t group_id,int dt,int tid);
    void write_to_hdf5_slices(const char *filename,int dt,int tid);
    void get_Density_matrix(int tstp);
    void step_W(int tstp,int SolverOrder,lattice_1d_1b &latt);
    void step_W2b(int tstp,int SolverOrder,lattice_1d_1b &latt);
    void step_D(int tstp,int SolverOrder,lattice_1d_1b &latt);
    void step_dyson(int tstp,int iter,int SolverOrder,lattice_1d_1b &latt);
    double step_dyson_with_error(int tstp,int iter,int SolverOrder,lattice_1d_1b &latt);
    double step_W_with_error(int tstp,int iter, int n,int SolverOrder,lattice_1d_1b &latt);

    double beta_;
    double h_;
    int nt_;
    int ntau_;
    int nrpa_; // Dimension of the Green's function

    CFUNC rho_,SHartree_,SFock_,hk_,hkeff_,vertex_,hkeff_eigen_;
    GREEN G_,G0_,G0Sigma_,SigmaG0_,Sigma_,P_,VP_,PV_,W_;
    // For Hartree-Fock we need info of density/interactions for all k
    std::vector<CFUNC> density_k_;
    CFUNC rho_loc_;


    double mu_; 
    double kk_;
  };
}