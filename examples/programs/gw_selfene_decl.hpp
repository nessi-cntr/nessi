#pragma once

#include <sys/stat.h>
#include <complex>
#include "cntr/cntr.hpp"
#include "gw_kpoints_decl.hpp"
#include "gw_latt_decl.hpp"

using namespace cntr;

#define CFUNC cntr::function<double>
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
#define DIST_TIMESTEP cntr::distributed_timestep_array<double>


namespace diag {
	void init_G_mat_nointeraction(int Nk_rank,std::vector<gw::kpoint> &corrK_rank,lattice_1d_1b &lattice,int SolverOrder);
	void gather_gk_timestep(int tstp,int Nk_rank,DIST_TIMESTEP &gk_all_timesteps, std::vector<gw::kpoint> &corrK_rank,
		std::vector<int> &kindex_rank);
	void gather_wk_timestep(int tstp,int Nk_rank,DIST_TIMESTEP &wk_all_timesteps, std::vector<gw::kpoint> &corrK_rank,
		std::vector<int> &kindex_rank);
	void set_density_k(int tstp, int Norb,DIST_TIMESTEP &gk_all_timesteps,lattice_1d_1b &lattice,std::vector<CFUNC> & density_k,std::vector<int> &kindex_rank,CFUNC &rho_loc);
	void get_loc(int tstp,int Ntau,int Norb,lattice_1d_1b &lattice,GREEN &Gloc,DIST_TIMESTEP &gk_all_timesteps);
	void sigma_Hartree(int tstp,int Norb,CFUNC &S,lattice_1d_1b &lattice,std::vector<CFUNC> &density_k,std::vector<CFUNC> &vertex,CFUNC &Ut);
	void sigma_Fock(int tstp,int Norb,int kk,CFUNC &S,lattice_1d_1b &lattice,std::vector<CFUNC> &density_k,std::vector<CFUNC> &vertex,CFUNC &Ut);
	void sigma_GW(int tstp,int kk,GREEN &S,DIST_TIMESTEP &gk_all_timesteps,DIST_TIMESTEP &wk_all_timesteps,lattice_1d_1b &lattice,int Ntau,int Norb);
	void get_Polarization_Bubble(int tstp,int Norb,int Ntau,int qq,GREEN &P,DIST_TIMESTEP &gk_all_timesteps,lattice_1d_1b &lattice);
	void symmetrise_mat(GREEN &G,int Ntau);
	void extrapolate_timestep_W(int tstp,int Nk_rank,int SolverOrder,int Nt,std::vector<gw::kpoint> &corrK_rank);
	void extrapolate_timestep_G(int tstp,int Nk_rank,int SolverOrder,int Nt,std::vector<gw::kpoint> &corrK_rank);
	void extrapolate_rho(int tstp,int Nk,std::vector<CFUNC> &density_k);
	double KineticEnergy(int tstp,lattice_1d_1b &lattice,std::vector<CFUNC> & density_k);
	double current(int tstp,lattice_1d_1b &lattice,std::vector<CFUNC> & density_k);
	double CorrelationEnergy(int tstp,int Nk_rank,int SolverOrder,double beta, double h,std::vector<gw::kpoint> &corrK_rank,std::vector<int> &kindex_rank,lattice_1d_1b &lattice);
}
