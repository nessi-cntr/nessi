#pragma once

#include <sys/stat.h>
#include <complex>
#include "cntr/cntr.hpp"
#include "cntr/hdf5/hdf5_interface.hpp"
#include "cntr/hdf5/hdf5_interface_cntr.hpp"

#define CFUNC cntr::function<double>
#define GREEN cntr::herm_matrix<double>
#define GREEN_TSTP cntr::herm_matrix_timestep<double>
 


class lattice_1d_1b{
/*///////////////////////////////////////////////////////////////////////////

1d-lattice, 1 band


setup of kpoints: 

** The model stays inversion symmetric (k <-> k-1); (there is no electric 
field);
   
** each point in the BZ has one representatives under the point-group symmetry;
   - kpoints_[0...nk_-1] is the list of all representatives (we choose kk<=0)
   - later, Green functions are only stored for representatives


G=(nk-1),  (X]: point X included, [X: point X excluded,).

       kk             -PI   0     PI
 	
index k (full BZ)      0    G    [2G
index k (stored)       0    G ]

mappings: 
BZ -> stored:
representative_kk  (kbz<=G ? kbz : 2G-kbz) 

///////////////////////////////////////////////////////////////////////////*/
public:
	int nt_;    
	int nk_; // the RBZ  
	int G_; // = index of Gamma point in BZ and stored sector

	std::vector<double>  kpoints_;  // 0...nk_-1
	std::vector<double>  kweight_;  // 0...nk_-1, normalized to 1
	//////////////////////////////////////////////////////////////////////
	// MODEL PARAMETERS
	int Norb_;
	double mu_;
	CFUNC   tt_; // time-dependent hopping
	CFUNC   U_,V_,Apulse_;

	lattice_1d_1b(void);
	lattice_1d_1b(int nk,int nt,CFUNC &tt,CFUNC &Ut,CFUNC &Vt,std::vector<double> &Epulse,double mu,int Norb,int SolverOrder,double h);
	int representative_kk(int kbz);
	int add_kpoints(int k1,int s1,int k2,int s2);
	void init_kk(int nk);
	// the dispersion (particle-hole symmetric for A=0)
	void hk(cdmatrix &hkmatrix,int tstp,double kk,int iter);
	void hkfree(cdmatrix &hkmatrix,int tstp,double kk);
	void efield_to_afield(int nt,double h,std::vector<double> &efield,int SolverOrder);
	//Velocity is still 1d vector
	void vk(cdmatrix &vkmatrix,int tstp,double kk);
	void V(int tstp,double qq,cdmatrix &V);
};
