#ifndef CNTR_MPITOOLS_DECL_H
#define CNTR_MPITOOLS_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

#if CNTR_USE_MPI == 1
	template <typename T> class herm_matrix;
	template <typename T> class herm_matrix_timestep;
	template <typename T> class herm_matrix_timestep_view;

	template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix_timestep<T> &Gred, 
		herm_matrix_timestep<T> &G);
	template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix_timestep_view<T> &Gred, 
		herm_matrix_timestep_view<T> &G);
	template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix<T> &Gred, 
		herm_matrix_timestep<T> &G);
	template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix<T> &Gred, 
		herm_matrix_timestep_view<T> &G);
	template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix_timestep<T> &Gred, 
		herm_matrix<T> &G);
	template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix_timestep_view<T> &Gred, 
		herm_matrix<T> &G);
	template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix<T> &Gred, 
		herm_matrix<T> &G);
#endif 
}

#endif // CNTR_MPITOOLS_DECL_H