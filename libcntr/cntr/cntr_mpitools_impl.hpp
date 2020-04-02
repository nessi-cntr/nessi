#ifndef CNTR_MPITOOLS_IMPL_H
#define CNTR_MPITOOLS_IMPL_H

#include "cntr_mpitools_decl.hpp"
#include "cntr_herm_matrix_decl.hpp"
#include "cntr_herm_matrix_timestep_decl.hpp"
#include "cntr_herm_matrix_timestep_view_decl.hpp"
#include "cntr_herm_matrix_timestep_view_impl.hpp"

namespace cntr{

/** \brief <b> MPI reduce for the `herm_matrix_timestep` to a `herm_matrix_timestep` </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > MPI reduce for the `herm_matrix_timestep` to the `root`
* > Works for scalar or square-matrix contour objects.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > time step
* @param root
* > Index of root
* @param Gred
* > The reduced `herm_matrix_timestep` on rank `root`.
* @param G
* > The `herm_matrix_timestep` on the individual ranks.
*/
template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix_timestep<T> &Gred, 
	herm_matrix_timestep<T> &G){
	assert(tstp == G.tstp()); 
	int taskid;
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	if (taskid == root) {
		assert(tstp == Gred.tstp());
		assert(G.ntau() == Gred.ntau());
		assert(G.size1() == Gred.size1());
		assert(G.size2() == Gred.size2());
	}

	int len = 2 * (2 * (tstp + 1) + G.ntau() + 1) * G.size1() * G.size2();

	if (sizeof(T) == sizeof(double)) {
		MPI_Reduce((double *)G.data_, (double *)Gred.data_, len, MPI_DOUBLE_PRECISION, MPI_SUM, root,
           MPI_COMM_WORLD);
   } else {
      if (taskid == root) std::cerr << "herm_matrix_timestep<T>::MPI_Reduce only for double " << std::endl;
      MPI_Finalize();
      exit(0);
   }

}


/** \brief <b> MPI reduce for the `herm_matrix_timestep_view` to a `herm_matrix_timestep_view` </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > MPI reduce for the `herm_matrix_timestep_view` to the `root`
* > Works for scalar or square-matrix contour objects.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > time step
* @param root
* > Index of root
* @param Gred
* > The reduced `herm_matrix_timestep_view` on rank `root`.
* @param G
* > The `herm_matrix_timestep_view` on the individual ranks.
*/
template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix_timestep_view<T> &Gred, 
	herm_matrix_timestep_view<T> &G){
	assert(tstp == G.tstp());
	int taskid;
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	if (taskid == root) {
		assert(tstp == Gred.tstp());
		assert(G.ntau() == Gred.ntau());
		assert(G.size1() == Gred.size1());
		assert(G.size2() == Gred.size2());
	}

	int len_rt = 2 * (tstp + 1) * G.size1() * G.size2();
	int len_it = 2 * (G.ntau() + 1) * G.size1() * G.size2();

	if (sizeof(T) == sizeof(double)) {
		if(tstp == -1){
			MPI_Reduce((double *)G.mat_, (double *)Gred.mat_, len_it, MPI_DOUBLE_PRECISION, MPI_SUM, root,
            	MPI_COMM_WORLD);
		} else{
			MPI_Reduce((double *)G.les_, (double *)Gred.les_, len_rt, MPI_DOUBLE_PRECISION, MPI_SUM, root,
            	MPI_COMM_WORLD);
			MPI_Reduce((double *)G.ret_, (double *)Gred.ret_, len_rt, MPI_DOUBLE_PRECISION, MPI_SUM, root,
            	MPI_COMM_WORLD);
			MPI_Reduce((double *)G.tv_, (double *)Gred.tv_, len_it, MPI_DOUBLE_PRECISION, MPI_SUM, root,
            	MPI_COMM_WORLD);
		}
		
   } else {
   	  if (taskid == root) std::cerr << "herm_matrix_timestep_view<T>::MPI_Reduce only for double " << std::endl;
   	  MPI_Finalize();
      exit(0);
   }

}


/** \brief <b> MPI reduce for the `herm_matrix_timestep` to a `herm_matrix` </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > MPI reduce for the `herm_matrix_timestep` to a `herm_matrix` on rank `root`.
* > Works for scalar or square-matrix contour objects.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > time step
* @param root
* > Index of root
* @param Gred
* > The reduced `herm_matrix` on rank `root`.
* @param G
* > The `herm_matrix_timestep` on the individual ranks.
*/
template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix<T> &Gred, 
	herm_matrix_timestep<T> &G){
	assert(tstp == G.tstp());
	int taskid;
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	if (taskid == root) {
		assert(tstp <= Gred.nt());
		assert(G.ntau() == Gred.ntau());
		assert(G.size1() == Gred.size1());
		assert(G.size2() == Gred.size2());
	}

    herm_matrix_timestep<T> Gtemp;
	if (taskid == root){
	    Gtemp.resize(tstp, G.ntau(), G.size1());
	}

	Reduce_timestep(tstp, root, Gtemp, G);

	if (taskid == root){
		Gred.set_timestep(tstp, Gtemp);
	}
}


/** \brief <b> MPI reduce for the `herm_matrix_timestep_view` to a `herm_matrix` </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > MPI reduce for the `herm_matrix_timestep_view` to a `herm_matrix` on rank `root`.
* > Works for scalar or square-matrix contour objects.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > time step
* @param root
* > Index of root
* @param Gred
* > The reduced `herm_matrix` on rank `root`.
* @param G
* > The `herm_matrix_timestep_view` on the individual ranks.
*/
template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix<T> &Gred, 
	herm_matrix_timestep_view<T> &G){
	assert(tstp == G.tstp());
	int taskid;
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	if (taskid == root) {
		assert(tstp <= Gred.nt());
		assert(G.ntau() == Gred.ntau());
		assert(G.size1() == Gred.size1());
		assert(G.size2() == Gred.size2());
	}

	bool check_tstp = (taskid == root);
	herm_matrix_timestep_view<T> Gred_tmp(tstp, Gred, check_tstp);
	Reduce_timestep(tstp, root, Gred_tmp, G);

}



/** \brief <b> MPI reduce for the `herm_matrix` to a `herm_matrix_timestep` </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > MPI reduce for the `herm_matrix` to a `herm_matrix_timestep` on rank `root`.
* > Works for scalar or square-matrix contour objects.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > time step
* @param root
* > Index of root
* @param Gred
* > The reduced `herm_matrix_timestep` on rank `root`.
* @param G
* > The `herm_matrix` on the individual ranks.
*/
template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix_timestep<T> &Gred, 
	herm_matrix<T> &G){
	assert(tstp <= G.nt());
	int taskid;
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	if (taskid == root) {
		assert(tstp == Gred.tstp());
		assert(G.ntau() == Gred.ntau());
		assert(G.size1() == Gred.size1());
		assert(G.size2() == Gred.size2());
	}

	herm_matrix_timestep<T> Gtemp;
	Gtemp.resize(tstp, G.ntau(), G.size1());
	G.get_timestep(tstp, Gtemp);

	Reduce_timestep(tstp, root, Gred, Gtemp);

}

/** \brief <b> MPI reduce for the `herm_matrix` to a `herm_matrix_timestep_view` </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > MPI reduce for the `herm_matrix` to a `herm_matrix_timestep_view` on rank `root`.
* > Works for scalar or square-matrix contour objects.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > time step
* @param root
* > Index of root
* @param Gred
* > The reduced `herm_matrix_timestep_view` on rank `root`.
* @param G
* > The `herm_matrix` on the individual ranks.
*/
template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix_timestep_view<T> &Gred, 
	herm_matrix<T> &G){
	assert(tstp <= G.nt());
	int taskid;
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	if (taskid == root) {
		assert(tstp == Gred.tstp());
		assert(G.ntau() == Gred.ntau());
		assert(G.size1() == Gred.size1());
		assert(G.size2() == Gred.size2());
	}

	herm_matrix_timestep_view<T> Gtemp(tstp, G);
	Reduce_timestep(tstp, root, Gred, Gtemp);

}


/** \brief <b> MPI reduce for the `herm_matrix` to a `herm_matrix` </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* > MPI reduce for the `herm_matrix` to a `herm_matrix` on rank `root`.
* > Works for scalar or square-matrix contour objects.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > time step
* @param root
* > Index of root
* @param Gred
* > The reduced `herm_matrix` on rank `root`.
* @param G
* > The `herm_matrix` on the individual ranks.
*/
template <typename T> void Reduce_timestep(int tstp, int root, herm_matrix<T> &Gred, 
	herm_matrix<T> &G){
	assert(tstp <= G.nt());
	int taskid;
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	if (taskid == root) {
		assert(tstp <= Gred.nt());
		assert(G.ntau() == Gred.ntau());
		assert(G.size1() == Gred.size1());
		assert(G.size2() == Gred.size2());
	}

	bool check_tstp = (taskid == root);
	herm_matrix_timestep_view<T> Gred_tmp(tstp, Gred, check_tstp);
	herm_matrix_timestep_view<T> Garr_tmp(tstp, G);

	Reduce_timestep(tstp, root, Gred_tmp, Garr_tmp);
}



} // namespace cntr


#endif // CNTR_MPITOOLS_IMPL_H