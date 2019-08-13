#ifndef CNTR_DISTRIBUTED_TIMESTEP_ARRAY_IMPL_H
#define CNTR_DISTRIBUTED_TIMESTEP_ARRAY_IMPL_H

#include "cntr_herm_matrix_timestep_view_decl.hpp"
#include "cntr_distributed_array_decl.hpp"
#include "cntr_distributed_timestep_array_decl.hpp"

namespace cntr {

/* #######################################################################################
#
#   CONSTRUCTION/DESTRUCTION
#
########################################################################################*/
template <typename T> distributed_timestep_array<T>::distributed_timestep_array(){
	n_=data_.n();
	tid_=data_.tid();
	ntasks_=data_.ntasks();
	G_=std::vector<cntr::herm_matrix_timestep_view<T> >(n_);
	tstp_=-2;
	nt_=-2;
	ntau_=-1;
	size_=0;
	sig_=-1;	
}
template <typename T> distributed_timestep_array<T>::~distributed_timestep_array(){
	// donothing
}
template <typename T> distributed_timestep_array<T>::distributed_timestep_array(const distributed_timestep_array &a){
	// default
	n_=a.n();
	tid_=a.tid();
	ntasks_=a.ntasks();
	data_=a.data(); // note: this copies the data of a
	tstp_=a.tstp();
	nt_=a.nt();
	ntau_=a.ntau();
	size_=a.size();
	sig_=a.sig();
	G_=a.G(); // this copies only the pointers of G_[j], they still point to the data of a
	// reset G_[j] to point to the new data
	for(int j=0;j<n_;j++) G_[j].set_to_data(data_.block(j),tstp_,ntau_,size_,sig_);
}
template <typename T> distributed_timestep_array<T>& distributed_timestep_array<T>::operator=(const distributed_timestep_array &a){
	if(this==&a) return *this;
	n_=a.n();
	tid_=a.tid();
	ntasks_=a.ntasks();
	data_=a.data(); // note: this copies the data of a
	tstp_=a.tstp();
	nt_=a.nt();
	ntau_=a.ntau();
	size_=a.size();
	sig_=a.sig();
	G_=a.G(); // this copies only the poiters of G_[j], they still point to the data of a
	// reset G_[j] to point to the new data
	for(int j=0;j<n_;j++) G_[j].set_to_data(data_.block(j),tstp_,ntau_,size_,sig_);
	return *this;
}

/** \brief <b> Initializes the `distributed_timestep_array` class.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Initializes the `distributed_timestep_array` class including n blocks of the class herm_matrix_ti* mestep_view, which serves as an interface to herm_matrix_timestep without copying the data. At   * the initialization the pointers don't point anywhere yet.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > Number of blocks
* @param nt
* > The maximum allowed number of timestep (set by size of data_)
* @param ntau
* > Number of points on Matsubara axis
* @param size
* > Matrix rank of the herm_matrix_timestep function
* @param sig
* > Set `sig = -1` for fermions or `sig = +1` for bosons
* @param mpi
* > If 'true' use MPI, otherwise one task with tid_=0 
*/
  
template <typename T> distributed_timestep_array<T>::distributed_timestep_array(int n,int nt,int ntau,int size,int sig,bool mpi){
	assert(-1<=nt && 0<=ntau &&  sig*sig==1 && 1<=size);
	int size_tstp=(ntau+1+2*(nt+1))*size*size;
	int maxlen=size_tstp;
	data_=cntr::distributed_array<std::complex<T> >(n,maxlen,mpi);
	n_=data_.n();	
	tid_=data_.tid();
	ntasks_=data_.ntasks();	
	tstp_=-2;
	nt_=nt;
	ntau_=ntau;
	size_=size;
	sig_=sig;
	// they all point nowhere, since tstp=-2
	G_=std::vector<cntr::herm_matrix_timestep_view<T> >(n_);
}

/** \brief <b> Reset the data to new timestep of herm_matrix_timestep  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Reset the data to new timestep which includes the change in the size and setting to the data at t* he given timestep 
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param tstp
* > New timestep to which the data is reset 
*/
  
template <typename T> void distributed_timestep_array<T>::reset_tstp(int tstp){
	assert(-1<=tstp && -1<=nt_);
	int blocksize=(ntau_+1+2*(tstp+1))*size_*size_;
	tstp_=tstp;
	data_.reset_blocksize(blocksize);
	for(int j=0;j<n_;j++) G_[j].set_to_data(data_.block(j),tstp_,ntau_,size_,sig_);
}

/** \brief <b> Clear the data   </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Clear the data  
*
* <!-- ARGUMENTS
*      ========= -->
*
*/
template <typename T> void distributed_timestep_array<T>::clear(void){
	data_.clear();
}

/** \brief <b> Get the pointer to the herm_matrix_timestep for the j-th block  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Return the pointer to the to the herm_matrix_timestep for the j-th block (=mpi rank)
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Block number (=mpi rank number) 
*/
  
template <typename T> cntr::herm_matrix_timestep_view<T>& distributed_timestep_array<T>::G(int j){
	assert(0<=j && 0<=n_-1);
	return G_[j];
}

/** \brief <b> MPI broadcast equivalent for the j-th block  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* MPI broadcast  equivalent of the distributed timestep array for the j-th block (=mpi rank)
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Block number (=mpi rank number) 
*/
  
template <typename T> void distributed_timestep_array<T>::mpi_bcast_block(int j){
	data_.mpi_bcast_block(j);
}

/** \brief <b> MPI allgather equivalent  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* MPI allgather  equivalent of the distributed timestep array 
* <!-- ARGUMENTS
*      ========= -->
*
*/
  
template <typename T> void distributed_timestep_array<T>::mpi_bcast_all(void){
	data_.mpi_bcast_all();
}

}
#endif
