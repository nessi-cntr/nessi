/* 
 * Copyright: D. Golez (2018), The Simons Foundation (2019)
 * Authors: D. Golez, H. UR Strand
 */

#ifndef CNTR_DISTRIBUTED_ARRAY_IMPL_H
#define CNTR_DISTRIBUTED_ARRAY_IMPL_H

#include "cntr_herm_matrix_timestep_view_decl.hpp"
#include "cntr_distributed_array_decl.hpp"


namespace cntr {
/* #######################################################################################
#
#   CONSTRUCTION/DESTRUCTION
#
########################################################################################*/
template <typename T> 
distributed_array<T>::distributed_array(){
	data_=0;
	maxlen_=0;
	blocksize_=0;
	n_=0;
	tid_=0;
	ntasks_=1;
	tid_map_=std::vector<int>(0);
}
template <typename T> distributed_array<T>::~distributed_array(){ 
   if(data_!=0) delete [] data_;
}
template <typename T> distributed_array<T>::distributed_array(const distributed_array &g){
	size_t len;
	blocksize_=g.blocksize_;
	n_=g.n();
	tid_=g.tid();
	ntasks_=g.ntasks();
	tid_map_=g.tid_map();
	maxlen_=g.maxlen();
	len=maxlen_*n_;
	if(len>0){
	  data_ = new T [len];
	  memcpy(data_, g.data(), sizeof(T)*len);
	}else{
	   data_=0;
	}
}
template <typename T> distributed_array<T> & distributed_array<T>::operator=(const distributed_array &g){
 	size_t len;
	if(this==&g) return *this;
	if(data_!=0) delete [] data_;
	blocksize_=g.blocksize();
	n_=g.n();
	tid_=g.tid();
	ntasks_=g.ntasks();
	tid_map_=g.tid_map();
	maxlen_=g.maxlen();
	len=maxlen_*n_;
	if(len>0){
		data_ = new T [len];
		memcpy(data_, g.data(), sizeof(T)*len);
	}else{
		data_=0;
	}
	return *this;
}

/** \brief <b> Initializes the `distributed_array` class.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Initializes the `distributed_array` class, where template T is a type of the basic unit of data,  with or without the MPI support.
*
* <!-- ARGUMENTS
*      ========= -->
*
* @param n
* > Number of blocks
* @param maxlen
* > Maximum size of block 
* @param mpi
* > If 'true' use MPI, otherwise one task with tid_=0 
*/

  template <typename T> distributed_array<T>::distributed_array(int n,int maxlen,bool mpi){
	assert(0<=maxlen && 0<=n);
	size_t len;
	maxlen_=maxlen;
	blocksize_ = maxlen_;
	n_=n;
	tid_map_.resize(n_);
	len=maxlen_*n_;
	if(len==0){
		data_=0;
	}else{ 
		data_ = new T [len];
		memset(data_, 0, sizeof(T)*len);
	}
	
	if(mpi){
		#if CNTR_USE_MPI==1
			ntasks_=MPI::COMM_WORLD.Get_size();
			tid_=MPI::COMM_WORLD.Get_rank();
		#else
			tid_=0;
			ntasks_=1;		
		#endif	
	}else{
		tid_=0;
		ntasks_=1;		
	}
	{
		// This distribution tries to spread evenly number of tasks
		int nr_alloced = 0;
		int remainer,buckets;
		for(int i=0;i<ntasks_;i++){
			remainer = n - nr_alloced;
			buckets = (ntasks_ - i);
			int size=remainer/ buckets;
			for(int k=0;k<size;k++){
				tid_map_[nr_alloced+k]=i;
			}
			nr_alloced +=size; //Ceiling division
		}


		// int rank_totblock=0,rank_firstblock=0;
		// // Set position of first block on a given rank
		// for(int i=0;i<n_;i++){
		// 	if(tid_map_[i]==tid_){
		// 		rank_firstblock=i;
		// 		break;
		// 	}
		// }
		// // Set number of all blocks on a given rank
		// for(int i=0;i<n_;i++){
		// 	if(tid_map_[i]==tid_){
		// 		rank_totblock+=1;
		// 	}
		// }

		// for(int i=0;i<n_;i++){
		// 	std::cout << "Tid map " << tid_ << " " << i << " " << tid_map_[i] << " " <<  rank_totblock << " " << rank_firstblock << std::endl;	
		// }

		// std::cout << " -------------------------------------- " << std::endl;

	}
}

  /** \brief <b> Returns the pointer to block j.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Returns the pointer to block j.
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Block index .
*/
  
template <typename T> T* distributed_array<T>::block(int j){
	assert(j<=n_-1);
	return data_+j*blocksize_;
}

  /** \brief <b> Clear all data  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Clear all data
*/
  
template <typename T> void distributed_array<T>::clear(void){
	if(data_!=0) memset(data_, 0, sizeof(T)*maxlen_*n_);
}

    /** \brief <b> Reset the block size for each block.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* > Clear the previous data and reset the block size for each block.
* <!-- ARGUMENTS
*      ========= -->
*
* @param blocksize
* > New size of the block
*/
  
template <typename T> void distributed_array<T>::reset_blocksize(int blocksize){
	assert(0<=blocksize && 0<=blocksize_ && 0<= (int) maxlen_);
	clear();
	blocksize_=blocksize;
}

/** \brief <b> Return number of blocks on rank.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* >  Return number of blocks on rank. 
*/
  
template <typename T> int distributed_array<T>::numblock_rank(void){
	int rank_totblock=0;

	for(int i=0;i<n_;i++){
		if(tid_map_[i]==tid_){
			rank_totblock+=1;
		}
	}

	return rank_totblock;
}


/** \brief <b> Return first blocks on rank.  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
*
* >  Return first blocks on rank. 
*/
  
template <typename T> int distributed_array<T>::firstblock_rank(void){
	int rank_firstblock=0;

	for(int i=0;i<n_;i++){
		if(tid_map_[i]==tid_){
			rank_firstblock=i;
			break;
		}
	}
	return rank_firstblock;
}

/* #######################################################################################
#
#   MPI UTILS
#
########################################################################################*/
// In case you don't use MPI all routines are given "trivial" no MPI versions, assume ntasks=1 and do nothing
#if CNTR_USE_MPI==0
template <typename T> void distributed_array<T>::mpi_bcast_block(int j){
    // donothing
}
template <typename T> void distributed_array<T>::mpi_bcast_all(void){
    // donothing
}
template <typename T> void distributed_array<T>::mpi_gather(void){
    // donothing
}
#else  // CNTR_USE_MPI==1
/** \brief <b> Sends the j-th block to the MPI rank dest  </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* Sends the j-th block to the MPI rank (=tid_ of the block) dest
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Index of the block send
* @param dest
* > Index of the MPI rank (=tid_ of the block) receiving the data
*/
template <typename T> void distributed_array<T>::mpi_send_block(int j,int dest){
	assert(0<=dest && 0<=ntasks_-1);
	int root=tid_map_[j];
	if(dest==root){
		return;
	}else{
		int tag=100;
		size_t int_per_t=sizeof(T)/sizeof(int);
		size_t len=int_per_t*blocksize_;
		if(tid_map_[j]==tid_){
			MPI_Send((int*) block(j), len, MPI_INT, dest, tag, MPI_COMM_WORLD);
			// MPI::COMM_WORLD.Send((int*) block(j),len,MPI::INT,dest,tag);
		}
		if(tid_==dest){
			// MPI::COMM_WORLD.Recv((int*) block(j),len,MPI::INT,root,tag);
			MPI_Recv((int*) block(j), len, MPI_INT, root, tag, MPI_COMM_WORLD, 
				MPI_STATUS_IGNORE);
		}
	}
}
/** \brief <b> MPI gather for the distributed array   </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* MPI gather equivalent for the distributed array, where we gather the data on rank dest 
* <!-- ARGUMENTS
*      ========= -->
*
* @param dest
* > Index of the MPI rank (=tid_ of the block) receiving the data
*/
template <typename T> void distributed_array<T>::mpi_gather(int dest){
	assert(0<=dest && 0<=ntasks_-1);
	for(int j=0;j<n_;j++) mpi_send_block(j,dest);
}

/** \brief <b> MPI broadcast of j-th block   </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* MPI broadcast equivalent for the distributed array, where we broadcast the j-th block 
* <!-- ARGUMENTS
*      ========= -->
*
* @param j
* > Index of the MPI rank (=tid_ of the block), which we want to broadcast
*/
template <typename T> void distributed_array<T>::mpi_bcast_block(int j){
	int root=tid_map_[j];
	size_t int_per_t=sizeof(T)/sizeof(int);
	assert(int_per_t*sizeof(int)==sizeof(T));
	size_t len=int_per_t*blocksize_;
	// MPI::COMM_WORLD.Bcast((int*)block(j),len,MPI::INT,root);
	MPI_Bcast((int*)block(j), len, MPI_INT, root, MPI_COMM_WORLD);
}

//* in a global allgather operation, the data are send from the root to all 
//  other processes

/** \brief <b> MPI Allgather equivalent for the distributed array   </b>
*
* <!-- ====== DOCUMENTATION ====== -->
*
*  \par Purpose
* <!-- ========= -->
*
* MPI Allgather equivalent for the distributed array, which is used when all
* processes needs to aggregate the data.
* <!-- ARGUMENTS
*      ========= -->
*
*/
template <typename T> void distributed_array<T>::mpi_bcast_all(void){

  size_t int_per_t = sizeof(T) / sizeof(int);
  assert(int_per_t * sizeof(int) == sizeof(T));
  int element_size = blocksize_ * int_per_t;
	
  std::vector<int> recvcount(ntasks_, 0);
  for(int j=0;j<n_;j++) {
    int tid_ = tid_map_[j];
    recvcount[tid_] += 1;
  }
  
  std::vector<int> displs(ntasks_, 0);
  for(int rank = 1; rank < ntasks_; rank++) {
    displs[rank] = displs[rank - 1] + recvcount[rank - 1];
  }
  
  for(int rank = 0; rank < ntasks_; rank++) {
    recvcount[rank] *= element_size;
    displs[rank] *= element_size;
  }
	
  // MPI::COMM_WORLD.Allgatherv(MPI::IN_PLACE, 0, MPI::DATATYPE_NULL,
		// 	     data_, recvcount.data(), displs.data(), MPI::INT);

  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, data_, 
  	recvcount.data(), displs.data(), MPI_INT, MPI_COMM_WORLD);
	
}
#endif // CNTR_USE_MPI


}
#endif
