#ifndef CNTR_DISTRIBUTED_ARRAY_DECL_H
#define CNTR_DISTRIBUTED_ARRAY_DECL_H

#include "cntr_herm_matrix_timestep_view_decl.hpp"


namespace cntr{

/** \brief <b> Auxiliary data structure for handling 
set of data blocks and includes usual MPI processes on them.</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * Auxiliary data structure for handling of data blocks (total number is n_) ,
 * which are stored in  * contiguous form in member data_ and
 * includes usual MPI routines. The class identity is marked by * tid_  \f$ \in (0,\ldots,ntasks_-1) \f$
 * and the value of the tid_ is just the MPI rank or 0 if MPI is not defined. 
 * Each data block j is owned by precisely one process, which is given by \f$ tid * _map(j) = tid_\f$.
 * The member maxlen marks the maximum size reserved for the block. 
 * NOTE: even if the block is not owned by the process, the space for 
 * the data is allocated; the ownership plays a role when the data are 
 * manipulated
 */


template <typename T> 
class distributed_array{
public:
  /* construction, destruction */
  distributed_array();
  ~distributed_array();
  distributed_array(const distributed_array &g);
  distributed_array &operator=(const distributed_array &g);
  // #if __cplusplus >= 201103L
  // distributed_array(distributed_array &&g) noexcept;
  // distributed_array &operator=(distributed_array &&g) noexcept;
  // #endif
  distributed_array(int n,int maxlen,bool mpi);
  T* block(int j);
  void clear(void);
  void reset_blocksize(int blocksize);
  T* data(void) const {return data_;}
  int n(void) const {return n_;}
  int numblock_rank(void);
  int firstblock_rank(void);
  int blocksize(void) const {return blocksize_;}
  int maxlen(void) const {return maxlen_;}
  std::vector<int> tid_map(void) const {return tid_map_;}
  int tid(void) const {return tid_;}
  int ntasks(void) const {return ntasks_;}
  bool rank_owns(int k) const {return tid_map()[k] == tid();}
  
  // MPI UTILS
  // all MPI routines are given "trivial" no MPI versions, which basically assume ntasks=1 and do nothing
#if CNTR_USE_MPI==0
  void mpi_bcast_block(int j);
  void mpi_bcast_all(void);
  void mpi_gather(void);
#else  // CNTR_USE_MPI==1
  void mpi_send_block(int j,int dest);
  void mpi_gather(int dest);
  void mpi_bcast_block(int j);
  void mpi_bcast_all(void);
#endif  // CNTR_USE_MPI
private:
  T *data_;                  /*!< Pointer to the contiguous data */
  int n_;                    /*!< Number of data blocks */ 
  int blocksize_;            /*!< Size of one element in the block. In combination with distributed_timestep_array this number is changed with increasing timesteps. */ 
  int maxlen_;               /*!< Size of one element in the block. In combination with distributed_timestep_array this number is fixed to maximal length */ 
  std::vector<int> tid_map_; /*!< Returns the rank owning block j \f$ tid_map(j) = tid_\f$. */
  int tid_;                  /*!< mpi rank if MPI is defined, else 0 */
  int ntasks_;               /*!< mpi size if MPI is defined, else 1 */
};

} //namespace cntr
#endif
