#ifndef CNTR_DISTRIBUTED_TIMESTEP_ARRAY_MOVING_DECL_H
#define CNTR_DISTRIBUTED_TIMESTEP_ARRAY_MOVING_DECL_H

#include "cntr_distributed_array_decl.hpp"
#include "cntr_herm_matrix_timestep_moving_view_decl.hpp"

namespace cntr{

/** \brief <b> Specialization of the distributed_array in which data-blocks are associated with the herm_matrix_timestep_moving objects </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * Specialization of the distributed_array for the herm_matrix_timestep_moving objects, which is
 * used for problems, where all ranks need to have the full information [for instance spatial] about the system for a given timestep.
 * In practice the time stepping procedure is used and the last
 * timestep needs to be communicated between all  MPI processes.
 */

  template <typename T>
  class distributed_timestep_array_moving{
  public:

    /* construction,destruction  */
    distributed_timestep_array_moving();
    ~distributed_timestep_array_moving();
    distributed_timestep_array_moving(const distributed_timestep_array_moving &a);
    distributed_timestep_array_moving<T> & operator=(const distributed_timestep_array_moving &a);
    distributed_timestep_array_moving(int n,int tc,int t0,int size,int sig,bool mpi);
    // #if __cplusplus >= 201103L
    //   		distributed_timestep_array_moving(distributed_timestep_array_moving &&a) noexcept;
    //   		distributed_timestep_array_moving &operator=(distributed_timestep_array_moving &&a) noexcept;
    // #endif
    ////////////////////////////////////////////
    void reset_tstp(int tstp);
    void clear(void); // set all data to zero
    ////////////////////////////////////////////
    // access:
    cntr::herm_matrix_timestep_moving_view<T> &G(int j);
    ////////////////////////////////////////////
    // MPI collective distribution routines
    distributed_array<std::complex<T> > data(void) const {return data_;}
    int n(void) const {return n_;}
    int tid(void) const {return tid_;}
    int ntasks(void) const {return ntasks_;}
    std::vector<cntr::herm_matrix_timestep_moving_view<T> > G(void) const {return G_;}
    int t0(void) const {return t0_;}
    int tc(void) const {return tc_;}
    int size(void) const {return size_;}
    int sig(void) const {return sig_;}
    void set_t0(int t0) {t0_=t0;}
    void mpi_bcast_block(int j);
    void mpi_bcast_all(void);

  private:
    distributed_array<std::complex<T> > data_;          /*!< Pointer to the contiguous data */
    int n_;                                             /*!< Number of blocks */
    int tid_;                                           /*!< MPI rank if MPI is defined, else 0 */
    int ntasks_;                                        /*!< MPI size if MPI is defined, else 1 */
    std::vector<cntr::herm_matrix_timestep_moving_view<T> > G_;/*!< Vector of views to the herm_matrix_timestep */
    int t0_;                                            /*!< Current timestep */
    int tc_;                                            /*!< The maximum allowed timestep (set by size of data_) */
    int size_;                                          /*!< Size of green's function*/
    int sig_;                                           /*!< Fermion(boson)=-1(1)*/
  };

} //namespace cntr
#endif
