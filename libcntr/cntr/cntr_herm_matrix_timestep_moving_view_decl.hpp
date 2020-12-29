#ifndef CNTR_HERM_TIMESTEP_MOVING_VIEW_DECL_H
#define CNTR_HERM_TIMESTEP_MOVING_VIEW_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

  template <typename T> class herm_matrix_timestep;
  template <typename T> class herm_matrix;
  template <typename T> class herm_matrix_timestep_moving;
  template <typename T> class herm_matrix_moving;
  template <typename T>
  /** \brief <b> Class `herm_matrix_timestep_trunc_view` serves for interfacing with class herm_matrix_timestep_trunc</b>
   * without copying the data.
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *  Occasionally you may have a pre-defined object of herm_matrix_trunc or herm_matrix_timestep_trunc that you want to use.
   *  While one option is to make a copy of the data, most commonly you probably want to re-use this memory as a libcntr type,
   *  and this class serves as an interface. Mark similar concept of a map class in Eigen.
   *  To construct this class you need: a pointers to the region of memory defining different components (`ret`,`tv`,`les`,`mat`)
   *  and usual defining constants for the herm_matrix_timestep.
   *  The contour function can be of scalar type or matrix-valued.
   *
   *  Square-matrix-type contour functions are fully supported, while non-square-
   *  matrix contour functions are under development.
   *
   *
   */
  /// !! NOTE: MIND DANGLING POINTERS

  class herm_matrix_timestep_moving_view {
  public:
    typedef std::complex<T> cplx;
    typedef T scalar_type;
    herm_matrix_timestep_moving_view();
    ~herm_matrix_timestep_moving_view();
    herm_matrix_timestep_moving_view(const herm_matrix_timestep_moving_view &g);
    herm_matrix_timestep_moving_view(herm_matrix_timestep_moving<T> &g);
    herm_matrix_timestep_moving_view &operator=(const herm_matrix_timestep_moving_view &g);
    herm_matrix_timestep_moving_view(int tstp, herm_matrix_moving<T> &g);
    herm_matrix_timestep_moving_view(int t0, herm_matrix_timestep_moving<T> &g);
    herm_matrix_timestep_moving_view(int t0, herm_matrix_timestep_moving_view<T> &g);
    herm_matrix_timestep_moving_view(int tc, int t0, int size1, int size2, int sig);
    ///////////////////////////////////////////
    int size1(void) const { return size1_; }
    int size2(void) const { return size2_; }
    int sig(void) const { return sig_; }
    int tc(void) const { return tc_; }
    int t0(void) const { return t0_; }
    inline std::complex<T> *retptr(int i) { return ret_ + i * element_size_; }
    inline std::complex<T> *lesptr(int i) {
      return les_ + i * element_size_;
    };
    ///////////////////////////////////////////
    // the following routines set the pointers:
    void set_to_data(std::complex<T> *data, int tc, int t0, int size, int sig);
    void set_to_data(int tstp, herm_matrix_moving<T> &g);
    void set_to_data(herm_matrix_timestep_moving<T> &g);
    //TO-DO
    // HDF5 I/O
    // #if CNTR_USE_HDF5 == 1
    //     void write_to_hdf5(hid_t group_id);
    //     void write_to_hdf5(hid_t group_id, const char *groupname);
    //     void write_to_hdf5(const char *filename, const char *groupname);
    //     // READ: Read only data, no resize!
    //     void read_from_hdf5(hid_t group_id, const char *groupname);
    //     void read_from_hdf5(hid_t group_id);
    //     void read_from_hdf5(const char *filename, const char *groupname);
    // #endif

    /*//////////////////////////////////////////////////////////////////////////////////
      OPERATE ON THE DATA POINTED TO BY THE VIEW OBJECT:

      the argument G can be anything on which one can create a temporary
      herm_matrix_timestep_moving_view object by
      herm_matrix_timestep_moving_view<T>(tstp_,G)

    //////////////////////////////////////////////////////////////////////////////////*/
    // ** GET_DATA: copy data into the array pointed to by the view:
    void get_data(herm_matrix_timestep_moving_view<T> &g);
    template <class GG>
    void get_data(GG &g);
    template <class GG>
    void get_data(int tstp,GG &g);
    //void get_timestep(int tstp,herm_matrix_timestep<T> &timestep);
    // ** SET MATRIX_ELEMENT:  this(i1,i2) <= g(j1,j2)
    void set_matrixelement(int i1, int i2, herm_matrix_timestep_moving_view<T> &g,
                           int j1, int j2);
    template <class GG>
    void set_matrixelement(int i1, int i2, GG &g, int j1, int j2);
    // ** INCR_TIMESTEP  this += alpha * g
    // TO-DO
    // ** Would be more more consistent with description, when read only... **
    // void incr_timestep(herm_matrix_timestep_moving_view<T> &g,
    //                    std::complex<T> alpha);
    // void incr_timestep(herm_matrix_timestep_moving_view<T> &g, T alpha);
    // template <class GG>
    // void incr_timestep(GG &g, std::complex<T> alpha);
    // template <class GG>
    // void incr_timestep(GG &g, T alpha);
    // // ** SMUL  this *= alpha
    // void smul(T alpha);
#if CNTR_USE_MPI == 1
    void MPI_Reduce(int root);
#endif

    /////////////////////////////////////////////////////////////////////////
    std::complex<T> density_matrix(int tstp);
    template <class Matrix>
    void density_matrix(int tstp, Matrix &M);
    /////////////////////////////////////////////////////////////////////////
    // private:
    void get_data(std::complex<T> *ret, std::complex<T> *les);
    /// @private
    /** \brief <b> Pointer to the retarded component. 'ret_+ \f$t \;*\f$ element\_size'  corresponds to \f$(0,0)\f$-component of \f$ G^R(tstp,t) \f$ </b> */
    std::complex<T> *ret_;
    /// @private
    /** \brief <b> Pointer to the lesser component. 'less_+ \f$t \;*\f$ element\_size'  corresponds to \f$(0,0)\f$-component of \f$ G^<(t,tstp) \f$ </b> */
    std::complex<T> *les_;
    /// @private
    /** \brief <b> time cut-off.</b> */
    int tc_;
    /// @private
    /** \brief <b> time step.</b> */
    int t0_;
    /// @private
    /** \brief <b> Number of the colums in the Matrix form.</b> */
    int size1_;
    /// @private
    /** \brief <b> Number of the rows in the Matrix form.</b> */
    int size2_;
    /// @private
    /** \brief <b> Size of the Matrix form; size1*size2. </b> */
    int element_size_;
    /// @private
    /** \brief <b> Bose = +1, Fermi =-1. </b> */
    int sig_;
  };

}  // namespace cntr

#endif  // CNTR_HERM_TIMESTEP_TRUNC_VIEW_DECL_H
