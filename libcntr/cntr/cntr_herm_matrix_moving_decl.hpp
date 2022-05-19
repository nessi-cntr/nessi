#ifndef CNTR_HERM_MATRIX_MOVING_DECL_H
#define CNTR_HERM_MATRIX_MOVING_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {
  /*//////////////////////////////////////////////////////////

    Internal storage order such that
    G.ret_[t1] + t2*element_size_  -> Gret(t0-t1,t0-t1-t2)
    G.les_[t1] + t2*element_size_  -> Gles(t0-t1,t0-t1-t2)

  ///////////////////////////////////////////////////////////*/

  template <typename T> class function;
  template <typename T> class herm_matrix_timestep;
  template <typename T> class herm_matrix_timestep_moving;
  template <typename T> class herm_matrix_timestep_view;
  template <typename T> class herm_matrix_timestep_moving;
  template <typename T> class function_moving;
  template <typename T> class herm_pseudo;
  
/** \brief <b> Class `herm_matrix_moving` for two-time contour objects \f$ C(t,t') \f$
   * stored for a reduced range of times.</b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   *  This class is a modification of `herm_matrix` where data is stored in the time
   *  range \f$ t_0\geq t>=t_0-t_c \f$. 
   *  The class `herm_matrix_moving` stores two non-redundant Keldysh components
   *   - retarded component \f$ C^\mathrm{R}(t_i,t_i-s) \f$ for \f$ s=0,\ldots,t_c\f$
   *   - lesser component \f$ C^<(t_i,t_i-s) \f$ for \f$ s=0,\ldots,t_c.\f$
   *
   *  The contour function \f$ C(t,t') \f$ can be of scalar type or matrix-valued.
   *  Square-matrix-type contour functions are fully supported, while non-square-
   *  matrix contour functions are under development.
   *
   *
   */
  template <typename T> class herm_matrix_moving{
  public:
    typedef std::complex<T> cplx;
    typedef T scalar_type;
    /* construction, destruction */
    herm_matrix_moving();
    ~herm_matrix_moving();
    herm_matrix_moving(int tc,int t0,int size1=1,int sig=-1);
    herm_matrix_moving(const herm_matrix_moving &g);
    herm_matrix_moving & operator=(const herm_matrix_moving &g);
    void clear(void);
    void resize(int tc,int size1);
    /* access size etc ... */
    /// @private
    int element_size(void) const{ return element_size_;}
    int size1(void) const{ return size1_;}
    int size2(void) const{ return size2_;}
    int tc(void) const{ return tc_;}
    int t0(void) const{ return t0_;}
    int sig(void) const{ return sig_;}
    void set_t0(int t0);
    // raw pointer to elements ... to be used with care
    /// @private
    inline cplx * lesptr(int i,int j){return les_[i] + j*element_size_;}  // points to Gret(t0-i,t0-i-j)
    /// @private
    inline cplx * retptr(int i,int j){return ret_[i] + j*element_size_;}  // points to Gles(t0-i,t0-i-j)
    // reading basic and derived elements to any Matrix type:
    // the get_... address time-arguments "relative to t0"  (i,j) == (t0-i,t0-i-j)
    // and work only for 0 <= i,j <= tc
    template<class Matrix> void get_les(int i,int j,Matrix &M) const;
    template<class Matrix> void get_gtr(int i,int j,Matrix &M) const;
    template<class Matrix> void get_ret(int i,int j,Matrix &M) const;
    // these will adress only (0,0) element for dim>1:
    inline void get_les(int i,int j,cplx &x) const;
    inline void get_gtr(int i,int j,cplx &x) const;
    inline void get_ret(int i,int j,cplx &x) const;
    cplx density_matrix(int i);  //  -sig*ii*Gles(t0-i,t0-i)
    template<class Matrix> void density_matrix(int tstp,Matrix &M);
    // writing basic elements (also relative to t0)
    /// @private
    template<class Matrix> void set_les(int i,int j,Matrix &M);
    /// @private
    template<class Matrix> void set_ret(int i,int j,Matrix &M);
    /// @private
    inline void set_les(int i,int j,cplx x);
    /// @private
    inline void set_ret(int i,int j,cplx x);
    // INPUT/OUTPUT
    void print_to_file(const char *file,int precision=16);
    void read_from_file(const char *file);
    // DATA exchange with HERM_MATRIX
    // read data to slice i (relative to t0)
    void clear_timestep(int i);
    void set_timestep(int i,herm_matrix_timestep_view<T> &g,herm_matrix_timestep_view<T> &gcc);
    void set_timestep(int i,int tstp,herm_matrix<T> &g,herm_matrix<T> &gcc);
    void set_from_G_backward(int tstp,herm_matrix<T> &g,herm_matrix<T> &gcc);
    void set_from_G_backward(int tstp,herm_matrix<T> &g);
    void set_from_G_backward(int tstp, const herm_pseudo<T> &g);
    void get_timestep(int i,herm_matrix_timestep_moving<T> &g);

    void set_timestep(int i,herm_matrix_timestep_moving<T> &g);
    void set_timestep(int i,herm_matrix_moving<T> &g,int j);
    void incr_timestep(int i,herm_matrix_timestep_moving<T> &g,cplx alpha);
    void incr_timestep(int i,herm_matrix_moving<T> &g,int j,cplx alpha);
    // only for timestep 0
    void left_multiply(function_moving<T> &g,T weight, int i=0);
    void right_multiply(function_moving<T> &g,T weight, int i=0);
    ///
    void smul(int tstp, T weight);
    void smul(int tstp, cplx weight);

    void forward(void);
#if CNTR_USE_HDF5 == 1
    void write_to_hdf5(hid_t group_id);
    void write_to_hdf5(hid_t group_id, const char *groupname);
    void write_to_hdf5(const char *filename, const char *groupname);
    void read_from_hdf5(hid_t group_id);
    void read_from_hdf5(hid_t group_id, const char *groupname);
    // TODO 
    // void read_from_hdf5(const char *filename, const char *groupname);
    // read up to timestep nt1: NO RESIZE
    // void read_from_hdf5(int nt1, hid_t group_id);
    // void read_from_hdf5(int nt1, hid_t group_id, const char *groupname);
    // void read_from_hdf5(int nt1, const char *filename, const char *groupname);
    // writing some other compressed format: just a collection of timeslices
    // -1,0,dt,2*dt,...
    // void write_to_hdf5_slices(hid_t group_id, int dt);
    // void write_to_hdf5_slices(hid_t group_id, const char *groupname, int dt);
    // void write_to_hdf5_slices(const char *filename, const char *groupname,
    // int dt);
    // void write_to_hdf5_tavtrel(hid_t group_id, int dt);
    // void write_to_hdf5_tavtrel(hid_t group_id, const char *groupname, int dt);
    // void write_to_hdf5_tavtrel(const char *filename, const char *groupname,
    //                            int dt);
#endif


  private:
    /** \brief <b> Pointer to the data array stoirng lesser and retarded component.</b> */
    cplx* data_;
    /** \brief <b> Pointer to a pointer storing the lesser component.</b> */
    cplx** les_;
    /** \brief <b> Pointer to a pointer storing the retarded component.</b> */
    cplx** ret_;
    /** \brief <b> Cutoff time.</b> */
    int tc_;
    /** \brief <b> Physical time.</b> */
    int t0_;
    /** \brief <b> Number of the colums in the Matrix form.</b> */
    int size1_;
    /** \brief <b> Number of the rows in the Matrix form.</b> */
    int size2_;
    /** \brief <b> Size of the Matrix form; size1*size2. </b> */
    int element_size_;
    /** \brief <b> Bose = +1, Fermi =-1. </b> */
    int sig_; // Bose = +1, Fermi =-1
  };


}  // namespace cntr

#endif  // CNTR_HERM_MATRIX_MOVING_DECL_H
