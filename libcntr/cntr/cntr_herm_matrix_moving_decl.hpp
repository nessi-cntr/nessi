#ifndef CNTR_HERM_MATRIX_MOVING_DECL_H
#define CNTR_HERM_MATRIX_MOVING_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

  template <typename T> class function;
  template <typename T> class herm_matrix_timestep;
  template <typename T> class herm_matrix_timestep_moving;
  template <typename T> class herm_matrix_timestep_view;
  template <typename T> class herm_matrix_timestep_moving;
  template <typename T> class herm_matrix_timestep_moving_view;
  template <typename T> class function_moving;
  template <typename T> class herm_pseudo;
  
/**
 * @brief <b> Class `herm_matrix_moving` for two-time contour objects \f$ C(t,t') \f$ on a memory-truncated window.</b>
 *
 *  This class can represent a Green's function in a time range \f$ t_0 \geq t_i \geq t_0 - t_c \f$.
 *  of `tc` timesteps backward from some leading timestep `t0`. The value `t0` is not stored in the class.
 * 
 *  The class `herm_matrix_moving` stores two non-redundant Keldysh components:
 *  - Retarded component \f$ C^\mathrm{R}(t_0-i,t_0-i-s) \f$ for \f$ i,s=0,\ldots,t_c\f$
 *  - Lesser component \f$ C^<(t0-i,t_0-i-s) \f$ for \f$ i,s=0,\ldots,t_c \f$.
 *
 *  The \f$ C(t,t') \f$ can be scalar type or matrix-valued.  Only square-matrix-type contour functions are fully supported.
 *
 */
  template <typename T> class herm_matrix_moving{
  public:
    typedef std::complex<T> cplx;
    typedef T scalar_type;
    /* construction, destruction */
    herm_matrix_moving();
    ~herm_matrix_moving();
    herm_matrix_moving(int tc,int size1=1,int sig=-1);
    herm_matrix_moving(const herm_matrix_moving &g);
    herm_matrix_moving & operator=(const herm_matrix_moving &g);
    void clear(void);
    void resize(int tc,int size1);
    int element_size(void) const{ return element_size_;}
    int size1(void) const{ return size1_;}
    int size2(void) const{ return size2_;}
    int tc(void) const{ return tc_;}
    int sig(void) const{ return sig_;}
    inline cplx * lesptr(int delti,int deltj){return les_[delti] + deltj*element_size_;}
    inline cplx * retptr(int delti,int deltj){return ret_[delti] + deltj*element_size_;}
    template<class Matrix> void get_les(int delti,int deltj,Matrix &M) ;
    template<class Matrix> void get_gtr(int delti,int deltj,Matrix &M);
    template<class Matrix> void get_ret(int delti,int deltj,Matrix &M);
    inline void get_les(int delti,int deltj,cplx &x);
    inline void get_gtr(int delti,int deltj,cplx &x);
    inline void get_ret(int delti,int deltj,cplx &x);
    cplx density_matrix(int delti);
    template<class Matrix> void density_matrix(int i,Matrix &M);
    template<class Matrix> void set_les(int delti,int deltj,Matrix &M);
    template<class Matrix> void set_ret(int delti,int deltj,Matrix &M);
    inline void set_les(int delti,int deltj,cplx x);
    inline void set_ret(int delti,int deltj,cplx x);
    void clear_timestep(int delti);
    void set_timestep_zero(int delti);
    void smul(int delti, T alpha);
    void smul(int delti, cplx alpha);
    void set_timestep(int delti,herm_matrix_timestep_moving_view<T> &g);
    void set_timestep(int delti,herm_matrix_timestep_moving<T> &g);
    void set_timestep(int delti,herm_matrix_moving<T> &g,int delti_g);
    void set_matrixelement(int delti, int i1, int i2, herm_matrix_timestep_moving_view<T> &g, int j1, int j2);
    void set_matrixelement(int delti, int i1, int i2, herm_matrix_timestep_moving<T> &g, int j1, int j2);
    void set_matrixelement(int delti, int i1, int i2, herm_matrix_moving<T> &g, int delti_g, int j1, int j2);
    void incr_timestep(int delti,herm_matrix_timestep_moving_view<T> &g,cplx alpha=1.0);
    void incr_timestep(int delti,herm_matrix_timestep_moving<T> &g,cplx alpha=1.0);
    void incr_timestep(int delti,herm_matrix_moving<T> &g,int delti_g,cplx alpha=1.0);
    // --------------------- DATA exchange with HERM_MATRIX ---------------------
    void set_timestep(int delti,herm_matrix_timestep_view<T> &g,herm_matrix_timestep_view<T> &gcc);
    void set_timestep(int delti,herm_matrix_timestep<T> &g,herm_matrix_timestep<T> &gcc);
    void set_timestep(int delti,herm_matrix<T> &g,herm_matrix<T> &gcc,int tstp_g);
    void set_from_G_backward(herm_matrix<T> &g,herm_matrix<T> &gcc,int tstp_g);
    void set_from_G_backward(herm_matrix<T> &g,int tstp_g);
    /** @brief Sets the data from a `herm_matrix` backward for the given timestep index `tstp_g`. */
    void set_from_G_backward(const herm_pseudo<T> &g,int tstp_g);
    /** @brief Sets the data from a `herm_matrix` backward for the given timestep index `tstp_g`. */
    void set_from_G_backward(int tstp_g,const herm_pseudo<T> &g);
    // --------------------- text file INPUT/OUTPUT ---------------------
    void print_to_file(const char *file,int precision=16);
    void read_from_file(const char *file);        
    // --------------------- function mult (only for leading timestep) ---------------------
    void left_multiply(function_moving<T> &ft);
    void right_multiply(function_moving<T> &ft);
    void left_multiply_hermconj(function_moving<T> &ft);
    void right_multiply_hermconj(function_moving<T> &ft);
    // --------------------- forward nove ---------------------
    void forward(void);
    // --------------------- hdf5 ---------------------
#if CNTR_USE_HDF5 == 1
    /** @brief Writes data to an HDF5 file for the specified group ID. */
    void write_to_hdf5(hid_t group_id);
      /** @brief Writes data to an HDF5 file for the specified group ID and creating  a subgroup given a name. */
    void write_to_hdf5(hid_t group_id, const char *groupname);
    /** @brief Writes data to an HDF5 file with a specified filename and group name. */
    void write_to_hdf5(const char *filename, const char *groupname,
		       h5_mode mode = h5_mode::create_truncate);
    /** @brief Writes data in timeslice delti to an HDF5 file for the specified group ID in the form a timeslice*/
    void write_timestep_to_hdf5(int delti,hid_t group_id);
    void write_timestep_to_hdf5(int delti,hid_t group_id, const char *groupname);
    /** @brief Writes data to an HDF5 file with a specified filename and group name. */
    void write_timestep_to_hdf5(int delti,const char *filename, const char *groupname,
		       h5_mode mode = h5_mode::create_truncate);

    /** @brief Reads data from an HDF5 file for the specified group ID. */
    void read_from_hdf5(hid_t group_id);
    /** @brief Reads data from an HDF5 file for the specified group ID and groupname. */
    void read_from_hdf5(hid_t group_id, const char *groupname);
    /** @brief Reads data from an HDF5 file with a specified filename and group name. */
    void read_from_hdf5(const char *filename, const char *groupname);   
#endif
    // TODO: MPI routines!
    
  private:
    /** \brief <b> Pointer to the data array stoirng lesser and retarded component.</b> */
    cplx* data_;
    /** \brief <b> Pointer to a pointer storing the lesser component.</b> */
    cplx** les_;
    /** \brief <b> Pointer to a pointer storing the retarded component.</b> */
    cplx** ret_;
    /** \brief <b> Cutoff time.</b> */
    int tc_;
    /** \brief <b> Number of the colums in the Matrix form.</b> */
    int size1_;
    /** \brief <b> Number of the rows in the Matrix form (currently: always `size1=size2)</b> */
    int size2_;
    /** \brief <b> Size of the Matrix form; `size1*size2`. </b> */
    int element_size_;
    /** \brief <b> Bose = +1, Fermi =-1. </b> */
    int sig_; // Bose = +1, Fermi =-1
  };


}  // namespace cntr

#endif  // CNTR_HERM_MATRIX_MOVING_DECL_H
