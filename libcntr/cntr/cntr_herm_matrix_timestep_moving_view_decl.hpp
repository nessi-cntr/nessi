#ifndef CNTR_HERM_TIMESTEP_MOVING_VIEW_DECL_H
#define CNTR_HERM_TIMESTEP_MOVING_VIEW_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {
  
  template <typename T> class herm_matrix_timestep;
  template <typename T> class herm_matrix;
  template <typename T> class herm_matrix_timestep_moving;
  template <typename T> class herm_matrix_moving;
  
  /** 
   * @brief Class for interfacing with herm_matrix_timestep_moving without copying data. 
   * 
   * This class provides an interface to `herm_matrix_timestep_moving` or `herm_matrix_moving`
   * storing only pointers to the data  instead of owning it. This avoids unnecessary memory 
   * duplication but requires careful handling to prevent memory errors due to raw pointer usage.
   * 
   * The herm_matrix_timestep_moving represent a Greens function at physical timestep t_lead.
   * Arguments of the moving timestep are understood relative to t_lead:
   * this->retptr(delt) points to G^ret(t_lead,t_lead-delt),  
   * this->lesptr(delt) points to G^les(t_lead,t_lead-delt) for delt=0,...,tc_.
   * The pysical timestep t_lead is not stored
   *
   * Each element of the Greens function is a size1 x size2 matrix. Currently only  size1 = size2 is supported
   * @tparam T Data type (typically a numerical type like double or complex<double>).
   */
  template <typename T>
  class herm_matrix_timestep_moving_view {
  public:
    typedef std::complex<T> cplx; ///< Complex number type.
    typedef T scalar_type; ///< Scalar type.
    
    herm_matrix_timestep_moving_view();
    ~herm_matrix_timestep_moving_view();
    herm_matrix_timestep_moving_view(const herm_matrix_timestep_moving_view &g);
    herm_matrix_timestep_moving_view(herm_matrix_timestep_moving<T> &g);
    herm_matrix_timestep_moving_view(herm_matrix_moving<T> &g,int delt_g);
    herm_matrix_timestep_moving_view(herm_matrix_moving<T> &g);
    herm_matrix_timestep_moving_view &operator=(const herm_matrix_timestep_moving_view &g);
    
    /** @brief Returns the number of columns. */
    int size1(void) const { return size1_; }    
    /** @brief Returns the number of rows (CurrentlyL always the same as number of columns!) */
    int size2(void) const { return size2_; }    
    /** @brief Returns the statistical sign (Bose = +1, Fermi = -1). */
    int sig(void) const { return sig_; }    
    /** @brief Returns the time cut-off. */
    int tc(void) const { return tc_; }
    int element_size(void) const{ return element_size_;}
    
    
    /** @brief Pointer to retarded component at a given time, representing G^R(t_lead,t_lead-deltj) */
    inline std::complex<T> *retptr(int deltj) { return ret_ + deltj * element_size_; }    
    /** @brief Pointer to lesser component at a given time, representing G^<(t_lead,t_lead-deltj) */
    inline std::complex<T> *lesptr(int deltj) { return les_ + deltj * element_size_; }
    
    template<class Matrix> void set_les(int deltj,Matrix &M);
    template<class Matrix> void set_ret(int deltj,Matrix &M);
    template<class Matrix> void get_les(int deltj,Matrix &M);
    template<class Matrix> void get_gtr(int deltj,Matrix &M);
    template<class Matrix> void get_ret(int deltj,Matrix &M);
    inline void set_les(int deltj,cplx x);
    inline void set_ret(int deltj,cplx x);
    inline void get_les(int deltj,cplx &x);
    inline void get_gtr(int deltj,cplx &x);
    inline void get_ret(int deltj,cplx &x);
    std::complex<T> density_matrix(void);
    template <class Matrix> void density_matrix(Matrix &M);

    void set_timestep(herm_matrix_timestep_moving_view<T> &g1);
    void set_timestep(herm_matrix_timestep_moving<T> &g1);
    void set_timestep(herm_matrix_moving<T> &g1,int delt_g);
    void set_timestep(herm_matrix_timestep_view<T> &G,herm_matrix_timestep_view<T> &Gcc);
    void set_timestep(herm_matrix_timestep<T> &G,herm_matrix_timestep<T> &Gcc);
    void set_timestep(herm_matrix<T> &G,herm_matrix<T> &Gcc,int tstp_g);
    
    void clear_timestep(void);
      void set_timestep_zero(void);
    void smul(T alpha);
    void smul(cplx alpha);
    
    void set_matrixelement(int i1, int i2, herm_matrix_timestep_moving_view<T> &G, int j1, int j2);
    void set_matrixelement(int i1, int i2, herm_matrix_timestep_moving<T> &G, int j1, int j2);
    void set_matrixelement(int i1, int i2, herm_matrix_moving<T> &G, int delt_g, int j1, int j2);
    void left_multiply(function_moving<T> &ft);
    void left_multiply_hermconj(function_moving<T> &ft);
    void right_multiply(function_moving<T> &ft);
    void right_multiply_hermconj(function_moving<T> &ft);
    void incr_timestep(herm_matrix_timestep_moving_view<T> &g, cplx alpha);
    void incr_timestep(herm_matrix_timestep_moving<T> &g, cplx alpha);
    void incr_timestep(herm_matrix_moving<T> &g,int delt_g, cplx alpha);
    
    
    #if CNTR_USE_MPI == 1
    /** @brief Performs MPI reduction. */
    void MPI_Reduce(int root);
    #endif

#if CNTR_USE_HDF5 == 1
    /** @brief Writes data to an HDF5 file for the specified group ID. */
    void write_to_hdf5(hid_t group_id);
      /** @brief Writes data to an HDF5 file for the specified group ID into a subgroup with specified name. */
    void write_to_hdf5(hid_t group_id, const char *groupname);
    /** @brief Writes data to an HDF5 file with a specified filename and group name. */
    void write_to_hdf5(const char *filename, const char *groupname,
		       h5_mode mode = h5_mode::create_truncate);
#endif

    
    /// @private
    std::complex<T> *ret_; ///< Pointer to retarded component.
    /// @private
    std::complex<T> *les_; ///< Pointer to lesser component.
    /// @private
    int tc_; ///< Time cut-off.
    /// @private
    int size1_; ///< Number of columns.
    /// @private
    int size2_; ///< Number of rows.
    /// @private
    int element_size_; ///< Total size of the matrix.
    /// @private
    int sig_; ///< Bose/Fermi statistical sign.
  };

}  // namespace cntr

#endif  // CNTR_HERM_TIMESTEP_MOVING_VIEW_DECL_H
