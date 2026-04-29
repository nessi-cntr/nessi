#ifndef CNTR_HERM_MATRIX_TIMESTEP_MOVING_DECL_H
#define CNTR_HERM_MATRIX_TIMESTEP_MOVING_DECL_H

#include "cntr_herm_matrix_timestep_decl.hpp"
//#include "cntr_exception.hpp"
#include "cntr_elements.hpp"
#include "cntr_function_decl.hpp"
#include "cntr_herm_matrix_decl.hpp"
#include "cntr_herm_matrix_timestep_view_impl.hpp"


namespace cntr {

   /** 
   * @brief Stored tha data at one memory-truncated timeslice of a greens function
   * 
   * The herm_matrix_timestep_moving represent a Greens function at physical timestep t0.
   * Arguments of the moving timestep are understood relative to t0:
   * this->retptr(trelj) points to G^ret(t0,t0-trelj),  this->lesptr(trelj) points to G^les(t0,t0-trelj) for trelj=0,...,tc_
   * the pysical timestep t0 is not stored
   * 
   * Each element of the Greens function is a size1 x size2 matrix. Currently only  size1 = size2 is supported
   * @tparam T Data type (typically a numerical type like double or complex<double>).
   */
  template <typename T> class herm_matrix_timestep_moving{
  public:
    typedef std::complex<T> cplx;
    typedef T scalar_type;
    /* construction, destruction */
    herm_matrix_timestep_moving();
    ~herm_matrix_timestep_moving();
    herm_matrix_timestep_moving(int tc,int size1=1,int sig=-1); 
    herm_matrix_timestep_moving(const herm_matrix_timestep_moving &g);
    herm_matrix_timestep_moving(herm_matrix_moving<T> &g,int delti_g);
    herm_matrix_timestep_moving & operator=(const herm_matrix_timestep_moving &g);
    void resize(int tc,int size1);

    /** @brief Returns the number of columns. */
    int size1(void) const { return size1_; }    
    /** @brief Returns the number of rows (CurrentlyL always the same as number of columns!) */
    int size2(void) const { return size2_; }    
    /** @brief Returns the statistical sign (Bose = +1, Fermi = -1). */
    int sig(void) const { return sig_; }    
    /** @brief Returns the time cut-off. */
    int tc(void) const { return tc_; }
    /** @brief Returns the time elemet size = size1*size2. */
    int element_size(void) const{ return element_size_;}
    
    
    /** @brief Pointer to retarded component at a given time, representing G^R(t0,t0-deltj) */
    inline std::complex<T> *retptr(int deltj) { return ret_ + deltj * element_size_; }    
    /** @brief Pointer to lesser component at a given time, representing G^<(t0,t0-deltj) */
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

    void set_timestep(herm_matrix_moving<T> &g1,int delti_g);

    void set_timestep(herm_matrix_timestep_view<T> &G,herm_matrix_timestep_view<T> &Gcc);  

    void set_timestep(herm_matrix_timestep<T> &G,herm_matrix_timestep<T> &Gcc); 

    void set_timestep(herm_matrix<T> &G,herm_matrix<T> &Gcc,int tstp_g);     
    
    /** @brief Clears the stored timestep data (set to zero) */
    void clear_timestep(void);
    /** @brief Clears the stored timestep data. (set to zero) */
    void set_timestep_zero(void);
    /** @brief Scales the stored values by a scalar. */
    void smul(T alpha);
    void smul(cplx alpha);
    

    void set_matrixelement(int i1, int i2, herm_matrix_timestep_moving_view<T> &G, int j1, int j2);

    void set_matrixelement(int i1, int i2, herm_matrix_timestep_moving<T> &G, int j1, int j2);

    void set_matrixelement(int i1, int i2, herm_matrix_moving<T> &G, int delti_g, int j1, int j2);    

    void left_multiply(function_moving<T> &ft);

    void left_multiply_hermconj(function_moving<T> &ft);

    void right_multiply(function_moving<T> &ft);

    void right_multiply_hermconj(function_moving<T> &ft);


    void incr_timestep(herm_matrix_timestep_moving_view<T> &g, cplx alpha);

    void incr_timestep(herm_matrix_timestep_moving<T> &g, cplx alpha);

    void incr_timestep(herm_matrix_moving<T> &g,int delti_g, cplx alpha);
    
    
    #if CNTR_USE_MPI == 1
    /** @brief Performs MPI reduction. */
    void MPI_Reduce(int root);
    //void Bcast_timestep(int tstp, int ntau, int size1, int root);
    //void Send_timestep(int tstp, int ntau, int size1, int dest, int tag);
    //void Recv_timestep(int tstp, int ntau, int size1, int root, int tag);
    #endif

    void print_to_file(const char *file,int precision=16);
    void read_from_file(const char *file);        

#if CNTR_USE_HDF5 == 1
    /** @brief Writes data to an HDF5 file for the specified group ID. */
    void write_to_hdf5(hid_t group_id);
    void write_to_hdf5(hid_t group_id, const char *groupname);
    /** @brief Writes data to an HDF5 file with a specified filename and group name. */
    void write_to_hdf5(const char *filename, const char *groupname,
		       h5_mode mode = h5_mode::create_truncate);
    /** @brief Reads data from an HDF5 file for the specified group ID. */
    void read_from_hdf5(hid_t group_id);
    void read_from_hdf5(hid_t group_id, const char *groupname);
    /** @brief Reads data from an HDF5 file with a specified filename and group name. */
    void read_from_hdf5(const char *filename, const char *groupname);   
#endif
    
    
  private:
    /** \brief <b> Pointer to the data for the time step (t0). </b> */
    cplx* data_;
    /** \brief <b> Pointer to the lesser part for the time step (t0). 'data_+\f$(tc+1+tc)\f$ * element_size' correponds to (0,0)-component of \f$ G^<(t0,t0-tc)\f$ and 'les_(tc)'. </b> */
    cplx* les_;
    /** \brief <b> Pointer to the retarded part for the time step (t0). 'data_+\f$tc\f$ * element_size' correponds to (0,0)-component of \f$ G^<(t0,t0-tc)\f$ and 'ret_(tc)'. </b> */
    cplx* ret_;
    /** \brief <b> Given cutoff time tc </b> */
    int tc_;
    /** \brief <b> Number of the colums in the Matrix form.</b> */
    int size1_;
    /** \brief <b> Number of the rows in the Matrix form.</b> */
    int size2_;
    /** \brief <b> Size of the Matrix form; size1*size2. </b> */
    int element_size_;
    /** \brief <b> Bose = +1, Fermi =-1. </b> */
    int sig_;
  };

} // namespace cntr

#endif  // CNTR_HERM_MATRIX_TIMESTEP_DECL_H
