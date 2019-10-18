#ifndef CNTR_HERM_TIMESTEP_VIEW_DECL_H
#define CNTR_HERM_TIMESTEP_VIEW_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

template <typename T> class herm_matrix_timestep;
template <typename T> class herm_matrix;
template <typename T> class function;

template <typename T>
/** \brief <b> Class `herm_matrix_timestep_view` serves for interfacing with class herm_matrix_timestep</b>
 * without copying the data.
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *  Occasionally you may have a pre-defined object of herm_matrix or herm_matrix_timestep that you want to use.
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

class herm_matrix_timestep_view {
  public:
    typedef std::complex<T> cplx;
    typedef T scalar_type;
    herm_matrix_timestep_view();
    ~herm_matrix_timestep_view();
    herm_matrix_timestep_view(const herm_matrix_timestep_view &g);
    herm_matrix_timestep_view &operator=(const herm_matrix_timestep_view &g);
    herm_matrix_timestep_view(herm_matrix_timestep<T> &g);
    herm_matrix_timestep_view(int tstp, herm_matrix<T> &g);
    herm_matrix_timestep_view(int tstp, int ntau, int size1, int size2, int sig);
    herm_matrix_timestep_view(int tstp, herm_matrix_timestep<T> &g);
    herm_matrix_timestep_view(int tstp, herm_matrix_timestep_view<T> &g);
    ///////////////////////////////////////////
    int size1(void) const { return size1_; }
    int size2(void) const { return size2_; }
    int ntau(void) const { return ntau_; }
    int tstp(void) const { return tstp_; }
    int sig(void) const { return sig_; }
    /// @private
    inline std::complex<T> *retptr(int i) { return ret_ + i * element_size_; }
    /// @private
    inline std::complex<T> *tvptr(int i) { return tv_ + i * element_size_; }
    /// @private
    inline std::complex<T> *lesptr(int i) {
        return les_ + i * element_size_;
    };
    /// @private
    inline std::complex<T> *matptr(int i) { return mat_ + i * element_size_; }
    ///////////////////////////////////////////
    // the following routines set the ponters:
    void set_to_data(std::complex<T> *data, int tstp, int ntau, int size,
                     int sig);
    void set_to_data(int tstp, herm_matrix<T> &g);
    void set_to_data(herm_matrix_timestep<T> &g);
// HDF5 I/O
#if CNTR_USE_HDF5 == 1
    void write_to_hdf5(hid_t group_id);
    void write_to_hdf5(hid_t group_id, const char *groupname);
    void write_to_hdf5(const char *filename, const char *groupname);
    // READ: Read only data, no resize!
    void read_from_hdf5(hid_t group_id, const char *groupname);
    void read_from_hdf5(hid_t group_id);
    void read_from_hdf5(const char *filename, const char *groupname);
#endif

    /*//////////////////////////////////////////////////////////////////////////////////
     OPERATE ON THE DATA POINTED TO BY THE VIEW OBJECT:

     the argument G can be anything on which one can create a temporary
     herm_matrix_timestep_view object by
     herm_matrix_timestep_view<T>(tstp_,G)

    //////////////////////////////////////////////////////////////////////////////////*/
    // ** GET_DATA: copy data into the array pointed to by the view:
    void get_data(herm_matrix_timestep_view<T> &g);
    template <class GG>
    void get_data(GG &g);

    /// @private
    template <class Matrix> void set_ret(int i, int j, Matrix &M);
    /// @private
    template <class Matrix> void set_les(int i, int j, Matrix &M);
    /// @private
    template <class Matrix> void set_tv(int i, int j, Matrix &M);
    /// @private
    template <class Matrix> void set_mat(int i, Matrix &M);

    /// @private
    template <class Matrix> void get_les(int i, int j, Matrix &M);
    /// @private
    template <class Matrix> void get_ret(int i, int j, Matrix &M);
    /// @private
    template <class Matrix> void get_tv(int i, int j, Matrix &M);
    /// @private
    template <class Matrix> void get_mat(int i, Matrix &M);
    /// @private
    template <class Matrix> void get_matminus(int i, Matrix &M);


    /// @private
    void left_multiply(int tstp, std::complex<T> *f0, std::complex<T> *ft, T weight = 1.0);
    void left_multiply(int tstp, function<T> &ft, T weight = 1.0);
    /// @private
    void right_multiply(int tstp, std::complex<T> *f0, std::complex<T> *ft, T weight = 1.0);
    void right_multiply(int tstp, function<T> &ft, T weight = 1.0);
    void left_multiply_hermconj(int tstp, function<T> &ft, T weight = 1.0);
    void right_multiply_hermconj(int tstp, function<T> &ft, T weight = 1.0);

    void set_timestep_zero(int tstp);

    void set_timestep(int tstp, herm_matrix<T> &g1);
    void set_timestep(int tstp, herm_matrix_timestep<T> &g1);

    /// @private
    std::complex<T> density_matrix(int tstp);

    void density_matrix(int tstp, std::complex<T> &M);
    template <class Matrix> void density_matrix(int tstp, Matrix &M);


    void get_timestep(int tstp,herm_matrix_timestep<T> &timestep);
    // ** SET MATRIX_ELEMENT:  this(i1,i2) <= g(j1,j2)
    void set_matrixelement(int i1, int i2, herm_matrix_timestep_view<T> &g,
                           int j1, int j2);
    template <class GG>
    void set_matrixelement(int i1, int i2, GG &g, int j1, int j2);
    // ** INCR_TIMESTEP  this += alpha * g
    void incr_timestep(herm_matrix_timestep_view<T> &g,
                       std::complex<T> alpha);
    void incr_timestep(herm_matrix_timestep_view<T> &g, T alpha);
    template <class GG>
    void incr_timestep(GG &g, std::complex<T> alpha);
    template <class GG>
    void incr_timestep(GG &g, T alpha);
    // ** SMUL  this *= alpha
    void smul(T alpha);
#if CNTR_USE_MPI == 1
    /// @private
    void MPI_Reduce(int root);
    // void Reduce_timestep(int tstp, int root); 
    void Bcast_timestep(int tstp, int root);
    void Send_timestep(int tstp, int dest, int tag);
    void Recv_timestep(int tstp, int root, int tag);
#endif
    
    /////////////////////////////////////////////////////////////////////////
    // private:
    void get_data(std::complex<T> *ret, std::complex<T> *les,
                  std::complex<T> *tv, std::complex<T> *mat);
    /// @private
    /** \brief <b> Pointer to the retarded component. 'ret_+ \f$t \;*\f$ element\_size'  corresponds to \f$(0,0)\f$-component of \f$ G^R(tstp,t) \f$ </b> */
    std::complex<T> *ret_;
    /// @private
    /** \brief <b> Pointer to the lesser component. 'less_+ \f$t \;*\f$ element\_size'  corresponds to \f$(0,0)\f$-component of \f$ G^<(t,tstp) \f$ </b> */
    std::complex<T> *les_;
    /// @private
    /** \brief <b> Pointer to the left-mixing component. 'tv_+ \f$\tau \;*\f$ element\_size'  corresponds to \f$(0,0)\f$-component of \f$ G^\rceil(tstp,\tau) \f$ </b> */
    std::complex<T> *tv_;
    /// @private
    /** \brief <b> Pointer to the Matsubara component. 'mat_+ \f$\tau \;*\f$ element\_size'  corresponds to \f$(0,0)\f$-component of \f$ G^M(\tau) \f$ </b> */
    std::complex<T> *mat_;
    /// @private
    /** \brief <b> time step.</b> */
    int tstp_;
    /// @private
    /** \brief <b> Number of the time grids on the Matsubara axis.</b> */
    int ntau_;
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

#endif  // CNTR_HERM_TIMESTEP_VIEW_DECL_H
