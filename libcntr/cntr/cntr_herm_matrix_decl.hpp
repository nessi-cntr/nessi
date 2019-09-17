#ifndef CNTR_HERM_MATRIX_DECL_H
#define CNTR_HERM_MATRIX_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

template <typename T> class function;
template <typename T> class herm_matrix_timestep;

template <typename T>
/** \brief <b> Class `herm_matrix` for two-time contour objects \f$ C(t,t') \f$
 * with hermitian symmetry.</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 *  This class contains the data structures for representing bosonic (\f$\eta = +1 \f$)
 *  or fermionic (\f$\eta = -1 \f$) or two-time contour
 *  functions \f$ C(t,t') \f$ and provides a number of tools to perform operations
 *  on such functions.
 *  The class `herm_matrix` stores the non-redundant Keldysh components
 *   - Matsubara component \f$ C^\mathrm{M}(\tau_k) \f$ for i=0,...,`nt`,
 *   - retarded component \f$ C^\mathrm{R}(t_i,t_j) \f$ for i=0,...,`nt`, j=0,...,i ,
 *   - lesser component \f$ C^<(t_j,t_i) \f$ for i=0,...,`nt`, j=0,...,i ,
 *   - left-mixing component \f$ C^\rceil(t_i,\tau_k) \f$ for i=0,...,`nt`, i=0,...,`nt`.
 *
 *  All other Keldysh components can be expressed by the hermitian symmetries
 > \f{align*}{
 *   C^<(t_i,t_j) &= - [C^<(t_j,t_i)]^\ddagger \ , \\
 *   C^\mathrm{R}(t_i,t_j) &= [C^\mathrm{A}(t_j,t_i)]^\ddagger \ , \\
 *   C^\lceil(\tau_k,t_i) &= -\eta [C^\rceil(t_i,\beta - \tau_k)]^\ddagger \ .
 * \f}
 *
 * The contour function \f$ C(t,t') \f$ can be of scalar type or matrix-valued.
 * Square-matrix-type contour functions are fully supported, while non-square-
 * matrix contour functions are under development.
 *
 *  If `nt = 0`, only the Matsubara component is stored.
 *
 */
class herm_matrix {
  public:
    typedef std::complex<T> cplx;
    typedef T scalar_type;

    /* construction, destruction */
    herm_matrix();
    ~herm_matrix();
    herm_matrix(int nt, int ntau, int size1 = 1, int sig = -1);
    herm_matrix(int nt, int ntau, int size1, int size2, int sig);
    herm_matrix(const herm_matrix &g);
    herm_matrix &operator=(const herm_matrix &g);
#if __cplusplus >= 201103L
    herm_matrix(herm_matrix &&g) noexcept;
    herm_matrix &operator=(herm_matrix &&g) noexcept;
#endif
    /* resize */
    void resize_discard(int nt, int ntau, int size1);
    void resize_nt(int nt);
    void resize(int nt, int ntau, int size1);
    void clear(void);
    /* access size etc ... */
    /// @private
    int element_size(void) const { return element_size_; }
    int size1(void) const { return size1_; }
    int size2(void) const { return size2_; }
    int ntau(void) const { return ntau_; }
    int nt(void) const { return nt_; }
    int sig(void) const { return sig_; }
    void set_sig(int sig) { sig_ = sig; }
    /* conversion from other types */
    // herm_matrix<T> & herm_matrix(const matrix<T> &g1);
    // herm_matrix<T> & herm_matrix(const scalar<T> &g1);
    // raw pointer to elements ... to be used with care
    inline cplx *lesptr(int i, int j);
    inline cplx *retptr(int i, int j);
    inline cplx *tvptr(int i, int j);
    inline cplx *matptr(int i);
    inline const cplx *lesptr(int i, int j) const;
    inline const cplx *retptr(int i, int j) const;
    inline const cplx *tvptr(int i, int j) const;
    inline const cplx *matptr(int i) const;
    // reading basic and derived elements to any Matrix type
    template <class Matrix>
    void get_les(int i, int j, Matrix &M) const;
    template <class Matrix>
    void get_gtr(int i, int j, Matrix &M) const;
    template <class Matrix>
    void get_ret(int i, int j, Matrix &M) const;
    template <class Matrix>
    void get_vt(int i, int j, Matrix &M) const;
    template <class Matrix>
    void get_tv(int i, int j, Matrix &M) const;
    template <class Matrix>
    void get_mat(int i, Matrix &M) const;
    template <class Matrix>
    void get_matminus(int i, Matrix &M) const;
    // reading complex numbers:
    // these will adress only (0,0) element for dim>1:
    inline void get_les(int i, int j, cplx &x) const;
    inline void get_gtr(int i, int j, cplx &x) const;
    inline void get_ret(int i, int j, cplx &x) const;
    inline void get_vt(int i, int j, cplx &x) const;
    inline void get_tv(int i, int j, cplx &x) const;
    inline void get_mat(int i, cplx &x) const;
    inline void get_matminus(int i, cplx &x) const;
    cplx density_matrix(int i) const;
    // should work with eigen
    template <class Matrix>
    void density_matrix(int tstp, Matrix &M) const;

    // writing basic elements:
    template <class Matrix>
    void set_les(int i, int j, Matrix &M);
    template <class Matrix>
    void set_ret(int i, int j, Matrix &M);
    template <class Matrix>
    void set_tv(int i, int j, Matrix &M);
    template <class Matrix>
    void set_mat(int i, Matrix &M);
    template<class Matrix> void set_mat_real(int i,Matrix &M);
    void set_mat_herm(int i);
    void set_mat_herm(void);
    inline void set_les(int i, int j, cplx x);
    inline void set_ret(int i, int j, cplx x);
    inline void set_tv(int i, int j, cplx x);
    inline void set_mat(int i, cplx x);
    // INPUT/OUTPUT
    void print_to_file(const char *file, int precision = 16) const;
    void read_from_file(const char *file);
    // read up to timestep nt1: NO RESIZE
    void read_from_file(int nt1, const char *file);
// HDF5 I/O
#if CNTR_USE_HDF5 == 1
    void write_to_hdf5(hid_t group_id);
    void write_to_hdf5(hid_t group_id, const char *groupname);
    void write_to_hdf5(const char *filename, const char *groupname);
    void read_from_hdf5(hid_t group_id);
    void read_from_hdf5(hid_t group_id, const char *groupname);
    void read_from_hdf5(const char *filename, const char *groupname);
    // read up to timestep nt1: NO RESIZE
    void read_from_hdf5(int nt1, hid_t group_id);
    void read_from_hdf5(int nt1, hid_t group_id, const char *groupname);
    void read_from_hdf5(int nt1, const char *filename, const char *groupname);
    // writing some other compressed format: just a collection of timeslices
    // -1,0,dt,2*dt,...
    void write_to_hdf5_slices(hid_t group_id, int dt);
    void write_to_hdf5_slices(hid_t group_id, const char *groupname, int dt);
    void write_to_hdf5_slices(const char *filename, const char *groupname,
                              int dt);
    void write_to_hdf5_tavtrel(hid_t group_id, int dt);
    void write_to_hdf5_tavtrel(hid_t group_id, const char *groupname, int dt);
    void write_to_hdf5_tavtrel(const char *filename, const char *groupname,
                               int dt);
#endif
    // simple operations on the timesteps (n<0: Matsubara)
    void set_timestep_zero(int tstp);
    void set_timestep(int tstp, herm_matrix &g1);
    void set_timestep(int tstp, herm_matrix_timestep<T> &timestep);
    void get_timestep(int tstp, herm_matrix_timestep<T> &timestep) const;
    void get_timestep(int tstp, herm_matrix<T> &timestep) const;
    void incr_timestep(int tstp, herm_matrix_timestep<T> &timestep,
                       cplx alpha = 1.0);
    void incr_timestep(int tstp, herm_matrix<T> &g, cplx alpha = 1.0);
    // check what this is actually doing
    void incr_timestep(herm_matrix<T> &g, cplx alpha = 1.0);

    void set_matrixelement(int tstp, int i1, int i2,
                           herm_matrix_timestep<T> &timestep);
    void set_matrixelement(int tstp, int i1, int i2, herm_matrix<T> &g);
    // this(i1,i2) <- g(j1,j2)
    void set_matrixelement(int tstp, int i1, int i2,
                           herm_matrix_timestep<T> &g, int j1, int j2);
    void set_matrixelement(int tstp, int i1, int i2, herm_matrix<T> &g,
                           int j1, int j2);
    void set_matrixelement(int i1,int i2,herm_matrix<T> &g);
    void set_matrixelement(int i1,int i2,herm_matrix<T> &g,int j1,int j2);

    void set_submatrix(std::vector<int> &i1,std::vector<int> &i2,
      herm_matrix<T> &g,std::vector<int> &j1,std::vector<int> &j2);
    void set_submatrix(int tstp, std::vector<int> &i1,std::vector<int> &i2,
      herm_matrix<T> &g,std::vector<int> &j1,std::vector<int> &j2);

    // multiplication with function

    // check if multiplication with pointer can be removed (there might be
    //  internal dependencies)
    void
    left_multiply(int tstp, cplx *f0, cplx *ft,
                  T weight = 1.0); // version with function object is safer
    void
    right_multiply(int tstp, cplx *f0, cplx *ft,
                   T weight = 1.0); // version with function object is safer
    void left_multiply(int tstp, function<T> &ft, T weight = 1.0);
    // - check why underscore (private routine?)
    // - check what std::integral_constant is doing
    // - add interface to multiply with a matrix from left/right
    template <int SIZE>
    void left_multiply_(int tstp, function<T> &ft, T weight,
                        std::integral_constant<int, SIZE>);
    void right_multiply(int tstp, function<T> &ft, T weight = 1.0);
    template <int SIZE>
    void right_multiply_(int tstp, function<T> &ft, T weight,
                         std::integral_constant<int, SIZE>);
    void left_multiply_hermconj(
        int tstp, function<T> &ft,
        T weight = 1.0); // G(t,t') -> weight*f(t)^dag * G(t,t')
    void right_multiply_hermconj(
        int tstp, function<T> &ft,
        T weight = 1.0); // G(t,t') -> weight*G(t,t') * f(t')^dag
    void smul(int tstp, T weight);
    void smul(int tstp, cplx weight);
// MPI UTILS
#if CNTR_USE_MPI == 1
    void Reduce_timestep(int tstp, int root);
    void Bcast_timestep(int tstp, int root);
    void Send_timestep(int tstp, int dest, int tag);
    void Recv_timestep(int tstp, int root, int tag);
#endif
  private:
    int les_offset(int t, int t1) const;
    int ret_offset(int t, int t1) const;
    int tv_offset(int t, int tau) const;
    int mat_offset(int tau) const;

  private:
    /// @private
    /** \brief <b> Pointer to the lesser component. For \f$ t \leq t_1 \f$,'les_+ \f$((t_1 * (t_1 + 1)) / 2 + t) \; *\f$ element\_size'  corresponds to \f$(0,0)\f$-component of \f$ G^<(t,t_1) \f$ </b> */
    cplx *les_;
    /// @private
    /** \brief <b> Pointer to the retarded component. For \f$ t \geq t_1 \f$,'ret_+ \f$((t * (t + 1)) / 2 + t_1) \;*\f$ element\_size'  corresponds to \f$(0,0)\f$-component of \f$ G^R(t,t_1) \f$ </b> */
    cplx *ret_;
    /// @private
    /** \brief <b> Pointer to the left-mixing component. 'tv_+ \f$(t * (n_{\rm tau} + 1) + \tau) \;*\f$ element\_size'  corresponds to \f$(0,0)\f$-component of \f$ G^\rceil(t,\tau) \f$ </b> */
    cplx *tv_; 
    /// @private
    /** \brief <b> Pointer to the Matsubra component. 'mat_+ \f$ \tau \;*\f$ element\_size'  corresponds to \f$(0,0)\f$-component of \f$ G^M(\tau) \f$ </b> */
    cplx *mat_;;
    /// @private
    /** \brief <b> Maximum number of the time steps.</b> */
    int nt_;
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
    int sig_; // Bose = +1, Fermi =-1
};

}  // namespace cntr

#endif  // CNTR_HERM_MATRIX_DECL_H
