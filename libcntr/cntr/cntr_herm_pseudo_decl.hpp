#ifndef CNTR_HERM_PSEUDE_DECL_H
#define CNTR_HERM_PSEUDE_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

template <typename T> class function;
template <typename T> class herm_matrix;
template <typename T> class herm_matrix_timestep;

/*#########################################################################################
#
#   ...   a pseudoparticle Greenfunction (herm_pseudo) is identical to a usual
#         herm_matrix, but one has ret=gtr instead of ret=gtr-les, such that the
#         conversion from a full_matrix and its timesteps is different
#
#########################################################################################*/
/// @private
template <typename T> class herm_pseudo {
  public:
    typedef std::complex<T> cplx;
    /* construction, destruction */
    herm_pseudo();
    ~herm_pseudo();
    herm_pseudo(int nt, int ntau, int size1 = 1, int sig = -1);
    herm_pseudo(const herm_pseudo &g);
    herm_pseudo &operator=(const herm_pseudo &g);
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
    /* conversion from other types: pseudoparticle GF,
    matrix GF with different storage structure etc...*/

    // herm_pseudo<T> & herm_pseudo(const matrix<T> &g1);
    // herm_pseudo<T> & herm_pseudo(const scalar<T> &g1);

    // raw pointer to elements ... to be used with care
    inline cplx *lesptr(int i, int j);
    inline cplx *retptr(int i, int j);
    inline cplx *tvptr(int i, int j);
    inline cplx *matptr(int i);
    // reading basic and derived elements to any Matrix type
    template <class Matrix> void get_les(int i, int j, Matrix &M);
    template <class Matrix> void get_gtr(int i, int j, Matrix &M);
    template <class Matrix> void get_ret(int i, int j, Matrix &M);
    template <class Matrix> void get_vt(int i, int j, Matrix &M);
    template <class Matrix> void get_tv(int i, int j, Matrix &M);
    template <class Matrix> void get_mat(int i, Matrix &M);
    template <class Matrix> void get_matminus(int i, Matrix &M);
    // reading complex numbers:
    // these will adress only (0,0) element for dim>1:
    inline void get_les(int i, int j, cplx &x);
    inline void get_gtr(int i, int j, cplx &x);
    inline void get_ret(int i, int j, cplx &x);
    inline void get_vt(int i, int j, cplx &x);
    inline void get_tv(int i, int j, cplx &x);
    inline void get_mat(int i, cplx &x);
    inline void get_matminus(int i, cplx &x);
    // writing basic elements:
    template <class Matrix> void set_les(int i, int j, Matrix &M);
    template <class Matrix> void set_ret(int i, int j, Matrix &M);
    template <class Matrix> void set_tv(int i, int j, Matrix &M);
    template <class Matrix> void set_mat(int i, Matrix &M);
    inline void set_les(int i, int j, cplx x);
    inline void set_ret(int i, int j, cplx x);
    inline void set_tv(int i, int j, cplx x);
    inline void set_mat(int i, cplx x);
    // density matrix
    cplx density_matrix(int i);
    template <class Matrix> void density_matrix(int tstp, Matrix &M);
    // INPUT/OUTPUT
    void print_to_file(const char *file, int precision = 16);
    void read_from_file(const char *file);
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
#endif

    // simple operations on the timesteps (n<0: Matsubara)
    void set_timestep_zero(int tstp);
    void set_timestep(int tstp, herm_pseudo &g1);
    void set_timestep(int tstp, herm_matrix_timestep<T> &timestep);
     void get_timestep(int tstp, herm_matrix_timestep<T> &timestep);
    void incr_timestep(int tstp, herm_matrix_timestep<T> &timestep, cplx alpha = 1.0);
    void left_multiply(int tstp, cplx *f0, cplx *ft,
                       T weight = 1.0); // version with function object is safer
    void right_multiply(int tstp, cplx *f0, cplx *ft,
                        T weight = 1.0); // version with function object is safer
    void left_multiply(int tstp, function<T> &ft, T weight = 1.0);
    void right_multiply(int tstp, function<T> &ft, T weight = 1.0);
    //
    void smul(int tstp, T alpha);
    void smul(int tstp, cplx alpha);

  private:
    /// @private
    cplx *les_;
    /// @private
    cplx *ret_;
    /// @private
    cplx *tv_;
    /// @private
    cplx *mat_;
    /// @private
    int nt_;
    /// @private
    int ntau_;
    /// @private
    int size1_;
    /// @private
    int size2_;
    /// @private
    int element_size_;
    /// @private
    int sig_; // Bose = +1, Fermi =-1
};

} // namespace cntr

#endif  // CNTR_HERM_PSEUDE_DECL_H
