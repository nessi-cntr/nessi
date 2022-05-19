#ifndef CNTR_FUNCTION_DECL_H
#define CNTR_FUNCTION_DECL_H

#include "cntr_global_settings.hpp"

namespace cntr {

  template <typename T>
  /** \brief <b> Class `function` for objects \f$ f(t) \f$ with time on real axis.</b>
   *
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   *  This class contains the data of a function of time on real axis. The function can be
   *  scalar/matrix valued. The class also contains various operations on functions.
   *  t = -1 means the function on the imaginary(Matsubara) axis.
   *
   *
   *
   */


  class function {
  public:
    // this is like a vector with pointer access
    typedef std::complex<T> cplx;
    function();
    ~function();
    function(int nt, int size1 = 1);
    function(int nt,int size1,int size2);
    function(const function &f);
    function &operator=(const function &f);
#if __cplusplus >= 201103L
    function(function &&f) noexcept;
    function &operator=(function &&f) noexcept;
#endif
    /// @private
    int element_size(void) const { return element_size_; }
    int size1(void) const { return size1_; }
    int size2(void) const { return size2_; }
    int nt(void) const { return nt_; }
    inline cplx *ptr(int t) { return data_ + (t + 1) * element_size_; }
    inline const cplx *ptr(int t) const {
      return data_ + (t + 1) * element_size_;
    }
    void resize(int nt, int size1);
    void set_zero(void);
    void set_constant(cplx *f0); // f0 must be size*size
    template<class EigenMatrix>
    void set_constant(EigenMatrix &M);
    template <class EigenMatrix>
    void set_value(int tstp, EigenMatrix &M);
    template <class EigenMatrix>
    void set_value(int tstp, cplx x);
    template <class EigenMatrix>
    void get_value(int tstp, EigenMatrix &M) const;
    void smul(T weight);
    template<class EigenMatrix>
    void set_matrixelement(int tstp,int i1,int i2,EigenMatrix &M,int j1,int j2);
    void set_matrixelement(int i1,int i2,function &g,int j1,int j2);
    void left_multiply(cplx *f, T weight = 1.0);
    void right_multiply(cplx *f, T weight = 1.0);
    void left_multiply(function<T> &f, T weight = 1.0);
    void right_multiply(function<T> &f, T weight = 1.0);
    void get_matrixelement(int i1, int i2, function<T> &g);
    void incr(function<T> &g, T weight = 1.0);

    cplx &operator[](int t) { return *ptr(t); } // useful only for size=1
    const cplx &operator[](int t) const {
      return *ptr(t);
    } // useful only for size=1
    void print_to_file(const char *file, int precision = 16) const;
    void read_from_file(const char *file);
    // READ UP TO NT1, NO RESIZE
    void read_from_file(int nt1, const char *file);
#if CNTR_USE_HDF5 == 1
    // HDF5 I/O
    void write_to_hdf5(hid_t group_id) const;
    void write_to_hdf5(hid_t group_id, const char *groupname) const;
    void write_to_hdf5(const char *filename, const char *groupname) const;
    void read_from_hdf5(hid_t group_id);
    void read_from_hdf5(hid_t group_id, const char *groupname);
    void read_from_hdf5(const char *filename, const char *groupname);
    void read_from_hdf5(int nt1, hid_t group_id);
    void read_from_hdf5(int nt1, hid_t group_id, const char *groupname);
    void read_from_hdf5(int nt1, const char *filename, const char *groupname);
#endif
#if CNTR_USE_MPI==1
    void Bcast_timestep(int tstp,int root);
#endif
    /** \brief <b> Pointer to the function in the Matrix form on the real-time axis (\f$f(t)\f$) ; 'data_+\f$ (t+1)\;*\f$element_size' corresponds to (0,0)-component of \f$f(t)\f$. </b> */
    cplx *data_;
    /** \brief <b> Maximum number of the time steps.</b> */
    int nt_;
    /** \brief <b> Number of the colums in the Matrix form.</b> */
    int size1_;
    /** \brief <b> Number of the rows in the Matrix form.</b> */
    int size2_;
    /// @private
    /** \brief <b> Size of the Matrix form; size1*size2. </b> */
    int element_size_;
    /** \brief <b> Size of the data stored for the function on the real-time axis including \f$ t=-1\f$; \f$ (n_t + 2)\f$ * size1 * size2 . </b> */
    int total_size_;
  };

  template <typename T>
  /** \brief <b> Class `function_moving` for objects \f$ f(t) \f$ with time on real axis.</b>
   *
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   *  This class contains the data of a function of time on real axis within a given cutoff 
   *  time. The function can be scalar/matrix valued. The class also contains various 
   *  operations on funtions.
   *
   */
  class function_moving {
  public:
    typedef std::complex<T> cplx;
    /* construction, destruction */
    function_moving();
    ~function_moving();
    function_moving(int tc,int size1=1);
    function_moving(const function_moving &g);
    function_moving & operator=(const function_moving &g);
    void clear(void);
    void resize(int tc,int size1);
    /* access size etc ... */
    /// @private
    int element_size(void) const{ return element_size_;}
    int size1(void) const{ return size1_;}
    int size2(void) const{ return size2_;}
    int tc(void) const{ return tc_;}
    void set_t0(int t0);
    // raw pointer to elements ... to be used with care
    inline cplx * ptr(int i){return value_[i];}  // points to Gret(t0-i,t0-i-j)
    template<class EigenMatrix>
    void set_value(int i,EigenMatrix &M);
    template<class EigenMatrix>
    void get_value(int i,EigenMatrix &M) const;
    cplx & operator[](int i){return *ptr(i);} // useful only for size=1
    cplx & operator[](int i) const{return *ptr(i);} // useful only for size=1
    // INPUT/OUTPUT
    void print_to_file(const char *file,int precision=16);
    void read_from_file(const char *file);
    // DATA exchange with HERM_MATRIX
    void forward(void);
    /* set from function only complex double so far... */
    //internal translation to complex double can cause memory issues think about passing correct Matrixtype
    void set_from_function_backward(int tstp, function<T>& f);

  private:
    cplx* data_;
    cplx** value_;
    /** \brief <b> Cutoff time</b> */
    int tc_;
    /** \brief <b> Current physical timestep.</b> */
    int t0_;
    /** \brief <b> Number of the colums in the Matrix form.</b> */
    int size1_;
    /** \brief <b> Number of the rows in the Matrix form.</b> */
    int size2_;
    /// @private
    /** \brief <b> Size of the Matrix form; size1*size2. </b> */
    int element_size_;
  };

}  // namespace cntr

#endif  // CNTR_FUNCTION_DECL_H
