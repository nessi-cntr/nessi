#ifndef NESS2_HERM_MATRIX_NESS_DECL_HPP
#define NESS2_HERM_MATRIX_NESS_DECL_HPP

#include "cntr/ness2_global_settings.hpp"
#include "cntr/ness2_fft_decl.hpp"
#include "cntr/fourier.hpp"
#include "cntr/cntr.hpp"

namespace ness2 {
  
/**
 * @brief Class representing a two-component Green's function: retarded (ret) and lesser (les).
 * Each component is an `fft_array` of size `(Nft) × (size1_ × size1_)`.
 */
class herm_matrix_ness {
public:
  /**
   * @brief Constructor.
   * 
   * @param Nft: Length of arrays
   * @param size Logical matrix size (`size × size`).
   * @param FFTW_FLAG Planning flag for FFTW (default: `FFTW_ESTIMATE`).
   */
  herm_matrix_ness(int Nft=4, int size=1,unsigned FFTW_FLAG = FFTW_ESTIMATE);

  //herm_matrix_ness();
  
  /// Destructor.
  ~herm_matrix_ness();
  
  /// Copy constructor.
  herm_matrix_ness(const herm_matrix_ness& other);
  
  /// Copy assignment operator.
  herm_matrix_ness& operator=(const herm_matrix_ness& other);
  
  /**
   * @brief Transform both components to the frequency domain.
   */
  void fft_to_freq();
  
  /**
   * @brief Transform both components back to the time domain.
   */
  void fft_to_time();
  
  
  /**
   * @brief Increment this `herm_matrix_ness` \f$ A \f$ by another scaled `herm_matrix_ness` \f$ B \f$: \f$ A^{R,<} \to A^{R,<} + \alpha \cdot B^{R,<} \f$.
   *
   * @param B Another herm_matrix_ness with matching dimensions.
   * @param alpha Complex scaling factor.
   * @param type Whether to operate on time or frequency arrays (`fft_domain::time` or `fft_domain::freq`).
   */
  void incr(const herm_matrix_ness& B, cplx alpha, fft_domain type);
  
  /**
   * @brief Scale both ret and les of this `herm_matrix_ness` \f$ A \f$ in-place by a complex scalar: \f$ A^{R,<} \to  \alpha \cdot A^{R,<} \f$.
   *
   * @param alpha Complex scaling factor.
   * @param type Whether to operate on time or frequency arrays (`fft_domain::time` or `fft_domain::freq`).
   */
  void smul(cplx alpha, fft_domain type);
  
  /**
   * @brief Set one matrix element (both ret and les) of this `herm_matrix_ness` \f$ A \f$  from another `herm_matrix_ness` \f$ B \f$.
   *
   * This sets \f$ A^{R,<}(i_1, i_2)(t) = B^{R,<}(j_1, j_2)(t) \f$ across all times or frequencies.
   *
   * @param i1 Row index in \f$ A(i_1, i_2) \f$.
   * @param i2 Column index in \f$ A(i_1, i_2) \f$.
   * @param B Source fft_array.
   * @param j1 Row index in \f$ B(j_1, j_2) \f$.
   * @param j2 Column index in \f$ B(j_1, j_2) \f$.
   * @param type Whether to operate on time or frequency arrays (`fft_domain::time` or `fft_domain::freq`).
   */
  void set_matrixelement(int i1, int i2, const herm_matrix_ness& B, int j1, int j2, fft_domain type);
  
  /**
   * @brief Set all elements of both ret and les to zero.
   * 
   * @param type Whether to clear time or frequency arrays (`fft_domain::time` or `fft_domain::freq`).
   */
  void set_zero(fft_domain type);
  
  /**
   * @brief Multiply each time/frequency matrix for both ret and les of this `herm_matrix_ness` \f$ A \f$ from the left by a complex matrix \f$ M \f$: \f$ A^{R,<}(t) \to M \cdot A^{R,<}(t) \f$.
   *
   * @tparam MatrixType Any Eigen-compatible complex matrix type.
   * @param M Matrix to left-multiply with.
   * @param type Whether to apply to time or frequency domain data (`fft_domain::time` or `fft_domain::freq`).
   */
  template <typename MatrixType>
  void left_multiply(const MatrixType& M, fft_domain type);
  
  /**
   * @brief Multiply each time/frequency matrix for both ret and les of this `herm_matrix_ness` \f$ A \f$ from the right by a complex matrix \f$ M \f$: \f$ A^{R,<}(t) \to A^{R,<}(t) \cdot M \f$.
   *
   * @tparam MatrixType Any Eigen-compatible complex matrix type.
   * @param M Matrix to right-multiply with.
   * @param type Whether to apply to time or frequency domain data (`fft_domain::time` or `fft_domain::freq`).
   */
  template <typename MatrixType>
  void right_multiply(const MatrixType& M, fft_domain type);
  
  /**
   * @brief Multiply each time/frequency matrix for both ret and les of this `herm_matrix_ness` \f$ A \f$ from the left by the hermitian conjugate of a complex matrix \f$ M \f$: \f$ A^{R,<}(t) \to M^\dagger \cdot A^{R,<}(t) \f$.
   *
   * @tparam MatrixType Any Eigen-compatible complex matrix type.
   * @param M Matrix to left-multiply with.
   * @param type Whether to apply to time or frequency domain data (`fft_domain::time` or `fft_domain::freq`).
   */
  template <typename MatrixType>
  void left_multiply_hermconj(const MatrixType& M, fft_domain type);
  
  /**
   * @brief Multiply each time/frequency matrix for both ret and les of this `herm_matrix_ness` \f$ A \f$ from the right by the hermitian conjugate of a complex matrix \f$ M \f$: \f$ A^{R,<}(t) \to A^{R,<}(t) \cdot M^\dagger \f$.
   *
   * @tparam MatrixType Any Eigen-compatible complex matrix type.
   * @param M Matrix to left-multiply with.
   * @param type Whether to apply to time or frequency domain data (`fft_domain::time` or `fft_domain::freq`).
   */
  template <typename MatrixType>
  void right_multiply_hermconj(const MatrixType& M, fft_domain type);
  
  /**
   * @brief Set the retarded component of this `herm_matrix_ness` \f$ A \f$ at time index \f$ t \f$ (with periodic wrapping) from a matrix \f$ M \f$: \f$ A^R(t) =  M \f$.
   *
   * @tparam MatrixType Any Eigen-compatible complex matrix type.
   * @param t Time/frequency index (wrapped periodically).
   * @param M Matrix to assign.
   * @param type Whether to apply to time or frequency array (`fft_domain::time` or `fft_domain::freq`).
   */
  template <typename MatrixType>
  void set_ret(int t, const MatrixType& M, fft_domain type);
  
  /**
   * @brief Get the retarded component of this `herm_matrix_ness` \f$ A \f$ at time index \f$ t \f$ (with periodic wrapping) and write to a matrix \f$ M \f$:  \f$ M =  A^R(t) \f$.
   *
   * @tparam MatrixType Any Eigen-compatible complex matrix type.
   * @param t Time/freq index (wrapped periodically).
   * @param M Matrix to fill with values.
   * @param type Whether to read from time or frequency array (`fft_domain::time` or `fft_domain::freq`).
   */
  template <typename MatrixType>
  void get_ret(int t, MatrixType& M, fft_domain type) const;
  
  /**
   * @brief Set the lesser component of this `herm_matrix_ness` \f$ A \f$ at time index \f$ t \f$ (with periodic wrapping) from a matrix \f$ M \f$:  \f$ A^<(t) =  M \f$.
   *
   * @tparam MatrixType Any Eigen-compatible complex matrix type.
   * @param t Time index (wrapped periodically).
   * @param M Matrix to assign.
   * @param type Whether to apply to time or frequency array (`fft_domain::time` or `fft_domain::freq`).
   */
  template <typename MatrixType>
  void set_les(int t, const MatrixType& M, fft_domain type);
  
  /**
   * @brief Get the lesser component of this `herm_matrix_ness` \f$ A \f$ at time index \f$ t \f$ (with periodic wrapping) and write to a matrix \f$ M \f$:  \f$ M =  A^<(t) \f$.
   *
   * @tparam MatrixType Any Eigen-compatible complex matrix type.
   * @param t Time/freq index (wrapped periodically).
   * @param M Matrix to fill with values.
   * @param type Whether to read from time or frequency array (`fft_domain::time` or `fft_domain::freq`).
   */
  template <typename MatrixType>
  void get_les(int t, MatrixType& M, fft_domain type) const;
  
  /// Access the retarded component.
  fft_array& retarded() { return ret_; }
  
  /// Access the lesser component.
  fft_array& lesser() { return les_; }
  
  /// Const access to the retarded component.
  const fft_array& retarded() const { return ret_; }
  
  /// Const access to the lesser component.
  const fft_array& lesser() const { return les_; }


  /**
   * \brief Compute the frequency-domain Fourier integral of the retarded Green's function.
   *
   * \details
   * This function computes the approximate Fourier integral
   * \f[
   * G^R(\omega) = \int_0^{T} dt\, e^{i \omega t} G^R(t),
   * \f]
   * and
   * \f[
   * G^<(\omega) = \int_{-T}^{T} dt\, e^{i \omega t} G^<(t),
   * \f]
   * where \f$ T = (N_{ft}/2- 1) \cdot h \f$, using either:
   *
   * - a simple discrete approximation if `method=FFT_TRAPEZ}` (default), and
   * - a cubically corrected Fourier transform if `method=FFT_CUBIC`, based on Numerical Recipes (Press et al., Chapter 13.9).
   *
   * The integral is computed for all frequencies
   * \f$ -(N_{ft}/2-1),\dots,(N_{ft}/2-1) \f$ , only the Nyqusit point \f$ N_{ft}/2 \f$ is set to zero.
   *
   * \param h timestep
   * \param method Method for calculating the integral transform (`FFT_TRAPEZ` default)
   */
  void integral_transform_to_freq(double h, fft_integral_method method=FFT_TRAPEZ);
  
  /**
   * \brief Compute the time-domain Fourier integral of the retarded Green's function.
   *
   * \details
   * This function computes the approximate Fourier integral
   * \f[
   * G^R(t) = -i\theta(t) \int_{-W}^{W} dw\, e^{-i \omega t} (-1/\pi)*ImG^R(w),
   * \f]
   * and 
   * \f[
   * G^<(t) = \int_{-W}^{W} dw\, e^{-i \omega t} G^<(w)/(2\pi)
   * \f]
   * where \f$ W = (N_{ft}/2 - 1) \cdot d\omega \f$, and \f$ d\omega=\pi/(N_{ft}/2 \cdot h) \f$ , and \f$ h \f$ is the timestep, using either:
   *
   * - a simple discrete approximation if `method=FFT_TRAPEZ` (default), and
   * - a cubically corrected Fourier transform if `method=FFT_CUBIC`, based on Numerical Recipes (Press et al., Chapter 13.9).
   *
   * Only times \f$ 0,...,(N_{ft}/2-1) \f$ are set, all others are set to zero
   * \param h timestep
   * \param method Method for calculating the DFT (`FFT_TRAPEZ` default)
   */
  void integral_transform_to_time(double h, fft_integral_method method=FFT_TRAPEZ);


  /**
   * @brief Enforce thermal equilibrium for the lesser Green's function.
   *
   * Computes \f$ G^<(t) \f$ in thermal equilibrium for either bosons or fermions,
   * using the fluctuation-dissipation relation and the retarded Green's function \f$ G^R(t) \f$.
   *
   * This method performs:
   * - Fourier transform of \f$ G^R(t) \f$ to \f$ G^R(\omega) \f$ using FFT with step size `h` 
   *   (equivalent to `integral_transform_to_freq` with `FFT_TRAPEZ`),
   * - Constructs \f$ G^<(\omega) \f$ via:
   *   \f[
   *   G^<(\omega) =  \mp 2i\, f(\omega - \mu)\, \mathrm{Im}\, G^R(\omega),
   *   \f]
   *   where \f$ f(x) = 1 / (e^{\beta x} \pm 1) \f$ is the Bose-Einstein (+) or Fermi-Dirac (–) distribution,
   * - Inverse FFT to obtain \f$ G^<(t) \f$.
   *
   * @param bosefermi Sign indicator: +1 for bosons, -1 for fermions
   * @param beta Inverse temperature \f$ \beta \f$
   * @param mu Chemical potential \f$ \mu \f$
   * @param h Time step size used in the FFT
   */
  void force_equilibrium(int bosefermi, double beta, double mu, double h);

  /**
   * @brief Print the Green's function (ret and les, time and frequency) to a text file.
   *
   * Outputs:
   * - Header line: # `Nft` `size1` `precision
   * - Blocks for ret_time, les_time, ret_freq, les_freq
   *
   * @param filename Path to the text file.
   * @param precision Number of digits for real/imaginary part (default: 12).
   */
  void print_to_file(const std::string& filename, int precision = 12) const;
  
  /**
   * @brief Read the Green's function data from a text file.
   *
   * Expects the same format as produced by `print_to_file`.
   * Resizes internal arrays if needed using the provided `FFTW_FLAG`.
   *
   * @param filename Path to the text file.
   * @param FFTW_FLAG FFTW planning flag to use if resizing (default: `FFTW_ESTIMATE).
   */
  void read_from_file(const std::string& filename, unsigned FFTW_FLAG = FFTW_ESTIMATE);
  
  
  
#if CNTR_USE_HDF5 == 1
  
  /**
   * @brief Helper function: Write the `herm_matrix_ness` data (ret and les components) to an hdf5 group.
   *
   * Stores attributes for `Nft_` and `size1_` and delegates writing
   * to the internal `fft_array` instances (`ret_` and `les_`).
   *
   * @param group_id hdf5 group identifier.
   */
  void write_to_hdf5(hid_t group_id) const;
  
  /**
   * @brief Helper function: Read the `herm_matrix_ness` data (ret and les components) from an hdf5 group.
   *
   * If the stored dimensions differ from the current object,
   * automatically resizes the object using the provided `FFTW_FLAG`.
   *
   * @param group_id hdf5 group identifier.
   * @param FFTW_FLAG FFTW planning flag to use if resizing (default: `FFTW_ESTIMATE`).
   */
  void read_from_hdf5(hid_t group_id, unsigned FFTW_FLAG = FFTW_ESTIMATE);


  /**
   * @brief Helper function: Write to a named subgroup under the given hdf5 group.
   *
   * Creates a subgroup named groupname, writes this object into it, and closes it.
   *
   * @param group_id hdf5 parent group identifier.
   * @param groupname Name of the subgroup to create.
   */
  void write_to_hdf5(hid_t group_id, const char* groupname) const;

  /**
   * @brief Write Green's function to a hdf5 file under the given group name.
   *
   * Opens or creates the hdf5 file if not exisitng, writes the`ret_`, `les_`,  `Nft` and `size` data of the `herm_matrix_ness` under the specified group name and closes the file.
   *
   * @param filename Name of the hdf5 file.
   * @param groupname Name of the top-level group to create.
   */
  void write_to_hdf5(const char* filename, const char* groupname) const;
  
  /**
   * @brief Helper function: Read from a named subgroup under the given hdf5 group.
   *
   * Opens the specified subgroup, reads the object, and closes the subgroup.
   *
   * @param group_id hdf5 parent group identifier.
   * @param groupname Name of the subgroup to open.
   * @param FFTW_FLAG FFTW planning flag to use if resizing (default: `FFTW_ESTIMATE`).
   */
  void read_from_hdf5(hid_t group_id, const char* groupname, unsigned FFTW_FLAG = FFTW_ESTIMATE);
  
  /**
   * @brief Read Green's function data from a given hdf5 file and group.
   *
   * Opens the desired hdf5 file, reads the `ret_`, `les_`,  `Nft` and `size` data from the specified group into the `herm_matrix_ness` object and closes the file.
   *
   * @param filename Name of the hdf5 file.
   * @param groupname Name of the top-level group.
   * @param FFTW_FLAG FFTW planning flag to use if resizing (default: `FFTW_ESTIMATE`).
   */
  void read_from_hdf5(const char* filename, const char* groupname, unsigned FFTW_FLAG = FFTW_ESTIMATE);

  
  /*
   * @brief Read data from `cntr::herm_matrix`  object `Gcntr` of matching size into this `ness2::herm_matrix_ness` \f$ A \f$.   *
   * set \f$ A^R(t)=G_{\mathrm {cntr}}^R({\tt tstp},{\tt tstp}-t) t=0...N_{ft}/2-1 \f$
   * and \f$ A^<(t)=G_{\mathrm {cntr}}^<({\tt tstp},{\tt tstp}-t) t=-(N_{ft}/2-1),...,N_{ft}/2-1 \f$
   * \f$ {\tt tstp} >= N_{ft}/2-1 \f$ required.
   *
   * @tparam T Scalar type of the input matrix (e.g., std::complex<double>).
   * @param Gcntr The source Green's function in `cntr::herm_matrix<T>` format.
   * @param tstp The reference time step (must satisfy \f$ tstp \leq Gcntr.nt() \f$).
   */
    
    
#endif
    
  template<typename T>
  void read_from_cntr(cntr::herm_matrix<T> &Gcntr,int tstp);
  
  //private:
  int Nft_;       ///< FFT length.
  int size1_;     ///< logical matrix size
  fft_array ret_; ///< Retarded component.
  fft_array les_; ///< Lesser component.
};

  // template definitions 


  template <typename MatrixType>
  void herm_matrix_ness::left_multiply(const MatrixType& M, fft_domain type) {
    ret_.left_multiply(M, type);
    les_.left_multiply(M, type);
  }
  
  template <typename MatrixType>
  void herm_matrix_ness::right_multiply(const MatrixType& M, fft_domain type) {
    ret_.right_multiply(M, type);
    les_.right_multiply(M, type);
  }
  
  template <typename MatrixType>
  void herm_matrix_ness::left_multiply_hermconj(const MatrixType& M, fft_domain type) {
    ret_.left_multiply_hermconj(M, type);
    les_.left_multiply_hermconj(M, type);
  }
  
  template <typename MatrixType>
  void herm_matrix_ness::right_multiply_hermconj(const MatrixType& M, fft_domain type) {
    ret_.right_multiply_hermconj(M, type);
    les_.right_multiply_hermconj(M, type);
  }
  
  template <typename MatrixType>
  void herm_matrix_ness::set_ret(int t, const MatrixType& M, fft_domain type) {
    ret_.set_element(t, M, type);
  }
  
  template <typename MatrixType>
  void herm_matrix_ness::get_ret(int t, MatrixType& M, fft_domain type) const {
    ret_.get_element(t, M, type);
  }
  
  template <typename MatrixType>
  void herm_matrix_ness::set_les(int t, const MatrixType& M, fft_domain type) {
    les_.set_element(t, M, type);
  }
  
  template <typename MatrixType>
  void herm_matrix_ness::get_les(int t, MatrixType& M, fft_domain type) const {
    les_.get_element(t, M, type);
  }
  
/**
 * @brief Import Green's function data from a `cntr::herm_matrix` at a specific time step.
 *
 * This method initializes the retarded and lesser components of the Green's function 
 * from a `cntr::herm_matrix<T>` at a given time step `tstp`. It extracts values 
 * \f$ G^R(t_{\text{stp}}, t_{\text{stp}} - i) \f$ and 
 * \f$ G^<(t_{\text{stp}}, t_{\text{stp}} - i) \f$ 
 * for all \f$ i \leq \min(t_{\text{stp}}, N_{\text{ft}} / 2 - 1) \f$,
 * and stores them in the time-domain representation of this `ness2::herm_matrix_ness` object.
 *
 * The lesser component is symmetrized to enforce Hermitian structure: 
 * \f$ G^<(-t) = -[G^<(t)]^\dagger \f$.
 *
 * @tparam T Scalar type of the input matrix (e.g., std::complex<double>).
 * @param Gcntr The source Green's function in `cntr::herm_matrix<T>` format.
 * @param tstp The reference time step (must satisfy \f$ {\tt tstp} \leq \tt Gcntr.nt() \f$).
 * @throws std::runtime_error If matrix dimensions mismatch or `tstp` is out of range.
 */
  template<typename T>
  void  herm_matrix_ness::read_from_cntr(cntr::herm_matrix<T> &Gcntr,int tstp){
    if(size1_!=Gcntr.size1() || tstp> Gcntr.nt()){
      throw std::runtime_error("dyson_ret read_ness_from_cntr:input param mismatch");
    }
    cdmatrix M;
    this->set_zero(fft_domain::time);
    for(int i=0;i<Nft_/2;i++){
      if(i<=tstp){
	Gcntr.get_ret(tstp,tstp-i,M);
	this->set_ret(i,M,fft_domain::time);
	Gcntr.get_les(tstp,tstp-i,M);
	this->set_les(i,M,fft_domain::time);
	if(i>0) this->set_les(-i,-M.adjoint(),fft_domain::time);
      }
    }
  }


  
} // namespace ness2

#endif // NESS2_HERM_MATRIX_NESS_DECL_HPP
