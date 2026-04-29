#ifndef FFT_ARRAY_DECL_HPP
#define FFT_ARRAY_DECL_HPP

#include <fftw3.h>
#include "ness2_global_settings.hpp"
#include <vector>
#include <stdexcept>

namespace ness2 {

  /**
   * \brief Grid information for time and frequency axes associated with FFT arrays.
   *
   * This class describes the time and frequency grids with steps \f$ h \f$ and \f$ d\omega \f$ and length \f$ N_{ft} \f$ used in `fft_array` and `herm_matrix_ness`.
   *
   * \details
   * Given:
   * - `h_`: time step,
   * - `Nft_`: total FFT grid points,
   * - \f$ {\tt dw\_} = \frac{2\pi}{{\tt Nft\_} \cdot {\tt h\_}} \f$: frequency step.
   *
   */
  class grid_info {
  public:
    int Nft_;   //!< total number of FFT grid points
    double h_;  //!< time step
    double dw_; //!< frequency step
    
    /**
     * \brief Constructor to initialize grid with time step \f$ h \f$ and length \f$ N_{ft} \f$.
     *
     * \param Nft number of time points (must be positive)
     * \param h  time step (must be positive)
     */
    grid_info(int Nft, double h);
    
    /**
     * \brief Update time step \f$ h \f$ and adjust frequency step \f$ d\omega \f$ accordingly.
     *
     * \param new_h new time step (must be positive)
     */
    void set_h(double new_h);
    
    /**
     * \brief Update frequency step \f$ d\omega \f$ and adjust time step \f$ h \f$ accordingly.
     *
     * \param new_dw new frequency step (must be positive)
     */
    void set_dw(double new_dw);
    
    /**
     * \brief Get time value \f$ t \f$ at grid index \f$ n \f$:
     * \f$ t = n \cdot h \quad \text{for} \quad n=0 \dots \frac{N_{ft}}{2}-1 \\
     * t = (n - N_{ft}) \cdot h \quad \text{for} \quad n = \frac{N_{ft}}{2} \dots N_{ft} - 1 \f$
     * \par further indices are periodically wrapped
     *
     * \param n time grid index
     */
    double time_at(int n) const;
    
      /**
       * \brief Get frequency value \f$ \omega \f$ at grid index \f$ k \f$:
       * \f$ w = k \cdot dw \quad \text{for} \quad k=0\dots \frac{N_{ft}}{2}-1 \\
       * w = (k - N_{ft}) \cdot dw \quad \text{for} \quad k=\frac{N_{ft}}{2} \dots N_{ft} - 1 \f$
       * \par further indices are periodically wrapped
       *
       * \param k freqeuncy  grid index
       */
    double freq_at(int k) const;
    
    /**
     * \brief Generate the full wrapped time grid.
     * \return vector of time values (length \f$ N_{ft} \f$)
     */
    std::vector<double> time_grid() const;
    
    /**
     * \brief Generate the full wrapped frequency grid.
     * \return vector of frequency values (length \f$ N_{ft} \f$)
     */
    std::vector<double> freq_grid() const;
  };
  

  /**
   * @brief  Class to manage FFTW threading setup and cleanup.
   * Use this to initialize FFTW with threads once at program start.
   */
  class FFT_OMP_Manager {
  public:
    /// Initialize FFTW threads (pthreads backend). \f$ n\f$ is used for plans created afterwards.
    static void initialize(void);
    
    /// Change the thread count used by *new* plans (does not affect existing plans) to `n`.
    static void set_threads_for_new_plans(int n);
    
    /// Clean up threading resources (call only after all plans are destroyed).
    static void finalize();
    
    FFT_OMP_Manager() = delete;
    FFT_OMP_Manager(const FFT_OMP_Manager&) = delete;
    FFT_OMP_Manager& operator=(const FFT_OMP_Manager&) = delete;
    
  private:
    static bool initialized_;
    static int  current_threads_;
  };
  

  /**
   * @brief A class that manages multiple strided FFTs using FFTW.
   *
   * The data is arranged such that each FFT has stride `element_size_ = size1_ * size1_`,
   * and is of length `Nft`. Forward and backward plans are created once and reused.
   */
  class fft_array {
  public:
    /**
     * @brief Constructor
     * @param Nft FFT length
     * @param size Logical size (`stride = size x size`)
     * @param FFTW_FLAG FFTW planning flag (e.g., `FFTW_ESTIMATE`)
     */
    fft_array(int Nft, int size, unsigned FFTW_FLAG = FFTW_ESTIMATE);
    // fft_array();
    
    /// Destructor
    ~fft_array();
    
    /// Copy constructor
    fft_array(const fft_array& other);
    
    /// Copy assignment operator
    fft_array& operator=(const fft_array& other);

    /**
     * @brief Resize the `fft_array` and rebuild FFTW plans.
     *
     * Frees any existing memory and FFTW plans, then allocates
     * new arrays and creates new FFTW plans for the given sizes.
     * After calling, the array is reset and ready for use.
     *
     * @param new_Nft New number of time samples per transform.
     * @param new_size1 New matrix dimension (`new_size1 × new_size1`).
     * @param FFTW_FLAG FFTW planning flag (default: `FFTW_ESTIMATE`).
     */
    void resize(int new_Nft, int new_size1, unsigned FFTW_FLAG = FFTW_ESTIMATE);

      /**
       * Performs a forward FFT (time domain → frequency domain).
       * This is because `plan_back_` was created with `FFTW_FORWARD`.
       *
       * \f$ X[k] = \sum_{n=0}^{N_{ft}-1} x[n] \cdot \exp(-2\pi i \cdot k \cdot n / N_{ft}) \f$
       *
       * Note: This is the standard DFT definition (no scaling).
       */
    void fft_to_freq();

      /**
       * Performs an inverse FFT (frequency domain → time domain).
       * This is because `plan_for_` was created with `FFTW_BACKWARD`.
       *
       * \f$ x[n] = \sum_{k=0}^{N_{ft}-1} X[k] \cdot \exp(+2\pi i \cdot k \cdot n / N_{ft}) \f$
       *
       * Note: FFTW does NOT normalize this, so the result will be scaled by \f$ N_{ft} \f$ .
       * You may need to manually divide by \f$ N_{ft} \f$ if you want the exact inverse.
      */
    void fft_to_time();
       
    
    /**
     * @brief Increment this `fft_array` \f$ A \f$  by another scaled array \f$ B \f$:  \f$ A  \to A +  \alpha \cdot B \f$
     * @param B Another fft_array with matching dimensions.
     * @param alpha Complex scaling factor.
     * @param type Whether to operate on time or frequency arrays (`fft_domain::time` or `fft_domain::freq`).
     */
    void incr(const fft_array& B, std::complex<double> alpha, fft_domain type);
      
    /**
     * @brief Scale this `fft_array` \f$ A \f$  in-place by a complex scalar \f$ \alpha \f$:  \f$ A  \to  \alpha \cdot A \f$
     * @param alpha Complex scaling factor.
     * @param type Indicates whether to scale time or frequency array (`fft_domain::time` or `fft_domain::freq`).
     */
    void smul(std::complex<double> alpha, fft_domain type);
    
    /**
     * @brief Set one matrix element's time or frequency series of this `fft_array` \f$ A \f$  from another `fft_array` \f$ B \f$.
     *
     * This sets \f$ A(i_1, i_2)(t) = B(j_1, j_2)(t) \f$ across all times or frequencies.
     *
     * @param i1 Row index in \f$ A(i_1, i_2) \f$.
     * @param i2 Column index in \f$ A(i_1, i_2) \f$.
     * @param B Source fft_array.
     * @param j1 Row index in \f$ B(j_1, j_2) \f$.
     * @param j2 Column index in \f$ B(j_1, j_2) \f$.
     * @param type Whether to operate on time or frequency arrays  (`fft_domain::time` or `fft_domain::freq`).
     */
    void set_matrixelement(int i1, int i2, const fft_array& B, int j1, int j2, fft_domain type);
    
    /**
     * @brief Set all elements of either the time or frequency array to zero.
     * @param type Indicates whether to clear the time or frequency array.
     */
    void set_zero(fft_domain type);
        
      /**
       * @brief Multiply each time/frequency matrix of this `fft_array` \f$ A \f$ from the left by a complex matrix \f$ M \f$: \f$ A(t) \to M \cdot A(t) \f$.
       *
       * @tparam MatrixType Any Eigen-compatible complex matrix type.
       * @param M The matrix to left-multiply with (must be `size1_ × size1_`)
       * @param type Whether to apply to time or frequency domain data (`fft_domain::time` or `fft_domain::freq`).
       */
    template <typename MatrixType>
    void left_multiply(const MatrixType& M, fft_domain type);

      /**
       * @brief Multiply each time/frequency matrix of this `fft_array` \f$ A \f$ from the right by a complex matrix \f$ M \f$: \f$ A(t) \to A(t) \cdot  M \f$.
       *
       * @tparam MatrixType Any Eigen-compatible complex matrix type.
       * @param M The matrix to right-multiply with (must be `size1_ xsize1_`)
       * @param type Whether to apply to time or frequency domain data (`fft_domain::time` or `fft_domain::freq`).
       */
    template <typename MatrixType>
    void right_multiply(const MatrixType& M, fft_domain type);

      /**
       * @brief Multiply each time/frequency matrix of this `fft_array` \f$ A \f$ from the left by the hermitian conjugate of a complex matrix \f$ M \f$: \f$ A(t) \to M^\dagger \cdot A(t) \f$.
       *
       * @tparam MatrixType Any Eigen-compatible complex matrix type.
       * @param M The matrix to left-multiply with (must be `size1_ × size1_`)
       * @param type Whether to apply to time or frequency domain data (`fft_domain::time` or `fft_domain::freq`).
       */
    template <typename MatrixType>
    void left_multiply_hermconj(const MatrixType& M, fft_domain type);
    
      /**
       * @brief Multiply each time/frequency matrix of this `fft_array` \f$ A \f$ from the right by the hermitian conjugate of a complex matrix \f$ M \f$: \f$ A(t) \to A(t)  \cdot  M^\dagger  \f$.
       *
       * @tparam MatrixType Any Eigen-compatible complex matrix type.
       * @param M The matrix to left-multiply with (must be `size1_ × size1_`)
       * @param type Whether to apply to time or frequency domain data (`fft_domain::time` or `fft_domain::freq`).
       */
    template <typename MatrixType>
    void right_multiply_hermconj(const MatrixType& M, fft_domain type);
    
      /**
       * @brief Set this `fft_array`\f$ A \f$ to \f$ A(t) \f$ = \f$ M \f$ at time or frequency index \f$ t \f$ (with periodic wrapping).
       *
       * If \f$ t < 0 \f$ or \f$ t \geq N_{ft} \f$, it is automatically wrapped periodically:
       * \f$ t \equiv t \bmod N_{ft} \f$ in \f$ [0, N_{ft} - 1] \f$.
       *
       * @tparam MatrixType Any Eigen-compatible complex matrix type.
       * @param t Time index (wrapped periodically).
       * @param M Matrix to assign.
       * @param type Whether to apply to time or frequency array (`fft_domain::time` or `fft_domain::freq`).
       */
    template <typename MatrixType>
    void set_element(int t, const MatrixType& M, fft_domain type);
    
      /**
       * @brief Get \f$ M = A(t) \f$ from this `fft_array`\f$ A \f$ at or frequency time index \f$ t \f$ (with periodic wrapping).
       *
       * If \f$ t < 0 \f$ or \f$ t \geq N_{ft} \f$, it is automatically wrapped periodically:
       * \f$ t \equiv t \bmod N_{ft} \f$ in \f$ [0, N_{ft} - 1] \f$.
       *
       * @tparam MatrixType Any Eigen-compatible complex matrix type.
       * @param t Time index (wrapped periodically).
       * @param M Matrix to fill with the values at time index \f$ t \f$.
       * @param type Whether to read from time or frequency array (`fft_domain::time` or `fft_domain::freq`).
       */
    template <typename MatrixType>
    void get_element(int t, MatrixType& M, fft_domain type) const;

      /**
       * @brief Set this `fft_array`\f$ A \f$ to \f$ A(t) = \tt{val} \f$ at time or frequency index \f$ t \f$ (with periodic wrapping).
       *
       * Only valid when \f$ \text{size1}\_ = 1 \f$ (scalar case).
       * If \f$ t < 0 \f$ or \f$ t \geq N_{ft} \f$, it is automatically wrapped periodically:
       * \f$ t \equiv t \bmod N_{ft} \f$ in \f$ [0, N_{ft} - 1] \f$.
       *
       * @param t Time index (wrapped periodically).
       * @param val Complex scalar to assign.
       * @param type Whether to apply to time or frequency array (`fft_domain::time` or `fft_domain::freq`).
       */
    void set_element(int t, cplx val, fft_domain type);
    
    /**
     * @brief Get \f$ {\tt val} = A(t) \f$ from this `fft_array`\f$ A \f$ at time or frequency index \f$ t \f$ (with periodic wrapping).
     *
     * Only valid when \f$ \text{size1}\_ = 1 \f$ (scalar case).
     * If \f$ t < 0 \f$ or \f$ t \geq N_{ft} \f$, it is automatically wrapped periodically:
     * \f$ t \equiv t \bmod N_{ft} \f$ in \f$ [0, N_{ft} - 1] \f$.
     *
     * @param t Time index (wrapped periodically).
     * @param val Complex scalar to retrieve.
     * @param type Whether to read from time or frequency array (`fft_domain::time` or `fft_domain::freq`).
     */
    void get_element(int t, cplx& val, fft_domain type) const;


    /**
     * @brief Print the time and frequency arrays to a plain text file.
     *
     * Outputs:
     * - Header: # `Nft` `size1_` `precision`
     * - Blocks: time and frequency (each with one line per \f$ t \f$).
     *
     * @param filename Path to the text file.
     * @param precision Number of digits for real/imaginary part (default: 12).
     */
    void print_to_file(const std::string& filename, int precision = 12) const;
    
    /**
     * @brief Read the time and frequency arrays from a plain text file written with `print_to_file`
     *
     * Resizes internal arrays if needed using the provided `FFTW_FLAG`.
     *
     * @param filename Path to the text file.
     * @param FFTW_FLAG FFTW planning flag to use if resizing (default: `FFTW_ESTIMATE`).
     */
    void read_from_file(const std::string& filename, unsigned FFTW_FLAG = FFTW_ESTIMATE);

    
#if CNTR_USE_HDF5 == 1

    /**
     * @brief Write the `fft_array` data (time and frequency arrays) to an HDF5 group.
     *
     * Stores attributes for `Nft_`, `size1_`, and `element_size_`,
     * and datasets for the `time_` and `freq_` arrays as 3D arrays of shape `[Nft_, size1_, size1_]`.
     *
     * @param group_id HDF5 group identifier.
     */
    void write_to_hdf5(hid_t group_id)  const;
    
    /**
     * @brief Read the `fft_array` data (time and frequency arrays) from an HDF5 group.
     *
     * If the stored array dimensions differ from the current object,
     * automatically resizes the object using the provided `FFTW_FLAG`.
     *
     * @param group_id HDF5 group identifier.
     * @param FFTW_FLAG FFTW planning flag to use if resizing (default: `FFTW_ESTIMATE`).
     */
    void read_from_hdf5(hid_t group_id, unsigned FFTW_FLAG = FFTW_ESTIMATE);
    
#endif

    fftw_complex* time_;  ///< Pointer to time-domain data
    fftw_complex* freq_;  ///< Pointer to frequency-domain data
    int Nft_;             ///< Number of time samples per transform
    int size1_;            ///< Logical size (`stride = size1_ * size1_`)
    int element_size_;    ///< Stride for each transform (`size1_ * size1_)

  private:
    fftw_plan plan_for_;   ///< FFTW forward plan
    fftw_plan plan_back_;  ///< FFTW backward plan
  };
  


  /// definitions of templated member functions:
  
  template <typename MatrixType>
  void fft_array::left_multiply(const MatrixType& M, fft_domain type) {
    if (M.rows() != size1_ || M.cols() != size1_) {
      throw std::invalid_argument("left_multiply: M must be size1_ x size1_");
    }
    
    if (size1_ == 1) {
      smul(M(0, 0), type);
      return;
    }
    
    int N = size1_;
    fftw_complex* data = (type == fft_domain::time) ? time_ : freq_;
        
    for (int t = 0; t < Nft_; ++t) {
      fftw_complex* frame_ptr = &data[t * element_size_];
      Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        A_frame(reinterpret_cast<std::complex<double>*>(frame_ptr), N, N);      
      A_frame = M * A_frame;
    }
  }
  
  template <typename MatrixType>
  void fft_array::right_multiply(const MatrixType& M, fft_domain type) {
    if (M.rows() != size1_ || M.cols() != size1_) {
      throw std::invalid_argument("right_multiply: M must be size1_ x size1_");
    } 
    if (size1_ == 1) {
      smul(M(0, 0), type);
      return;
    }
    int N = size1_;
    fftw_complex* data = (type == fft_domain::time) ? time_ : freq_;
    
    
    for (int t = 0; t < Nft_; ++t) {
      fftw_complex* frame_ptr = &data[t * element_size_];
      Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        A_frame(reinterpret_cast<std::complex<double>*>(frame_ptr), N, N);      
      A_frame = A_frame * M;
    }
  }
  
  template <typename MatrixType>
  void fft_array::left_multiply_hermconj(const MatrixType& M, fft_domain type) {
    if (M.rows() != size1_ || M.cols() != size1_) {
      throw std::invalid_argument("left_multiply_hermconj: M must be size1_ x size1_");
    }
    
    if (size1_ == 1) {
        smul(std::conj(M(0, 0)), type);
        return;
    }
    
    int N = size1_;
    fftw_complex* data = (type == fft_domain::time) ? time_ : freq_;
    
    for (int t = 0; t < Nft_; ++t) {
      fftw_complex* frame_ptr = &data[t * element_size_];
      Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        A_frame(reinterpret_cast<std::complex<double>*>(frame_ptr), N, N);      
      A_frame = M.adjoint() * A_frame;
    }
  }
  
  template <typename MatrixType>
  void fft_array::right_multiply_hermconj(const MatrixType& M, fft_domain type) {
    if (M.rows() != size1_ || M.cols() != size1_) {
      throw std::invalid_argument("right_multiply_hermconj: M must be size1_ x size1_");
    }    
    if (size1_ == 1) {
      smul(std::conj(M(0, 0)), type);
      return;
    }    
    int N = size1_;
    fftw_complex* data = (type == fft_domain::time) ? time_ : freq_;      
    for (int t = 0; t < Nft_; ++t) {
      fftw_complex* frame_ptr = &data[t * element_size_];
      Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        A_frame(reinterpret_cast<std::complex<double>*>(frame_ptr), N, N);      
      A_frame = A_frame * M.adjoint();
    }
  }


 
  template <typename MatrixType>
  void fft_array::set_element(int t, const MatrixType& M, fft_domain type) {
    // Wrap t periodically
    t = ((t % Nft_) + Nft_) % Nft_;    
    fftw_complex* data = (type == fft_domain::time) ? time_ : freq_;
    fftw_complex* frame_ptr = &data[t * element_size_];
    
    if (size1_ == 1) {
      frame_ptr[0][0] = M(0, 0).real();
      frame_ptr[0][1] = M(0, 0).imag();
    } else {
      if (M.rows() != size1_ || M.cols() != size1_) {
	throw std::invalid_argument("set_element: matrix size mismatch");
      }
      
      Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
	A_frame(reinterpret_cast<std::complex<double>*>(frame_ptr), size1_, size1_);
      A_frame = M;
    }
  }
  
  template <typename MatrixType>
  void fft_array::get_element(int t, MatrixType& M, fft_domain type) const {
    // Wrap t periodically
    t = ((t % Nft_) + Nft_) % Nft_;    
    const fftw_complex* data = (type == fft_domain::time) ? time_ : freq_;
    const fftw_complex* frame_ptr = &data[t * element_size_];
    
    if (size1_ == 1) {
      M.resize(1, 1);
      M(0, 0) = std::complex<double>(frame_ptr[0][0], frame_ptr[0][1]);
    } else {
      M.resize(size1_, size1_);
      
      Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
            A_frame(reinterpret_cast<const std::complex<double>*>(frame_ptr), size1_, size1_);
      M = A_frame;
    }
  }
  
  
  
} // namespace ness2


#endif // FFT_ARRAY_HPP
