#ifndef NESS_FFT_H
#define NESS_FFT_H
#include "ness_GF_decl.hpp"
#define SET_DW 1
#define SET_DT 1
namespace ness
{

  /** \brief <b> Class `fft_solver` contains the FFTW3 data structures for Fourier transformation of Green's functions. </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *  \par Purpose
   * <!-- ========= -->
   *
   * > This class contains FFT as member functions to transform Green's function objects.
   *
   * \details
   * Given:
   * - `nt_`: Number of time points (of the `GF` object).
   * - `nfreq_`: Number of frequency points. The frequency `GF` needs \f$ {\tt nfreq} = 2/3 ({\tt nt} - 1) \f$ frequency points for the FFT to work.
   * - `Nft_`: Total number of FFT points, \f$ {\tt Nft\_} = 2 ({\tt nt\_}-1)\f$.
   * - `in`: Internal FFT input array, GF data are read in. Has length `Nft_`.
   * - `out`: Internal FFT output array,  data are read out of this array.. Has length `Nft_`.
   * - `plan_for`: FFTW plan for forward transform (here to time)
   * - `plan_back`: FFTW plan for backward transform (here to frequency)
   */
  class fft_solver
  {
  public:
    fft_solver(int nt, size_t FFTW_FLAG);
    
    ~fft_solver();
    void to_time(GF &outG, const GF &inG, int set_dt = 0);
    void to_freq(GF &outG, const GF &inG, int set_dw = 0);
    
    
    fftw_complex* in; //!< Internal FFT input array.
    fftw_complex* out; //!< Internal FFT output array.
    int Nft_; //!< Total number of FFT points.
    int nfreq_; //!<  Number of frequency points (of the frequency `GF` object).
    int nt_; //!< Number of time points (of the `GF` object).
    fftw_plan plan_for;  //!< FFTW forward plan.
    fftw_plan plan_back; //!< FFTW backward plan.
    
  };
}
#endif
