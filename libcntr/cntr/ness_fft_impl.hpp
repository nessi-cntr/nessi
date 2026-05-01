#ifndef FFT_IMPL
#define FFT_IMPL

#include <chrono>
#include "ness_fft_decl.hpp"
#define HALF_N (Nft_ / 2 + 1)
#define HALF_n (nfreq_ / 2 + 1)

#define DEBUG_GTR

namespace ness {


  /** \brief <b> Constructor for a `fft_solver` object.  Uses FFTW3.</b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *   \par Purpose
   * <!-- ========= -->
   *
   * > Initiaizes an object, which sets the necessary input and output array as well as the plan for the FFT using FFTW3 from a number of `nt` timesteps in the `GF` object to be transformed. Also the FFTW plans are set using the provided `GF` size `nt` and `FFTW_FLAG`.
   * > The Fourier transform is then a member function which takes the respective Green's functions as arguments.
   *
   * <!-- ARGUMENTS
   *      ========= -->
   *
   * @param nt
   *  [int] Number of time points (of the `GF` object).
   * @param FFTW_FLAG
   *  [size_t] Flag for FFTW3 on how to perform the FFT, usually set to `FFTW_MEASURE` or `FFTW_ESTIMATE`.
   */
  fft_solver::fft_solver(int nt, size_t FFTW_FLAG) : nt_(nt)
  {
    Nft_ = 2 * (nt_-1);
    nfreq_ = 2 * (nt_-1) / 3;
    in = new fftw_complex[Nft_];
    out = new fftw_complex[Nft_];
    plan_for = fftw_plan_dft_1d(Nft_, in, out, FFTW_FORWARD, FFTW_FLAG);
    plan_back = fftw_plan_dft_1d(Nft_, in, out, FFTW_BACKWARD, FFTW_FLAG);
  }
  /// Destructor
  fft_solver::~fft_solver()
  {
    delete[] in;
    delete[] out;
  }
  


  

/** \brief <b> Performs FFT to frequency domain of a Green's function \f$G\f$.  Uses FFTW3.</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Transforms a time `GF` Green's function object to frequency
 * > \f$ G(\omega) = \mathcal{F}[G(t)]\f$
 * > using FFTW3 and writes it to a frequency type `GF` Green's function.
 * \note Note that for `nfreq` frequency points, the time `GF` needs \f$ 3/2 {\tt nfreq}+1 \f$ time points for the FFT to work.
 * \details
 * This function computes the approximate Fourier integral
 * \f[
 * G^R(\omega) = \int_0^{T} dt\, e^{i \omega t} G^R(t),
 * \f]
 * and
 * \f[
 * G^<(\omega) = \int_{-T}^{T} dt\, e^{i \omega t} G^<(t),
 * \f]
 * where \f$ T = ({\tt Nft\_}/2- 1) \cdot h \f$ in asimple trapezoidal approximation.
 *
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param outG
 *  [GF] Frequency type output Green's function to which the result of the FFT is given.
 * @param inG
 *  [GF] Time Green's function, which is to be Fourier tranformed to frequency.
 * @param set_dw
 *  [int] Flag for renormalizing the grid of the output Green's function, 0: False, 1: True.
 */
void fft_solver::to_freq(GF &outG, const GF &inG, int set_dw)
{
	CNTR_ASSERT_EQ(NESS_ASSERT_0, inG.size1_ , outG.size1_, __PRETTY_FUNCTION__);
	CNTR_ASSERT_EQ(NESS_ASSERT_0, inG.size2_ , outG.size2_, __PRETTY_FUNCTION__);
	CNTR_ASSERT_EQ(NESS_ASSERT_0, outG.ngrid_ , this -> nfreq_, __PRETTY_FUNCTION__);
	CNTR_ASSERT_EQ(NESS_ASSERT_0, inG.ngrid_ , this -> Nft_ / 2 + 1, __PRETTY_FUNCTION__);
	CNTR_ASSERT_EQ(NESS_ASSERT_0, inG.gf_type_ , time_gf, __PRETTY_FUNCTION__);
	CNTR_ASSERT_EQ(NESS_ASSERT_0, outG.gf_type_ , freq_gf, __PRETTY_FUNCTION__);

	//int HALF_N = Nft_ / 2 + 1;
	//int HALF_n = nfreq_ / 2 + 1;
	for(int idx1 = 0; idx1 < inG.size1_; idx1++)
	{
		for(int idx2 = 0; idx2 < inG.size2_; idx2++)
		{
			//compute gret
			for(int w = 0; w < HALF_N; w++)
			{
				in[w][0] = (*inG.p_ret(w, idx1, idx2)).real();
				in[w][1] = (*inG.p_ret(w, idx1, idx2)).imag();
			}
			for(int w = HALF_N; w < Nft_; w++) {in[w][0] = 0.0; in[w][1] = 0.0;}
			fftw_execute(plan_back);
			for(int w = 0; w < HALF_n; w++)
			{
				(*outG.p_ret(w, idx1, idx2)) = out[w][0] + out[w][1] * NESS_II;
			}
			for(int w = HALF_n; w < nfreq_; w++)
			{
				(*outG.p_ret(w, idx1, idx2)) = out[Nft_ - nfreq_ + w][0] + out[Nft_ - nfreq_ + w][1] * NESS_II;
			}

			//compute glss
			for(int w = 0; w < HALF_N; w++)
			{
				in[w][0] = (*inG.p_les(w, idx1, idx2)).real();
				in[w][1] = (*inG.p_les(w, idx1, idx2)).imag();
			}

			for(int w = HALF_N; w < Nft_; w++)
			{
				in[w][0] = -(*inG.p_les(Nft_ - w, idx2, idx1)).real();
				in[w][1] = (*inG.p_les(Nft_ - w, idx2, idx1)).imag();
			}

			fftw_execute(plan_back);
			for(int w = 0; w < HALF_n; w++)
			{
				(*outG.p_les(w, idx1, idx2)) = out[w][0] + out[w][1] * NESS_II;
			}
			for(int w = HALF_n; w < nfreq_; w++)
			{
				(*outG.p_les(w, idx1, idx2)) = out[Nft_ - nfreq_ + w][0] + out[Nft_ - nfreq_ + w][1] * NESS_II;
			}


		}
	}
	outG.smul(inG.dgrid_);
	if(set_dw) outG.reset_grid(2.0 * M_PI / (inG.dgrid_ * (Nft_)));
}

/** \brief <b> Performs FFT to time domain of a Green's function \f$G\f$.  Uses FFTW3.</b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *   \par Purpose
 * <!-- ========= -->
 *
 * > Transforms a frequency Green's function object to time
 * > \f$ G(t) = \mathcal{F}[G(\omega)]\f$
 * > using FFTW3 and writes it to a time type `GF` Green's function `GF`.
 * \note Note that for `nfreq` frequency points, the time `GF` needs \f$ 3/2 {\tt nfreq}+1 \f$ time points for the FFT to work.
 * \details
 * This function computes the approximate Fourier integral
 * \f[
 * G^R(t) = -i\theta(t) \int_{-W}^{W} dw\, e^{-i \omega t} (-1/\pi)*ImG^R(w),
 * \f]
 * and
 * \f[
 * G^<(t) = \int_{-W}^{W} dw\, e^{-i \omega t} G^<(w)/(2\pi)
 * \f]
 * where \f$ W = ({\tt Nft\_}/2 - 1) \cdot d\omega \f$, and \f$ d\omega=\pi/({\tt Nft\_}/2 \cdot h) \f$ , and \f$ h \f$ is the timestep using a simple trapezoidal approximation.
 * <!-- ARGUMENTS
 *      ========= -->
 *
 * @param outG
 *  [GF] Time type output Green's function to which the result of the FFT is given.
 * @param inG
 *  [GF] Frequency Green's function, which is to be Fourier tranformed to time.
 * @param set_dt
 *  [int] Flag for renormalizing the grid of the output Green's function, 0: False, 1: True.
 */
void fft_solver::to_time(GF &outG, const GF &inG, int set_dt)
{
	CNTR_ASSERT_EQ(NESS_ASSERT_0, inG.size1_ , outG.size1_, __PRETTY_FUNCTION__);
	CNTR_ASSERT_EQ(NESS_ASSERT_0, inG.size2_ , outG.size2_, __PRETTY_FUNCTION__);
	CNTR_ASSERT_EQ(NESS_ASSERT_0, inG.ngrid_ , this -> nfreq_, __PRETTY_FUNCTION__);
	CNTR_ASSERT_EQ(NESS_ASSERT_0, outG.ngrid_ , this -> Nft_ / 2 + 1, __PRETTY_FUNCTION__);
	CNTR_ASSERT_EQ(NESS_ASSERT_0, outG.gf_type_ , time_gf, __PRETTY_FUNCTION__);
	CNTR_ASSERT_EQ(NESS_ASSERT_0, inG.gf_type_ , freq_gf, __PRETTY_FUNCTION__);

	//int HALF_N = Nft_ / 2 + 1;
	//int HALF_n = nfreq_ / 2 + 1;
	for(int idx1 = 0; idx1 < inG.size1_; idx1++)
	{
		for(int idx2 = 0; idx2 < inG.size2_; idx2++)
		{
			//compute gret
			for(int w = 0; w < Nft_; w++) {in[w][0] = 0.0; in[w][1] = 0.0;}
			for(int w = 0; w < HALF_n; w++)
			{
				in[w][0] = ((*inG.p_ret(w, idx1, idx2)).real() - (*inG.p_ret(w, idx2, idx1)).real()) + (*inG.p_les(w, idx1, idx2)).real();
				in[w][1] = ((*inG.p_ret(w, idx1, idx2)).imag() + (*inG.p_ret(w, idx2, idx1)).imag()) + (*inG.p_les(w, idx1, idx2)).imag();

			}
			for(int w = HALF_n; w < nfreq_; w++)
{
				in[Nft_ - nfreq_ + w][0] = ((*inG.p_ret(w, idx1, idx2)).real() - (*inG.p_ret(w, idx2, idx1)).real()) + (*inG.p_les(w, idx1, idx2)).real();
				in[Nft_ - nfreq_ + w][1] = ((*inG.p_ret(w, idx1, idx2)).imag() + (*inG.p_ret(w, idx2, idx1)).imag()) + (*inG.p_les(w, idx1, idx2)).imag();
}

			fftw_execute(plan_for);
			for(int w = 0; w < HALF_N; w++)
			{
				(*outG.p_ret(w, idx1, idx2)) = out[w][0] + NESS_II * out[w][1];
			}

			//compute glss
			for(int w = 0; w < Nft_; w++) {in[w][0] = 0.0; in[w][1] = 0.0;}
			for(int w = 0; w < HALF_n; w++)
			{
				in[w][0] = (*inG.p_les(w, idx1, idx2)).real();
				in[w][1] = (*inG.p_les(w, idx1, idx2)).imag();

			}
			for(int w = HALF_n; w < nfreq_; w++)
			{
				in[Nft_ - nfreq_ + w][0] = (*inG.p_les(w, idx1, idx2)).real();
				in[Nft_ - nfreq_ + w][1] = (*inG.p_les(w, idx1, idx2)).imag();
			}
			fftw_execute(plan_for);
			for(int w = 0; w < HALF_N; w++)
			{
				(*outG.p_les(w, idx1, idx2)) = out[w][0] + out[w][1] * NESS_II;
				(*outG.p_ret(w, idx1, idx2)) -= (*outG.p_les(w, idx1, idx2));
			}
			(*outG.p_ret(0, idx1, idx2)) *= 0.5;
		}
	}
	outG.smul(inG.dgrid_ / (2.0 * M_PI));
	if(set_dt) outG.reset_grid(2.0 * M_PI / (inG.dgrid_ * Nft_));
}
};
#endif
