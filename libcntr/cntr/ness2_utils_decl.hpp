#ifndef NESS2_UTILS_H
#define NESS2_UTILS_H

#include "ness2.hpp"
#include <fstream>


namespace ness2
{

  /** \brief <b> Compute Euclidean distance between two `fft_array` objects in time or frequency domain </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *   \par Purpose
   * <!-- ========= -->
   *
   * > Computes the L2 norm (Euclidean distance)
   * > \f$ \|A - B\|_2 = \sqrt{ \sum_i |A_i - B_i|^2 } \f$
   * > between two `fft_array` objects in either time or frequency domain.
   *
   * <!-- ARGUMENTS
   *      ========= -->
   *
   * @param &A  
   *  First `fft_array` object  
   * @param &B  
   *  Second `fft_array` object  
   * @param type  
   *  `fft_domain::time` or `fft_domain::freq` – determines which internal data to compare  
   *
   * @returns  
   *  Euclidean norm of the difference
   *
   * @throws std::invalid_argument  
   *  If array sizes or shapes do not match
   *
   */
  double  distance_norm2(fft_array &A,fft_array &B,fft_domain type);
  

  /** \brief <b> Compute Euclidean distance between two `herm_matrix_ness` objects </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *   \par Purpose
   * <!-- ========= -->
   *
   * > Computes the combined L2 norm for retarded and lesser components:
   * > \f$ \|A - B\|_2 = \sqrt{ \|A^{ret} - B^{ret}\|_2^2 + \|A^{les} - B^{les}\|_2^2 } \f$
   *
   * <!-- ARGUMENTS
   *      ========= -->
   *
   * @param &A  
   *  First `herm_matrix_ness` object  
   * @param &B  
   *  Second `herm_matrix_ness` object  
   * @param type  
   *  `fft_domain::time` or `fft_domain::freq` – domain to compare  
   *
   * @returns  
   *  Total Euclidean norm of the difference across components
   *
   */
  double  distance_norm2(herm_matrix_ness &A,herm_matrix_ness &B,fft_domain type);


  /** \brief <b> Calculate the bubble diagram of first type: \f$ C_{c_1,c_2} (t)= i A_{a_1,a_2} (t) B_{b_2,b_1} (-t) \f$ </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *   \par Purpose
   * <!-- ========= -->
   *
   * > Evaluate Green's function C given by the bubble diagram equation
   * > \f$ C_{c_1,c_2} (t)= i A_{a_1,a_2} (t) B_{b_2,b_1} (-t) \f$.
   * > The evaluation is done for all time points and for the lesser and retarded component.
   *
   * <!-- ARGUMENTS
   *      ========= -->
   *
   * @param &C
   *   Green's function to store result in
   * @param &c1
   *  First orbital index of C
   * @param &c2
   *  Second orbital index of C
   * @param &A
   *   Time Green's function
   * @param &a1
   *  First orbital index of A
   * @param &a2
   *  Second orbital index of A
   * @param &B
   *   Time Green's function
   * @param &b1
   *  First orbital index of B
   * @param &b2
   *  Second orbital index of B
   */
  void Bubble1_ness(herm_matrix_ness &C,int c1,int c2,herm_matrix_ness &A,int a1,int a2,herm_matrix_ness &B,int b1,int b2);
  void Bubble1_ness(herm_matrix_ness &C, herm_matrix_ness &A, herm_matrix_ness &B);
  
  /** \brief <b> Calculate the bubble diagram of second type:  \f$ C_{c_1,c_2} (t)= i A_{a_1,a_2} (t) B_{b_1,b_2} (t) \f$ </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *   \par Purpose
   * <!-- ========= -->
   *
   * > Evaluate Green's function C given by the bubble diagram equation
   * > \f$ C_{c_1,c_2} (t)= i A_{a_1,a_2} (t) B_{b_1,b_2} (t) \f$.
   * > The evaluation is done for all time points and for the lesser ad retarded component.
   *
   * <!-- ARGUMENTS
   *      ========= -->
   *
   * @param &C
   *   Time Green's function to store result in
   * @param &c1
   *  First orbital index of C
   * @param &c2
   *  Second orbital index of C
   * @param &A
   *   Time Green's function
   * @param &a1
   *  First orbital index of A
   * @param &a2
   *  Second orbital index of A
   * @param &B
   *   Time Green's function
   * @param &b1
   *  First orbital index of B
   * @param &b2
   *  Second orbital index of B
   */
  void Bubble2_ness(herm_matrix_ness &C,int c1,int c2,herm_matrix_ness &A,int a1,int a2,herm_matrix_ness &B,int b1,int b2);
  void Bubble2_ness(herm_matrix_ness &C, herm_matrix_ness &A, herm_matrix_ness &B);   


  /**
   * @brief Convert CNTR Green's function `cntr::herm_matrix` to NESS `ness2::herm_matrix_ness` by sampling a given time slice.
   *
   * Read \f$ G_{\rm ness}^{R|<}(t) = G{\rm cntr}^{R|<}(tstp, tstp - t) \f$ from CNTR object at time slice `tstp`,
   * and write to NESS object.
   *
   * Values at negative times \f$ t < 0 \f$ are set by Hermitian symmetry:
   * \f$ G^<(-t) = -[G^<(t)]^\dagger \f$
   *
   * If `tstp` is smaller than the maximum time in Gness, the tail is left zero.
   *
   * @tparam GcntrType Type of the contour Green's function (cntr::herm_matrix,cntr::herm_matrix_timestep, cntr::herm_matrix_moving, cntr::herm_matrix_timestep_moving)
   * @param[out] Gness Target steady-state Green's function (time domain filled)
   * @param[in] Gcntr Source two-time Green's function
   * @param[in] tstp Time slice in `Gcntr` to sample from (defaults to max time)
   * @throws std::runtime_error If dimensions mismatch
   */
  template <typename T>
  void cntr2ness(herm_matrix_ness& Gness, cntr::herm_matrix<T>& Gcntr, int tstp = -1) {
    int Nft_half = Gness.Nft_ / 2;
    int size = Gness.size1_;
    int nt=Gcntr.nt(); 
    if (size != Gcntr.size1()) {
      throw std::runtime_error("cntr2ness: matrix size mismatch");
    }
    if (tstp>nt || nt<0) {
      throw std::runtime_error("cntr2ness: wrong dimensiuons in Gcntr");
    }    
    if (tstp < 0) tstp = Gcntr.nt(); // default to max physical time    
    Gness.set_zero(fft_domain::time);
    cdmatrix M;
    for (int t = 0; t <= std::min(tstp, Nft_half - 1); ++t) {
      Gcntr.get_ret(tstp, tstp - t, M);
      Gness.set_ret(t, M, fft_domain::time);
    }
    for (int t = -std::min(tstp, Nft_half - 1); t <= 0; ++t) {    
      Gcntr.get_les(tstp+t, tstp, M);
      Gness.set_les(t, M, fft_domain::time);
      if (t < 0) Gness.set_les(-t, -M.adjoint(), fft_domain::time);
    }
  }
  template <typename T>
  void cntr2ness(herm_matrix_ness& Gness, cntr::herm_matrix_timestep<T>& Gcntr, int tstp = -1) {
    int Nft_half = Gness.Nft_ / 2;
    int size = Gness.size1_;
    tstp=Gcntr.tstp(); 
    if (size != Gcntr.size1()) {
      throw std::runtime_error("cntr2ness: matrix size mismatch");
    }
    if (tstp<0) {
      throw std::runtime_error("cntr2ness: wrong dimensiuons in Gcntr");
    }    
    Gness.set_zero(fft_domain::time);
    cdmatrix M;
    for (int t = 0; t <= std::min(tstp, Nft_half - 1); ++t) {
      Gcntr.get_ret(tstp, tstp - t, M);
      Gness.set_ret(t, M, fft_domain::time);
    }
    for (int t = -std::min(tstp, Nft_half - 1); t <= 0; ++t) {    
      Gcntr.get_les(tstp+t, tstp, M);
      Gness.set_les(t, M, fft_domain::time);
      if (t < 0) Gness.set_les(-t, -M.adjoint(), fft_domain::time);
    }
  }

  //@private
  template <typename T>
  void cntr2ness(herm_matrix_ness& Gness, cntr::herm_matrix_timestep_moving_view<T>& Gcntr, int tstp = -1) {
    int Nft_half = Gness.Nft_ / 2;
    int size = Gness.size1_;
    tstp=Gcntr.tc(); 
    if (size != Gcntr.size1()) {
      throw std::runtime_error("cntr2ness: matrix size mismatch");
    }
    Gness.set_zero(fft_domain::time);
    cdmatrix M;
    for (int t = 0; t <= std::min(tstp, Nft_half - 1); ++t) {
      Gcntr.get_ret(t, M);
      Gness.set_ret(t, M, fft_domain::time);
      Gcntr.get_les(t, M);
      Gness.set_les(t, M, fft_domain::time);
      if (t > 0) Gness.set_les(-t, -M.adjoint(), fft_domain::time);
    }
  }
  template <typename T>
  void cntr2ness(herm_matrix_ness& Gness, cntr::herm_matrix_timestep_moving<T>& Gcntr, int tstp = -1) {
    cntr:: herm_matrix_timestep_moving_view<T> view(Gcntr);
    cntr2ness(Gness,view,tstp);
  }
  template <typename T>
  void cntr2ness(herm_matrix_ness& Gness, cntr::herm_matrix_moving<T>& Gcntr, int tstp = -1) {
    cntr::herm_matrix_timestep_moving_view<T> view(Gcntr,0);
    cntr2ness(Gness,view,tstp);
  }

  /**
   * @brief Convert a steady-state Green's function into a CNTR Green's function, assuming tranlstional invariance in time.
   *
   * For each pair \f$ (t, t') \f$, if the relative time \f$ t - t' \f$ exists in the NESS domain,
   * use \f$ G_{cntr}(t,t') =  G_{ness}(t - t') \f$. Else, set to zero.
   *
   * @tparam GcntrType Any writable two-time Green's function (cntr::herm_matrix,cntr::herm_matrix_timestep, cntr::herm_matrix_moving, cntr::herm_matrix_timestep_moving)
   * @param[out] Gcntr Target two-time function
   * @param[in] Gness Source steady-state function
   */
  template <typename T>
  void ness2cntr(cntr::herm_matrix<T>& Gcntr, const herm_matrix_ness& Gness) {
    const int nt = Gcntr.nt();
    const int size = Gcntr.size1();
    const int Nft_half = Gness.Nft_ / 2;
    if (size != Gness.size1_) {
      throw std::runtime_error("ness2cntr: matrix size mismatch");
    }
    cdmatrix M;
    for (int t = 0; t <= nt; ++t) {
      Gcntr.set_timestep_zero(t);
      for (int t1 = 0; t1 <= t; ++t1) {
	int dt = t - t1;
	if (dt < Nft_half) {
	  Gness.get_ret(dt, M, fft_domain::time);
	  Gcntr.set_ret(t, t1, M);
	  Gness.get_les(dt, M, fft_domain::time);
	  Gcntr.set_les(t1, t, -M.adjoint());
	}
      }
    }
  }
  template <typename T>
  void ness2cntr(cntr::herm_matrix_timestep<T>& Gcntr, const herm_matrix_ness& Gness) {
    const int t = Gcntr.tstp();
    const int size = Gcntr.size1();
    const int Nft_half = Gness.Nft_ / 2;
    if (size != Gness.size1_) {
      throw std::runtime_error("ness2cntr: matrix size mismatch");
    }
    cdmatrix M;
    Gcntr.set_timestep_zero(t);
    for (int t1 = 0; t1 <= t; ++t1) {
      int dt = t - t1;
      if (dt < Nft_half) {
	Gness.get_ret(dt, M, fft_domain::time);
	Gcntr.set_ret(t, t1, M);
	Gness.get_les(dt, M, fft_domain::time);
	Gcntr.set_les(t1, t, -M.adjoint());
      }
    }
  }
  template <typename T>
  void ness2cntr(cntr::herm_matrix_timestep_moving_view<T>& Gcntr, const herm_matrix_ness& Gness) {
    const int tc = Gcntr.tc();
    const int size = Gcntr.size1();
    const int Nft_half = Gness.Nft_ / 2;
    if (size != Gness.size1_) {
      throw std::runtime_error("ness2cntr: matrix size mismatch");
    }
    int dtmax=(tc<Nft_half-1 ? tc : Nft_half-1 );
    Gcntr.set_timestep_zero();
    cdmatrix M;
    for (int dt = 0; dt <= dtmax ; ++dt) {
      Gness.get_ret(dt, M, fft_domain::time);
      Gcntr.set_ret(dt, M);
      Gness.get_les(dt, M, fft_domain::time);
      Gcntr.set_les(dt,M);
    }
  }
  template <typename T>
  void ness2cntr(cntr::herm_matrix_moving<T>& Gcntr, const herm_matrix_ness& Gness) {
    const int tc = Gcntr.tc();
    for(int t=0;t<=tc;t++){
      cntr::herm_matrix_timestep_moving_view<T> view(Gcntr,t);
      ness2cntr(view,Gness);
    }
  }
  template <typename T>
  void ness2cntr(cntr::herm_matrix_timestep_moving<T>& Gcntr, const herm_matrix_ness& Gness) {
    cntr::herm_matrix_timestep_moving_view<T> view(Gcntr);
    ness2cntr(view,Gness);
  }
  
  /**
   * @brief Compute the equal-time lesser component of the product \f$ i (A \cdot B)^< \f$
   *        for two steady-state Green's functions in NESS representation.
   *
   * This routine evaluates
   * \f[
   *    M = bosefermi* i (A \cdot B)^<(t=0)
   *      = bosefermi*i \, \mathcal{F}^{-1} \left[ A^<(\omega) B^A(\omega) + A^R(\omega) B^<(\omega) \right]_{t=0}
   * \f]
   *
   * The function transforms  A and B to frequency (overwriing eisting values)
   * on the frequency grid of the NESS object, limited to the range
   *
   * @param[out] result
   *     Complex matrix of dimension `(size1 x size1)` containing the equal-time
   *     value \f$ i(A \cdot B)^<(t=0) \f$.
   * @param[in] bosefermi
   *     (+1) for bosonic functions, (-1) for fermionic functions
   * @param[in] h
   *     Time step used for the Fourier transforms.
   * @param[in] A
   *     First steady-state Green's function (`herm_matrix_ness`)
   * @param[in] B
   *     Second steady-state Green's function (`herm_matrix_ness`)
   *
   */
  void convolution_density_matrix(cdmatrix &result, int bosefermi,herm_matrix_ness &A, herm_matrix_ness &B,double h);
  
  /**
   * @brief Compute the equal-time density matrix from the lesser component of a NESS Green's function.
   *
   * This routine evaluates
   * \f[
   *    \rho = {\tt bosefermi} \cdot i \, A^<(t=0)
   * \f]
   * where \f$ A^<(t) \f$ is the lesser Green's function in the time domain, and \f$t=0\f$
   * denotes the equal-time value. 
   *
   * @param[out] result
   *     Complex matrix of dimension `(size1 x size1)` containing
   *     \f$ {\tt bosefermi} \cdot i A^<(t=0) \f$.
   * @param[in] bosefermi
   *     (+1) for bosonic functions, (-1) for fermionic functions
   * @param[in] A
   *     Steady-state Green's function (`herm_matrix_ness`) in the time domain.
   */
  void density_matrix( cdmatrix &result, int bosefermi,herm_matrix_ness &A);

  
   /**
    * @brief Downsampling of a `fft_array` to a coarser time grid, taking only every `factor` entry of the original object.
    *
    * The new length is \f$ N_{ft1} = \frac{N_{ft}}{\tt factor} \f$.
    * In frequency space, only values \f$ 0, \ldots, \frac{N_{ft1}}{2} - 1, N_{ft} - \frac{N_{ft1}}{2}, \ldots, N_{ft} - 1 \f$ of the input are written to the output object.
    * Returns a downsampled `fft_array` object.
    * It is required that \f$ N_{ft} \f$ is an integer multiple of `factor`, and \f$ N_{ft1} \f$ (the new length) is still an integer multiple of 2.
    *
    * @param in
    *     Input array (`fft_array`).
    * @param factor
    *     Coarse-graining factor, larger than one.
    */
  fft_array downsample(fft_array &in, int factor);

/**
 * @brief Downsampling of a `herm_matrix_ness` to a coarser time grid, taking only every `factor` entry of the original object.
 *
 * The new length is \f$ N_{ft1} = \frac{N_{ft}}{\tt factor} \f$.
 * In frequency space, only values \f$ 0, \ldots, \frac{N_{ft1}}{2} - 1, N_{ft} - \frac{N_{ft1}}{2}, \ldots, N_{ft} - 1 \f$ of the input are written to the output object.
 * Returns a downsampled `herm_matrix_ness` object.
 * It is required that \f$ N_{ft} \f$ is an integer multiple of `factor`, and \f$ N_{ft1} \f$ (the new length) is still an integer multiple of 2.
 *
 * @param in
 *     Input steady-state Green's function (`herm_matrix_ness`).
 * @param factor
 *     Coarse-graining factor, larger than one.
 */
  herm_matrix_ness downsample(herm_matrix_ness &in, int factor);

 /**
  * @brief Upsampling of a `herm_matrix_ness` to a finer time grid by transforming to frequency space,
  *        appending zeros until the new length \f$ N_{ft1} = N_{ft} \times factor \f$, and transforming back.
  *
  * @param h_in
  *     Initial time step.
  * @param in
  *     Input steady-state Green's function (`herm_matrix_ness`).
  * @param factor
  *     Fine-graining factor, larger than one.
  */
  herm_matrix_ness upsample(double h_in,herm_matrix_ness &in, int factor);

/**
 * @brief Automatic upsampling of a `herm_matrix_ness` by a factor of '-factor' if a negative value is passed and an automatic downsampling ba a factor of 'factor' if a positive value is passed.
 *
 * @param in
 *     Input steady-state Green's function (`herm_matrix_ness`).
 * @param factor
 *     Coarse-graining factor if positive, if negative '-factor' is fine-graining factor.
 */
  herm_matrix_ness resample(double h_in,herm_matrix_ness &in, int factor);


  
};//end of namespace

#endif // NESS2_UTILS_H
