#include "ness_GF_pair_decl.hpp"


namespace ness
{

 
  /** \brief <b> Standard Constructor of a Green's function pair </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *   \par Purpose
   * <!-- ========= -->
   * > Constructs a time `GF` and its frequency counterpart with all grid, shape and type parameters specified manually,
   * > and the Fourier solver object, enabling a transform between the `GF` objects.
   *
   * <!-- ARGUMENTS
   *      ========= -->
   *
   * @param h
   * > [double] Time grid spacing
   * @param nt
   * > [int] Number of time grid points
   * @param size
   * > [int] First orbital dimension, `GF` matrix rows
   * @param sign
   * > [int] Sign of particle type, +1 for boson,-1 for fermion.
   * @param  fft
   * > [fft_solver] Fourier tranform solver object, contains Fourier transform parameters
   */
  GF_pair::GF_pair(double h, int nt, int size, int sign, fft_solver & fft):
  time(h, nt, size, sign, time_gf),
  freq( M_PI / ((nt-1) * h) , 2 * (nt-1) / 3, size, sign, freq_gf),
    fft_(&fft)
  {
      assert(nt == fft_->nt_);
      assert(time.ngrid_ == fft_->nt_);
      assert(freq.ngrid_ == fft_->nfreq_);
  }

  /** \brief <b> Default constructor </b>*/
  GF_pair::GF_pair()
  {
    fft_ = nullptr;
  }
  
  /** \brief <b> Destructor </b> */
  GF_pair::~GF_pair() {}

  /** \brief <b> Copy constructor </b>
   *
   * <!-- ====== DOCUMENTATION ====== -->
   *
   *   \par Purpose
   * <!-- ========= -->
   *
   * > Constructs an empty `GF_pair` with all parameters given by another `GF_pair`,
   * > so that the new object has the same shape and type as the parameter `GF_pair`.
   *
   * <!-- ARGUMENTS
   *      ========= -->
   *
   * @param  gp
   * > [GF_pair ] Green's function pair to copy parameters from
   */
  GF_pair::GF_pair(const GF_pair & gp):
    newest(gp.newest), freq(gp.freq), time(gp.time), fft_(gp.fft_)
  {}
  
  GF_pair& GF_pair::operator=(const GF_pair & gp)
  {
    newest = gp.newest;
    freq = gp.freq;
    time = gp.time;
    fft_ = gp.fft_;
    
    return *this;
  }
  
  /** \brief <b> Clear time and frequency Green's function in the `GF_pair` </b> */
  void GF_pair::clear()
  {
    time.clear();
    freq.clear();
  }
  
  /** \brief <b> Transform frequency `GF` to time and write in the corresponding time `GF` member in the `GF_pair` </b> */
  void GF_pair::to_time() {fft_ -> to_time(time, freq); newest = time_gf;}
  /** \brief <b> Transform time `GF` to frequency and write in the corresponding frequency `GF` member in the `GF_pair` </b> */
  void GF_pair::to_freq() {fft_ -> to_freq(freq, time); newest = freq_gf;}

}
