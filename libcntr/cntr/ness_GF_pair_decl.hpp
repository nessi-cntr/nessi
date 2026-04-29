#include "ness_GF_decl.hpp"
#include "ness_fft_decl.hpp"

#ifndef GF_PAIR
#define GF_PAIR
namespace ness
{
/** \brief <b> Class `GF_pair` represents Green's functions in both time or frequency domain. </b>
 *
 * <!-- ====== DOCUMENTATION ====== -->
 *
 *  \par Purpose
 * <!-- ========= -->
 *
 * > This class contains two `GF` data structures for representing complex matrix valued Green's functions
 * > on a time and frequency grid (retarded and lesser component), and provides a FFT function given by the helper class `fft_solver` for calculating Fourier integral transforms.
 *
 */
class GF_pair
{

public:
/** \brief <b> Freqeuncy  `GF` object </b> */
GF freq, time; /**< \brief <b> Time  `GF` object </b> */
/** \brief <b> Pointer to Fourier solver object </b> */
fft_solver * fft_;
/** \brief <b> Green's function type flag, 0: time, 1frequency </b> */
int newest;

  GF_pair(double h, int nt, int size, int sign, fft_solver & fft);
/** \brief <b> Default constructor </b>*/
  GF_pair();

/** \brief <b> Destructor </b> */
  ~GF_pair();

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
* @param & gp
* > [GF_pair ] Green's function pair to copy parameters from
*/
  GF_pair(const GF_pair & gp);

  GF_pair& operator=(const GF_pair & gp);

/** \brief <b> Clear time and frequency Green's function in the `GF_pair` </b> */
  void clear();

/** \brief <b> Tranform frequency `GF` to time and write in the corresponding time `GF` member in the `GF_pair` </b> */
  void to_time();
/** \brief <b> Tranform time `GF` to frequency and write in the correspondingfrequency `GF` member in the `GF_pair` </b> */
  void to_freq();

};
/** \brief <b> Define a vector of `GF_pair` objects </b> */
typedef std::vector<GF_pair> GF_pairs;
};

#endif
