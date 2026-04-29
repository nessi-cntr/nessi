#ifndef NESS2_GLOBAL_SETTINGS_H
#define NESS2_GLOBAL_SETTINGS_H

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <stdexcept>

#include "eigen3/Eigen/Dense"
#include <fftw3.h>

#include "cntr_global_settings.hpp"

namespace ness2 {

#define NESS_ASSERT_0 1
  
  /** \brief Define complex data type ... same as in cntr:: */
  typedef std::complex<double> cplx;
  typedef std::complex<double> cdouble;
    
  const cplx NESS_II(0.0, 1.0);
  enum class fft_domain { time, freq };  
  //enum class fft_integral_method {  trapez,  cubic,  adaptive};
  enum fft_integral_method {  FFT_TRAPEZ,  FFT_CUBIC, FFT_ADAPTIVE};
  
} // namespace ness2

#define CNTR_ASSERT_EQ(d1, a, b, d2) (assert((a) == (b)))
#define CNTR_ASSERT_LESEQ(d1, a, b, d2) (assert((a) <= (b)))
#define CNTR_ASSERT_LESEQ_3(d1, a, b, c, d2) (assert((a) <= (b) && (b) <= (c)))

#endif // NESS_GLOBAL_SETTINGS_H


