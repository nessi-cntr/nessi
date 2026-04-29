#ifndef NESS_GLOBAL_SETTINGS_H
#define NESS_GLOBAL_SETTINGS_H
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstring>

#include "eigen3/Eigen/Dense"
#include <fftw3.h>


namespace ness{
#define NESS_ASSERT_0 1

/** \brief <b> Define complex data type </b> */

typedef std::complex<double> cplx;

/** \brief <b> define fermion as sign -1 and boson as sign +1 </b> */
enum {fermion = -1, boson = 1};

using Eigen::MatrixXcd;
using Eigen::VectorXcd;
/** \brief <b> Own complex eigen matric type </b> */
typedef MatrixXcd cdmatrix;
/** \brief <b> Own complex eigen vector type </b> */
typedef VectorXcd cdvector;
// awkward, to be conflicting with cntr
/** \brief <b> Imaginary unit </b> */
const cplx NESS_II(0.0, 1.0);

} //end of namespace

#define CNTR_ASSERT_EQ(d1, a,b, d2) (assert((a)==(b)))
#define CNTR_ASSERT_LESEQ(d1, a,b, d2) (assert((a)<=(b)))
#define CNTR_ASSERT_LESEQ_3(d1, a,b,c, d2) (assert((a)<=(b) && (b)<=(c)));



#endif
