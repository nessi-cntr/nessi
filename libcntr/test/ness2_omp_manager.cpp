#include "catch.hpp"
#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>
#include "ness2.hpp"

using namespace ness2;

/*///////////////////////////////////////////////////////////////////////////////////////

  This is a minimal test program to test the initialization of the OMP manager and that the omp related functions work (only the compatibility).

///////////////////////////////////////////////////////////////////////////////////////*/

//------------------------------------------------------------------------------

TEST_CASE("OMP manager initialization","[FFT_OMP_Manager]"){
    const int Nt = 1024;
    const int size = 1;
    const int nomp = 2;
    const double tol=1.0e-10;
    double err = 0.0;
    std::vector<std::complex<double>> MA(Nt);
    std::vector<std::complex<double>> MB(Nt);
    
    // FFT of B without OMP
    {
        fft_array B(Nt,size,FFTW_ESTIMATE);
        for (int t=0;t<Nt;t++){
            cdmatrix m(size,size);
            for(int i=0;i<size;i++){
                for(int j=0;j<size;j++){
                    double w=(0.33*i+0.242*j)*0.01;
                    m(i,j) = cdouble(cos(w*t),sin(w*t));
                }
            }
            B.set_element(t,m,fft_domain::time);
        }
        B.fft_to_freq();
        for (int t=0;t<Nt;t++){
            B.get_element(t,MB[t], fft_domain::freq);
        }
    }
    
    // FFT of A with OMP
    FFT_OMP_Manager::initialize();
    FFT_OMP_Manager::set_threads_for_new_plans(nomp);
    {
        fft_array A(Nt,size,FFTW_ESTIMATE);
        // fill some data:
        for (int t=0;t<Nt;t++){
            cdmatrix m(size,size);
            for(int i=0;i<size;i++){
                for(int j=0;j<size;j++){
                    double w=(0.33*i+0.242*j)*0.01;
                    m(i,j) = cdouble(cos(w*t),sin(w*t));
                }
            }
            A.set_element(t,m,fft_domain::time);
        }
        A.fft_to_freq();
        for (int t=0;t<Nt;t++){
            A.get_element(t,MA[t], fft_domain::freq);
        }
    }
    FFT_OMP_Manager::finalize();
    for (int i = 0; i < Nt; i++) {
        auto diff = MA[i] - MB[i];
        err += std::norm(diff);
    }
    REQUIRE(err/Nt<tol);
}
