#include <sys/stat.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cmath>
#include <cstring>
#include <chrono>

// contour library headers
#include "cntr/cntr.hpp"
#include "cntr/utils/read_inputfile.hpp"

// local headers to include
#include "formats.hpp"

using namespace std;

#define cplx std::complex<double>
// -----------------------------------------------------------------------
complex<double> GregoryIntegral(int n, double h, int k, integration::Integrator<double> &I){
  cplx iu(0.0,1.0);
  int nmax;
  double w;
  cplx exp_integ;

  exp_integ = 0.0*iu;

  nmax = (n > k ? n : k);

  if (n <= 2*k + 2){
    for(int j=0; j<=nmax; j++){
      w = I.gregory_weights(n, j);
      exp_integ += w*(cos(j*h) + iu*sin(j*h));
    }
  } else{
    for(int j=0; j<=k; j++){
      w = I.gregory_omega(j);
      exp_integ += w*(cos(j*h) + iu*sin(j*h));
    }
    for(int j=k+1; j<n-k; j++){
      exp_integ += (cos(j*h) + iu*sin(j*h));
    }
    for(int j=n-k; j<=n; j++){
      w = I.gregory_omega(n-j);
      exp_integ += w*(cos(j*h) + iu*sin(j*h));
    }
  }

  exp_integ = h*exp_integ;
  return exp_integ;

}
// -----------------------------------------------------------------------
cplx GregoryIntegral(int n, double h, vector<cplx> &f, integration::Integrator<double> &I){
  return h * I.integrate(f, n);
}

//==============================================================================
//         main program
//==============================================================================
int main(int argc,char *argv[]){
  //..................................................
  //                input
  //..................................................
  int npts,k;
  double h;
  std::complex<double> iu(0.0,1.0);
  //..................................................
  cout << endl;
  cout << " reading input file ..." << endl;
  cout << endl;
  try{
    //============================================================================
    //                          (II) READ INPUT
    //============================================================================
    {
      if(argc<2) throw("COMMAND LINE ARGUMENT MISSING");

      find_param(argv[1],"__npts=",npts);
      find_param(argv[1],"__h=",h);
      find_param(argv[1],"__k=",k);

      if (argc < 3) {
	 // Tell the user how to run the program
        std::cerr << " Please provide a prefix for the output files. Exiting ..." << std::endl;
        return 1;
      }

    }

    {
      int n;
      complex<double> integ_approx,integ_exact;
      double err;
      ofstream fout;
      fout.open(argv[2]);

      vector<cplx> fx(npts+1);
      for (n=0; n<= npts; n++) fx[n] = cos(n*h) + iu*sin(n*h);

      for (n=0; n<= npts; n++){
       integ_approx = GregoryIntegral(n, h, k, integration::I<double>(k));
       // integ_approx = GregoryIntegral(n, h, fx, integration::I<double>(k));
       integ_exact = (cos(n*h) + iu*sin(n*h) - 1.0)/iu;
       err = abs(integ_approx-integ_exact);
       fout << n*h << "  " << err << endl;
     }

     fout.close();
   }


  } // try
  catch(char *message){
    cerr << "exception\n**** " << message << " ****" << endl;
    cerr << " No input file found. Exiting ... " << endl;
  }
  catch(...){
    cerr << " No input file found. Exiting ... " << endl;
  }
  return 0;
}
//==============================================================================
