#include <sys/stat.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <cstring>

// contour library headers
#include "cntr/cntr.hpp"
#include "cntr/utils/read_inputfile.hpp"

// local headers to include
// #include "formats.hpp"

using namespace std;

//==============================================================================
//         main program
//==============================================================================
int main(int argc,char *argv[]){
  
  cdmatrix A;
  A.setRandom(10,10);
  
  cntr::hodlr_box<double> Abox(A,1);

  std::cout  << " Singular values " << Abox.singular() << " " << Abox.svdtol() << std::endl;

  std::cout  << " U " << Abox.U()  << std::endl;
  std::cout  << " U " << Abox.U().cols() << " " << Abox.U().rows()  << std::endl;

  std::cout  << " V "  << Abox.V().cols() << " " << Abox.V().rows()  << std::endl;

  cntr::herm_matrix_hodlr<double> G(64,3,0.1);



  return 0;
}
