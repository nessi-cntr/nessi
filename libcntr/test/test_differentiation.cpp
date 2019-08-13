#include "catch.hpp"
#include "cntr.hpp"
// #include "nca_impurity_model.hpp"


#define CPLX std::complex<double>  
#define GREEN cntr::herm_matrix<double>


TEST_CASE("Test differentiation","[Differetiation]"){
  
   int nt=100;
   int ntau=400;
   int kt=5;
   int tstp;
   
   double h=0.02;
   double beta=10;
   double err=0.0,eps=1e-7;
   
   cntr::herm_matrix<double> Gbethe;
   Gbethe=cntr::herm_matrix<double>(nt,ntau,1,-1);
   cntr::herm_matrix<double> Conv;
   Conv=cntr::herm_matrix<double>(nt,ntau,1,-1);
   
   cntr::herm_matrix<double> dG1;
   dG1=cntr::herm_matrix<double>(nt,ntau,1,-1);
   cntr::herm_matrix<double> dG2;
   dG2=cntr::herm_matrix<double>(nt,ntau,1,-1);
   
   cntr::herm_matrix<double> dGnew1;
   dGnew1=cntr::herm_matrix<double>(nt,ntau,1,-1);
   cntr::herm_matrix<double> dGnew2;
   dGnew2=cntr::herm_matrix<double>(nt,ntau,1,-1);
   
   cntr::green_equilibrium_bethe(Gbethe,beta, h);
   
   
   for(tstp=-1; tstp<=nt;tstp++)cntr::deriv1_timestep(tstp,dG1,Gbethe,Gbethe, integration::I<double>(kt), beta,h);
   for(tstp=-1; tstp<=nt;tstp++)cntr::deriv2_timestep(tstp,dG2,Gbethe,Gbethe, integration::I<double>(kt), beta,h);
   cntr::convolution(Conv,Gbethe,Gbethe,Gbethe,Gbethe,integration::I<double>(kt),beta,h);
   
	for(tstp=-1; tstp<=nt;tstp++){
		      cntr::deriv1_timestep(tstp,dGnew1,Gbethe,Gbethe, integration::I<double>(kt), beta,h);
		      // cout << "tstp: " << tstp<<endl;
		      err=cntr::distance_norm2(tstp,dGnew1,Conv);
		      REQUIRE(err<eps);
		      // cout << " |Conv- dGnew1| " << err;
		      err=cntr::distance_norm2_ret(tstp,dGnew1,Conv);
		      REQUIRE(err<eps);
		      // cout << " |Conv- dGnew1|ret " << err;
		      err=cntr::distance_norm2_tv(tstp,dGnew1,Conv);
		      REQUIRE(err<eps);		
		      // cout << " |Conv- dGnew1|tv " << err;
		      err=cntr::distance_norm2_les(tstp,dGnew1,Conv);
		      REQUIRE(err<eps);		
		      // cout << " |Conv- dGnew1|les " << err;
		      // cout << endl;
		      
		      
		      cntr::deriv2_timestep(tstp,dGnew2,Gbethe,Gbethe, integration::I<double>(kt), beta,h);
		      err=cntr::distance_norm2(tstp,dGnew2,Conv);
		      REQUIRE(err<eps);	
		      // cout << " |Conv- dGnew2| " << err;
		      err=cntr::distance_norm2_ret(tstp,dGnew2,Conv);
		      REQUIRE(err<eps);		
		      // cout << " |Conv- dGnew2|ret " << err;
		      err=cntr::distance_norm2_tv(tstp,dGnew2,Conv);
		      REQUIRE(err<eps);
		      // cout << " |Conv- dGnew2|tv " << err;
		      err=cntr::distance_norm2_les(tstp,dGnew2,Conv);
		      REQUIRE(err<eps);		
		      // cout << " |Conv- dGnew2|les " << err;
		      // cout << endl;
	}
	// cout<<"*------------------------------------*"<<endl;
	// dGnew1.print_to_file("dGnew1.out");
	// dGnew2.print_to_file("dGnew2.out");

	// dG1.print_to_file("dG1.out");
	// dG2.print_to_file("dG2.out");
	
	// Conv.print_to_file("Conv.out");
  // return 0;
}