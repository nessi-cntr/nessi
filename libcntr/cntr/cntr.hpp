/*********************************************************
*
*  Martin Eckstein, 2010
*  Contour Green funtion class
*
*********************************************************/
#ifndef GREEN_CNTR_FULL
#define GREEN_CNTR_FULL

/* #######################################################################################
#
#                single-particle Greenfunctions G(t,t')
#
##########################################################################################

     currently there exist three closely related Greenfunction types ...
	 
	 -  type T is double or float
     -  timesteps 0... ntau on imag contour,
	    0...nt on upper/lower real contour.
		for nt=-1, only the matsubara sector is stored
     -  each element has element_size_=size2_*size1_ complex<T> numbers
	    (only square matices allowed (size1_=size2_)
     -  strored elements:
        gtr(t,t'),les(t,t'),tv(t,t'),vt(t,t')
		mat(tau),mat(-tau)
		ret(t,t')=gtr-les  (stored for all t,t')
		
####################################################################################### */



#include "cntr_decl.hpp"
#include "cntr_impl.hpp"
#ifndef CNTR_NO_EXTERN_TEMPLATES
#include "cntr_extern_templates.hpp"
#endif

#endif // GREEN_CNTR_FULL
