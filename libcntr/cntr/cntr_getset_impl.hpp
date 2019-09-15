#ifndef CNTR_GETSET_IMPL_H
#define CNTR_GETSET_IMPL_H

#include "eigen_typedef.h"
#include "eigen_map.hpp"
#include "cntr_herm_matrix_decl.hpp"
#include "cntr_herm_matrix_timestep_decl.hpp"
#include "cntr_herm_matrix_timestep_view_decl.hpp"
#include "cntr_getset_decl.hpp"

namespace cntr {

/// @private	
template <typename T>
void map_component(int size1, int size2, std::complex<T> *ptr, cdmatrix &M){
	assert(size1 == size2);
	// M.resize(size1, size2);
	switch (size1){
		case 1:
			M = element_map<1, 1>(size1, size2, ptr);
			break;
		case 2:
			M = element_map<2, 2>(size1, size2, ptr);
			break;
		case 3:
			M = element_map<3, 3>(size1, size2, ptr);
			break;
		case 4:
			M = element_map<4, 4>(size1, size2, ptr);
			break;
		case 5:
			M = element_map<5, 5>(size1, size2, ptr);
			break;
		case 6:
			M = element_map<6, 6>(size1, size2, ptr);
			break;
		case 7:
			M = element_map<7, 7>(size1, size2, ptr);
			break;
		case 8:
			M = element_map<8, 8>(size1, size2, ptr);
			break;
		default:
			M = element_map<-1, -1>(size1, size2, ptr);
	}
}


//---------------------------------------------------------------------
//-------                  herm_matrix                          -------
//---------------------------------------------------------------------
#include "cntr_getset_herm_matrix_inc.hpp"
//---------------------------------------------------------------------
//-------             herm_matrix_timestep                      -------
//---------------------------------------------------------------------
#include "cntr_getset_herm_matrix_timestep_inc.hpp"
//---------------------------------------------------------------------
//-------           herm_matrix_timestep_view                   -------
//---------------------------------------------------------------------
#include "cntr_getset_herm_matrix_timestep_view_inc.hpp"




} // namespace cntr


#endif  // CNTR_GETSET_IMPL_H