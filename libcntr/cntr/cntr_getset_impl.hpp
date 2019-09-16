#ifndef CNTR_GETSET_IMPL_H
#define CNTR_GETSET_IMPL_H

#include "eigen_typedef.h"
#include "eigen_map.hpp"
#include "cntr_herm_matrix_decl.hpp"
#include "cntr_herm_matrix_timestep_decl.hpp"
#include "cntr_herm_matrix_timestep_view_decl.hpp"
#include "cntr_getset_decl.hpp"



namespace cntr {


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