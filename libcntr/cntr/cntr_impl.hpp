#ifndef CNTR_IMPL_H
#define CNTR_IMPL_H

#include "cntr_global_settings.hpp"

#include "integration.hpp"
#include "linalg.hpp"
#include "fourier.hpp"
#include "eigen_map.hpp"
#include "cntr_elements.hpp"

#include "cntr_matsubara_impl.hpp"

#include "cntr_function_impl.hpp"

#include "cntr_herm_matrix_timestep_view_impl.hpp"
#include "cntr_herm_matrix_timestep_impl.hpp"
#include "cntr_herm_matrix_impl.hpp"
#include "cntr_herm_matrix_hodlr_impl.hpp"

#include "cntr_herm_pseudo_impl.hpp"

#include "cntr_utilities_impl.hpp"
#include "cntr_differentiation_impl.hpp"
#include "cntr_convolution_impl.hpp"
#include "cntr_pseudo_convolution_impl.hpp"
#include "cntr_response_convolution_impl.hpp"

#include "cntr_equilibrium_impl.hpp"
#include "cntr_vie2_impl.hpp"
#include "cntr_pseudo_vie2_impl.hpp"
#include "cntr_dyson_impl.hpp"
#include "cntr_pseudodyson_impl.hpp"

#include "cntr_bubble_impl.hpp"
#include "cntr_distributed_array_impl.hpp"
#include "cntr_distributed_timestep_array_impl.hpp"

#include "cntr_getset_impl.hpp"
#if CNTR_USE_MPI == 1
#include "cntr_mpitools_impl.hpp"
#endif

#endif  // CNTR_IMPL_H
