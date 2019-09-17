#include "cntr_mpitools_extern_templates.hpp"
#include "cntr_mpitools_impl.hpp"

namespace cntr{
#if CNTR_USE_MPI == 1
	///@private
	template void Reduce_timestep<double>(int tstp, int root, herm_matrix_timestep<double> &Gred, 
		herm_matrix_timestep<double> G);
	///@private
	template void Reduce_timestep<double>(int tstp, int root, herm_matrix<double> &Gred, 
		herm_matrix_timestep<double> G);
	///@private
	template void Reduce_timestep<double>(int tstp, int root, herm_matrix_timestep<double> &Gred, 
		herm_matrix<double> G);
	///@private
	template void Reduce_timestep<double>(int tstp, int root, herm_matrix<double> &Gred, 
		herm_matrix<double> G);
#endif

} // namespace cntr