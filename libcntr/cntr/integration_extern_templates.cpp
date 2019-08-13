#include "integration.hpp"

namespace integration {

template class Integrator<double>;
template Integrator<double> &I<double>(int k);

}  // namespace integration
