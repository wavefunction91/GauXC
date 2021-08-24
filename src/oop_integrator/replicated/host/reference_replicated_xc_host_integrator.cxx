#include "reference_replicated_xc_host_integrator_exc_vxc.hpp"
#include "reference_replicated_xc_host_integrator_exc_grad.hpp"
 
namespace GauXC  {
namespace detail {

template <typename ValueType>
ReferenceReplicatedXCHostIntegrator<ValueType>::~ReferenceReplicatedXCHostIntegrator() noexcept = default;

template class ReferenceReplicatedXCHostIntegrator<double>;

}
}
