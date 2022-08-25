#include "incore_replicated_xc_device_integrator_integrate_den.hpp"
#include "incore_replicated_xc_device_integrator_exc_vxc.hpp"
#include "incore_replicated_xc_device_integrator_exc_grad.hpp"
#include "incore_replicated_xc_device_integrator_exx.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
IncoreReplicatedXCDeviceIntegrator<ValueType>::~IncoreReplicatedXCDeviceIntegrator() noexcept = default;



template class IncoreReplicatedXCDeviceIntegrator<double>;

}
}
