#include "shellbatched_replicated_xc_device_integrator_integrate_den.hpp"
#include "shellbatched_replicated_xc_device_integrator_exc_vxc.hpp"
#include "shellbatched_replicated_xc_device_integrator_exc_grad.hpp"
#include "shellbatched_replicated_xc_device_integrator_exx.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
ShellBatchedReplicatedXCDeviceIntegrator<ValueType>::~ShellBatchedReplicatedXCDeviceIntegrator() noexcept = default;

template class ShellBatchedReplicatedXCDeviceIntegrator<double>;

}
}
