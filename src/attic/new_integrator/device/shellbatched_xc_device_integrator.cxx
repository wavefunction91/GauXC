#include <gauxc/new_xc_integrator/replicated/shellbatched_xc_device_integrator.hpp>

#include "device/shellbatched_xc_device_exc_vxc.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
ShellBatchedXCDeviceIntegrator<ValueType>::
  ShellBatchedXCDeviceIntegrator( const ShellBatchedXCDeviceIntegrator& ) = default;

template <typename ValueType>
ShellBatchedXCDeviceIntegrator<ValueType>::
  ShellBatchedXCDeviceIntegrator( ShellBatchedXCDeviceIntegrator&& ) noexcept = default;

template <typename ValueType>
ShellBatchedXCDeviceIntegrator<ValueType>::
  ~ShellBatchedXCDeviceIntegrator() noexcept = default;





template class ShellBatchedXCDeviceIntegrator<double>;

}
}
