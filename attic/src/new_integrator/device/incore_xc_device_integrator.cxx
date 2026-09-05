#include <gauxc/new_xc_integrator/replicated/incore_xc_device_integrator.hpp>

#include "device/incore_xc_device_exc_vxc.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
IncoreXCDeviceIntegrator<ValueType>::
  IncoreXCDeviceIntegrator( const IncoreXCDeviceIntegrator& ) = default;

template <typename ValueType>
IncoreXCDeviceIntegrator<ValueType>::
  IncoreXCDeviceIntegrator( IncoreXCDeviceIntegrator&& ) noexcept = default;

template <typename ValueType>
IncoreXCDeviceIntegrator<ValueType>::
  ~IncoreXCDeviceIntegrator() noexcept = default;





template class IncoreXCDeviceIntegrator<double>;

}
}
