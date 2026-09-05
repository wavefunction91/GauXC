#include <gauxc/new_xc_integrator/replicated/reference_xc_host_integrator.hpp>

#include "host/reference_xc_host_exc_vxc.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
ReferenceXCHostIntegrator<ValueType>::
  ReferenceXCHostIntegrator( const ReferenceXCHostIntegrator& ) = default;

template <typename ValueType>
ReferenceXCHostIntegrator<ValueType>::
  ReferenceXCHostIntegrator( ReferenceXCHostIntegrator&& ) noexcept = default;

template <typename ValueType>
ReferenceXCHostIntegrator<ValueType>::
  ~ReferenceXCHostIntegrator() noexcept = default;





template class ReferenceXCHostIntegrator<double>;

}
}
