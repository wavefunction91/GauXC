#include "incore_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"
#include <stdexcept>
#include "device/xc_device_aos_data.hpp"
#include <fstream>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_grad_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* EXC_GRAD ) { 
                 
  throw std::runtime_error(std::string(__PRETTY_FUNCTION__) + " NYI" );                 
}

}
}
