#include "shellbatched_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"
#include "device/xc_device_aos_data.hpp"
#include "integrator_util/integrator_common.hpp"
#include "host/util.hpp"
#include <gauxc/util/misc.hpp>

#include <stdexcept>
#include <fstream>
#include <queue>
#include <mutex>
#include <future>
#include <set>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void ShellBatchedReplicatedXCDeviceIntegrator<ValueType>::
  eval_exx_( int64_t m, int64_t n, const value_type* P,
             int64_t ldp, value_type* K, int64_t ldk, 
             const IntegratorSettingsEXX& settings ) { 
             
  throw std::runtime_error(std::string(__PRETTY_FUNCTION__) + " NYI" );                 
}

}
}
