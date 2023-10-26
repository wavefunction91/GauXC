/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include <gauxc/xc_integrator/pgas_dist/pgas_dist_xc_device_integrator.hpp>
#include "shellbatched_pgas_dist_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
PGASDistributedXCDeviceIntegrator<ValueType>::~PGASDistributedXCDeviceIntegrator() noexcept = default;

template class PGASDistributedXCDeviceIntegrator<double>;


template <typename ValueType>
typename PGASDistributedXCDeviceIntegratorFactory<ValueType>::ptr_return_t
  PGASDistributedXCDeviceIntegratorFactory<ValueType>::make_integrator_impl(
    std::string integrator_kernel,
    std::shared_ptr<functional_type> func,
    std::shared_ptr<LoadBalancer> lb, 
    std::unique_ptr<LocalWorkDriver>&& lwd
    ) {

  // Make sure that the LWD is a valid LocalDeviceWorkDriver
  if(not dynamic_cast<LocalDeviceWorkDriver*>(lwd.get())) {
    GAUXC_GENERIC_EXCEPTION("Passed LWD Not valid for Device ExSpace");
  }

  std::transform(integrator_kernel.begin(), integrator_kernel.end(), 
    integrator_kernel.begin(), ::toupper );

  if( integrator_kernel == "DEFAULT" ) integrator_kernel = "SHELLBATCHED";

  if( integrator_kernel == "SHELLBATCHED" )
    return std::make_unique<ShellBatchedPGASDistributedXCDeviceIntegrator<ValueType>>(
      func, lb, std::move(lwd)
    );

  else
    GAUXC_GENERIC_EXCEPTION("INTEGRATOR KERNEL " + integrator_kernel + " NOT RECOGNIZED");


}

template struct PGASDistributedXCDeviceIntegratorFactory<double>;


}
}

