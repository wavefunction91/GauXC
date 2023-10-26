/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include <gauxc/xc_integrator/pgas_dist/pgas_dist_xc_host_integrator.hpp>
#include "shell_batched_pgas_dist_xc_host_integrator.hpp"
#include "host/local_host_work_driver.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
PGASDistributedXCHostIntegrator<ValueType>::~PGASDistributedXCHostIntegrator() noexcept = default;

template class PGASDistributedXCHostIntegrator<double>;


template <typename ValueType>
typename PGASDistributedXCHostIntegratorFactory<ValueType>::ptr_return_t
  PGASDistributedXCHostIntegratorFactory<ValueType>::make_integrator_impl(
    std::string integrator_kernel,
    std::shared_ptr<functional_type> func,
    std::shared_ptr<LoadBalancer> lb, 
    std::unique_ptr<LocalWorkDriver>&& lwd
    ) {

  // Make sure that the LWD is a valid LocalHostWorkDriver
  if(not dynamic_cast<LocalHostWorkDriver*>(lwd.get())) {
    GAUXC_GENERIC_EXCEPTION("Passed LWD Not valid for Host ExSpace");
  }

  std::transform(integrator_kernel.begin(), integrator_kernel.end(), 
    integrator_kernel.begin(), ::toupper );

  if( integrator_kernel == "DEFAULT" ) integrator_kernel = "SHELLBATCHED";


  if( integrator_kernel == "SHELLBATCHED" )
    return std::make_unique<ShellBatchedPGASDistributedXCHostIntegrator<ValueType>>(
      func, lb, std::move(lwd)
    );

  else
    GAUXC_GENERIC_EXCEPTION("INTEGRATOR KERNEL: " + integrator_kernel + " NOT RECOGNIZED");

  return nullptr;


}

template class PGASDistributedXCHostIntegratorFactory<double>;


}
}

