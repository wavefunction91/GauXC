/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/xc_integrator/pgas_dist/pgas_dist_xc_host_integrator.hpp>
#ifdef GAUXC_ENABLE_DEVICE
#include <gauxc/xc_integrator/pgas_dist/pgas_dist_xc_device_integrator.hpp>
#endif
#include <gauxc/xc_integrator/pgas_dist/impl.hpp>
#include <gauxc/exceptions.hpp>

namespace GauXC {

/// Factory to generate PGASDistributedXCIntegrator instances
template <typename MatrixType>
struct PGASDistributedXCIntegratorFactory {

  using integrator_type = detail::PGASDistributedXCIntegrator<MatrixType>;
  using value_type      = typename integrator_type::value_type;
  using ptr_return_t    = std::unique_ptr<integrator_type>;

  
  /** Generate a PGASDistributedXCIntegrator instance
   *
   *  @param[in]  ex                 Execution space for integrator instance
   *  @param[in]  integration_kernel Name of integration scaffold to load ("Default", "Reference", etc)
   *  @param[in]  func               XC functional to integrate
   *  @param[in]  lb                 Pregenerated LoadBalancer instance
   *  @param[in]  lwd                Local Work Driver
   */
  static ptr_return_t make_integrator_impl( 
    ExecutionSpace ex,
    std::string integrator_kernel,
    std::shared_ptr<functional_type>   func,
    std::shared_ptr<LoadBalancer>      lb,
    std::unique_ptr<LocalWorkDriver>&& lwd
    ) {



    switch(ex) {

      using host_factory = 
        detail::PGASDistributedXCHostIntegratorFactory<value_type>;
      case ExecutionSpace::Host:
        return std::make_unique<integrator_type>( 
          host_factory::make_integrator_impl(
            integrator_kernel, func, lb, std::move(lwd)
          )
        );

      #ifdef GAUXC_ENABLE_DEVICE
      using device_factory = 
        detail::PGASDistributedXCDeviceIntegratorFactory<value_type>;
      case ExecutionSpace::Device:
        return std::make_unique<integrator_type>( 
          device_factory::make_integrator_impl(
            integrator_kernel, func, lb, std::move(lwd)
          )
        );
      #endif

      default:
        GAUXC_GENERIC_EXCEPTION("PGASDistributedXCIntegrator ExecutionSpace Not Supported");
    }

    return nullptr;

  }

 
};


}
