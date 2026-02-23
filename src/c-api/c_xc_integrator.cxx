/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */

#include <gauxc/xc_integrator.h>
#include <gauxc/xc_integrator/replicated/replicated_xc_host_integrator.hpp>
#ifdef GAUXC_HAS_DEVICE
#include <gauxc/xc_integrator/replicated/replicated_xc_device_integrator.hpp>
#endif
#include <gauxc/xc_integrator/replicated/impl.hpp>
#include <gauxc/exceptions.hpp>

#include "c_load_balancer.hpp"
#include "c_functional.hpp"
#include "c_xc_integrator.hpp"
#include "c_status.hpp"

namespace GauXC::detail {

static inline std::unique_ptr<LocalWorkDriver>
get_local_work_driver(const char* local_work_kernel_name, ExecutionSpace ex) {
  return LocalWorkDriverFactory::make_local_work_driver( ex, 
    std::string(local_work_kernel_name), LocalWorkSettings() );
}

static inline std::shared_ptr<ReductionDriver>
get_reduction_driver(std::shared_ptr<LoadBalancer> lb, const char* reduction_kernel_name) {
  return ReductionDriverFactory::get_shared_instance( 
      lb->runtime(), std::string(reduction_kernel_name) );
}

static inline std::unique_ptr<ReplicatedXCIntegratorImpl<double>>
get_replicated_integrator( const char* integrator_kernel_name,
                          std::shared_ptr<functional_type> func,
                          std::shared_ptr<LoadBalancer> lb,
                          std::unique_ptr<LocalWorkDriver> lwd,
                          std::shared_ptr<ReductionDriver> rd,
                          ExecutionSpace ex ) {
  switch(ex) {
  case ExecutionSpace::Host:
    return ReplicatedXCHostIntegratorFactory<double>::make_integrator_impl(
          std::string(integrator_kernel_name), func, lb, std::move(lwd), rd );

  #ifdef GAUXC_HAS_DEVICE
  case ExecutionSpace::Device:
    return ReplicatedXCDeviceIntegratorFactory<double>::make_integrator_impl(
          std::string(integrator_kernel_name), func, lb, std::move(lwd), rd );

  #endif

  default:
    GAUXC_GENERIC_EXCEPTION("ReplicatedXCIntegrator ExecutionSpace Not Supported");
  }
  return nullptr;
}
}


namespace GauXC::C {
extern "C" {

void gauxc_integrator_delete(
  GauXCStatus* status,
  GauXCIntegrator* integrator
) {
  detail::gauxc_status_init(status);
  if (integrator == nullptr) return;
  if (integrator->ptr != nullptr) {
    delete detail::get_xc_integrator_ptr(*integrator);
  }
  integrator->ptr = nullptr;
}

GauXCIntegrator gauxc_integrator_new(
  GauXCStatus* status,
  const GauXCFunctional functional,
  const GauXCLoadBalancer lb,
  enum GauXC_ExecutionSpace execution_space,
  const char* integrator_input_type,
  const char* integrator_kernel_name,
  const char* local_work_kernel_name,
  const char* reduction_kernel_name
) {
  detail::gauxc_status_init(status);
  GauXCIntegrator integrator{};
  integrator.hdr = GauXCHeader{GauXC_Type_Integrator};
  integrator.ptr = nullptr;

  try {
    auto ex = static_cast<ExecutionSpace>(execution_space);
    auto func = std::make_shared<functional_type>(*detail::get_functional_ptr(functional));
    auto lb_ptr = *detail::get_load_balancer_ptr(lb);
    // Create Local Work Driver
    auto lwd = detail::get_local_work_driver(local_work_kernel_name, ex);
    // Create Reduction Driver
    auto rd = detail::get_reduction_driver(lb_ptr, reduction_kernel_name);

    // Create Integrator instance
    auto input_type_ = std::string(integrator_input_type);
    std::transform( input_type_.begin(), input_type_.end(), input_type_.begin(), ::toupper );

    if( input_type_ != "REPLICATED" ) GAUXC_GENERIC_EXCEPTION("INTEGRATOR TYPE NOT RECOGNIZED");

    integrator.ptr = 
      detail::get_replicated_integrator(integrator_kernel_name, func, lb_ptr, std::move(lwd), rd, ex).release();
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return integrator;
}


void gauxc_integrator_integrate_den(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix,
  int64_t ldp,
  double* den
) {
  detail::gauxc_status_init(status);
  try {
    detail::get_xc_integrator_ptr(integrator)->integrate_den(
      m, n,
      density_matrix, ldp,
      den );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

void gauxc_integrator_eval_exc_rks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix,
  int64_t ldp,
  double* exc
) {
  detail::gauxc_status_init(status);
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_exc(
      m, n,
      density_matrix, ldp,
      exc,
      IntegratorSettingsXC{} );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

void gauxc_integrator_eval_exc_uks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix_s,
  int64_t ldp_s,
  const double* density_matrix_z,
  int64_t ldp_z,
  double* exc
) {
  detail::gauxc_status_init(status);
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_exc(
      m, n,
      density_matrix_s, ldp_s,
      density_matrix_z, ldp_z,
      exc,
      IntegratorSettingsXC{} );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

void gauxc_integrator_eval_exc_gks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix_s,
  int64_t ldp_s,
  const double* density_matrix_z,
  int64_t ldp_z,
  const double* density_matrix_y,
  int64_t ldp_y,
  const double* density_matrix_x,
  int64_t ldp_x,
  double* exc
) {
  detail::gauxc_status_init(status);
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_exc(
      m, n,
      density_matrix_s, ldp_s,
      density_matrix_z, ldp_z,
      density_matrix_y, ldp_y,
      density_matrix_x, ldp_x,
      exc,
      IntegratorSettingsXC{} );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

void gauxc_integrator_eval_exc_vxc_rks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix,
  int64_t ldp,
  double* exc,
  double* vxc_matrix,
  int64_t vxc_ld
) {
  detail::gauxc_status_init(status);
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_exc_vxc(
      m, n,
      density_matrix, ldp,
      vxc_matrix, vxc_ld,
      exc,
      IntegratorSettingsXC{} );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

void gauxc_integrator_eval_exc_vxc_uks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix_s,
  int64_t ldp_s,
  const double* density_matrix_z,
  int64_t ldp_z,
  double* exc,
  double* vxc_matrix_s,
  int64_t vxc_ld_s,
  double* vxc_matrix_z,
  int64_t vxc_ld_z
) {
  detail::gauxc_status_init(status);
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_exc_vxc(
      m, n,
      density_matrix_s, ldp_s,
      density_matrix_z, ldp_z,
      vxc_matrix_s, vxc_ld_s,
      vxc_matrix_z, vxc_ld_z,
      exc,
      IntegratorSettingsXC{} );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

void gauxc_integrator_eval_exc_vxc_gks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix_s,
  int64_t ldp_s,
  const double* density_matrix_z,
  int64_t ldp_z,
  const double* density_matrix_y,
  int64_t ldp_y,
  const double* density_matrix_x,
  int64_t ldp_x,
  double* exc,
  double* vxc_matrix_s,
  int64_t vxc_ld_s,
  double* vxc_matrix_z,
  int64_t vxc_ld_z,
  double* vxc_matrix_y,
  int64_t vxc_ld_y,
  double* vxc_matrix_x,
  int64_t vxc_ld_x
) {
  detail::gauxc_status_init(status);
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_exc_vxc(
      m, n,
      density_matrix_s, ldp_s,
      density_matrix_z, ldp_z,
      density_matrix_y, ldp_y,
      density_matrix_x, ldp_x,
      vxc_matrix_s, vxc_ld_s,
      vxc_matrix_z, vxc_ld_z,
      vxc_matrix_y, vxc_ld_y,
      vxc_matrix_x, vxc_ld_x,
      exc,
      IntegratorSettingsXC{} );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

} // extern "C"
} // namespace GauXC::C
