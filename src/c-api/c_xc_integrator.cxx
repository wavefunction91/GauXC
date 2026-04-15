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

#include <algorithm>
#include <cctype>

#include <gauxc/c/xc_integrator.h>

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
  if (functional.ptr == nullptr || functional.hdr.type != GauXC_Type_Functional) {
    detail::gauxc_status_handle(status, 1, "Invalid Functional handle");
    return integrator;
  }
  if (lb.ptr == nullptr || lb.hdr.type != GauXC_Type_LoadBalancer) {
    detail::gauxc_status_handle(status, 1, "Invalid LoadBalancer handle");
    return integrator;
  }
  if (integrator_input_type == nullptr) {
    detail::gauxc_status_handle(status, 1, "Integrator input type string cannot be null");
    return integrator;
  }
  if (integrator_kernel_name == nullptr) {
    detail::gauxc_status_handle(status, 1, "Integrator kernel name string cannot be null");
    return integrator;
  }
  if (local_work_kernel_name == nullptr) {
    detail::gauxc_status_handle(status, 1, "Local work kernel name string cannot be null");
    return integrator;
  }
  if (reduction_kernel_name == nullptr) {
    detail::gauxc_status_handle(status, 1, "Reduction kernel name string cannot be null");
    return integrator;
  }

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
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix pointer cannot be null");
    return;
  }
  if (den == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density output pointer cannot be null");
    return;
  }
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
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix pointer cannot be null");
    return;
  }
  if (exc == nullptr) {
    detail::gauxc_status_handle(status, 1, "Exc output pointer cannot be null");
    return;
  }
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
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix S pointer cannot be null");
    return;
  }
  if (density_matrix_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix Z pointer cannot be null");
    return;
  }
  if (exc == nullptr) {
    detail::gauxc_status_handle(status, 1, "Exc output pointer cannot be null");
    return;
  }
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
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix S pointer cannot be null");
    return;
  }
  if (density_matrix_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix Z pointer cannot be null");
    return;
  }
  if (density_matrix_y == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix Y pointer cannot be null");
    return;
  }
  if (density_matrix_x == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix X pointer cannot be null");
    return;
  }
  if (exc == nullptr) {
    detail::gauxc_status_handle(status, 1, "Exc output pointer cannot be null");
    return;
  }
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
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix pointer cannot be null");
    return;
  }
  if (exc == nullptr) {
    detail::gauxc_status_handle(status, 1, "Exc output pointer cannot be null");
    return;
  }
  if (vxc_matrix == nullptr) {
    detail::gauxc_status_handle(status, 1, "VXC matrix pointer cannot be null");
    return;
  }
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
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix S pointer cannot be null");
    return;
  }
  if (density_matrix_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix Z pointer cannot be null");
    return;
  }
  if (exc == nullptr) {
    detail::gauxc_status_handle(status, 1, "Exc output pointer cannot be null");
    return;
  }
  if (vxc_matrix_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "VXC matrix S pointer cannot be null");
    return;
  }
  if (vxc_matrix_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "VXC matrix Z pointer cannot be null");
    return;
  }
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

void gauxc_integrator_eval_exc_vxc_onedft_uks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix_s,
  int64_t ldp_s,
  const double* density_matrix_z,
  int64_t ldp_z,
  const char* model,
  double* exc,
  double* vxc_matrix_s,
  int64_t vxc_ld_s,
  double* vxc_matrix_z,
  int64_t vxc_ld_z
) {
  detail::gauxc_status_init(status);
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix S pointer cannot be null");
    return;
  }
  if (density_matrix_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix Z pointer cannot be null");
    return;
  }
  if (model == nullptr) {
    detail::gauxc_status_handle(status, 1, "Model string cannot be null");
    return;
  }
  if (exc == nullptr) {
    detail::gauxc_status_handle(status, 1, "Exc output pointer cannot be null");
    return;
  }
  if (vxc_matrix_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "VXC matrix S pointer cannot be null");
    return;
  }
  if (vxc_matrix_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "VXC matrix Z pointer cannot be null");
    return;
  }
  OneDFTSettings settings{};
  settings.model = std::string(model);
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_exc_vxc_onedft(
      m, n,
      density_matrix_s, ldp_s,
      density_matrix_z, ldp_z,
      vxc_matrix_s, vxc_ld_s,
      vxc_matrix_z, vxc_ld_z,
      exc,
      settings );
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
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix S pointer cannot be null");
    return;
  }
  if (density_matrix_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix Z pointer cannot be null");
    return;
  }
  if (density_matrix_y == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix Y pointer cannot be null");
    return;
  }
  if (density_matrix_x == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix X pointer cannot be null");
    return;
  }
  if (exc == nullptr) {
    detail::gauxc_status_handle(status, 1, "Exc output pointer cannot be null");
    return;
  }
  if (vxc_matrix_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "VXC matrix S pointer cannot be null");
    return;
  }
  if (vxc_matrix_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "VXC matrix Z pointer cannot be null");
    return;
  }
  if (vxc_matrix_y == nullptr) {
    detail::gauxc_status_handle(status, 1, "VXC matrix Y pointer cannot be null");
    return;
  }
  if (vxc_matrix_x == nullptr) {
    detail::gauxc_status_handle(status, 1, "VXC matrix X pointer cannot be null");
    return;
  }
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

void gauxc_integrator_eval_exc_grad_rks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix,
  int64_t ldp,
  double* exc_grad
) {
  detail::gauxc_status_init(status);
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix pointer cannot be null");
    return;
  }
  if (exc_grad == nullptr) {
    detail::gauxc_status_handle(status, 1, "Exc gradient output pointer cannot be null");
    return;
  }
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_exc_grad(
      m, n,
      density_matrix, ldp,
      exc_grad,
      IntegratorSettingsXC{} );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

void gauxc_integrator_eval_exc_grad_uks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix_s,
  int64_t ldp_s,
  const double* density_matrix_z,
  int64_t ldp_z,
  double* exc_grad
) {
  detail::gauxc_status_init(status);
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix S pointer cannot be null");
    return;
  }
  if (density_matrix_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix Z pointer cannot be null");
    return;
  }
  if (exc_grad == nullptr) {
    detail::gauxc_status_handle(status, 1, "Exc gradient output pointer cannot be null");
    return;
  }
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_exc_grad(
      m, n,
      density_matrix_s, ldp_s,
      density_matrix_z, ldp_z,
      exc_grad,
      IntegratorSettingsXC{} );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

void gauxc_integrator_eval_exc_grad_onedft_uks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix_s,
  int64_t ldp_s,
  const double* density_matrix_z,
  int64_t ldp_z,
  const char* model,
  double* exc_grad
) {
  detail::gauxc_status_init(status);
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix S pointer cannot be null");
    return;
  }
  if (density_matrix_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix Z pointer cannot be null");
    return;
  }
  if (model == nullptr) {
    detail::gauxc_status_handle(status, 1, "Model string cannot be null");
    return;
  }
  if (exc_grad == nullptr) {
    detail::gauxc_status_handle(status, 1, "Exc gradient output pointer cannot be null");
    return;
  }
  OneDFTSettings settings{};
  settings.model = std::string(model);
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_exc_grad_onedft(
      m, n,
      density_matrix_s, ldp_s,
      density_matrix_z, ldp_z,
      exc_grad,
      settings );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

void gauxc_integrator_eval_exx_rks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix,
  int64_t ldp,
  double* K,
  int64_t ldk
) {
  detail::gauxc_status_init(status);
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix pointer cannot be null");
    return;
  }
  if (K == nullptr) {
    detail::gauxc_status_handle(status, 1, "K output pointer cannot be null");
    return;
  }
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_exx(
      m, n,
      density_matrix, ldp,
      K, ldk,
      IntegratorSettingsEXX{} );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

void gauxc_integrator_eval_fxc_contraction_rks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix,
  int64_t ldp,
  const double* t_density_matrix,
  int64_t ldtp,
  double* fxc,
  int64_t ldfxc
) {
  detail::gauxc_status_init(status);
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix pointer cannot be null");
    return;
  }
  if (t_density_matrix == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix pointer cannot be null");
    return;
  }
  if (fxc == nullptr) {
    detail::gauxc_status_handle(status, 1, "FXC output pointer cannot be null");
    return;
  }
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_fxc_contraction(
      m, n,
      density_matrix, ldp,
      t_density_matrix, ldtp,
      fxc, ldfxc,
      IntegratorSettingsXC{} );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

void gauxc_integrator_eval_fxc_contraction_uks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  int64_t m,
  int64_t n,
  const double* density_matrix_s,
  int64_t ldp_s,
  const double* density_matrix_z,
  int64_t ldp_z,
  const double* t_density_matrix_s,
  int64_t ldtp_s,
  const double* t_density_matrix_z,
  int64_t ldtp_z,
  double* fxc_s,
  int64_t ldfxc_s,
  double* fxc_z,
  int64_t ldfxc_z
) {
  detail::gauxc_status_init(status);
  if (integrator.ptr == nullptr || integrator.hdr.type != GauXC_Type_Integrator) {
    detail::gauxc_status_handle(status, 1, "Invalid Integrator handle");
    return;
  }
  if (density_matrix_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix S pointer cannot be null");
    return;
  }
  if (density_matrix_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix Z pointer cannot be null");
    return;
  }
  if (t_density_matrix_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix S pointer cannot be null");
    return;
  }
  if (t_density_matrix_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "Density matrix Z pointer cannot be null");
    return;
  }
  if (fxc_s == nullptr) {
    detail::gauxc_status_handle(status, 1, "FXC S output pointer cannot be null");
    return;
  }
  if (fxc_z == nullptr) {
    detail::gauxc_status_handle(status, 1, "FXC Z output pointer cannot be null");
    return;
  }
  try {
    detail::get_xc_integrator_ptr(integrator)->eval_fxc_contraction(
      m, n,
      density_matrix_s, ldp_s,
      density_matrix_z, ldp_z,
      t_density_matrix_s, ldtp_s,
      t_density_matrix_z, ldtp_z,
      fxc_s, ldfxc_s,
      fxc_z, ldfxc_z,
      IntegratorSettingsXC{} );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

} // extern "C"
} // namespace GauXC::C
