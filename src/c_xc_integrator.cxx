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
#include <gauxc/xc_integrator.hpp>
#include <gauxc/util/c_load_balancer.hpp>
#include <gauxc/util/c_functional.hpp>
#include <gauxc/util/c_xc_integrator.hpp>
#include <gauxc/util/c_status.hpp>

namespace GauXC::C {
extern "C" {

void gauxc_integrator_delete(
  GauXCStatus* status,
  GauXCIntegrator* integrator
) {
  status->code = 0;
  if (integrator == nullptr) return;
  if (integrator->ptr != nullptr) {
    if (integrator->owned)
      delete detail::get_xc_integrator_ptr(*integrator);
    else
      delete detail::get_xc_integrator_shared(*integrator);
  }
  integrator->ptr = nullptr;
}

void gauxc_integrator_factory_delete(
  GauXCStatus* status,
  GauXCIntegratorFactory* factory
) {
  status->code = 0;
  if (factory == nullptr) return;
  if (factory->ptr != nullptr)
    delete detail::get_xc_integrator_factory_ptr(*factory);
  factory->ptr = nullptr;
}

GauXCIntegratorFactory gauxc_integrator_factory_new(
  GauXCStatus* status,
  enum GauXC_ExecutionSpace execution_space,
  const char* integrator_input_type,
  const char* integrator_kernel_name,
  const char* local_work_kernel_name,
  const char* reduction_kernel_name
) {
  GauXCIntegratorFactory factory{};
  factory.hdr = GauXCHeader{GauXC_Type_IntegratorFactory};
  factory.ptr = nullptr;

  try {
    factory.ptr = new XCIntegratorFactory<detail::CMatrix>(
      static_cast<ExecutionSpace>(execution_space),
      std::string(integrator_input_type),
      std::string(integrator_kernel_name),
      std::string(local_work_kernel_name),
      std::string(reduction_kernel_name)
    );
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return factory;
}


GauXCIntegrator gauxc_integrator_factory_get_instance(
  GauXCStatus* status,
  const GauXCIntegratorFactory factory,
  const GauXCFunctional functional,
  const GauXCLoadBalancer lb
) {
  GauXCIntegrator integrator{};
  integrator.hdr = GauXCHeader{GauXC_Type_Integrator};
  integrator.ptr = nullptr;
  integrator.owned = true;

  try {
    auto integrator_instance = detail::get_integrator_instance(factory, functional, lb);
    integrator.ptr = new XCIntegrator<detail::CMatrix>(std::move(integrator_instance));
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return integrator;
}

GauXCIntegrator gauxc_integrator_factory_get_shared_instance(
  GauXCStatus* status,
  const GauXCIntegratorFactory factory,
  const GauXCFunctional functional,
  const GauXCLoadBalancer lb
) {
  GauXCIntegrator integrator{};
  integrator.hdr = GauXCHeader{GauXC_Type_Integrator};
  integrator.ptr = nullptr;
  integrator.owned = false;

  try {
    auto integrator_instance = detail::get_shared_integrator_instance(factory, functional, lb);
    integrator.ptr = new std::shared_ptr<XCIntegrator<detail::CMatrix>>(std::move(integrator_instance));
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return integrator;
}

void gauxc_integrator_integrate_den(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const GauXCMatrix density_matrix,
  double* den_out
) {
  try {
    auto& dm = *detail::get_matrix_ptr(density_matrix);

    auto den = integrator.owned
      ? detail::get_xc_integrator_ptr(integrator)->integrate_den(dm)
      : detail::get_xc_integrator_shared(integrator)->get()->integrate_den(dm);

    *den_out = den;
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
}

void gauxc_integrator_eval_exc_rks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const GauXCMatrix density_matrix,
  double* exc_out
) {
  try {
    auto& dm = *detail::get_matrix_ptr(density_matrix);

    auto exc = integrator.owned
      ? detail::get_xc_integrator_ptr(integrator)->eval_exc(dm)
      : detail::get_xc_integrator_shared(integrator)->get()->eval_exc(dm);

    *exc_out = exc;
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
}

void gauxc_integrator_eval_exc_uks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const GauXCMatrix density_matrix_s,
  const GauXCMatrix density_matrix_z,
  double* exc_out
) {
  try {
    auto& dm_s = *detail::get_matrix_ptr(density_matrix_s);
    auto& dm_z = *detail::get_matrix_ptr(density_matrix_z);

    auto exc = integrator.owned
      ? detail::get_xc_integrator_ptr(integrator)->eval_exc(dm_s, dm_z)
      : detail::get_xc_integrator_shared(integrator)->get()->eval_exc(dm_s, dm_z);

    *exc_out = exc;
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
}

void gauxc_integrator_eval_exc_gks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const GauXCMatrix density_matrix_s,
  const GauXCMatrix density_matrix_z,
  const GauXCMatrix density_matrix_x,
  const GauXCMatrix density_matrix_y,
  double* exc_out
) {
  try {
    auto& dm_s = *detail::get_matrix_ptr(density_matrix_s);
    auto& dm_z = *detail::get_matrix_ptr(density_matrix_z);
    auto& dm_x = *detail::get_matrix_ptr(density_matrix_x);
    auto& dm_y = *detail::get_matrix_ptr(density_matrix_y);

    auto exc = integrator.owned
      ? detail::get_xc_integrator_ptr(integrator)->eval_exc(dm_s, dm_z, dm_x, dm_y)
      : detail::get_xc_integrator_shared(integrator)->get()->eval_exc(dm_s, dm_z, dm_x, dm_y);

    *exc_out = exc;
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
}

void gauxc_integrator_eval_exc_vxc_rks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const GauXCMatrix density_matrix,
  double* exc_out,
  GauXCMatrix* vxc_matrix
) {
  vxc_matrix->ptr = nullptr;

  try {
    auto& dm = *detail::get_matrix_ptr(density_matrix);

    auto [exc, vxc] = integrator.owned
      ? detail::get_xc_integrator_ptr(integrator)->eval_exc_vxc(dm)
      : detail::get_xc_integrator_shared(integrator)->get()->eval_exc_vxc(dm);

    *exc_out = exc;
    vxc_matrix->ptr = new detail::CMatrix(std::move(vxc));
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
}

void gauxc_integrator_eval_exc_vxc_uks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const GauXCMatrix density_matrix_s,
  const GauXCMatrix density_matrix_z,
  double* exc_out,
  GauXCMatrix* vxc_matrix_s,
  GauXCMatrix* vxc_matrix_z
) {
  vxc_matrix_s->ptr = nullptr;
  vxc_matrix_z->ptr = nullptr;

  try {
    auto& dm_s = *detail::get_matrix_ptr(density_matrix_s);
    auto& dm_z = *detail::get_matrix_ptr(density_matrix_z);

    auto [exc, vxc_s, vxc_z] = integrator.owned
      ? detail::get_xc_integrator_ptr(integrator)->eval_exc_vxc(dm_s, dm_z)
      : detail::get_xc_integrator_shared(integrator)->get()->eval_exc_vxc(dm_s, dm_z);

    *exc_out = exc;
    vxc_matrix_s->ptr = new detail::CMatrix(std::move(vxc_s));
    vxc_matrix_z->ptr = new detail::CMatrix(std::move(vxc_z));
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
}

void gauxc_integrator_eval_exc_vxc_onedft_uks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const GauXCMatrix density_matrix_s,
  const GauXCMatrix density_matrix_z,
  const char* model,
  double* exc_out,
  GauXCMatrix* vxc_matrix_s,
  GauXCMatrix* vxc_matrix_z
) {
  vxc_matrix_s->ptr = nullptr;
  vxc_matrix_z->ptr = nullptr;
  OneDFTSettings onedft_settings;
  onedft_settings.model = std::string(model);

  try {
    auto& dm_s = *detail::get_matrix_ptr(density_matrix_s);
    auto& dm_z = *detail::get_matrix_ptr(density_matrix_z);

    auto [exc, vxc_s, vxc_z] = integrator.owned
      ? detail::get_xc_integrator_ptr(integrator)->eval_exc_vxc_onedft(dm_s, dm_z, onedft_settings)
      : detail::get_xc_integrator_shared(integrator)->get()->eval_exc_vxc_onedft(dm_s, dm_z, onedft_settings);

    *exc_out = exc;
    vxc_matrix_s->ptr = new detail::CMatrix(std::move(vxc_s));
    vxc_matrix_z->ptr = new detail::CMatrix(std::move(vxc_z));
    status->code = 0;
  } catch (...) {
    status->code = 1;
  }
}

void gauxc_integrator_eval_exc_vxc_gks(
  GauXCStatus* status,
  const GauXCIntegrator integrator,
  const GauXCMatrix density_matrix_s,
  const GauXCMatrix density_matrix_z,
  const GauXCMatrix density_matrix_x,
  const GauXCMatrix density_matrix_y,
  double* exc_out,
  GauXCMatrix* vxc_matrix_s,
  GauXCMatrix* vxc_matrix_z,
  GauXCMatrix* vxc_matrix_x,
  GauXCMatrix* vxc_matrix_y
) {
  vxc_matrix_s->ptr = nullptr;
  vxc_matrix_z->ptr = nullptr;
  vxc_matrix_x->ptr = nullptr;
  vxc_matrix_y->ptr = nullptr;

  try {
    auto& dm_s = *detail::get_matrix_ptr(density_matrix_s);
    auto& dm_z = *detail::get_matrix_ptr(density_matrix_z);
    auto& dm_x = *detail::get_matrix_ptr(density_matrix_x);
    auto& dm_y = *detail::get_matrix_ptr(density_matrix_y);

    auto [exc, vxc_s, vxc_z, vxc_x, vxc_y] = integrator.owned
      ? detail::get_xc_integrator_ptr(integrator)->eval_exc_vxc(dm_s, dm_z, dm_x, dm_y)
      : detail::get_xc_integrator_shared(integrator)->get()->eval_exc_vxc(dm_s, dm_z, dm_x, dm_y);

    *exc_out = exc;
    vxc_matrix_s->ptr = new detail::CMatrix(std::move(vxc_s));
    vxc_matrix_z->ptr = new detail::CMatrix(std::move(vxc_z));
    vxc_matrix_x->ptr = new detail::CMatrix(std::move(vxc_x));
    vxc_matrix_y->ptr = new detail::CMatrix(std::move(vxc_y));
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
}

} // extern "C"
} // namespace GauXC::C
