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
#pragma once

#include <gauxc/xc_integrator.h>
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>
#include <gauxc/util/c_matrix.hpp>
#include <gauxc/util/c_load_balancer.hpp>
#include <gauxc/util/c_functional.hpp>

namespace GauXC::detail {

static inline XCIntegrator<CMatrix>* get_xc_integrator_ptr(C::GauXCIntegrator integrator) noexcept {
  return static_cast<XCIntegrator<CMatrix>*>(integrator.ptr);
}
static inline std::shared_ptr<XCIntegrator<CMatrix>>* get_xc_integrator_shared(C::GauXCIntegrator integrator) noexcept {
  return static_cast<std::shared_ptr<XCIntegrator<CMatrix>>*>(integrator.ptr);
}
static inline XCIntegratorFactory<CMatrix>* get_xc_integrator_factory_ptr(C::GauXCIntegratorFactory factory) noexcept {
  return static_cast<XCIntegratorFactory<CMatrix>*>(factory.ptr);
}
static inline XCIntegrator<CMatrix> get_integrator_instance(
  C::GauXCIntegratorFactory factory,
  C::GauXCFunctional functional,
  C::GauXCLoadBalancer lb
) {
  if (lb.owned)
    return get_xc_integrator_factory_ptr(factory)->get_instance(
      *get_functional_ptr(functional),
      *get_load_balancer_ptr(lb)
    );
  else
    return get_xc_integrator_factory_ptr(factory)->get_instance(
      *get_functional_ptr(functional),
      **get_load_balancer_shared(lb)
    );
}
static inline std::shared_ptr<XCIntegrator<CMatrix>> get_shared_integrator_instance(
  C::GauXCIntegratorFactory factory,
  C::GauXCFunctional functional,
  C::GauXCLoadBalancer lb
) {
  if (lb.owned)
    return get_xc_integrator_factory_ptr(factory)->get_shared_instance(
      *get_functional_ptr(functional),
      *get_load_balancer_ptr(lb)
    );
  else
    return get_xc_integrator_factory_ptr(factory)->get_shared_instance(
      *get_functional_ptr(functional),
      **get_load_balancer_shared(lb)
    );
}

} // namespace GauXC::detail