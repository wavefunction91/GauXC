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

#include <gauxc/load_balancer.h>
#include <gauxc/load_balancer.hpp>

namespace GauXC::detail {
static inline std::shared_ptr<LoadBalancer>* get_load_balancer_ptr(C::GauXCLoadBalancer lb) noexcept {
  return static_cast<std::shared_ptr<LoadBalancer>*>(lb.ptr);
}
static inline LoadBalancerFactory* get_load_balancer_factory_ptr(C::GauXCLoadBalancerFactory lbf) noexcept {
  return static_cast<LoadBalancerFactory*>(lbf.ptr);
}
}