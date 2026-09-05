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

#include <memory>
#include <vector>

#include <gauxc/basisset.hpp>
#include <gauxc/basisset_map.hpp>
#include <gauxc/orbital_evaluator.hpp>
#include <gauxc/xc_integrator/local_work_driver.hpp>

#include "xc_integrator/local_work_driver/host/local_host_work_driver.hpp"

namespace GauXC::detail {

/// Host implementation state for OrbitalEvaluator
class OrbitalEvaluatorImpl {

public:

  BasisSet<double> basis;                     ///< Basis set to evaluate
  std::unique_ptr<LocalWorkDriver> driver_owner;
  LocalHostWorkDriver* host_driver = nullptr; ///< Non-owning view of the driver
  std::unique_ptr<BasisSetMap> basis_map;     ///< shell -> AO range lookups
  std::vector<double> shell_cutoff_r2;        ///< Per-shell squared cutoff radius
  int32_t nbf_ = 0;

  OrbitalEvaluatorImpl( BasisSet<double> bs, double screening_tolerance );

}; // class OrbitalEvaluatorImpl

}  // namespace GauXC::detail
