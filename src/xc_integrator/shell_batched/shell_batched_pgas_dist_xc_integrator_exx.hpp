/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "shell_batched_pgas_dist_xc_integrator.hpp"
#include <gauxc/util/misc.hpp>
#include <gauxc/util/unused.hpp>

namespace GauXC  {
namespace detail {

template <typename BaseIntegratorType, typename IncoreIntegratorType>
void ShellBatchedPGASDistributedXCIntegrator<BaseIntegratorType, IncoreIntegratorType>::
  eval_exx_( const matrix_type& P, matrix_type& K, 
            const IntegratorSettingsEXX& settings ) {
  GAUXC_GENERIC_EXCEPTION("NYI" );                 
  util::unused(P,K,settings);
}

}
}
