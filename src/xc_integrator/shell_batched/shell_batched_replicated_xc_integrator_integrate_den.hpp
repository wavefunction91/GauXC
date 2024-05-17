/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "shell_batched_replicated_xc_integrator.hpp"
#include <gauxc/util/misc.hpp>
#include <gauxc/util/unused.hpp>

namespace GauXC  {
namespace detail {

template <typename BaseIntegratorType, typename IncoreIntegratorType>
void ShellBatchedReplicatedXCIntegrator<BaseIntegratorType, IncoreIntegratorType>::
  integrate_den_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* N_EL ) { 
                 
  GAUXC_GENERIC_EXCEPTION("NYI" );                 
  util::unused(m,n,P,ldp,N_EL);
}

}
}
