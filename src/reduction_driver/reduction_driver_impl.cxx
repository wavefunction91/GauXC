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
#include "reduction_driver_impl.hpp"
#include <gauxc/exceptions.hpp>

namespace GauXC::detail {


ReductionDriverImpl::ReductionDriverImpl( const RuntimeEnvironment& rt ) 
  : runtime_(rt){}

ReductionDriverImpl::~ReductionDriverImpl() noexcept = default;
ReductionDriverImpl::ReductionDriverImpl(const ReductionDriverImpl& ) = default;

void ReductionDriverImpl::allgather_v_typeerased( const void*, size_t,
  std::vector<std::byte>&, std::type_index, std::any ) {
  GAUXC_GENERIC_EXCEPTION("Variable-size allgather is not supported by this ReductionDriver");
}

int ReductionDriverImpl::comm_rank() const { return runtime_.comm_rank(); }
int ReductionDriverImpl::comm_size() const { return runtime_.comm_size(); }

}
