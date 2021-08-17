#pragma once
#include <gauxc/xc_task.hpp>
#include <vector>
#include <gauxc/basisset_map.hpp>

namespace GauXC {

// Base class for all XCDeviceData types
struct XCDeviceData {

  using host_task_type        = XCTask;
  using host_task_container   = std::vector<host_task_type>;
  using host_task_iterator    = host_task_container::iterator;

  virtual ~XCDeviceData() noexcept = default;

  virtual void allocate_static_data( int32_t natoms,
    int32_t nbf, int32_t nshells ) = 0;

  virtual host_task_iterator generate_buffers( 
    const BasisSetMap& basis_map, host_task_iterator task_begin,
    host_task_iterator task_end ) = 0;

};

}
