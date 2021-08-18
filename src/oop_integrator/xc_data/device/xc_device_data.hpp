#pragma once
#include <gauxc/xc_task.hpp>
#include <vector>
#include <gauxc/basisset_map.hpp>
#include <gauxc/molmeta.hpp>

namespace GauXC {

// Base class for all XCDeviceData types
struct XCDeviceData {

  using host_task_type        = XCTask;
  using host_task_container   = std::vector<host_task_type>;
  using host_task_iterator    = host_task_container::iterator;

  virtual ~XCDeviceData() noexcept = default;

  virtual void allocate_static_data( int32_t natoms,
    int32_t nbf, int32_t nshells ) = 0;

  virtual void send_static_data( const double* P, int32_t ldp,
    const BasisSet<double>& basis, const Molecule& mol,
    const MolMeta& meta ) = 0;

  virtual void zero_integrands() = 0;

  virtual host_task_iterator generate_buffers( 
    const BasisSetMap& basis_map, host_task_iterator task_begin,
    host_task_iterator task_end ) = 0;

  virtual void retrieve_xc_integrands( double* EXC, double* N_EL,
    double* VXC, int32_t ldvxc ) = 0;

};

}
