#pragma once
#include <gauxc/xc_task.hpp>
#include <vector>
#include <gauxc/basisset_map.hpp>
#include <gauxc/molmeta.hpp>
//#include <gauxc/reduction_driver.hpp>
#include <any>
#include <cstring>
#include "type_erased_queue.hpp"

namespace GauXC {

struct integrator_term_tracker {
  bool weights = false;
  bool exc_vxc = false;
  inline void reset() {
    std::memset( this, 0, sizeof(integrator_term_tracker) );
  }
};

/** Base class for all XCDeviceData types
 *
 *  Exposes virtual API to manage device memory and batch XC
 *  integration tasks.
 */
struct XCDeviceData {

  using host_task_type        = XCTask;
  using host_task_container   = std::vector<host_task_type>;
  using host_task_iterator    = host_task_container::iterator;

  virtual ~XCDeviceData() noexcept = default;

  /// Allocate device memory for data that will persist on the device.
  virtual void reset_allocations() = 0;
  virtual void allocate_static_data_weights( int32_t natoms ) = 0;
  virtual void allocate_static_data_exc_vxc( int32_t nbf, int32_t nshells ) = 0;

  // Send persistent data from host to device
  virtual void send_static_data_weights( const Molecule& mol, const MolMeta& meta ) = 0;
  virtual void send_static_data_exc_vxc( const double* P, int32_t ldp, const BasisSet<double>& basis ) = 0;

  /** Zero out the integrands in device memory
   *
   *  TODO: this will depend on the integrand, we should refactor this
   *  to only allocate what is needed
   */
  virtual void zero_integrands() = 0;

  /** Generate task batch to execute on device
   *
   *  Generate a batch of XC tasks to execute on the device and 
   *  populate device memory for said batch.
   *
   *  TODO: this will depend on the integrand, we should refactor this
   *  to only allocate what is needed
   *
   *  @param[in] basis_map  Basis set map instance for passed basis object
   *                        (TODO, this should probably persist to avoid clashes)
   *  @param[in] task_begin Start iterator for XC task queue
   *  @param[in] task_end   End iterator for XC task queue
   *
   *  @returns iterator to last XC task queue which was not kept in the
   *           allocated batch (!= task_end)
   */
  virtual host_task_iterator generate_buffers( integrator_term_tracker terms,
    const BasisSetMap& basis_map, host_task_iterator task_begin,
    host_task_iterator task_end ) = 0;

  /** Retreive XC integrands from device memory
   *
   *  TODO: this will depend on the integrand, we should refactor this
   *  to only allocate what is needed
   *
   *  TODO: this might be merged with reduction to allow for e.g. NCCL
   *
   *  @param[out] EXC  Integrated XC energy (host) for XC task
   *  @param[out] N_EL Integrated # electrons (host) for XC queue (accuracy metric)
   *  @param[out[ VXC  Integrated XC potential (host) for XC queue
   */
  virtual void retrieve_xc_integrands( double* EXC, double* N_EL,
    double* VXC, int32_t ldvxc ) = 0;


  virtual void copy_weights_to_tasks( host_task_iterator task_begin, host_task_iterator task_end ) = 0;
  virtual void populate_submat_maps ( size_t, host_task_iterator begin, host_task_iterator end, const BasisSetMap& ) = 0;

  virtual double* vxc_device_data() = 0;
  virtual double* exc_device_data() = 0;
  virtual double* nel_device_data() = 0;
  virtual type_erased_queue queue() = 0;


};

}
