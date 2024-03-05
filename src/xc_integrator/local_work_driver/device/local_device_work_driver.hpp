/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/xc_integrator/local_work_driver.hpp>

#include <memory>
#include <gauxc/molmeta.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/basisset_map.hpp>
#include <gauxc/xc_task.hpp>

#include "device/xc_device_data.hpp"
#include <gauxc/runtime_environment/fwd.hpp>

namespace GauXC {
namespace detail {

struct LocalDeviceWorkDriverPIMPL;

}

/// Base class for local work drivers in Device execution spaces 
class LocalDeviceWorkDriver : public LocalWorkDriver {

  using pimpl_type = std::unique_ptr<detail::LocalDeviceWorkDriverPIMPL>;

public:

  using host_task_iterator = std::vector<XCTask>::iterator;

  /// Construct LocalDeviceWorkDriver instance in invalid state
  LocalDeviceWorkDriver();

  /** Construct LocalDeviceWorkDriver instance given implementation pointer
   *  @param[in] ptr Pointer to implementation
   */
  LocalDeviceWorkDriver( pimpl_type&& ptr );

  /// Destructor (default)
  ~LocalDeviceWorkDriver() noexcept;

  // Remove copy ctor
  LocalDeviceWorkDriver( const LocalDeviceWorkDriver& ) = delete;

  /** Construct LocalDeviceWorkDriver by transferring ownership
   *  @param[in] other LocalDeviceWorkDriver instance to take ownership
   */
  LocalDeviceWorkDriver( LocalDeviceWorkDriver&& other ) noexcept;


  // Public APIs

  void partition_weights( XCDeviceData* );

  void eval_collocation( XCDeviceData* );
  void eval_collocation_gradient( XCDeviceData* );
  void eval_collocation_hessian( XCDeviceData* );

  void eval_xmat( double fac, XCDeviceData*, bool do_grad, density_id den );
  
  void eval_uvars_lda( XCDeviceData*, integrator_ks_scheme );
  void eval_uvars_gga( XCDeviceData*, integrator_ks_scheme );
  void eval_vvar( XCDeviceData*, bool, density_id );

  void eval_kern_exc_vxc_lda( const functional_type&, XCDeviceData* );
  void eval_kern_exc_vxc_gga( const functional_type&, XCDeviceData* );


  void eval_zmat_lda_vxc( XCDeviceData*, integrator_ks_scheme,  density_id );
  void eval_zmat_gga_vxc( XCDeviceData*, integrator_ks_scheme,  density_id );


  void eval_exx_fmat( XCDeviceData* );
  void eval_exx_gmat( XCDeviceData*, const BasisSetMap& );

  void inc_exc( XCDeviceData* );
  void inc_nel( XCDeviceData* );
  void inc_vxc( XCDeviceData*, density_id );
  void inc_exc_grad_lda( XCDeviceData* );
  void inc_exc_grad_gga( XCDeviceData* );
  void inc_exx_k( XCDeviceData* );

  void eval_exx_ek_screening_bfn_stats( XCDeviceData* );
  void exx_ek_shellpair_collision( double eps_E, double eps_K, XCDeviceData*, 
    host_task_iterator, host_task_iterator, const ShellPairCollection<double>& );

  void symmetrize_vxc( XCDeviceData*, density_id );
  void symmetrize_exx_k( XCDeviceData* );

  std::unique_ptr<XCDeviceData> create_device_data(const DeviceRuntimeEnvironment&);

private: 

  pimpl_type pimpl_; ///< Implementation

};

}

