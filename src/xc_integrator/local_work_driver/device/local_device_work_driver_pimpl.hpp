/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "local_device_work_driver.hpp"


namespace GauXC {
namespace detail {

struct LocalDeviceWorkDriverPIMPL;


/// Base class for local work drivers in Device execution spaces 
struct LocalDeviceWorkDriverPIMPL {

  using host_task_iterator = LocalDeviceWorkDriver::host_task_iterator;

  LocalDeviceWorkDriverPIMPL();
  virtual ~LocalDeviceWorkDriverPIMPL() noexcept;

  LocalDeviceWorkDriverPIMPL( const LocalDeviceWorkDriverPIMPL& )     = delete;
  LocalDeviceWorkDriverPIMPL( LocalDeviceWorkDriverPIMPL&& ) noexcept = delete;


  // Public APIs

  virtual void partition_weights( XCDeviceData* ) = 0;
  virtual void eval_collocation( XCDeviceData* ) = 0;
  virtual void eval_collocation_gradient( XCDeviceData* ) = 0;
  virtual void eval_collocation_hessian( XCDeviceData* ) = 0;
  virtual void eval_xmat( double fac, XCDeviceData*, bool do_grad, density_id den ) = 0;
  virtual void eval_den( XCDeviceData*, density_id ) = 0;
  virtual void eval_exx_fmat( XCDeviceData* ) = 0;
  //virtual void eval_exx_gmat( XCDeviceData* ) = 0;
  virtual void eval_exx_gmat( XCDeviceData*, const BasisSetMap& ) = 0;
  virtual void eval_uvvar_lda( XCDeviceData*, integrator_term_tracker ) = 0;
  virtual void eval_uvvar_gga( XCDeviceData*, integrator_term_tracker ) = 0;
  virtual void eval_kern_exc_vxc_lda( const functional_type&, XCDeviceData* ) = 0;
  virtual void eval_kern_exc_vxc_gga( const functional_type&, XCDeviceData* ) = 0;
  virtual void eval_zmat_lda_vxc_rks( XCDeviceData* ) = 0;
  virtual void eval_zmat_gga_vxc_rks( XCDeviceData* ) = 0;
  virtual void eval_zmat_lda_vxc_uks( XCDeviceData*, density_id ) = 0;
  virtual void eval_zmat_gga_vxc_uks( XCDeviceData* ) = 0;
  virtual void eval_zmat_lda_vxc_gks( XCDeviceData* ) = 0;
  virtual void eval_zmat_gga_vxc_gks( XCDeviceData* ) = 0;
  virtual void inc_exc( XCDeviceData* ) = 0;
  virtual void inc_nel( XCDeviceData* ) = 0;
  virtual void inc_vxc( XCDeviceData* , density_id) = 0;
  virtual void inc_vxc( XCDeviceData* ) = 0;
  virtual void inc_exc_grad_lda( XCDeviceData* ) = 0;
  virtual void inc_exc_grad_gga( XCDeviceData* ) = 0;
  virtual void inc_exx_k( XCDeviceData* ) = 0;
  virtual void symmetrize_vxc( XCDeviceData*, density_id ) = 0;
  virtual void symmetrize_exx_k( XCDeviceData* ) = 0;

  virtual void eval_exx_ek_screening_bfn_stats( XCDeviceData* ) = 0;
  virtual void exx_ek_shellpair_collision( double eps_E, double eps_K, 
    XCDeviceData*, host_task_iterator, host_task_iterator, 
    const ShellPairCollection<double>&) = 0;

  virtual std::unique_ptr<XCDeviceData> create_device_data(const DeviceRuntimeEnvironment&) = 0;

};

}
}

