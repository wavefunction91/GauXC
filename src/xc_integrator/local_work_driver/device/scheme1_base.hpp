/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "device/local_device_work_driver_pimpl.hpp"
#include "device/scheme1_data_base.hpp"

namespace GauXC {

struct AoSScheme1Base : public detail::LocalDeviceWorkDriverPIMPL {

  // Device Common APIs (final overrides)
  void eval_collocation( XCDeviceData* ) override final;
  void eval_collocation_gradient( XCDeviceData* ) override final;
  void eval_collocation_hessian( XCDeviceData* ) override final;
  void eval_uvvar_lda_rks( XCDeviceData* ) override final;
  void eval_uvvar_gga_rks( XCDeviceData* ) override final;
  void eval_zmat_lda_vxc_rks( XCDeviceData* ) override final;
  void eval_zmat_gga_vxc_rks( XCDeviceData* ) override final;

  void eval_uvvar_lda_uks( XCDeviceData* ) override final;
  void eval_uvvar_gga_uks( XCDeviceData* ) override final;
  void eval_zmat_lda_vxc_uks( XCDeviceData*, density_id ) override final;
  void eval_zmat_gga_vxc_uks( XCDeviceData* ) override final;

  void eval_uvvar_lda_gks( XCDeviceData* ) override final;
  void eval_uvvar_gga_gks( XCDeviceData* ) override final;
  void eval_zmat_lda_vxc_gks( XCDeviceData* ) override final;
  void eval_zmat_gga_vxc_gks( XCDeviceData* ) override final;

  void eval_kern_exc_vxc_lda( const functional_type&, XCDeviceData* ) override final;
  void eval_kern_exc_vxc_gga( const functional_type&, XCDeviceData* ) override final;
  void inc_exc( XCDeviceData* ) override final;
  void inc_nel( XCDeviceData* ) override final;
  void inc_exc_grad_lda( XCDeviceData* ) override final;
  void inc_exc_grad_gga( XCDeviceData* ) override final;
  void symmetrize_vxc( XCDeviceData* ) override final;
  void symmetrize_vxc( XCDeviceData* , density_id) override final;
  void symmetrize_exx_k( XCDeviceData* ) override final;
  //void eval_exx_gmat( XCDeviceData* ) override final;
  void eval_exx_gmat( XCDeviceData*, const BasisSetMap& ) override final;

  void eval_exx_ek_screening_bfn_stats( XCDeviceData* ) override final;
  void exx_ek_shellpair_collision( double eps_E, double eps_K, 
    XCDeviceData*, host_task_iterator, host_task_iterator,
    const ShellPairCollection<double>& ) override final;


  // Overridable APIs
  virtual void eval_xmat( double fac, XCDeviceData*, bool do_grad, density_id den_sel ) override;
  virtual void eval_exx_fmat( XCDeviceData* ) override;
  virtual void inc_vxc( XCDeviceData*, density_id) override;
  virtual void inc_vxc( XCDeviceData* ) override; 
  virtual void inc_exx_k( XCDeviceData* ) override;

  void eval_den( XCDeviceData*, density_id) override final;

  using Data = Scheme1DataBase;

  AoSScheme1Base();
  virtual ~AoSScheme1Base() noexcept;

  double* dev_boys_table = nullptr;
};

}
