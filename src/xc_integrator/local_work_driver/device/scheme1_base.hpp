#pragma once
#include "device/local_device_work_driver_pimpl.hpp"
#include "device/scheme1_data_base.hpp"

namespace GauXC {

struct AoSScheme1Base : public detail::LocalDeviceWorkDriverPIMPL {

  // Device Common APIs (final overrides)
  void eval_collocation( XCDeviceData* ) override final;
  void eval_collocation_gradient( XCDeviceData* ) override final;
  void eval_uvvar_lda( XCDeviceData* ) override final;
  void eval_uvvar_gga( XCDeviceData* ) override final;
  void eval_zmat_lda_vxc( XCDeviceData* ) override final;
  void eval_zmat_gga_vxc( XCDeviceData* ) override final;
  void eval_kern_exc_vxc_lda( const functional_type&, XCDeviceData* ) override final;
  void eval_kern_exc_vxc_gga( const functional_type&, XCDeviceData* ) override final;
  void inc_exc( XCDeviceData* ) override final;
  void inc_nel( XCDeviceData* ) override final;
  void symmetrize_vxc( XCDeviceData* ) override final;


  // Overridable APIs
  virtual void eval_xmat( XCDeviceData* ) override;
  virtual void inc_vxc( XCDeviceData* ) override;

  using Data = Scheme1DataBase;

  virtual ~AoSScheme1Base() = default;
};

}
