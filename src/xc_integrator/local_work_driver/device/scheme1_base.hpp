#pragma once
#include "device/local_device_work_driver_pimpl.hpp"
#include "device/scheme1_data_base.hpp"

namespace GauXC {

struct AoSScheme1Base : public detail::LocalDeviceWorkDriverPIMPL {

  // Device Common APIs
  void eval_zmat_lda_vxc( XCDeviceData* ) override final;
  void eval_zmat_gga_vxc( XCDeviceData* ) override final;

  using Data = Scheme1DataBase;

};

}
