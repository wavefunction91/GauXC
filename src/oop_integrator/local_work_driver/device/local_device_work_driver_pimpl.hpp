#pragma once
#include "local_device_work_driver.hpp"


namespace GauXC {
namespace detail {

class LocalDeviceWorkDriverPIMPL;


/// Base class for local work drivers in Device execution spaces 
struct LocalDeviceWorkDriverPIMPL {

  LocalDeviceWorkDriverPIMPL();
  virtual ~LocalDeviceWorkDriverPIMPL() noexcept;

  LocalDeviceWorkDriverPIMPL( const LocalDeviceWorkDriverPIMPL& )     = delete;
  LocalDeviceWorkDriverPIMPL( LocalDeviceWorkDriverPIMPL&& ) noexcept = delete;


  // Public APIs

  virtual void partition_weights( XCDeviceData* ) = 0;
  virtual void eval_collocation( XCDeviceData* ) = 0;
  virtual void eval_collocation_gradient( XCDeviceData* ) = 0;
  virtual void eval_xmat( XCDeviceData* ) = 0;
  virtual void eval_uvvar_lda( XCDeviceData* ) = 0;
  virtual void eval_uvvar_gga( XCDeviceData* ) = 0;
  virtual void eval_zmat_lda_vxc( XCDeviceData* ) = 0;
  virtual void eval_zmat_gga_vxc( XCDeviceData* ) = 0;
  virtual void inc_vxc( XCDeviceData* ) = 0;

  virtual std::unique_ptr<XCDeviceData> create_device_data() = 0;

};

}
}

