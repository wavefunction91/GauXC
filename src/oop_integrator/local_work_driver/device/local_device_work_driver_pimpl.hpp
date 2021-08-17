#pragma once
#include "local_device_work_driver.hpp"


namespace GauXC {
namespace detail {

class LocalDeviceWorkDriverPIMPL;


/// Base class for local work drivers in Device execution spaces 
struct LocalDeviceWorkDriverPIMPL {

  using task_t         = LocalDeviceWorkDriver::task_t;         
  using task_container = LocalDeviceWorkDriver::task_container; 
  using task_iterator  = LocalDeviceWorkDriver::task_iterator;  

  /// Construct LocalDeviceWorkDriver instance in invalid state
  LocalDeviceWorkDriverPIMPL();

  /// Destructor (default)
  virtual ~LocalDeviceWorkDriverPIMPL() noexcept;

  LocalDeviceWorkDriverPIMPL( const LocalDeviceWorkDriverPIMPL& )     = delete;
  LocalDeviceWorkDriverPIMPL( LocalDeviceWorkDriverPIMPL&& ) noexcept = delete;


  // Public APIs

  virtual void partition_weights( task_t*, XCDeviceData* ) = 0;
  virtual void eval_collocation( task_t*, XCDeviceData* ) = 0;
  virtual void eval_collocation_gradient( task_t*, XCDeviceData* ) = 0;
  virtual void eval_xmat( task_t*, XCDeviceData* ) = 0;
  virtual void eval_uvvar_lda( task_t*, XCDeviceData* ) = 0;
  virtual void eval_uvvar_gga( task_t*, XCDeviceData* ) = 0;
  virtual void eval_zmat_lda_vxc( task_t*, XCDeviceData* ) = 0;
  virtual void eval_zmat_gga_vxc( task_t*, XCDeviceData* ) = 0;
  virtual void inc_vxc( task_t*, XCDeviceData* ) = 0;

};

}
}

