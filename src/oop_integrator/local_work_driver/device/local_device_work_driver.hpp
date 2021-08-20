#pragma once
#include <gauxc/oop_xc_integrator/local_work_driver.hpp>

#include <memory>
#include <gauxc/molmeta.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/xc_task.hpp>

#include "device/xc_device_data.hpp"

namespace GauXC {
namespace detail {

class LocalDeviceWorkDriverPIMPL;

}

/// Base class for local work drivers in Device execution spaces 
class LocalDeviceWorkDriver : public LocalWorkDriver {

  using pimpl_type = std::unique_ptr<detail::LocalDeviceWorkDriverPIMPL>;

public:

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

  void eval_xmat( XCDeviceData* );

  void eval_uvvar_lda( XCDeviceData* );
  void eval_uvvar_gga( XCDeviceData* );

  void eval_kern_exc_vxc_lda( const functional_type&, XCDeviceData* );
  void eval_kern_exc_vxc_gga( const functional_type&, XCDeviceData* );

  void eval_zmat_lda_vxc( XCDeviceData* );
  void eval_zmat_gga_vxc( XCDeviceData* );

  void inc_exc( XCDeviceData* );
  void inc_nel( XCDeviceData* );
  void inc_vxc( XCDeviceData* );

  std::unique_ptr<XCDeviceData> create_device_data();

private: 

  pimpl_type pimpl_; ///< Implementation

};

}

