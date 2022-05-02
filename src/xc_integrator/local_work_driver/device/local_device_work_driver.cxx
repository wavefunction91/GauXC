#include "local_device_work_driver_pimpl.hpp"
#include <stdexcept>

namespace GauXC {

LocalDeviceWorkDriver::LocalDeviceWorkDriver() : 
  pimpl_(nullptr) { }
LocalDeviceWorkDriver::LocalDeviceWorkDriver(pimpl_type&& ptr) :
  pimpl_( std::move(ptr) ){ }

LocalDeviceWorkDriver::~LocalDeviceWorkDriver() noexcept = default;

LocalDeviceWorkDriver::LocalDeviceWorkDriver( LocalDeviceWorkDriver&& other ) noexcept :
  pimpl_(std::move(other.pimpl_)) { }

#define throw_if_invalid_pimpl(ptr) \
  if(not ptr) GAUXC_PIMPL_NOT_INITIALIZED()



#define FWD_TO_PIMPL(NAME) \
void LocalDeviceWorkDriver::NAME( XCDeviceData* device_data ) { \
  throw_if_invalid_pimpl(pimpl_);                               \
  pimpl_->NAME(device_data);                                    \
}


FWD_TO_PIMPL(partition_weights)         // Partition weights

FWD_TO_PIMPL(eval_collocation)          // Collocation
FWD_TO_PIMPL(eval_collocation_gradient) // Collocation Gradient
FWD_TO_PIMPL(eval_collocation_hessian)  // Collocation Hessian

FWD_TO_PIMPL(eval_uvvar_lda)            // U/VVar LDA (density)
FWD_TO_PIMPL(eval_uvvar_gga)            // U/VVar GGA (density + grad, gamma)

FWD_TO_PIMPL(eval_zmat_lda_vxc)         // Eval Z Matrix LDA VXC
FWD_TO_PIMPL(eval_zmat_gga_vxc)         // Eval Z Matrix GGA VXC

FWD_TO_PIMPL(eval_exx_fmat)             // Eval EXX F Matrix
//FWD_TO_PIMPL(eval_exx_gmat)             // Eval EXX G Matrix

FWD_TO_PIMPL(inc_exc)
FWD_TO_PIMPL(inc_nel)
FWD_TO_PIMPL(inc_vxc)                   // Increment VXC by Z 
FWD_TO_PIMPL(inc_exx_k)     
FWD_TO_PIMPL(inc_exc_grad_lda)
FWD_TO_PIMPL(inc_exc_grad_gga)

FWD_TO_PIMPL(symmetrize_vxc)
FWD_TO_PIMPL(symmetrize_exx_k)

// X     = P * B
// dX/dx = P * dB/dx (do_grad)
void LocalDeviceWorkDriver::eval_xmat( XCDeviceData* device_data, bool do_grad ) {
  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_xmat(device_data, do_grad);
}

void LocalDeviceWorkDriver::eval_exx_gmat( XCDeviceData* device_data, 
  const BasisSetMap& basis_map) {
  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_exx_gmat(device_data, basis_map);
}

void LocalDeviceWorkDriver::eval_kern_exc_vxc_lda( const functional_type& func,
  XCDeviceData* data) {
  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_kern_exc_vxc_lda(func,data);
}

void LocalDeviceWorkDriver::eval_kern_exc_vxc_gga( const functional_type& func,
  XCDeviceData* data) {
  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_kern_exc_vxc_gga(func,data);
}


std::unique_ptr<XCDeviceData> LocalDeviceWorkDriver::create_device_data() {
  throw_if_invalid_pimpl(pimpl_);
  return pimpl_->create_device_data();
}

}
