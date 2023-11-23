/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
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

FWD_TO_PIMPL(eval_uvvar_lda_rks)            // U/VVar LDA (density)
FWD_TO_PIMPL(eval_uvvar_gga_rks)            // U/VVar GGA (density + grad, gamma)

FWD_TO_PIMPL(eval_uvvar_lda_uks)            // U/VVar LDA (density)
FWD_TO_PIMPL(eval_uvvar_gga_uks)            // U/VVar GGA (density + grad, gamma)

FWD_TO_PIMPL(eval_uvvar_lda_gks)            // U/VVar LDA (density)
FWD_TO_PIMPL(eval_uvvar_gga_gks)            // U/VVar GGA (density + grad, gamma)

FWD_TO_PIMPL(eval_zmat_lda_vxc_rks)         // Eval Z Matrix LDA VXC
FWD_TO_PIMPL(eval_zmat_gga_vxc_rks)         // Eval Z Matrix GGA VXC

FWD_TO_PIMPL(eval_zmat_lda_vxc_uks)         // Eval Z Matrix LDA VXC
FWD_TO_PIMPL(eval_zmat_gga_vxc_uks)         // Eval Z Matrix GGA VXC

FWD_TO_PIMPL(eval_zmat_lda_vxc_gks)         // Eval Z Matrix LDA VXC
FWD_TO_PIMPL(eval_zmat_gga_vxc_gks)         // Eval Z Matrix GGA VXC

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
FWD_TO_PIMPL(eval_exx_ek_screening_bfn_stats)

// X     = fac * P * B
// dX/dx = fac * P * dB/dx (do_grad)
void LocalDeviceWorkDriver::eval_xmat( double fac, XCDeviceData* device_data, bool do_grad ) {
  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_xmat(fac, device_data, do_grad);
}

void LocalDeviceWorkDriver::eval_xmat( double fac, XCDeviceData* device_data, bool do_grad, density_id den ) {
  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_xmat(fac, device_data, do_grad, den);
}

void LocalDeviceWorkDriver::eval_den(  XCDeviceData* device_data, density_id den ) {
  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_den( device_data, den);
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


std::unique_ptr<XCDeviceData> LocalDeviceWorkDriver::create_device_data(const DeviceRuntimeEnvironment& rt) {
  throw_if_invalid_pimpl(pimpl_);
  return pimpl_->create_device_data(rt);
}

void LocalDeviceWorkDriver::exx_ek_shellpair_collision( double eps_E, double eps_K, 
  XCDeviceData* device_data, host_task_iterator tb, host_task_iterator te,
  const ShellPairCollection<double>& shpairs ) {
  throw_if_invalid_pimpl(pimpl_);
  pimpl_->exx_ek_shellpair_collision( eps_E, eps_K, device_data, tb, te, shpairs );
}

}
