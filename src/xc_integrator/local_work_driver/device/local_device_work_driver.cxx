/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
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
#define FWD_TO_PIMPL_BOOL(NAME) \
void LocalDeviceWorkDriver::NAME( XCDeviceData* device_data, bool b ) { \
  throw_if_invalid_pimpl(pimpl_);                                       \
  pimpl_->NAME(device_data, b);                                         \
}

#define FWD_TO_PIMPL_DEN_ID(NAME) \
void LocalDeviceWorkDriver::NAME( XCDeviceData* device_data, density_id den ) { \
  throw_if_invalid_pimpl(pimpl_);                               \
  pimpl_->NAME(device_data, den);                               \
}

#define FWD_TO_PIMPL_DEN_ID_BOOL(NAME) \
void LocalDeviceWorkDriver::NAME( XCDeviceData* device_data, density_id den, bool b ) { \
  throw_if_invalid_pimpl(pimpl_);                               \
  pimpl_->NAME(device_data, den, b);                               \
}

#define FWD_TO_PIMPL_KS_SCHEME(NAME) \
void LocalDeviceWorkDriver::NAME( XCDeviceData* device_data, integrator_ks_scheme track ) { \
  throw_if_invalid_pimpl(pimpl_);                               \
  pimpl_->NAME(device_data, track);                               \
}
#define FWD_TO_PIMPL_KS_SCHEME_DEN_ID(NAME) \
void LocalDeviceWorkDriver::NAME( XCDeviceData* device_data, integrator_ks_scheme track, density_id den ) { \
  throw_if_invalid_pimpl(pimpl_);                               \
  pimpl_->NAME(device_data, track, den);                               \
}

FWD_TO_PIMPL(partition_weights)         // Partition weights

FWD_TO_PIMPL(eval_collocation)          // Collocation
FWD_TO_PIMPL(eval_collocation_gradient) // Collocation Gradient
FWD_TO_PIMPL(eval_collocation_hessian)  // Collocation Hessian
FWD_TO_PIMPL(eval_collocation_laplacian)  // Collocation Laplacian
FWD_TO_PIMPL(eval_collocation_lapgrad)  // Collocation Laplacian gradient


FWD_TO_PIMPL_KS_SCHEME(eval_uvars_lda)            // U variables LDA (rho)
FWD_TO_PIMPL_KS_SCHEME(eval_uvars_gga)            // U variables GGA (gamma)
FWD_TO_PIMPL_BOOL(eval_uvars_mgga)                // U variables MGGA (tau, lapl)
FWD_TO_PIMPL_DEN_ID_BOOL(eval_vvar)               // V variable (density + grad)

FWD_TO_PIMPL_KS_SCHEME_DEN_ID(eval_zmat_lda_vxc)         // Eval Z Matrix LDA VXC
FWD_TO_PIMPL_KS_SCHEME_DEN_ID(eval_zmat_gga_vxc)         // Eval Z Matrix GGA VXC
FWD_TO_PIMPL_BOOL(eval_zmat_mgga_vxc)                    // Eval Z Matrix mGGA VXC
FWD_TO_PIMPL_BOOL(eval_mmat_mgga_vxc)                    // Eval M Matrix mGGA VXC

FWD_TO_PIMPL(eval_exx_fmat)             // Eval EXX F Matrix
//FWD_TO_PIMPL(eval_exx_gmat)             // Eval EXX G Matrix


FWD_TO_PIMPL(inc_exc)
FWD_TO_PIMPL(inc_nel)
FWD_TO_PIMPL_DEN_ID_BOOL(inc_vxc)            // Increment VXC_I by Z

FWD_TO_PIMPL(inc_exx_k)     
FWD_TO_PIMPL(inc_exc_grad_lda)
FWD_TO_PIMPL(inc_exc_grad_gga)
FWD_TO_PIMPL_BOOL(inc_exc_grad_mgga)

FWD_TO_PIMPL_DEN_ID(symmetrize_vxc)
FWD_TO_PIMPL(symmetrize_exx_k)
FWD_TO_PIMPL(eval_exx_ek_screening_bfn_stats)

// X     = fac * P * B
// dX/dx = fac * P * dB/dx (do_grad)

void LocalDeviceWorkDriver::eval_xmat( double fac, XCDeviceData* device_data, bool do_grad, density_id den ) {
  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_xmat(fac, device_data, do_grad, den);
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

void LocalDeviceWorkDriver::eval_kern_exc_vxc_mgga( const functional_type& func,
  XCDeviceData* data) {
  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_kern_exc_vxc_mgga(func,data);
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
