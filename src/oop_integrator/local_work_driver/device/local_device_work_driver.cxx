#include "local_host_work_driver_pimpl.hpp"
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
  if(not ptr) throw std::runtime_error(std::string("INVALID LocalDeviceWorkDriver PIMPL: ") + std::string(__PRETTY_FUNCTION__) );





// Partition weights
void LocalDeviceWorkDriver::partition_weights( task_t* tasks, XCDeviceData* device_data ) {
  

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->partition_weights(tasks,device_data);

}


// Collocation
void LocalDeviceWorkDriver::eval_collocation( task_t* tasks, XCDeviceData* device_data ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_collocation(tasks,device_data);

}


// Collocation Gradient
void LocalDeviceWorkDriver::eval_collocation_gradient( task_t* tasks, XCDeviceData* device_data ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_collocation_gradient(tasks,device_data);

}


// X matrix (P * B)
void LocalDeviceWorkDriver::eval_xmat(task_t* tasks, XCDeviceData* device_data ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_xmat(tasks,device_data);

}


// U/VVar LDA (density)
void LocalDeviceWorkDriver::eval_uvvar_lda(task_t* tasks, XCDeviceData* device_data ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_uvvar_lda(tasks,device_data);

}


// U/VVar GGA (density + grad, gamma)
void LocalDeviceWorkDriver::eval_uvvar_gga(task_t* tasks, XCDeviceData* device_data ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_uvvar_gga(tasks,device_data);

}

// Eval Z Matrix LDA VXC
void LocalDeviceWorkDriver::eval_zmat_lda_vxc(task_t* tasks, XCDeviceData* device_data ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_zmat_lda_vxc(tasks,device_data);

}

// Eval Z Matrix GGA VXC
void LocalDeviceWorkDriver::eval_zmat_gga_vxc(task_t* tasks, XCDeviceData* device_data ) { 

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_zmat_gga_vxc(tasks,device_data);

}


// Increment VXC by Z
void LocalDeviceWorkDriver::inc_vxc(task_t* tasks, XCDeviceData* device_data ) { 

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->inc_vxc(tasks,device_data);

}



}
