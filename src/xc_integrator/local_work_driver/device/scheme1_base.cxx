#include "scheme1_base.hpp"
#include "device/common/zmat_vxc.hpp"
#include "device/common/collocation_device.hpp"
#include "device/common/device_blas.hpp"

namespace GauXC {

void AoSScheme1Base::eval_zmat_lda_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  if( not data->device_backend_ ) throw std::runtime_error("INVALID DEVICE BACKEND");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  auto aos_stack     = data->aos_stack;
  zmat_lda_vxc( ntasks, nbe_max, npts_max, aos_stack.device_tasks,
    data->device_backend_->queue() );

}

void AoSScheme1Base::eval_zmat_gga_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  if( not data->device_backend_ ) throw std::runtime_error("INVALID DEVICE BACKEND");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  auto aos_stack     = data->aos_stack;
  zmat_gga_vxc( ntasks, nbe_max, npts_max, aos_stack.device_tasks,
    data->device_backend_->queue() );

}



void AoSScheme1Base::eval_collocation( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  if( not data->device_backend_ ) throw std::runtime_error("INVALID DEVICE BACKEND");

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  size_t npts_max = 0, nshells_max = 0;
  for( auto& task : tasks ) {
    npts_max    = std::max( npts_max, task.npts );
    nshells_max = std::max( nshells_max, task.nshells );
  }

  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  eval_collocation_masked_combined( ntasks, npts_max, nshells_max,
    static_stack.shells_device, aos_stack.device_tasks, 
    data->device_backend_->queue() );

}

void AoSScheme1Base::eval_collocation_gradient( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  if( not data->device_backend_ ) throw std::runtime_error("INVALID DEVICE BACKEND");

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  size_t npts_max = 0, nshells_max = 0;
  for( auto& task : tasks ) {
    npts_max    = std::max( npts_max, task.npts );
    nshells_max = std::max( nshells_max, task.nshells );
  }

  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  eval_collocation_masked_combined_deriv1( ntasks, npts_max, nshells_max,
    static_stack.shells_device, aos_stack.device_tasks, 
    data->device_backend_->queue() );
  
}





void AoSScheme1Base::inc_exc( XCDeviceData* _data ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  if( not data->device_backend_ ) throw std::runtime_error("INVALID DEVICE BACKEND");

  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;
  gdot( data->device_backend_->master_blas_handle(), data->total_npts_task_batch,
    base_stack.eps_eval_device, 1, base_stack.den_eval_device, 1, 
    static_stack.acc_scr_device, static_stack.exc_device );

}
void AoSScheme1Base::inc_nel( XCDeviceData* _data ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  if( not data->device_backend_ ) throw std::runtime_error("INVALID DEVICE BACKEND");

  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;
  gdot( data->device_backend_->master_blas_handle(), data->total_npts_task_batch,
    base_stack.weights_device, 1, base_stack.den_eval_device, 1, 
    static_stack.acc_scr_device, static_stack.nel_device );

}

}
