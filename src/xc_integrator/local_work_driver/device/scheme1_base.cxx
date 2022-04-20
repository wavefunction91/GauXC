#include "scheme1_base.hpp"
#include "device/common/zmat_vxc.hpp"
#include "device/common/collocation_device.hpp"
#include "device/common/device_blas.hpp"
#include "device/common/xc_functional_eval_wrapper.hpp"
#include "device/common/uvvars.hpp"
#include "device/common/pack_submat.hpp"
#include "device/common/inc_potential.hpp"
#include "device/common/symmetrize_mat.hpp"
#include "device/common/increment_exc_grad.hpp"

namespace GauXC {

void AoSScheme1Base::eval_zmat_lda_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

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
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

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
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

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
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

#ifdef GAUXC_ENABLE_HIP
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
#else
  auto aos_stack     = data->aos_stack;

  auto max_l = data->l_batched_shell_to_task.size() - 1;
  eval_collocation_shell_to_task_gradient( max_l, 
    data->l_batched_shell_to_task.data(), aos_stack.device_tasks,
    data->device_backend_->queue() );
#endif
  
}

void AoSScheme1Base::eval_collocation_hessian( XCDeviceData* _data ) {
#ifdef GAUXC_ENABLE_HIP
  GAUXC_GENERIC_EXCEPTION("Hessian NYI for HIP Backends");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto aos_stack     = data->aos_stack;

  auto max_l = data->l_batched_shell_to_task.size() - 1;
  eval_collocation_shell_to_task_hessian( max_l, 
    data->l_batched_shell_to_task.data(), aos_stack.device_tasks,
    data->device_backend_->queue() );
#endif
}





void AoSScheme1Base::inc_exc( XCDeviceData* _data ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;
  gdot( data->device_backend_->master_blas_handle(), data->total_npts_task_batch,
    base_stack.eps_eval_device, 1, base_stack.den_eval_device, 1, 
    static_stack.acc_scr_device, static_stack.exc_device );

}
void AoSScheme1Base::inc_nel( XCDeviceData* _data ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;
  gdot( data->device_backend_->master_blas_handle(), data->total_npts_task_batch,
    base_stack.weights_device, 1, base_stack.den_eval_device, 1, 
    static_stack.acc_scr_device, static_stack.nel_device );

}















void AoSScheme1Base::eval_uvvar_lda( XCDeviceData* _data ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  // Zero density
  auto base_stack    = data->base_stack;
  data->device_backend_->set_zero_async_master_queue( data->total_npts_task_batch, base_stack.den_eval_device, "Den Zero" );
    

  // Evaluate U variables
  auto aos_stack     = data->aos_stack;
  eval_uvvars_lda( ntasks, nbe_max, npts_max, 
    aos_stack.device_tasks, data->device_backend_->queue() );

}




void AoSScheme1Base::eval_uvvar_gga( XCDeviceData* _data ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  // Zero density + gradient
  auto base_stack    = data->base_stack;
  data->device_backend_->set_zero_async_master_queue( data->total_npts_task_batch, base_stack.den_eval_device, "Den Zero" );
  data->device_backend_->set_zero_async_master_queue( data->total_npts_task_batch, base_stack.den_x_eval_device, "DenX Zero" );
  data->device_backend_->set_zero_async_master_queue( data->total_npts_task_batch, base_stack.den_y_eval_device, "DenY Zero" );
  data->device_backend_->set_zero_async_master_queue( data->total_npts_task_batch, base_stack.den_z_eval_device, "DenZ Zero" );

  // Evaluate U variables
  auto aos_stack     = data->aos_stack;
  eval_uvvars_gga( ntasks, data->total_npts_task_batch, nbe_max, npts_max, 
    aos_stack.device_tasks, base_stack.den_x_eval_device, base_stack.den_y_eval_device,
    base_stack.den_z_eval_device, base_stack.gamma_eval_device, 
    data->device_backend_->queue() );

}









void AoSScheme1Base::eval_kern_exc_vxc_lda( const functional_type& func, 
  XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  if( !func.is_lda() ) GAUXC_GENERIC_EXCEPTION("XC Kernel not LDA!");

  auto base_stack    = data->base_stack;
  GauXC::eval_kern_exc_vxc_lda( func, data->total_npts_task_batch, 
    base_stack.den_eval_device, base_stack.eps_eval_device, 
    base_stack.vrho_eval_device, data->device_backend_->queue() );

  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.eps_eval_device, 1 );
  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.vrho_eval_device, 1 );

}


void AoSScheme1Base::eval_kern_exc_vxc_gga( const functional_type& func, 
  XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  if( !func.is_gga() ) GAUXC_GENERIC_EXCEPTION("XC Kernel not GGA!");

  auto base_stack    = data->base_stack;
  GauXC::eval_kern_exc_vxc_gga( func, data->total_npts_task_batch, 
    base_stack.den_eval_device, base_stack.gamma_eval_device, 
    base_stack.eps_eval_device, base_stack.vrho_eval_device, 
    base_stack.vgamma_eval_device, data->device_backend_->queue() );

  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.eps_eval_device, 1 );
  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.vrho_eval_device, 1 );
  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.vgamma_eval_device, 1 );

}










void AoSScheme1Base::eval_xmat( XCDeviceData* _data, bool do_grad ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  // Pack density matrix 
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  pack_submat( ntasks, aos_stack.device_tasks, static_stack.dmat_device, 
    nbf, submat_block_size, data->device_backend_->queue() );


  // Sync blas streams with master stream
  data->device_backend_->sync_blas_pool_with_master();

  auto do_gemm = [&]( auto& handle, size_t npts, size_t nbe, auto* bf_ptr, auto* den_ptr, int ldden, auto* x_ptr ) {
    gemm( handle, DeviceBlasOp::NoTrans, DeviceBlasOp::NoTrans, npts, nbe, nbe, 1., bf_ptr, npts,
      den_ptr, ldden, 0., x_ptr, npts ); 
  };

  // Launch GEMM in round-robin
  const auto n_blas_streams = data->device_backend_->blas_pool_size();
  //size_t nsingle = 0;
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];
      auto den_ptr = task.ncut > 1 ? task.nbe_scr : static_stack.dmat_device + task.ibf_begin*(nbf+1);
      int  ldden   = task.ncut > 1 ? task.nbe : nbf;
      auto handle = data->device_backend_->blas_pool_handle( iT % n_blas_streams );
      do_gemm( handle, task.npts, task.nbe, task.bf, den_ptr, ldden, task.zmat );
      if( do_grad ) {
        do_gemm( handle, task.npts, task.nbe, task.dbfx, den_ptr, ldden, task.xmat_x );
        do_gemm( handle, task.npts, task.nbe, task.dbfy, den_ptr, ldden, task.xmat_y );
        do_gemm( handle, task.npts, task.nbe, task.dbfz, den_ptr, ldden, task.xmat_z );
      }
  }

  // Record completion of BLAS ops on master stream
  data->device_backend_->sync_master_with_blas_pool();

}















void AoSScheme1Base::inc_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  // Sync blas streams with master stream
  data->device_backend_->sync_blas_pool_with_master();

  // Launch SYR2K in round robin
  const auto n_blas_streams = data->device_backend_->blas_pool_size();
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];
    syr2k( data->device_backend_->blas_pool_handle(iT % n_blas_streams), 
      DeviceBlasUplo::Lower, DeviceBlasOp::Trans, task.nbe, task.npts, 1.,
      task.bf, task.npts, task.zmat, task.npts, 0., task.nbe_scr,
      task.nbe );
  }

  // Record completion of BLAS ops on master stream
  data->device_backend_->sync_master_with_blas_pool();

  // Increment global VXC
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  task_inc_potential( ntasks, aos_stack.device_tasks, 
    static_stack.vxc_device, nbf, submat_block_size, 
    data->device_backend_->queue() );
}













void AoSScheme1Base::symmetrize_vxc( XCDeviceData* _data) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  const auto nbf = data->global_dims.nbf;
  auto static_stack  = data->static_stack;
  symmetrize_matrix( nbf, static_stack.vxc_device, nbf, 
    data->device_backend_->queue() ); 

}




void AoSScheme1Base::inc_exc_grad_lda( XCDeviceData* _data ) {
#ifdef GAUXC_ENABLE_HIP
  GAUXC_GENERIC_EXCEPTION("LDA Grad NYI for HIP Backends");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  const auto nshell = data->global_dims.nshells;
  increment_exc_grad_lda( nshell, 
    data->shell_to_task_stack.shell_to_task_device, 
    data->aos_stack.device_tasks,
    data->static_stack.exc_grad_device,
    data->device_backend_->queue() ); 
#endif
}

void AoSScheme1Base::inc_exc_grad_gga( XCDeviceData* _data ) {
#ifdef GAUXC_ENABLE_HIP
  GAUXC_GENERIC_EXCEPTION("GGA Grad NYI for HIP Backends");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  const auto nshell = data->global_dims.nshells;
  increment_exc_grad_gga( nshell, 
    data->shell_to_task_stack.shell_to_task_device, 
    data->aos_stack.device_tasks,
    data->static_stack.exc_grad_device,
    data->device_backend_->queue() ); 
#endif
}


void AoSScheme1Base::inc_exx_k( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();
  GAUXC_GENERIC_EXCEPTION("NYI");

}

void AoSScheme1Base::symmetrize_exx_k( XCDeviceData* _data) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  GAUXC_GENERIC_EXCEPTION("NYI");

}

void AoSScheme1Base::eval_exx_fmat( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  GAUXC_GENERIC_EXCEPTION("NYI");

}


void AoSScheme1Base::eval_exx_gmat( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  GAUXC_GENERIC_EXCEPTION("NYI");

}

}
