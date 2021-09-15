#include "hip_aos_scheme1.hpp"
#include "kernels/collocation_device.hpp"
#include "device/hip/hip_backend.hpp"
#include "kernels/grid_to_center.hpp"
#include "kernels/hip_ssf_1d.hpp"
//#include "hip_aos_scheme1_weights.hpp"
#include "kernels/hipblas_extensions.hpp"
//#include "kernels/uvvars.hpp"
#include "kernels/zmat_vxc.hpp"
//#include "kernels/pack_submat.hpp"
//#include "kernels/hip_inc_potential.hpp"

namespace GauXC {

std::unique_ptr<XCDeviceData> HipAoSScheme1::create_device_data() {
  return std::make_unique<Data>();
}

 
void HipAoSScheme1::partition_weights( XCDeviceData* _data ) {
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  const auto ldatoms = data->get_ldatoms();
  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;
  auto scheme1_stack = data->scheme1_stack;

  // Compute distances from grid to atomic centers
  compute_grid_to_center_dist( data->total_npts_task_batch, data->global_dims.natoms,
    static_stack.coords_device, base_stack.points_device, 
    scheme1_stack.dist_scratch_device, ldatoms, *device_backend->master_stream );

  // Modify weights
  partition_weights_ssf_1d( data->total_npts_task_batch, data->global_dims.natoms,
    static_stack.rab_device, ldatoms, static_stack.coords_device, 
    scheme1_stack.dist_scratch_device, ldatoms, scheme1_stack.iparent_device, 
    scheme1_stack.dist_nearest_device, base_stack.weights_device,
    *device_backend->master_stream );

}


void HipAoSScheme1::eval_collocation( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

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
    *device_backend->master_stream );

}

void HipAoSScheme1::eval_collocation_gradient( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

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
    *device_backend->master_stream );
  
}

#if 0
void HipAoSScheme1::eval_xmat( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  // Pack density matrix 
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  pack_submat( ntasks, aos_stack.device_tasks, data->dmat_device, nbf, 
    submat_block_size, *device_backend->master_stream );


  // Sync blas streams with master stream
  util::hip_event master_event;
  master_event.record( *device_backend->master_stream );
  for( auto& stream : device_backend->blas_streams ) stream.wait( master_event );

  // Launch GEMM in round-robin
  const auto n_blas_streams = device_backend->blas_streams.size();
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];
    gemm( device_backend->blas_handles[iT % n_blas_streams], HIPBLAS_OP_N,
      HIPBLAS_OP_N, task.npts, task.nbe, task.nbe, 1., task.bf, task.npts,
      task.nbe_scr, task.nbe, 0., task.zmat, task.npts );
  }

  // Record completion of BLAS ops on master stream
  std::vector< util::hip_event > blas_events( n_blas_streams );
  for( size_t iS = 0; iS < n_blas_streams; ++iS )
    blas_events[iS].record( device_backend->blas_streams[iS] );

  for( auto& event : blas_events )
    device_backend->master_stream->wait( event );

}






void HipAoSScheme1::eval_uvvar_lda( XCDeviceData* _data ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  // Zero density
  util::hip_set_zero_async( data->total_npts_task_batch, data->den_eval_device,
    *device_backend->master_stream, "Den Zero" );

  // Evaluate U variables
  eval_uvvars_lda_hip( ntasks, nbe_max, npts_max, aos_stack.device_tasks,
    *device_backend->master_stream );

}




void HipAoSScheme1::eval_uvvar_gga( XCDeviceData* _data ){


  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  // Zero density + gradient
  util::hip_set_zero_async( data->total_npts_task_batch, data->den_eval_device,
    *device_backend->master_stream, "Den Zero" );
  util::hip_set_zero_async( data->total_npts_task_batch, data->den_x_eval_device,
    *device_backend->master_stream, "DenX Zero" );
  util::hip_set_zero_async( data->total_npts_task_batch, data->den_y_eval_device,
    *device_backend->master_stream, "DenY Zero" );
  util::hip_set_zero_async( data->total_npts_task_batch, data->den_z_eval_device,
    *device_backend->master_stream, "DenZ Zero" );

  // Evaluate U variables
  eval_uvvars_gga_hip( ntasks, data->total_npts_task_batch, nbe_max, npts_max, 
    aos_stack.device_tasks, data->den_x_eval_device, data->den_y_eval_device,
    data->den_z_eval_device, data->gamma_eval_device, 
    *device_backend->master_stream );

}
#endif

void HipAoSScheme1::eval_kern_exc_vxc_lda( const functional_type& func, 
  XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  if( !func.is_lda() ) throw std::runtime_error("XC Kernel not LDA!");

  auto base_stack    = data->base_stack;
  func.eval_exc_vxc_device( data->total_npts_task_batch, 
    base_stack.den_eval_device, base_stack.eps_eval_device, 
    base_stack.vrho_eval_device, *device_backend->master_stream );

  hadamard_product( *device_backend->master_handle, 
    data->total_npts_task_batch, 1, base_stack.weights_device, 1, 
    base_stack.eps_eval_device, 1 );
  hadamard_product( *device_backend->master_handle, 
    data->total_npts_task_batch, 1, base_stack.weights_device, 1, 
    base_stack.vrho_eval_device, 1 );

}


void HipAoSScheme1::eval_kern_exc_vxc_gga( const functional_type& func, 
  XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  if( !func.is_gga() ) throw std::runtime_error("XC Kernel not GGA!");

  auto base_stack    = data->base_stack;
  func.eval_exc_vxc_device( data->total_npts_task_batch, 
    base_stack.den_eval_device, base_stack.gamma_eval_device, 
    base_stack.eps_eval_device, base_stack.vrho_eval_device, 
    base_stack.vgamma_eval_device, *device_backend->master_stream );

  hadamard_product( *device_backend->master_handle, 
    data->total_npts_task_batch, 1, base_stack.weights_device, 1, 
    base_stack.eps_eval_device, 1 );
  hadamard_product( *device_backend->master_handle, 
    data->total_npts_task_batch, 1, base_stack.weights_device, 1, 
    base_stack.vrho_eval_device, 1 );
  hadamard_product( *device_backend->master_handle, 
    data->total_npts_task_batch, 1, base_stack.weights_device, 1, 
    base_stack.vgamma_eval_device, 1 );

}


void HipAoSScheme1::eval_zmat_lda_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  auto aos_stack     = data->aos_stack;
  zmat_lda_vxc_hip( ntasks, nbe_max, npts_max, aos_stack.device_tasks,
    *device_backend->master_stream );

}

void HipAoSScheme1::eval_zmat_gga_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  auto aos_stack     = data->aos_stack;
  zmat_gga_vxc_hip( ntasks, nbe_max, npts_max, aos_stack.device_tasks,
    *device_backend->master_stream );

}

void HipAoSScheme1::inc_exc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;
  gdot( *device_backend->master_handle, data->total_npts_task_batch,
    base_stack.eps_eval_device, 1, base_stack.den_eval_device, 1, 
    static_stack.acc_scr_device, static_stack.exc_device );

}

void HipAoSScheme1::inc_nel( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;
  gdot( *device_backend->master_handle, data->total_npts_task_batch,
    base_stack.weights_device, 1, base_stack.den_eval_device, 1, 
    static_stack.acc_scr_device, static_stack.nel_device );

}

#if 0
void HipAoSScheme1::inc_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  // Sync blas streams with master stream
  util::hip_event master_event;
  master_event.record( *device_backend->master_stream );
  for( auto& stream : device_backend->blas_streams ) stream.wait( master_event );

  // Launch SYR2K in round robin
  const auto n_blas_streams = device_backend->blas_streams.size();
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];
    syr2k( device_backend->blas_handles[iT % n_blas_streams], 
      HIPBLAS_FILL_MODE_LOWER, HIPBLAS_OP_T, task.nbe, task.npts, 1.,
      task.bf, task.npts, task.zmat, task.npts, 0., task.nbe_scr,
      task.nbe );
  }

  // Record completion of BLAS ops on master stream
  std::vector< util::hip_event > blas_events( n_blas_streams );
  for( size_t iS = 0; iS < n_blas_streams; ++iS )
    blas_events[iS].record( device_backend->blas_streams[iS] );

  for( auto& event : blas_events )
    device_backend->master_stream->wait( event );

  // Increment global VXC
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  task_inc_potential( ntasks, aos_stack.device_tasks, data->vxc_device, nbf,
    submat_block_size, *device_backend->master_stream );
}
#endif

}
