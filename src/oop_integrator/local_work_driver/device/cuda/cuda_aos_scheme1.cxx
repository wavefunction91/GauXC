#include "cuda_aos_scheme1.hpp"
#include "kernels/collocation_device.hpp"
#include "device/cuda/cuda_backend.hpp"
#include "cuda_aos_scheme1_weights.hpp"
#include "kernels/cublas_extensions.hpp"
#include "kernels/uvvars.hpp"
#include "kernels/zmat_vxc.hpp"
#include "kernels/pack_submat.hpp"
#include "kernels/cuda_inc_potential.hpp"

namespace GauXC {

std::unique_ptr<XCDeviceData> CudaAoSScheme1::create_device_data() {
  return std::make_unique<Data>();
}

 
void CudaAoSScheme1::partition_weights( XCDeviceData* _data ) {
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");


  // Compute distances from grid to atomic centers
  const auto ldatoms = data->get_ldatoms();
  cuda_aos_scheme1_weights_wrapper( data->total_npts_task_batch, data->global_dims.natoms,
    data->points_device, data->rab_device, ldatoms, data->coords_device, 
    data->dist_scratch_device, ldatoms, data->iparent_device, 
    data->dist_nearest_device, data->weights_device,
    *device_backend->master_stream );
}


void CudaAoSScheme1::eval_collocation( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  size_t npts_max = 0, nshells_max = 0;
  for( auto& task : tasks ) {
    npts_max    = std::max( npts_max, task.npts );
    nshells_max = std::max( nshells_max, task.nshells );
  }

  eval_collocation_masked_combined( ntasks, npts_max, nshells_max,
    data->shells_device, data->device_tasks, *device_backend->master_stream );

}

void CudaAoSScheme1::eval_collocation_gradient( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  size_t npts_max = 0, nshells_max = 0;
  for( auto& task : tasks ) {
    npts_max    = std::max( npts_max, task.npts );
    nshells_max = std::max( nshells_max, task.nshells );
  }

  eval_collocation_masked_combined_deriv1( ntasks, npts_max, nshells_max,
    data->shells_device, data->device_tasks, *device_backend->master_stream );
  
}

void CudaAoSScheme1::eval_xmat( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  // Pack density matrix 
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  pack_submat( ntasks, data->device_tasks, data->dmat_device, nbf, 
    submat_block_size, *device_backend->master_stream );


  // Sync blas streams with master stream
  util::cuda_event master_event;
  master_event.record( *device_backend->master_stream );
  for( auto& stream : device_backend->blas_streams ) stream.wait( master_event );

  // Launch GEMM in round-robin
  const auto n_blas_streams = device_backend->blas_streams.size();
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];
    gemm( device_backend->blas_handles[iT % n_blas_streams], CUBLAS_OP_N,
      CUBLAS_OP_N, task.npts, task.nbe, task.nbe, 1., task.bf, task.npts,
      task.nbe_scr, task.nbe, 0., task.zmat, task.npts );
  }

  // Record completion of BLAS ops on master stream
  std::vector< util::cuda_event > blas_events( n_blas_streams );
  for( size_t iS = 0; iS < n_blas_streams; ++iS )
    blas_events[iS].record( device_backend->blas_streams[iS] );

  for( auto& event : blas_events )
    device_backend->master_stream->wait( event );

}






void CudaAoSScheme1::eval_uvvar_lda( XCDeviceData* _data ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  // Zero density
  util::cuda_set_zero_async( data->total_npts_task_batch, data->den_eval_device,
    *device_backend->master_stream, "Den Zero" );

  // Evaluate U variables
  eval_uvvars_lda_cuda( ntasks, nbe_max, npts_max, data->device_tasks,
    *device_backend->master_stream );

}




void CudaAoSScheme1::eval_uvvar_gga( XCDeviceData* _data ){


  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  // Zero density + gradient
  util::cuda_set_zero_async( data->total_npts_task_batch, data->den_eval_device,
    *device_backend->master_stream, "Den Zero" );
  util::cuda_set_zero_async( data->total_npts_task_batch, data->den_x_eval_device,
    *device_backend->master_stream, "DenX Zero" );
  util::cuda_set_zero_async( data->total_npts_task_batch, data->den_y_eval_device,
    *device_backend->master_stream, "DenY Zero" );
  util::cuda_set_zero_async( data->total_npts_task_batch, data->den_z_eval_device,
    *device_backend->master_stream, "DenZ Zero" );

  // Evaluate U variables
  eval_uvvars_gga_cuda( ntasks, data->total_npts_task_batch, nbe_max, npts_max, 
    data->device_tasks, data->den_x_eval_device, data->den_y_eval_device,
    data->den_z_eval_device, data->gamma_eval_device, 
    *device_backend->master_stream );

}

void CudaAoSScheme1::eval_kern_exc_vxc_lda( const functional_type& func, 
  XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  if( !func.is_lda() ) throw std::runtime_error("XC Kernel not LDA!");

  func.eval_exc_vxc_device( data->total_npts_task_batch, data->den_eval_device,
    data->eps_eval_device, data->vrho_eval_device, *device_backend->master_stream );

  hadamard_product( *device_backend->master_handle, data->total_npts_task_batch, 1, 
                    data->weights_device, 1, data->eps_eval_device, 1 );
  hadamard_product( *device_backend->master_handle, data->total_npts_task_batch, 1, 
                    data->weights_device, 1, data->vrho_eval_device, 1 );

}


void CudaAoSScheme1::eval_kern_exc_vxc_gga( const functional_type& func, 
  XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  if( !func.is_gga() ) throw std::runtime_error("XC Kernel not GGA!");

  func.eval_exc_vxc_device( data->total_npts_task_batch, data->den_eval_device,
    data->gamma_eval_device, data->eps_eval_device, data->vrho_eval_device, 
    data->vgamma_eval_device, *device_backend->master_stream );

  hadamard_product( *device_backend->master_handle, data->total_npts_task_batch, 1, 
                    data->weights_device, 1, data->eps_eval_device, 1 );
  hadamard_product( *device_backend->master_handle, data->total_npts_task_batch, 1, 
                    data->weights_device, 1, data->vrho_eval_device, 1 );
  hadamard_product( *device_backend->master_handle, data->total_npts_task_batch, 1, 
                    data->weights_device, 1, data->vgamma_eval_device, 1 );

}




void CudaAoSScheme1::eval_zmat_lda_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  zmat_lda_vxc_cuda( ntasks, nbe_max, npts_max, data->device_tasks,
    *device_backend->master_stream );

}

void CudaAoSScheme1::eval_zmat_gga_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  zmat_gga_vxc_cuda( ntasks, nbe_max, npts_max, data->device_tasks,
    *device_backend->master_stream );

}
void CudaAoSScheme1::inc_exc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  gdot( *device_backend->master_handle, data->total_npts_task_batch,
    data->eps_eval_device, 1, data->den_eval_device, 1, data->acc_scr_device,
    data->exc_device );

}
void CudaAoSScheme1::inc_nel( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  gdot( *device_backend->master_handle, data->total_npts_task_batch,
    data->weights_device, 1, data->den_eval_device, 1, data->acc_scr_device,
    data->nel_device );

}



void CudaAoSScheme1::inc_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  // Sync blas streams with master stream
  util::cuda_event master_event;
  master_event.record( *device_backend->master_stream );
  for( auto& stream : device_backend->blas_streams ) stream.wait( master_event );

  // Launch SYR2K in round robin
  const auto n_blas_streams = device_backend->blas_streams.size();
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];
    syr2k( device_backend->blas_handles[iT % n_blas_streams], 
      CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, task.nbe, task.npts, 1.,
      task.bf, task.npts, task.zmat, task.npts, 0., task.nbe_scr,
      task.nbe );
  }

  // Record completion of BLAS ops on master stream
  std::vector< util::cuda_event > blas_events( n_blas_streams );
  for( size_t iS = 0; iS < n_blas_streams; ++iS )
    blas_events[iS].record( device_backend->blas_streams[iS] );

  for( auto& event : blas_events )
    device_backend->master_stream->wait( event );

  // Increment global VXC
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  task_inc_potential( ntasks, data->device_tasks, data->vxc_device, nbf,
    submat_block_size, *device_backend->master_stream );
}

}
