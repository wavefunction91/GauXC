#include "cuda_aos_scheme1.hpp"
#include "device/cuda/cuda_backend.hpp"
#include "cuda_aos_scheme1_weights.hpp"
//#include "kernels/cublas_extensions.hpp"
#include "device/common/device_blas.hpp"
#include "kernels/uvvars.hpp"
#include "kernels/pack_submat.hpp"
#include "kernels/cuda_inc_potential.hpp"
#include "kernels/symmetrize_mat.hpp"

//#include "device/common/zmat_vxc.hpp"

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
  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;
  auto scheme1_stack = data->scheme1_stack;
  cuda_aos_scheme1_weights_wrapper( data->total_npts_task_batch, data->global_dims.natoms,
    base_stack.points_device, static_stack.rab_device, ldatoms, 
    static_stack.coords_device, scheme1_stack.dist_scratch_device, ldatoms, 
    scheme1_stack.iparent_device, scheme1_stack.dist_nearest_device, 
    base_stack.weights_device, *device_backend->master_stream );
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
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  pack_submat( ntasks, aos_stack.device_tasks, static_stack.dmat_device, 
    nbf, submat_block_size, *device_backend->master_stream );


  // Sync blas streams with master stream
  data->device_backend_->sync_blas_pool_with_master();

  // Launch GEMM in round-robin
  const auto n_blas_streams = data->device_backend_->blas_pool_size();
  size_t nsingle = 0;
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];
    if(task.ncut > 1)
      gemm( data->device_backend_->blas_pool_handle(iT % n_blas_streams), 
        DeviceBlasOp::NoTrans, DeviceBlasOp::NoTrans, 
        task.npts, task.nbe, task.nbe, 1., task.bf, task.npts, task.nbe_scr, 
        task.nbe, 0., task.zmat, task.npts );
    else {
      gemm( data->device_backend_->blas_pool_handle(iT % n_blas_streams), 
        DeviceBlasOp::NoTrans, DeviceBlasOp::NoTrans,
        task.npts, task.nbe, task.nbe, 1., task.bf, task.npts,
        static_stack.dmat_device + task.ibf_begin*(nbf+1), nbf, 
        0., task.zmat, task.npts );
      nsingle++;
    }
  }

  // Record completion of BLAS ops on master stream
  data->device_backend_->sync_master_with_blas_pool();

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
  auto base_stack    = data->base_stack;
  util::cuda_set_zero_async( data->total_npts_task_batch, 
    base_stack.den_eval_device, *device_backend->master_stream, "Den Zero" );

  // Evaluate U variables
  auto aos_stack     = data->aos_stack;
  eval_uvvars_lda_cuda( ntasks, nbe_max, npts_max, 
    aos_stack.device_tasks, *device_backend->master_stream );

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
  auto base_stack    = data->base_stack;
  util::cuda_set_zero_async( data->total_npts_task_batch, base_stack.den_eval_device,
    *device_backend->master_stream, "Den Zero" );
  util::cuda_set_zero_async( data->total_npts_task_batch, base_stack.den_x_eval_device,
    *device_backend->master_stream, "DenX Zero" );
  util::cuda_set_zero_async( data->total_npts_task_batch, base_stack.den_y_eval_device,
    *device_backend->master_stream, "DenY Zero" );
  util::cuda_set_zero_async( data->total_npts_task_batch, base_stack.den_z_eval_device,
    *device_backend->master_stream, "DenZ Zero" );

  // Evaluate U variables
  auto aos_stack     = data->aos_stack;
  eval_uvvars_gga_cuda( ntasks, data->total_npts_task_batch, nbe_max, npts_max, 
    aos_stack.device_tasks, base_stack.den_x_eval_device, base_stack.den_y_eval_device,
    base_stack.den_z_eval_device, base_stack.gamma_eval_device, 
    *device_backend->master_stream );

}

void CudaAoSScheme1::eval_kern_exc_vxc_lda( const functional_type& func, 
  XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  if( !func.is_lda() ) throw std::runtime_error("XC Kernel not LDA!");

  auto base_stack    = data->base_stack;
  func.eval_exc_vxc_device( data->total_npts_task_batch, base_stack.den_eval_device,
    base_stack.eps_eval_device, base_stack.vrho_eval_device, 
    *device_backend->master_stream );

  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.eps_eval_device, 1 );
  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.vrho_eval_device, 1 );

}


void CudaAoSScheme1::eval_kern_exc_vxc_gga( const functional_type& func, 
  XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  if( !func.is_gga() ) throw std::runtime_error("XC Kernel not GGA!");

  auto base_stack    = data->base_stack;
  func.eval_exc_vxc_device( data->total_npts_task_batch, base_stack.den_eval_device,
    base_stack.gamma_eval_device, base_stack.eps_eval_device, base_stack.vrho_eval_device, 
    base_stack.vgamma_eval_device, *device_backend->master_stream );

  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.eps_eval_device, 1 );
  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.vrho_eval_device, 1 );
  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.vgamma_eval_device, 1 );

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
    syr2k( data->device_backend_->blas_pool_handle(iT % n_blas_streams), 
      DeviceBlasUplo::Lower, DeviceBlasOp::Trans, task.nbe, task.npts, 1.,
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
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  task_inc_potential( ntasks, aos_stack.device_tasks, 
    static_stack.vxc_device, nbf, submat_block_size, 
    *device_backend->master_stream );
}

void CudaAoSScheme1::symmetrize_vxc( XCDeviceData* _data) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  const auto nbf = data->global_dims.nbf;
  auto static_stack  = data->static_stack;
  symmetrize_matrix( nbf, static_stack.vxc_device, nbf, 
    *device_backend->master_stream ); 

}

}
