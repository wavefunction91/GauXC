#include "scheme1_magma_base.hpp"
#include "device/common/pack_submat.hpp"
#include "device/common/inc_potential.hpp"
#include "device/common/device_blas.hpp"
#include "device/cuda/cuda_backend.hpp"

#include "device_specific/magma_util.hpp"

namespace GauXC {

void AoSScheme1MAGMABase::eval_xmat( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  if( not data->device_backend_ ) throw std::runtime_error("INVALID DEVICE BACKEND");

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  // Pack density matrix 
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  pack_submat( ntasks, aos_stack.device_tasks, static_stack.dmat_device, 
    nbf, submat_block_size, data->device_backend_->queue() );


  auto cuda_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  util::magma_queue master_queue( 0, *cuda_backend->master_stream, *cuda_backend->master_handle );

  auto magma_stack = data->magma_stack;
  magmablas_dgemm_vbatched( MagmaNoTrans, MagmaNoTrans,
    magma_stack.m_array_device, magma_stack.n_array_device, 
    magma_stack.k_array_device, 
    1., magma_stack.bf_array_device,   magma_stack.ld_bf_array_device,
        magma_stack.dmat_array_device, magma_stack.ld_dmat_array_device,
    0., magma_stack.zmat_array_device, magma_stack.ld_zmat_array_device,
    ntasks, master_queue );
  data->device_backend_->master_queue_synchronize(); 

}

void AoSScheme1MAGMABase::inc_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  if( not data->device_backend_ ) throw std::runtime_error("INVALID DEVICE BACKEND");

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  auto cuda_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  util::magma_queue master_queue( 0, *cuda_backend->master_stream, *cuda_backend->master_handle );

  auto magma_stack = data->magma_stack;
  magmablas_dsyr2k_vbatched( MagmaLower, MagmaTrans, 
    magma_stack.n_array_device, magma_stack.m_array_device,
    1., magma_stack.bf_array_device,   magma_stack.ld_bf_array_device, 
        magma_stack.zmat_array_device, magma_stack.ld_zmat_array_device,
    0., magma_stack.vmat_array_device, magma_stack.ld_vmat_array_device, 
    ntasks, master_queue );
  data->device_backend_->master_queue_synchronize(); 

  // Increment global VXC
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  task_inc_potential( ntasks, aos_stack.device_tasks, 
    static_stack.vxc_device, nbf, submat_block_size, 
    data->device_backend_->queue() );
}

}
