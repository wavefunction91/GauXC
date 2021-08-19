#include "cuda_aos_scheme1.hpp"
#include "device/cuda/cuda_backend.hpp"
#include "kernels/grid_to_center.hpp"
#include "kernels/cuda_ssf_2d.hu"
#include "cuda_aos_scheme1_weights.hpp"
 
namespace GauXC {

void cuda_aos_scheme1_weights_wrapper( int32_t npts, int32_t natoms,
  const double* points, const double* RAB, int32_t ldRAB, const double* coords, 
  double* dist, int32_t lddist, const int32_t* iparent,
  const double* dist_nearest, double* weights, cudaStream_t stream ) {

  constexpr auto weight_unroll = detail::CudaAoSScheme1::weight_unroll;
  constexpr auto weight_thread_block = detail::CudaAoSScheme1::weight_thread_block;
  constexpr auto weight_thread_block_per_sm = 
    detail::CudaAoSScheme1::weight_thread_block_per_sm;

  // Compute distances from grid to atomic centers
  compute_grid_to_center_dist( npts, natoms, coords, points, dist, lddist, stream );

  // Get the number of SM's on the device
  int num_sm;
  int dev_id = 0;
  cudaDeviceGetAttribute(&num_sm, cudaDevAttrMultiProcessorCount, dev_id);

  // Modify weights
  dim3 threads( cuda::warp_size, weight_thread_block / cuda::warp_size );
  dim3 blocks ( 1, num_sm * weight_thread_block_per_sm );
  modify_weights_ssf_kernel_2d
    <weight_unroll, weight_thread_block, weight_thread_block_per_sm>
    <<< blocks, threads, 0, stream >>> (
      npts, natoms, RAB, ldRAB, coords, dist, lddist, iparent, dist_nearest,
      weights
    );

}

namespace detail {

 
void CudaAoSScheme1::partition_weights( XCDeviceData* _data ) {
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");


  // Compute distances from grid to atomic centers
  const auto ldatoms = data->get_ldatoms();
  cuda_aos_scheme1_weights_wrapper( data->total_npts_task_batch, data->natoms,
    data->points_device, data->rab_device, ldatoms, data->coords_device, 
    data->dist_scratch_device, ldatoms, data->iparent_device, 
    data->dist_nearest_device, data->weights_device,
    *device_backend->master_stream );
}

}

}
