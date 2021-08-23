#include "hip/hip_runtime.h"
#include <gauxc/util/div_ceil.hpp>
#include "grid_to_center.hpp"
#include "device_specific/hip_device_constants.hpp"

namespace GauXC {

#if 0
__global__ void compute_grid_to_center_dist(
        int32_t npts,
        int32_t natoms,
  const double* coords,
  const double* points,
        double* dist,
        int32_t lddist
) {

  __shared__ double3 point_buffer[hip::warp_size];
  double3 coord_reg;

  const int natoms_block = (natoms + hip::warp_size-1) / hip::warp_size;
  const int coords_block = (npts + hip::warp_size-1) / hip::warp_size;

  const double3* coords_vec = (double3*) coords;
  const double3* points_vec = (double3*) points;

  for (int j = blockIdx.x; j < natoms_block; j += gridDim.x) {
    const int iAtom = j * hip::warp_size + threadIdx.x;
    // Load blocks into registers/shared memory
    if (iAtom < natoms) {
      coord_reg = coords_vec[iAtom];
    }
    for (int i = blockIdx.y; i < coords_block; i += gridDim.y) {
      const int iPt_load = i * hip::warp_size + threadIdx.x;
      if (iPt_load < npts) {
        point_buffer[threadIdx.x] = points_vec[iPt_load];
      }
      __syncthreads();

      // do the computation
      #pragma unroll 2
      for (int k = threadIdx.y; k < hip::warp_size; k+=hip::warp_size/2) {
        const int iPt_sm = k;
        const int iPt = i * hip::warp_size + iPt_sm;
        const double rx = point_buffer[iPt_sm].x - coord_reg.x;
        const double ry = point_buffer[iPt_sm].y - coord_reg.y;
        const double rz = point_buffer[iPt_sm].z - coord_reg.z;

        if (iAtom < natoms and iPt < npts) {
          dist[ iAtom + iPt * lddist ] = std::sqrt( rx*rx + ry*ry + rz*rz );
        }
      }
      __syncthreads();
    }
  }
}

void compute_grid_to_center_dist( int32_t npts, int32_t natoms,
  const double* coords, const double* points,  double* dist,
  int32_t lddist, hipStream_t stream ) {

    const int distance_thread_y = hip::max_warps_per_thread_block / 2;
    dim3 threads( hip::warp_size, distance_thread_y );
    dim3 blocks( util::div_ceil( natoms,   threads.x), 
                 util::div_ceil( npts, threads.y * distance_thread_y) );

    hipLaunchKernelGGL(compute_grid_to_center_dist, dim3(blocks), dim3(threads), 0, stream, 
      npts, natoms, coords, points, dist, lddist
    );

}

#else

__global__ void compute_grid_to_center_dist(
        int32_t npts,
        int32_t natoms,
  const double* coords,
  const double* points,
        double* dist,
        int32_t lddist
) {

  const auto tid_x = threadIdx.x + blockIdx.x*blockDim.x;
  const auto tid_y = threadIdx.y + blockIdx.y*blockDim.y;

  if( tid_x < natoms && tid_y < npts ) {
    const double3* coords_vec = (double3*) coords;
    const double3* points_vec = (double3*) points;
   
    auto RA = coords_vec[tid_x];
    auto ri = points_vec[tid_y];

    const auto rx = RA.x - ri.x;
    const auto ry = RA.y - ri.y;
    const auto rz = RA.z - ri.z;

    dist[ tid_x + tid_y * lddist ] = std::sqrt(rx*rx + ry*ry + rz*rz);
  }

}

void compute_grid_to_center_dist( int32_t npts, int32_t natoms,
  const double* coords, const double* points,  double* dist,
  int32_t lddist, hipStream_t stream ) {


  dim3 threads( hip::warp_size, hip::max_warps_per_thread_block );
  dim3 blocks( util::div_ceil( natoms, threads.x ),
               util::div_ceil( npts, threads.y ) );

  hipLaunchKernelGGL( compute_grid_to_center_dist, blocks, threads, 0, stream,
    npts, natoms, coords, points, dist, lddist );

}
#endif

}
