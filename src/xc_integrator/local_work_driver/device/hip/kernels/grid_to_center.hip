#include "hip/hip_runtime.h"
#include <gauxc/util/div_ceil.hpp>
#include "grid_to_center.hpp"
#include "device_specific/hip_device_constants.hpp"

namespace GauXC {

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

}
