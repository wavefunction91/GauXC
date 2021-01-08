#include <iostream>

#include <thrust/device_vector.h>

#include "cuda_collision_detection.hpp"


namespace GauXC      {
namespace integrator {
namespace cuda       {


__device__ __inline__ 
int cube_sphere_intersect( 
  const double3 lo, 
  const double3 up,
  const double3 center,
  const double                rad
) {

  double dist = rad * rad;

  if( center.x < lo.x ) {
    const double r_lo = center.x - lo.x;
    const double dist_lo = r_lo * r_lo;
    dist -= dist_lo;
  } else if( center.x > up.x ) {
    const double r_up = center.x - up.x;
    const double dist_up = r_up * r_up;
    dist -= dist_up;
  }

  if( dist < 0. ) return false;

  if( center.y < lo.y ) {
    const double r_lo = center.y - lo.y;
    const double dist_lo = r_lo * r_lo;
    dist -= dist_lo;
  } else if( center.y > up.y ) {
    const double r_up = center.y - up.y;
    const double dist_up = r_up * r_up;
    dist -= dist_up;
  }

  if( dist < 0. ) return false;


  if( center.z < lo.z ) {
    const double r_lo = center.z - lo.z;
    const double dist_lo = r_lo * r_lo;
    dist -= dist_lo;
  } else if( center.z > up.z ) {
    const double r_up = center.z - up.z;
    const double dist_up = r_up * r_up;
    dist -= dist_up;
  }

  return dist > 0.;

}


__global__ void collision_detection_gpu( size_t ncubes,
                          size_t nspheres,
                          size_t LD_bit,
                          const double* low_points,
                          const double* high_points,
                          const double* centers,
                          const double* radii,
                               int32_t* collisions,
                               int32_t* counts) {
  const size_t nspheres_block = (nspheres + 31) / 32;
  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < ncubes; i += blockDim.x * gridDim.x) {
    double3 low_point;
    double3 high_point;
    low_point.x = low_points[3*i+0];
    low_point.y = low_points[3*i+1];
    low_point.z = low_points[3*i+2];

    high_point.x = high_points[3*i+0];
    high_point.y = high_points[3*i+1];
    high_point.z = high_points[3*i+2];


    for (int j_block = 0; j_block < nspheres_block; j_block++) {
      int temp_collisions = 0;
      for (int j_inner = 0; j_inner < 32; j_inner++) {
        int j = j_block * 32 + j_inner;
        if (j < nspheres) {
          double3 center;
          double radius; 
          center.x = centers[3*j+0];
          center.y = centers[3*j+1];
          center.z = centers[3*j+2];

          radius = radii[j];
          temp_collisions |= (cube_sphere_intersect(low_point, high_point, center, radius) ? 1 << (j_inner) : 0);
        }
      }
      collisions[i * LD_bit + j_block] = temp_collisions;
      counts[i] += __popc(temp_collisions);
    }
  }
}

__global__ void bitvector_to_position_list( size_t ncubes, size_t nspheres, size_t LD_bit, const int32_t* collisions, const int32_t* counts, int32_t* position_list) {

  for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < ncubes; i += blockDim.x * gridDim.x) {
    int32_t* out = position_list;
    if (i != 0) {
      out += counts[i-1];
    } 

    int current = 0;
    for (int j = 0; j < nspheres; j++) {
      if( collisions[i * LD_bit + j/32] & (1 << (j%32)) ) {
        out[current++] = j;
      }
    }
  }

}


void collision_detection( size_t ncubes,
                          size_t nspheres,
                          size_t LD_bit,
                          const double* low_points,
                          const double* high_points,
                          const double* centers,
                          const double* radii,
                               int32_t* counts,
                               int32_t** position_list) {

    double *low_points_device, *high_points_device, *centers_device, *radii_device;
    int32_t *collisions_device, *counts_device, *position_list_device;

    // Allocate
    cudaMalloc(&low_points_device, sizeof(double) * 3 * ncubes);
    cudaMalloc(&high_points_device, sizeof(double) * 3 * ncubes);
    cudaMalloc(&centers_device, sizeof(double) * 3 * nspheres);
    cudaMalloc(&radii_device, sizeof(double) * nspheres);
    cudaMalloc(&collisions_device, sizeof(int32_t) * LD_bit * ncubes);

    thrust::device_vector<int32_t> counts_device_vec(ncubes);
    counts_device = thrust::raw_pointer_cast(counts_device_vec.data());

    // Transfer
    cudaMemcpy(low_points_device, low_points, sizeof(double) * 3 * ncubes, cudaMemcpyHostToDevice);
    cudaMemcpy(high_points_device, high_points, sizeof(double) * 3 * ncubes, cudaMemcpyHostToDevice);
    cudaMemcpy(centers_device, centers, sizeof(double) * 3 * nspheres, cudaMemcpyHostToDevice);
    cudaMemcpy(radii_device, radii, sizeof(double) * nspheres, cudaMemcpyHostToDevice);

    // Compute bitvector
    collision_detection_gpu<<<(ncubes + 1023) / 1024, 1024>>>(ncubes, nspheres, LD_bit, low_points_device, high_points_device, centers_device, radii_device, collisions_device, counts_device);
    cudaDeviceSynchronize();

    // Compute scan of count vector and allocate memory
    thrust::inclusive_scan(counts_device_vec.begin(), counts_device_vec.end(), counts_device_vec.begin());
    int32_t total;
    cudaMemcpy(&total, counts_device + ncubes - 1, sizeof(int32_t), cudaMemcpyDeviceToHost);
    cudaMalloc(&position_list_device, sizeof(int32_t) * total);
    (*position_list) = (int32_t*) malloc(sizeof(int32_t) * total);

    // convert from bitvector to position list
    bitvector_to_position_list<<<(ncubes + 1023) / 1024, 1024>>>(ncubes, nspheres, LD_bit, collisions_device, counts_device, position_list_device);
    cudaDeviceSynchronize();

    // Transfer
    cudaMemcpy(*position_list, position_list_device, sizeof(int32_t) * total, cudaMemcpyDeviceToHost);
    cudaMemcpy(counts, counts_device, sizeof(int32_t) * ncubes, cudaMemcpyDeviceToHost);

    // Free
    cudaFree(low_points_device);
    cudaFree(high_points_device);
    cudaFree(centers_device);
    cudaFree(radii_device);
    cudaFree(collisions_device);
    cudaFree(position_list_device);
}

}
}
}

