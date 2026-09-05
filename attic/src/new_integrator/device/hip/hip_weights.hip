#include "hip/hip_runtime.h"
#include <gauxc/util/div_ceil.hpp>

#include "device/hip/hip_weights.hpp"
#include "common/integrator_constants.hpp"
#include "device/hip/hip_extensions.hpp"
#include "device/hip/hip_device_properties.hpp"

constexpr double eps_d = std::numeric_limits<double>::epsilon();


namespace GauXC      {
namespace integrator {
namespace hip       {

using namespace GauXC::hip;

__global__ void reciprocal_kernel(size_t length, double* vec) {
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < length; i += blockDim.x * gridDim.x) {
     vec[i] = 1. / vec[i];
   }
}

__global__ void compute_point_center_dist(
        size_t      npts,
        size_t      LDatoms,
        size_t      natoms,
  const double*     coords,
  const double*     points,
        double*     dist
) {

  __shared__ double3 point_buffer[warp_size];
  double3 coord_reg;

  const int natoms_block = (natoms + warp_size-1) / warp_size;
  const int coords_block = (npts + warp_size-1) / warp_size;

  const double3* coords_vec = (double3*) coords;
  const double3* points_vec = (double3*) points;

  for (int j = blockIdx.x; j < natoms_block; j += gridDim.x) {
    const int iAtom = j * warp_size + threadIdx.x;
    // Load blocks into registers/shared memory
    if (iAtom < natoms) {
      coord_reg = coords_vec[iAtom];
    }
    for (int i = blockIdx.y; i < coords_block; i += gridDim.y) {
      const int iPt_load = i * warp_size + threadIdx.x;
      if (iPt_load < npts) {
        point_buffer[threadIdx.x] = points_vec[iPt_load];
      }
      __syncthreads();

      // do the computation
      #pragma unroll 2
      for (int k = threadIdx.y; k < warp_size; k+=warp_size/2) {
        const int iPt_sm = k;
        const int iPt = i * warp_size + iPt_sm;
        const double rx = point_buffer[iPt_sm].x - coord_reg.x;
        const double ry = point_buffer[iPt_sm].y - coord_reg.y;
        const double rz = point_buffer[iPt_sm].z - coord_reg.z;

        if (iAtom < natoms and iPt < npts) {
          dist[ iAtom + iPt * LDatoms ] = std::sqrt( rx*rx + ry*ry + rz*rz );
        }
      }
      __syncthreads();
    }
  }
}

#if 0
__global__ void modify_weights_becke_kernel(
        size_t                            npts,
        size_t                            natoms,
  const double*                           RAB,
  const double*                           coords,
  const double*                           dist_scratch,
  const int32_t*                          iparent_device,
        double*                           weights_device
) {

  // Becke partition functions
  auto hBecke = [](double x) {return 1.5 * x - 0.5 * x * x * x;}; // Eq. 19
  auto gBecke = [&](double x) {return hBecke(hBecke(hBecke(x)));}; // Eq. 20 f_3


  __shared__ double shared[2048]; 
  for( int ipt = blockIdx.x; ipt < npts; ipt += gridDim.x ) {

    const auto iParent = iparent_device[ipt];

    double sum = 0.; 
    double parent_weight = 0.;

    const double* const local_dist_scratch = dist_scratch + ipt * natoms;
    for( int iCenter = threadIdx.y; iCenter < natoms; iCenter += blockDim.y ) {

      const double ri = local_dist_scratch[ iCenter ];

      const double* const local_rab = RAB + iCenter * natoms;

      double ps = 1.;
      for( int jCenter = threadIdx.x; jCenter < natoms; jCenter += blockDim.x ) {
      
        const double rj = local_dist_scratch[ jCenter ];

        const double mu = (ri - rj) / local_rab[ jCenter ]; // XXX: RAB is symmetric
        const double s  = 0.5 * ( 1. - gBecke( mu ) );

        ps *= (iCenter == jCenter) ? 1. : s ;

      }

      ps = warp_prod_reduce( ps ); // XXX: Assumes blockDim.x == 32

      if( iCenter == iParent ) parent_weight = ps;

      sum += ps;

    }

    // XXX: Assumes blockDim.x == blockDim.y == 32
    if( threadIdx.x == 0 ) {
      shared[ threadIdx.y ]        = sum;
      shared[ threadIdx.y + 1024]  = parent_weight;
    }

    __syncthreads();
    sum = shared[ threadIdx.x ];
    sum = warpReduceSum( sum );

    __syncthreads();
    parent_weight = shared[ threadIdx.x + 1024];
    parent_weight = __shfl_sync(0xffffffff, parent_weight, iParent % 32, 32 );

    if( threadIdx.x == 0 and threadIdx.y == 0 )
      weights_device[ipt] *= parent_weight / sum;
    

  }


}



__global__ void modify_weights_ssf_kernel(
        size_t                            npts,
        size_t                            natoms,
  const double*                           RAB,
  const double*                           coords,
  const double*                           dist_scratch,
  const int32_t*                          iparent_device,
  const double*                           dist_nearest_device,
        double*                           weights_device
) {

  // Frisch partition functions
  auto gFrisch = [](double x) {

    const double s_x  = x / magic_ssf_factor<>;
    const double s_x2 = s_x  * s_x;
    const double s_x3 = s_x  * s_x2;
    const double s_x5 = s_x3 * s_x2;
    const double s_x7 = s_x5 * s_x2;

    return (35.*(s_x - s_x3) + 21.*s_x5 - 5.*s_x7) / 16.;
  };
  
  auto sFrisch = [&] (double x) {
    const double g = 0.5 * (1. - gFrisch(x));
    return (x >= magic_ssf_factor<>) ? 0. : (x <= -magic_ssf_factor<>) ? 1. : g;
  };

  constexpr double weight_tol = 1e-10;

  __shared__ double shared[2048]; 
  for( int ipt = blockIdx.x; ipt < npts; ipt += gridDim.x ) {

    const auto iParent = iparent_device[ipt];

    double sum = 0.; 
    double parent_weight = 0.;

    const double* const local_dist_scratch = dist_scratch + ipt * natoms;
    const double dist_cutoff = 0.5 * (1 - magic_ssf_factor<> ) * 
      dist_nearest_device[ipt];
    if( local_dist_scratch[iParent] < dist_cutoff ) continue;

    for( int iCenter = threadIdx.y; iCenter < natoms; iCenter += blockDim.y ) {

      const double ri = local_dist_scratch[ iCenter ];

      const double* const local_rab = RAB + iCenter * natoms;

      double ps = 1.;
      for( int jCenter = threadIdx.x; jCenter < natoms; jCenter += blockDim.x ) 
      if( fabs(ps) > weight_tol ) {
      
        const double rj = local_dist_scratch[ jCenter ];

        const double mu = (ri - rj) / local_rab[ jCenter ]; // XXX: RAB is symmetric
        const double s  = sFrisch( mu );
        ps *= (iCenter == jCenter) ? 1. : s ;

      }

      ps = warp_prod_reduce( ps ); // XXX: Assumes blockDim.x == 32

      if( iCenter == iParent ) parent_weight = ps;

      sum += ps;

    }

    // XXX: Assumes blockDim.x == blockDim.y == 32
    if( threadIdx.x == 0 ) {
      shared[ threadIdx.y ]        = sum;
      shared[ threadIdx.y + 1024]  = parent_weight;
    }

    __syncthreads();
    sum = shared[ threadIdx.x ];
    sum = warpReduceSum( sum );

    __syncthreads();
    parent_weight = shared[ threadIdx.x + 1024];
    parent_weight = __shfl_sync(0xffffffff, parent_weight, iParent % 32, 32 );

    if( threadIdx.x == 0 and threadIdx.y == 0 )
      weights_device[ipt] *= parent_weight / sum;
    

  }


}
#endif

// SIMT over points: 1D kernel
__global__ void modify_weights_ssf_kernel_1d(
        size_t                            npts,
        size_t                            natoms,
  const double*                           RAB,
  const double*                           coords,
  const double*                           dist_scratch,
  const int32_t*                          iparent_device,
  const double*                           dist_nearest_device,
        double*                           weights_device
) {

  // Frisch partition functions
  auto gFrisch = [](double x) {

    const double s_x  = x / magic_ssf_factor<>;
    const double s_x2 = s_x  * s_x;
    const double s_x3 = s_x  * s_x2;
    const double s_x5 = s_x3 * s_x2;
    const double s_x7 = s_x5 * s_x2;

    return (35.*(s_x - s_x3) + 21.*s_x5 - 5.*s_x7) / 16.;
  };
  
#if 0
  auto sFrisch = [&] (double x) {
    const double g = 0.5 * (1. - gFrisch(x));
    return (x >= magic_ssf_factor<>) ? 0. : (x <= -magic_ssf_factor<>) ? 1. : g;
  };
#else
  auto sFrisch = [&] (double x) {
    if( fabs(x) < magic_ssf_factor<> ) return 0.5 * (1. - gFrisch(x));
    else if( x >= magic_ssf_factor<> ) return 0.;
    else                               return 1.;
  };
#endif

  constexpr double weight_tol = 1e-10;

  const int tid_x = threadIdx.x + blockIdx.x * blockDim.x;
  const int nt_x  = blockDim.x  * gridDim.x;

  //__shared__ double shared[2048]; 
  for( int ipt = tid_x; ipt < npts; ipt += nt_x ) {

    const auto iParent = iparent_device[ipt];

    double sum = 0.; 
    double parent_weight = 0.;

    const double* const local_dist_scratch = dist_scratch + ipt * natoms;
    const double dist_cutoff = 0.5 * (1 - magic_ssf_factor<> ) * 
      dist_nearest_device[ipt];
    if( local_dist_scratch[iParent] < dist_cutoff ) continue;

#if 0
    for( int iCenter = 0; iCenter < natoms; iCenter++ ) {

      const double ri = local_dist_scratch[ iCenter ];

      const double* const local_rab = RAB + iCenter * natoms;

      double ps = 1.;
      for( int jCenter = 0; jCenter < natoms; jCenter++ ) 
      if( fabs(ps) > weight_tol ) {
      if( iCenter != jCenter ) {
      
        const double rj = local_dist_scratch[ jCenter ];

        const double mu = (ri - rj) / local_rab[ jCenter ]; // XXX: RAB is symmetric
        ps *= sFrisch( mu );

      }
      } else break;

      //__syncwarp();

      if( iCenter == iParent ) parent_weight = ps;

      sum += ps;

    }
#else

    // Do iParent First
    {

      const double ri = local_dist_scratch[ iParent ];
      const double* const local_rab = RAB + iParent * natoms;

      parent_weight = 1.;
      for( int jCenter = 0; jCenter < natoms; jCenter++ ) 
      if( parent_weight > weight_tol ) {
      if( iParent != jCenter ) {
      
        const double rj = local_dist_scratch[ jCenter ];

        const double mu = (ri - rj) / local_rab[ jCenter ]; // XXX: RAB is symmetric
        parent_weight *= sFrisch( mu );

      }
      } else break;

      //__syncwarp();
      sum += parent_weight;

    }

    if( parent_weight < eps_d ) {
      weights_device[ipt] = 0.;
      continue;
    }

    for( int iCenter = 0; iCenter < natoms; iCenter++ ) 
    if( iParent != iCenter ) {

      const double ri = local_dist_scratch[ iCenter ];

      const double* const local_rab = RAB + iCenter * natoms;

      double ps = 1.;
      for( int jCenter = 0; jCenter < natoms; jCenter++ ) 
      if( ps > weight_tol ) {
      if( iCenter != jCenter ) {
      
        const double rj = local_dist_scratch[ jCenter ];

        const double mu = (ri - rj) / local_rab[ jCenter ]; // XXX: RAB is symmetric
        ps *= sFrisch( mu );

      }
      } else break;

      //__syncwarp();
      sum += ps;

    }

#endif

    weights_device[ipt] *= parent_weight / sum;
    

  }


}

__device__ __inline__ double gFrisch(double x) {
  // Frisch partition functions
//  const double s_x  = x / magic_ssf_factor<>;
  const double s_x  = x * 1.5625;
  const double s_x2 = s_x  * s_x;
  const double s_x3 = s_x  * s_x2;
  const double s_x5 = s_x3 * s_x2;
  const double s_x7 = s_x5 * s_x2;

  return ((35.) *(s_x - s_x3) + (21.) *s_x5 - (5.) *s_x7);
}


__device__ __inline__ double sFrisch(double x) {
    //double frisch_val = (0.5 - (0.5/ 16.0) * gFrisch(x));

    if( fabs(x) < magic_ssf_factor<> ) return (0.5 - (0.5/ 16.0) * gFrisch(x));
    else if( x >= magic_ssf_factor<> ) return 0.;
    else                               return 1.;
}

__global__ __launch_bounds__(weight_thread_block, weight_thread_block_per_sm)
void modify_weights_ssf_kernel_2d(
        size_t                            npts,
        size_t                            LDatoms,
        size_t                            natoms,
  const double*                           RAB,
  const double*                           coords,
  const double*                           dist_scratch,
  const int32_t*                          iparent_device,
  const double*                           dist_nearest_device,
        double*                           weights_device
) {
  constexpr double weight_tol = 1e-10;
  int natom_block = ((natoms + blockDim.x - 1) / blockDim.x) * blockDim.x;

  const int tid_x = threadIdx.y + blockIdx.y * blockDim.y;
  const int nt_x  = blockDim.y  * gridDim.y;

  __shared__ int jCounter_sm[max_warps_per_thread_block];
  int* jCounter = reinterpret_cast<int *>(jCounter_sm) + threadIdx.y;

  // Each warp will work together on a point
  for( int ipt = tid_x; ipt < npts; ipt += nt_x ) {

    const auto iParent = iparent_device[ipt];

    double sum = 0.; 
    double parent_weight = 0.;

    const double* const local_dist_scratch = dist_scratch + ipt * LDatoms;
    const double dist_cutoff = 0.5 * (1 - magic_ssf_factor<> ) * 
      dist_nearest_device[ipt];
    if( local_dist_scratch[iParent] < dist_cutoff ) continue;

    // Do iParent First
    {

      const double ri = local_dist_scratch[ iParent ];
      const double* const local_rab = RAB + iParent * LDatoms;

      parent_weight = 1.;
      for( int jCenter = threadIdx.x; jCenter < natom_block; jCenter+=blockDim.x ) {
        double contribution = 1.0;
        if (jCenter < natoms && iParent != jCenter) {
          const double rj = local_dist_scratch[ jCenter ];
          const double mu = (ri - rj) * local_rab[ jCenter ]; // XXX: RAB is symmetric
          contribution = sFrisch( mu );
        }
        contribution = warpReduceProd(contribution);
        parent_weight *= contribution;

        if (parent_weight < weight_tol) break;
      }
    }

    if( parent_weight < eps_d ) {
      if (threadIdx.x == 0)
        weights_device[ipt] = 0.;
      __syncwarp();
      continue;
    }

    // Initialize each counter to 0
    if (threadIdx.x == 0) {
      jCounter[0] = 0;
    }
    __syncwarp();

    // Each thread will process an iCenter. Atomic operations are used to assign
    // an iCenter value to each thread.
    int iCenter = atomicAdd(jCounter, 1);
    if (iCenter >= iParent) iCenter++; // iCenter == iParent is skipped

    // The entire warp processes the same jCenter value at the same time
    int jCenter = 0;

    const double* local_rab = RAB + iCenter * LDatoms;
    double ri = local_dist_scratch[ iCenter ];
    double ps = 1.;
    int iCount = 0; 
    int cont = (iCenter < natoms);

    // We will continue iterating until all of the threads have cont set to 0
    while (__any_sync(0xffffffff, cont)) {
      if (cont) {
        double2 rj[weight_unroll/2];
        double2 rab_val[weight_unroll/2];
        double mu[weight_unroll];
        iCount += weight_unroll;

        #pragma unroll
        for (int k = 0; k < weight_unroll/2; k++) {
          rj[k]      = *((double2*)(local_dist_scratch + jCenter) + k);
          rab_val[k] = *((double2*)(local_rab          + jCenter) + k); 
        }

        #pragma unroll
        for (int k = 0; k < weight_unroll/2; k++) {
          mu[2*k+0] = (ri - rj[k].x) * rab_val[k].x; // XXX: RAB is symmetric
          mu[2*k+1] = (ri - rj[k].y) * rab_val[k].y; 
        }

        #pragma unroll
        for (int k = 0; k < weight_unroll; k++) {
          if((iCenter != jCenter + k) && (jCenter + k < natoms)) {
            mu[k] = sFrisch( mu[k] );
            ps *= mu[k];
          }
        }

        // A thread is done with a iCenter based on 2 conditions. Weight tolerance
        // Or if it has seen all of the jCenters
        if( !(ps > weight_tol && iCount < LDatoms )) {
          // In the case were the thread is done, it begins processing another iCenter
          sum += ps;
          iCenter = atomicAdd(jCounter, 1);
          if (iCenter >= iParent) iCenter++;

          // If there are no more iCenters left to process, it signals it is ready to exit
          cont = (iCenter < natoms);
          ri = local_dist_scratch[ iCenter ];
          local_rab = RAB + iCenter * LDatoms;
          ps = 1.;
          iCount = 0;
        }
      }
      // Wraps jCenter around. This was faster than modulo
      jCenter += weight_unroll;
      jCenter = (jCenter < LDatoms) ? jCenter : 0;
    }

    // All of the threads then sum their contributions. Only thread 0 needs to add the parent
    // contribution.
    __syncwarp();
    sum = warpReduceSum(sum);
    if (threadIdx.x == 0) {
      sum += parent_weight;
      weights_device[ipt] *= parent_weight / sum;
    }

    __syncwarp();

  }
}


void hip_reciprocal(size_t length, double* vec, hipStream_t stream) {
  dim3 threads(max_threads_per_thread_block);
  dim3 blocks( get_device_sm_count(0) ); 
  hipLaunchKernelGGL(reciprocal_kernel, dim3(threads), dim3(blocks), 0, stream, length, vec);
}


template <typename F>
void partition_weights_hip_SoA( XCWeightAlg    weight_alg,
                                 size_t         npts,
                                 size_t         LDatoms,
                                 size_t         natoms,
                                 const F*       points_device,
                                 const int32_t* iparent_device,
                                 const F*       dist_nearest_device,
                                 const F*       rab_device,
                                 const F*       atomic_coords_device,
                                       F*       weights_device,
                                       F*       dist_scratch_device,
                                 hipStream_t   stream ) {



  // Evaluate point-to-atom collocation
  {
    const int distance_thread_y = max_warps_per_thread_block / 2;
    dim3 threads(  warp_size, distance_thread_y );
    dim3 blocks( util::div_ceil( natoms,   threads.x), 
                 util::div_ceil( npts, threads.y * distance_thread_y) );

    hipLaunchKernelGGL(compute_point_center_dist, dim3(blocks), dim3(threads), 0, stream, 
      npts, LDatoms, natoms, atomic_coords_device, points_device, dist_scratch_device
    );

  }
  const bool partition_weights_1d_kernel = true;

  if( partition_weights_1d_kernel ) {

    dim3 threads( warp_size, weight_thread_block / warp_size );
    dim3 blocks(  1, get_device_sm_count(0) * weight_thread_block_per_sm); 
    hipLaunchKernelGGL(modify_weights_ssf_kernel_2d, dim3(blocks), dim3(threads), 0, stream , 
      npts, LDatoms, natoms, rab_device, atomic_coords_device, dist_scratch_device, 
      iparent_device, dist_nearest_device, weights_device
    );

  } else {

#if 0
    dim3 threads( 32, 32 );
    dim3 blocks ( npts, 1 );

    if( weight_alg == XCWeightAlg::SSF ) 
      hipLaunchKernelGGL(modify_weights_ssf_kernel, dim3(blocks), dim3(threads), 0, stream , 
        npts, natoms, rab_device, atomic_coords_device, dist_scratch_device, 
        iparent_device, dist_nearest_device, weights_device
      );
    else
      hipLaunchKernelGGL(modify_weights_becke_kernel, dim3(blocks), dim3(threads), 0, stream , 
        npts, natoms, rab_device, atomic_coords_device, dist_scratch_device, 
        iparent_device, weights_device
      );
#endif

  }


}

template
void partition_weights_hip_SoA( XCWeightAlg    weight_alg,
                                 size_t         npts,
                                 size_t         LDatoms,
                                 size_t         natoms,
                                 const double*  points_device,
                                 const int32_t* iparent_device,
                                 const double*  dist_nearest_device,
                                 const double*  rab_device,
                                 const double*  atomic_coords_device,
                                       double*  weights_device,
                                       double*  dist_scratch_device,
                                 hipStream_t   stream );

}
}
}
