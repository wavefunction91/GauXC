#include "cuda/cuda_inc_potential.hpp"
#include "cuda/cuda_device_properties.hpp"
#include <gauxc/util/div_ceil.hpp>

#include "cuda/cuda_device_properties.hpp"

namespace GauXC      {
namespace integrator {
namespace cuda       {

using namespace GauXC::cuda;


#define WARP_X 16
#define WARP_Y 1
#define UNROLL_FACTOR 4
#define EFF_UNROLL 4
#define CUT_X 8
#define CUT_Y 8


template <typename T>
__global__ __launch_bounds__(1024, 1)
void inc_by_submat_combined_kernel( size_t           ntasks,
                                    XCTaskDevice<T>* device_tasks,
                                    T*               A,
                                    size_t           LDA, 
				    const int block_y,
				    const int block_x ) {

  const int batch_id = blockIdx.z;
  auto& task = device_tasks[ batch_id ];

  const auto* submat_cut_device = task.submat_cut;
  const auto* submat_block_device = task.submat_block;
  const auto  LDAS              = task.nbe;
        auto* ASmall_device     = task.nbe_scr;

  //if( LDAS == LDAB ) return;
  const int tid_xx = threadIdx.x % WARP_X;
  const int tid_xy = threadIdx.x / WARP_X;

  const int tid_yx = threadIdx.y % CUT_X;
  const int tid_yy = threadIdx.y / CUT_X;

  const int start_cut_y = submat_block_device[block_y];
  const int end_cut_y   = submat_block_device[block_y+1];
  const int start_cut_x = submat_block_device[block_x];
  const int end_cut_x   = submat_block_device[block_x+1];

  for( int i_cut = tid_yy + start_cut_y; i_cut < end_cut_y; i_cut += CUT_Y ) {
    const int3 i_data = *((int3*)(submat_cut_device + 3*i_cut));
    const int i_cut_first  = i_data.x;
    const int delta_i      = i_data.y;
    const int i_cut_small  = i_data.z;

  for( int j_cut = tid_yx + start_cut_x; j_cut < end_cut_x; j_cut += CUT_X ) {
    const int3 j_data = *((int3*)(submat_cut_device + 3*j_cut));
    const int j_cut_first  = j_data.x; 
    const int delta_j      = j_data.y;
    const int j_cut_small  = j_data.z;

    auto* ASmall_begin = ASmall_device + i_cut_small + j_cut_small*LDAS;
    auto* ABig_begin   = A   + i_cut_first + j_cut_first*LDA;

    int J;
    for( J = tid_xy; J < (delta_j / EFF_UNROLL) * EFF_UNROLL; J += EFF_UNROLL ) {
      for( int I = tid_xx; I < delta_i; I += WARP_X ) {

        double val[UNROLL_FACTOR];
        double* address[UNROLL_FACTOR];
#pragma unroll
        for (int k = 0; k < UNROLL_FACTOR; k++) {
          val[k] = ASmall_begin[I + (J+k*WARP_Y)*LDAS];
          address[k] = ABig_begin + I + (J+k*WARP_Y)*LDA;
        }
#pragma unroll
        for (int k = 0; k < UNROLL_FACTOR; k++) {
          atomicAdd(address[k], val[k] );
        }
      }
    }

    for ( ; J < delta_j; J += WARP_Y) {
      for( int I = tid_xx; I < delta_i; I += WARP_X ) {
        atomicAdd(ABig_begin + I + J*LDA, ASmall_begin[I + J*LDAS] );
      }
    }

  }
  }
}


template <typename T>
void task_inc_potential( size_t           ntasks,
                         XCTaskDevice<T>* device_tasks,
                         T*               V_device,
                         size_t           LDV,
                         cudaStream_t     stream ) {
  dim3 threads(warp_size / 2, max_warps_per_thread_block * 2, 1), blocks(1,1,ntasks);

  const int submat_block_size = get_submat_cut_block(LDV, 0);
  for (int i = 0; i < util::div_ceil(LDV, submat_block_size); i++) {
    for (int j = 0; j < util::div_ceil(LDV, submat_block_size); j++) {
      inc_by_submat_combined_kernel<<< blocks, threads, 0, stream >>>(
        ntasks, device_tasks, V_device, LDV, i, j
      );
    }
  }
}

template 
void task_inc_potential( size_t                ntasks,
                         XCTaskDevice<double>* device_tasks,
                         double*               V_device,
                         size_t                LDV,
                         cudaStream_t          stream );

template <typename T>
__global__ void symmetrize_matrix_device( size_t nbf, size_t LDA, T* A ) {
  const size_t block_size = warp_size;

  __shared__ T buffer[block_size][block_size+1];  // Pad shared memory to resolve shared memory

  const size_t num_blocks = ((nbf + block_size - 1) / block_size);

  for (int i = blockIdx.x; i < num_blocks; i += gridDim.x) {
    // TODO This could be load balanced if need be
    const int i_coord = i * block_size;
    for (int j = i; j < num_blocks; j++) {
      const int j_coord = j * block_size;

      // Read in block to buffer
      // TODO These could be vector reads/writes if this becomes significant
      if (i_coord + threadIdx.y < nbf && j_coord + threadIdx.x < nbf) {
        buffer[threadIdx.y][threadIdx.x] = A[(i_coord + threadIdx.y) * LDA + j_coord + threadIdx.x];
      }
      __syncthreads();

      // Write buffer
      if (j_coord + threadIdx.y < nbf && i_coord + threadIdx.x < nbf) {
        if ((j_coord != i_coord || threadIdx.x < threadIdx.y)) { // handles the diagonal block
          A[(j_coord + threadIdx.y) * LDA + i_coord + threadIdx.x] = buffer[threadIdx.x][threadIdx.y];
        }
      }
      __syncthreads();
    }
  }
}

template <typename T>
void symmetrize_matrix( size_t nbf, size_t LDV, T* V_device, cudaStream_t stream) {
  const size_t num_blocks = ((LDV + warp_size - 1) / warp_size);
  // Warp size must equal max_warps_per_thread_block must equal 32
  dim3 threads(warp_size, max_warps_per_thread_block), blocks(num_blocks);
  symmetrize_matrix_device<<<blocks, threads, 0, stream>>>(nbf, LDV, V_device);
}

template
void symmetrize_matrix( size_t nbf, size_t LDV, double* V_device, cudaStream_t stream );


}
}
}

