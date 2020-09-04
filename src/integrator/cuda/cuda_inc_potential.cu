#include "cuda_inc_potential.hpp"

namespace GauXC      {
namespace integrator {
namespace cuda       {


template <typename T>
__global__ void inc_by_submat_combined_kernel( size_t           ntasks,
                                               XCTaskDevice<T>* device_tasks,
                                               T*               A,
                                               size_t           LDA ) {

  const int batch_id = blockIdx.z;

  if( batch_id < ntasks ) {

  auto& task = device_tasks[ batch_id ];

  const auto  ncut              = task.ncut;
  const auto* submat_cut_device = task.submat_cut;
  const auto  LDAS              = task.nbe;
        auto* ASmall_device     = task.nbe_scr;

  //if( LDAS == LDAB ) return;


  const int tid_x = blockDim.x * blockIdx.x + threadIdx.x;
  const int tid_y = blockDim.y * blockIdx.y + threadIdx.y;

  int64_t i(0);
  for( size_t i_cut = 0; i_cut < ncut; ++i_cut ) {
    const int64_t i_cut_first  = submat_cut_device[ 2*i_cut ];
    const int64_t i_cut_second = submat_cut_device[ 2*i_cut + 1 ];
    const int64_t delta_i      = i_cut_second - i_cut_first;

    int64_t j(0);
  for( size_t j_cut = 0; j_cut < ncut; ++j_cut ) {
    const int64_t j_cut_first  = submat_cut_device[ 2*j_cut ];
    const int64_t j_cut_second = submat_cut_device[ 2*j_cut + 1 ];
    const int64_t delta_j      = j_cut_second - j_cut_first;

    auto* ASmall_begin = ASmall_device + i           + j          *LDAS;
    auto* ABig_begin   = A             + i_cut_first + j_cut_first*LDA ;
    
    for( size_t J = tid_y; J < delta_j; J += blockDim.y )      
    for( size_t I = tid_x; I < delta_i; I += blockDim.x )
      //ABig_begin[I + J*LDA] += ASmall_begin[I + J*LDAS];
      atomicAdd( ABig_begin + I + J*LDA, ASmall_begin[I+J*LDAS] );

    j += delta_j;
  }
    i += delta_i;
  }

  } // batch_id check
}


template <typename T>
void task_inc_potential( size_t           ntasks,
                         XCTaskDevice<T>* device_tasks,
                         T*               V_device,
                         size_t           LDV,
                         cudaStream_t     stream ) {

  dim3 threads(32,32,1), blocks(1,1,ntasks);
  inc_by_submat_combined_kernel<<< blocks, threads, 0, stream >>>(
    ntasks, device_tasks, V_device, LDV
  );


}

template 
void task_inc_potential( size_t                ntasks,
                         XCTaskDevice<double>* device_tasks,
                         double*               V_device,
                         size_t                LDV,
                         cudaStream_t          stream );

}
}
}

