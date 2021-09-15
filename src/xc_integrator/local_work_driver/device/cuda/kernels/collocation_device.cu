#include <gauxc/util/div_ceil.hpp>
#include <gauxc/util/cuda_util.hpp>
#include "exceptions/cuda_exception.hpp"
#include <gauxc/xc_task.hpp>

#include "device/cuda/kernels/collocation_masked_kernels.hpp"
#include "device/cuda/kernels/collocation_masked_combined_kernels.hpp"
#include "device/cuda/kernels/collocation_shell_to_task_kernels.hpp"

#include "device_specific/cuda_device_constants.hpp"

namespace GauXC {

 
template <typename T>
void eval_collocation_masked(
  size_t          nshells,
  size_t          nbf,
  size_t          npts,
  const Shell<T>* shells_device,
  const size_t*   mask_device,
  const size_t*   offs_device,
  const T*        pts_device,
  T*              eval_device,
  cudaStream_t    stream
) {

  auto nmax_threads = util::cuda_kernel_max_threads_per_block( 
    collocation_device_masked_kernel<T>
  );
  auto max_warps_per_thread_block = nmax_threads / cuda::warp_size;

  dim3 threads(cuda::warp_size, max_warps_per_thread_block, 1);
  dim3 blocks( util::div_ceil( npts,    threads.x ),
               util::div_ceil( nshells, threads.y ) );

  collocation_device_masked_kernel<T>
    <<<blocks, threads, 0, stream>>>
    ( nshells, nbf, npts, shells_device, mask_device,
      offs_device, pts_device, eval_device );

}
 
template             
void eval_collocation_masked(
  size_t               nshells,
  size_t               nbf,
  size_t               npts,
  const Shell<double>* shells_device,
  const size_t*        mask_device,
  const size_t*        offs_device,
  const double*        pts_device,
  double*              eval_device,
  cudaStream_t         stream
);




template <typename T>
void eval_collocation_masked_combined(
  size_t        ntasks,
  size_t        npts_max,
  size_t        nshells_max,
  Shell<T>*     shells_device,
  XCDeviceTask* device_tasks,
  cudaStream_t  stream
) {

  auto nmax_threads = util::cuda_kernel_max_threads_per_block( 
    collocation_device_masked_combined_kernel<T>
  );

  auto max_warps_per_thread_block = nmax_threads / cuda::warp_size;
  dim3 threads(cuda::warp_size, max_warps_per_thread_block, 1);
  dim3 blocks( util::div_ceil( npts_max,    threads.x ),
               util::div_ceil( nshells_max, threads.y ),
               ntasks );

  collocation_device_masked_combined_kernel<T>
    <<<blocks, threads, 0, stream>>>
    ( ntasks, shells_device, device_tasks );
     
}

template
void eval_collocation_masked_combined(
  size_t                ntasks,
  size_t                npts_max,
  size_t                nshells_max,
  Shell<double>*        shells_device,
  XCDeviceTask* device_tasks,
  cudaStream_t          stream
);











template <typename T>
void eval_collocation_masked_deriv1(
  size_t          nshells,
  size_t          nbf,
  size_t          npts,
  const Shell<T>* shells_device,
  const size_t*   mask_device,
  const size_t*   offs_device,
  const T*        pts_device,
  T*              eval_device,
  T*              deval_device_x,
  T*              deval_device_y,
  T*              deval_device_z,
  cudaStream_t    stream
) {

  auto nmax_threads = util::cuda_kernel_max_threads_per_block( 
    collocation_device_masked_combined_kernel<T>
  );

  auto max_warps_per_thread_block = nmax_threads / cuda::warp_size;
  dim3 threads(cuda::warp_size, max_warps_per_thread_block, 1);
  dim3 blocks( util::div_ceil( npts,    threads.x ),
               util::div_ceil( nshells, threads.y ) );

  collocation_device_masked_kernel_deriv1<T>
    <<<blocks, threads, 0, stream>>>
    ( nshells, nbf, npts, shells_device, mask_device, offs_device,
      pts_device, eval_device, deval_device_x, deval_device_y,
      deval_device_z );

}

template
void eval_collocation_masked_deriv1(
  size_t               nshells,
  size_t               nbf,
  size_t               npts,
  const Shell<double>* shells_device,
  const size_t*        mask_device,
  const size_t*        offs_device,
  const double*        pts_device,
  double*              eval_device,
  double*              deval_device_x,
  double*              deval_device_y,
  double*              deval_device_z,
  cudaStream_t         stream
);
















template <typename T>
void eval_collocation_masked_combined_deriv1(
  size_t        ntasks,
  size_t        npts_max,
  size_t        nshells_max,
  Shell<T>*     shells_device,
  XCDeviceTask* device_tasks,
  cudaStream_t  stream
) {

  auto nmax_threads = util::cuda_kernel_max_threads_per_block( 
    collocation_device_masked_combined_kernel_deriv1<T>
  );

  dim3 threads(cuda::warp_size, nmax_threads/cuda::warp_size, 1);
  dim3 blocks( util::div_ceil( npts_max,    threads.x ),
               util::div_ceil( nshells_max, threads.y ),
               ntasks );

  collocation_device_masked_combined_kernel_deriv1<T>
    <<<blocks, threads, 0, stream>>>
    ( ntasks, shells_device, device_tasks );
     
}

template
void eval_collocation_masked_combined_deriv1(
  size_t                ntasks,
  size_t                npts_max,
  size_t                nshells_max,
  Shell<double>*        shells_device,
  XCDeviceTask* device_tasks,
  cudaStream_t          stream
);




uint32_t max_threads_shell_to_task_collocation( int32_t l, bool pure ) {
  if( pure ) {
    switch(l) {
      case 0: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_cartesian_0 );
      case 1: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_spherical_1 );
      case 2: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_spherical_2 );
      case 3: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_spherical_3 );
    }
  } else {
    switch(l) {
      case 0: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_cartesian_0 );
      case 1: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_cartesian_1 );
      case 2: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_cartesian_2 );
      case 3: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_cartesian_3 );
    }
  }
  return 0;
}

template <typename... Args>
void dispatch_shell_to_task_collocation( cudaStream_t stream, int32_t l, bool pure, Args&&... args ) {
  dim3 block(10);
  dim3 threads = max_threads_shell_to_task_collocation(l,pure);
  if( pure ) {
    switch(l) {
      case 0:
        collocation_device_shell_to_task_kernel_cartesian_0<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
      case 1:
        collocation_device_shell_to_task_kernel_spherical_1<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
      case 2:
        collocation_device_shell_to_task_kernel_spherical_2<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
      case 3:
        collocation_device_shell_to_task_kernel_spherical_3<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
    }
  } else {
    switch(l) {
      case 0:
        collocation_device_shell_to_task_kernel_cartesian_0<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
      case 1:
        collocation_device_shell_to_task_kernel_cartesian_1<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
      case 2:
        collocation_device_shell_to_task_kernel_cartesian_2<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
      case 3:
        collocation_device_shell_to_task_kernel_cartesian_3<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
    }
  }
}


void eval_collocation_shell_to_task(
  size_t          nshells,
  Shell<double>*  shells_device,
  const int32_t** shell_to_task_idx,
  const int32_t** shell_to_task_off,
  const int32_t*  shell_to_task_ntask,
  const int32_t*  shell_to_task_ls,
  const int32_t*  shell_to_task_pure,
  XCDeviceTask*   device_tasks,
  cudaStream_t    stream
) {

  // Loop over shells
  for( auto i = 0; i < nshells; ++i ) {

    const auto* shell    = shells_device + i;
    const auto* task_idx = shell_to_task_idx[i];
    const auto* task_off = shell_to_task_off[i];
    const auto  ntask    = shell_to_task_ntask[i];
    const auto  l        = shell_to_task_ls[i];
    const auto  pure     = shell_to_task_pure[i];

    if( ntask ) {
      dispatch_shell_to_task_collocation( stream, l, pure, 
        ntask, shell, task_idx, task_off, device_tasks
      );
    }


  }

}


uint32_t max_threads_shell_to_task_collocation_gradient( int32_t l, bool pure ) {
  if( pure ) {
    switch(l) {
      case 0: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_cartesian_gradient_0 );
      case 1: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_spherical_gradient_1 );
      case 2: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_spherical_gradient_2 );
      case 3: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_spherical_gradient_3 );
    }
  } else {
    switch(l) {
      case 0: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_cartesian_gradient_0 );
      case 1: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_cartesian_gradient_1 );
      case 2: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_cartesian_gradient_2 );
      case 3: return util::cuda_kernel_max_threads_per_block( collocation_device_shell_to_task_kernel_cartesian_gradient_3 );
    }
  }
  return 0;
}

template <typename... Args>
void dispatch_shell_to_task_collocation_gradient( cudaStream_t stream, int32_t l, bool pure, Args&&... args ) {
  dim3 block(10);
  dim3 threads = max_threads_shell_to_task_collocation(l,pure);
  std::cout << threads.x << ", " << block.x << std::endl;
  if( pure ) {
    switch(l) {
      case 0:
        collocation_device_shell_to_task_kernel_cartesian_gradient_0<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
      case 1:
        collocation_device_shell_to_task_kernel_spherical_gradient_1<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
      case 2:
        collocation_device_shell_to_task_kernel_spherical_gradient_2<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
      case 3:
        collocation_device_shell_to_task_kernel_spherical_gradient_3<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
    }
  } else {
    switch(l) {
      case 0:
        collocation_device_shell_to_task_kernel_cartesian_gradient_0<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
      case 1:
        collocation_device_shell_to_task_kernel_cartesian_gradient_1<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
      case 2:
        collocation_device_shell_to_task_kernel_cartesian_gradient_2<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
      case 3:
        collocation_device_shell_to_task_kernel_cartesian_gradient_3<<<block,threads,0,stream>>>(
          std::forward<Args>(args)... );
        break;
    }
  }

  auto err = cudaDeviceSynchronize();
  GAUXC_CUDA_ERROR( "COLLOCATION KERNEL", err );
}


void eval_collocation_shell_to_task_gradient(
  size_t          nshells,
  Shell<double>*  shells_device,
  const int32_t** shell_to_task_idx,
  const int32_t** shell_to_task_off,
  const int32_t*  shell_to_task_ntask,
  const int32_t*  shell_to_task_ls,
  const int32_t*  shell_to_task_pure,
  XCDeviceTask*   device_tasks,
  cudaStream_t    stream
) {

  // Loop over shells
  for( auto i = 0; i < nshells; ++i ) {

    const auto* shell    = shells_device + i;
    const auto* task_idx = shell_to_task_idx[i];
    const auto* task_off = shell_to_task_off[i];
    const auto  ntask    = shell_to_task_ntask[i];
    const auto  l        = shell_to_task_ls[i];
    const auto  pure     = shell_to_task_pure[i];

    if( ntask ) {
    //if(l == 0)
      dispatch_shell_to_task_collocation_gradient( stream, l, pure, 
        ntask, shell, task_idx, task_off, device_tasks
      );
    //else
    //  dispatch_shell_to_task_collocation( stream, l, pure, 
    //    ntask, shell, task_idx, task_off, device_tasks
    //  );
    }


  }

}





} // namespace GauXC
