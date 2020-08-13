#include <gauxc/util/div_ceil.hpp>

#include "collocation_petite_kernels.hpp"
#include "collocation_masked_kernels.hpp"

namespace GauXC      {
namespace integrator {
namespace cuda       {

template <typename T>
void eval_collocation_petite(
  size_t          nshells,
  size_t          nbf,
  size_t          npts,
  const Shell<T>* shells_device,
  const size_t*   offs_device,
  const T*        pts_device,
  T*              eval_device,
  cudaStream_t    stream
) {


  dim3 threads(32, 32, 1);
  dim3 blocks( util::div_ceil( npts,    threads.x ),
               util::div_ceil( nshells, threads.y ) );

  collocation_device_petite_kernel<T>
    <<<blocks, threads, 0, stream>>>
    ( nshells, nbf, npts, shells_device, offs_device,
      pts_device, eval_device );

}
 
template             
void eval_collocation_petite(
  size_t               nshells,
  size_t               nbf,
  size_t               npts,
  const Shell<double>* shells_device,
  const size_t*        offs_device,
  const double*        pts_device,
  double*              eval_device,
  cudaStream_t         stream
);









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

  dim3 threads(32, 32, 1);
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
void eval_collocation_petite_deriv1(
  size_t          nshells,
  size_t          nbf,
  size_t          npts,
  const Shell<T>* shells_device,
  const size_t*   offs_device,
  const T*        pts_device,
  T*              eval_device,
  T*              deval_device_x,
  T*              deval_device_y,
  T*              deval_device_z,
  cudaStream_t    stream
) {

  dim3 threads(32, 32, 1);
  dim3 blocks( util::div_ceil( npts,    threads.x ),
               util::div_ceil( nshells, threads.y ) );

  collocation_device_petite_kernel_deriv1<T>
    <<<blocks, threads, 0, stream>>>
    ( nshells, nbf, npts, shells_device, offs_device,
      pts_device, eval_device, deval_device_x, deval_device_y,
      deval_device_z );

}

template
void eval_collocation_petite_deriv1(
  size_t               nshells,
  size_t               nbf,
  size_t               npts,
  const Shell<double>* shells_device,
  const size_t*        offs_device,
  const double*        pts_device,
  double*              eval_device,
  double*              deval_device_x,
  double*              deval_device_y,
  double*              deval_device_z,
  cudaStream_t         stream
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

  dim3 threads(32, 32, 1);
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


} // namespace cuda
} // namespace integrator
} // namespace GauXC
