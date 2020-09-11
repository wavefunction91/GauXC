#pragma once
#include <gauxc/shell.hpp>
#include <gauxc/xc_task.hpp>

namespace GauXC      {
namespace integrator {
namespace cuda       {

using namespace GauXC::cuda;

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
);

template <typename T>
void eval_collocation_petite_combined(
  size_t           ntasks,
  size_t           npts_max,
  size_t           nshells_max,
  XCTaskDevice<T>* device_tasks,
  cudaStream_t     stream
);

template <typename T>
void eval_collocation_masked_combined(
  size_t           ntasks,
  size_t           npts_max,
  size_t           nshells_max,
  Shell<T>*        shells_device,
  XCTaskDevice<T>* device_tasks,
  cudaStream_t     stream
);



template <typename T>
void eval_collocation_petite_combined_deriv1(
  size_t           ntasks,
  size_t           npts_max,
  size_t           nshells_max,
  XCTaskDevice<T>* device_tasks,
  cudaStream_t     stream
);

template <typename T>
void eval_collocation_masked_combined_deriv1(
  size_t           ntasks,
  size_t           npts_max,
  size_t           nshells_max,
  Shell<T>*        shells_device,
  XCTaskDevice<T>* device_tasks,
  cudaStream_t     stream
);

} // namespace cuda
} // namespace integrator
} // namespace GauXC
