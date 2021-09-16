#pragma once
#include <gauxc/shell.hpp>
#include "device/xc_device_task.hpp"
#include "device/type_erased_queue.hpp"

namespace GauXC {


template <typename T>
void eval_collocation_masked(
  size_t            nshells,
  size_t            nbf,
  size_t            npts,
  const Shell<T>*   shells_device,
  const size_t*     mask_device,
  const size_t*     offs_device,
  const T*          pts_device,
  T*                eval_device,
  type_erased_queue queue );


template <typename T>
void eval_collocation_masked_deriv1(
  size_t            nshells,
  size_t            nbf,
  size_t            npts,
  const Shell<T>*   shells_device,
  const size_t*     mask_device,
  const size_t*     offs_device,
  const T*          pts_device,
  T*                eval_device,
  T*                deval_device_x,
  T*                deval_device_y,
  T*                deval_device_z,
  type_erased_queue queue );


template <typename T>
void eval_collocation_masked_combined(
  size_t            ntasks,
  size_t            npts_max,
  size_t            nshells_max,
  Shell<T>*         shells_device,
  XCDeviceTask*     device_tasks,
  type_erased_queue queue );




template <typename T>
void eval_collocation_masked_combined_deriv1(
  size_t        ntasks,
  size_t        npts_max,
  size_t        nshells_max,
  Shell<T>*     shells_device,
  XCDeviceTask* device_tasks,
  type_erased_queue queue );








void eval_collocation_shell_to_task(
  size_t            nshells,
  Shell<double>*    shells_device,
  const int32_t**   shell_to_task_idx,
  const int32_t**   shell_to_task_off,
  const int32_t*    shell_to_task_ntask,
  const int32_t*    shell_to_task_ls,
  const int32_t*    shell_to_task_pure,
  XCDeviceTask*     device_tasks,
  type_erased_queue queue );

void eval_collocation_shell_to_task_gradient(
  size_t            nshells,
  Shell<double>*    shells_device,
  const int32_t**   shell_to_task_idx,
  const int32_t**   shell_to_task_off,
  const int32_t*    shell_to_task_ntask,
  const int32_t*    shell_to_task_ls,
  const int32_t*    shell_to_task_pure,
  XCDeviceTask*     device_tasks,
  type_erased_queue queue );

} // namespace GauXC
