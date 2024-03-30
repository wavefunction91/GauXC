/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include <stdexcept>

#include "scheme1_cutlass_base.hpp"
#include "device/common/pack_submat.hpp"
#include "device/common/inc_potential.hpp"
#include "device/common/device_blas.hpp"

#include "device/cuda/kernels/cutlass_wrapper.hpp"

namespace GauXC {

void AoSScheme1CUTLASSBase::eval_xmat(double fac, XCDeviceData* _data, bool do_grad ){

  if( do_grad ) GAUXC_GENERIC_EXCEPTION("CUTLASS + X Gradient NYI");

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  // Pack density matrix 
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  sym_pack_submat( ntasks, aos_stack.device_tasks, static_stack.dmat_device, 
    nbf, submat_block_size, data->device_backend_->queue() );

  auto cutlass_stack = data->cutlass_stack;
  cutlass_gemm(
    cutlass_stack.problem_sizes_device,
    data->problem_sizes_host.data(),
    ntasks,
    cutlass_stack.bf_array_device, cutlass_stack.dmat_array_device,
    cutlass_stack.zmat_array_device, cutlass_stack.zmat_array_device,
    cutlass_stack.ld64_bf_array_device, cutlass_stack.ld64_dmat_array_device,
    cutlass_stack.ld64_zmat_array_device, cutlass_stack.ld64_zmat_array_device,
    fac, 0.0,
    data->device_backend_->queue()
  );
}

void AoSScheme1CUTLASSBase::inc_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  auto cutlass_stack = data->cutlass_stack;
  cutlass_syr2k(
    cutlass_stack.syr2k_sizes_device,
    data->syr2k_sizes_host.data(),
    ntasks,
    cutlass_stack.bf_array_device, cutlass_stack.zmat_array_device,
    cutlass_stack.vmat_array_device, cutlass_stack.vmat_array_device,
    cutlass_stack.ld64_bf_array_device, cutlass_stack.ld64_zmat_array_device,
    cutlass_stack.ld64_vmat_array_device, cutlass_stack.ld64_vmat_array_device,
    1.0, 0.0,
    data->device_backend_->queue()
  );

  // Increment global VXC
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  sym_task_inc_potential( ntasks, aos_stack.device_tasks, 
    static_stack.vxc_device, nbf, submat_block_size, 
    data->device_backend_->queue() );
}

}
