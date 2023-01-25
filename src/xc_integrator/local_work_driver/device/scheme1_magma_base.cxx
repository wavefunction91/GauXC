/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#include "scheme1_magma_base.hpp"
#include "device/common/pack_submat.hpp"
#include "device/common/inc_potential.hpp"
#include "device/common/device_blas.hpp"

namespace GauXC {

void AoSScheme1MAGMABase::eval_xmat( XCDeviceData* _data, bool do_grad ){

  if( do_grad ) GAUXC_GENERIC_EXCEPTION("MAGMA + X Gradient NYI");

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


  auto master_queue = data->device_backend_->master_magma_queue();
  auto magma_stack = data->magma_stack;
  magmablas_dgemm_vbatched( MagmaNoTrans, MagmaNoTrans,
    magma_stack.m_array_device, magma_stack.n_array_device, 
    magma_stack.k_array_device, 
    1., magma_stack.bf_array_device,   magma_stack.ld_bf_array_device,
        magma_stack.dmat_array_device, magma_stack.ld_dmat_array_device,
    0., magma_stack.zmat_array_device, magma_stack.ld_zmat_array_device,
    ntasks, *master_queue );

}

void AoSScheme1MAGMABase::inc_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  auto master_queue = data->device_backend_->master_magma_queue();
  auto magma_stack = data->magma_stack;
  magmablas_dsyr2k_vbatched( MagmaLower, MagmaTrans, 
    magma_stack.n_array_device, magma_stack.m_array_device,
    1., magma_stack.bf_array_device,   magma_stack.ld_bf_array_device, 
        magma_stack.zmat_array_device, magma_stack.ld_zmat_array_device,
    0., magma_stack.vmat_array_device, magma_stack.ld_vmat_array_device, 
    ntasks, *master_queue );

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
