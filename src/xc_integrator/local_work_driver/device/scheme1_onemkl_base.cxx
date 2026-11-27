/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "scheme1_onemkl_base.hpp"
#include "device/common/pack_submat.hpp"
#include "device/common/inc_potential.hpp"
#include "device/common/device_blas.hpp"
#include "device/sycl/sycl_backend.hpp"

// Only the BLAS domain is needed. The umbrella <oneapi/mkl.hpp> also pulls in
// the C LAPACK prototypes, whose parameter named "csl" collides with the csl
// macro from buffer_adaptor.hpp
#include <oneapi/mkl/blas.hpp>

#ifdef GAUXC_HAS_CUDA
#define GAUXC_ENABLE_EXX
#endif

namespace GauXC {

namespace {
  // onemkl_data stores transpose flags as plain char (see scheme1_onemkl_base.hpp
  // for why); reinterpret to the real oneMKL type at the call site.
  inline oneapi::mkl::transpose* as_transpose( char* p ) {
    return reinterpret_cast<oneapi::mkl::transpose*>(p);
  }
  inline const oneapi::mkl::transpose* as_transpose( const char* p ) {
    return reinterpret_cast<const oneapi::mkl::transpose*>(p);
  }
}

void AoSScheme1OneMKLBase::eval_xmat( double fac, XCDeviceData* _data, bool do_grad, density_id den ){

  if( do_grad ) GAUXC_GENERIC_EXCEPTION("oneMKL + X Gradient NYI");

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
  auto* backend = dynamic_cast<SYCLBackend*>(data->device_backend_);
  if( !backend ) GAUXC_BAD_BACKEND_CAST();
  ::sycl::queue& master_queue = backend->master_stream->queue;
  auto& onemkl_stack = data->onemkl_stack;
  double* dmat_ptr   = nullptr;
  switch( den ) {
    case DEN_S:
      dmat_ptr = static_stack.dmat_s_device;
      break;
    case DEN_Z:
      dmat_ptr = static_stack.dmat_z_device;
      break;
    case DEN_Y:
      dmat_ptr = static_stack.dmat_y_device;
      break;
    case DEN_X:
      dmat_ptr = static_stack.dmat_x_device;
      break;
    default:
      GAUXC_GENERIC_EXCEPTION( "eval_xmat called with invalid density specifier" );
  }
  sym_pack_submat( ntasks, aos_stack.device_tasks, dmat_ptr,
    nbf, submat_block_size, data->device_backend_->queue() );

  // Update dmat pointers on onemkl_stack. xdmat_array_shared is USM-shared,
  // so this can be written directly from the host instead of routed through
  // device_backend_->copy_async as MAGMA's device-only array requires.
  for( decltype(tasks.size()) i = 0; i < ntasks; i++ ) {
    auto& task = tasks[i];
    if( task.bfn_screening.ncut > 1 ) {
      onemkl_stack.xdmat_array_shared[i] = task.nbe_scr;
    } else {
      onemkl_stack.xdmat_array_shared[i] = dmat_ptr + task.bfn_screening.ibf_begin*(nbf+1);
    }
  }

  // fac is only known at this call site (e.g. 2.0 for RKS), unlike beta which
  // is always 0 here, so alpha is refilled on every call; groupsize, m/n/k,
  // leading dimensions, and the (constant) transpose flags were already
  // filled once in pack_and_send_xmat.
  for( decltype(tasks.size()) i = 0; i < ntasks; i++ )
    onemkl_stack.xmat_alpha_array_shared[i] = fac;

  // MAGMA's vbatched GEMM takes one m/n/k per matrix. oneMKL's group-form
  // gemm_batch instead takes one m/n/k/transpose per GROUP plus a groupsize
  // array; "every matrix has its own size" is expressed by making each task
  // its own group of size 1 (group_count == ntasks, groupsize[i] == 1).
  GAUXC_ONEMKL_ERROR( "ONEMKL DGEMM_BATCH (XMAT) FAILED",
    oneapi::mkl::blas::column_major::gemm_batch( master_queue,
      as_transpose(onemkl_stack.xmat_transa_array_shared),
      as_transpose(onemkl_stack.xmat_transb_array_shared),
      onemkl_stack.xmat_m_array_shared, onemkl_stack.xmat_n_array_shared,
      onemkl_stack.xmat_k_array_shared, onemkl_stack.xmat_alpha_array_shared,
      (const double**)onemkl_stack.bf_array_shared,    onemkl_stack.ld_bf_array_shared,
      (const double**)onemkl_stack.xdmat_array_shared, onemkl_stack.ld_xdmat_array_shared,
      onemkl_stack.xmat_beta_array_shared, onemkl_stack.zmat_array_shared, onemkl_stack.ld_zmat_array_shared,
      ntasks, onemkl_stack.groupsize_array_shared ).wait() );

}

void AoSScheme1OneMKLBase::eval_exx_fmat( XCDeviceData* _data ) {
#ifndef GAUXC_ENABLE_EXX
  GAUXC_GENERIC_EXCEPTION("EXX + non-CUDA NYI");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  // Pack the density matrix into (bfn, cou) shape
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  asym_pack_submat( ntasks, aos_stack.device_tasks, static_stack.dmat_s_device,
    nbf, submat_block_size, data->device_backend_->queue() );

  auto* backend = dynamic_cast<SYCLBackend*>(data->device_backend_);
  if( !backend ) GAUXC_BAD_BACKEND_CAST();
  ::sycl::queue& master_queue = backend->master_stream->queue;
  auto& onemkl_stack = data->onemkl_stack;
  GAUXC_ONEMKL_ERROR( "ONEMKL DGEMM_BATCH (FMAT) FAILED",
    oneapi::mkl::blas::column_major::gemm_batch( master_queue,
      as_transpose(onemkl_stack.fmat_transa_array_shared),
      as_transpose(onemkl_stack.fmat_transb_array_shared),
      onemkl_stack.fmat_m_array_shared, onemkl_stack.fmat_n_array_shared,
      onemkl_stack.fmat_k_array_shared, onemkl_stack.fmat_alpha_array_shared,
      (const double**)onemkl_stack.bf_array_shared,    onemkl_stack.ld_bf_array_shared,
      (const double**)onemkl_stack.fdmat_array_shared, onemkl_stack.ld_fdmat_array_shared,
      onemkl_stack.fmat_beta_array_shared, onemkl_stack.fmat_array_shared, onemkl_stack.ld_fmat_array_shared,
      ntasks, onemkl_stack.groupsize_array_shared ).wait() );
#endif
}

// inc_vxc is intentionally NOT overridden here: there is no batched syr2k in
// oneMKL (only syrk_batch, see usm_decls.hpp), and the correct decomposition
// of syr2k C := alpha*(A^T B + B^T A) + beta*C into two GEMMs is exactly what
// AoSScheme1Base::inc_vxc already does via inc_potential_impl's round-robin
// do_syr2k lambda, which itself calls the device_blas syr2k() entry point
// already specialized for oneMKL in kernels/onemkl_extensions.cxx (over the
// BLAS stream pool). Substituting anything else here would be a numerically
// different, and still unbatched (no group syr2k exists to batch it with),
// implementation for no benefit, so the base class implementation stands.

void AoSScheme1OneMKLBase::inc_exx_k( XCDeviceData* _data){
#ifndef GAUXC_ENABLE_EXX
  GAUXC_GENERIC_EXCEPTION("EXX + non-CUDA NYI");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  auto* backend = dynamic_cast<SYCLBackend*>(data->device_backend_);
  if( !backend ) GAUXC_BAD_BACKEND_CAST();
  ::sycl::queue& master_queue = backend->master_stream->queue;
  auto& onemkl_stack = data->onemkl_stack;
  GAUXC_ONEMKL_ERROR( "ONEMKL DGEMM_BATCH (EXX_K) FAILED",
    oneapi::mkl::blas::column_major::gemm_batch( master_queue,
      as_transpose(onemkl_stack.kmat_transa_array_shared),
      as_transpose(onemkl_stack.kmat_transb_array_shared),
      onemkl_stack.fmat_k_array_shared, onemkl_stack.fmat_n_array_shared,
      onemkl_stack.fmat_m_array_shared, onemkl_stack.fmat_alpha_array_shared,
      (const double**)onemkl_stack.bf_array_shared,   onemkl_stack.ld_bf_array_shared,
      (const double**)onemkl_stack.gmat_array_shared, onemkl_stack.ld_fmat_array_shared,
      onemkl_stack.fmat_beta_array_shared, onemkl_stack.kmat_array_shared, onemkl_stack.ld_fdmat_array_shared,
      ntasks, onemkl_stack.groupsize_array_shared ).wait() );

  // Increment EXX_K
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  asym_task_inc_potential( ntasks, aos_stack.device_tasks,
    static_stack.exx_k_device, nbf, submat_block_size,
    data->device_backend_->queue() );
#endif
}

}
