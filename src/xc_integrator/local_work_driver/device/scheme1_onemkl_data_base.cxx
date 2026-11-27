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
#include "buffer_adaptor.hpp"
#include "device/sycl/sycl_backend.hpp"
#include "device_specific/sycl_util.hpp"
#include <type_traits>

// Only the BLAS domain is needed for oneapi::mkl::transpose -- see
// scheme1_onemkl_base.cxx for why the umbrella <oneapi/mkl.hpp> is avoided.
#include <oneapi/mkl/blas.hpp>

namespace GauXC {

AoSScheme1OneMKLBase::Data::~Data() {
  free_shared_arrays();
}

// onemkl_stack holds no device-only allocations of its own (unlike MAGMA's
// magma_stack, whose double**/int32_t* arrays live in the buffer_adaptor
// dynamic stack) -- every field here is USM-shared memory allocated and
// freed independently in allocate_dynamic_stack/free_shared_arrays. There is
// therefore nothing for reset_allocations to hand back to the dynamic stack;
// it only needs to release the shared allocations from the previous task
// batch before the next one is sized, exactly mirroring what
// base_type::reset_allocations() does for the (device-only) base stacks.
void AoSScheme1OneMKLBase::Data::reset_allocations() {
  base_type::reset_allocations();
  free_shared_arrays();
  onemkl_stack.reset();
}

void AoSScheme1OneMKLBase::Data::free_shared_arrays() {
  auto* backend = dynamic_cast<SYCLBackend*>(device_backend_);
  if( !backend ) return; // Not yet bound to a backend, or already torn down
  auto& queue = backend->master_stream->queue;
  auto& s = onemkl_stack;

  auto free_if = [&]( auto*& ptr ) {
    if( ptr ) { util::sycl_free( queue, ptr ); }
  };

  free_if(s.xdmat_array_shared); free_if(s.fdmat_array_shared);
  free_if(s.vmat_array_shared);  free_if(s.kmat_array_shared);
  free_if(s.zmat_array_shared);  free_if(s.fmat_array_shared);
  free_if(s.gmat_array_shared);  free_if(s.bf_array_shared);

  free_if(s.xmat_m_array_shared); free_if(s.xmat_n_array_shared);
  free_if(s.xmat_k_array_shared); free_if(s.ld_xdmat_array_shared);

  free_if(s.fmat_m_array_shared); free_if(s.fmat_n_array_shared);
  free_if(s.fmat_k_array_shared); free_if(s.ld_fdmat_array_shared);

  free_if(s.ld_vmat_array_shared); free_if(s.ld_zmat_array_shared);
  free_if(s.ld_fmat_array_shared); free_if(s.ld_bf_array_shared);

  free_if(s.groupsize_array_shared);

  free_if(s.xmat_transa_array_shared); free_if(s.xmat_transb_array_shared);
  free_if(s.fmat_transa_array_shared); free_if(s.fmat_transb_array_shared);
  free_if(s.kmat_transa_array_shared); free_if(s.kmat_transb_array_shared);

  free_if(s.xmat_alpha_array_shared); free_if(s.xmat_beta_array_shared);
  free_if(s.fmat_alpha_array_shared); free_if(s.fmat_beta_array_shared);
}

// allocate_dynamic_stack is otherwise the point at which per-batch device
// memory is carved out of the fixed device_buffer_t handed down from
// XCDeviceStackData::generate_buffers. onemkl_stack's arrays are excluded
// from that accounting (see the .hpp) because the group-form gemm_batch
// requires them to be host-readable USM-shared memory, which buffer_adaptor
// -- built to slice a single device-only allocation -- cannot provide.
// Instead they are allocated directly here via sycl::malloc_shared, sized to
// the same ntask this task batch was already sized to by the base class.
AoSScheme1OneMKLBase::Data::device_buffer_t
  AoSScheme1OneMKLBase::Data::allocate_dynamic_stack(
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end,
  device_buffer_t buf ){

  // Allocate base info on the stack
  buf = base_type::allocate_dynamic_stack( terms, task_begin, task_end,
    buf );

  required_term_storage reqt(terms);

  auto* backend = dynamic_cast<SYCLBackend*>(device_backend_);
  if( !backend ) GAUXC_BAD_BACKEND_CAST();
  auto& queue = backend->master_stream->queue;

  // A single operation may need more than one task batch to fit the fixed
  // device dynamic-stack budget, meaning allocate_dynamic_stack can be called
  // more than once between two reset_allocations() calls. Unlike MAGMA's
  // buffer_adaptor slices (harmlessly re-carved from the same fixed device
  // buffer every call), onemkl_stack's fields are independent
  // sycl::malloc_shared allocations, so the previous batch's arrays must be
  // freed here -- not only in reset_allocations() -- or every batch after the
  // first in a multi-batch operation leaks.
  free_shared_arrays();

  const auto ntask = std::distance( task_begin, task_end );
  auto& s = onemkl_stack;

  // bf_array_shared and groupsize_array_shared are shared between the xmat
  // and fmat call sites (same bf pointers, same "every task is its own
  // group of 1" groupsize semantics), so -- unlike MAGMA's buffer_adaptor
  // slices, where allocating twice just wastes stack space -- they must only
  // be sycl::malloc_shared'd once or the first allocation leaks.
  if(reqt.task_xmat or reqt.task_fmat) {
    s.bf_array_shared        = util::sycl_malloc_shared<double*>(ntask, queue);
    s.ld_bf_array_shared     = util::sycl_malloc_shared<int32_t>(ntask, queue);
    s.groupsize_array_shared = util::sycl_malloc_shared<int32_t>(ntask, queue);
  }

  if(reqt.task_xmat) {
    s.xdmat_array_shared = util::sycl_malloc_shared<double*>(ntask, queue);
    s.zmat_array_shared  = util::sycl_malloc_shared<double*>(ntask, queue);

    s.xmat_m_array_shared   = util::sycl_malloc_shared<int32_t>(ntask, queue);
    s.xmat_n_array_shared   = util::sycl_malloc_shared<int32_t>(ntask, queue);
    s.xmat_k_array_shared   = util::sycl_malloc_shared<int32_t>(ntask, queue);
    s.ld_xdmat_array_shared = util::sycl_malloc_shared<int32_t>(ntask, queue);
    s.ld_zmat_array_shared  = util::sycl_malloc_shared<int32_t>(ntask, queue);

    s.xmat_transa_array_shared = util::sycl_malloc_shared<char>(ntask, queue);
    s.xmat_transb_array_shared = util::sycl_malloc_shared<char>(ntask, queue);
    s.xmat_alpha_array_shared  = util::sycl_malloc_shared<double>(ntask, queue);
    s.xmat_beta_array_shared   = util::sycl_malloc_shared<double>(ntask, queue);
  }

  if(reqt.task_fmat) {
    s.fdmat_array_shared = util::sycl_malloc_shared<double*>(ntask, queue);
    s.fmat_array_shared  = util::sycl_malloc_shared<double*>(ntask, queue);
    s.gmat_array_shared  = util::sycl_malloc_shared<double*>(ntask, queue);
    s.kmat_array_shared  = util::sycl_malloc_shared<double*>(ntask, queue);

    s.fmat_m_array_shared   = util::sycl_malloc_shared<int32_t>(ntask, queue);
    s.fmat_n_array_shared   = util::sycl_malloc_shared<int32_t>(ntask, queue);
    s.fmat_k_array_shared   = util::sycl_malloc_shared<int32_t>(ntask, queue);
    s.ld_fdmat_array_shared = util::sycl_malloc_shared<int32_t>(ntask, queue);
    s.ld_fmat_array_shared  = util::sycl_malloc_shared<int32_t>(ntask, queue);

    s.fmat_transa_array_shared = util::sycl_malloc_shared<char>(ntask, queue);
    s.fmat_transb_array_shared = util::sycl_malloc_shared<char>(ntask, queue);
    s.kmat_transa_array_shared = util::sycl_malloc_shared<char>(ntask, queue);
    s.kmat_transb_array_shared = util::sycl_malloc_shared<char>(ntask, queue);
    s.fmat_alpha_array_shared  = util::sycl_malloc_shared<double>(ntask, queue);
    s.fmat_beta_array_shared   = util::sycl_malloc_shared<double>(ntask, queue);
  }

  // Nothing was taken from the device-only dynamic stack, so hand it back
  // untouched for further derived-class allocations.
  return buf;
}

void AoSScheme1OneMKLBase::Data::pack_and_send(
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end,
  const BasisSetMap& basis_map ) {


  base_type::pack_and_send( terms, task_begin, task_end, basis_map );
  required_term_storage reqt(terms);

  if(reqt.task_xmat) pack_and_send_xmat(task_begin, task_end);
  if(reqt.task_fmat) pack_and_send_fmat(task_begin, task_end);

}


void AoSScheme1OneMKLBase::Data::pack_and_send_xmat(
  host_task_iterator task_begin, host_task_iterator task_end
) {

  const auto ntask = std::distance( task_begin, task_end );
  auto& s = onemkl_stack;

  const auto nbf = global_dims.nbf;

  // Every field written here is USM-shared, so it can be filled directly
  // from the host without a device_backend_->copy_async round trip -- unlike
  // MAGMA's device-only magma_stack, which requires exactly that.
  //
  // xdmat_array_shared is intentionally NOT filled here: the density pointer
  // depends on which of DEN_S/DEN_Z/DEN_Y/DEN_X eval_xmat is called for, so
  // eval_xmat fills it (and xmat_alpha_array_shared, since alpha == fac is
  // also only known at that call site) itself, on every call.
  for( std::decay_t<decltype(ntask)> i = 0; i < ntask; ++i ) {
    auto& task = host_device_tasks[i];
    s.zmat_array_shared[i] = task.zmat;   s.ld_zmat_array_shared[i]  = task.npts;
    s.bf_array_shared[i]   = task.bf;     s.ld_bf_array_shared[i]    = task.npts;

    if( task.bfn_screening.ncut > 1 ) {
      s.ld_xdmat_array_shared[i] = task.bfn_screening.nbe;
    } else {
      s.ld_xdmat_array_shared[i] = nbf;
    }

    s.xmat_m_array_shared[i] = task.npts;
    s.xmat_n_array_shared[i] = task.bfn_screening.nbe;
    s.xmat_k_array_shared[i] = task.bfn_screening.nbe;

    s.groupsize_array_shared[i]   = 1;
    s.xmat_transa_array_shared[i] = static_cast<char>(oneapi::mkl::transpose::nontrans);
    s.xmat_transb_array_shared[i] = static_cast<char>(oneapi::mkl::transpose::nontrans);
    s.xmat_beta_array_shared[i]   = 0.;
  }

}



void AoSScheme1OneMKLBase::Data::pack_and_send_fmat(
  host_task_iterator task_begin, host_task_iterator task_end
) {

  const auto ntask = std::distance( task_begin, task_end );
  auto& s = onemkl_stack;

  // host_device_tasks should be populated by parent impl called at top
  for( std::decay_t<decltype(ntask)> i = 0; i < ntask; ++i ) {
    auto& task = host_device_tasks[i];
    s.fmat_array_shared[i] = task.fmat;    s.ld_fmat_array_shared[i] = task.npts;
    s.gmat_array_shared[i] = task.gmat;
    s.bf_array_shared[i]   = task.bf;      s.ld_bf_array_shared[i]   = task.npts;
    s.kmat_array_shared[i] = task.nbe_scr;

    s.fdmat_array_shared[i]    = task.nbe_scr;
    s.ld_fdmat_array_shared[i] = task.bfn_screening.nbe;

    s.fmat_m_array_shared[i] = task.npts;
    s.fmat_n_array_shared[i] = task.cou_screening.nbe;
    s.fmat_k_array_shared[i] = task.bfn_screening.nbe;

    s.groupsize_array_shared[i]   = 1;
    s.fmat_transa_array_shared[i] = static_cast<char>(oneapi::mkl::transpose::nontrans);
    s.fmat_transb_array_shared[i] = static_cast<char>(oneapi::mkl::transpose::nontrans);
    s.fmat_beta_array_shared[i]   = 0.;

    // inc_exx_k reuses these same groupsize/beta arrays with its own
    // (transposed) transa/transb pair.
    s.kmat_transa_array_shared[i] = static_cast<char>(oneapi::mkl::transpose::trans);
    s.kmat_transb_array_shared[i] = static_cast<char>(oneapi::mkl::transpose::nontrans);
  }

}
}
