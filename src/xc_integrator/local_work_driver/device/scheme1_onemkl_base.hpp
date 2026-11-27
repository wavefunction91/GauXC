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
#pragma once
#include "scheme1_base.hpp"

namespace GauXC {

struct AoSScheme1OneMKLBase : public AoSScheme1Base {

  void eval_xmat( double fac, XCDeviceData*, bool do_grad, density_id den ) override final;
  void eval_exx_fmat( XCDeviceData* ) override final;
  void inc_exx_k( XCDeviceData* ) override final;

  struct Data;

  virtual ~AoSScheme1OneMKLBase() = default;
};

struct AoSScheme1OneMKLBase::Data : public AoSScheme1Base::Data {

  using base_type = AoSScheme1Base::Data;
  using base_type::host_task_type;
  using base_type::device_buffer_t;

  // oneMKL only exposes a variable-size batched GEMM through the "group" form
  // of gemm_batch: it takes one m/n/k/transpose/alpha/beta PER GROUP plus a
  // groupsize array, not one m/n/k per matrix the way MAGMA's vbatched
  // routines do. To express "every matrix in the batch has its own size"
  // (exactly what MAGMA's per-matrix arrays give for free), each task is
  // modeled as its own group of size 1: group_count == ntasks and
  // groupsize[i] == 1 for all i. That promotes transa/transb/alpha/beta from
  // MAGMA's single scalars-for-the-whole-batch to per-task arrays here.
  //
  // gemm_batch's group form is a host-side dispatch call: it walks
  // transa/transb/m/n/k/alpha/beta/lda/ldb/ldc/groupsize -- and the arrays of
  // matrix pointers themselves (a/b/c) -- on the host to build the batched
  // command lists before anything is enqueued to the device. The reference
  // oneMKL SYCL example (oneMKL's examples_sycl.tgz, sycl/blas/source/gemm_batch.cpp)
  // allocates every one of those parameter arrays, INCLUDING the
  // pointer-to-pointer matrix arrays a/b/c, with sycl::malloc_shared for
  // exactly this reason. Device-only (malloc_device) memory, which is what
  // the rest of this dynamic stack uses via buffer_adaptor, is not
  // host-readable and would make that host-side walk undefined behavior.
  // Every array passed directly to gemm_batch is therefore USM-shared here
  // (suffix "_shared"); only the matrix DATA the pointers point to (task.bf,
  // task.zmat, static_stack.dmat_s_device, ...) stays device-resident exactly
  // as it already is elsewhere in GauXC, since only the batched kernel itself
  // (running on-device) ever dereferences that data.
  struct onemkl_data {
    // Arrays of matrix pointers passed as gemm_batch's a/b/c. Shared USM: the
    // pointer VALUES must be host-readable for gemm_batch's dispatch, even
    // though each pointer refers to device-resident matrix data.
    double** xdmat_array_shared = nullptr;
    double** fdmat_array_shared = nullptr;
    double** vmat_array_shared  = nullptr;
    double** kmat_array_shared  = nullptr;
    double** zmat_array_shared  = nullptr;
    double** fmat_array_shared  = nullptr;
    double** gmat_array_shared  = nullptr;
    double** bf_array_shared    = nullptr;

    // Per-group (== per-task) sizes and leading dimensions.
    int32_t* xmat_m_array_shared   = nullptr;
    int32_t* xmat_n_array_shared   = nullptr;
    int32_t* xmat_k_array_shared   = nullptr;
    int32_t* ld_xdmat_array_shared = nullptr;

    int32_t* fmat_m_array_shared   = nullptr;
    int32_t* fmat_n_array_shared   = nullptr;
    int32_t* fmat_k_array_shared   = nullptr;
    int32_t* ld_fdmat_array_shared = nullptr;

    int32_t* ld_vmat_array_shared  = nullptr;
    int32_t* ld_zmat_array_shared  = nullptr;
    int32_t* ld_fmat_array_shared  = nullptr;
    int32_t* ld_bf_array_shared    = nullptr;

    // groupsize[i] == 1 for all i (each task is its own group); allocated and
    // filled once since it never changes across calls or task batches.
    int32_t* groupsize_array_shared = nullptr;

    // Constant per-group transpose flags. eval_xmat/eval_exx_fmat/inc_exx_k
    // each use one fixed (transa,transb) pair across the whole batch, but the
    // group API still requires a per-group array, so one array per call site
    // is filled once and reused on every call.
    // Stored as plain `char` (the underlying type of oneapi::mkl::transpose,
    // an enum class : char) rather than the oneMKL enum itself, so this
    // header -- like scheme1_magma_base.hpp -- stays free of any
    // backend-specific #include; the .cxx reinterpret_casts these to
    // oneapi::mkl::transpose* at the gemm_batch call site.
    char* xmat_transa_array_shared = nullptr;
    char* xmat_transb_array_shared = nullptr;
    char* fmat_transa_array_shared = nullptr;
    char* fmat_transb_array_shared = nullptr;
    char* kmat_transa_array_shared = nullptr;
    char* kmat_transb_array_shared = nullptr;

    // alpha varies per call (eval_xmat's fac is a call-site argument), so it
    // is refilled in eval_xmat itself on every call; beta is always 0 for
    // these calls but still needs a per-group array to satisfy the API, so
    // it is filled once in pack_and_send, like groupsize and the transposes.
    double* xmat_alpha_array_shared = nullptr;
    double* xmat_beta_array_shared  = nullptr;
    double* fmat_alpha_array_shared = nullptr;
    double* fmat_beta_array_shared  = nullptr;

    inline void reset(){ std::memset(this,0,sizeof(onemkl_data)); }
  };

  onemkl_data onemkl_stack;

  template <typename... Args>
  Data( Args&&... args ) : base_type( std::forward<Args>(args)... ) { }

  virtual ~Data();

  // get_mem_req/get_static_mem_requirement are deliberately NOT overridden:
  // unlike MAGMA's arrays, every onemkl_stack field is USM-shared memory
  // allocated directly with sycl::malloc_shared in allocate_dynamic_stack
  // (see .cxx), not carved out of the buffer_adaptor dynamic stack that
  // get_mem_req/get_static_mem_requirement account for. It therefore does not
  // affect, and must not be added to, the task-batch-size accounting those
  // two methods perform against the fixed device buffer.
  void reset_allocations() override final;
  device_buffer_t allocate_dynamic_stack( integrator_term_tracker terms,
    host_task_iterator begin, host_task_iterator end, device_buffer_t buf )
    override final;
  void pack_and_send( integrator_term_tracker terms,
    host_task_iterator begin, host_task_iterator end,
    const BasisSetMap& basis_map ) override final;


  void pack_and_send_xmat( host_task_iterator, host_task_iterator );
  void pack_and_send_fmat( host_task_iterator, host_task_iterator );

private:
  // All of onemkl_data's fields are USM-shared allocations made outside the
  // device-only buffer_adaptor dynamic stack (see allocate_dynamic_stack),
  // so they need their own explicit sycl::free at teardown; done from the
  // .cxx (the only place a SYCL queue is available to free them).
  void free_shared_arrays();
};

}
