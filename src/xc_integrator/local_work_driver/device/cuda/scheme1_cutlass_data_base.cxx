#include "scheme1_cutlass_base.hpp"
#include "buffer_adaptor.hpp"

#include "cutlass/cutlass.h"
#include "cutlass/gemm/gemm.h"

namespace GauXC {


void AoSScheme1CUTLASSBase::Data::reset_allocations() {
  base_type::reset_allocations();
  cutlass_stack.reset();
}

size_t AoSScheme1CUTLASSBase::Data::get_static_mem_requirement() {
  return base_type::get_static_mem_requirement() + 8 * sizeof(int32_t); 
    // Extra elements in CUTLASS dimension arrays
}

size_t AoSScheme1CUTLASSBase::Data::get_mem_req( integrator_term_tracker terms, 
  const host_task_type& task ) {

  
  size_t base_size = base_type::get_mem_req(terms, task);

  // All CUTLASS memory is related to exc/vxc +
  if( not terms.exc_vxc ) return base_size;

  // TODO Update this for cutlass
  return base_size + 
         4*sizeof(double*) + // batch device pointers
         8*sizeof(int32_t) +
         8*sizeof(int64_t);  // Dimensions + leading dimensions 
                             // (extra handled by get_static_mem_requirement)


}
AoSScheme1CUTLASSBase::Data::device_buffer_t 
  AoSScheme1CUTLASSBase::Data::allocate_dynamic_stack( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end, 
  device_buffer_t buf ){

  // Allocate base info on the stack
  buf = base_type::allocate_dynamic_stack( terms, task_begin, task_end,
    buf );

  // All CUTLASS memory is related to exc/vxc +
  if( not terms.exc_vxc ) return buf;

  // Allocate additional device memory 
  auto [ ptr, sz ] = buf;
  buffer_adaptor mem( ptr, sz );

  const auto ntask = std::distance( task_begin, task_end );
  cutlass_stack.dmat_array_device = mem.aligned_alloc<double*>( ntask, csl );
  cutlass_stack.vmat_array_device = mem.aligned_alloc<double*>( ntask, csl );
  cutlass_stack.zmat_array_device = mem.aligned_alloc<double*>( ntask, csl );
  cutlass_stack.bf_array_device   = mem.aligned_alloc<double*>( ntask, csl );

  cutlass_stack.ld64_dmat_array_device = mem.aligned_alloc<int64_t>( ntask + 1, csl );
  cutlass_stack.ld64_zmat_array_device = mem.aligned_alloc<int64_t>( ntask + 1, csl );
  cutlass_stack.ld64_vmat_array_device = mem.aligned_alloc<int64_t>( ntask + 1, csl );
  cutlass_stack.ld64_bf_array_device   = mem.aligned_alloc<int64_t>( ntask + 1, csl );
  
  cutlass_stack.problem_sizes_device = mem.aligned_alloc<cutlass::gemm::GemmCoord>(  ntask + 1, csl );
  cutlass_stack.syr2k_sizes_device   = mem.aligned_alloc<cutlass::gemm::GemmCoord>(  ntask + 1, csl );

  // Update dynmem data for derived impls
  return device_buffer_t{ mem.stack(), mem.nleft() };
}

void AoSScheme1CUTLASSBase::Data::pack_and_send( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end,
  const BasisSetMap& basis_map ) {

  base_type::pack_and_send( terms, task_begin, task_end, basis_map );
  if( not terms.exc_vxc ) return;

  const auto ntask = std::distance( task_begin, task_end );
  std::vector<double*> dmat_host( ntask ), zmat_host( ntask ), bf_host( ntask ),
                       vmat_host( ntask );
  std::vector<cutlass::gemm::GemmCoord> problem_sizes( ntask );
  std::vector<cutlass::gemm::GemmCoord> sym2k_sizes( ntask );
  std::vector<int64_t> ld64_dmat_host( ntask ), ld64_zmat_host( ntask ), 
                       ld64_vmat_host( ntask ), ld64_bf_host( ntask );

  double* static_dmat = static_stack.dmat_device;
  const auto nbf = global_dims.nbf;

  // host_device_tasks should be populated by parent impl called at top
  for( auto i = 0; i < ntask; ++i ) {
    auto& task = host_device_tasks[i];
    zmat_host[i] = task.zmat;    ld64_zmat_host[i] = task.npts;
    bf_host[i]   = task.bf;      ld64_bf_host[i]   = task.npts;
    vmat_host[i] = task.nbe_scr; ld64_vmat_host[i] = task.nbe;
    if( task.ncut > 1 ) {
      dmat_host[i]    = task.nbe_scr;
      ld64_dmat_host[i] = task.nbe;
    } else {
      dmat_host[i]    = static_dmat + task.ibf_begin*(nbf+1);
      ld64_dmat_host[i] = nbf;
    }

    cutlass::gemm::GemmCoord problem(task.npts, task.nbe, task.nbe);
    problem_sizes[i] = problem;

    cutlass::gemm::GemmCoord problem2(task.nbe, task.nbe, task.npts);
    sym2k_sizes[i] = problem2;

  }

  // Send to device
  device_backend_->copy_async( ntask, dmat_host.data(), 
    cutlass_stack.dmat_array_device, "send dmat array" );
  device_backend_->copy_async( ntask, zmat_host.data(), 
    cutlass_stack.zmat_array_device, "send zmat array" );
  device_backend_->copy_async( ntask, vmat_host.data(), 
    cutlass_stack.vmat_array_device, "send vmat array" );
  device_backend_->copy_async( ntask, bf_host.data(), 
    cutlass_stack.bf_array_device, "send bf array" );

  device_backend_->copy_async( ntask, problem_sizes.data(), 
    cutlass_stack.problem_sizes_device, "send problemsize array" );
  device_backend_->copy_async( ntask, sym2k_sizes.data(), 
    cutlass_stack.syr2k_sizes_device, "send problemsize array" );
  device_backend_->copy_async( ntask, ld64_dmat_host.data(), 
    cutlass_stack.ld64_dmat_array_device, "send ld dmat array" );
  device_backend_->copy_async( ntask, ld64_zmat_host.data(), 
    cutlass_stack.ld64_zmat_array_device, "send ld zmat array" );
  device_backend_->copy_async( ntask, ld64_vmat_host.data(), 
    cutlass_stack.ld64_vmat_array_device, "send ld vmat array" );
  device_backend_->copy_async( ntask, ld64_bf_host.data(), 
    cutlass_stack.ld64_bf_array_device, "send ld bf array" );

  device_backend_->master_queue_synchronize(); 

}
}
