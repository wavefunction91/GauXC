#include "scheme1_magma_base.hpp"
#include "buffer_adaptor.hpp"

namespace GauXC {


void AoSScheme1MAGMABase::Data::reset_allocations() {
  base_type::reset_allocations();
  magma_stack.reset();
}

size_t AoSScheme1MAGMABase::Data::get_static_mem_requirement() {
  return base_type::get_static_mem_requirement() + 8 * sizeof(int32_t); 
    // Extra elements in MAGMA dimension arrays
}

size_t AoSScheme1MAGMABase::Data::get_mem_req( integrator_term_tracker terms, 
  const host_task_type& task ) {

  
  size_t base_size = base_type::get_mem_req(terms, task);

  // TODO: There is probably a better way to check this
  required_term_storage reqt(terms);
  if( reqt.task_nbe_scr ) {
    base_size += 
      4*sizeof(double*) + // batch device pointers
      8*sizeof(int32_t);  // Dimensions + leading dimensions 
                          // (extra handled by get_static_mem_requirement)
  }
  return base_size;

}
AoSScheme1MAGMABase::Data::device_buffer_t 
  AoSScheme1MAGMABase::Data::allocate_dynamic_stack( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end, 
  device_buffer_t buf ){

  // Allocate base info on the stack
  buf = base_type::allocate_dynamic_stack( terms, task_begin, task_end,
    buf );

  required_term_storage reqt(terms);
  if( not reqt.task_nbe_scr ) return buf;

  // Allocate additional device memory 
  auto [ ptr, sz ] = buf;
  buffer_adaptor mem( ptr, sz );

  const auto ntask = std::distance( task_begin, task_end );
  magma_stack.dmat_array_device = mem.aligned_alloc<double*>( ntask, csl );
  magma_stack.vmat_array_device = mem.aligned_alloc<double*>( ntask, csl );
  magma_stack.zmat_array_device = mem.aligned_alloc<double*>( ntask, csl );
  magma_stack.bf_array_device   = mem.aligned_alloc<double*>( ntask, csl );

  magma_stack.m_array_device   = mem.aligned_alloc<int32_t>( ntask + 1, csl );
  magma_stack.n_array_device   = mem.aligned_alloc<int32_t>( ntask + 1, csl );
  magma_stack.k_array_device   = mem.aligned_alloc<int32_t>( ntask + 1, csl );
  magma_stack.ld_dmat_array_device = mem.aligned_alloc<int32_t>( ntask + 1, csl );
  magma_stack.ld_zmat_array_device = mem.aligned_alloc<int32_t>( ntask + 1, csl );
  magma_stack.ld_vmat_array_device = mem.aligned_alloc<int32_t>( ntask + 1, csl );
  magma_stack.ld_bf_array_device   = mem.aligned_alloc<int32_t>( ntask + 1, csl );

  // Update dynmem data for derived impls
  return device_buffer_t{ mem.stack(), mem.nleft() };
}

void AoSScheme1MAGMABase::Data::pack_and_send( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end,
  const BasisSetMap& basis_map ) {

  base_type::pack_and_send( terms, task_begin, task_end, basis_map );
  required_term_storage reqt(terms);
  if( not reqt.task_nbe_scr ) return;

  const auto ntask = std::distance( task_begin, task_end );
  std::vector<double*> dmat_host( ntask ), zmat_host( ntask ), bf_host( ntask ),
                       vmat_host( ntask );
  std::vector<int32_t> m_host( ntask ), n_host( ntask ), k_host( ntask ),
                       ld_dmat_host( ntask ), ld_zmat_host( ntask ), 
                       ld_vmat_host( ntask ), ld_bf_host( ntask );

  double* static_dmat = static_stack.dmat_device;
  const auto nbf = global_dims.nbf;

  // host_device_tasks should be populated by parent impl called at top
  for( auto i = 0; i < ntask; ++i ) {
    auto& task = host_device_tasks[i];
    zmat_host[i] = task.zmat;    ld_zmat_host[i] = task.npts;
    bf_host[i]   = task.bf;      ld_bf_host[i]   = task.npts;
    vmat_host[i] = task.nbe_scr; ld_vmat_host[i] = task.bfn_screening.nbe;
    if( task.bfn_screening.ncut > 1 ) {
      dmat_host[i]    = task.nbe_scr;
      ld_dmat_host[i] = task.bfn_screening.nbe;
    } else {
      dmat_host[i]    = static_dmat + task.bfn_screening.ibf_begin*(nbf+1);
      ld_dmat_host[i] = nbf;
    }

    m_host[i] = task.npts;
    n_host[i] = task.bfn_screening.nbe;
    k_host[i] = task.bfn_screening.nbe;
  }

  // Send to device
  device_backend_->copy_async( ntask, dmat_host.data(), 
    magma_stack.dmat_array_device, "send dmat array" );
  device_backend_->copy_async( ntask, zmat_host.data(), 
    magma_stack.zmat_array_device, "send zmat array" );
  device_backend_->copy_async( ntask, vmat_host.data(), 
    magma_stack.vmat_array_device, "send vmat array" );
  device_backend_->copy_async( ntask, bf_host.data(), 
    magma_stack.bf_array_device, "send bf array" );

  device_backend_->copy_async( ntask, m_host.data(), magma_stack.m_array_device,
    "send m array" );
  device_backend_->copy_async( ntask, n_host.data(), magma_stack.n_array_device,
    "send n array" );
  device_backend_->copy_async( ntask, k_host.data(), magma_stack.k_array_device,
    "send k array" );
  device_backend_->copy_async( ntask, ld_dmat_host.data(), 
    magma_stack.ld_dmat_array_device, "send ld dmat array" );
  device_backend_->copy_async( ntask, ld_zmat_host.data(), 
    magma_stack.ld_zmat_array_device, "send ld zmat array" );
  device_backend_->copy_async( ntask, ld_vmat_host.data(), 
    magma_stack.ld_vmat_array_device, "send ld vmat array" );
  device_backend_->copy_async( ntask, ld_bf_host.data(), 
    magma_stack.ld_bf_array_device, "send ld bf array" );
  device_backend_->master_queue_synchronize(); 

}
}
