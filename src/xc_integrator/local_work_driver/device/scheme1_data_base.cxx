#include "scheme1_data_base.hpp"
#include "buffer_adaptor.hpp"

namespace GauXC {

Scheme1DataBase::~Scheme1DataBase() noexcept = default;

Scheme1DataBase::Scheme1DataBase(std::unique_ptr<DeviceBackend>&& ptr, 
  bool batch_l3_blas) : XCDeviceAoSData(std::move(ptr)) {

  if( batch_l3_blas && device_backend_ ) 
    device_backend_->create_blas_queue_pool(4);

}

void Scheme1DataBase::reset_allocations() {
  XCDeviceAoSData::reset_allocations();
  scheme1_stack.reset();
}

size_t Scheme1DataBase::get_static_mem_requirement() {
  return 0;
}









size_t Scheme1DataBase::get_mem_req( integrator_term_tracker terms, 
  const host_task_type& task ) {

  // All local memory is weights related
  size_t base_size = base_type::get_mem_req(terms, task);
  if( not terms.weights ) return base_size;

  const auto ldatoms = get_ldatoms();
  const auto mem_dist_scr = ldatoms * task.npts;
  const auto mem_dist_ner = task.npts;
  const auto mem_iparent  = task.npts;

  return base_size + 
    (mem_dist_scr + mem_dist_ner) * sizeof(double) + 
    mem_iparent * sizeof(int32_t);
}










Scheme1DataBase::device_buffer_t Scheme1DataBase::allocate_dynamic_stack( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end, 
  device_buffer_t buf ){

  // Allocate base info on the stack
  buf = XCDeviceAoSData::allocate_dynamic_stack( terms, task_begin, task_end,
    buf );

  // All local memory is weights related
  if( not terms.weights ) return buf;

  // Allocate additional device memory 
  auto [ ptr, sz ] = buf;
  buffer_adaptor mem( ptr, sz );

  const auto ldatoms = get_ldatoms();
  scheme1_stack.dist_scratch_device = 
    mem.aligned_alloc<double>( ldatoms * total_npts_task_batch, sizeof(double2) );
  scheme1_stack.dist_nearest_device = mem.aligned_alloc<double>( total_npts_task_batch );
  scheme1_stack.iparent_device = mem.aligned_alloc<int32_t>( total_npts_task_batch );

  // Update dynmem data for derived impls
  return device_buffer_t{ mem.stack(), mem.nleft() };
}









void Scheme1DataBase::pack_and_send( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end,
  const BasisSetMap& basis_map ) {

  // Pack and send base data
  XCDeviceAoSData::pack_and_send( terms, task_begin, task_end, basis_map );

  // All local memory is weights related
  if( not terms.weights ) return;

  // Host Packing Arrays
  std::vector< int32_t > iparent_pack;
  std::vector< double >  dist_nearest_pack;

  iparent_pack.reserve( total_npts_task_batch );
  dist_nearest_pack.reserve( total_npts_task_batch );

  // Pack additional host data and send
  for( auto it = task_begin; it != task_end; ++it ) {
    iparent_pack.insert( iparent_pack.end(), it->points.size(), it->iParent );
    dist_nearest_pack.insert( dist_nearest_pack.end(), it->points.size(), 
      it->dist_nearest );
  }

  device_backend_->copy_async( iparent_pack.size(), iparent_pack.data(), 
              scheme1_stack.iparent_device, "send iparent"  );
  device_backend_->copy_async( dist_nearest_pack.size(), dist_nearest_pack.data(), 
              scheme1_stack.dist_nearest_device, "send dist_nearest" );

  device_backend_->master_queue_synchronize(); 
}








void Scheme1DataBase::add_extra_to_indirection( 
  integrator_term_tracker terms, std::vector<XCDeviceTask>& tasks  ) {

  // All local memory is weights related
  if( not terms.weights ) return;

  const auto ldatoms = get_ldatoms();
  buffer_adaptor dist_scratch_mem( scheme1_stack.dist_scratch_device, 
    ldatoms * total_npts_task_batch * sizeof(double) );

  // Extra indirection for dist scratch
  for( auto& task : tasks ) {
    task.dist_scratch  = dist_scratch_mem.aligned_alloc<double>( 
      ldatoms * task.npts, sizeof(double2) );
  }

}

}
