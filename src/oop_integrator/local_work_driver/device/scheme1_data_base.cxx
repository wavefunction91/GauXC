#include "scheme1_data_base.hpp"
#include "buffer_adaptor.hpp"

namespace GauXC {

Scheme1DataBase::~Scheme1DataBase() noexcept = default;

Scheme1DataBase::Scheme1DataBase(std::unique_ptr<DeviceBackend>&& ptr, 
  bool batch_l3_blas) : XCDeviceAoSData(std::move(ptr)) {

  if( batch_l3_blas && device_backend_ ) 
    device_backend_->create_blas_queue_pool(4);

}

size_t Scheme1DataBase::get_static_mem_requirement() {
  return 0;
}

size_t Scheme1DataBase::get_mem_req( const host_task_type& task,
  const BasisSetMap& basis_map ) {

  const auto ldatoms = get_ldatoms();
  const auto mem_dist_scr = ldatoms * task.npts;
  const auto mem_dist_ner = task.npts;
  const auto mem_iparent  = task.npts;

  size_t base_size = base_type::get_mem_req(task, basis_map);
  return base_size + 
    (mem_dist_scr + mem_dist_ner) * sizeof(double) + 
    mem_iparent * sizeof(int32_t);
}



Scheme1DataBase::device_buffer_t Scheme1DataBase::add_extra_to_indirection( 
  std::vector<XCDeviceTask>& tasks, device_buffer_t buf ) {

  if( not device_backend_ ) throw std::runtime_error("Invalid Device Backend");
  
  // Host Packing Arrays
  std::vector< int32_t > iparent_pack;
  std::vector< double >  dist_nearest_pack;

  // Allocate additional device memory 
  auto [ ptr, sz ] = buf;
  buffer_adaptor mem( ptr, sz );


  const auto ldatoms = get_ldatoms();
  dist_scratch_device = 
    mem.aligned_alloc<double>( ldatoms * total_npts_task_batch, sizeof(double2) );
  dist_nearest_device = mem.aligned_alloc<double>( total_npts_task_batch );
  iparent_device = mem.aligned_alloc<int32_t>( total_npts_task_batch );

  double* dist_scratch_ptr = dist_scratch_device;
  // Pack additional host data and send
  for( auto& task : tasks ) {
    iparent_pack.insert( iparent_pack.end(), task.npts, task.iParent );
    dist_nearest_pack.insert( dist_nearest_pack.end(), task.npts, 
      task.dist_nearest );

    // Extra indirection for dist scratch
    task.dist_scratch  = dist_scratch_ptr;
    dist_scratch_ptr   += ldatoms * task.npts;
  }

  device_backend_->copy_async( iparent_pack.size(), iparent_pack.data(), 
              iparent_device, "send iparent"  );
  device_backend_->copy_async( dist_nearest_pack.size(), dist_nearest_pack.data(), 
              dist_nearest_device, "send dist_nearest" );

  device_backend_->master_queue_synchronize(); 

  return device_buffer_t{ mem.stack(), mem.nleft() };
}

}
