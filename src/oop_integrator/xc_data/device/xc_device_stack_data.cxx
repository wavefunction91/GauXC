#include "xc_device_stack_data.hpp"
#include "buffer_adaptor.hpp"

namespace GauXC {

void XCDeviceStackData::allocate_static_data( int32_t _natoms, int32_t _nbf,
  int32_t _nshells ) {

  // Save state
  nshells = _nshells;
  nbf     = _nbf; 
  natoms  = _natoms;

  // Allocate static memory with proper alignment
  buffer_adaptor mem( device_ptr, devmem_sz );

  shells_device     = mem.aligned_alloc<Shell<double>>( nshells );
  exc_device        = mem.aligned_alloc<double>( 1 );
  nel_device        = mem.aligned_alloc<double>( 1 );
  acc_scr_device    = mem.aligned_alloc<double>( 1 );
  coords_device     = mem.aligned_alloc<double>( 3 * natoms );

  vxc_device  = mem.aligned_alloc<double>( nbf * nbf );
  dmat_device = mem.aligned_alloc<double>( nbf * nbf );

  // Get current stack location
  dynmem_ptr = mem.stack();
  dynmem_sz  = mem.nleft(); 

  // XXX: RAB Device Storage can be strided for performance, must be handled
  // in derived implementation
  allocate_rab();

}

void XCDeviceStackData::send_static_data( const double* P, int32_t ldp,
  const BasisSet<double>& basis, const Molecule& mol,
  const MolMeta& meta ) {

  if( ldp != nbf ) throw std::runtime_error("LDP must bf NBF");

  // Copy Density
  copy_async( nbf*nbf, P, dmat_device, "P H2D" );

  // Copy Basis Set
  copy_async( basis.nshells(), basis.data(), shells_device,
    "Shells H2D" );

  // Copy Atomic Coordinates
  std::vector<double> coords( 3*natoms );
  for( auto i = 0ul; i < natoms; ++i ) {
    coords[ 3*i + 0 ] = mol[i].x;
    coords[ 3*i + 1 ] = mol[i].y;
    coords[ 3*i + 2 ] = mol[i].z;
  }
  copy_async( 3*natoms, coords.data(), coords_device, "Coords H2D" );

  // Implementation dependent send of RAB
  send_rab(meta);

  master_queue_synchronize(); 
}


void XCDeviceStackData::zero_integrands() {

  set_zero( nbf*nbf, vxc_device, "VXC Zero" );
  set_zero( 1,       nel_device, "NEL Zero" );
  set_zero( 1,       exc_device, "EXC Zero" );

}

void XCDeviceStackData::retrieve_xc_integrands( double* EXC, double* N_EL,
  double* VXC, int32_t ldvxc ) {

  if( ldvxc != nbf ) throw std::runtime_error("LDVXC must bf NBF");
  
  copy_async( nbf*nbf, vxc_device, VXC,  "VXC D2H" );
  copy_async( 1,       nel_device, N_EL, "NEL D2H" );
  copy_async( 1,       exc_device, EXC,  "EXC D2H" );

}





XCDeviceStackData::host_task_iterator XCDeviceStackData::generate_buffers(
  const BasisSetMap& basis_map,
  host_task_iterator task_begin,
  host_task_iterator task_end
) {

  size_t mem_left = dynmem_sz - get_static_mem_requirement();

  // Determine the number of batches that will fit into device memory
  host_task_iterator task_it = task_begin;
  while( task_it != task_end ) {

    // Get memory requirement for batch
    size_t mem_req_batch = get_mem_req( *task_it, basis_map );

    // Break out of loop if we can't allocate for this batch
    if( mem_req_batch > mem_left ) break;

    // Update remaining memory and increment task iterator
    mem_left -= mem_req_batch;
    task_it++;

  }


  // Pack host data and send to device
  alloc_pack_and_send( task_begin, task_it, 
    device_buffer_t{dynmem_ptr, dynmem_sz}, basis_map );


  return task_it;
}




size_t XCDeviceStackData::get_mem_req( const host_task_type& task,
  const BasisSetMap& ) {

  const auto& points = task.points;
  const size_t npts  = points.size();

  // Grid
  const auto mem_points = 3 * npts;
  const auto mem_weights = npts;

  // U variables
  const auto mem_den     = npts;
  const auto mem_den_x   = npts;
  const auto mem_den_y   = npts;
  const auto mem_den_z   = npts;

  // V variables
  const auto mem_gamma = npts;

  // XC output
  const auto mem_eps    = npts;
  const auto mem_vrho   = npts;
  const auto mem_vgamma = npts;

  // Everything is double here
  return (mem_points + mem_weights + 
          mem_den + mem_den_x + mem_den_y + mem_den_z +
          mem_gamma +
          mem_eps + mem_vrho + mem_vgamma) * sizeof(double);
}



XCDeviceStackData::device_buffer_t XCDeviceStackData::alloc_pack_and_send( 
  host_task_iterator task_begin, host_task_iterator task_end, device_buffer_t buf,
  const BasisSetMap& ) {

  // Host data packing arrays
  std::vector< std::array<double,3> > points_pack;
  std::vector< double > weights_pack;

  // Contatenation utility
  auto concat_iterable = []( auto& a, const auto& b ) {
    a.insert( a.end(), b.begin(), b.end() );
  };

  // Pack points / weights
  for( auto it = task_begin; it != task_end; ++it ) {

    const auto& points  = it->points;
    const auto& weights = it->weights;

    concat_iterable( points_pack,  points  );
    concat_iterable( weights_pack, weights );
    
  } // Loop over tasks

  total_npts_task_batch = weights_pack.size();

  // Allocate device memory
  auto [ ptr, sz ] = buf;
  buffer_adaptor mem( ptr, sz );

  // Grid
  points_device  = mem.aligned_alloc<double>( 3 * total_npts_task_batch );
  weights_device = mem.aligned_alloc<double>( total_npts_task_batch );

  // U Variables
  den_eval_device   = mem.aligned_alloc<double>( total_npts_task_batch );
  den_x_eval_device = mem.aligned_alloc<double>( total_npts_task_batch );
  den_y_eval_device = mem.aligned_alloc<double>( total_npts_task_batch );
  den_z_eval_device = mem.aligned_alloc<double>( total_npts_task_batch );

  // V Variables
  gamma_eval_device = mem.aligned_alloc<double>( total_npts_task_batch );

  // XC output
  eps_eval_device    = mem.aligned_alloc<double>( total_npts_task_batch );
  vrho_eval_device   = mem.aligned_alloc<double>( total_npts_task_batch );
  vgamma_eval_device = mem.aligned_alloc<double>( total_npts_task_batch );

  // Send grid data
  copy_async( 3*points_pack.size(), points_pack.data()->data(),
              points_device, "send points buffer" );
  copy_async( weights_pack.size(), weights_pack.data(),
              weights_device, "send weights buffer" );

  // Synchronize on the copy stream to keep host vecs in scope
  master_queue_synchronize(); 

  // Update dynmem data for derived impls
  return device_buffer_t{ mem.stack(), mem.nleft() };
}



}
