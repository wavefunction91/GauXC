#include "xc_device_stack_data.hpp"
#include "buffer_adaptor.hpp"

namespace GauXC {

XCDeviceStackData::XCDeviceStackData( std::unique_ptr<DeviceBackend>&& ptr ) :
  device_backend_(std::move(ptr)) { 

  // Allocate Device memory
  if( device_backend_ ) {
    auto avail = device_backend_->get_available_mem();
    std::tie( device_ptr, devmem_sz ) = 
      device_backend_->allocate_device_buffer(0.9 * avail);
    reset_allocations();
  }

}





XCDeviceStackData::~XCDeviceStackData() noexcept {
  // Free memory if allocated
  if( device_backend_ and devmem_sz and device_ptr )
    device_backend_->free_device_buffer(device_ptr);
}


double* XCDeviceStackData::vxc_device_data() { return static_stack.vxc_device; }
double* XCDeviceStackData::exc_device_data() { return static_stack.exc_device; }
double* XCDeviceStackData::nel_device_data() { return static_stack.nel_device; }
device_queue XCDeviceStackData::queue() { 
  if( not device_backend_ ) GAUXC_GENERIC_EXCEPTION("Invalid Device Backend");
  return device_backend_->queue();
}




void XCDeviceStackData::reset_allocations() {
  dynmem_ptr = device_ptr;
  dynmem_sz  = devmem_sz;
  allocated_terms.reset();
  static_stack.reset();
  base_stack.reset();
}

void XCDeviceStackData::allocate_static_data_weights( int32_t natoms ) {

  if( allocated_terms.weights ) 
    GAUXC_GENERIC_EXCEPTION("Attempting to reallocate Stack Weights");

  // Save state
  global_dims.natoms  = natoms;

  // Allocate static memory with proper alignment
  buffer_adaptor mem( dynmem_ptr, dynmem_sz );

  static_stack.coords_device = mem.aligned_alloc<double>( 3 * natoms, csl );

  // Allow for RAB to be strided and properly aligned
  const auto ldatoms   = get_ldatoms();
  const auto rab_align = get_rab_align();
  static_stack.rab_device = mem.aligned_alloc<double>( natoms * ldatoms, rab_align, csl );

  // Get current stack location
  dynmem_ptr = mem.stack();
  dynmem_sz  = mem.nleft(); 

  allocated_terms.weights = true;
}

void XCDeviceStackData::allocate_static_data_exc_vxc( int32_t nbf, int32_t nshells ) {

  if( allocated_terms.exc_vxc ) 
    GAUXC_GENERIC_EXCEPTION("Attempting to reallocate Stack EXC VXC");

  // Save state
  global_dims.nshells = nshells;
  global_dims.nbf     = nbf; 

  // Allocate static memory with proper alignment
  buffer_adaptor mem( dynmem_ptr, dynmem_sz );

  static_stack.shells_device     = mem.aligned_alloc<Shell<double>>( nshells , csl);
  static_stack.exc_device        = mem.aligned_alloc<double>( 1 , csl);
  static_stack.nel_device        = mem.aligned_alloc<double>( 1 , csl);
  static_stack.acc_scr_device    = mem.aligned_alloc<double>( 1 , csl);

  static_stack.vxc_device  = mem.aligned_alloc<double>( nbf * nbf , csl);
  static_stack.dmat_device = mem.aligned_alloc<double>( nbf * nbf , csl);

  // Get current stack location
  dynmem_ptr = mem.stack();
  dynmem_sz  = mem.nleft(); 

  allocated_terms.exc_vxc = true;
}







void XCDeviceStackData::send_static_data_weights( const Molecule& mol, const MolMeta& meta ) {

  if( not allocated_terms.weights ) 
    GAUXC_GENERIC_EXCEPTION("Weights Not Stack Allocated");

  const auto natoms = global_dims.natoms;
  if( not device_backend_ ) GAUXC_GENERIC_EXCEPTION("Invalid Device Backend");

  // Copy Atomic Coordinates
  std::vector<double> coords( 3*natoms );
  for( auto i = 0ul; i < natoms; ++i ) {
    coords[ 3*i + 0 ] = mol[i].x;
    coords[ 3*i + 1 ] = mol[i].y;
    coords[ 3*i + 2 ] = mol[i].z;
  }
  device_backend_->copy_async( 3*natoms, coords.data(), static_stack.coords_device, 
    "Coords H2D" );

  // Invert and send RAB
  const auto ldatoms = get_ldatoms();
  std::vector<double> rab_inv(natoms*natoms);
  for( auto i = 0ul; i < (natoms*natoms); ++i) rab_inv[i] = 1./meta.rab().data()[i];
  device_backend_->copy_async_2d( natoms, natoms, rab_inv.data(), natoms,
    static_stack.rab_device, ldatoms, "RAB H2D" );

  device_backend_->master_queue_synchronize(); 
}

void XCDeviceStackData::send_static_data_exc_vxc( const double* P, int32_t ldp,
  const BasisSet<double>& basis ) {

  if( not allocated_terms.exc_vxc ) 
    GAUXC_GENERIC_EXCEPTION("EXC_VXC Not Stack Allocated");

  const auto nbf    = global_dims.nbf;
  if( ldp != (int)nbf ) GAUXC_GENERIC_EXCEPTION("LDP must bf NBF");
  if( not device_backend_ ) GAUXC_GENERIC_EXCEPTION("Invalid Device Backend");

  // Copy Density
  device_backend_->copy_async( nbf*nbf, P, static_stack.dmat_device, "P H2D" );

  // Copy Basis Set
  device_backend_->copy_async( basis.nshells(), basis.data(), static_stack.shells_device,
    "Shells H2D" );


  device_backend_->master_queue_synchronize(); 
}






void XCDeviceStackData::zero_integrands() {

  if( not device_backend_ ) GAUXC_GENERIC_EXCEPTION("Invalid Device Backend");

  const auto nbf = global_dims.nbf;
  device_backend_->set_zero( nbf*nbf, static_stack.vxc_device, "VXC Zero" );
  device_backend_->set_zero( 1,       static_stack.nel_device, "NEL Zero" );
  device_backend_->set_zero( 1,       static_stack.exc_device, "EXC Zero" );

}





void XCDeviceStackData::retrieve_xc_integrands( double* EXC, double* N_EL,
  double* VXC, int32_t ldvxc ) {

  const auto nbf = global_dims.nbf;
  if( ldvxc != (int)nbf ) GAUXC_GENERIC_EXCEPTION("LDVXC must bf NBF");
  if( not device_backend_ ) GAUXC_GENERIC_EXCEPTION("Invalid Device Backend");
  
  device_backend_->copy_async( nbf*nbf, static_stack.vxc_device, VXC,  "VXC D2H" );
  device_backend_->copy_async( 1,       static_stack.nel_device, N_EL, "NEL D2H" );
  device_backend_->copy_async( 1,       static_stack.exc_device, EXC,  "EXC D2H" );

}





XCDeviceStackData::host_task_iterator XCDeviceStackData::generate_buffers(
  integrator_term_tracker terms,
  const BasisSetMap& basis_map,
  host_task_iterator task_begin,
  host_task_iterator task_end
) {

  if( get_static_mem_requirement() > dynmem_sz )
    GAUXC_GENERIC_EXCEPTION("Insufficient memory to even start!");

  size_t mem_left = dynmem_sz - get_static_mem_requirement();

  // Determine the number of batches that will fit into device memory
  host_task_iterator task_it = task_begin;
  while( task_it != task_end ) {

    // Get memory requirement for batch
    size_t mem_req_batch = get_mem_req( terms, *task_it );

    // Break out of loop if we can't allocate for this batch
    if( mem_req_batch > mem_left ) break;

    // Update remaining memory and increment task iterator
    mem_left -= mem_req_batch;
    task_it++;

  }

  // TODO: print this if verbose
  //std::cout << "XCDeviceStackData will allocate for " << std::distance(task_begin, task_it) << " Tasks" << std::endl;

  // Pack host data and send to device
  allocate_dynamic_stack( terms, task_begin, task_it,
    device_buffer_t{dynmem_ptr, dynmem_sz} );

  pack_and_send( terms, task_begin, task_it, basis_map );

  return task_it;
}





size_t XCDeviceStackData::get_mem_req( 
  integrator_term_tracker terms,
  const host_task_type& task
) {

  const auto& points = task.points;
  const size_t npts  = points.size();

  // Grid
  const auto mem_points = 3 * npts * sizeof(double);
  const auto mem_weights = npts * sizeof(double);

  // All terms require grid
  size_t mem_req = mem_points + mem_weights;

  // XC Specific terms
  if( terms.exc_vxc or terms.exc_grad ) {

    if( terms.xc_approx == _UNDEFINED ) GAUXC_GENERIC_EXCEPTION("NO XC APPROX SET");
    const bool is_lda = terms.xc_approx == LDA;
    const bool is_gga = terms.xc_approx == GGA;
  
    const bool need_grad = is_gga or terms.exc_grad;

    // U variables
    const auto mem_den      = npts * sizeof(double);
    const auto mem_den_grad = need_grad ? 3*mem_den : 0;

    // V variables
    const auto mem_gamma = (!is_lda) ? npts * sizeof(double) : 0;

    // XC output
    const auto mem_eps    = npts * sizeof(double);
    const auto mem_vrho   = npts * sizeof(double);
    const auto mem_vgamma = (!is_lda) ? npts * sizeof(double) : 0;

    mem_req += mem_den + mem_den_grad + mem_gamma + mem_eps + mem_vrho + mem_vgamma;
   
  }

  return mem_req;
}









XCDeviceStackData::device_buffer_t XCDeviceStackData::allocate_dynamic_stack( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end, 
  device_buffer_t buf ) {


  // Get total npts
  total_npts_task_batch = std::accumulate( task_begin, task_end, 0ul,
    [](const auto& a, const auto& b){ return a + b.points.size(); } );

  // Allocate device memory
  auto [ ptr, sz ] = buf;
  buffer_adaptor mem( ptr, sz );

  // Grid
  //base_stack.points_device = 
  //  mem.aligned_alloc<double>( 3 * total_npts_task_batch , csl);
  base_stack.points_x_device = mem.aligned_alloc<double>( total_npts_task_batch, 256, csl);
  base_stack.points_y_device = mem.aligned_alloc<double>( total_npts_task_batch, 256, csl);
  base_stack.points_z_device = mem.aligned_alloc<double>( total_npts_task_batch, 256, csl);

  base_stack.weights_device = 
    mem.aligned_alloc<double>( total_npts_task_batch , csl);

  // XC Specific terms
  if( terms.exc_vxc or terms.exc_grad ) {

    if( terms.xc_approx == _UNDEFINED ) GAUXC_GENERIC_EXCEPTION("NO XC APPROX SET");
    const bool is_lda = terms.xc_approx == LDA;
    const bool is_gga = terms.xc_approx == GGA;
  
    const bool need_grad = is_gga or terms.exc_grad;

    // U Variables
    base_stack.den_eval_device   = 
      mem.aligned_alloc<double>( total_npts_task_batch, 256, csl);

    if( need_grad ) {
      base_stack.den_x_eval_device = 
        mem.aligned_alloc<double>( total_npts_task_batch, 256, csl);
      base_stack.den_y_eval_device = 
        mem.aligned_alloc<double>( total_npts_task_batch, 256, csl);
      base_stack.den_z_eval_device = 
        mem.aligned_alloc<double>( total_npts_task_batch, 256, csl);
    }

    // V Variables
    if( !is_lda ) {
      base_stack.gamma_eval_device = 
        mem.aligned_alloc<double>( total_npts_task_batch, csl);
    }

    // XC output
    base_stack.eps_eval_device    = 
      mem.aligned_alloc<double>( total_npts_task_batch, csl);
    base_stack.vrho_eval_device   = 
      mem.aligned_alloc<double>( total_npts_task_batch, csl);
    if( !is_lda ) {
      base_stack.vgamma_eval_device = 
        mem.aligned_alloc<double>( total_npts_task_batch, csl);
    }
  }

  // Update dynmem data for derived impls
  return device_buffer_t{ mem.stack(), mem.nleft() };
}

void XCDeviceStackData::pack_and_send( integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end, const BasisSetMap& ) {

  if( not device_backend_ ) GAUXC_GENERIC_EXCEPTION("Invalid Device Backend");

  // Host data packing arrays
  //std::vector< std::array<double,3> > points_pack;
  std::vector<double> points_x_pack, points_y_pack, points_z_pack;
  std::vector< double > weights_pack;

  // Contatenation utility
  auto concat_iterable = []( auto& a, const auto& b ) {
    a.insert( a.end(), b.begin(), b.end() );
  };

  // Pack points / weights
  for( auto it = task_begin; it != task_end; ++it ) {

    const auto& points  = it->points;
    const auto& weights = it->weights;

    //concat_iterable( points_pack,  points  );
    std::vector<double> pts_x, pts_y, pts_z;
    for( auto pt : points ) {
      pts_x.emplace_back( pt[0] );
      pts_y.emplace_back( pt[1] );
      pts_z.emplace_back( pt[2] );
    }
    concat_iterable( points_x_pack, pts_x );
    concat_iterable( points_y_pack, pts_y );
    concat_iterable( points_z_pack, pts_z );

    concat_iterable( weights_pack, weights );
    
  } // Loop over tasks


  // Send grid data
  //device_backend_->copy_async( 3*points_pack.size(), points_pack.data()->data(),
  //            base_stack.points_device, "send points buffer" );
  device_backend_->copy_async( points_x_pack.size(), points_x_pack.data(),
              base_stack.points_x_device, "send points_x buffer" );
  device_backend_->copy_async( points_y_pack.size(), points_y_pack.data(),
              base_stack.points_y_device, "send points_y buffer" );
  device_backend_->copy_async( points_z_pack.size(), points_z_pack.data(),
              base_stack.points_z_device, "send points_z buffer" );
  device_backend_->copy_async( weights_pack.size(), weights_pack.data(),
              base_stack.weights_device, "send weights buffer" );


  // Synchronize on the copy stream to keep host vecs in scope
  device_backend_->master_queue_synchronize(); 

}


void XCDeviceStackData::copy_weights_to_tasks( host_task_iterator task_begin, host_task_iterator task_end ) {

  if( not device_backend_ ) GAUXC_GENERIC_EXCEPTION("Invalid Device Backend");

  // Sanity check that npts is consistent
  size_t local_npts = std::accumulate( task_begin, task_end, 0ul, 
    []( const auto& a, const auto& b ) { return a + b.points.size(); } );

  if( local_npts != total_npts_task_batch )
    GAUXC_GENERIC_EXCEPTION("NPTS MISMATCH");

  // Copy weights into contiguous host data
  std::vector<double> weights_host(local_npts);
  device_backend_->copy_async( local_npts, base_stack.weights_device, weights_host.data(),
    "Weights D2H" );
  device_backend_->master_queue_synchronize(); 

  // Place into host memory 
  auto* weights_ptr = weights_host.data();
  for( auto it = task_begin; it != task_end; ++it ) {
    const auto npts = it->points.size();
    std::copy_n( weights_ptr, npts, it->weights.data() );
    weights_ptr += npts;
  }

}


}
