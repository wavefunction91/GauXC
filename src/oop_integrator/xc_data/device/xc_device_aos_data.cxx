#include "xc_device_aos_data.hpp"
#include "buffer_adaptor.hpp"
#include "integrator_util/integrator_common.hpp"

namespace GauXC {

void XCDeviceAoSData::reset_allocations() {
  XCDeviceStackData::reset_allocations(); // Base implementation
  aos_stack.reset();
}

size_t XCDeviceAoSData::get_mem_req( integrator_term_tracker terms,
  const host_task_type& task, const BasisSetMap& basis_map ) {

  size_t base_size = XCDeviceStackData::get_mem_req(terms, task, basis_map);
  
  // Everything in AoS is not required for current implementations of
  // the weights kernel
  if( not terms.exc_vxc ) return base_size;



  const auto& points     = task.points;
  const auto& shell_list = task.shell_list;

  // Get packing size 
  const size_t submat_chunk_size = this->get_submat_chunk_size(global_dims.nbf,0);

  // Generate basis map
  // TODO: This should happen once and get stored
  auto [submat_cut, submat_block] = 
    gen_compressed_submat_map( basis_map, shell_list, global_dims.nbf, submat_chunk_size );


  // Dimensions
  const size_t npts     = points.size();
  const size_t nbe      = task.nbe;
  const size_t ncut     = submat_cut.size();
  const size_t nblock   = submat_block.size();
  const size_t nshells  = shell_list.size();

  // Collocation + derivatives
  // TODO: this is dependent on integrand
  const size_t mem_bf   = nbe * npts * sizeof(double);
  const size_t mem_dbfx = mem_bf;
  const size_t mem_dbfy = mem_bf;
  const size_t mem_dbfz = mem_bf;

  // LDA/GGA Z Matrix 
  // TODO: this is dependent on integrand
  const size_t mem_zmat_lda_gga = nbe * npts * sizeof(double);

  // nbe * nbe scratch
  const size_t mem_nbe_scr = nbe * nbe * sizeof(double);

  // Shell index packing 
  const size_t mem_shell_list = nshells * sizeof(size_t);
  const size_t mem_shell_offs = nshells * sizeof(size_t);
  const size_t mem_submat_cut = 3 * ncut * sizeof(int32_t);
  const size_t mem_submat_block = nblock * sizeof(int32_t);

  // Memroty associated with added a task to the indirection
  const size_t mem_task = sizeof(XCDeviceTask);


  return base_size + 
    ( mem_bf + mem_dbfx + mem_dbfy + mem_dbfz + mem_nbe_scr ) +
    ( mem_zmat_lda_gga )                                      +
    ( mem_shell_list + mem_shell_offs )                       +
    ( mem_submat_cut + mem_submat_block )                     +
    ( mem_task );
}



#if 0
XCDeviceAoSData::device_buffer_t XCDeviceAoSData::alloc_pack_and_send( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end, device_buffer_t buf,
  const BasisSetMap& basis_map ) {

  if( not device_backend_ ) throw std::runtime_error("Invalid Device Backend");

  // Allocate base info in the stack
  buf = XCDeviceStackData::alloc_pack_and_send( terms, task_begin, task_end, buf, 
    basis_map );

  // Reset AoS
  host_device_tasks.clear();

  // Current Stack
  auto [ ptr, sz ] = buf;
  buffer_adaptor mem( ptr, sz );

  // We only need AoS for XC integrations 
  if( terms.exc_vxc ) {

  // Host Packing Arrays
  std::vector< size_t > shell_list_pack;
  std::vector< size_t > shell_offs_pack;
  std::vector< std::array<int32_t, 3> > submat_cut_pack;
  std::vector< int32_t > submat_block_pack;

  // Contatenation utility
  auto concat_iterable = []( auto& a, const auto& b ) {
    a.insert( a.end(), b.begin(), b.end() );
  };

  const size_t submat_chunk_size = this->get_submat_chunk_size(global_dims.nbf,0);

  // Pack AoS data and construct indirections
  total_nbe_sq_task_batch   = 0;
  total_nbe_npts_task_batch = 0; 
  total_nshells_task_batch  = 0; 
  total_ncut_task_batch     = 0; 
  total_nblock_task_batch   = 0; 

  for( auto it = task_begin; it != task_end; ++it ) {

    const auto  iAtom       = it->iParent;
    const auto& points      = it->points;
    const auto& shell_list  = it->shell_list;
    const auto dist_nearest = it->dist_nearest;

    // Generate basis map
    // TODO: This should happen once and get stored
    auto [submat_cut, submat_block] = 
      gen_compressed_submat_map( basis_map, shell_list, global_dims.nbf, submat_chunk_size );

    // Dimensions
    const size_t ncut     = submat_cut.size();
    const size_t nblock   = submat_block.size();
    const size_t nshells  = shell_list.size();
    const size_t npts     = points.size();
    const auto nbe        = it->nbe;

    // Compute offsets
    // TODO: this should happen once and get stored
    std::vector< size_t > shell_offs( nshells );
    shell_offs.at(0) = 0;
    for( auto i = 1ul; i < nshells; ++i )
      shell_offs.at(i) = shell_offs.at(i-1) + 
                           basis_map.shell_size( shell_list.at(i-1) );

    // Pack Shell indexing
    concat_iterable( shell_list_pack, shell_list );
    concat_iterable( shell_offs_pack, shell_offs );
    concat_iterable( submat_cut_pack, submat_cut );
    concat_iterable( submat_block_pack, submat_block );

    // Increment running counters
    total_nbe_sq_task_batch   += nbe * nbe;
    total_nbe_npts_task_batch += nbe * npts;
    total_nshells_task_batch  += nshells;
    total_ncut_task_batch     += ncut;
    total_nblock_task_batch   += nblock;

    // Add task to device indirection
    host_device_tasks.emplace_back();

    // Populate indirection with dimensions
    host_device_tasks.back().nbe          = nbe;
    host_device_tasks.back().npts         = npts;
    host_device_tasks.back().ncut         = ncut;
    host_device_tasks.back().nblock       = nblock;
    host_device_tasks.back().nshells      = nshells;
    host_device_tasks.back().iParent      = iAtom;
    host_device_tasks.back().dist_nearest = dist_nearest;

  }


  // TODO: Print this if verbose
  //std::cout << "XCDeviceAoSData buf = " << ptr << ", " << sz << std::endl;


  // Device task indirection
  const size_t ntask = std::distance( task_begin, task_end );
  aos_stack.device_tasks = mem.aligned_alloc<XCDeviceTask>( ntask , csl);

  // Collocation 
  aos_stack.bf_eval_device     = mem.aligned_alloc<double>( total_nbe_npts_task_batch , csl);
  aos_stack.dbf_x_eval_device  = mem.aligned_alloc<double>( total_nbe_npts_task_batch , csl);
  aos_stack.dbf_y_eval_device  = mem.aligned_alloc<double>( total_nbe_npts_task_batch , csl);
  aos_stack.dbf_z_eval_device  = mem.aligned_alloc<double>( total_nbe_npts_task_batch , csl);

  // VXC Z Matrix
  aos_stack.zmat_vxc_lda_gga_device = mem.aligned_alloc<double>( total_nbe_npts_task_batch , csl);

  // Scratch buffer
  aos_stack.nbe_scr_device = mem.aligned_alloc<double>( total_nbe_sq_task_batch , csl);

  // AoS buffers
  aos_stack.shell_list_device   = mem.aligned_alloc<size_t>( total_nshells_task_batch , csl);
  aos_stack.shell_offs_device   = mem.aligned_alloc<size_t>( total_nshells_task_batch , csl);
  aos_stack.submat_cut_device   = mem.aligned_alloc<int32_t>( 3 * total_ncut_task_batch , csl);
  aos_stack.submat_block_device = mem.aligned_alloc<int32_t>( total_nblock_task_batch , csl);

  // Send AoS information early to overlap with indirection construction
  device_backend_->copy_async( shell_list_pack.size(), shell_list_pack.data(), 
    aos_stack.shell_list_device, "send_shell_list" );
  device_backend_->copy_async( shell_offs_pack.size(), shell_offs_pack.data(), 
    aos_stack.shell_offs_device, "send_shell_offs" );
  device_backend_->copy_async( 3 * submat_cut_pack.size(), 
    submat_cut_pack.data()->data(), aos_stack.submat_cut_device, "send_submat_cut"  ); 
  device_backend_->copy_async( submat_block_pack.size(), submat_block_pack.data(), 
    aos_stack.submat_block_device, "send_submat_block"  ); 

  // Construct full indirection
  {

  const size_t total_npts    = total_npts_task_batch * sizeof(double);
  buffer_adaptor points_mem ( base_stack.points_device,  3*total_npts );
  buffer_adaptor weights_mem( base_stack.weights_device,   total_npts );

  const size_t total_nshells = total_nshells_task_batch * sizeof(size_t);
  const size_t total_ncut    = total_ncut_task_batch * sizeof(int32_t);
  const size_t total_nblock  = total_nblock_task_batch * sizeof(int32_t);
  buffer_adaptor shell_list_mem( aos_stack.shell_list_device, total_nshells );
  buffer_adaptor shell_offs_mem( aos_stack.shell_offs_device, total_nshells );
  buffer_adaptor submat_cut_mem( aos_stack.submat_cut_device, 3*total_ncut  );
  buffer_adaptor submat_block_mem( aos_stack.submat_block_device, total_nblock);

  const size_t total_nbe_sq = total_nbe_sq_task_batch * sizeof(double);
  const size_t total_nbe_npts = total_nbe_npts_task_batch * sizeof(double);
  buffer_adaptor nbe_mem( aos_stack.nbe_scr_device, total_nbe_sq );
  buffer_adaptor zmat_mem( aos_stack.zmat_vxc_lda_gga_device, total_nbe_npts );

  buffer_adaptor bf_mem   ( aos_stack.bf_eval_device,    total_nbe_npts );
  buffer_adaptor dbf_x_mem( aos_stack.dbf_x_eval_device, total_nbe_npts );
  buffer_adaptor dbf_y_mem( aos_stack.dbf_y_eval_device, total_nbe_npts );
  buffer_adaptor dbf_z_mem( aos_stack.dbf_z_eval_device, total_nbe_npts );

  buffer_adaptor den_mem   ( base_stack.den_eval_device,   total_npts );
  buffer_adaptor dden_x_mem( base_stack.den_x_eval_device, total_npts );
  buffer_adaptor dden_y_mem( base_stack.den_y_eval_device, total_npts );
  buffer_adaptor dden_z_mem( base_stack.den_z_eval_device, total_npts );

  buffer_adaptor eps_mem( base_stack.eps_eval_device, total_npts );
  buffer_adaptor gamma_mem( base_stack.gamma_eval_device, total_npts );
  buffer_adaptor vrho_mem( base_stack.vrho_eval_device, total_npts );
  buffer_adaptor vgamma_mem( base_stack.vgamma_eval_device, total_npts );

  size_t i = 0;
  for( auto& task : host_device_tasks ) {
    const auto npts    = task.npts;
    const auto nbe     = task.nbe;
    const auto ncut    = task.ncut;
    const auto nblock  = task.nblock;
    const auto nshells = task.nshells;

    task.points       = points_mem .aligned_alloc<double>(3*npts, csl);
    task.weights      = weights_mem.aligned_alloc<double>(npts, csl); 
    task.shell_list   = shell_list_mem.aligned_alloc<size_t>( nshells , csl); 
    task.shell_offs   = shell_offs_mem.aligned_alloc<size_t>( nshells , csl); 
    task.submat_cut   = submat_cut_mem.aligned_alloc<int32_t>( 3*ncut , csl);
    task.submat_block = submat_block_mem.aligned_alloc<int32_t>(nblock, csl);

    task.nbe_scr = nbe_mem .aligned_alloc<double>( nbe * nbe  , csl);
    task.zmat    = zmat_mem.aligned_alloc<double>( nbe * npts , csl);

    task.bf   = bf_mem   .aligned_alloc<double>( nbe * npts , csl);
    task.dbfx = dbf_x_mem.aligned_alloc<double>( nbe * npts , csl);
    task.dbfy = dbf_y_mem.aligned_alloc<double>( nbe * npts , csl);
    task.dbfz = dbf_z_mem.aligned_alloc<double>( nbe * npts , csl);

    task.den    = den_mem   .aligned_alloc<double>(npts, csl);
    task.ddenx  = dden_x_mem.aligned_alloc<double>(npts, csl);
    task.ddeny  = dden_y_mem.aligned_alloc<double>(npts, csl);
    task.ddenz  = dden_z_mem.aligned_alloc<double>(npts, csl);

    task.eps    = eps_mem   .aligned_alloc<double>(npts, csl);
    task.gamma  = gamma_mem .aligned_alloc<double>(npts, csl);
    task.vrho   = vrho_mem  .aligned_alloc<double>(npts, csl);
    task.vgamma = vgamma_mem.aligned_alloc<double>(npts, csl);

    //std::cout << i << ", " << nshells << std::endl;

    ++i;
  } // Loop over device tasks

  } // Setup indirection

  } // end AoS setup (XC Integrands only)

  // Setup extra pieces to indirection which are algorithm specific
  // XXX: this also allocates additional memory
  device_buffer_t buf_left{ mem.stack(), mem.nleft() };
  buf_left = add_extra_to_indirection(terms, host_device_tasks, buf_left);


  // Send indirection (XC only)
  if( terms.exc_vxc ) {
    device_backend_->copy_async( host_device_tasks.size(), host_device_tasks.data(), 
      aos_stack.device_tasks, "send_tasks_device" );
  }


  // Synchronize on the copy stream to keep host vecs in scope
  device_backend_->master_queue_synchronize(); 

  // Update dynmem data for derived impls
  return buf_left;
}
#else


XCDeviceAoSData::device_buffer_t XCDeviceAoSData::allocate_dynamic_stack( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end, device_buffer_t buf,
  const BasisSetMap& basis_map ) {

  // Allocate base info in the stack
  buf = XCDeviceStackData::allocate_dynamic_stack( terms, task_begin, task_end, 
    buf, basis_map );

  // All data that currently resides in AoS is XC related and can be skipped
  // for weights
  if( not terms.exc_vxc ) return buf; 


  // Current Stack
  auto [ ptr, sz ] = buf;
  buffer_adaptor mem( ptr, sz );

  const size_t submat_chunk_size = this->get_submat_chunk_size(global_dims.nbf,0);


  // Get dimensions
  total_nbe_sq_task_batch   = 0;
  total_nbe_npts_task_batch = 0; 
  total_nshells_task_batch  = 0; 
  total_ncut_task_batch     = 0; 
  total_nblock_task_batch   = 0; 
  for( auto it = task_begin; it != task_end; ++it ) {

    const auto& points      = it->points;
    const auto& shell_list  = it->shell_list;

    // Generate basis map
    // TODO: This should happen once and get stored
    auto [submat_cut, submat_block] = 
      gen_compressed_submat_map( basis_map, shell_list, global_dims.nbf, 
        submat_chunk_size );

    const size_t ncut     = submat_cut.size();
    const size_t nblock   = submat_block.size();
    const size_t nshells  = shell_list.size();
    const size_t npts     = points.size();
    const auto nbe        = it->nbe;

    total_nbe_sq_task_batch   += nbe * nbe;
    total_nbe_npts_task_batch += nbe * npts;
    total_nshells_task_batch  += nshells;
    total_ncut_task_batch     += ncut;
    total_nblock_task_batch   += nblock;

  }

  // Device task indirection
  const size_t ntask = std::distance( task_begin, task_end );
  aos_stack.device_tasks = mem.aligned_alloc<XCDeviceTask>( ntask , csl);

  // Collocation 
  aos_stack.bf_eval_device     = mem.aligned_alloc<double>( total_nbe_npts_task_batch , csl);
  aos_stack.dbf_x_eval_device  = mem.aligned_alloc<double>( total_nbe_npts_task_batch , csl);
  aos_stack.dbf_y_eval_device  = mem.aligned_alloc<double>( total_nbe_npts_task_batch , csl);
  aos_stack.dbf_z_eval_device  = mem.aligned_alloc<double>( total_nbe_npts_task_batch , csl);

  // VXC Z Matrix
  aos_stack.zmat_vxc_lda_gga_device = mem.aligned_alloc<double>( total_nbe_npts_task_batch , csl);

  // Scratch buffer
  aos_stack.nbe_scr_device = mem.aligned_alloc<double>( total_nbe_sq_task_batch , csl);

  // AoS buffers
  aos_stack.shell_list_device   = mem.aligned_alloc<size_t>( total_nshells_task_batch , csl);
  aos_stack.shell_offs_device   = mem.aligned_alloc<size_t>( total_nshells_task_batch , csl);
  aos_stack.submat_cut_device   = mem.aligned_alloc<int32_t>( 3 * total_ncut_task_batch , csl);
  aos_stack.submat_block_device = mem.aligned_alloc<int32_t>( total_nblock_task_batch , csl);

  // Update dynmem data for derived impls
  return device_buffer_t{ mem.stack(), mem.nleft() };
}


void XCDeviceAoSData::pack_and_send( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end,
  const BasisSetMap& basis_map ) {

  // Pack and send base data
  XCDeviceStackData::pack_and_send( terms, task_begin, task_end, basis_map );

  if( not device_backend_ ) throw std::runtime_error("Invalid Device Backend");

  // All data that currently resides in AoS is XC related and can be skipped
  // for weights
  if( not terms.exc_vxc ) return; 

  // Reset AoS
  host_device_tasks.clear();

  // Host Packing Arrays
  std::vector< size_t > shell_list_pack;
  std::vector< size_t > shell_offs_pack;
  std::vector< std::array<int32_t, 3> > submat_cut_pack;
  std::vector< int32_t > submat_block_pack;

  // Contatenation utility
  auto concat_iterable = []( auto& a, const auto& b ) {
    a.insert( a.end(), b.begin(), b.end() );
  };

  const size_t submat_chunk_size = this->get_submat_chunk_size(global_dims.nbf,0);

  // Pack AoS data and construct indirections
  for( auto it = task_begin; it != task_end; ++it ) {

    const auto  iAtom       = it->iParent;
    const auto& points      = it->points;
    const auto& shell_list  = it->shell_list;
    const auto dist_nearest = it->dist_nearest;

    // Generate basis map
    // TODO: This should happen once and get stored
    auto [submat_cut, submat_block] = 
      gen_compressed_submat_map( basis_map, shell_list, global_dims.nbf, submat_chunk_size );

    // Dimensions
    const size_t ncut     = submat_cut.size();
    const size_t nblock   = submat_block.size();
    const size_t nshells  = shell_list.size();
    const size_t npts     = points.size();
    const auto nbe        = it->nbe;

    // Compute offsets
    // TODO: this should happen once and get stored
    std::vector< size_t > shell_offs( nshells );
    shell_offs.at(0) = 0;
    for( auto i = 1ul; i < nshells; ++i )
      shell_offs.at(i) = shell_offs.at(i-1) + 
                           basis_map.shell_size( shell_list.at(i-1) );

    // Pack Shell indexing
    concat_iterable( shell_list_pack, shell_list );
    concat_iterable( shell_offs_pack, shell_offs );
    concat_iterable( submat_cut_pack, submat_cut );
    concat_iterable( submat_block_pack, submat_block );

    // Add task to device indirection
    host_device_tasks.emplace_back();

    // Populate indirection with dimensions
    host_device_tasks.back().nbe          = nbe;
    host_device_tasks.back().npts         = npts;
    host_device_tasks.back().ncut         = ncut;
    host_device_tasks.back().nblock       = nblock;
    host_device_tasks.back().nshells      = nshells;
    host_device_tasks.back().iParent      = iAtom;
    host_device_tasks.back().dist_nearest = dist_nearest;

  }


  // TODO: Print this if verbose
  //std::cout << "XCDeviceAoSData buf = " << ptr << ", " << sz << std::endl;



  // Send AoS information early to overlap with indirection construction
  device_backend_->copy_async( shell_list_pack.size(), shell_list_pack.data(), 
    aos_stack.shell_list_device, "send_shell_list" );
  device_backend_->copy_async( shell_offs_pack.size(), shell_offs_pack.data(), 
    aos_stack.shell_offs_device, "send_shell_offs" );
  device_backend_->copy_async( 3 * submat_cut_pack.size(), 
    submat_cut_pack.data()->data(), aos_stack.submat_cut_device, "send_submat_cut"  ); 
  device_backend_->copy_async( submat_block_pack.size(), submat_block_pack.data(), 
    aos_stack.submat_block_device, "send_submat_block"  ); 

  // Construct full indirection
  {

  const size_t total_npts    = total_npts_task_batch * sizeof(double);
  buffer_adaptor points_mem ( base_stack.points_device,  3*total_npts );
  buffer_adaptor weights_mem( base_stack.weights_device,   total_npts );

  const size_t total_nshells = total_nshells_task_batch * sizeof(size_t);
  const size_t total_ncut    = total_ncut_task_batch * sizeof(int32_t);
  const size_t total_nblock  = total_nblock_task_batch * sizeof(int32_t);
  buffer_adaptor shell_list_mem( aos_stack.shell_list_device, total_nshells );
  buffer_adaptor shell_offs_mem( aos_stack.shell_offs_device, total_nshells );
  buffer_adaptor submat_cut_mem( aos_stack.submat_cut_device, 3*total_ncut  );
  buffer_adaptor submat_block_mem( aos_stack.submat_block_device, total_nblock);

  const size_t total_nbe_sq = total_nbe_sq_task_batch * sizeof(double);
  const size_t total_nbe_npts = total_nbe_npts_task_batch * sizeof(double);
  buffer_adaptor nbe_mem( aos_stack.nbe_scr_device, total_nbe_sq );
  buffer_adaptor zmat_mem( aos_stack.zmat_vxc_lda_gga_device, total_nbe_npts );

  buffer_adaptor bf_mem   ( aos_stack.bf_eval_device,    total_nbe_npts );
  buffer_adaptor dbf_x_mem( aos_stack.dbf_x_eval_device, total_nbe_npts );
  buffer_adaptor dbf_y_mem( aos_stack.dbf_y_eval_device, total_nbe_npts );
  buffer_adaptor dbf_z_mem( aos_stack.dbf_z_eval_device, total_nbe_npts );

  buffer_adaptor den_mem   ( base_stack.den_eval_device,   total_npts );
  buffer_adaptor dden_x_mem( base_stack.den_x_eval_device, total_npts );
  buffer_adaptor dden_y_mem( base_stack.den_y_eval_device, total_npts );
  buffer_adaptor dden_z_mem( base_stack.den_z_eval_device, total_npts );

  buffer_adaptor eps_mem( base_stack.eps_eval_device, total_npts );
  buffer_adaptor gamma_mem( base_stack.gamma_eval_device, total_npts );
  buffer_adaptor vrho_mem( base_stack.vrho_eval_device, total_npts );
  buffer_adaptor vgamma_mem( base_stack.vgamma_eval_device, total_npts );

  size_t i = 0;
  for( auto& task : host_device_tasks ) {
    const auto npts    = task.npts;
    const auto nbe     = task.nbe;
    const auto ncut    = task.ncut;
    const auto nblock  = task.nblock;
    const auto nshells = task.nshells;

    task.points       = points_mem .aligned_alloc<double>(3*npts, csl);
    task.weights      = weights_mem.aligned_alloc<double>(npts, csl); 
    task.shell_list   = shell_list_mem.aligned_alloc<size_t>( nshells , csl); 
    task.shell_offs   = shell_offs_mem.aligned_alloc<size_t>( nshells , csl); 
    task.submat_cut   = submat_cut_mem.aligned_alloc<int32_t>( 3*ncut , csl);
    task.submat_block = submat_block_mem.aligned_alloc<int32_t>(nblock, csl);

    task.nbe_scr = nbe_mem .aligned_alloc<double>( nbe * nbe  , csl);
    task.zmat    = zmat_mem.aligned_alloc<double>( nbe * npts , csl);

    task.bf   = bf_mem   .aligned_alloc<double>( nbe * npts , csl);
    task.dbfx = dbf_x_mem.aligned_alloc<double>( nbe * npts , csl);
    task.dbfy = dbf_y_mem.aligned_alloc<double>( nbe * npts , csl);
    task.dbfz = dbf_z_mem.aligned_alloc<double>( nbe * npts , csl);

    task.den    = den_mem   .aligned_alloc<double>(npts, csl);
    task.ddenx  = dden_x_mem.aligned_alloc<double>(npts, csl);
    task.ddeny  = dden_y_mem.aligned_alloc<double>(npts, csl);
    task.ddenz  = dden_z_mem.aligned_alloc<double>(npts, csl);

    task.eps    = eps_mem   .aligned_alloc<double>(npts, csl);
    task.gamma  = gamma_mem .aligned_alloc<double>(npts, csl);
    task.vrho   = vrho_mem  .aligned_alloc<double>(npts, csl);
    task.vgamma = vgamma_mem.aligned_alloc<double>(npts, csl);

    //std::cout << i << ", " << nshells << std::endl;

    ++i;
  } // Loop over device tasks

  } // Setup indirection


  // Setup extra pieces to indirection which are algorithm specific
  add_extra_to_indirection(terms, host_device_tasks);

  // Send indirection 
  device_backend_->copy_async( host_device_tasks.size(), host_device_tasks.data(), 
    aos_stack.device_tasks, "send_tasks_device" );

  // Synchronize on the copy stream to keep host vecs in scope
  device_backend_->master_queue_synchronize(); 


}
#endif

}
