#include "xc_device_aos_data.hpp"
#include "buffer_adaptor.hpp"
#include "integrator_util/integrator_common.hpp"

namespace GauXC {

void XCDeviceAoSData::reset_allocations() {
  XCDeviceStackData::reset_allocations(); // Base implementation
  aos_stack.reset();
}

size_t XCDeviceAoSData::get_mem_req( integrator_term_tracker terms,
  const host_task_type& task ) {

  size_t base_size = XCDeviceStackData::get_mem_req(terms, task);
  
  // Everything in AoS is not required for current implementations of
  // the weights kernel
  if( not terms.exc_vxc ) return base_size;

  const auto& points     = task.points;
  const auto& shell_list = task.shell_list;
  const auto& submat_cut = task.submat_map;
  const auto& submat_block = task.submat_block;
  if( !submat_cut.size() or !submat_block.size() )
    throw std::runtime_error("Must Populate Submat Maps");



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

  // Shell -> Task maps
  const size_t mem_shell_to_task_idx  = nshells * sizeof(int32_t);
  const size_t mem_shell_to_task_off  = nshells * sizeof(int32_t);
  const size_t mem_shell_to_task_data = nshells * sizeof(ShellToTaskDevice);

  // Memroty associated with added a task to the indirection
  const size_t mem_task = sizeof(XCDeviceTask);


  return base_size + 
    ( mem_bf + mem_dbfx + mem_dbfy + mem_dbfz + mem_nbe_scr ) +
    ( mem_zmat_lda_gga )                                      +
    ( mem_shell_list + mem_shell_offs )                       +
    ( mem_submat_cut + mem_submat_block )                     +
    ( mem_shell_to_task_idx + mem_shell_to_task_off +
      mem_shell_to_task_data )                                +
    ( mem_task );
}





XCDeviceAoSData::device_buffer_t XCDeviceAoSData::allocate_dynamic_stack( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end, 
  device_buffer_t buf ) {

  // Allocate base info in the stack
  buf = XCDeviceStackData::allocate_dynamic_stack( terms, task_begin, task_end, 
    buf );

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
    const auto& submat_cut = it->submat_map;
    const auto& submat_block = it->submat_block;
    if( !submat_cut.size() or !submat_block.size() )
      throw std::runtime_error("Must Populate Submat Maps");

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


  // Shell -> Task buffers
  shell_to_task_stack.shell_to_task_idx_device = 
    mem.aligned_alloc<int32_t>( total_nshells_task_batch );
  shell_to_task_stack.shell_to_task_off_device = 
    mem.aligned_alloc<int32_t>( total_nshells_task_batch );
  shell_to_task_stack.shell_to_task_device =
    mem.aligned_alloc<ShellToTaskDevice>( global_dims.nshells );


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

  std::vector< std::vector<int32_t> > shell_to_task_idx( global_dims.nshells ),
                                      shell_to_task_off( global_dims.nshells );

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
    const auto& submat_cut = it->submat_map;
    const auto& submat_block = it->submat_block;
    if( !submat_cut.size() or !submat_block.size() )
      throw std::runtime_error("Must Populate Submat Maps");
    const auto dist_nearest = it->dist_nearest;

    // Dimensions
    const size_t ncut     = submat_cut.size();
    const size_t nblock   = submat_block.size();
    const size_t nshells  = shell_list.size();
    const size_t npts     = points.size();
    const auto nbe        = it->nbe;

    // Compute offsets
    std::vector< size_t > shell_offs( nshells );
    shell_offs.at(0) = 0;
    for( auto i = 1ul; i < nshells; ++i )
      shell_offs.at(i) = shell_offs.at(i-1) + 
                           basis_map.shell_size( shell_list.at(i-1) );

    // Setup Shell -> Task
    const auto itask = std::distance( task_begin, it );
    for( auto i = 0; i < nshells; ++i ) {
      shell_to_task_idx[ shell_list.at(i) ].emplace_back(itask);
      shell_to_task_off[ shell_list.at(i) ].emplace_back(shell_offs[i]);
    }

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

    host_device_tasks.back().ibf_begin = 
      basis_map.shell_to_first_ao(shell_list[0]);
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

  // Construct Shell -> Task
  std::vector<int32_t> concat_shell_to_task_idx, concat_shell_to_task_off;
  {
    const size_t total_nshells = total_nshells_task_batch * sizeof(int32_t);
    buffer_adaptor 
      shell_idx_mem( shell_to_task_stack.shell_to_task_idx_device, total_nshells );
    buffer_adaptor 
      shell_off_mem( shell_to_task_stack.shell_to_task_off_device, total_nshells );


    std::vector<ShellToTaskDevice> host_shell_to_task(global_dims.nshells);
    for( auto ish = 0; ish < global_dims.nshells; ++ish ) {
      const auto ntask = shell_to_task_idx[ish].size();
      auto& bck = host_shell_to_task[ish];
      bck.ntask = ntask;
      bck.shell_device = static_stack.shells_device + ish;

      bck.task_idx_device        = shell_idx_mem.aligned_alloc<int32_t>( ntask );
      bck.task_shell_offs_device = shell_off_mem.aligned_alloc<int32_t>( ntask );

      // Pack data
      concat_iterable( concat_shell_to_task_idx, shell_to_task_idx[ish] );
      concat_iterable( concat_shell_to_task_off, shell_to_task_off[ish] );
    }

    // Send Data
    device_backend_->copy_async( concat_shell_to_task_idx.size(),
      concat_shell_to_task_idx.data(), 
      shell_to_task_stack.shell_to_task_idx_device, "shell_to_task_idx_device" );
    device_backend_->copy_async( concat_shell_to_task_off.size(),
      concat_shell_to_task_off.data(), 
      shell_to_task_stack.shell_to_task_off_device, "shell_to_task_off_device" );


    // Sort shell indices by L
    std::vector<uint32_t> shell_idx( global_dims.nshells );
    std::iota( shell_idx.begin(), shell_idx.end(), 0 );
    std::stable_sort( shell_idx.begin(), shell_idx.end(),
      [&]( auto i, auto j ){ 
        return basis_map.shell_l(i) < basis_map.shell_l(j); 
      } );

    {
    std::vector<ShellToTaskDevice> shell_to_task_sorted( global_dims.nshells );
    for( auto i = 0; i < global_dims.nshells; ++i ) 
      shell_to_task_sorted[i] = host_shell_to_task[shell_idx[i]];
    host_shell_to_task = std::move(shell_to_task_sorted);
    }

    // Send shell to task maps
    device_backend_->copy_async( global_dims.nshells, 
      host_shell_to_task.data(), shell_to_task_stack.shell_to_task_device,
      "shell_to_task_device" );


    // Form angular momenta batches
    auto max_l = basis_map.max_l();
    shell_to_task_stack.l_batched_shell_to_task.resize(max_l + 1);
    {
    auto* p = shell_to_task_stack.shell_to_task_device;
    for( auto l = 0; l <= max_l; ++l ) {
      auto nsh  = basis_map.nshells_with_l(l);
      auto pure = basis_map.l_purity(l);
      shell_to_task_stack.l_batched_shell_to_task[l].nshells_in_batch     = nsh;
      shell_to_task_stack.l_batched_shell_to_task[l].pure                 = pure;
      shell_to_task_stack.l_batched_shell_to_task[l].shell_to_task_device = p;
      p += nsh;
    }
    }
  } // Shell -> Task data

  // Synchronize on the copy stream to keep host vecs in scope
  device_backend_->master_queue_synchronize(); 


}



void XCDeviceAoSData::populate_submat_maps( 
  size_t N,
  host_task_iterator task_begin, host_task_iterator task_end, 
  const BasisSetMap& basis_map ) {


  // Get packing size 
  const size_t submat_chunk_size = this->get_submat_chunk_size(N,0);

  for( auto it = task_begin; it != task_end; ++it ) {

    const auto& shell_list = it->shell_list;
    std::tie( it->submat_map, it->submat_block ) = 
      gen_compressed_submat_map( basis_map, shell_list, N, submat_chunk_size );

  }

}

}
