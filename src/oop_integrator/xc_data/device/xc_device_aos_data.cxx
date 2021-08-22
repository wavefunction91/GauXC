#include "xc_device_aos_data.hpp"
#include "buffer_adaptor.hpp"
#include "integrator_util/integrator_common.hpp"

namespace GauXC {

size_t XCDeviceAoSData::get_mem_req( const host_task_type& task,
  const BasisSetMap& basis_map ) {

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
  const size_t mem_bf   = nbe * npts;
  const size_t mem_dbfx = mem_bf;
  const size_t mem_dbfy = mem_bf;
  const size_t mem_dbfz = mem_bf;

  // LDA/GGA Z Matrix 
  // TODO: this is dependent on integrand
  const size_t mem_zmat_lda_gga = nbe * npts;

  // nbe * nbe scratch
  const size_t mem_nbe_scr = nbe * nbe;

  // Shell index packing 
  const size_t mem_shell_list = nshells;
  const size_t mem_shell_offs = nshells;
  const size_t mem_submat_cut = 3 * ncut;
  const size_t mem_submat_block = nblock;

  // Memroty associated with added a task to the indirection
  const size_t mem_task = 1;


  size_t base_size = XCDeviceStackData::get_mem_req(task, basis_map);
  return base_size + 
    ( mem_bf + mem_dbfx + mem_dbfy + mem_dbfz + mem_nbe_scr ) * sizeof(double) +
    ( mem_zmat_lda_gga )                                      * sizeof(double) +
    ( mem_shell_list + mem_shell_offs )                       * sizeof(size_t) +
    ( mem_submat_cut + mem_submat_block )                     * sizeof(int32_t) +
    ( mem_task ) * sizeof(XCDeviceTask);
}



XCDeviceAoSData::device_buffer_t XCDeviceAoSData::alloc_pack_and_send( 
  host_task_iterator task_begin, host_task_iterator task_end, device_buffer_t buf,
  const BasisSetMap& basis_map ) {

  if( not device_backend_ ) throw std::runtime_error("Invalid Device Backend");

  // Allocate base info in the stack
  buf = XCDeviceStackData::alloc_pack_and_send( task_begin, task_end, buf, 
    basis_map );

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
  host_device_tasks.clear();

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

  // Allocate device memory
  auto [ ptr, sz ] = buf;
  buffer_adaptor mem( ptr, sz );

  // TODO: Print this if verbose
  //std::cout << "XCDeviceAoSData buf = " << ptr << ", " << sz << std::endl;


  // Device task indirection
  const size_t ntask = std::distance( task_begin, task_end );
  device_tasks = mem.aligned_alloc<XCDeviceTask>( ntask );

  // Collocation 
  bf_eval_device     = mem.aligned_alloc<double>( total_nbe_npts_task_batch );
  dbf_x_eval_device  = mem.aligned_alloc<double>( total_nbe_npts_task_batch );
  dbf_y_eval_device  = mem.aligned_alloc<double>( total_nbe_npts_task_batch );
  dbf_z_eval_device  = mem.aligned_alloc<double>( total_nbe_npts_task_batch );

  // VXC Z Matrix
  zmat_vxc_lda_gga_device = mem.aligned_alloc<double>( total_nbe_npts_task_batch );

  // Scratch buffer
  nbe_scr_device = mem.aligned_alloc<double>( total_nbe_sq_task_batch );

  // AoS buffers
  shell_list_device   = mem.aligned_alloc<size_t>( total_nshells_task_batch );
  shell_offs_device   = mem.aligned_alloc<size_t>( total_nshells_task_batch );
  submat_cut_device   = mem.aligned_alloc<int32_t>( 3 * total_ncut_task_batch );
  submat_block_device = mem.aligned_alloc<int32_t>( total_nblock_task_batch );

  // Send AoS information early to overlap with indirection construction
  device_backend_->copy_async( shell_list_pack.size(), shell_list_pack.data(), 
    shell_list_device, "send_shell_list" );
  device_backend_->copy_async( shell_offs_pack.size(), shell_offs_pack.data(), 
    shell_offs_device, "send_shell_offs" );
  device_backend_->copy_async( 3 * submat_cut_pack.size(), 
    submat_cut_pack.data()->data(), submat_cut_device, "send_submat_cut"  ); 
  device_backend_->copy_async( submat_block_pack.size(), submat_block_pack.data(), 
    submat_block_device, "send_submat_block"  ); 

  // Construct full indirection
  {
  double* points_ptr  = this->points_device;
  double* weights_ptr = this->weights_device;

  size_t* shell_list_ptr  = shell_list_device;
  size_t* shell_offs_ptr  = shell_offs_device;
  int32_t* submat_cut_ptr = submat_cut_device;
  int32_t* submat_block_ptr = submat_block_device;

  double*      nbe_ptr    = nbe_scr_device;
  double*      zmat_ptr   = zmat_vxc_lda_gga_device;

  double*      bf_ptr     = bf_eval_device;
  double*      dbfx_ptr   = dbf_x_eval_device;
  double*      dbfy_ptr   = dbf_y_eval_device;
  double*      dbfz_ptr   = dbf_z_eval_device;
  
  double*      den_ptr    = den_eval_device;
  double*      ddenx_ptr  = den_x_eval_device;
  double*      ddeny_ptr  = den_y_eval_device;
  double*      ddenz_ptr  = den_z_eval_device;

  double*      eps_ptr     = eps_eval_device;
  double*      gamma_ptr   = gamma_eval_device;
  double*      vrho_ptr    = vrho_eval_device;
  double*      vgamma_ptr  = vgamma_eval_device;

  for( auto& task : host_device_tasks ) {

    task.points     = points_ptr;
    task.weights    = weights_ptr;
    task.shell_list = shell_list_ptr;
    task.shell_offs = shell_offs_ptr;
    task.submat_cut = submat_cut_ptr;
    task.submat_block = submat_block_ptr;
    
    task.nbe_scr = nbe_ptr;
    task.zmat    = zmat_ptr;
    task.bf      = bf_ptr;
    task.dbfx    = dbfx_ptr;
    task.dbfy    = dbfy_ptr;
    task.dbfz    = dbfz_ptr;
    task.den     = den_ptr;
    task.ddenx   = ddenx_ptr;
    task.ddeny   = ddeny_ptr;
    task.ddenz   = ddenz_ptr;

    task.eps    = eps_ptr;
    task.gamma  = gamma_ptr;
    task.vrho   = vrho_ptr;
    task.vgamma = vgamma_ptr;

    //task.dist_scratch      = dist_scratch_ptr;

    auto npts    = task.npts;
    auto nbe     = task.nbe;
    auto ncut    = task.ncut;
    auto nblock  = task.nblock;
    auto nshells = task.nshells;

    points_ptr     += 3 * npts;
    weights_ptr    += npts;
    shell_list_ptr += nshells;
    shell_offs_ptr += nshells;
    submat_cut_ptr += 3 * ncut;
    submat_block_ptr += nblock;
    
    nbe_ptr    += nbe * nbe;
    zmat_ptr   += nbe * npts;

    bf_ptr     += nbe * npts;
    dbfx_ptr   += nbe * npts;
    dbfy_ptr   += nbe * npts;
    dbfz_ptr   += nbe * npts;

    den_ptr    += npts;
    ddenx_ptr  += npts;
    ddeny_ptr  += npts;
    ddenz_ptr  += npts;

    eps_ptr    += npts;
    gamma_ptr  += npts;
    vrho_ptr   += npts;
    vgamma_ptr += npts;
    //dist_scratch_ptr += LDatoms * npts;

  } // Loop over device tasks

  } // End task setup

  // Setup extra pieces to indirection which are algorithm specific
  device_buffer_t buf_left{ mem.stack(), mem.nleft() };
  buf_left = add_extra_to_indirection(host_device_tasks, buf_left);


  // Send indirection
  device_backend_->copy_async( host_device_tasks.size(), host_device_tasks.data(), 
    device_tasks, "send_tasks_device" );


  // Synchronize on the copy stream to keep host vecs in scope
  device_backend_->master_queue_synchronize(); 

  // Update dynmem data for derived impls
  return buf_left;
}

}
