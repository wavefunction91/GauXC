/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#include "scheme1_data_base.hpp"
#include "buffer_adaptor.hpp"

namespace GauXC {

Scheme1DataBase::~Scheme1DataBase() noexcept = default;

Scheme1DataBase::Scheme1DataBase(const DeviceRuntimeEnvironment& rt) : 
  base_type(rt) {

  if( device_backend_ ) 
    device_backend_->create_blas_queue_pool(4);

}

void Scheme1DataBase::reset_allocations() {
  base_type::reset_allocations();
  scheme1_stack.reset();
  collocation_stack.reset();
  //coulomb_stack.reset();
  shell_to_task_stack.reset();
  shell_pair_to_task_stack.reset();
  
  l_batched_shell_to_task.clear();

  task_to_shell_pair_stack.reset();
  subtask.clear();
  nprim_pairs_host.clear();
  sp_ptr_host.clear();

  sp_X_AB_host.clear();
  sp_Y_AB_host.clear();
  sp_Z_AB_host.clear();

  task_to_shell_pair.clear();
}

size_t Scheme1DataBase::get_static_mem_requirement() {
  size_t size = 0;

  const size_t nsp = global_dims.nshell_pairs;
  size += 
    // Shell Pair map
    global_dims.nshells * sizeof(ShellToTaskDevice) +
    nsp * sizeof(ShellPairToTaskDevice) +
    // Task Map
    nsp * sizeof(int32_t) +      // nprim_pairs
    nsp * sizeof(shell_pair*) +  // shell_pair pointer
    nsp * 3 * sizeof(double);    // X_AB, Y_AB, Z_AB

  return size;
}









size_t Scheme1DataBase::get_mem_req( integrator_term_tracker terms, 
  const host_task_type& task ) {

  // All local memory is weights related
  size_t base_size = base_type::get_mem_req(terms, task);

#if 0
  if( terms.weights ) {
    const auto ldatoms = get_ldatoms();
    const auto mem_dist_scr = ldatoms * task.npts;
    const auto mem_dist_ner = task.npts;
    const auto mem_iparent  = task.npts;

    base_size += (mem_dist_scr + mem_dist_ner) * sizeof(double) + 
                 mem_iparent * sizeof(int32_t);
  }

  if( terms.exx or terms.exc_vxc or terms.exc_grad or terms.den ) {
    const auto& shell_list_bfn = task.bfn_screening.shell_list;
    const size_t nshells_bfn  = shell_list_bfn.size();
    const auto& shell_list_cou = task.cou_screening.shell_list;
    const size_t nshells_cou  = shell_list_cou.size();

    const size_t mem_shell_list_bfn = nshells_bfn * sizeof(size_t);
    const size_t mem_shell_offs_bfn = nshells_bfn * sizeof(size_t);

    const size_t mem_shell_list_cou = /* terms.exx ? nshells_cou * sizeof(size_t) : */ 0;
    const size_t mem_shell_offs_cou = /* terms.exx ? nshells_cou * sizeof(size_t) : */ 0;

    // Shell -> Task maps
    const size_t mem_shell_to_task_idx  = nshells_bfn * sizeof(int32_t);
    const size_t mem_shell_to_task_off  = nshells_bfn * sizeof(int32_t);

    // ShellPair -> Task maps
    const size_t nshells_cou_sqlt = (nshells_cou*(nshells_cou+1))/2;
    const size_t mem_shell_pair_to_task_idx     = terms.exx ? nshells_cou_sqlt * sizeof(int32_t) : 0;
    const size_t mem_shell_pair_to_task_row_off = terms.exx ? nshells_cou_sqlt * sizeof(int32_t) : 0;
    const size_t mem_shell_pair_to_task_col_off = terms.exx ? nshells_cou_sqlt * sizeof(int32_t) : 0;

    base_size +=  mem_shell_list_bfn + mem_shell_offs_bfn +
      mem_shell_list_cou + mem_shell_offs_cou +
      mem_shell_to_task_idx + mem_shell_to_task_off +
      mem_shell_pair_to_task_idx + mem_shell_pair_to_task_row_off + mem_shell_pair_to_task_col_off;
  }
#else
  required_term_storage reqt(terms);
  const auto ldatoms = get_ldatoms();
  const auto npts = task.npts;
  const auto& shell_list_bfn = task.bfn_screening.shell_list;
  const auto& shell_list_cou = task.cou_screening.shell_list;
  const size_t nshells_bfn  = shell_list_bfn.size();
  const size_t nshells_cou  = shell_list_cou.size();

  const int max_l = global_dims.max_l;
  const size_t n_sp_types = (max_l+1) * (max_l+1); 
  const size_t n_sp_types_with_diag = n_sp_types + (max_l+1); 

  base_size += 
    // Weights specific memory
    reqt.grid_to_center_dist_scr_size(ldatoms, npts) * sizeof(double)  +
    reqt.grid_to_center_dist_nearest_size(npts)      * sizeof(double)  +
    reqt.grid_to_parent_center_size(npts)            * sizeof(int32_t) +

    // Shell / Shell Pair lists + indirection
    reqt.task_shell_list_bfn_size(nshells_bfn)            * sizeof(size_t)  +
    reqt.task_shell_offs_bfn_size(nshells_bfn)            * sizeof(size_t)  +
    reqt.shell_to_task_idx_bfn_size(nshells_bfn)          * sizeof(int32_t) +
    reqt.shell_to_task_off_bfn_size(nshells_bfn)          * sizeof(int32_t) +
    reqt.shell_pair_to_task_idx_cou_size(nshells_cou)     * sizeof(int32_t) +
    reqt.shell_pair_to_task_row_off_cou_size(nshells_cou) * sizeof(int32_t) +
    reqt.shell_pair_to_task_col_off_cou_size(nshells_cou) * sizeof(int32_t) +

    // Task to shell pair map
    reqt.task_to_shell_pair_cou_size() * n_sp_types_with_diag * sizeof(TaskToShellPairDevice) +
    reqt.task_to_shell_pair_col_off_cou_size(nshells_cou) * sizeof(int32_t) +
    reqt.task_to_shell_pair_row_off_cou_size(nshells_cou) * sizeof(int32_t) +
    reqt.task_to_shell_pair_idx_cou_size(nshells_cou) * sizeof(int32_t) + 
    reqt.task_to_shell_pair_cou_subtask_size(npts, 256) * sizeof(std::array<int32_t, 4>);

#endif

  //std::cout << "MEM REQ: " << base_size << std::endl;
  return base_size;
}










Scheme1DataBase::device_buffer_t Scheme1DataBase::allocate_dynamic_stack( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end, 
  device_buffer_t buf ){

  // Allocate base info on the stack
  buf = base_type::allocate_dynamic_stack( terms, task_begin, task_end,
    buf );

  // Allocate additional device memory 
  auto [ ptr, sz ] = buf;
  buffer_adaptor mem( ptr, sz );

#if 0
  // Weights related memory
  if( terms.weights ) { 
    const auto ldatoms = get_ldatoms();
    scheme1_stack.dist_scratch_device = 
      mem.aligned_alloc<double>( ldatoms * total_npts_task_batch, sizeof(double2), csl );
    scheme1_stack.dist_nearest_device = 
      mem.aligned_alloc<double>( total_npts_task_batch, csl );
    scheme1_stack.iparent_device = 
      mem.aligned_alloc<int32_t>( total_npts_task_batch, csl );
  }

  if( terms.exx or terms.exc_vxc or terms.exc_grad or terms.den ) {
    total_nshells_bfn_task_batch  = 0; 
    //total_nshells_cou_task_batch  = 0; 
    total_nshells_cou_sqlt_task_batch  = 0; 
    for( auto it = task_begin; it != task_end; ++it ) {
      const auto& shell_list_bfn  = it->bfn_screening.shell_list;
      const size_t nshells_bfn  = shell_list_bfn.size();
      total_nshells_bfn_task_batch  += nshells_bfn;

      const auto& shell_list_cou  = it->cou_screening.shell_list;
      const size_t nshells_cou  = shell_list_cou.size();
      const size_t nshells_cou_sqlt = (nshells_cou*(nshells_cou+1))/2;
      //total_nshells_cou_task_batch  += nshells_cou;
      total_nshells_cou_sqlt_task_batch  += nshells_cou_sqlt;
    }

    collocation_stack.shell_list_device = 
      mem.aligned_alloc<size_t>( total_nshells_bfn_task_batch , csl);
    collocation_stack.shell_offs_device = 
      mem.aligned_alloc<size_t>( total_nshells_bfn_task_batch , csl);

    // Shell -> Task buffers
    shell_to_task_stack.shell_to_task_idx_device = 
      mem.aligned_alloc<int32_t>( total_nshells_bfn_task_batch, csl );
    shell_to_task_stack.shell_to_task_off_device = 
      mem.aligned_alloc<int32_t>( total_nshells_bfn_task_batch, csl );
    shell_to_task_stack.shell_to_task_device =
      mem.aligned_alloc<ShellToTaskDevice>( global_dims.nshells, csl );

    if(terms.exx) {
    #if 0
      coulomb_stack.shell_list_device = 
        mem.aligned_alloc<size_t>( total_nshells_cou_task_batch , csl);
      coulomb_stack.shell_offs_device = 
        mem.aligned_alloc<size_t>( total_nshells_cou_task_batch , csl);
    #endif
      
      // ShellPair -> Task buffers
      const size_t nsp = (global_dims.nshells*(global_dims.nshells+1))/2;
      shell_pair_to_task_stack.shell_pair_to_task_idx_device = 
        mem.aligned_alloc<int32_t>( total_nshells_cou_sqlt_task_batch, csl );
      shell_pair_to_task_stack.shell_pair_to_task_row_off_device = 
        mem.aligned_alloc<int32_t>( total_nshells_cou_sqlt_task_batch, csl );
      shell_pair_to_task_stack.shell_pair_to_task_col_off_device = 
        mem.aligned_alloc<int32_t>( total_nshells_cou_sqlt_task_batch, csl );
      shell_pair_to_task_stack.shell_pair_to_task_device =
        mem.aligned_alloc<ShellPairToTaskDevice>( nsp, csl );
    }

  }
#else

  required_term_storage reqt(terms);

  // Weights related memory
  if(reqt.grid_to_center_dist_scr) {
    const auto ldatoms = get_ldatoms();
    scheme1_stack.dist_scratch_device = mem.aligned_alloc<double>( 
      ldatoms * total_npts_task_batch, alignof(double2), csl );
  }
  if(reqt.grid_to_center_dist_nearest) {
    scheme1_stack.dist_nearest_device = 
      mem.aligned_alloc<double>( total_npts_task_batch, csl );
  }
  if(reqt.grid_to_parent_center) {
    scheme1_stack.iparent_device = 
      mem.aligned_alloc<int32_t>( total_npts_task_batch, csl );
  }

  // Compute total dimensions for shell(pair) lists
  total_nshells_bfn_task_batch       = 0; 
  total_nshells_cou_sqlt_task_batch  = 0; 
  size_t num_subtasks = 0;
  const int points_per_subtask = get_points_per_subtask();
  for( auto it = task_begin; it != task_end; ++it ) {
    const auto& shell_list_bfn  = it->bfn_screening.shell_list;
    const size_t nshells_bfn  = shell_list_bfn.size();
    total_nshells_bfn_task_batch  += nshells_bfn;

    const auto& shell_list_cou  = it->cou_screening.shell_list;
    const size_t nshells_cou  = shell_list_cou.size();
    const size_t nshells_cou_sqlt = (nshells_cou*(nshells_cou+1))/2;
    total_nshells_cou_sqlt_task_batch  += nshells_cou_sqlt;

    num_subtasks += util::div_ceil(it->npts, points_per_subtask);
  }

  // Shell lists and offs (bfn)
  if(reqt.task_shell_list_bfn) {
    collocation_stack.shell_list_device = 
      mem.aligned_alloc<size_t>( total_nshells_bfn_task_batch , csl);
  }
  if(reqt.task_shell_offs_bfn) {
    collocation_stack.shell_offs_device = 
      mem.aligned_alloc<size_t>( total_nshells_bfn_task_batch , csl);
  }

  // Shell -> Task buffers
  if(reqt.shell_to_task_bfn) {
    shell_to_task_stack.shell_to_task_idx_device = 
      mem.aligned_alloc<int32_t>( total_nshells_bfn_task_batch, csl );
    shell_to_task_stack.shell_to_task_off_device = 
      mem.aligned_alloc<int32_t>( total_nshells_bfn_task_batch, csl );
    shell_to_task_stack.shell_to_task_device =
      mem.aligned_alloc<ShellToTaskDevice>( global_dims.nshells, csl );
  }

  const size_t nsp = global_dims.nshell_pairs;

  // ShellPair -> Task buffer (cou)
  if(reqt.shell_pair_to_task_cou) {
    throw std::runtime_error("SPARSE + SP2TASK NYI");
    shell_pair_to_task_stack.shell_pair_to_task_idx_device = 
      mem.aligned_alloc<int32_t>( total_nshells_cou_sqlt_task_batch, csl );
    shell_pair_to_task_stack.shell_pair_to_task_row_off_device = 
      mem.aligned_alloc<int32_t>( total_nshells_cou_sqlt_task_batch, csl );
    shell_pair_to_task_stack.shell_pair_to_task_col_off_device = 
      mem.aligned_alloc<int32_t>( total_nshells_cou_sqlt_task_batch, csl );

    shell_pair_to_task_stack.shell_pair_to_task_device =
      mem.aligned_alloc<ShellPairToTaskDevice>( nsp, csl );
  }

  // Task -> ShellPair (cou)
  if(reqt.task_to_shell_pair_cou) { 
    const size_t ntasks = std::distance(task_begin, task_end);
    const int max_l = global_dims.max_l;
    const size_t n_sp_types = (max_l+1) * (max_l+1); 
    const size_t n_sp_types_with_diag = n_sp_types + (max_l+1); 
    task_to_shell_pair_stack.task_to_shell_pair_device = 
      mem.aligned_alloc<TaskToShellPairDevice>( ntasks * n_sp_types_with_diag, csl );

    task_to_shell_pair_stack.task_shell_linear_idx_device = 
      mem.aligned_alloc<int32_t>( total_nshells_cou_sqlt_task_batch, csl);
    task_to_shell_pair_stack.task_shell_off_row_device = 
      mem.aligned_alloc<int32_t>( total_nshells_cou_sqlt_task_batch, csl);
    task_to_shell_pair_stack.task_shell_off_col_device = 
      mem.aligned_alloc<int32_t>( total_nshells_cou_sqlt_task_batch, csl);

    task_to_shell_pair_stack.subtask_device = 
      mem.aligned_alloc<std::array<int32_t, 4>>( num_subtasks, 16, csl );

    task_to_shell_pair_stack.nprim_pairs_device = 
      mem.aligned_alloc<int32_t>( nsp, 16, csl );
    task_to_shell_pair_stack.sp_ptr_device = 
      mem.aligned_alloc<shell_pair*>( nsp, 16, csl );
    task_to_shell_pair_stack.sp_X_AB_device = 
      mem.aligned_alloc<double>( nsp, 16, csl );
    task_to_shell_pair_stack.sp_Y_AB_device = 
      mem.aligned_alloc<double>( nsp, 16, csl );
    task_to_shell_pair_stack.sp_Z_AB_device = 
      mem.aligned_alloc<double>( nsp, 16, csl );

  }


#endif

  // Update dynmem data for derived impls
  return device_buffer_t{ mem.stack(), mem.nleft() };
}









void Scheme1DataBase::pack_and_send( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end,
  const BasisSetMap& basis_map ) {

  // Pack and send base data
  base_type::pack_and_send( terms, task_begin, task_end, basis_map );

#if 0
  // All local memory is weights related
  if( terms.weights ) { 
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
    device_backend_->copy_async( dist_nearest_pack.size(), 
      dist_nearest_pack.data(), scheme1_stack.dist_nearest_device, 
      "send dist_nearest" );
  }

  if( terms.exx or terms.exc_vxc or terms.exc_grad or terms.den ) {
    // Contatenation utility
    auto concat_iterable = []( auto& a, const auto& b ) {
      a.insert( a.end(), b.begin(), b.end() );
    };

    std::vector< size_t > shell_list_bfn_pack;
    std::vector< size_t > shell_offs_bfn_pack;

    std::vector< std::vector<int32_t> > shell_to_task_idx( global_dims.nshells ),
                                        shell_to_task_off( global_dims.nshells );


    // Shell Pair -> Task
    if( terms.exx ) {
      shell_pair_to_task.clear();
      const size_t nsp = (global_dims.nshells*(global_dims.nshells+1))/2;
      shell_pair_to_task.resize(nsp);
      for( auto i = 0ul; i < nsp; ++i ) {
        auto& sptt = shell_pair_to_task[i];
        sptt.shell_pair_device =
          this->shell_pair_soa.shell_pair_dev_ptr[i];
        std::tie(sptt.lA, sptt.lB) =
          this->shell_pair_soa.shell_pair_ls[i];
        std::tie(sptt.rA, sptt.rB) =
          this->shell_pair_soa.shell_pair_centers[i];
        std::tie(sptt.idx_bra, sptt.idx_ket) =
          this->shell_pair_soa.shell_pair_shidx[i];
      }
    }

    for( auto it = task_begin; it != task_end; ++it ) {
      const auto& shell_list_bfn  = it->bfn_screening.shell_list;
      const size_t nshells_bfn  = shell_list_bfn.size();

      // Compute offsets (bfn)
#if 0 
      std::vector< size_t > shell_offs_bfn( nshells_bfn );
      shell_offs_bfn.at(0) = 0;
      for( auto i = 1ul; i < nshells_bfn; ++i )
        shell_offs_bfn.at(i) = shell_offs_bfn.at(i-1) + 
                             basis_map.shell_size( shell_list_bfn.at(i-1) );
#else
      auto shell_offs_bfn = basis_map.shell_offs<size_t>( 
        shell_list_bfn.begin(), shell_list_bfn.end() );
#endif

      concat_iterable( shell_list_bfn_pack, shell_list_bfn );
      concat_iterable( shell_offs_bfn_pack, shell_offs_bfn );



      // Setup Shell -> Task
      const auto itask = std::distance( task_begin, it );
      for( auto i = 0ul; i < nshells_bfn; ++i ) {
        shell_to_task_idx[ shell_list_bfn.at(i) ].emplace_back(itask);
        shell_to_task_off[ shell_list_bfn.at(i) ].emplace_back(shell_offs_bfn[i]);
      }


      if(terms.exx) {
        const auto& shell_list_cou  = it->cou_screening.shell_list;
        const size_t nshells_cou  = shell_list_cou.size();

        // Compute offsets (cou)
#if 0
        std::vector< size_t > shell_offs_cou( nshells_cou );
        shell_offs_cou.at(0) = 0;
        for( auto i = 1ul; i < nshells_cou; ++i )
          shell_offs_cou.at(i) = shell_offs_cou.at(i-1) + 
                               basis_map.shell_size( shell_list_cou.at(i-1) );
#else
        auto shell_offs_cou = basis_map.shell_offs<size_t>( 
          shell_list_cou.begin(), shell_list_cou.end() );
#endif

        // Setup Shell Pair -> Task
        for( auto j = 0ul; j < nshells_cou; ++j )
        for( auto i = j;   i < nshells_cou; ++i ) {
          const auto ish = shell_list_cou[i];
          const auto jsh = shell_list_cou[j];
          const auto idx = detail::packed_lt_index(ish,jsh, global_dims.nshells);

          auto& sptt = shell_pair_to_task[idx];
          sptt.task_idx.emplace_back(itask);
          sptt.task_shell_off_row.emplace_back(shell_offs_cou[i]);
          sptt.task_shell_off_col.emplace_back(shell_offs_cou[j]);
        }
      }


    }

    // Send Shell List data (bfn)
    device_backend_->copy_async( shell_list_bfn_pack.size(), 
      shell_list_bfn_pack.data(), collocation_stack.shell_list_device, 
      "send_shell_list_bfn" );
    device_backend_->copy_async( shell_offs_bfn_pack.size(), 
      shell_offs_bfn_pack.data(), collocation_stack.shell_offs_device, 
      "send_shell_offs_bfn" );

    // Construct Shell -> Task (bfn)
    std::vector<int32_t> concat_shell_to_task_idx, concat_shell_to_task_off;
    {
      const size_t total_nshells_bfn = total_nshells_bfn_task_batch * sizeof(int32_t);
      buffer_adaptor 
        shell_idx_mem( shell_to_task_stack.shell_to_task_idx_device, total_nshells_bfn );
      buffer_adaptor                                                                  
        shell_off_mem( shell_to_task_stack.shell_to_task_off_device, total_nshells_bfn );


      std::vector<ShellToTaskDevice> host_shell_to_task(global_dims.nshells);
      for( auto ish = 0ul; ish < global_dims.nshells; ++ish ) {
        const auto ntask = shell_to_task_idx[ish].size();
        auto& bck = host_shell_to_task[ish];
        bck.ntask = ntask;
        bck.center_idx = basis_map.shell_to_center( ish );
        bck.true_idx   = ish;
        bck.shell_device = static_stack.shells_device + ish;

        bck.task_idx_device        = shell_idx_mem.aligned_alloc<int32_t>( ntask, csl );
        bck.task_shell_offs_device = shell_off_mem.aligned_alloc<int32_t>( ntask, csl );

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
      for( auto i = 0ul; i < global_dims.nshells; ++i ) 
        shell_to_task_sorted[i] = host_shell_to_task[shell_idx[i]];
      host_shell_to_task = std::move(shell_to_task_sorted);
      }

      // Send shell to task maps
      device_backend_->copy_async( global_dims.nshells, 
        host_shell_to_task.data(), shell_to_task_stack.shell_to_task_device,
        "shell_to_task_device" );


      // Form angular momenta batches
      auto max_l = basis_map.max_l();
      l_batched_shell_to_task.resize(max_l + 1);
      {
      auto* p = shell_to_task_stack.shell_to_task_device;
      auto* h = host_shell_to_task.data();
      for( auto l = 0ul; l <= max_l; ++l ) {
        auto nsh  = basis_map.nshells_with_l(l);
        auto pure = basis_map.l_purity(l);
        l_batched_shell_to_task[l].nshells_in_batch     = nsh;
        l_batched_shell_to_task[l].pure                 = pure;
        l_batched_shell_to_task[l].shell_to_task_device = p;
                            
        size_t max_ntask = std::max_element( h, h+nsh,
          [](auto& a, auto& b){ return a.ntask < b.ntask; } )->ntask;

        l_batched_shell_to_task[l].ntask_average = max_ntask;
        l_batched_shell_to_task[l].npts_average  = 0;

        p += nsh;
        h += nsh;
      }
      }
    } // Shell -> Task data





    // Construct ShellPair -> Task
    std::vector<int32_t> concat_shell_pair_to_task_idx, concat_shell_pair_to_task_off_row, concat_shell_pair_to_task_off_col;
    if(terms.exx) {
      const size_t total_nshells_cou_sqlt = total_nshells_cou_sqlt_task_batch * sizeof(int32_t);
      buffer_adaptor 
        shell_idx_mem( shell_pair_to_task_stack.shell_pair_to_task_idx_device, total_nshells_cou_sqlt );
      buffer_adaptor                                                                            
        shell_row_off_mem( shell_pair_to_task_stack.shell_pair_to_task_row_off_device, total_nshells_cou_sqlt );
      buffer_adaptor                                                                            
        shell_col_off_mem( shell_pair_to_task_stack.shell_pair_to_task_col_off_device, total_nshells_cou_sqlt );

      // Sort Shell Pairs by diag + L pairs
      const size_t nsp = (global_dims.nshells*(global_dims.nshells+1))/2;
      {
        std::vector<int> sp_idx(nsp);

        // First partition the shell pairs into diagonal and off diagonal
        std::iota( sp_idx.begin(), sp_idx.end(), 0 );
        auto diag_end = std::stable_partition( sp_idx.begin(), sp_idx.end(),
          [&](auto idx) {
            auto [i,j] = detail::from_packed_lt_index(idx, global_dims.nshells);
            return i == j;
          });

        // Next sort diag and off diag by shell pair
        auto lp_functor = [&](int i, int j){
          const auto& sp_i = shell_pair_to_task[i];
          const auto& sp_j = shell_pair_to_task[j];

          auto l_i = std::make_pair( sp_i.lA, sp_i.lB );
          auto l_j = std::make_pair( sp_j.lA, sp_j.lB );
          return l_i < l_j;
        };

        std::stable_sort( sp_idx.begin(), diag_end, lp_functor );
        std::stable_sort( diag_end, sp_idx.end(),   lp_functor );

        // Reorder shell pairs
        std::vector<ShellPairToTaskHost> sorted_shell_pair_to_task(nsp);
        for( auto i = 0ul; i < nsp; ++i ) {
          sorted_shell_pair_to_task[i] = shell_pair_to_task[sp_idx[i]];
        }
        shell_pair_to_task = std::move(sorted_shell_pair_to_task);
      }


      std::vector<ShellPairToTaskDevice> host_shell_pair_to_task(nsp);
      for( auto isp = 0ul; isp < nsp; ++isp ) {
        auto& sptt = shell_pair_to_task[isp];
        auto& bck = host_shell_pair_to_task[isp];

        const auto ntask = sptt.task_idx.size();
        bck.ntask = ntask;
        bck.X_AB = sptt.rA.x - sptt.rB.x;
        bck.Y_AB = sptt.rA.y - sptt.rB.y;
        bck.Z_AB = sptt.rA.z - sptt.rB.z;
        bck.shell_pair_device = sptt.shell_pair_device;

        bck.task_idx_device           = shell_idx_mem.aligned_alloc<int32_t>( ntask, csl );
        bck.task_shell_off_row_device = shell_row_off_mem.aligned_alloc<int32_t>( ntask, csl );
        bck.task_shell_off_col_device = shell_col_off_mem.aligned_alloc<int32_t>( ntask, csl );

        // Pack data
        concat_iterable( concat_shell_pair_to_task_idx,     sptt.task_idx           );
        concat_iterable( concat_shell_pair_to_task_off_row, sptt.task_shell_off_row );
        concat_iterable( concat_shell_pair_to_task_off_col, sptt.task_shell_off_col );
      }

      // Send Data
      device_backend_->copy_async( concat_shell_pair_to_task_idx.size(),
        concat_shell_pair_to_task_idx.data(), 
        shell_pair_to_task_stack.shell_pair_to_task_idx_device, "shell_pair_to_task_idx_device" );

      device_backend_->copy_async( concat_shell_pair_to_task_off_row.size(),
        concat_shell_pair_to_task_off_row.data(), 
        shell_pair_to_task_stack.shell_pair_to_task_row_off_device, "shell_pair_to_task_row_off_device" );

      device_backend_->copy_async( concat_shell_pair_to_task_off_col.size(),
        concat_shell_pair_to_task_off_col.data(), 
        shell_pair_to_task_stack.shell_pair_to_task_col_off_device, "shell_pair_to_task_col_off_device" );

      // Send shell to task maps
      device_backend_->copy_async( nsp, 
        host_shell_pair_to_task.data(), shell_pair_to_task_stack.shell_pair_to_task_device,
        "shell_pair_to_task_device" );

      
      const int max_l = basis_map.max_l();

      // Setup DIAGONAL SP AM batches
      l_batched_shell_pair_to_task_diag.clear();
      l_batched_shell_pair_to_task_diag.resize(max_l+1);
      {
      auto* p = shell_pair_to_task_stack.shell_pair_to_task_device;
      auto h  = shell_pair_to_task.begin();
      for( auto l = 0; l <= max_l; ++l ) {
        auto& batch_map = l_batched_shell_pair_to_task_diag[l];
        const auto nsh = basis_map.nshells_with_l(l);

        size_t max_npts = 0;
        size_t max_ntask = std::max_element( h, h + nsh,
          [](auto& a, auto& b){
            return a.task_idx.size() < b.task_idx.size();
          })->task_idx.size();

        batch_map.ntask_average = max_ntask;
        batch_map.npts_average  = max_npts;
        batch_map.nshells_in_batch = nsh;
        batch_map.shell_pair_to_task_device = p;
        batch_map.lA = l;
        batch_map.lB = l;

        p += nsh;
        h += nsh;
      }
      }

      // Setup OFFDIAGONAL SP AM batches
      l_batched_shell_pair_to_task_off_diag.clear();
      l_batched_shell_pair_to_task_off_diag.resize((max_l+1)*(max_l+1));
      {
      auto* p = shell_pair_to_task_stack.shell_pair_to_task_device +
        global_dims.nshells;
      auto h  = shell_pair_to_task.begin() + global_dims.nshells;
      for( auto l_i = 0, l_ij = 0; l_i <= max_l; ++l_i )
      for( auto l_j = 0; l_j <= max_l; ++l_j, ++l_ij ) {
        auto& batch_map = l_batched_shell_pair_to_task_off_diag[l_ij];
        size_t nsh = std::count_if(
          h, shell_pair_to_task.end(),
          [=](const auto& sp) {
            return sp.lA == l_i and sp.lB == l_j;
          });

        size_t max_npts = 0;
        size_t max_ntask = std::max_element( h, h + nsh,
          [](auto& a, auto& b){
            return a.task_idx.size() < b.task_idx.size();
          })->task_idx.size();

        batch_map.ntask_average = max_ntask;
        batch_map.npts_average  = max_npts;
        batch_map.nshells_in_batch = nsh;
        batch_map.shell_pair_to_task_device = p;
        batch_map.lA = l_i;
        batch_map.lB = l_j;

        p += nsh;
        h += nsh;
      }
      }


      std::cout << "DIAG SP Batches" << std::endl;
      for( auto& b : l_batched_shell_pair_to_task_diag ) {
        std::cout << b.lA << " " << b.lB << " " << b.nshells_in_batch << " " 
          << b.ntask_average << std::endl;
      }

      std::cout << "OFFDIAG SP Batches" << std::endl;
      for( auto& b : l_batched_shell_pair_to_task_off_diag ) {
        std::cout << b.lA << " " << b.lB << " " << b.nshells_in_batch << " " 
          << b.ntask_average << std::endl;
      }
    } // ShellPair -> Task data
    device_backend_->master_queue_synchronize(); 
  }
#else

  required_term_storage reqt(terms);

  // Host Packing Arrays
  std::vector<int32_t> iparent_pack;
  std::vector<double>  dist_nearest_pack;
  std::vector<size_t> shell_list_bfn_pack, shell_offs_bfn_pack;
  std::vector< std::vector<int32_t> > 
    shell_to_task_idx_bfn, shell_to_task_off_bfn;
  std::vector<int32_t>
    concat_shell_to_task_idx_bfn, concat_shell_to_task_off_bfn,
    concat_shell_pair_to_task_idx_cou, concat_shell_pair_to_task_off_row_cou, 
    concat_shell_pair_to_task_off_col_cou;
  std::vector<ShellToTaskDevice> host_shell_to_task_bfn;
  std::vector<ShellPairToTaskDevice> host_shell_pair_to_task_cou;

  // Contatenation utility
  auto concat_iterable = []( auto& a, const auto& b ) {
    a.insert( a.end(), b.begin(), b.end() );
  };

  using hrt_t = std::chrono::high_resolution_clock;
  using dur_t = std::chrono::duration<double,std::milli>;

  /*******************************************
   *         WEIGHTS RELATED MEMORY          *
   *******************************************/

  auto w_mem_st = hrt_t::now();
  // Nearest Distance Array
  if(reqt.grid_to_center_dist_nearest) {

    // Pack on host
    dist_nearest_pack.reserve( total_npts_task_batch );
    for( auto it = task_begin; it != task_end; ++it ) {
      dist_nearest_pack.insert( dist_nearest_pack.end(), it->points.size(), 
        it->dist_nearest );
    }

    // Send to device
    device_backend_->copy_async( dist_nearest_pack.size(), 
      dist_nearest_pack.data(), scheme1_stack.dist_nearest_device, 
      "send dist_nearest" );

  }

  // IParent Array
  if(reqt.grid_to_parent_center) {

    // Pack on host
    iparent_pack.reserve( total_npts_task_batch );
    for( auto it = task_begin; it != task_end; ++it ) {
      iparent_pack.insert( iparent_pack.end(), it->points.size(), it->iParent );
    }

    // Send to device
    device_backend_->copy_async( iparent_pack.size(), iparent_pack.data(), 
      scheme1_stack.iparent_device, "send iparent"  );

  }
  auto w_mem_en = hrt_t::now();

  /************************************************
   * SHELL LIST, OFFSET and TASK MAP MEMORY (bfn) *
   ************************************************/

  auto sl_mem_st = hrt_t::now();
  // Resize host arrays for Shell -> Task
  if(reqt.shell_to_task_bfn) {
    shell_to_task_idx_bfn.resize( global_dims.nshells );
    shell_to_task_off_bfn.resize( global_dims.nshells );
  }

  // Shell list, offsets + task map (bfn)
  for( auto it = task_begin; it != task_end; ++it ) {
    const auto& shell_list_bfn  = it->bfn_screening.shell_list;
    const size_t nshells_bfn  = shell_list_bfn.size();

    // Pack shell list (bfn)
    if(reqt.task_shell_list_bfn) {
      concat_iterable( shell_list_bfn_pack, shell_list_bfn );
    }
    
    // Generate and pack shell offsets (bfn)
    std::vector<size_t> shell_offs_bfn;
    if(reqt.task_shell_offs_bfn) {
      shell_offs_bfn = basis_map.shell_offs<size_t>( 
        shell_list_bfn.begin(), shell_list_bfn.end() );
      concat_iterable( shell_offs_bfn_pack, shell_offs_bfn );
    }

    // Setup meta data for Shell -> Task (bfn)
    if(reqt.shell_to_task_bfn) {
      const auto itask = std::distance( task_begin, it );
      for( auto i = 0ul; i < nshells_bfn; ++i ) {
        const auto sh_idx = shell_list_bfn.at(i);
        shell_to_task_idx_bfn[sh_idx].emplace_back(itask);
        shell_to_task_off_bfn[sh_idx].emplace_back(shell_offs_bfn.at(i));
      }
    }
  }

  // Send Shell list and offsets (bfn) to device
  if(reqt.task_shell_list_bfn) {
    device_backend_->copy_async( shell_list_bfn_pack.size(), 
      shell_list_bfn_pack.data(), collocation_stack.shell_list_device, 
      "send_shell_list_bfn" );
  } 
  if(reqt.task_shell_offs_bfn) {
    device_backend_->copy_async( shell_offs_bfn_pack.size(), 
      shell_offs_bfn_pack.data(), collocation_stack.shell_offs_device, 
      "send_shell_offs_bfn" );
  }

  auto sl_mem_en = hrt_t::now();
  /*****************************************
   *     GENERATE SHELL -> TASK (bfn)      *
   *****************************************/
  auto s2t_mem_st = hrt_t::now();
  if(reqt.shell_to_task_bfn) {
      
    // Set up buffer allocations from preallocated device segments
    const size_t total_nshells_bfn = 
      total_nshells_bfn_task_batch * sizeof(int32_t);
    buffer_adaptor shell_idx_mem( shell_to_task_stack.shell_to_task_idx_device, 
      total_nshells_bfn );
    buffer_adaptor shell_off_mem( shell_to_task_stack.shell_to_task_off_device, 
      total_nshells_bfn );

    // Reserve memory 
    host_shell_to_task_bfn.resize(global_dims.nshells);

    for( auto ish = 0ul; ish < global_dims.nshells; ++ish ) {
      const auto ntask = shell_to_task_idx_bfn[ish].size();
      auto& bck = host_shell_to_task_bfn[ish];

      // Unpack meta data
      bck.ntask = ntask;
      bck.center_idx = basis_map.shell_to_center( ish );
      bck.true_idx   = ish;
      bck.shell_device = static_stack.shells_device + ish;

      // Allocate device memory 
      bck.task_idx_device = 
        shell_idx_mem.aligned_alloc<int32_t>( ntask, csl );
      bck.task_shell_offs_device = 
        shell_off_mem.aligned_alloc<int32_t>( ntask, csl );

      // Pack host data
      concat_iterable(concat_shell_to_task_idx_bfn, shell_to_task_idx_bfn[ish]);
      concat_iterable(concat_shell_to_task_off_bfn, shell_to_task_off_bfn[ish]);
    }

    // Send data to device
    device_backend_->copy_async( concat_shell_to_task_idx_bfn.size(),
      concat_shell_to_task_idx_bfn.data(), 
      shell_to_task_stack.shell_to_task_idx_device, "shell_to_task_idx_device" );
    device_backend_->copy_async( concat_shell_to_task_off_bfn.size(),
      concat_shell_to_task_off_bfn.data(), 
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
    for( auto i = 0ul; i < global_dims.nshells; ++i ) 
      shell_to_task_sorted[i] = host_shell_to_task_bfn[shell_idx[i]];
    host_shell_to_task_bfn = std::move(shell_to_task_sorted);
    }

    // Send Shell -> Task (bfn) map to device
    device_backend_->copy_async( global_dims.nshells, 
      host_shell_to_task_bfn.data(), shell_to_task_stack.shell_to_task_device,
      "shell_to_task_device" );


    // Form angular momenta batches
    auto max_l = basis_map.max_l();
    l_batched_shell_to_task.resize(max_l + 1);
    auto* p = shell_to_task_stack.shell_to_task_device;
    auto* h = host_shell_to_task_bfn.data();
    for( auto l = 0ul; l <= max_l; ++l ) {
      auto nsh  = basis_map.nshells_with_l(l);
      auto pure = basis_map.l_purity(l);
      l_batched_shell_to_task[l].nshells_in_batch     = nsh;
      l_batched_shell_to_task[l].pure                 = pure;
      l_batched_shell_to_task[l].shell_to_task_device = p;
                          
      size_t max_ntask = std::max_element( h, h+nsh,
        [](auto& a, auto& b){ return a.ntask < b.ntask; } )->ntask;

      l_batched_shell_to_task[l].ntask_average = max_ntask;
      l_batched_shell_to_task[l].npts_average  = 0;

      p += nsh;
      h += nsh;
    }
  
  } // Generate Shell -> Task (bfn)
  auto s2t_mem_en = hrt_t::now();

  /*****************************************
   *   GENERATE SHELLPAIR TO TASK (cou)    *
   *****************************************/
  auto sp2t_mem_st = hrt_t::now();
  if(reqt.shell_pair_to_task_cou or reqt.task_to_shell_pair_cou) {

    const size_t nsp = global_dims.nshell_pairs;
#define USE_TASK_MAP 1

#if !USE_TASK_MAP
    throw std::runtime_error("SP2T + Sparse NYI");
    // Unpack ShellPair SoA (populated in static allocation)
    auto sp2t_1 = hrt_t::now();
    shell_pair_to_task.clear();
    shell_pair_to_task.resize(nsp);
    for( auto i = 0ul; i < nsp; ++i ) {
      auto& sptt = shell_pair_to_task[i];
      sptt.shell_pair_device =
        this->shell_pair_soa.shell_pair_dev_ptr[i];
      std::tie(sptt.lA, sptt.lB) =
        this->shell_pair_soa.shell_pair_ls[i];
      std::tie(sptt.rA, sptt.rB) =
        this->shell_pair_soa.shell_pair_centers[i];
      std::tie(sptt.idx_bra, sptt.idx_ket) =
        this->shell_pair_soa.shell_pair_shidx[i];
    }

    
    // Unpack meta data for cou shell pairs
    auto sp2t_2 = hrt_t::now();
    {
    for( auto it = task_begin; it != task_end; ++it ) {
      const auto& shell_list_cou  = it->cou_screening.shell_list;
      const size_t nshells_cou  = shell_list_cou.size();

      // Compute shell offsets (cou)
      auto shell_offs_cou = basis_map.shell_offs<size_t>( 
        shell_list_cou.begin(), shell_list_cou.end() );
    
      // Setup ShellPair -> Task meta data (cou)
      const auto itask = std::distance( task_begin, it );

      std::map<int,int> sh_off;
      for( auto i = 0ul; i < nshells_cou; ++i )
        sh_off[shell_list_cou[i]] = shell_offs_cou[i];

      for( auto [ish,jsh] : it->cou_screening.shell_pair_list ) {
        const auto idx = detail::packed_lt_index(ish,jsh, global_dims.nshells);

        auto& sptt = shell_pair_to_task[idx];
        sptt.task_idx.emplace_back(itask);
        sptt.task_shell_off_row.emplace_back(sh_off[ish]);
        sptt.task_shell_off_col.emplace_back(sh_off[jsh]);
      }
    }
    }


    // Sort Shell Pairs by diag + L pairs
    auto sp2t_3 = hrt_t::now();
    {
      std::vector<int> sp_idx(nsp);

      // First partition the shell pairs into diagonal and off diagonal
      #if 0
      std::iota( sp_idx.begin(), sp_idx.end(), 0 );
      auto diag_end = std::partition( sp_idx.begin(), sp_idx.end(),
        [=](auto idx) {
          auto [i,j] = detail::from_packed_lt_index(idx, global_dims.nshells);
          return i == j;
        });

      std::sort(sp_idx.begin(), diag_end);
      size_t ioff = 0;
      for( auto i = 0; i < global_dims.nshells; ++i) {
        std::cout << sp_idx[i] << " " << ioff << std::endl;
        ioff += global_dims.nshells - i;
      }
      #else
      auto diag_end = sp_idx.begin() + global_dims.nshells;
      const size_t ns = global_dims.nshells;
      for(size_t j = 0, ioff = 0, odoff = ns; j < ns; ++j) {
        sp_idx[j] = ioff; ++ioff;
        for(size_t i = j + 1; i < ns; ++i) {
          sp_idx[odoff] = ioff; ioff++; odoff++;
        }
      }
      #endif

      // Next sort diag and off diag by shell pair angular momenta
      auto lp_functor = [&](int i, int j){
        const auto& sp_i = shell_pair_to_task[i];
        const auto& sp_j = shell_pair_to_task[j];

        auto l_i = std::make_pair( sp_i.lA, sp_i.lB );
        auto l_j = std::make_pair( sp_j.lA, sp_j.lB );
        return l_i < l_j;
      };

      std::sort( sp_idx.begin(), diag_end, lp_functor );
      std::sort( diag_end, sp_idx.end(),   lp_functor );

      // Reorder shell pairs
      std::vector<ShellPairToTaskHost> sorted_shell_pair_to_task(nsp);
      for( auto i = 0ul; i < nsp; ++i ) {
        sorted_shell_pair_to_task[i] = std::move(shell_pair_to_task[sp_idx[i]]);
      }
      shell_pair_to_task = std::move(sorted_shell_pair_to_task);
    }


    // Set up buffer allocations from preallocated device segments
    const size_t total_nshells_cou_sqlt = 
      total_nshells_cou_sqlt_task_batch * sizeof(int32_t);
    buffer_adaptor shell_idx_mem( 
      shell_pair_to_task_stack.shell_pair_to_task_idx_device, 
      total_nshells_cou_sqlt );
    buffer_adaptor shell_row_off_mem( 
      shell_pair_to_task_stack.shell_pair_to_task_row_off_device, 
      total_nshells_cou_sqlt );
    buffer_adaptor shell_col_off_mem( 
      shell_pair_to_task_stack.shell_pair_to_task_col_off_device, 
      total_nshells_cou_sqlt );

    // Reserve memory 
    auto sp2t_4 = hrt_t::now();
    host_shell_pair_to_task_cou.resize(nsp);
    size_t sptt_sz = std::accumulate(shell_pair_to_task.begin(),shell_pair_to_task.end(),0ul,[](const auto& a, const auto& b){ return a + b.task_idx.size(); });
    concat_shell_pair_to_task_idx_cou.reserve(sptt_sz);
    concat_shell_pair_to_task_off_row_cou.reserve(sptt_sz);
    concat_shell_pair_to_task_off_col_cou.reserve(sptt_sz);

    // Setup ShellPair -> Task (cou) data structures
    for( auto isp = 0ul; isp < nsp; ++isp ) {
      auto& sptt = shell_pair_to_task[isp];
      auto& bck = host_shell_pair_to_task_cou[isp];

      // Unpack meta data
      const auto ntask = sptt.task_idx.size();
      bck.ntask = ntask;
      bck.X_AB = sptt.rA.x - sptt.rB.x;
      bck.Y_AB = sptt.rA.y - sptt.rB.y;
      bck.Z_AB = sptt.rA.z - sptt.rB.z;
      bck.shell_pair_device = sptt.shell_pair_device;

      // Allocate device memory 
      bck.task_idx_device = 
        shell_idx_mem.aligned_alloc<int32_t>( ntask, csl );
      bck.task_shell_off_row_device = 
        shell_row_off_mem.aligned_alloc<int32_t>( ntask, csl );
      bck.task_shell_off_col_device = 
        shell_col_off_mem.aligned_alloc<int32_t>( ntask, csl );

      // Pack host data
      concat_iterable(concat_shell_pair_to_task_idx_cou, sptt.task_idx);
      concat_iterable(
        concat_shell_pair_to_task_off_row_cou, sptt.task_shell_off_row);
      concat_iterable(
        concat_shell_pair_to_task_off_col_cou, sptt.task_shell_off_col);
    }

    // Send packed data to device
    device_backend_->copy_async( concat_shell_pair_to_task_idx_cou.size(),
      concat_shell_pair_to_task_idx_cou.data(), 
      shell_pair_to_task_stack.shell_pair_to_task_idx_device, 
      "shell_pair_to_task_idx_device" );

    device_backend_->copy_async( concat_shell_pair_to_task_off_row_cou.size(),
      concat_shell_pair_to_task_off_row_cou.data(), 
      shell_pair_to_task_stack.shell_pair_to_task_row_off_device, 
      "shell_pair_to_task_row_off_device" );

    device_backend_->copy_async( concat_shell_pair_to_task_off_col_cou.size(),
      concat_shell_pair_to_task_off_col_cou.data(), 
      shell_pair_to_task_stack.shell_pair_to_task_col_off_device, 
      "shell_pair_to_task_col_off_device" );

    // Send ShellPair -> Task (cou) to device 
    device_backend_->copy_async( nsp, 
      host_shell_pair_to_task_cou.data(), 
      shell_pair_to_task_stack.shell_pair_to_task_device,
      "shell_pair_to_task_device" );

    
    const int max_l = basis_map.max_l();

    // Setup DIAGONAL SP AM batches
    auto sp2t_5 = hrt_t::now();
    l_batched_shell_pair_to_task_diag.clear();
    l_batched_shell_pair_to_task_diag.resize(max_l+1);
    {
      auto* p = shell_pair_to_task_stack.shell_pair_to_task_device;
      auto h  = shell_pair_to_task.begin();
      for( auto l = 0; l <= max_l; ++l ) {
        auto& batch_map = l_batched_shell_pair_to_task_diag[l];
        const auto nsh = basis_map.nshells_with_l(l);

        size_t max_npts = 0;
        size_t max_ntask = std::max_element( h, h + nsh,
          [](auto& a, auto& b){
            return a.task_idx.size() < b.task_idx.size();
          })->task_idx.size();

        batch_map.ntask_average = max_ntask;
        batch_map.npts_average  = max_npts;
        batch_map.nshells_in_batch = nsh;
        batch_map.shell_pair_to_task_device = p;
        batch_map.lA = l;
        batch_map.lB = l;

        p += nsh;
        h += nsh;
      }
    } // DIAGONAL SP AM Scope

    // Setup OFFDIAGONAL SP AM batches
    auto sp2t_6 = hrt_t::now();
    l_batched_shell_pair_to_task_off_diag.clear();
    l_batched_shell_pair_to_task_off_diag.resize((max_l+1)*(max_l+1));
    {
      // XXX: NSHELLS offset b/c there are exactly NSHELLS diagonal blocks
      auto* p = shell_pair_to_task_stack.shell_pair_to_task_device +
        global_dims.nshells;
      auto h  = shell_pair_to_task.begin() + global_dims.nshells;
      for( auto l_i = 0, l_ij = 0; l_i <= max_l; ++l_i )
      for( auto l_j = 0; l_j <= max_l; ++l_j, ++l_ij ) {
        auto& batch_map = l_batched_shell_pair_to_task_off_diag[l_ij];
        size_t nsh = std::count_if(
          h, shell_pair_to_task.end(),
          [=](const auto& sp) {
            return sp.lA == l_i and sp.lB == l_j;
          });

        size_t max_npts = 0;
        size_t max_ntask = std::max_element( h, h + nsh,
          [](auto& a, auto& b){
            return a.task_idx.size() < b.task_idx.size();
          })->task_idx.size();

        batch_map.ntask_average = max_ntask;
        batch_map.npts_average  = max_npts;
        batch_map.nshells_in_batch = nsh;
        batch_map.shell_pair_to_task_device = p;
        batch_map.lA = l_i;
        batch_map.lB = l_j;

        p += nsh;
        h += nsh;
      }
    } // OFFDIAGONAL SP AM Scope

    auto sp2t_7 = hrt_t::now();

#if 0
    std::cout << "DIAG SP Batches" << std::endl;
    for( auto& b : l_batched_shell_pair_to_task_diag ) {
      std::cout << b.lA << " " << b.lB << " " << b.nshells_in_batch << " " 
        << b.ntask_average << std::endl;
    }

    std::cout << "OFFDIAG SP Batches" << std::endl;
    for( auto& b : l_batched_shell_pair_to_task_off_diag ) {
      std::cout << b.lA << " " << b.lB << " " << b.nshells_in_batch << " " 
        << b.ntask_average << std::endl;
    }
#endif

  dur_t sp2t_dur_1 = sp2t_2 - sp2t_1;
  dur_t sp2t_dur_2 = sp2t_3 - sp2t_2;
  dur_t sp2t_dur_3 = sp2t_4 - sp2t_3;
  dur_t sp2t_dur_4 = sp2t_5 - sp2t_4;
  dur_t sp2t_dur_5 = sp2t_6 - sp2t_5;
  dur_t sp2t_dur_6 = sp2t_7 - sp2t_6;
  //std::cout << "SP2T 1 = " << sp2t_dur_1.count() << std::endl;
  //std::cout << "SP2T 2 = " << sp2t_dur_2.count() << std::endl;
  //std::cout << "SP2T 3 = " << sp2t_dur_3.count() << std::endl;
  //std::cout << "SP2T 4 = " << sp2t_dur_4.count() << std::endl;
  //std::cout << "SP2T 5 = " << sp2t_dur_5.count() << std::endl;
  //std::cout << "SP2T 6 = " << sp2t_dur_6.count() << std::endl;


#else
  /*****************************************
   *   GENERATE TASK TO SHELLPAIR (cou)    *
   *****************************************/

    auto t2sp_start = hrt_t::now();

    hrt_t::time_point t2sp_1, t2sp_2, t2sp_3, t2sp_4, t2sp_5;

    subtask.clear();
    nprim_pairs_host.clear();
    sp_ptr_host.clear();
    sp_X_AB_host.clear();
    sp_Y_AB_host.clear();
    sp_Z_AB_host.clear();
    task_to_shell_pair.clear();
    l_batch_task_to_shell_pair.clear();
    l_batch_diag_task_to_shell_pair.clear();

    {
      //using point = detail::cartesian_point;
      const int max_l = basis_map.max_l();
      const size_t ntasks = std::distance(task_begin, task_end);

      // Set up task maps for the AM
      for( auto l_i = 0, l_ij = 0; l_i <= max_l; ++l_i )
      for( auto l_j = 0; l_j <= max_l; ++l_j, ++l_ij ) {
        l_batch_task_to_shell_pair.emplace_back();
        auto& batch = l_batch_task_to_shell_pair[l_ij];
        batch.task_to_shell_pair.resize(ntasks);
        batch.lA = l_i;
        batch.lB = l_j;
        batch.max_prim_pairs = 0;
      }

      // Diag terms
      for( auto l_i = 0; l_i <= max_l; ++l_i ) {
        l_batch_diag_task_to_shell_pair.emplace_back();
        auto& batch = l_batch_diag_task_to_shell_pair[l_i];
        batch.task_to_shell_pair.resize(ntasks);
        batch.lA = l_i;
        batch.lB = l_i;
        batch.max_prim_pairs = 0;
      }

      // Generate shell pair device buffer
      nprim_pairs_host = this->shell_pair_soa.shell_pair_nprim_pairs;
      sp_ptr_host = this->shell_pair_soa.shell_pair_dev_ptr;
      for( auto i = 0ul; i < nsp; ++i ) {
        //nprim_pairs_host.push_back(
        //  this->shell_pair_soa.shell_pair_nprim_pairs[i]
        //);
        //sp_ptr_host.push_back(
        //  this->shell_pair_soa.shell_pair_dev_ptr[i]
        //);
        //point rA, rB;
        const auto& [rA, rB] = this->shell_pair_soa.shell_pair_centers[i];

        sp_X_AB_host.push_back(rA.x - rB.x);
        sp_Y_AB_host.push_back(rA.y - rB.y);
        sp_Z_AB_host.push_back(rA.z - rB.z);
      }
    }

    t2sp_1 = hrt_t::now();

    // Total length of the concatenated task map buffers
    size_t task_map_aggregate_length = 0;
    double interim_dur = 0.0;

    {
    const int max_l = basis_map.max_l();

    std::vector<int> sh_off_flat(nsp);
    //std::vector<size_t> sp_idx(nsp);


    //const auto& sp_row_ptr = this->shell_pair_soa.sp_row_ptr;
    //const auto& sp_col_ind = this->shell_pair_soa.sp_col_ind;
    const size_t ntask = std::distance(task_begin,task_end);
    //for( auto it = task_begin; it != task_end; ++it ) {
    //const auto itask = std::distance( task_begin, it );
    for( size_t itask = 0; itask < ntask; ++itask ) {
      auto it = task_begin + itask;

      // Construct the subtasks
      const int points_per_subtask = get_points_per_subtask();
      for (int subtask_i = 0; subtask_i < it->npts; subtask_i += points_per_subtask) {
        subtask.push_back({itask, subtask_i, std::min(it->npts, subtask_i+points_per_subtask), 0});
      }

      // Setup ShellPair offset data
      const auto& shell_list_cou  = it->cou_screening.shell_list;
      const size_t nshells_cou  = shell_list_cou.size();

      // Compute shell offsets (cou)
      auto shell_offs_cou = basis_map.shell_offs<size_t>( 
        shell_list_cou.begin(), shell_list_cou.end() );

      for( auto i = 0ul; i < nshells_cou; ++i )
        sh_off_flat[shell_list_cou[i]] = shell_offs_cou[i];

      // Calculate indices
      //for(auto i = 0ul; i < it->cou_screening.shell_pair_list.size(); ++i) {
      //  auto [ish, jsh] = it->cou_screening.shell_pair_list[i];
      //  sp_idx[i] = detail::csr_index(ish, jsh, global_dims.nshells, sp_row_ptr.data(), sp_col_ind.data()); 
      //}


      // Count the number of shell pairs per task
      for(auto i = 0ul; i < it->cou_screening.shell_pair_list.size(); ++i) {
        auto [ish, jsh] = it->cou_screening.shell_pair_list[i];
        //const auto idx = detail::packed_lt_index(ish,jsh, global_dims.nshells);
        //const auto idx = detail::csr_index(ish, jsh, global_dims.nshells, sp_row_ptr.data(), sp_col_ind.data()); 
        //const auto idx = sp_idx[i];
        const auto idx = it->cou_screening.shell_pair_idx_list[i];;

        int32_t lA, lB;
        std::tie(lA, lB) = this->shell_pair_soa.shell_pair_ls[idx];

        // Filter out diag shell pairs
        if (ish != jsh) {
          const int type_index = lA * (max_l+1) + lB;
          auto& ttsp = l_batch_task_to_shell_pair[type_index].task_to_shell_pair[itask];
          ttsp.nsp++;
          task_map_aggregate_length++;

          l_batch_task_to_shell_pair[type_index].max_prim_pairs = std::max(
            l_batch_task_to_shell_pair[type_index].max_prim_pairs,
            nprim_pairs_host[idx]);

        } else {
          const int type_index = lA;
          auto& ttsp = l_batch_diag_task_to_shell_pair[type_index].task_to_shell_pair[itask];
          ttsp.nsp++;
          task_map_aggregate_length++;

          l_batch_diag_task_to_shell_pair[type_index].max_prim_pairs = std::max(
            l_batch_diag_task_to_shell_pair[type_index].max_prim_pairs,
            nprim_pairs_host[idx]);

        }
      }

      auto _st = hrt_t::now();
      // Allocate space for the shell pair data
      for (auto& batch : l_batch_task_to_shell_pair) {
        for (auto& ttsp : batch.task_to_shell_pair) {
          ttsp.shell_pair_linear_idx.resize(ttsp.nsp);
          ttsp.task_shell_off_row.resize(ttsp.nsp);
          ttsp.task_shell_off_col.resize(ttsp.nsp);
          ttsp.nsp_filled = 0;
        }
      }
      for (auto& batch : l_batch_diag_task_to_shell_pair) {
        for (auto& ttsp : batch.task_to_shell_pair) {
          ttsp.shell_pair_linear_idx.resize(ttsp.nsp);
          ttsp.task_shell_off_row.resize(ttsp.nsp);
          ttsp.task_shell_off_col.resize(ttsp.nsp);
          ttsp.nsp_filled = 0;
        }
      }
      interim_dur += dur_t(hrt_t::now() - _st).count();

      // Iterate over shell pairs adding to tasks
      for(auto i = 0ul; i < it->cou_screening.shell_pair_list.size(); ++i) {
        auto [ish, jsh] = it->cou_screening.shell_pair_list[i];
        //const auto idx = detail::packed_lt_index(ish,jsh, global_dims.nshells);
        //const auto idx = detail::csr_index(ish, jsh, global_dims.nshells, sp_row_ptr.data(), sp_col_ind.data()); 
        //const auto idx = sp_idx[i];
        const auto idx = it->cou_screening.shell_pair_idx_list[i];;

        int32_t lA, lB;
        std::tie(lA, lB) = this->shell_pair_soa.shell_pair_ls[idx];

        // Filter out diag shell pairs
        if (ish != jsh) {
          const int type_index = lA * (max_l+1) + lB;
          auto& ttsp = l_batch_task_to_shell_pair[type_index].task_to_shell_pair[itask];

          const int index = ttsp.nsp_filled++;
          ttsp.shell_pair_linear_idx[index] = idx;
          ttsp.task_shell_off_row[index] = (sh_off_flat[ish] * it->npts);
          ttsp.task_shell_off_col[index] = (sh_off_flat[jsh] * it->npts);
        } else {
          const int type_index = lA;
          auto& ttsp = l_batch_diag_task_to_shell_pair[type_index].task_to_shell_pair[itask];

          const int index = ttsp.nsp_filled++;
          ttsp.shell_pair_linear_idx[index] = idx;
          ttsp.task_shell_off_row[index] = (sh_off_flat[ish] * it->npts);
        }
      }
    }

    }

    t2sp_2 = hrt_t::now();

    // Concat host buffers and copy to device
    buffer_adaptor task_sp_mem( 
      task_to_shell_pair_stack.task_shell_linear_idx_device, 
      total_nshells_cou_sqlt_task_batch * sizeof(int32_t) );
    buffer_adaptor task_row_off_mem( 
      task_to_shell_pair_stack.task_shell_off_row_device, 
      total_nshells_cou_sqlt_task_batch * sizeof(int32_t) );
    buffer_adaptor task_col_off_mem( 
      task_to_shell_pair_stack.task_shell_off_col_device, 
      total_nshells_cou_sqlt_task_batch  * sizeof(int32_t) );

    {
    const size_t ntasks = std::distance(task_begin, task_end);
    const int max_l = basis_map.max_l();
    const int num_sp_types = (max_l+1)*(max_l+1);
    const int num_sp_types_with_diag = num_sp_types + (max_l+1);

    std::vector<TaskToShellPairDevice> host_task_to_shell_pair_task(ntasks * num_sp_types_with_diag);

    std::vector<int32_t> concat_task_to_shell_pair_idx;
    std::vector<int32_t> concat_task_to_shell_pair_off_row;
    std::vector<int32_t> concat_task_to_shell_pair_off_col;

    concat_task_to_shell_pair_idx.reserve(task_map_aggregate_length);
    concat_task_to_shell_pair_off_row.reserve(task_map_aggregate_length);
    concat_task_to_shell_pair_off_col.reserve(task_map_aggregate_length);

    t2sp_3 = hrt_t::now();

    for( auto l_i = 0, l_ij = 0; l_i <= max_l; ++l_i )
    for( auto l_j = 0; l_j <= max_l; ++l_j, ++l_ij ) {
      for( auto itask = 0ul; itask < ntasks; ++itask ) {
        auto& ttsp = l_batch_task_to_shell_pair[l_ij].task_to_shell_pair[itask];

        auto& bck = host_task_to_shell_pair_task[l_ij * ntasks + itask];
        const auto nsp = ttsp.nsp;
        bck.nsp = nsp;

        bck.shell_pair_linear_idx_device = task_sp_mem.aligned_alloc<int32_t>(nsp, csl);
        bck.task_shell_off_row_device = task_row_off_mem.aligned_alloc<int32_t>(nsp, csl);
        bck.task_shell_off_col_device = task_col_off_mem.aligned_alloc<int32_t>(nsp, csl);

        concat_iterable( concat_task_to_shell_pair_idx, ttsp.shell_pair_linear_idx );
        concat_iterable( concat_task_to_shell_pair_off_row, ttsp.task_shell_off_row );
        concat_iterable( concat_task_to_shell_pair_off_col, ttsp.task_shell_off_col );
      }
    }
    for( auto l_i = 0; l_i <= max_l; ++l_i ) {
      for( auto itask = 0ul; itask < ntasks; ++itask ) {
        auto& ttsp = l_batch_diag_task_to_shell_pair[l_i].task_to_shell_pair[itask];

        auto& bck = host_task_to_shell_pair_task[(num_sp_types + l_i) * ntasks + itask];
        const auto nsp = ttsp.nsp;
        bck.nsp = nsp;

        bck.shell_pair_linear_idx_device = task_sp_mem.aligned_alloc<int32_t>(nsp, csl);
        bck.task_shell_off_row_device = task_row_off_mem.aligned_alloc<int32_t>(nsp, csl);
        bck.task_shell_off_col_device = task_col_off_mem.aligned_alloc<int32_t>(nsp, csl);

        concat_iterable( concat_task_to_shell_pair_idx, ttsp.shell_pair_linear_idx );
        concat_iterable( concat_task_to_shell_pair_off_row, ttsp.task_shell_off_row );
        concat_iterable( concat_task_to_shell_pair_off_col, ttsp.task_shell_off_col );
      }
    }
    t2sp_4 = hrt_t::now();

    device_backend_->copy_async( concat_task_to_shell_pair_idx.size(),
      concat_task_to_shell_pair_idx.data(),
      task_to_shell_pair_stack.task_shell_linear_idx_device,
      "task_shell_linear_idx_device");
    device_backend_->copy_async( concat_task_to_shell_pair_off_row.size(),
      concat_task_to_shell_pair_off_row.data(),
      task_to_shell_pair_stack.task_shell_off_row_device,
      "task_shell_off_row_device");
    device_backend_->copy_async( concat_task_to_shell_pair_off_col.size(),
      concat_task_to_shell_pair_off_col.data(),
      task_to_shell_pair_stack.task_shell_off_col_device,
      "task_shell_off_col_device");

    device_backend_->copy_async(host_task_to_shell_pair_task.size(),
      host_task_to_shell_pair_task.data(),
      task_to_shell_pair_stack.task_to_shell_pair_device,
      "task_to_shell_pair_device");

    device_backend_->copy_async(subtask.size(),
      subtask.data(),
      task_to_shell_pair_stack.subtask_device,
      "subtask_device");

    device_backend_->copy_async(nprim_pairs_host.size(),
      nprim_pairs_host.data(),
      task_to_shell_pair_stack.nprim_pairs_device,
      "nprim_pairs_device");

    device_backend_->copy_async(sp_ptr_host.size(),
      sp_ptr_host.data(),
      task_to_shell_pair_stack.sp_ptr_device,
      "sp_ptr_device");
      
    device_backend_->copy_async(sp_X_AB_host.size(),
      sp_X_AB_host.data(),
      task_to_shell_pair_stack.sp_X_AB_device,
      "sp_X_AB_device");

    device_backend_->copy_async(sp_Y_AB_host.size(),
      sp_Y_AB_host.data(),
      task_to_shell_pair_stack.sp_Y_AB_device,
      "sp_Y_AB_device");

    device_backend_->copy_async(sp_Z_AB_host.size(),
      sp_Z_AB_host.data(),
      task_to_shell_pair_stack.sp_Z_AB_device,
      "sp_Z_AB_device");

    t2sp_5 = hrt_t::now();

    l_batch_task_to_shell_pair_device.clear();
    l_batch_task_to_shell_pair_device.resize(num_sp_types);
    for( auto l_i = 0, l_ij = 0; l_i <= max_l; ++l_i )
    for( auto l_j = 0; l_j <= max_l; ++l_j, ++l_ij ) {
      auto& map = l_batch_task_to_shell_pair_device[l_ij];
      map.task_to_shell_pair_device = task_to_shell_pair_stack.task_to_shell_pair_device + l_ij * ntasks;
      map.lA = l_i;
      map.lB = l_j;
      map.max_prim_pairs = l_batch_task_to_shell_pair[l_ij].max_prim_pairs;
    }

    l_batch_diag_task_to_shell_pair_device.clear();
    l_batch_diag_task_to_shell_pair_device.resize(max_l+1);
    for( auto l_i = 0; l_i <= max_l; ++l_i ) {
      auto& map = l_batch_diag_task_to_shell_pair_device[l_i];
      const int offset = (l_i + num_sp_types) * ntasks;
      map.task_to_shell_pair_device = task_to_shell_pair_stack.task_to_shell_pair_device + offset;
      map.lA = l_i;
      map.lB = l_i;
      map.max_prim_pairs = l_batch_diag_task_to_shell_pair[l_i].max_prim_pairs;
    }
    }

    auto t2sp_end = hrt_t::now();

  dur_t t2sp_dur_total = t2sp_end - t2sp_start;
  dur_t t2sp_dur_1 = t2sp_1 - t2sp_start;
  dur_t t2sp_dur_2 = t2sp_2 - t2sp_1;
  dur_t t2sp_dur_3 = t2sp_3 - t2sp_2;
  dur_t t2sp_dur_4 = t2sp_4 - t2sp_3;
  dur_t t2sp_dur_5 = t2sp_5 - t2sp_4;
  dur_t t2sp_dur_6 = t2sp_end - t2sp_5;
  std::cout << "T2SP TOTAL  = " << t2sp_dur_total.count() << std::endl;
  std::cout << "T2SP 1 = " << t2sp_dur_1.count() << std::endl;
  std::cout << "T2SP 2 = " << t2sp_dur_2.count() << std::endl;
  std::cout << "T2SP 3 = " << t2sp_dur_3.count() << std::endl;
  std::cout << "T2SP 4 = " << t2sp_dur_4.count() << std::endl;
  std::cout << "T2SP 5 = " << t2sp_dur_5.count() << std::endl;
  std::cout << "T2SP 6 = " << t2sp_dur_6.count() << std::endl;
  std::cout << "INTERIM = " << interim_dur << std::endl;

#endif

  } // Generate ShellPair -> Task (cou)
  auto sp2t_mem_en = hrt_t::now();
#endif

  dur_t w_mem_dur = w_mem_en - w_mem_st;
  dur_t sl_mem_dur = sl_mem_en - sl_mem_st;
  dur_t s2t_mem_dur = s2t_mem_en - s2t_mem_st;
  dur_t sp2t_mem_dur = sp2t_mem_en - sp2t_mem_st;

  //std::cout << "W DUR    = " << w_mem_dur.count() << std::endl;
  //std::cout << "SL DUR   = " << sl_mem_dur.count() << std::endl;
  //std::cout << "S2T DUR  = " << s2t_mem_dur.count() << std::endl;
  //std::cout << "SP2T DUR = " << sp2t_mem_dur.count() << std::endl;

  auto snd_st = hrt_t::now();
  device_backend_->master_queue_synchronize(); 
  auto snd_en = hrt_t::now();
  //std::cout << "SND_WAIT = " << dur_t(snd_en-snd_st).count() << std::endl;
}








void Scheme1DataBase::add_extra_to_indirection( 
  integrator_term_tracker terms, std::vector<XCDeviceTask>& tasks  ) {

  // Weights Specific
  if( terms.weights ) {
    const auto ldatoms = get_ldatoms();
    buffer_adaptor dist_scratch_mem( scheme1_stack.dist_scratch_device, 
      ldatoms * total_npts_task_batch * sizeof(double) );

    // Extra indirection for dist scratch
    for( auto& task : tasks ) {
      task.dist_scratch  = dist_scratch_mem.aligned_alloc<double>( 
        ldatoms * task.npts, sizeof(double2), csl );
    }
  }

  if( terms.exx or terms.exc_vxc or terms.exc_grad or terms.den or terms.exx_ek_screening ) {
    const size_t total_nshells_bfn = total_nshells_bfn_task_batch * sizeof(size_t);
    buffer_adaptor 
      shell_list_bfn_mem( collocation_stack.shell_list_device, total_nshells_bfn );
    buffer_adaptor 
      shell_offs_bfn_mem( collocation_stack.shell_offs_device, total_nshells_bfn );

    for( auto& task : tasks ) {
      const auto nshells_bfn = task.bfn_screening.nshells;
      task.bfn_screening.shell_list = 
        shell_list_bfn_mem.aligned_alloc<size_t>( nshells_bfn , csl); 
      task.bfn_screening.shell_offs = 
        shell_offs_bfn_mem.aligned_alloc<size_t>( nshells_bfn , csl); 
    }
  }

}

}
