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
}

size_t Scheme1DataBase::get_static_mem_requirement() {
  const size_t nsp = (global_dims.nshells*(global_dims.nshells+1))/2;
  return global_dims.nshells * sizeof(ShellToTaskDevice) +
    nsp * sizeof(ShellPairToTaskDevice);
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
    reqt.shell_pair_to_task_col_off_cou_size(nshells_cou) * sizeof(int32_t) ;

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
  for( auto it = task_begin; it != task_end; ++it ) {
    const auto& shell_list_bfn  = it->bfn_screening.shell_list;
    const size_t nshells_bfn  = shell_list_bfn.size();
    total_nshells_bfn_task_batch  += nshells_bfn;

    const auto& shell_list_cou  = it->cou_screening.shell_list;
    const size_t nshells_cou  = shell_list_cou.size();
    const size_t nshells_cou_sqlt = (nshells_cou*(nshells_cou+1))/2;
    total_nshells_cou_sqlt_task_batch  += nshells_cou_sqlt;
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

  // ShellPair -> Task buffer (cou)
  if(reqt.shell_pair_to_task_cou) {
    shell_pair_to_task_stack.shell_pair_to_task_idx_device = 
      mem.aligned_alloc<int32_t>( total_nshells_cou_sqlt_task_batch, csl );
    shell_pair_to_task_stack.shell_pair_to_task_row_off_device = 
      mem.aligned_alloc<int32_t>( total_nshells_cou_sqlt_task_batch, csl );
    shell_pair_to_task_stack.shell_pair_to_task_col_off_device = 
      mem.aligned_alloc<int32_t>( total_nshells_cou_sqlt_task_batch, csl );

    const size_t nsp = (global_dims.nshells*(global_dims.nshells+1))/2;
    shell_pair_to_task_stack.shell_pair_to_task_device =
      mem.aligned_alloc<ShellPairToTaskDevice>( nsp, csl );
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

  /*******************************************
   *         WEIGHTS RELATED MEMORY          *
   *******************************************/

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

  /************************************************
   * SHELL LIST, OFFSET and TASK MAP MEMORY (bfn) *
   ************************************************/

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

  /*****************************************
   *     GENERATE SHELL -> TASK (bfn)      *
   *****************************************/
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

  /*****************************************
   *   GENERATE SHELLPAIR TO TASK (cou)    *
   *****************************************/
  if(reqt.shell_pair_to_task_cou) {

    // Unpack ShellPair SoA (populated in static allocation)
    const size_t nsp = (global_dims.nshells*(global_dims.nshells+1))/2;
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
    for( auto it = task_begin; it != task_end; ++it ) {
      const auto& shell_list_cou  = it->cou_screening.shell_list;
      const size_t nshells_cou  = shell_list_cou.size();

      // Compute shell offsets (cou)
      auto shell_offs_cou = basis_map.shell_offs<size_t>( 
        shell_list_cou.begin(), shell_list_cou.end() );
    
      // Setup ShellPair -> Task meta data (cou)
      const auto itask = std::distance( task_begin, it );
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


    // Sort Shell Pairs by diag + L pairs
    {
      std::vector<int> sp_idx(nsp);

      // First partition the shell pairs into diagonal and off diagonal
      std::iota( sp_idx.begin(), sp_idx.end(), 0 );
      auto diag_end = std::stable_partition( sp_idx.begin(), sp_idx.end(),
        [&](auto idx) {
          auto [i,j] = detail::from_packed_lt_index(idx, global_dims.nshells);
          return i == j;
        });

      // Next sort diag and off diag by shell pair angular momenta
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
    host_shell_pair_to_task_cou.resize(nsp);

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

  } // Generate ShellPair -> Task (cou)
#endif

  device_backend_->master_queue_synchronize(); 
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

  if( terms.exx or terms.exc_vxc or terms.exc_grad or terms.den ) {
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
