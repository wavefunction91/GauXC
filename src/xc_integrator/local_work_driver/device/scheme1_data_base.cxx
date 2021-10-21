#include "scheme1_data_base.hpp"
#include "buffer_adaptor.hpp"

namespace GauXC {

Scheme1DataBase::~Scheme1DataBase() noexcept = default;

Scheme1DataBase::Scheme1DataBase(std::unique_ptr<DeviceBackend>&& ptr ) : 
  base_type(std::move(ptr)) {

  if( device_backend_ ) 
    device_backend_->create_blas_queue_pool(4);

}

void Scheme1DataBase::reset_allocations() {
  base_type::reset_allocations();
  scheme1_stack.reset();
}

size_t Scheme1DataBase::get_static_mem_requirement() {
  return global_dims.nshells * sizeof(ShellToTaskDevice);
}









size_t Scheme1DataBase::get_mem_req( integrator_term_tracker terms, 
  const host_task_type& task ) {

  // All local memory is weights related
  size_t base_size = base_type::get_mem_req(terms, task);

  if( terms.weights ) {
    const auto ldatoms = get_ldatoms();
    const auto mem_dist_scr = ldatoms * task.npts;
    const auto mem_dist_ner = task.npts;
    const auto mem_iparent  = task.npts;

    base_size += (mem_dist_scr + mem_dist_ner) * sizeof(double) + 
                 mem_iparent * sizeof(int32_t);
  }

  if( terms.exc_vxc or terms.exc_grad ) {
    const auto& shell_list = task.shell_list;
    const size_t nshells  = shell_list.size();

    const size_t mem_shell_list = nshells * sizeof(size_t);
    const size_t mem_shell_offs = nshells * sizeof(size_t);

    // Shell -> Task maps
    const size_t mem_shell_to_task_idx  = nshells * sizeof(int32_t);
    const size_t mem_shell_to_task_off  = nshells * sizeof(int32_t);

    base_size +=  mem_shell_list + mem_shell_offs +
      mem_shell_to_task_idx + mem_shell_to_task_off;
  }

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

  if( terms.exc_vxc or terms.exc_grad ) {
    total_nshells_task_batch  = 0; 
    for( auto it = task_begin; it != task_end; ++it ) {
      const auto& shell_list  = it->shell_list;
      const size_t nshells  = shell_list.size();
      total_nshells_task_batch  += nshells;
    }
    collocation_stack.shell_list_device = 
      mem.aligned_alloc<size_t>( total_nshells_task_batch , csl);
    collocation_stack.shell_offs_device = 
      mem.aligned_alloc<size_t>( total_nshells_task_batch , csl);

    // Shell -> Task buffers
    shell_to_task_stack.shell_to_task_idx_device = 
      mem.aligned_alloc<int32_t>( total_nshells_task_batch, csl );
    shell_to_task_stack.shell_to_task_off_device = 
      mem.aligned_alloc<int32_t>( total_nshells_task_batch, csl );
    shell_to_task_stack.shell_to_task_device =
      mem.aligned_alloc<ShellToTaskDevice>( global_dims.nshells, csl );

  }

  // Update dynmem data for derived impls
  return device_buffer_t{ mem.stack(), mem.nleft() };
}









void Scheme1DataBase::pack_and_send( 
  integrator_term_tracker terms,
  host_task_iterator task_begin, host_task_iterator task_end,
  const BasisSetMap& basis_map ) {

  // Pack and send base data
  base_type::pack_and_send( terms, task_begin, task_end, basis_map );

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

  if( terms.exc_vxc or terms.exc_grad ) {
    // Contatenation utility
    auto concat_iterable = []( auto& a, const auto& b ) {
      a.insert( a.end(), b.begin(), b.end() );
    };

    std::vector< size_t > shell_list_pack;
    std::vector< size_t > shell_offs_pack;
    std::vector< std::vector<int32_t> > shell_to_task_idx( global_dims.nshells ),
                                        shell_to_task_off( global_dims.nshells );

    for( auto it = task_begin; it != task_end; ++it ) {
      const auto& shell_list  = it->shell_list;
      const size_t nshells  = shell_list.size();

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

      concat_iterable( shell_list_pack, shell_list );
      concat_iterable( shell_offs_pack, shell_offs );
    }

    device_backend_->copy_async( shell_list_pack.size(), shell_list_pack.data(), 
      collocation_stack.shell_list_device, "send_shell_list" );
    device_backend_->copy_async( shell_offs_pack.size(), shell_offs_pack.data(), 
      collocation_stack.shell_offs_device, "send_shell_offs" );

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
      l_batched_shell_to_task.resize(max_l + 1);
      {
      auto* p = shell_to_task_stack.shell_to_task_device;
      auto* h = host_shell_to_task.data();
      for( auto l = 0; l <= max_l; ++l ) {
        auto nsh  = basis_map.nshells_with_l(l);
        auto pure = basis_map.l_purity(l);
        l_batched_shell_to_task[l].nshells_in_batch     = nsh;
        l_batched_shell_to_task[l].pure                 = pure;
        l_batched_shell_to_task[l].shell_to_task_device = p;
        size_t max_npts = 0;
        //for( auto ish = 0; ish < global_dims.nshells; ++ish ) 
        //if( basis_map.shell_l(ish) == l ) {
        //  if( not shell_to_task_idx[ish].size() ) continue;
        //  auto task_idx = 
        //    *std::max_element( shell_to_task_idx[ish].begin(),
        //                       shell_to_task_idx[ish].end(),
        //                       [&](auto i, auto j) {
        //                         return (task_begin+i)->points.size() <
        //                                (task_begin+j)->points.size();
        //                       } );
        //  max_npts = std::max( max_npts, (task_begin+task_idx)->points.size() );
       
        //}

        size_t total_ntask = std::accumulate( h, h + nsh, 0ul,
          [](auto& a, auto& b){ return a + b.ntask; } );
        size_t max_ntask = std::max_element( h, h+nsh,
          [](auto& a, auto& b){ return a.ntask < b.ntask; } )->ntask;
        std::cout << "L = " << l << " AVG = " << (total_ntask/nsh) << 
          " MAX = " << max_ntask << ", " << max_npts << std::endl;
        //l_batched_shell_to_task[l].ntask_average = total_ntask / nsh;
        l_batched_shell_to_task[l].ntask_average = max_ntask;
        l_batched_shell_to_task[l].npts_average  = max_npts;

        p += nsh;
        h += nsh;
      }
      }
    } // Shell -> Task data
  }

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

  if( terms.exc_vxc or terms.exc_grad ) {
    const size_t total_nshells = total_nshells_task_batch * sizeof(size_t);
    buffer_adaptor 
      shell_list_mem( collocation_stack.shell_list_device, total_nshells );
    buffer_adaptor 
      shell_offs_mem( collocation_stack.shell_offs_device, total_nshells );
    for( auto& task : tasks ) {
      const auto nshells = task.nshells;
      task.shell_list   = shell_list_mem.aligned_alloc<size_t>( nshells , csl); 
      task.shell_offs   = shell_offs_mem.aligned_alloc<size_t>( nshells , csl); 
    }
  }

}

}
