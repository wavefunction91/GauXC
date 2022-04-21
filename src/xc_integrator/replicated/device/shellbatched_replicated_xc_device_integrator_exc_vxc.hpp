#include "shellbatched_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"
#include "device/xc_device_aos_data.hpp"
#include "integrator_util/integrator_common.hpp"
#include "host/util.hpp"
#include <gauxc/util/misc.hpp>

#include <stdexcept>
#include <fstream>
#include <queue>
#include <mutex>
#include <future>
#include <set>

namespace GauXC  {
namespace detail {


template <typename ValueType>
void ShellBatchedReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* VXC, int64_t ldvxc,
                 value_type* EXC ) {


  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
  const int64_t nbf = basis.nbf();
  if( m != n ) 
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Be Square");
  if( m != nbf ) 
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Have Same Dimension as Basis");
  if( ldp < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");
  if( ldvxc < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXC");


  // Get Tasks
  auto& tasks = this->load_balancer_->get_tasks();

  // Allocate Device memory
  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
  auto device_data_ptr = 
    this->timer_.time_op("XCIntegrator.DeviceAlloc",
      [=](){ return lwd->create_device_data(); });

  // Generate incore integrator instance, transfer ownership of LWD
  incore_integrator_type incore_integrator( this->func_, this->load_balancer_,
    this->release_local_work_driver(), this->reduction_driver_ );

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  // Compute local contributions to EXC/VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    exc_vxc_local_work_( basis, P, ldp, VXC, ldvxc, EXC, 
      &N_EL, tasks.begin(), tasks.end(), incore_integrator,
      *device_data_ptr );
  });

  // Release ownership of LWD back to this integrator instance
  this->local_work_driver_ = std::move( incore_integrator.release_local_work_driver() );


  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){
    this->reduction_driver_->allreduce_inplace( VXC, nbf*nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( EXC,   1    , ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_EL, 1    , ReductionOp::Sum );
  });


}











template <typename ValueType>
void ShellBatchedReplicatedXCDeviceIntegrator<ValueType>::
  exc_vxc_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                       value_type* VXC, int64_t ldvxc, value_type* EXC, value_type *N_EL,
                       host_task_iterator task_begin, host_task_iterator task_end,
                       incore_integrator_type& incore_integrator, XCDeviceData& device_data ) {


  //incore_integrator.exc_vxc_local_work( basis, P, ldp, VXC, ldvxc, EXC, N_EL, task_begin, task_end, device_data );
  //return;


  const auto nbf = basis.nbf();
  const uint32_t nbf_threshold = 8000;
  const auto& mol = this->load_balancer_->molecule();
  // Zero out integrands on host
  this->timer_.time_op("XCIntegrator.ZeroHost", [&](){
    *EXC  = 0.;
    *N_EL = 0.;
    for( auto j = 0; j < nbf; ++j )
    for( auto i = 0; i < nbf; ++i )
      VXC[i + j*ldvxc] = 0.;
  });


  // Task queue
  std::queue< incore_device_task > incore_device_task_queue;

  // Task queue modification mutex
  std::mutex queue_mod_ex;

  // Lambda for the execution of incore tasks on the device
  auto execute_incore_device_task = [&]() {

    // Early return if there is no task to execute
    if( incore_device_task_queue.empty() ) return;

    incore_device_task next_task;
    {
      std::lock_guard<std::mutex> lock(queue_mod_ex);

      // Move the next task into local scope and remove
      // from queue
      next_task = std::move( incore_device_task_queue.front() );
      incore_device_task_queue.pop();
    }

    // Execute task
    execute_task_batch( next_task, basis, mol, P, ldp, VXC, ldvxc, EXC, N_EL,
      incore_integrator, device_data );
  };


  // Setup future to track execution of currently running
  // device task
  std::future<void> task_future;

  auto task_it = task_begin;
  while( task_it != task_end ) {

    // Generate and enqueue task
    incore_device_task_queue.emplace(
      generate_incore_device_task( nbf_threshold, basis, task_it, task_end )
    );

    // Update iterator for next task generation
    task_it = incore_device_task_queue.back().task_end;

    if( not task_future.valid() ) {
      // No device task to wait on
      task_future = std::async( std::launch::async, execute_incore_device_task );
    } else {
      // Check the status of current device task
      auto status = task_future.wait_for( std::chrono::milliseconds(5) );
      if( status == std::future_status::ready ) {
        // If the status is ready - execute the next task in queue
        task_future.get();
        task_future = std::async( std::launch::async, execute_incore_device_task ); 
      }
    }

  } // Loop until all tasks have been enqued 

  // TODO: Try to merge remaining tasks appropriately

  // Execute remaining tasks sequentially
  if( task_future.valid() ) {
    task_future.wait();
  }
  while( not incore_device_task_queue.empty() ) {
    execute_incore_device_task();
  }
}






template <typename ValueType>
typename ShellBatchedReplicatedXCDeviceIntegrator<ValueType>::incore_device_task 
  ShellBatchedReplicatedXCDeviceIntegrator<ValueType>::
    generate_incore_device_task( const uint32_t     nbf_threshold,
                                 const basis_type&  basis,
                                 host_task_iterator task_begin,
                                 host_task_iterator task_end ) {


  auto nbe_comparator = []( const auto& task_a, const auto& task_b ) {
    return task_a.bfn_screening.nbe < task_b.bfn_screening.nbe;
  };

  // Find task with largest NBE
  auto max_task = this->timer_.time_op_accumulate("XCIntegrator.MaxTask", [&]() {
    return std::max_element( task_begin, task_end, nbe_comparator );
  } );

  const auto max_shell_list = max_task->bfn_screening.shell_list; // copy for reset


  // Init union shell list to max shell list outside of loop
  std::set<int32_t> union_shell_set(max_shell_list.begin(), 
                                    max_shell_list.end());



  int n_overlap_pthresh     = 20;
  double overlap_pthresh_delta = 1. / n_overlap_pthresh;
  std::vector<double> overlap_pthresh;
  for( int i = 1; i < n_overlap_pthresh; ++i )
    overlap_pthresh.emplace_back( i*overlap_pthresh_delta );

  std::vector<int> overlap_pthresh_idx( overlap_pthresh.size() );
  std::iota( overlap_pthresh_idx.begin(), overlap_pthresh_idx.end(), 0 );

  std::map<int, std::pair<host_task_iterator, decltype(union_shell_set)>> 
    cached_task_ends;

  int cur_partition_pthresh_idx = -1;

  auto _it = std::partition_point( overlap_pthresh_idx.rbegin(), 
                                   overlap_pthresh_idx.rend(), 
  [&](int idx) {

    uint32_t overlap_threshold = 
      std::max(1., max_shell_list.size() * overlap_pthresh[idx] );


    host_task_iterator search_st = task_begin;
    host_task_iterator search_en = task_end;

    // Make a local copy of union list
    std::set<int32_t> local_union_shell_set;

    // Attempt to limit task search based on current partition
    if( cur_partition_pthresh_idx >= 0 ) {

      const auto& last_pthresh = 
        cached_task_ends.at(cur_partition_pthresh_idx);

      if( cur_partition_pthresh_idx > idx ) {
        search_st = last_pthresh.first;    
        local_union_shell_set = last_pthresh.second;
      } else {
        search_en = last_pthresh.first;    
        local_union_shell_set = union_shell_set;
      }

    } else {
      local_union_shell_set = union_shell_set;
    }


    // Partition tasks into those which overlap max_task up to
    // specified threshold
    auto local_task_end = 
    this->timer_.time_op_accumulate("XCIntegrator.TaskIntersection", [&]() {
      return std::partition( search_st, search_en, [&](const auto& t) {
        return util::integral_list_intersect( max_shell_list, t.bfn_screening.shell_list,
                                        overlap_threshold );
      } );
    } );



    // Take union of shell list for all overlapping tasks
    this->timer_.time_op_accumulate("XCIntegrator.ShellListUnion",[&]() {
      for( auto task_it = search_st; task_it != local_task_end; ++task_it ) {
        local_union_shell_set.insert( task_it->bfn_screening.shell_list.begin(), 
                                      task_it->bfn_screening.shell_list.end() );
      }
    } );

    auto cur_nbe = basis.nbf_subset( local_union_shell_set.begin(), 
                                     local_union_shell_set.end() );

    //std::cout << "  Threshold %       = " << std::setw(5)  << overlap_pthresh[idx] << ", ";
    //std::cout << "  Overlap Threshold = " << std::setw(8)  << overlap_threshold    << ", ";
    //std::cout << "  Current NBE       = " << std::setw(8)  << cur_nbe              << std::endl;

    // Cache the data
    cached_task_ends[idx] = std::make_pair( local_task_end, local_union_shell_set );

    // Update partitioned threshold
    cur_partition_pthresh_idx = idx;

    return (uint32_t)cur_nbe < nbf_threshold;

  } );

  host_task_iterator local_task_end;
  auto _idx_partition = (_it == overlap_pthresh_idx.rend()) ? 0 : *_it;
  std::tie( local_task_end, union_shell_set ) = cached_task_ends.at(_idx_partition);





  //std::cout << "FOUND " << std::distance( task_begin, local_task_end ) 
  //                      << " OVERLAPPING TASKS" << std::endl;


  std::vector<int32_t> union_shell_list( union_shell_set.begin(),
                                         union_shell_set.end() );

  // Try to add additional tasks given current union list
  local_task_end = this->timer_.time_op_accumulate("XCIntegrator.SubtaskGeneration", [&]() {
    return std::partition( local_task_end, task_end, [&]( const auto& t ) {
      return util::list_subset( union_shell_list, t.bfn_screening.shell_list );
    } );
  } );

  //std::cout << "FOUND " << std::distance( task_begin, local_task_end ) 
  //                      << " SUBTASKS" << std::endl;


  incore_device_task ex_task;
  ex_task.task_begin = task_begin;
  ex_task.task_end   = local_task_end;
  ex_task.shell_list = std::move( union_shell_list );

  return ex_task;

}












template <typename ValueType>
void ShellBatchedReplicatedXCDeviceIntegrator<ValueType>::
  execute_task_batch( incore_device_task& task, const basis_type& basis, const Molecule& mol, 
                      const value_type* P, int64_t ldp, value_type* VXC, int64_t ldvxc, 
                      value_type* EXC, value_type *N_EL, incore_integrator_type& incore_integrator, 
                      XCDeviceData& device_data ) {

  // Alias information
  auto task_begin  = task.task_begin;
  auto task_end    = task.task_end;
  auto& union_shell_list = task.shell_list;


  // Extract subbasis
  BasisSet<double> basis_subset; basis_subset.reserve(union_shell_list.size());
  this->timer_.time_op_accumulate("XCIntegrator.CopySubBasis",[&]() {
    for( auto i : union_shell_list ) {
      basis_subset.emplace_back( basis.at(i) );
    }
  });

  // Setup basis maps
  BasisSetMap basis_map( basis, mol );

  //const size_t nshells = basis_subset.nshells();
  const size_t nbe     = basis_subset.nbf();
  //std::cout << "TASK_UNION HAS:"   << std::endl
  //          << "  NSHELLS    = " <<  nshells << std::endl
  //          << "  NBE        = " <<  nbe     << std::endl;

  // Recalculate shell_list based on subbasis
  this->timer_.time_op_accumulate("XCIntegrator.RecalcShellList",[&]() {
    for( auto _it = task_begin; _it != task_end; ++_it ) {
      auto union_list_idx = 0;
      auto& cur_shell_list = _it->bfn_screening.shell_list;
      for( auto j = 0ul; j < cur_shell_list.size(); ++j ) {
        while( union_shell_list[union_list_idx] != cur_shell_list[j] )
          union_list_idx++;
        cur_shell_list[j] = union_list_idx;
      }
    }
  } );


  // Allocate host temporaries
  std::vector<double> P_submat_host(nbe*nbe), VXC_submat_host(nbe*nbe,0.);
  double EXC_tmp, NEL_tmp;
  double* P_submat   = P_submat_host.data();
  double* VXC_submat = VXC_submat_host.data();



  // Extract subdensity
  std::vector<std::array<int32_t,3>> union_submat_cut;
  std::vector<int32_t> foo;
  std::tie(union_submat_cut,foo) = 
    gen_compressed_submat_map( basis_map, union_shell_list, 
      basis.nbf(), basis.nbf() );

  this->timer_.time_op_accumulate("XCIntegrator.ExtractSubDensity",[&]() {
    detail::submat_set( basis.nbf(), basis.nbf(), nbe, nbe, P, ldp, 
                        P_submat, nbe, union_submat_cut );
  } );


  // Process selected task batch
  incore_integrator.exc_vxc_local_work( basis_subset, P_submat, nbe, VXC_submat, nbe,
    &EXC_tmp, &NEL_tmp, task_begin, task_end, device_data );


  // Update full quantities
  *EXC += EXC_tmp;
  *N_EL += NEL_tmp;
  this->timer_.time_op_accumulate("XCIntegrator.IncrementSubPotential",[&]() {
    detail::inc_by_submat( basis.nbf(), basis.nbf(), nbe, nbe, VXC, ldvxc, 
                           VXC_submat, nbe, union_submat_cut );
  });


  // Reset shell_list to be wrt full basis
  this->timer_.time_op_accumulate("XCIntegrator.ResetShellList",[&]() {
    for( auto _it = task_begin; _it != task_end; ++_it ) 
    for( auto j = 0ul; j < _it->bfn_screening.shell_list.size();  ++j  ) {
      _it->bfn_screening.shell_list[j] = union_shell_list[_it->bfn_screening.shell_list[j]];
    }
  });

}

}
}

