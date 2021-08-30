#include <set>
#include <queue>
#include <future>
#include <cstring>

#include <gauxc/util/unused.hpp>
#include <gauxc/util/misc.hpp>

#include "device/local_work_replicated_incore_exc_vxc.hpp"
#include "device/local_work_replicated_shellbatched_exc_vxc.hpp"
#include "host/util.hpp"
#include "common/integrator_common.hpp"

namespace GauXC  {
namespace integrator::device {

struct dev_ex_task {
  host_task_iterator   task_begin;
  host_task_iterator   task_end;
  std::vector<int32_t> shell_list;
};




dev_ex_task generate_dev_batch( const uint32_t nbf_threshold,
                                host_task_iterator task_begin,
                                host_task_iterator local_work_end,
                                const BasisSet<double>& basis,
                                util::Timer&            timer ) {


  auto nbe_comparator = []( const auto& task_a, const auto& task_b ) {
    return task_a.nbe < task_b.nbe;
  };

  // Find task with largest NBE
  auto max_task = timer.time_op_accumulate("XCIntegrator.MaxTask", [&]() {
    return std::max_element( task_begin, local_work_end, nbe_comparator );
  } );

  const auto max_shell_list = max_task->shell_list; // copy for reset

  // Init uniion shell list to max shell list outside of loop
  std::set<int32_t> union_shell_set(max_shell_list.begin(), 
                                    max_shell_list.end());



  size_t n_overlap_pthresh     = 20;
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
    host_task_iterator search_en = local_work_end;

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
    auto task_end = 
    timer.time_op_accumulate("XCIntegrator.TaskIntersection", [&]() {
      return std::partition( search_st, search_en, [&](const auto& t) {
        return util::integral_list_intersect( max_shell_list, t.shell_list,
                                        overlap_threshold );
      } );
    } );



    // Take union of shell list for all overlapping tasks
    timer.time_op_accumulate("XCIntegrator.ShellListUnion",[&]() {
      for( auto task_it = search_st; task_it != task_end; ++task_it ) {
        local_union_shell_set.insert( task_it->shell_list.begin(), 
                                      task_it->shell_list.end() );
      }
    } );

    auto cur_nbe = basis.nbf_subset( local_union_shell_set.begin(), 
                                     local_union_shell_set.end() );

    //std::cout << "  Threshold %       = " << std::setw(5)  << overlap_pthresh[idx] << ", ";
    //std::cout << "  Overlap Threshold = " << std::setw(8)  << overlap_threshold    << ", ";
    //std::cout << "  Current NBE       = " << std::setw(8)  << cur_nbe              << std::endl;

    // Cache the data
    cached_task_ends[idx] = std::make_pair( task_end, local_union_shell_set );

    // Update partitioned threshold
    cur_partition_pthresh_idx = idx;

    return cur_nbe < nbf_threshold;

  } );

  host_task_iterator task_end;
  auto _idx_partition = (_it == overlap_pthresh_idx.rend()) ? 0 : *_it;
  std::tie( task_end, union_shell_set ) = cached_task_ends.at(_idx_partition);





  //std::cout << "FOUND " << std::distance( task_begin, task_end ) 
  //                      << " OVERLAPPING TASKS" << std::endl;


  std::vector<int32_t> union_shell_list( union_shell_set.begin(),
                                         union_shell_set.end() );

  // Try to add additional tasks given current union list
  task_end = timer.time_op_accumulate("XCIntegrator.SubtaskGeneration", [&]() {
    return std::partition( task_end, local_work_end, [&]( const auto& t ) {
      return util::list_subset( union_shell_list, t.shell_list );
    } );
  } );

  //std::cout << "FOUND " << std::distance( task_begin, task_end ) 
  //                      << " SUBTASKS" << std::endl;


  dev_ex_task ex_task;
  ex_task.task_begin = task_begin;
  ex_task.task_end   = task_end;
  ex_task.shell_list = std::move( union_shell_list );

  return ex_task;

}

template <typename F, size_t n_deriv>
void device_execute_shellbatched(
  util::Timer&           timer,
  XCWeightAlg            weight_alg,
  XCIntegratorState      state,
  const functional_type& func,
  const BasisSet<F>&     basis,
  const Molecule   &     mol,
  const MolMeta    &     meta,
  XCDeviceData<F>  &     device_data,
  const F*               P,
  F*                     VXC,
  F*                     EXC,
  F*                     NEL,
  const dev_ex_task&     ex_task_obj
) {

  // Alias information
  auto task_begin  = ex_task_obj.task_begin;
  auto task_end    = ex_task_obj.task_end;
  auto& union_shell_list = ex_task_obj.shell_list;

  const auto natoms = mol.natoms();

  // Extract subbasis
  BasisSet<F> basis_subset; basis_subset.reserve(union_shell_list.size());
  timer.time_op_accumulate("XCIntegrator.CopySubBasis",[&]() {
    for( auto i : union_shell_list ) {
      basis_subset.emplace_back( basis.at(i) );
    }
    //basis_subset.generate_shell_to_ao();
  });

  // Setup basis maps
  BasisSetMap basis_map( basis );

  const size_t nshells = basis_subset.size();
  const size_t nbe     = basis_subset.nbf();
  std::cout << "TASK_UNION HAS:"   << std::endl
            << "  NSHELLS    = " <<  nshells << std::endl
            << "  NBE        = " <<  nbe     << std::endl;

  // Recalculate shell_list based on subbasis
  timer.time_op_accumulate("XCIntegrator.RecalcShellList",[&]() {
    for( auto _it = task_begin; _it != task_end; ++_it ) {
      auto union_list_idx = 0;
      auto& cur_shell_list = _it->shell_list;
      for( auto j = 0; j < cur_shell_list.size(); ++j ) {
        while( union_shell_list[union_list_idx] != cur_shell_list[j] )
          union_list_idx++;
        cur_shell_list[j] = union_list_idx;
      }
    }
  } );
  


  // Allocate host temporaries
  std::vector<F> P_submat_host(nbe*nbe), VXC_submat_host(nbe*nbe);
  F EXC_tmp, NEL_tmp;
  F* P_submat   = P_submat_host.data();
  F* VXC_submat = VXC_submat_host.data();

  // Extract subdensity
  std::vector<std::array<int32_t,3>> union_submat_cut;
  std::vector<int32_t> foo;
  //auto [union_submat_cut, foo] = 
  std::tie(union_submat_cut,foo) = 
    integrator::gen_compressed_submat_map( basis_map, union_shell_list, 
      basis.nbf(), basis.nbf() );

  timer.time_op_accumulate("XCIntegrator.ExtractSubDensity",[&]() {
    detail::submat_set( basis.nbf(), basis.nbf(), nbe, nbe, P, basis.nbf(), 
                        P_submat, nbe, union_submat_cut );
  } );
 

  // Allocate static quantities on device stack
  device_data.allocate_static_data( natoms, n_deriv, nbe, nshells );

  // Process batches on device with subobjects
  local_work_replicated_incore_exc_vxc_impl<F,n_deriv>(
    weight_alg, state, func, basis_subset, mol, meta, device_data, 
    task_begin, task_end, P_submat, VXC_submat, &EXC_tmp, &NEL_tmp
  );

  // Update full quantities
  *EXC += EXC_tmp;
  *NEL += NEL_tmp;
  timer.time_op_accumulate("XCIntegrator.IncrementSubPotential",[&]() {
    detail::inc_by_submat( basis.nbf(), basis.nbf(), nbe, nbe, VXC, basis.nbf(), 
                           VXC_submat, nbe, union_submat_cut );
  });


  // Reset shell_list to be wrt full basis
  timer.time_op_accumulate("XCIntegrator.ResetShellList",[&]() {
    for( auto _it = task_begin; _it != task_end; ++_it ) 
    for( auto j = 0; j < _it->shell_list.size();  ++j  ) {
      _it->shell_list[j] = union_shell_list[_it->shell_list[j]];
    }
  });

}





template <typename F, size_t n_deriv>
void local_work_replicated_shellbatched_exc_vxc_impl(
  util::Timer&           timer,
  XCWeightAlg            weight_alg,
  XCIntegratorState      state,
  const functional_type& func,
  const BasisSet<F>&     basis,
  const Molecule   &     mol,
  const MolMeta    &     meta,
  XCDeviceData<F>  &     device_data,
  host_task_iterator     local_work_begin,
  host_task_iterator     local_work_end,
  const F*               P,
  F*                     VXC,
  F*                     EXC,
  F*                     NEL
) {

  const uint32_t nbf_threshold = 8000;
  std::cout << "IN SHELL BATCHED\n" << std::flush;
  std::cout << "TOTAL NTASKS = " << std::distance( local_work_begin, local_work_end ) << std:: endl;
  std::cout << "TOTAL NBF    = " << basis.nbf() << std::endl;
  std::cout << "NBF THRESH   = " << nbf_threshold << std::endl;


  // Zero out final results
  timer.time_op( "XCIntegrator.ZeroHost", [&]() {
    *EXC = 0.;
    *NEL = 0.;
    std::memset( VXC, 0, basis.nbf()*basis.nbf()*sizeof(F) );
  });

#if 0
  size_t nbf     = basis.nbf();
  size_t nshells = basis.nshells();
  size_t natoms  = mol.size();

  // Allocate static quantities on device stack
  device_data.allocate_static_data( natoms, n_deriv, nbf, nshells );

  process_batches_cuda_replicated_density_incore_p<F,n_deriv>(
    weight_alg, func, basis, mol, meta, device_data, 
    local_work_begin, local_work_end, P, VXC, EXC, NEL
  );
#else

  auto nbe_comparator = []( const auto& task_a, const auto& task_b ) {
    return task_a.nbe < task_b.nbe;
  };


  size_t batch_iter = 0;
  auto task_begin = local_work_begin;

  const size_t natoms  = mol.size();

  //std::future<void> device_ex;

  std::cout << "MASTER THREAD ID = " << std::this_thread::get_id() << std::endl;
  std::queue< dev_ex_task > dev_tasks;

  auto execute_device_task = [&] () {

    if( dev_tasks.empty() ) return;

    std::cout << "Executing device tasks on thread " << std::this_thread::get_id() << std::endl;

    dev_ex_task batch_task = std::move( dev_tasks.front() ); // Move task to local scope
    dev_tasks.pop(); // Remove from queue
    
    // Execute task
    timer.time_op_accumulate( "XCIntegrator.DeviceWork", [&]() {
      device_execute_shellbatched<F,n_deriv>( timer, weight_alg, state, func, basis, mol,
                                              meta, device_data, P, VXC, EXC, NEL,
                                              batch_task );
    });


  };

  std::future<void> dev_future;
  while( task_begin != local_work_end ) {

    // Generate task
    dev_tasks.emplace( generate_dev_batch( nbf_threshold, task_begin, 
                                           local_work_end, basis, timer ) );

    if( not dev_future.valid() ) {
      dev_future = std::async( std::launch::async, execute_device_task );
    } else {
      auto status = dev_future.wait_for( std::chrono::milliseconds(5) );
      if( status == std::future_status::ready ) {
        dev_future.get();
        dev_future = std::async( std::launch::async, execute_device_task );
      }
    }

    // Update task iterator for next set of batches
    task_begin = dev_tasks.back().task_end;

  }


  if( dev_future.valid() ) dev_future.wait();

  // TODO: Try to merge tasks if possible
  //for( auto _task_it = dev_tasks.begin(); _task_it != dev_tasks.end()-1; ++_task_it ) {
  //  const auto& shell_list = _task_it->union_shell_list;
  //  auto task_nbe = basis.nbf_subset( shell_list.begin(), shell_list.end() );
  //  auto _merge_it = _task_it + 1;
  //  while( task_nbe <= nbf_threshold and _merge_it != dev_tasks.end() ) {
  //    _merge_it = std::find_if( _merge_it, dev_tasks.end(), [&]( const auto& t ) {
  //      const auto& local_shell_list
  //    } );
  //  }
  //}

  while( not dev_tasks.empty() ) {
    // Execute remaining tasks
    execute_device_task();
  }



#endif

}


#define CUDA_IMPL( F, ND ) \
template \
void local_work_replicated_shellbatched_exc_vxc_impl<F,ND>(\
  util::Timer&           timer,\
  XCWeightAlg            weight_alg,\
  XCIntegratorState      state,\
  const functional_type& func,\
  const BasisSet<F>&     basis,\
  const Molecule   &     mol,\
  const MolMeta    &     meta,\
  XCDeviceData<F>  &     device_data,\
  host_task_iterator     local_work_begin,\
  host_task_iterator     local_work_end,\
  const F*               P,\
  F*                     VXC,\
  F*                     exc,\
  F*                     n_el\
) 

CUDA_IMPL( double, 0 );
CUDA_IMPL( double, 1 );

}
}

