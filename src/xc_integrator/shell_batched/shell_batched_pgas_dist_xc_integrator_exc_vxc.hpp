/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "shell_batched_pgas_dist_xc_integrator.hpp"
#ifdef GAUXC_ENABLE_DEVICE
#include "device/local_device_work_driver.hpp"
#include "device/xc_device_aos_data.hpp"
#endif
#include "integrator_util/integrator_common.hpp"
#include "host/util.hpp"
#include <gauxc/util/misc.hpp>
#include <gauxc/util/unused.hpp>

#include <stdexcept>
#include <fstream>
#include <queue>
#include <mutex>
#include <future>
#include <set>

namespace GauXC  {
namespace detail {


template <typename BaseIntegratorType, typename IncoreIntegratorType>
void ShellBatchedPGASDistributedXCIntegrator<BaseIntegratorType, IncoreIntegratorType>::
  eval_exc_vxc_( const matrix_type& P, matrix_type& VXC, value_type* EXC ) {


  const auto& basis = this->load_balancer_->basis();

  // TODO: Check that P / VXC are sane

  // Get Tasks
  auto& tasks = this->load_balancer_->get_tasks();

#ifdef GAUXC_ENABLE_DEVICE
  // TODO: Allocate Device memory
  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
  auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
  if constexpr (IncoreIntegratorType::is_device) {
    device_data_ptr_ = 
      this->timer_.time_op("XCIntegrator.DeviceAlloc",
        [&](){ return lwd->create_device_data(rt); });
  }
#endif

  // Generate incore integrator instance, transfer ownership of LWD
  incore_integrator_type incore_integrator( this->func_, this->load_balancer_,
    this->release_local_work_driver(), nullptr );

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  // Compute local contributions to EXC/VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    exc_vxc_local_work_( basis, P, VXC, EXC, &N_EL, 
      tasks.begin(), tasks.end(), incore_integrator );
  });

  // Release ownership of LWD back to this integrator instance
  this->local_work_driver_ = std::move( incore_integrator.release_local_work_driver() );

#ifdef GAUXC_ENABLE_DEVICE
  device_data_ptr_.reset();
#endif

}



template <typename BaseIntegratorType, typename IncoreIntegratorType>
void ShellBatchedPGASDistributedXCIntegrator<BaseIntegratorType, IncoreIntegratorType>::
  eval_exc_vxc_( const matrix_type& Ps, const matrix_type& Pz,
                      matrix_type& VXCs, matrix_type& VXCz,
                      value_type* EXC ) {

  GauXC::util::unused(Ps,Pz,VXCs,VXCz,EXC);
  GAUXC_GENERIC_EXCEPTION("UKS NOT YET IMPLEMENTED FOR DEVICE");
}

template <typename BaseIntegratorType, typename IncoreIntegratorType>
void ShellBatchedPGASDistributedXCIntegrator<BaseIntegratorType, IncoreIntegratorType>::
  exc_vxc_local_work_( const basis_type& basis, const matrix_type& P, 
                       matrix_type& VXC, value_type* EXC, value_type *N_EL,
                       host_task_iterator task_begin, host_task_iterator task_end,
                       incore_integrator_type& incore_integrator ) {

  std::cout << "***** MADE IT IN PGAS DRIVER *****" << std::endl;
#if 0
  //incore_integrator.exc_vxc_local_work( basis, P, ldp, VXC, ldvxc, EXC, N_EL, task_begin, task_end, device_data );
  //return;

  std::cout << "TOP SB" << std::endl;

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
  std::queue< incore_task_data > incore_task_data_queue;

  // Task queue modification mutex
  std::mutex queue_mod_ex;

  // Lambda for the execution of incore tasks on the device
  auto execute_incore_task = [&]() {

    // Early return if there is no task to execute
    if( incore_task_data_queue.empty() ) return;

    incore_task_data next_task;
    {
      std::lock_guard<std::mutex> lock(queue_mod_ex);

      // Move the next task into local scope and remove
      // from queue
      next_task = std::move( incore_task_data_queue.front() );
      incore_task_data_queue.pop();
    }

    // Execute task
    execute_task_batch( next_task, basis, mol, P, ldp, VXC, ldvxc, EXC, N_EL,
      incore_integrator );
  };


  // Setup future to track execution of currently running
  // device task
  std::future<void> task_future;

  auto task_it = task_begin;
  while( task_it != task_end ) {

    // Generate and enqueue task
    incore_task_data_queue.emplace(
      generate_incore_task( nbf_threshold, basis, task_it, task_end )
    );

    // Update iterator for next task generation
    task_it = incore_task_data_queue.back().task_end;

    if( not task_future.valid() ) {
      // No device task to wait on
      task_future = std::async( std::launch::async, execute_incore_task );
    } else {
      // Check the status of current device task
      auto status = task_future.wait_for( std::chrono::milliseconds(5) );
      if( status == std::future_status::ready ) {
        // If the status is ready - execute the next task in queue
        task_future.get();
        task_future = std::async( std::launch::async, execute_incore_task ); 
      }
    }

  } // Loop until all tasks have been enqued 

  // TODO: Try to merge remaining tasks appropriately

  // Execute remaining tasks sequentially
  if( task_future.valid() ) {
    task_future.wait();
  }
  while( not incore_task_data_queue.empty() ) {
    execute_incore_task();
  }
#endif
}



template <typename BaseIntegratorType, typename IncoreIntegratorType>
void ShellBatchedPGASDistributedXCIntegrator<BaseIntegratorType, IncoreIntegratorType>::
  execute_task_batch( incore_task_data& task, const basis_type& basis, const Molecule& mol, 
                      const value_type* P, int64_t ldp, value_type* VXC, int64_t ldvxc, 
                      value_type* EXC, value_type *N_EL, incore_integrator_type& incore_integrator ) {

  std::cout << "IN TASK BATCH" << std::endl;

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
#ifdef GAUXC_ENABLE_DEVICE
  if constexpr (IncoreIntegratorType::is_device) {
    incore_integrator.exc_vxc_local_work( basis_subset, P_submat, nbe, VXC_submat, nbe,
      &EXC_tmp, &NEL_tmp, task_begin, task_end, *device_data_ptr_ );
  } else {
#endif
    incore_integrator.exc_vxc_local_work( basis_subset, P_submat, nbe, VXC_submat, nbe,
      &EXC_tmp, &NEL_tmp, task_begin, task_end );
#ifdef GAUXC_ENABLE_DEVICE
  }
#endif


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

