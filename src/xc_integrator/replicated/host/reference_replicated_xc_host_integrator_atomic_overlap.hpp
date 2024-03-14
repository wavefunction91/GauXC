/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include "reference_replicated_xc_host_integrator.hpp"
#include "host/local_host_work_driver.hpp"
#include "host/blas.hpp"
#include <stdexcept>
#include <set>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_atomic_overlap_( int64_t iAtom, int64_t m, int64_t n, 
             value_type* S, int64_t lds ) {

  const auto& basis = this->load_balancer_->basis();

  // Check that S is sane
  const int64_t nbf = basis.nbf();
  if( m != n ) 
    GAUXC_GENERIC_EXCEPTION(" S Must Be Square");
  if( m != nbf ) 
    GAUXC_GENERIC_EXCEPTION(" S Must Have Same Dimension as Basis");
  if( lds < nbf )
    GAUXC_GENERIC_EXCEPTION(" Invalid LDS");


  // Get Tasks
  this->load_balancer_->get_tasks();

  // Compute Local contributions to atomic overlap
  this->timer_.time_op("XCIntegrator.LocalWork_AtomicOverlap", [&](){
    atomic_overlap_local_work_(iAtom, S, lds);
  });

  #ifdef GAUXC_ENABLE_MPI
  this->timer_.time_op("XCIntegrator.LocalWait_AtomicOverlap", [&](){
    MPI_Barrier( this->load_balancer_->runtime().comm() );
  });
  #endif

  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce_AtomicOverlap", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");

    this->reduction_driver_->allreduce_inplace( S, nbf*nbf, ReductionOp::Sum );

  });

}










template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  atomic_overlap_local_work_( int64_t iAtom, value_type* S, int64_t lds ) {

  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  // Setup Aliases
  const auto& basis = this->load_balancer_->basis();
  const auto& mol   = this->load_balancer_->molecule();


  // Get basis map
  BasisSetMap basis_map(basis,mol);

  const int32_t nbf = basis.nbf();

  // Sort tasks on size (XXX: maybe doesnt matter?)
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.bfn_screening.nbe) > (b.points.size() * b.bfn_screening.nbe);
  };

  auto& tasks = this->load_balancer_->get_tasks();
  auto task_begin = tasks.begin();
  auto task_end   = tasks.end();
  if(iAtom >= 0) {
    task_end = std::partition(task_begin, task_end, [=](const auto& t){ return t.iParent == iAtom; });
  }
  std::sort( task_begin, task_end, task_comparator );


  // Check that Partition Weights have been calculated
  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored ) {
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Beed Modified"); 
  }

  // Zero out integrands
  for( auto j = 0; j < nbf; ++j )
  for( auto i = 0; i < nbf; ++i ) 
    S[i + j*lds] = 0.;


  const size_t ntasks = std::distance(task_begin, task_end);
  #pragma omp parallel
  {

  XCHostData<value_type> host_data; // Thread local host data

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {

    // Alias current task
    const auto& task = tasks[iT];

    // Get tasks constants
    const int32_t  npts    = task.points.size();
    const int32_t  nbe     = task.bfn_screening.nbe;
    const int32_t  nshells = task.bfn_screening.shell_list.size();

    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();
    const int32_t* shell_list = task.bfn_screening.shell_list.data();

    // Allocate enough memory for batch
    host_data.basis_eval.resize(npts * nbe);
    host_data.nbe_scr.resize(nbe  * nbe);
    host_data.zmat.resize(nbe  * nbe);

    // Alias/Partition out scratch memory
    auto* basis_eval = host_data.basis_eval.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* zmat       = host_data.zmat.data();

    // Get the submatrix map for batch
    std::vector< std::array<int32_t, 3> > submat_map;
    std::tie(submat_map, std::ignore) =
          gen_compressed_submat_map(basis_map, task.bfn_screening.shell_list, nbf, nbf);

    // Evaluate Collocation
    lwd->eval_collocation( npts, nshells, nbe, points, basis, shell_list,
      basis_eval );

    // Copy BFN -> Z
    blas::lacpy('A', npts, nbe, basis_eval, nbe, zmat, nbe);

    // Scale columns of Z by weights
    for(auto j = 0; j < npts; ++j) {
      blas::scal(nbe, 0.5 * weights[j], zmat + j*nbe, nbe); 
    }

    // Incremet LT of S 
    #pragma omp critical
    {
      lwd->inc_vxc( npts, nbf, nbe, basis_eval, submat_map, zmat, nbe, S, lds,
        nbe_scr );
    }

  } // Loop over tasks
  } // End OpenMP Scope

  // Symmetrize VXC
  for( int32_t j = 0;   j < nbf; ++j ) {
    for( int32_t i = j+1; i < nbf; ++i ) {
      S[ j + i*lds ] = S[ i + j*lds ];
    }
  }
}

}
}
