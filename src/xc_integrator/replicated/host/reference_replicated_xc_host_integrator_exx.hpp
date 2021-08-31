#pragma once

#include "reference_replicated_xc_host_integrator.hpp"
#include "integrator_util/integrator_common.hpp"
#include "host/local_host_work_driver.hpp"
#include "host/blas.hpp"
#include <stdexcept>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_exx_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* K, int64_t ldk ) {

  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
  const int64_t nbf = basis.nbf();
  std::string fun_name = __PRETTY_FUNCTION__;
  if( m != n ) 
    throw std::logic_error(fun_name + " P/VXC Must Be Square");
  if( m != nbf ) 
    throw std::logic_error(fun_name + " P/VXC Must Have Same Dimension as Basis");
  if( ldp < nbf )
    throw std::logic_error(fun_name + " Invalid LDP");
  if( ldk < nbf )
    throw std::logic_error(fun_name + " Invalid LDVXC");


  // Get Tasks
  this->load_balancer_->get_tasks();

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    exx_local_work_( P, ldp, K, ldk );
  });


  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      throw std::runtime_error("This Module Only Works With Host Reductions");

    this->reduction_driver_->allreduce_inplace( K, nbf*nbf, ReductionOp::Sum );

  });

}

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  exx_local_work_( const value_type* P, int64_t ldp, 
    value_type* K, int64_t ldk ) {

  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  // Setup Aliases
  const auto& basis = this->load_balancer_->basis();
  const auto& mol   = this->load_balancer_->molecule();
  const auto& meta  = this->load_balancer_->molmeta();

  // Get basis map
  BasisSetMap basis_map(basis,mol);

  const int32_t nbf = basis.nbf();

  // Sort tasks on size (XXX: maybe doesnt matter?)
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.nbe) > (b.points.size() * b.nbe);
  };

  auto& tasks = this->load_balancer_->get_tasks();
  std::sort( tasks.begin(), tasks.end(), task_comparator );


  // Compute Partition Weights
  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored ) {
    lwd->partition_weights( XCWeightAlg::SSF, mol, meta, 
      tasks.begin(), tasks.end() );
    lb_state.modified_weights_are_stored = true;
  }

  // TODO: Get max weight / task + screen

  // Zero out integrands
  for( auto j = 0; j < nbf; ++j )
  for( auto i = 0; i < nbf; ++i ) 
    K[i + j*ldk] = 0.;


  // Loop over tasks
  const size_t ntasks = tasks.size();


  #pragma omp parallel
  {

  XCHostData<value_type> host_data; // Thread local host data

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {

    std::cout << iT << "/" << ntasks << std::endl;
    // Alias current task
    const auto& task = tasks[iT];

    // Get tasks constants
    const int32_t  npts    = task.points.size();

    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();

    std::vector<int32_t> shell_list_( basis.nshells() );
    std::iota( shell_list_.begin(), shell_list_.end(), 0 ); // Don't screen for now
    int32_t* shell_list = shell_list_.data();
    size_t nshells = basis.nshells();
    size_t nbe     = nbf;

    // Allocate data
    host_data.basis_eval.resize( npts * nbe );
    host_data.zmat      .resize( npts * nbe );
    host_data.gmat      .resize( npts * nbe );
    host_data.nbe_scr   .resize( nbe * nbe  );

    // Alias/Partition out scratch memory
    auto* basis_eval = host_data.basis_eval.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* zmat       = host_data.zmat.data();
    auto* gmat       = host_data.gmat.data();


    // Evaluate collocation
    lwd->eval_collocation( npts, nshells, nbe, points, basis, shell_list, 
      basis_eval );

    // Get the submatrix map for batch
    auto [submat_map, foo] = 
      gen_compressed_submat_map( basis_map, shell_list_, nbf, nbf );

    // Evaluate X = P * B -> store in Z
    // XXX: This scales by 2 for total density
    lwd->eval_xmat( npts, nbf, nbe, submat_map, P, ldp, basis_eval, nbe,
      zmat, nbe, nbe_scr );

    // Revert factor of 2 in X
    for( auto i = 0; i < npts * nbe; ++i ) zmat[i] = zmat[i] / 2.;

    // Compute G(mu,i) = w(i) * A(mu,nu,i) * X(nu,i)
    lwd->eval_exx_gmat( npts, nbe, points, weights, basis, basis_map,
      zmat, nbe, gmat, nbe );

    // Increment K
    #pragma omp critical
    blas::gemm( 'N', 'T', nbf, nbf, npts, 1., basis_eval, nbf, 
      gmat, nbf, 1., K, ldk );

  } // Loop over tasks 

  } // End OpenMP region

  // Symmetrize VXC

}

}
}
