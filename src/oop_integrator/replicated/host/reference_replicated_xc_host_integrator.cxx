#include "reference_replicated_xc_host_integrator.hpp"
#include "integrator_util/integrator_common.hpp"
#include "host/local_host_work_driver.hpp"
#include <stdexcept>

namespace GauXC  {
namespace detail {

template <typename ValueType>
ReferenceReplicatedXCHostIntegrator<ValueType>::~ReferenceReplicatedXCHostIntegrator() noexcept = default;

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* VXC, int64_t ldvxc,
                 value_type* EXC ) {

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
  if( ldvxc < nbf )
    throw std::logic_error(fun_name + " Invalid LDVXC");


  // Get Tasks
  this->load_balancer_->get_tasks();

#if 0
  // Setup Memory
  auto host_data_ptr = this->timer_.time_op("XCIntegrator.HostAlloc",
    [=](){
      //size_t max_npts       = this->load_balancer_->max_npts();
      ////size_t max_nbe        = this->load_balancer_->max_nbe();
      //size_t max_npts_x_nbe = this->load_balancer_->max_npts_x_nbe();
      return std::make_unique<XCHostData<value_type>>(
        //n_deriv, nbf, max_npts, max_npts_x_nbe
      );
    });
#endif

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    exc_vxc_local_work_( P, ldp, VXC, ldvxc, EXC, &N_EL );
  });

#ifdef GAUXC_ENABLE_MPI

  int world_size;
  auto comm = this->load_balancer_->comm();
  MPI_Comm_size( comm, &world_size );

  if( world_size > 1 ) {

  this->timer_.time_op("XCIntegrator.Allreduce", [&](){
    // Test of communicator is an inter-communicator
    // XXX: Can't think of a case when this would be true, but who knows...
    int inter_flag;
    MPI_Comm_test_inter( comm, &inter_flag );

    // Is Intra-communicator, Allreduce can be done inplace
    if( not inter_flag ) {

      MPI_Allreduce( MPI_IN_PLACE, VXC, nbf*nbf, MPI_DOUBLE,
                     MPI_SUM, comm );
      MPI_Allreduce( MPI_IN_PLACE, EXC,  1, MPI_DOUBLE, MPI_SUM, comm );
      MPI_Allreduce( MPI_IN_PLACE, &N_EL, 1, MPI_DOUBLE, MPI_SUM, comm );

    // Isn't Intra-communicator (weird), Allreduce can't be done inplace
    } else {

      std::allocator<value_type> alloc;
      auto VXC_cpy = alloc.allocate( nbf*nbf );
      value_type EXC_cpy = *EXC, N_EL_cpy = N_EL;

      MPI_Allreduce( VXC_cpy, VXC, nbf*nbf, MPI_DOUBLE,
                     MPI_SUM, comm );
      MPI_Allreduce( &EXC_cpy,  EXC,  1, MPI_DOUBLE, MPI_SUM, comm );
      MPI_Allreduce( &N_EL_cpy, &N_EL, 1, MPI_DOUBLE, MPI_SUM, comm );
      

    }
  });

  }

#endif
}

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  exc_vxc_local_work_( const value_type* P, int64_t ldp, 
    value_type* VXC, int64_t ldvxc, value_type* EXC, 
    value_type* N_EL ) {

  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  // Setup Aliases
  const auto& func  = *this->func_;
  const auto& basis = this->load_balancer_->basis();
  const auto& mol   = this->load_balancer_->molecule();
  const auto& meta  = this->load_balancer_->molmeta();

  // Get basis map
  BasisSetMap basis_map(basis);

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
    VXC[i + j*ldvxc] = 0.;
  *EXC = 0.;


  // Loop over tasks
  const size_t ntasks = tasks.size();

  #pragma omp parallel
  {

  XCHostData<value_type> host_data; // Thread local host data

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {

    // Alias current task
    const auto& task = tasks[iT];

    // Get tasks constants
    const int32_t  npts    = task.points.size();
    const int32_t  nbe     = task.nbe;
    const int32_t  nshells = task.shell_list.size();

    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();
    const int32_t* shell_list = task.shell_list.data();

    // Allocate enough memory for batch

    // Things that every calc needs
    host_data.nbe_scr .resize( nbe * nbe  );
    host_data.zmat    .resize( npts * nbe );
    host_data.eps     .resize( npts );
    host_data.vrho    .resize( npts );

    // LDA data requirements
    if( func.is_lda() ){
      host_data.basis_eval .resize( npts * nbe );
      host_data.den_scr    .resize( npts );
    }

    // GGA data requirements
    if( func.is_gga() ){
      host_data.basis_eval .resize( 4 * npts * nbe );
      host_data.den_scr    .resize( 4 * npts );
      host_data.gamma      .resize( npts );
      host_data.vgamma     .resize( npts );
    }

    // Alias/Partition out scratch memory
    auto* basis_eval = host_data.basis_eval.data();
    auto* den_eval   = host_data.den_scr.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* zmat       = host_data.zmat.data();

    auto* eps        = host_data.eps.data();
    auto* gamma      = host_data.gamma.data();
    auto* vrho       = host_data.vrho.data();
    auto* vgamma     = host_data.vgamma.data();

    value_type* dbasis_x_eval = nullptr;
    value_type* dbasis_y_eval = nullptr;
    value_type* dbasis_z_eval = nullptr;
    value_type* dden_x_eval = nullptr;
    value_type* dden_y_eval = nullptr;
    value_type* dden_z_eval = nullptr;

    if( func.is_gga() ) {
      dbasis_x_eval = basis_eval    + npts * nbe;
      dbasis_y_eval = dbasis_x_eval + npts * nbe;
      dbasis_z_eval = dbasis_y_eval + npts * nbe;
      dden_x_eval   = den_eval    + npts;
      dden_y_eval   = dden_x_eval + npts;
      dden_z_eval   = dden_y_eval + npts;
    }


    // Get the submatrix map for batch
    auto [submat_map, foo] = 
      gen_compressed_submat_map( basis_map, task.shell_list, nbf, nbf );

    // Evaluate Collocation (+ Grad)
    if( func.is_gga() )
      lwd->eval_collocation_gradient( npts, nshells, nbe, points, basis, shell_list, 
        basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );
    else
      lwd->eval_collocation( npts, nshells, nbe, points, basis, shell_list, 
        basis_eval );


    // Evaluate X matrix (P * B) -> store in Z
    lwd->eval_xmat( npts, nbf, nbe, submat_map, P, ldp, basis_eval, nbe,
      zmat, nbe, nbe_scr );


    // Evaluate U and V variables
    if( func.is_gga() )
      lwd->eval_uvvar_gga( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
        dbasis_z_eval, zmat, nbe, den_eval, dden_x_eval, dden_y_eval, dden_z_eval,
        gamma );
     else
      lwd->eval_uvvar_lda( npts, nbe, basis_eval, zmat, nbe, den_eval );

    // Evaluate XC functional
    if( func.is_gga() )
      func.eval_exc_vxc( npts, den_eval, gamma, eps, vrho, vgamma );
    else
      func.eval_exc_vxc( npts, den_eval, eps, vrho );

    // Factor weights into XC results
    for( int32_t i = 0; i < npts; ++i ) {
      eps[i]  *= weights[i];
      vrho[i] *= weights[i];
    }

    if( func.is_gga() )
      for( int32_t i = 0; i < npts; ++i ) vgamma[i] *= weights[i];




    // Evaluate Z matrix for VXC
    if( func.is_gga() )
      lwd->eval_zmat_gga_vxc( npts, nbe, vrho, vgamma, basis_eval, dbasis_x_eval,
                              dbasis_y_eval, dbasis_z_eval, dden_x_eval, dden_y_eval,
                              dden_z_eval, zmat, nbe); 
    else
      lwd->eval_zmat_lda_vxc( npts, nbe, vrho, basis_eval, zmat, nbe ); 


    // Incremeta LT of VXC
    #pragma omp critical
    {
      // Scalar integrations
      for( int32_t i = 0; i < npts; ++i ) {
        *N_EL += weights[i] * den_eval[i];
        *EXC  += eps[i]     * den_eval[i];
      }

      // Increment VXC
      lwd->inc_vxc( npts, nbf, nbe, basis_eval, submat_map, zmat, nbe, VXC, ldvxc,
        nbe_scr );
    }

  } // Loop over tasks 

  } // End OpenMP region

  // Symmetrize VXC
  for( int32_t j = 0;   j < nbf; ++j )
  for( int32_t i = j+1; i < nbf; ++i )
    VXC[ j + i*nbf ] = VXC[ i + j*nbf ];

}

template class ReferenceReplicatedXCHostIntegrator<double>;
}
}
