#pragma once

#include "reference_replicated_xc_host_integrator.hpp"
#include "integrator_util/integrator_common.hpp"
#include "integrator_util/integral_bounds.hpp"
#include "host/local_host_work_driver.hpp"
#include "host/blas.hpp"
#include <stdexcept>
#include <set>
#include <atomic>

#include <gauxc/util/geometry.hpp>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_exx_( int64_t m, int64_t n, const value_type* P,
             int64_t ldp, value_type* K, int64_t ldk,
             const IntegratorSettingsEXX& settings ) {

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
    exx_local_work_( P, ldp, K, ldk, settings );
  });

  #ifdef GAUXC_ENABLE_MPI
  this->timer_.time_op("XCIntegrator.LocalWait", [&](){
    MPI_Barrier( this->load_balancer_->comm() );
  });
  #endif

  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      throw std::runtime_error("This Module Only Works With Host Reductions");

    this->reduction_driver_->allreduce_inplace( K, nbf*nbf, ReductionOp::Sum );

  });

}






// MBFS(i) = sqrt(W[i]) * sum_mu B(mu,i)
// return max_i MBFS(i)
double compute_max_bf_sum( size_t npts, size_t nbe_bfn, const double* weights,
  const double* basis_eval, size_t ldb ) {

  std::vector<double> bf_sums( npts );
  for( auto ipt = 0ul; ipt < npts; ++ipt ) {
    double tmp = 0.;
    for( auto ibf = 0ul; ibf < nbe_bfn; ++ibf ) 
      tmp += std::abs( basis_eval[ibf + ipt * ldb] );
    bf_sums[ipt] = std::sqrt(weights[ipt]) * tmp;
  }

  return *std::max_element( bf_sums.begin(), bf_sums.end() );

}


auto compute_approx_f_max( size_t npts, size_t nshells_bf, size_t nbf, 
  size_t nbe_bfn, const BasisSetMap& basis_map, const double* weights, 
  const std::vector<std::array<int32_t,3>>& submat_bfn, const double* basis_eval,
  size_t ldb, const double* P_abs, size_t ldp, LocalHostWorkDriver* lwd,
  double* nbe_scr) {

  // Get max value for each bfn over grid 
  std::vector<double> max_bf_grid( nbe_bfn );
  for( auto ibf = 0ul; ibf < nbe_bfn; ++ibf ) {
    double tmp = 0.;
    for( auto ipt = 0ul; ipt < npts; ++ipt )
      tmp = std::max( tmp,
        std::sqrt(weights[ipt]) *
        std::abs(basis_eval[ibf + ipt*ldb])
      );
    max_bf_grid[ibf] = tmp;
  }

  // Compute approximate F max over basis functions
  std::vector<double> max_F_approx_bfn( nbf );
  std::decay_t<decltype(submat_bfn)> submat_full = {
    {0, nbf, 0}
  };

  lwd->eval_exx_fmat( 1, nbf, nbf, nbe_bfn, submat_full, submat_bfn,
    P_abs, ldp, max_bf_grid.data(), nbe_bfn, max_F_approx_bfn.data(),
    nbf, nbe_scr );

  // Collapse approx F max over shells 
  std::vector<double> max_F_approx( nshells_bf );
  for( auto ish = 0ul; ish < nshells_bf; ++ish ) {
    const auto sh_st = basis_map.shell_to_first_ao(ish);
    const auto sh_sz = basis_map.shell_size(ish);
    double tmp = 0.;
    for( auto i = sh_st; i < sh_st + sh_sz; ++i )
      tmp = std::max( tmp, std::abs(max_F_approx_bfn[i]) );
    max_F_approx[ish] = tmp;
  }

  return max_F_approx;
}



auto compute_true_f_max( size_t npts, size_t nshells_bra, size_t nbe_bra,
  const BasisSetMap& basis_map, const std::vector<int32_t>& shell_list_bra,
  const double* weights, const double* F, size_t ldf ) {

  std::vector<double> max_F( nshells_bra );
  size_t sh_st = 0;
  for( auto i = 0ul; i < nshells_bra; ++i ) {
    const auto ish = shell_list_bra[i];
    const auto sh_sz = basis_map.shell_size(ish);

    double tmp_max = 0.;
    for( auto ipt = 0ul; ipt < npts; ++ipt ) 
    for( auto ii = 0ul;   ii < sh_sz;  ++ii  )
      tmp_max = std::max( tmp_max,
        std::sqrt(weights[ipt]) *
        std::abs( F[ sh_st + ii + ipt*ldf] )
      );
    max_F[i] = tmp_max;

    sh_st += sh_sz;
  }

  return max_F;

}



auto compute_sn_LinK_E_set( size_t nshells, const std::vector<int32_t>& shell_list,
  const double* V_max, size_t ldv, const double* max_F, double E_tol ) {

  std::set<int32_t> E_shells;
  for( auto i = 0ul; i < nshells; ++i )
  for( auto j = 0ul; j <= i;      ++j ) {

    const auto ish = shell_list[i];
    const auto jsh = shell_list[j];

    const auto V_ij = V_max[ish + jsh*ldv];
    const auto F_i  = max_F[i];
    const auto F_j  = max_F[j];

    const double eps_E_compare = F_i * F_j * V_ij;
    if( eps_E_compare > E_tol )  {
      E_shells.insert(ish); 
      E_shells.insert(jsh); 
    }

  }

  return E_shells;

}

auto compute_sn_LinK_K_set( size_t nshells, const std::vector<int32_t>& shell_list,
  const double* V_max, size_t ldv, const double* max_F, double max_bf_sum,
  double K_tol ) {

  std::set<int32_t> K_shells;
  for( auto i = 0ul; i < nshells; ++i )
  for( auto j = 0ul; j <= i;      ++j ) {

    const auto ish = shell_list[i];
    const auto jsh = shell_list[j];

    const auto V_ij = V_max[ish + jsh*ldv];
    const auto F_i  = max_F[i];
    const auto F_j  = max_F[j];

    const double eps_K_compare = std::max(F_i, F_j) * V_ij * max_bf_sum;
    if( eps_K_compare > K_tol )  {
      K_shells.insert(ish); 
      K_shells.insert(jsh); 
    }

  }

  return K_shells;

}

auto compute_sn_LinK_ek_set( size_t nshells, const std::vector<int32_t>& shell_list,
  const double* V_max, size_t ldv, const double* max_F, double max_bf_sum,
  double E_tol, double K_tol ) {

  std::set<int32_t> ek_shells;
  for( auto i = 0ul; i < nshells; ++i )
  for( auto j = 0ul; j <= i;      ++j ) {

    const auto ish = shell_list[i];
    const auto jsh = shell_list[j];

    const auto V_ij = V_max[ish + jsh*ldv];
    const auto F_i  = max_F[i];
    const auto F_j  = max_F[j];

    const double eps_E_compare = F_i * F_j * V_ij;
    const double eps_K_compare = std::max(F_i, F_j) * V_ij * max_bf_sum;
    if( eps_K_compare > K_tol or eps_E_compare > E_tol)  {
      ek_shells.insert(ish); 
      ek_shells.insert(jsh); 
    }

  }

  return ek_shells;

}














template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  exx_local_work_( const value_type* P, int64_t ldp, 
    value_type* K, int64_t ldk, const IntegratorSettingsEXX& settings ) {

  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  // Setup Aliases
  const auto& basis = this->load_balancer_->basis();
  const auto& mol   = this->load_balancer_->molecule();
  const auto& meta  = this->load_balancer_->molmeta();

  // XXX Check that basis is cartesian
  //for( auto& sh : basis ) if(sh.pure()) throw std::runtime_error("sn-LinK EXX Only Works With Cartesian Functions");

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

  // Zero out integrands
  for( auto j = 0; j < nbf; ++j )
  for( auto i = 0; i < nbf; ++i ) 
    K[i + j*ldk] = 0.;

   
  // Compute V upper bounds per shell pair
  const size_t nshells_bf = basis.size();
  std::vector<double> V_max( nshells_bf * nshells_bf );
  for( auto i = 0; i < nshells_bf; ++i )
  for( auto j = 0; j < nshells_bf; ++j ) {
    V_max[i + j*nshells_bf] = util::max_coulomb( basis.at(i), basis.at(j) );
  }

  // Absolute value of P
  std::vector<double> P_abs(nbf*nbf);
  for( auto i = 0; i < nbf*nbf; ++i ) P_abs[i] = std::abs(P[i]);

  // Loop over tasks
  const size_t ntasks = tasks.size();
  std::atomic_size_t nskip = 0;

  // Full shell list
  std::vector<int32_t> full_shell_list_( basis.nshells() );
  std::iota( full_shell_list_.begin(), full_shell_list_.end(), 0 );
  std::vector< std::array<int32_t,3> > full_submat_map = { {0, nbf, 0} };

  // Screening settings
  IntegratorSettingsSNLinK sn_link_settings;
  if( auto* tmp = dynamic_cast<const IntegratorSettingsSNLinK*>(&settings) ) {
    sn_link_settings = *tmp;
  }

  const bool screen_ek = sn_link_settings.screen_ek;
  const double eps_K   = sn_link_settings.k_tol;
  const double eps_E   = sn_link_settings.energy_tol;

  int world_rank = 0;
  #ifdef GAUXC_ENABLE_MPI
  auto comm = this->load_balancer_->comm();
  MPI_Comm_rank( comm, &world_rank );
  #endif
  if( !world_rank ) {
    std::cout << "sn-LinK Settings:" << std::endl
              << "  SCREEN_EK     = " << std::boolalpha << screen_ek << std::endl
              << "  EPS_E         = " << eps_E << std::endl
              << "  EPS_K         = " << eps_K << std::endl
              << std::endl;
  }

  #pragma omp parallel
  {

  XCHostData<value_type> host_data; // Thread local host data

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {

    //std::cout << iT << "/" << ntasks << std::endl;
    // Alias current task
    const auto& task = tasks[iT];

    // Get tasks constants
    const int32_t  npts    = task.points.size();

    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();

    // Basis function shell list
    auto shell_list_bfn_ = task.shell_list;
    int32_t* shell_list_bfn = shell_list_bfn_.data();
    size_t nshells_bfn = shell_list_bfn_.size();
    size_t nbe_bfn     = 
      basis.nbf_subset( shell_list_bfn_.begin(), shell_list_bfn_.end() );
    auto [submat_map_bfn, foo1] = 
      gen_compressed_submat_map( basis_map, shell_list_bfn_, nbf, nbf );
    


    // Allocate data screening independent data
    host_data.basis_eval.resize( npts * nbe_bfn );
    host_data.nbe_scr   .resize( nbe_bfn * nbf );
    auto* basis_eval = host_data.basis_eval.data();
    auto* nbe_scr    = host_data.nbe_scr.data();



    // Evaluate collocation B(mu,i)
    // mu ranges over the bfn shell list and i runs over all points
    lwd->eval_collocation( npts, nshells_bfn, nbe_bfn, points, basis, 
      shell_list_bfn, basis_eval );

    // Compute Max BF Sum
    auto max_bf_sum = 
      compute_max_bf_sum( npts, nbe_bfn, weights, basis_eval, nbe_bfn );

    // Compute Approximate max F
    auto max_F_approx =
      compute_approx_f_max( npts, nshells_bf, nbf, nbe_bfn, basis_map,
        weights, submat_map_bfn, basis_eval, nbe_bfn, P_abs.data(), nbf,
        lwd, nbe_scr );

    // Get shell pair screening for integrals
    auto ek_shell_set = 
      compute_sn_LinK_ek_set( nshells_bf, full_shell_list_,
        V_max.data(), nshells_bf, max_F_approx.data(), max_bf_sum,
        eps_E, eps_K );

    // Bail on task if no shells are needed after ek screening
    if( ek_shell_set.size() == 0 ) {
      nskip++; // this is atomic
      continue;
    }


    // Convert ek shell set -> vector
    std::vector<int32_t> ek_shell_list;
    std::vector< std::array<int32_t,3> > ek_submat_map;
    if( ek_shell_list.size() == nshells_bf or !screen_ek ) {
      ek_shell_list = full_shell_list_;
      ek_submat_map = full_submat_map;
    } else {
      ek_shell_list = decltype(ek_shell_list)( ek_shell_set.begin(), ek_shell_set.end() );
      std::tie( ek_submat_map, std::ignore ) =
        gen_compressed_submat_map( basis_map, ek_shell_list, nbf, nbf );
    }
    const auto nbe_ek = basis.nbf_subset( ek_shell_list.begin(), ek_shell_list.end() );
    const auto nshells_ek = ek_shell_list.size();


    // Allocate Screening Dependent Data
    host_data.zmat.resize( npts * nbe_ek );
    host_data.gmat.resize( npts * nbe_ek );
    auto* zmat = host_data.zmat.data();
    auto* gmat = host_data.gmat.data();

    // Evaluate F(mu,i) = P(mu,nu) * B(nu,i)
    // mu runs over significant ek shells
    // nu runs over the bfn shell list
    // i runs over all points
    lwd->eval_exx_fmat( npts, nbf, nbe_ek, nbe_bfn, ek_submat_map,
      submat_map_bfn, P, ldp, basis_eval, nbe_bfn, zmat, nbe_ek, nbe_scr );

    // Get True Max F for shell pairs
    //auto max_F = compute_true_f_max( npts, nshells_ek, nbe_ek, basis_map,
    //  ek_shell_list, weights, zmat, nbe_ek );


    // Compute G(mu,i) = w(i) * A(mu,nu,i) * F(nu,i)
    // mu/nu run over significant ek shells
    // i runs over all points
    lwd->eval_exx_gmat( npts, nshells_ek, nbe_ek, points, weights, 
      basis, basis_map, ek_shell_list.data(), zmat, nbe_ek, gmat, nbe_ek );

    // Increment K(mu,nu) += B(mu,i) * G(nu,i)
    // mu runs over bfn shell list
    // nu runs over ek shells
    // i runs over all points
    #pragma omp critical
    lwd->inc_exx_k( npts, nbf, nbe_bfn, nbe_ek, basis_eval, submat_map_bfn,
      ek_submat_map, gmat, nbe_ek, K, ldk, nbe_scr );

  } // Loop over tasks 

  } // End OpenMP region

  // Symmetrize K
  for( auto j = 0; j < nbf; ++j ) 
  for( auto i = 0; i < j;   ++i ) {
    const auto K_ij = K[i + j*ldk];
    const auto K_ji = K[j + i*ldk];
    const auto K_symm = 0.5 * (K_ij + K_ji);
    K[i + j*ldk] = K_symm;
    K[j + i*ldk] = K_symm;
  }

  std::cout << "NSKIP = " << nskip << " / " << ntasks << std::endl;
}

}
}
