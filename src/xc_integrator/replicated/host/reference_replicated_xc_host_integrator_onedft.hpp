/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include "reference_replicated_xc_host_integrator.hpp"
#include "integrator_util/integrator_common.hpp"
#include "integrator_util/onedft_util.hpp"
#include "host/local_host_work_driver.hpp"
#include "host/blas.hpp"

#include <stdexcept>
#include <string>

namespace GauXC  {
namespace detail {
  
FeatureDict prepare_onedft_features(const int ndm, std::vector<XCTask>& tasks, const Molecule& mol, 
  const std::vector<std::string> feature_keys, const RuntimeEnvironment& rt, std::vector<int>& sendcounts, 
  std::vector<int>& displs, std::vector<int64_t>& atom_reorder_inv_perm);

void send_buffer_onedft_outputs(const int ndm, const FeatureDict features_dict, std::vector<XCTask>& tasks, 
  const RuntimeEnvironment& rt, std::vector<int> sendcounts, std::vector<int> displs,
  const std::vector<int64_t>& atom_reorder_inv_perm);

void interleave_data(const double* a, const double* b, const size_t n, double* out);

void eval_zmat_gga_vxc_uks(size_t npts, size_t nbf, 
  const double* vdden_eval_a, const double* vdden_eval_b, 
  const double* vdden_x_eval_a, const double* vdden_x_eval_b, const double* vdden_y_eval_a, 
  const double* vdden_y_eval_b, const double* vdden_z_eval_a, const double* vdden_z_eval_b,
  const double* basis_eval, const double* dbasis_x_eval, const double* dbasis_y_eval, const double* dbasis_z_eval,
  double* Zs, size_t ldzs, double* Zz, size_t ldzz);

void eval_zmat_mgga_vxc_uks(size_t npts, size_t nbf, 
  const double* vdden_eval_a, const double* vdden_eval_b, 
  const double* vlapl_a, const double* vlapl_b,
  const double* vdden_x_eval_a, const double* vdden_x_eval_b, const double* vdden_y_eval_a, 
  const double* vdden_y_eval_b, const double* vdden_z_eval_a, const double* vdden_z_eval_b,
  const double* basis_eval, const double* dbasis_x_eval, const double* dbasis_y_eval, const double* dbasis_z_eval,
  const double* lbasis_eval, double* Zs, size_t ldzs, double* Zz, size_t ldzz);

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_exc_vxc_onedft_( int64_t m, int64_t n, 
    const value_type* Ps, int64_t ldps,
    const value_type* Pz, int64_t ldpz,
    value_type* VXCs, int64_t ldvxcs,
    value_type* VXCz, int64_t ldvxcz,
    value_type* EXC, const IntegratorSettingsXC& settings ) {

  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
  const int64_t nbf = basis.nbf();
  if( m != n ) 
    GAUXC_GENERIC_EXCEPTION("P Must Be Square");
  if( m != nbf ) 
    GAUXC_GENERIC_EXCEPTION("P Must Have Same Dimension as Basis");

  if( ldps and ldps < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");
  if( ldpz and ldpz < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDPZ");

  if( ldvxcs < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCS");
  if( ldvxcz and ldvxcz < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCZ");

  // const bool is_exc_only = (!VXCs) and (!VXCz) and (!VXCy) and (!VXCx);

  const bool is_uks = (Pz != nullptr);
  if (not is_uks) {
    // TODO: duplicate the density matrix ot duplicate the feature results?
    // Pz = Ps;
    GAUXC_GENERIC_EXCEPTION("RKS Not Yet Implemented");
  }
  // Get Tasks
  auto& tasks = this->load_balancer_->get_tasks();  
  auto rt = this->load_balancer_->runtime();
  int32_t world_rank = rt.comm_rank();
  // Temporary electron count to judge integrator accuracy
  value_type N_EL; 
  
  // load model from parameter/interator settings
  OneDFTSettings onedft_settings;
  if( auto* tmp = dynamic_cast<const OneDFTSettings*>(&settings) ) {
    onedft_settings = *tmp;
  }
  const auto model_path = onedft_settings.model;
  torch::DeviceType device = torch::kCPU;
  auto [exc_func, feature_keys] = load_model(model_path, device);
  
  // determine what feature we need based on the keys
  if (feature_keys.size() == 0) {
    GAUXC_GENERIC_EXCEPTION("No feature keys found in model");
  }
  bool is_gga = false;
  bool is_mgga = false;
  for (const auto& key : feature_keys) {
    if ( not valueExists(key) ) GAUXC_GENERIC_EXCEPTION("Feature Key Required Not Implemented: " + key);
    if (key == feat_map.at(ONEDFT_FEATURE::TAU)) is_mgga = true;
    if (key == feat_map.at(ONEDFT_FEATURE::DDEN)) is_gga = true;
  }
  if (is_mgga) is_gga = false;

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    pre_onedft_local_work_( basis, Ps, ldps, Pz, ldpz, &N_EL, is_gga, is_mgga, false /*needs_laplacian*/);
  });
  std::vector<int> sendcounts(rt.comm_size(), 0);
  std::vector<int> displs(rt.comm_size(), 0);
  std::vector<int64_t> atom_reorder_inv_perm;
  FeatureDict features_dict = prepare_onedft_features(2/*ndm*/, tasks, this->load_balancer_->molecule(), feature_keys, rt, 
    sendcounts, displs, atom_reorder_inv_perm);
  if (world_rank == 0) {
    auto exc_on_grid = get_exc(exc_func, features_dict);
    // check is_nan
    if (exc_on_grid.isnan().any().item<bool>()) {
      GAUXC_GENERIC_EXCEPTION("exc_on_grid has NaN");
    }
    auto exc = (exc_on_grid * features_dict.at(feat_map.at(ONEDFT_FEATURE::WEIGHTS))).sum();
    exc.backward();
    EXC[0] = exc.item().to<double>();
  } else {
    EXC[0] = 0.0;
  }
  // TODO: stop here if only exc

  send_buffer_onedft_outputs(2/*ndm*/, features_dict, tasks, rt, sendcounts, displs, atom_reorder_inv_perm);

  this->timer_.time_op("XCIntegrator.LocalWork2", [&](){
    post_onedft_local_work_( basis, Ps, ldps, Pz, ldpz, VXCs, n, VXCz, n, is_gga, is_mgga, false /*needs_laplacian*/);
  });

  this->timer_.time_op("XCIntegrator.Allreduce", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");

    this->reduction_driver_->allreduce_inplace( VXCs, nbf*nbf, ReductionOp::Sum );
    if(VXCz) this->reduction_driver_->allreduce_inplace( VXCz, nbf*nbf, ReductionOp::Sum );

    this->reduction_driver_->allreduce_inplace( EXC,   1    , ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_EL, 1    , ReductionOp::Sum );

  });

}

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  pre_onedft_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
    const value_type* Pz, int64_t ldpz, value_type *N_EL, 
    const bool is_gga, const bool is_mgga, const bool needs_laplacian) {

  const bool is_uks = (Pz != nullptr);
  const bool is_rks = not is_uks;
  const bool is_lda = not is_gga and not is_mgga;
  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());
  const auto& mol   = this->load_balancer_->molecule();
  BasisSetMap basis_map(basis,mol);
  const int32_t nbf = basis.nbf();

  auto& tasks = this->load_balancer_->get_tasks();
  const size_t ntasks = tasks.size();

  // Check that Partition Weights have been calculated
  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored ) {
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Been Modified");
  }

  double NEL_WORK = 0.0;

  #pragma omp parallel
  {

  XCHostData<value_type> host_data; // Thread local host data
  
  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];

    // Get tasks constants
    const int32_t  npts    = task.points.size();
    const int32_t  nbe     = task.bfn_screening.nbe;
    const int32_t  nshells = task.bfn_screening.shell_list.size();

    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();
    const int32_t* shell_list = task.bfn_screening.shell_list.data();
    const size_t spin_dim_scal = is_rks ? 1 : 2; 
    const size_t mgga_dim_scal = is_mgga ? 4 : 1; // basis + d1basis
 
    // Things that every calc needs
    host_data.nbe_scr .resize(nbe  * nbe);
    host_data.zmat    .resize(npts * nbe * spin_dim_scal * mgga_dim_scal); 

    // LDA data requirements
    if( is_lda ){
      host_data.basis_eval .resize( npts * nbe );
      task.feat.den_eval.resize(npts * spin_dim_scal);
    }
     
    // GGA data requirements
    const size_t gga_dim_scal = is_rks ? 1 : 3;
    if( is_gga ){
      host_data.basis_eval .resize( 4 * npts * nbe );
      host_data.gamma      .resize( gga_dim_scal * npts ); // TODO: delete gamma
      task.feat.den_eval.resize(npts * spin_dim_scal);
      task.feat.dden_x_eval.resize(npts * spin_dim_scal);
      task.feat.dden_y_eval.resize(npts * spin_dim_scal);
      task.feat.dden_z_eval.resize(npts * spin_dim_scal);
    }
    if( is_mgga ){
      host_data.basis_eval .resize( 4 * npts * nbe ); // basis + grad (3)
      host_data.gamma      .resize( gga_dim_scal * npts );
      task.feat.den_eval.resize(npts * spin_dim_scal);
      task.feat.dden_x_eval.resize(npts * spin_dim_scal);
      task.feat.dden_y_eval.resize(npts * spin_dim_scal);
      task.feat.dden_z_eval.resize(npts * spin_dim_scal);
      task.feat.tau.resize(npts * spin_dim_scal);
    }

    // Alias/Partition out scratch memory
    auto* basis_eval = host_data.basis_eval.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* zmat       = host_data.zmat.data();
    auto* gamma      = host_data.gamma.data();
    auto* lapl       = host_data.lapl.data();

    decltype(zmat) zmat_z = nullptr;
    decltype(zmat) zmat_x = nullptr;
    decltype(zmat) zmat_y = nullptr;
    if(!is_rks) {
      zmat_z = zmat + mgga_dim_scal * nbe * npts;
    }
    
    auto* den_eval   = task.feat.den_eval.data();
    auto* tau        = task.feat.tau.data();
    auto* dden_x_eval = task.feat.dden_x_eval.data();
    auto* dden_y_eval = task.feat.dden_y_eval.data();
    auto* dden_z_eval = task.feat.dden_z_eval.data();

    value_type* dbasis_x_eval = nullptr;
    value_type* dbasis_y_eval = nullptr;
    value_type* dbasis_z_eval = nullptr;
    value_type* lbasis_eval = nullptr;

    value_type* mmat_x      = nullptr;
    value_type* mmat_y      = nullptr;
    value_type* mmat_z      = nullptr;
    value_type* mmat_x_z    = nullptr;
    value_type* mmat_y_z    = nullptr;
    value_type* mmat_z_z    = nullptr;

    if( is_gga ) {
      dbasis_x_eval = basis_eval    + npts * nbe;
      dbasis_y_eval = dbasis_x_eval + npts * nbe;
      dbasis_z_eval = dbasis_y_eval + npts * nbe;
    }

    if ( is_mgga ) {
      dbasis_x_eval = basis_eval    + npts * nbe;
      dbasis_y_eval = dbasis_x_eval + npts * nbe;
      dbasis_z_eval = dbasis_y_eval + npts * nbe;
      mmat_x        = zmat + npts * nbe;
      mmat_y        = mmat_x + npts * nbe;
      mmat_z        = mmat_y + npts * nbe;
      if(is_uks) {
        mmat_x_z = zmat_z + npts * nbe;
        mmat_y_z = mmat_x_z + npts * nbe;
        mmat_z_z = mmat_y_z + npts * nbe;
      }
    }


    // Get the submatrix map for batch
    std::vector< std::array<int32_t, 3> > submat_map;
    std::tie(submat_map, std::ignore) =
          gen_compressed_submat_map(basis_map, task.bfn_screening.shell_list, nbf, nbf);

    // Evaluate Collocation
    if( is_mgga ) {
        lwd->eval_collocation_gradient( npts, nshells, nbe, points, basis, shell_list,
          basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );
    } else if( is_gga )
      lwd->eval_collocation_gradient( npts, nshells, nbe, points, basis, shell_list,
        basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );
    else
      lwd->eval_collocation( npts, nshells, nbe, points, basis, shell_list,
        basis_eval );

     
    // Evaluate X matrix (fac * P * B) -> store in Z
    const auto xmat_fac = is_rks ? 2.0 : 1.0; // TODO Fix for spinor RKS input
    lwd->eval_xmat( mgga_dim_scal * npts, nbf, nbe, submat_map, xmat_fac, Ps, ldps, basis_eval, nbe,
      zmat, nbe, nbe_scr );

    // X matrix for Pz
    if(not is_rks) {
      lwd->eval_xmat( mgga_dim_scal * npts, nbf, nbe, submat_map, 1.0, Pz, ldpz, basis_eval, nbe,
        zmat_z, nbe, nbe_scr);
    }

    // Evaluate U and V variables
    if( is_mgga ) {
      if (is_rks) {
        lwd->eval_uvvar_mgga_rks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
          dbasis_z_eval, lbasis_eval, zmat, nbe, mmat_x, mmat_y, mmat_z, 
          nbe, den_eval, dden_x_eval, dden_y_eval, dden_z_eval, gamma, tau, lapl);
      } else if (is_uks) {
        lwd->eval_uvvar_mgga_uks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
          dbasis_z_eval, lbasis_eval, zmat, nbe, zmat_z, nbe, 
          mmat_x, mmat_y, mmat_z, nbe, mmat_x_z, mmat_y_z, mmat_z_z, nbe, 
          den_eval, dden_x_eval, dden_y_eval, dden_z_eval, gamma, tau, lapl);
      }
    } else if ( is_gga ) {
      if(is_rks) {
        lwd->eval_uvvar_gga_rks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
          dbasis_z_eval, zmat, nbe, den_eval, dden_x_eval, dden_y_eval, dden_z_eval,
          gamma );
      } else if(is_uks) {
        lwd->eval_uvvar_gga_uks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
          dbasis_z_eval, zmat, nbe, zmat_z, nbe, den_eval, dden_x_eval, 
          dden_y_eval, dden_z_eval, gamma );
      }
     } else {
      if(is_rks) {
        lwd->eval_uvvar_lda_rks( npts, nbe, basis_eval, zmat, nbe, den_eval );
      } else if(is_uks) {
        lwd->eval_uvvar_lda_uks( npts, nbe, basis_eval, zmat, nbe, zmat_z, nbe,
          den_eval );
      }
    }
    
    // Scalar integrations
    double NEL_local = 0.0;
    for( int32_t i = 0; i < npts; ++i ) {
      const auto den = is_rks ? den_eval[i] : (den_eval[2*i] + den_eval[2*i+1]);
      NEL_local += weights[i] * den;
    }

    // Atomic updates
    #pragma omp atomic
    NEL_WORK += NEL_local;
  } // Loop over tasks
}  // End OpenMP region
*N_EL = NEL_WORK;
// std::cout << "N_EL: " << *N_EL << std::endl;
}

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  post_onedft_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
    const value_type* Pz, int64_t ldpz,
    value_type* VXCs, int64_t ldvxcs,
    value_type* VXCz, int64_t ldvxcz,
    const bool is_gga, const bool is_mgga, const bool needs_laplacian) {


    const bool is_uks = (Pz != nullptr);
    const bool is_rks = not is_uks;
    const bool is_lda = not is_gga and not is_mgga;
    auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

    const auto& mol   = this->load_balancer_->molecule();
    BasisSetMap basis_map(basis,mol);
    const int32_t nbf = basis.nbf();

    // Zero out integrands
    
    if(VXCs)
    for( auto j = 0; j < nbf; ++j ) {
      for( auto i = 0; i < nbf; ++i ) {
        VXCs[i + j*ldvxcs] = 0.;
      }
    }

    if(VXCz) {
      for( auto j = 0; j < nbf; ++j ) {
        for( auto i = 0; i < nbf; ++i ) {
          VXCz[i + j*ldvxcz] = 0.;
        }
      }
    }

  // Loop over tasks
  auto& tasks = this->load_balancer_->get_tasks();
  const size_t ntasks = tasks.size();

  #pragma omp parallel
  {

  XCHostData<value_type> host_data; // Thread local host data

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    const auto& task = tasks[iT];

    // Get tasks constants
    const int32_t  npts    = task.points.size();
    const int32_t  nbe     = task.bfn_screening.nbe;
    const int32_t  nshells = task.bfn_screening.shell_list.size();

    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();
    const int32_t* shell_list = task.bfn_screening.shell_list.data();
    const size_t spin_dim_scal = is_rks ? 1 : is_uks ? 2 : 4; // last case is_gks
    const size_t mgga_dim_scal = is_mgga ? 4 : 1; // basis + d1basis
 
    // Things that every calc needs
    host_data.nbe_scr .resize(nbe  * nbe);
    host_data.zmat    .resize(npts * nbe * spin_dim_scal * mgga_dim_scal); 

    // LDA data requirements
    if( is_lda ){
      host_data.basis_eval .resize( npts * nbe );
    }
    // GGA data requirements
    const size_t gga_dim_scal = is_rks ? 1 : 3;
    if( is_gga ){
      host_data.basis_eval .resize( 4 * npts * nbe );
    }
    if( is_mgga ){
      host_data.basis_eval .resize( 4 * npts * nbe ); // basis + grad (3)
    }

    // Alias/Partition out scratch memory
    auto* basis_eval = host_data.basis_eval.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* zmat       = host_data.zmat.data();
    auto* lapl       = host_data.lapl.data();

    decltype(zmat) zmat_z = nullptr;
    decltype(zmat) zmat_x = nullptr;
    decltype(zmat) zmat_y = nullptr;
    if(!is_rks) {
      zmat_z = zmat + mgga_dim_scal * nbe * npts;
    }

    value_type* dbasis_x_eval = nullptr;
    value_type* dbasis_y_eval = nullptr;
    value_type* dbasis_z_eval = nullptr;
    value_type* lbasis_eval = nullptr;
    value_type* mmat_x      = nullptr;
    value_type* mmat_y      = nullptr;
    value_type* mmat_z      = nullptr;
    value_type* mmat_x_z    = nullptr;
    value_type* mmat_y_z    = nullptr;
    value_type* mmat_z_z    = nullptr;

    if( is_gga ) {
      dbasis_x_eval = basis_eval    + npts * nbe;
      dbasis_y_eval = dbasis_x_eval + npts * nbe;
      dbasis_z_eval = dbasis_y_eval + npts * nbe;
    }

    if ( is_mgga ) {
      dbasis_x_eval = basis_eval    + npts * nbe;
      dbasis_y_eval = dbasis_x_eval + npts * nbe;
      dbasis_z_eval = dbasis_y_eval + npts * nbe;
      mmat_x        = zmat + npts * nbe;
      mmat_y        = mmat_x + npts * nbe;
      mmat_z        = mmat_y + npts * nbe;
      if(is_uks) {
        mmat_x_z = zmat_z + npts * nbe;
        mmat_y_z = mmat_x_z + npts * nbe;
        mmat_z_z = mmat_y_z + npts * nbe;
      }
    }

    // assume always uks
    const value_type* vdden_eval_a, *vdden_eval_b;
    const value_type* vdden_x_eval_a, *vdden_y_eval_a, *vdden_z_eval_a;
    const value_type* vdden_x_eval_b, *vdden_y_eval_b, *vdden_z_eval_b;
    const value_type* vtau;
    std::vector<value_type> vrho;

    vdden_eval_a = task.feat.vdden_eval_a.data();
    vdden_eval_b = task.feat.vdden_eval_b.data();
    if (is_gga || is_mgga) {
      vdden_x_eval_a = task.feat.vdden_x_eval_a.data();
      vdden_y_eval_a = task.feat.vdden_y_eval_a.data();
      vdden_z_eval_a = task.feat.vdden_z_eval_a.data();
      vdden_x_eval_b = task.feat.vdden_x_eval_b.data();
      vdden_y_eval_b = task.feat.vdden_y_eval_b.data();
      vdden_z_eval_b = task.feat.vdden_z_eval_b.data();
    } else { // lda
      vrho.resize(npts * spin_dim_scal);
      interleave_data(task.feat.vdden_eval_a.data(), task.feat.vdden_eval_b.data(), npts, vrho.data());
    }
    if (is_mgga) {
      vtau = task.feat.vtau.data();
    }
    // Get the submatrix map for batch
    std::vector< std::array<int32_t, 3> > submat_map;
    std::tie(submat_map, std::ignore) =
          gen_compressed_submat_map(basis_map, task.bfn_screening.shell_list, nbf, nbf);

    // Evaluate Collocation (+ Grad and Hessian)
    if( is_mgga ) {
      lwd->eval_collocation_gradient( npts, nshells, nbe, points, basis, shell_list,
        basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );
    } else if( is_gga )
      lwd->eval_collocation_gradient( npts, nshells, nbe, points, basis, shell_list,
        basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );
    else
      lwd->eval_collocation( npts, nshells, nbe, points, basis, shell_list,
        basis_eval );

    // Evaluate Z matrix for VXC
    if (is_gga){
        eval_zmat_gga_vxc_uks( npts, nbe, vdden_eval_a, vdden_eval_b, vdden_x_eval_a, vdden_x_eval_b, vdden_y_eval_a, 
                                vdden_y_eval_b, vdden_z_eval_a, vdden_z_eval_b, basis_eval,
                                dbasis_x_eval, dbasis_y_eval, dbasis_z_eval, zmat, nbe, zmat_z, nbe);
    } else if (is_mgga) {
        eval_zmat_mgga_vxc_uks( npts, nbe, vdden_eval_a, vdden_eval_b, 
                              nullptr, nullptr, /* vlapl_a, vlapl_b */
                              vdden_x_eval_a, vdden_x_eval_b, vdden_y_eval_a, 
                              vdden_y_eval_b, vdden_z_eval_a, vdden_z_eval_b, basis_eval,
                              dbasis_x_eval, dbasis_y_eval, dbasis_z_eval, 
                              lbasis_eval, zmat, nbe, zmat_z, nbe);
        lwd->eval_mmat_mgga_vxc_uks( npts, nbe, vtau, nullptr /*vlapl*/, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval,
                                mmat_x, mmat_y, mmat_z, nbe, mmat_x_z, mmat_y_z, mmat_z_z, nbe);
    } else {
      lwd->eval_zmat_lda_vxc_uks( npts, nbe, vrho.data(), basis_eval, zmat, nbe, zmat_z, nbe );
    }

    // Increment VXC
    lwd->inc_vxc( mgga_dim_scal * npts, nbf, nbe, basis_eval, submat_map, zmat, nbe, VXCs, ldvxcs, nbe_scr );
    if(is_uks) {
      lwd->inc_vxc( mgga_dim_scal * npts, nbf, nbe, basis_eval, submat_map, zmat_z, nbe,VXCz, ldvxcz, nbe_scr);
    }
  } // loop over tasks
  } // end OpenMP region

  for( int32_t j = 0;   j < nbf; ++j ) {
    for( int32_t i = j+1; i < nbf; ++i ) {
      VXCs[ j + i*ldvxcs ] = VXCs[ i + j*ldvxcs ];
    }
  }
  if(not is_rks) {
    for( int32_t j = 0;   j < nbf; ++j ) {
      for( int32_t i = j+1; i < nbf; ++i ) {
        VXCz[ j + i*ldvxcz ] = VXCz[ i + j*ldvxcz ];
      }
    }
  }
}

void eval_zmat_gga_vxc_uks(size_t npts, size_t nbf, 
  const double* vdden_eval_a, const double* vdden_eval_b, const double* vdden_x_eval_a, const double* vdden_x_eval_b,
  const double* vdden_y_eval_a, const double* vdden_y_eval_b, const double* vdden_z_eval_a, const double* vdden_z_eval_b,
  const double* basis_eval, const double* dbasis_x_eval, const double* dbasis_y_eval, const double* dbasis_z_eval,
  double* Zs, size_t ldzs, double* Zz, size_t ldzz) {
  if( ldzs != nbf ) GAUXC_GENERIC_EXCEPTION(std::string("Invalid Dims"));
  if( ldzz != nbf ) GAUXC_GENERIC_EXCEPTION(std::string("Invalid Dims"));
  blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Zs, ldzs);
  blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Zz, ldzz);

  for( int32_t i = 0; i < (int32_t)npts; ++i ) {

    const int32_t ioff = i * nbf;

    auto* zs_col = Zs + ioff;
    auto* zz_col = Zz + ioff;
    auto* bf_x_col = dbasis_x_eval + ioff;
    auto* bf_y_col = dbasis_y_eval + ioff;
    auto* bf_z_col = dbasis_z_eval + ioff;

    const double factp = 0.5 * vdden_eval_a[i];
    const double factm = 0.5 * vdden_eval_b[i];

    GauXC::blas::scal( nbf, 0.5*(factp + factm), zs_col, 1 );
    GauXC::blas::scal( nbf, 0.5*(factp - factm), zz_col, 1 );

    const double x_factp = 0.5*(vdden_x_eval_a[i] + vdden_x_eval_b[i]);
    const double y_factp = 0.5*(vdden_y_eval_a[i] + vdden_y_eval_b[i]);
    const double z_factp = 0.5*(vdden_z_eval_a[i] + vdden_z_eval_b[i]);
    const double x_factm = 0.5*(vdden_x_eval_a[i] - vdden_x_eval_b[i]);
    const double y_factm = 0.5*(vdden_y_eval_a[i] - vdden_y_eval_b[i]);
    const double z_factm = 0.5*(vdden_z_eval_a[i] - vdden_z_eval_b[i]);
    GauXC::blas::axpy( nbf, x_factp, bf_x_col, 1, zs_col, 1 );
    GauXC::blas::axpy( nbf, y_factp, bf_y_col, 1, zs_col, 1 );
    GauXC::blas::axpy( nbf, z_factp, bf_z_col, 1, zs_col, 1 );

    GauXC::blas::axpy( nbf, x_factm, bf_x_col, 1, zz_col, 1 );
    GauXC::blas::axpy( nbf, y_factm, bf_y_col, 1, zz_col, 1 );
    GauXC::blas::axpy( nbf, z_factm, bf_z_col, 1, zz_col, 1 );
  }
}

void eval_zmat_mgga_vxc_uks(size_t npts, size_t nbf, 
  const double* vdden_eval_a, const double* vdden_eval_b, 
  const double* vlapl_a, const double* vlapl_b,
  const double* vdden_x_eval_a, const double* vdden_x_eval_b, const double* vdden_y_eval_a, 
  const double* vdden_y_eval_b, const double* vdden_z_eval_a, const double* vdden_z_eval_b,
  const double* basis_eval, const double* dbasis_x_eval, const double* dbasis_y_eval, const double* dbasis_z_eval,
  const double* lbasis_eval, double* Zs, size_t ldzs, double* Zz, size_t ldzz){

  if( ldzs != nbf ) GAUXC_GENERIC_EXCEPTION(std::string("Invalid Dims"));
  if( ldzz != nbf ) GAUXC_GENERIC_EXCEPTION(std::string("Invalid Dims"));
  blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Zs, ldzs);
  blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Zz, ldzz);

  for( int32_t i = 0; i < (int32_t)npts; ++i ) {

    const int32_t ioff = i * nbf;

    auto* zs_col = Zs + ioff;
    auto* zz_col = Zz + ioff;
    auto* bf_x_col = dbasis_x_eval + ioff;
    auto* bf_y_col = dbasis_y_eval + ioff;
    auto* bf_z_col = dbasis_z_eval + ioff;

    const double factp = 0.5 * vdden_eval_a[i];
    const double factm = 0.5 * vdden_eval_b[i];

    GauXC::blas::scal( nbf, 0.5*(factp + factm), zs_col, 1 ); 
    GauXC::blas::scal( nbf, 0.5*(factp - factm), zz_col, 1 );

    const double x_factp = 0.5*(vdden_x_eval_a[i] + vdden_x_eval_b[i]);
    const double y_factp = 0.5*(vdden_y_eval_a[i] + vdden_y_eval_b[i]);
    const double z_factp = 0.5*(vdden_z_eval_a[i] + vdden_z_eval_b[i]);
    const double x_factm = 0.5*(vdden_x_eval_a[i] - vdden_x_eval_b[i]);
    const double y_factm = 0.5*(vdden_y_eval_a[i] - vdden_y_eval_b[i]);
    const double z_factm = 0.5*(vdden_z_eval_a[i] - vdden_z_eval_b[i]);
    
    GauXC::blas::axpy( nbf, x_factp, bf_x_col, 1, zs_col, 1 );
    GauXC::blas::axpy( nbf, y_factp, bf_y_col, 1, zs_col, 1 );
    GauXC::blas::axpy( nbf, z_factp, bf_z_col, 1, zs_col, 1 );

    GauXC::blas::axpy( nbf, x_factm, bf_x_col, 1, zz_col, 1 );
    GauXC::blas::axpy( nbf, y_factm, bf_y_col, 1, zz_col, 1 );
    GauXC::blas::axpy( nbf, z_factm, bf_z_col, 1, zz_col, 1 );

    if (vlapl_a != nullptr) {
      auto* lbf_col = lbasis_eval + ioff;
      const auto lfactp = vlapl_a[i];
      const auto lfactm = vlapl_b[i];
      blas::axpy( nbf, 0.5*(lfactp + lfactm), lbf_col, 1, zs_col, 1);
      blas::axpy( nbf, 0.5*(lfactp - lfactm), lbf_col, 1, zz_col, 1);
    }
  }
}

void sz_to_ab(std::vector<double>& dden){
  for (size_t i = 0; i < dden.size()/2; i++) {
    double s = dden[2*i];
    double z = dden[2*i+1];
    dden[2*i] = 0.5 * (s + z);
    dden[2*i+1] = 0.5 * (s - z);
  }
}

void interleave_data(const double* a, const double* b, const size_t n, double* result) {
  for (size_t i = 0; i < n; ++i) {
    result[2*i] = a[i];
    result[2*i+1] = b[i];
  }
}

FeatureDict prepare_onedft_features(const int ndm, std::vector<XCTask>& tasks, const Molecule& mol, 
                  const std::vector<std::string> feature_keys, const RuntimeEnvironment& rt,
                  std::vector<int>& sendcounts, std::vector<int>& displs,
                  std::vector<int64_t>& atom_reorder_inv_perm) {
  std::vector<double> den_eval, dden_eval, tau, grid_coords, grid_weights, raw_grid_weights;
  // Sort tasks by atom index so that grid points are grouped by atom.
  // build_atom_reorder_perm assumes this contiguous-by-atom layout.
  std::stable_sort(tasks.begin(), tasks.end(),
    [](const auto& a, const auto& b) { return a.iParent < b.iParent; });
  int total_npts = std::accumulate( tasks.begin(), tasks.end(), 0,
    [](const auto& a, const auto& b) { return a + b.npts; } );
  grid_coords.reserve(total_npts * 3);
  grid_weights.reserve(total_npts);
  raw_grid_weights.reserve(total_npts);
  den_eval.reserve(total_npts * ndm);
  dden_eval.resize(total_npts * 6);  // 2 values per point, 3 components
  tau.reserve(total_npts * ndm);

  int offset = 0;
  for (auto& task : tasks) {
    for (const auto& point : task.points) {
      grid_coords.push_back(point[0]);
      grid_coords.push_back(point[1]);
      grid_coords.push_back(point[2]);
    }
    std::copy(task.weights.begin(), task.weights.end(), std::back_inserter(grid_weights));
    std::copy(task.raw_weights.begin(), task.raw_weights.end(), std::back_inserter(raw_grid_weights));
    std::copy(task.feat.den_eval.begin(), task.feat.den_eval.end(), std::back_inserter(den_eval));

    if (task.feat.dden_x_eval.size() != 0){
      sz_to_ab(task.feat.dden_x_eval);
      sz_to_ab(task.feat.dden_y_eval);
      sz_to_ab(task.feat.dden_z_eval);
      for (size_t i = 0; i < task.points.size(); i++) {
        dden_eval[6 * i + 0 + 6 * offset] = task.feat.dden_x_eval[2 * i];
        dden_eval[6 * i + 1 + 6 * offset] = task.feat.dden_y_eval[2 * i];
        dden_eval[6 * i + 2 + 6 * offset] = task.feat.dden_z_eval[2 * i];
        dden_eval[6 * i + 3 + 6 * offset] = task.feat.dden_x_eval[2 * i + 1];  
        dden_eval[6 * i + 4 + 6 * offset] = task.feat.dden_y_eval[2 * i + 1];
        dden_eval[6 * i + 5 + 6 * offset] = task.feat.dden_z_eval[2 * i + 1];
      }
    }
    offset += task.points.size();
    std::copy(task.feat.tau.begin(), task.feat.tau.end(), std::back_inserter(tau));
  }

  // Compute per-atom grid sizes from local tasks
  int natoms = mol.size();
  std::vector<int64_t> atomic_grid_sizes_vec(natoms, 0);
  for (const auto& task : tasks) {
    if (task.iParent >= 0 && task.iParent < natoms) {
      atomic_grid_sizes_vec[task.iParent] += task.npts;
    }
  }

  // MPI gather all data to rank 0 and reorder from rank-order to atom-order
  int world_rank = rt.comm_rank();
  int local_npts = total_npts; // save before gather overwrites
  auto reorder_result = mpi_gather_and_reorder(
    den_eval, dden_eval, tau, grid_coords, grid_weights,
    atomic_grid_sizes_vec, total_npts, natoms, rt, sendcounts, displs);
  total_npts = reorder_result.total_npts;
  atom_reorder_inv_perm = std::move(reorder_result.inv_perm);
  auto& global_atomic_grid_sizes_vec = reorder_result.global_atomic_grid_sizes;

  // Gather and reorder raw_grid_weights using the same MPI layout
  GAUXC_MPI_CODE(
    if (rt.comm_size() > 1) {
      std::vector<double> recv_raw(world_rank == 0 ? total_npts : 0);
      MPI_Gatherv(raw_grid_weights.data(), local_npts, MPI_DOUBLE,
                  recv_raw.data(), sendcounts.data(), displs.data(),
                  MPI_DOUBLE, 0, rt.comm());
      if (world_rank == 0) {
        raw_grid_weights = std::move(recv_raw);
        // Reconstruct forward perm from inv_perm and reorder
        std::vector<int64_t> perm(total_npts);
        for (int64_t j = 0; j < total_npts; j++) perm[atom_reorder_inv_perm[j]] = j;
        std::vector<double> tmp(total_npts);
        for (int64_t i = 0; i < total_npts; i++) tmp[perm[i]] = raw_grid_weights[i];
        raw_grid_weights = std::move(tmp);
      }
    }
  )

  FeatureDict featmap;
  if (world_rank == 0) {
    int64_t max_grid_size = *std::max_element(
      global_atomic_grid_sizes_vec.begin(), global_atomic_grid_sizes_vec.end());
    std::vector<double> coarse_0_atomic_coords (natoms*3); 
    for (int i = 0; i < natoms; i++) {
      coarse_0_atomic_coords[3*i] = mol[i].x;
      coarse_0_atomic_coords[3*i+1] = mol[i].y;
      coarse_0_atomic_coords[3*i+2] = mol[i].z;
    }

    auto options = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);
    for (const auto& key : feature_keys) {
      auto enum_key = reverse_feat_map.at(key);
      at::Tensor tensor;
      switch (enum_key) {
      case ONEDFT_FEATURE::DEN: {
        auto flat_tensor = torch::from_blob(den_eval.data(), {ndm, total_npts}, {1, ndm}, options);
        tensor = flat_tensor.clone().requires_grad_(true);
        break;
      }
      case ONEDFT_FEATURE::DDEN: {
        auto flat_tensor = torch::from_blob(dden_eval.data(), {ndm, 3, total_npts}, {3, 1, 3*ndm}, options);
        tensor = flat_tensor.clone().requires_grad_(true);
        break;
      }
      case ONEDFT_FEATURE::TAU: {
        auto flat_tensor = torch::from_blob(tau.data(), {ndm, total_npts}, {1, ndm}, options);
        tensor = flat_tensor.clone().requires_grad_(true);
        break;
      }
      case ONEDFT_FEATURE::POINTS: {
        auto flat_tensor = torch::from_blob(grid_coords.data(), {total_npts, 3}, options);
        tensor = flat_tensor.clone();
        break;
      }
      case ONEDFT_FEATURE::WEIGHTS: {
        auto flat_tensor = torch::from_blob(grid_weights.data(), {total_npts}, options);
        tensor = flat_tensor.clone();
        break;
      }
      case ONEDFT_FEATURE::COORDS: {
        auto flat_tensor = torch::from_blob(coarse_0_atomic_coords.data(), {natoms, 3}, options);
        tensor = flat_tensor.clone();
        break;
      }
      case ONEDFT_FEATURE::ATOMIC_GRID_WEIGHTS: {
        auto flat_tensor = torch::from_blob(raw_grid_weights.data(), {total_npts}, options);
        tensor = flat_tensor.clone();
        break;
      }
      case ONEDFT_FEATURE::ATOMIC_GRID_SIZES: {
        auto sizes_options = torch::TensorOptions().dtype(torch::kInt64).device(torch::kCPU);
        tensor = torch::from_blob(global_atomic_grid_sizes_vec.data(), {natoms}, sizes_options).clone();
        break;
      }
      case ONEDFT_FEATURE::ATOMIC_GRID_SIZE_BOUND_SHAPE: {
        auto sizes_options = torch::TensorOptions().dtype(torch::kInt64).device(torch::kCPU);
        tensor = torch::zeros({max_grid_size, 0}, sizes_options);
        break;
      }
      default:
        GAUXC_GENERIC_EXCEPTION("Feature Key Not Implemented: " + key);
      }
      if (tensor.isnan().any().item<bool>()) {
        GAUXC_GENERIC_EXCEPTION("NaN detected in feature tensor: " + key);
      }
      featmap.insert(key, tensor);
    }
  }
  return featmap;
}

void send_buffer_onedft_outputs(const int ndm, const FeatureDict features_dict, std::vector<XCTask>& tasks, 
                                const RuntimeEnvironment& rt, std::vector<int> sendcounts, std::vector<int> displs,
                                const std::vector<int64_t>& atom_reorder_inv_perm) {

  std::vector<double> den_eval, dden_eval, tau;
  auto total_npts = mpi_scatter_onedft_outputs(features_dict, rt.comm_rank(), rt.comm_size(),
                                                sendcounts, displs, atom_reorder_inv_perm,
                                                den_eval, dden_eval, tau);

  size_t offset = 0;
  for (auto&task : tasks) {
    int64_t npts = task.points.size();
    task.feat.vdden_eval_a.resize(npts);
    task.feat.vdden_eval_b.resize(npts);
    auto den_a_slice = den_eval.data() + offset;
    auto den_b_slice = den_eval.data() + total_npts + offset;
    std::copy(den_a_slice, den_a_slice + npts, task.feat.vdden_eval_a.begin());
    std::copy(den_b_slice, den_b_slice + npts, task.feat.vdden_eval_b.begin());
    if (task.feat.dden_x_eval.size() != 0){
      task.feat.vdden_x_eval_a.resize(npts);
      task.feat.vdden_y_eval_a.resize(npts);
      task.feat.vdden_z_eval_a.resize(npts);
      task.feat.vdden_x_eval_b.resize(npts);
      task.feat.vdden_y_eval_b.resize(npts);
      task.feat.vdden_z_eval_b.resize(npts);

      auto dden_a_x_slice = dden_eval.data() + offset;
      auto dden_a_y_slice = dden_eval.data() + total_npts + offset;
      auto dden_a_z_slice = dden_eval.data() + total_npts * 2 + offset;
      auto dden_b_x_slice = dden_eval.data() + total_npts * 3 + offset;
      auto dden_b_y_slice = dden_eval.data() + total_npts * 4 + offset;
      auto dden_b_z_slice = dden_eval.data() + total_npts * 5 + offset;

      std::copy(dden_a_x_slice, dden_a_x_slice + npts, task.feat.vdden_x_eval_a.begin());
      std::copy(dden_a_y_slice, dden_a_y_slice + npts, task.feat.vdden_y_eval_a.begin());
      std::copy(dden_a_z_slice, dden_a_z_slice + npts, task.feat.vdden_z_eval_a.begin());
      std::copy(dden_b_x_slice, dden_b_x_slice + npts, task.feat.vdden_x_eval_b.begin());
      std::copy(dden_b_y_slice, dden_b_y_slice + npts, task.feat.vdden_y_eval_b.begin());
      std::copy(dden_b_z_slice, dden_b_z_slice + npts, task.feat.vdden_z_eval_b.begin());
    }

    if (task.feat.tau.size() != 0){
      task.feat.vtau.resize(npts * 2);
      auto vtau_a = tau.data() + offset;
      auto vtau_b = tau.data() + total_npts + offset;
      interleave_data(vtau_a, vtau_b, npts, task.feat.vtau.data());
    }
    offset += npts;
  }
  if (offset != total_npts) {
    GAUXC_GENERIC_EXCEPTION("Mismatch in number of points for onedft features.");
  }
}

// RKS OneDFT driver - delegates to generic GKS impl
// template <typename ValueType>
// void ReferenceReplicatedXCHostIntegrator<ValueType>::
//   eval_exc_vxc_onedft_( int64_t m, int64_t n, 
//                  const value_type* P, int64_t ldp,
//                  value_type* VXC, int64_t ldvxc,
//                  value_type* EXC, const IntegratorSettingsXC& ks_settings) {

//   eval_exc_vxc_onedft_(m, n, P, ldp, nullptr, 0, nullptr, 0, nullptr, 0,
//     VXC, ldvxc, nullptr, 0, nullptr, 0, nullptr, 0, EXC, ks_settings);

// }


// ============================================================================
// OneDFT EXC Gradient
// ============================================================================

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_exc_grad_onedft_( int64_t m, int64_t n, const value_type* Ps, int64_t ldps,
                         const value_type* Pz, int64_t ldpz,
                         value_type* EXC_GRAD, const IntegratorSettingsXC& settings ) {

  const auto& basis = this->load_balancer_->basis();
  const int64_t nbf = basis.nbf();
  if( m != n )    GAUXC_GENERIC_EXCEPTION("P Must Be Square");
  if( m != nbf )  GAUXC_GENERIC_EXCEPTION("P Must Have Same Dimension as Basis");
  if( ldps < nbf) GAUXC_GENERIC_EXCEPTION("Invalid LDPS");
  if( ldpz && ldpz < nbf ) GAUXC_GENERIC_EXCEPTION("Invalid LDPZ");

  const bool is_uks = (Pz != nullptr);
  if (not is_uks) {
    GAUXC_GENERIC_EXCEPTION("RKS OneDFT gradient Not Yet Implemented");
  }

  // Get Tasks
  auto& tasks = this->load_balancer_->get_tasks();
  auto rt = this->load_balancer_->runtime();
  int32_t world_rank = rt.comm_rank();

  // Load model
  OneDFTSettings onedft_settings;
  if( auto* tmp = dynamic_cast<const OneDFTSettings*>(&settings) ) {
    onedft_settings = *tmp;
  }
  const auto model_path = onedft_settings.model;
  torch::DeviceType device = torch::kCPU;
  auto [exc_func, feature_keys] = load_model(model_path, device);

  // Determine feature requirements
  bool is_gga = false;
  bool is_mgga = false;
  for (const auto& key : feature_keys) {
    if ( not valueExists(key) ) GAUXC_GENERIC_EXCEPTION("Feature Key Required Not Implemented: " + key);
    if (key == feat_map.at(ONEDFT_FEATURE::TAU)) is_mgga = true;
    if (key == feat_map.at(ONEDFT_FEATURE::DDEN)) is_gga = true;
  }
  if (is_mgga) is_gga = false;

  value_type N_EL;

  // Step 1: Pre-work (basis eval, density computation)
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    pre_onedft_local_work_( basis, Ps, ldps, Pz, ldpz, &N_EL, is_gga, is_mgga, false);
  });

  // Step 2: Gather features and build torch tensors
  std::vector<int> sendcounts(rt.comm_size(), 0);
  std::vector<int> displs(rt.comm_size(), 0);
  std::vector<int64_t> atom_reorder_inv_perm;
  FeatureDict features_dict = prepare_onedft_features(2/*ndm*/, tasks,
    this->load_balancer_->molecule(), feature_keys, rt,
    sendcounts, displs, atom_reorder_inv_perm);

  // Step 3: Forward + backward with grad on points and coords
  std::vector<double> eps_on_grid_global; // exc_on_grid values for weight derivative
  std::vector<double> points_grad_global; // [total_npts * 3]
  std::vector<double> coords_grad_global; // [natoms * 3]
  const int natoms = this->load_balancer_->molecule().natoms();

  if (world_rank == 0) {
    // Enable requires_grad on points and coords tensors
    if (features_dict.find(feat_map.at(ONEDFT_FEATURE::POINTS)) != features_dict.end()) {
      features_dict.at(feat_map.at(ONEDFT_FEATURE::POINTS)).requires_grad_(true);
    }
    if (features_dict.find(feat_map.at(ONEDFT_FEATURE::COORDS)) != features_dict.end()) {
      features_dict.at(feat_map.at(ONEDFT_FEATURE::COORDS)).requires_grad_(true);
    }

    auto exc_on_grid = get_exc(exc_func, features_dict);
    if (exc_on_grid.isnan().any().item<bool>()) {
      GAUXC_GENERIC_EXCEPTION("exc_on_grid has NaN");
    }
    auto exc = (exc_on_grid * features_dict.at(feat_map.at(ONEDFT_FEATURE::WEIGHTS))).sum();
    exc.backward();

    // Extract eps_on_grid for weight derivative term
    int total_npts = exc_on_grid.size(0);
    at::Tensor eps_cpu = exc_on_grid.detach().cpu().contiguous();
    eps_on_grid_global.resize(total_npts);
    std::memcpy(eps_on_grid_global.data(), eps_cpu.data_ptr<double>(), total_npts * sizeof(double));

    // Extract points.grad() -> per-grid-point forces
    if (features_dict.find(feat_map.at(ONEDFT_FEATURE::POINTS)) != features_dict.end()) {
      auto pg = features_dict.at(feat_map.at(ONEDFT_FEATURE::POINTS)).grad();
      if (pg.defined()) {
        at::Tensor pg_cpu = pg.cpu().contiguous();
        points_grad_global.resize(total_npts * 3);
        std::memcpy(points_grad_global.data(), pg_cpu.data_ptr<double>(), total_npts * 3 * sizeof(double));
      }
    }

    // Extract coords.grad() -> per-atom forces
    if (features_dict.find(feat_map.at(ONEDFT_FEATURE::COORDS)) != features_dict.end()) {
      auto cg = features_dict.at(feat_map.at(ONEDFT_FEATURE::COORDS)).grad();
      if (cg.defined()) {
        at::Tensor cg_cpu = cg.cpu().contiguous();
        coords_grad_global.resize(natoms * 3);
        std::memcpy(coords_grad_global.data(), cg_cpu.data_ptr<double>(), natoms * 3 * sizeof(double));
      }
    }

    // Reorder eps_on_grid from atom-order back to rank-order for scatter
    if (!atom_reorder_inv_perm.empty()) {
      std::vector<double> tmp(total_npts);
      for (int64_t i = 0; i < total_npts; i++) {
        tmp[atom_reorder_inv_perm[i]] = eps_on_grid_global[i];
      }
      eps_on_grid_global = std::move(tmp);

      // Also reorder points_grad
      if (!points_grad_global.empty()) {
        std::vector<double> tmp3(total_npts * 3);
        for (int64_t i = 0; i < total_npts; i++) {
          int64_t j = atom_reorder_inv_perm[i];
          tmp3[j*3+0] = points_grad_global[i*3+0];
          tmp3[j*3+1] = points_grad_global[i*3+1];
          tmp3[j*3+2] = points_grad_global[i*3+2];
        }
        points_grad_global = std::move(tmp3);
      }
    }
  }

  // Step 4: Scatter Vxc back to tasks
  send_buffer_onedft_outputs(2/*ndm*/, features_dict, tasks, rt, sendcounts, displs, atom_reorder_inv_perm);

  // Scatter eps_on_grid to local tasks (for single rank, just distribute)
  // For MPI, would need MPI_Scatterv — for now handle single rank
  std::vector<double> eps_on_grid_local;
  if (rt.comm_size() == 1) {
    eps_on_grid_local = std::move(eps_on_grid_global);
  } else {
    // TODO: MPI scatter of eps_on_grid
    GAUXC_GENERIC_EXCEPTION("OneDFT gradient with MPI not yet implemented");
  }

  // Zero out EXC_GRAD
  for (int i = 0; i < 3*natoms; ++i) EXC_GRAD[i] = 0.0;

  // Step 6: Add autograd forces BEFORE Pulay (which re-sorts tasks!)
  // points.grad gives ∂E/∂r_g. Since grid points move with their parent atom,
  // the force on atom A = Σ_{g∈A} points_grad[g].
  // NOTE: Must be done while tasks are still in iParent-sorted order
  //       (matching points_grad_global layout). exc_grad_local_work_onedft_
  //       re-sorts tasks by workload, breaking the correspondence.
  if (!points_grad_global.empty() && world_rank == 0) {
    size_t offset = 0;
    for (const auto& task : tasks) {
      int iParent = task.iParent;
      for (size_t ipt = 0; ipt < task.points.size(); ++ipt) {
        EXC_GRAD[3*iParent + 0] += points_grad_global[(offset + ipt)*3 + 0];
        EXC_GRAD[3*iParent + 1] += points_grad_global[(offset + ipt)*3 + 1];
        EXC_GRAD[3*iParent + 2] += points_grad_global[(offset + ipt)*3 + 2];
      }
      offset += task.points.size();
    }
  }

  // coords.grad gives ∂E/∂R_A directly (no task-order dependence)
  if (!coords_grad_global.empty() && world_rank == 0) {
    for (int a = 0; a < natoms; ++a) {
      EXC_GRAD[3*a + 0] += coords_grad_global[3*a + 0];
      EXC_GRAD[3*a + 1] += coords_grad_global[3*a + 1];
      EXC_GRAD[3*a + 2] += coords_grad_global[3*a + 2];
    }
  }

  // Step 5: Pulay + weight derivative term (re-sorts tasks internally!)
  this->timer_.time_op("XCIntegrator.LocalWork2", [&](){
    exc_grad_local_work_onedft_( Ps, ldps, Pz, ldpz, EXC_GRAD, eps_on_grid_local, is_gga, is_mgga);
  });

  // Step 7: Allreduce
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){
    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");
    this->reduction_driver_->allreduce_inplace( EXC_GRAD, 3*natoms, ReductionOp::Sum );
  });
}


// Pulay + weight derivative local work using OneDFT Vxc format
template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  exc_grad_local_work_onedft_( const value_type* Ps, int64_t ldps,
                               const value_type* Pz, int64_t ldpz,
                               value_type* EXC_GRAD,
                               const std::vector<double>& eps_on_grid,
                               const bool is_gga, const bool is_mgga) {

  const bool is_uks = Pz != nullptr;
  const bool is_rks = not is_uks;

  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  const auto& basis = this->load_balancer_->basis();
  const auto& mol   = this->load_balancer_->molecule();
  const auto& molmeta = this->load_balancer_->molmeta();

  // Weight derivative settings
  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored ) {
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Been Modified");
  }
  XCWeightAlg& weight_alg = lb_state.weight_alg;

  BasisSetMap basis_map(basis, mol);
  const int32_t nbf = basis.nbf();
  const int32_t natoms = mol.natoms();

  auto& tasks = this->load_balancer_->get_tasks();
  const size_t ntasks = tasks.size();

  // Sort tasks for load balancing
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.bfn_screening.nbe) > (b.points.size() * b.bfn_screening.nbe);
  };
  std::sort( tasks.begin(), tasks.end(), task_comparator );

  // Build global eps_on_grid offset map: since tasks may be re-sorted,
  // we need to distribute eps_on_grid to tasks. We use the task ordering
  // from send_buffer_onedft_outputs (which was sorted by iParent).
  // After re-sorting by task_comparator, we need a different approach.
  // Actually, eps_on_grid_local was already in the send_buffer_onedft_outputs
  // task order (sorted by iParent). The tasks are now re-sorted.
  // We need to store eps per-task before re-sorting.
  //
  // WORKAROUND: Store per-task eps in task.feat before sorting.
  // Actually, the simpler approach: don't re-sort. Use the current task order.
  // The tasks were already sorted by iParent from prepare_onedft_features.
  // Let's just rebuild the eps_per_task mapping.

  // Build task -> eps mapping from the eps_on_grid vector (in iParent-sorted order)
  // First, re-sort back to iParent order to match eps_on_grid
  std::stable_sort( tasks.begin(), tasks.end(),
    [](const auto& a, const auto& b) { return a.iParent < b.iParent; });

  // Distribute eps_on_grid to per-task storage
  {
    size_t offset = 0;
    for (auto& task : tasks) {
      int64_t npts = task.points.size();
      task.feat.eps.resize(npts);
      std::copy(eps_on_grid.data() + offset,
                eps_on_grid.data() + offset + npts,
                task.feat.eps.begin());
      offset += npts;
    }
  }

  // Now sort by workload for the Pulay loop
  std::sort( tasks.begin(), tasks.end(), task_comparator );

  #pragma omp parallel
  {

  XCHostData<value_type> host_data;

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {

    auto& task = tasks[iT];
    const int32_t  npts    = task.points.size();
    const int32_t  nbe     = task.bfn_screening.nbe;
    const int32_t  nshells = task.bfn_screening.shell_list.size();

    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();
    const int32_t* shell_list = task.bfn_screening.shell_list.data();

    // Allocate memory for basis evaluation (up to hessian for GGA/MGGA)
    if (is_gga || is_mgga) {
      host_data.basis_eval.resize(10 * npts * nbe); // B, dB_xyz, d2B_6
      host_data.zmat.resize(4 * 2 * npts * nbe);    // xN, xN_xyz, xZ, xZ_xyz
    } else {
      host_data.basis_eval.resize(4 * npts * nbe);  // B, dB_xyz
      host_data.zmat.resize(2 * npts * nbe);         // xN, xZ
    }
    host_data.nbe_scr.resize(nbe * nbe);
    host_data.eps.resize(npts);

    auto* basis_eval = host_data.basis_eval.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* eps_buf    = host_data.eps.data();

    auto* dbasis_x_eval = basis_eval    + npts * nbe;
    auto* dbasis_y_eval = dbasis_x_eval + npts * nbe;
    auto* dbasis_z_eval = dbasis_y_eval + npts * nbe;

    value_type* d2basis_xx_eval = nullptr;
    value_type* d2basis_xy_eval = nullptr;
    value_type* d2basis_xz_eval = nullptr;
    value_type* d2basis_yy_eval = nullptr;
    value_type* d2basis_yz_eval = nullptr;
    value_type* d2basis_zz_eval = nullptr;

    if (is_gga || is_mgga) {
      d2basis_xx_eval = dbasis_z_eval   + npts * nbe;
      d2basis_xy_eval = d2basis_xx_eval + npts * nbe;
      d2basis_xz_eval = d2basis_xy_eval + npts * nbe;
      d2basis_yy_eval = d2basis_xz_eval + npts * nbe;
      d2basis_yz_eval = d2basis_yy_eval + npts * nbe;
      d2basis_zz_eval = d2basis_yz_eval + npts * nbe;
    }

    // X-matrix pointers
    auto* xNmat = host_data.zmat.data();
    value_type* xNmat_x = nullptr;
    value_type* xNmat_y = nullptr;
    value_type* xNmat_z = nullptr;
    value_type* xZmat   = nullptr;
    value_type* xZmat_x = nullptr;
    value_type* xZmat_y = nullptr;
    value_type* xZmat_z = nullptr;

    if (is_gga || is_mgga) {
      xNmat_x = xNmat   + npts*nbe;
      xNmat_y = xNmat_x + npts*nbe;
      xNmat_z = xNmat_y + npts*nbe;
      xZmat   = xNmat_z + npts*nbe;
      xZmat_x = xZmat   + npts*nbe;
      xZmat_y = xZmat_x + npts*nbe;
      xZmat_z = xZmat_y + npts*nbe;
    } else {
      xZmat = xNmat + npts*nbe;
    }

    // Get submat map
    auto [submat_map, foo] =
      gen_compressed_submat_map( basis_map, task.bfn_screening.shell_list, nbf, nbf );

    // Evaluate collocation (gradient + hessian for GGA/MGGA)
    if (is_gga || is_mgga) {
      lwd->eval_collocation_hessian( npts, nshells, nbe, points, basis, shell_list,
        basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval,
        d2basis_xx_eval, d2basis_xy_eval, d2basis_xz_eval,
        d2basis_yy_eval, d2basis_yz_eval, d2basis_zz_eval );
    } else {
      lwd->eval_collocation_gradient( npts, nshells, nbe, points, basis, shell_list,
        basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );
    }

    // Evaluate X-matrices: xN = Ps * B, xZ = Pz * B
    const int xmat_len = (is_gga || is_mgga) ? 4 : 1;
    lwd->eval_xmat( xmat_len*npts, nbf, nbe, submat_map, 1.0, Ps, ldps, basis_eval, nbe,
                    xNmat, nbe, nbe_scr );
    if (is_uks) {
      lwd->eval_xmat( xmat_len*npts, nbf, nbe, submat_map, 1.0, Pz, ldpz, basis_eval, nbe,
                      xZmat, nbe, nbe_scr );
    }

    // Read OneDFT Vxc from task.feat (already includes grid weights from autograd)
    const value_type* vdden_a = task.feat.vdden_eval_a.data();
    const value_type* vdden_b = task.feat.vdden_eval_b.data();
    const value_type* vdden_x_a = nullptr;
    const value_type* vdden_y_a = nullptr;
    const value_type* vdden_z_a = nullptr;
    const value_type* vdden_x_b = nullptr;
    const value_type* vdden_y_b = nullptr;
    const value_type* vdden_z_b = nullptr;
    const value_type* vtau_data = nullptr;

    if (is_gga || is_mgga) {
      vdden_x_a = task.feat.vdden_x_eval_a.data();
      vdden_y_a = task.feat.vdden_y_eval_a.data();
      vdden_z_a = task.feat.vdden_z_eval_a.data();
      vdden_x_b = task.feat.vdden_x_eval_b.data();
      vdden_y_b = task.feat.vdden_y_eval_b.data();
      vdden_z_b = task.feat.vdden_z_eval_b.data();
    }
    if (is_mgga) {
      vtau_data = task.feat.vtau.data();
    }

    // --- Weight derivative term ---
    // eps_contracted[ipt] = exc_on_grid[ipt] * w[ipt]
    for (int ipt = 0; ipt < npts; ++ipt) {
      eps_buf[ipt] = task.feat.eps[ipt] * weights[ipt];
    }
    lwd->eval_weight_1st_deriv_contracted( weight_alg, mol, molmeta,
      task, eps_buf, EXC_GRAD);

    // --- Pulay gradient loop ---
    // Using OneDFT's native per-component Vxc (weights already included)
    size_t bf_off = 0;
    for (auto ish = 0; ish < nshells; ++ish) {
      const int sh_idx = shell_list[ish];
      const int sh_sz  = basis[sh_idx].size();
      const int iAt    = basis_map.shell_to_center( sh_idx );

      // Skip basis functions on the parent atom (handled by weight derivative)
      if (iAt == task.iParent) {
        bf_off += sh_sz;
        continue;
      }

      double g_acc_x(0), g_acc_y(0), g_acc_z(0);

      for (int ibf = 0, mu = bf_off; ibf < sh_sz; ++ibf, ++mu)
      for (int ipt = 0; ipt < npts; ++ipt) {

        const int32_t mu_i = mu + ipt*nbe;

        // OneDFT Vxc: vdden_a = w * ∂ε/∂ρ_α, vdden_b = w * ∂ε/∂ρ_β
        // vrho_s = vdden_a + vdden_b (total density derivative, weighted)
        // vrho_z = vdden_a - vdden_b (magnetization derivative, weighted)
        const double vrho_s = vdden_a[ipt] + vdden_b[ipt];
        const double vrho_z = vdden_a[ipt] - vdden_b[ipt];

        const double xN = xNmat[mu_i];
        const double xZ = is_uks ? xZmat[mu_i] : 0.0;

        const double dbx = dbasis_x_eval[mu_i];
        const double dby = dbasis_y_eval[mu_i];
        const double dbz = dbasis_z_eval[mu_i];

        // LDA contribution (no separate weight multiplication — already in Vxc)
        g_acc_x += 0.5 * vrho_s * xN * dbx;
        g_acc_y += 0.5 * vrho_s * xN * dby;
        g_acc_z += 0.5 * vrho_s * xN * dbz;

        if (is_uks) {
          g_acc_x += 0.5 * vrho_z * xZ * dbx;
          g_acc_y += 0.5 * vrho_z * xZ * dby;
          g_acc_z += 0.5 * vrho_z * xZ * dbz;
        }

        if (is_gga || is_mgga) {
          // GGA contribution using OneDFT per-component derivatives
          // vdden_d_a = w * ∂ε/∂(∂ρ_α/∂d)
          // Force = Σ_d vdden_d_s * (d2B_{cd} * xN + dBc * xN_d) + z terms
          const double vds_x = 0.5 * (vdden_x_a[ipt] + vdden_x_b[ipt]);
          const double vds_y = 0.5 * (vdden_y_a[ipt] + vdden_y_b[ipt]);
          const double vds_z = 0.5 * (vdden_z_a[ipt] + vdden_z_b[ipt]);
          const double vdz_x = 0.5 * (vdden_x_a[ipt] - vdden_x_b[ipt]);
          const double vdz_y = 0.5 * (vdden_y_a[ipt] - vdden_y_b[ipt]);
          const double vdz_z = 0.5 * (vdden_z_a[ipt] - vdden_z_b[ipt]);

          const double xNx = xNmat_x[mu_i];
          const double xNy = xNmat_y[mu_i];
          const double xNz = xNmat_z[mu_i];
          const double xZx = is_uks ? xZmat_x[mu_i] : 0.0;
          const double xZy = is_uks ? xZmat_y[mu_i] : 0.0;
          const double xZz = is_uks ? xZmat_z[mu_i] : 0.0;

          const double d2bxx = d2basis_xx_eval[mu_i];
          const double d2bxy = d2basis_xy_eval[mu_i];
          const double d2bxz = d2basis_xz_eval[mu_i];
          const double d2byy = d2basis_yy_eval[mu_i];
          const double d2byz = d2basis_yz_eval[mu_i];
          const double d2bzz = d2basis_zz_eval[mu_i];

          // s (total) contribution: Σ_d vds_d * (d2B_{c,d} * xN + dBc * xN_d)
          // x-component of force:
          const double d2_xN_x = d2bxx * xN + dbx * xNx;
          const double d2_xN_y = d2bxy * xN + dbx * xNy;
          const double d2_xN_z = d2bxz * xN + dbx * xNz;
          g_acc_x += vds_x * d2_xN_x + vds_y * d2_xN_y + vds_z * d2_xN_z;

          // y-component:
          const double d2_yN_x = d2bxy * xN + dby * xNx;
          const double d2_yN_y = d2byy * xN + dby * xNy;
          const double d2_yN_z = d2byz * xN + dby * xNz;
          g_acc_y += vds_x * d2_yN_x + vds_y * d2_yN_y + vds_z * d2_yN_z;

          // z-component:
          const double d2_zN_x = d2bxz * xN + dbz * xNx;
          const double d2_zN_y = d2byz * xN + dbz * xNy;
          const double d2_zN_z = d2bzz * xN + dbz * xNz;
          g_acc_z += vds_x * d2_zN_x + vds_y * d2_zN_y + vds_z * d2_zN_z;

          if (is_uks) {
            // z (magnetization) contribution: Σ_d vdz_d * (d2B_{c,d} * xZ + dBc * xZ_d)
            const double d2_xZ_x = d2bxx * xZ + dbx * xZx;
            const double d2_xZ_y = d2bxy * xZ + dbx * xZy;
            const double d2_xZ_z = d2bxz * xZ + dbx * xZz;
            g_acc_x += vdz_x * d2_xZ_x + vdz_y * d2_xZ_y + vdz_z * d2_xZ_z;

            const double d2_yZ_x = d2bxy * xZ + dby * xZx;
            const double d2_yZ_y = d2byy * xZ + dby * xZy;
            const double d2_yZ_z = d2byz * xZ + dby * xZz;
            g_acc_y += vdz_x * d2_yZ_x + vdz_y * d2_yZ_y + vdz_z * d2_yZ_z;

            const double d2_zZ_x = d2bxz * xZ + dbz * xZx;
            const double d2_zZ_y = d2byz * xZ + dbz * xZy;
            const double d2_zZ_z = d2bzz * xZ + dbz * xZz;
            g_acc_z += vdz_x * d2_zZ_x + vdz_y * d2_zZ_y + vdz_z * d2_zZ_z;
          }
        }

        if (is_mgga) {
          // MGGA τ contribution
          // vtau is interleaved [α₀, β₀, α₁, β₁, ...]
          const double vtaup = 0.5 * vtau_data[2*ipt];     // α, already weighted
          const double vtaum = 0.5 * vtau_data[2*ipt + 1];  // β, already weighted
          const double vtaun = vtaup + vtaum;
          const double vtauz = vtaup - vtaum;

          const double xNx = xNmat_x[mu_i];
          const double xNy = xNmat_y[mu_i];
          const double xNz = xNmat_z[mu_i];

          const double d2bxx = d2basis_xx_eval[mu_i];
          const double d2bxy = d2basis_xy_eval[mu_i];
          const double d2bxz = d2basis_xz_eval[mu_i];
          const double d2byy = d2basis_yy_eval[mu_i];
          const double d2byz = d2basis_yz_eval[mu_i];
          const double d2bzz = d2basis_zz_eval[mu_i];

          auto d2_term_x = d2bxx * xNx + d2bxy * xNy + d2bxz * xNz;
          auto d2_term_y = d2bxy * xNx + d2byy * xNy + d2byz * xNz;
          auto d2_term_z = d2bxz * xNx + d2byz * xNy + d2bzz * xNz;

          g_acc_x += 0.5 * vtaun * d2_term_x;
          g_acc_y += 0.5 * vtaun * d2_term_y;
          g_acc_z += 0.5 * vtaun * d2_term_z;

          if (is_uks) {
            const double xZx = xZmat_x[mu_i];
            const double xZy = xZmat_y[mu_i];
            const double xZz = xZmat_z[mu_i];

            d2_term_x = d2bxx * xZx + d2bxy * xZy + d2bxz * xZz;
            d2_term_y = d2bxy * xZx + d2byy * xZy + d2byz * xZz;
            d2_term_z = d2bxz * xZx + d2byz * xZy + d2bzz * xZz;

            g_acc_x += 0.5 * vtauz * d2_term_x;
            g_acc_y += 0.5 * vtauz * d2_term_y;
            g_acc_z += 0.5 * vtauz * d2_term_z;
          }
        }

      } // end loop over bfns + grid points

      #pragma omp atomic
      EXC_GRAD[3*iAt + 0] += -2 * g_acc_x;
      #pragma omp atomic
      EXC_GRAD[3*iAt + 1] += -2 * g_acc_y;
      #pragma omp atomic
      EXC_GRAD[3*iAt + 2] += -2 * g_acc_z;

      // Weight derivative counterpart for non-parent atoms
      #pragma omp atomic
      EXC_GRAD[3*task.iParent + 0] -= -2 * g_acc_x;
      #pragma omp atomic
      EXC_GRAD[3*task.iParent + 1] -= -2 * g_acc_y;
      #pragma omp atomic
      EXC_GRAD[3*task.iParent + 2] -= -2 * g_acc_z;

      bf_off += sh_sz;

    } // end loop over shells

  } // end loop over tasks

  } // end OpenMP region
}

} // namespace detail
} // namespace GauXC
