/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
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
  neo_eval_exc_vxc_( int64_t elec_m, int64_t elec_n, int64_t prot_m, int64_t prot_n, 
                     const value_type* elec_Ps, int64_t elec_ldps,
                     const value_type* prot_Ps, int64_t prot_ldps,
                     const value_type* prot_Pz, int64_t prot_ldpz,
                     value_type* elec_VXCs,     int64_t elec_ldvxcs,
                     value_type* prot_VXCs,     int64_t prot_ldvxcs,
                     value_type* prot_VXCz,     int64_t prot_ldvxcz,
                     value_type* elec_EXC,  value_type* prot_EXC, const IntegratorSettingsXC& settings){
  
  const auto& elec_basis = this->load_balancer_->basis();
  const auto& prot_basis = this->load_balancer_->protonic_basis();

  // Check that P / VXC are sane
  const int64_t elec_nbf = elec_basis.nbf();
  const int64_t prot_nbf = prot_basis.nbf();

  if( elec_m != elec_n   | prot_m != prot_n)
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Be Square");
  if( elec_m != elec_nbf | prot_m != prot_nbf)
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Have Same Dimension as Basis");
  if( elec_ldps < elec_nbf | prot_ldps < prot_nbf | prot_ldpz < prot_nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");
  if( elec_ldvxcs < elec_nbf | prot_ldvxcs < prot_nbf | prot_ldvxcz < prot_nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXC");

  // Get Tasks
  this->load_balancer_->get_tasks();

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;
  // Temporary proton count to judge integrator accuracy
  value_type N_PROT;

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    neo_exc_vxc_local_work_( elec_Ps,   elec_ldps,
                             nullptr,   0,
                             prot_Ps,   prot_ldps,
                             prot_Pz,   prot_ldpz,
                             elec_VXCs, elec_ldvxcs,
                             nullptr,   0,
                             prot_VXCs, prot_ldvxcs,
                             prot_VXCz, prot_ldvxcz,
                             elec_EXC,  prot_EXC, 
                             &N_EL, &N_PROT, settings  );
  });


  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");

    this->reduction_driver_->allreduce_inplace( elec_VXCs, elec_nbf*elec_nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( prot_VXCs, prot_nbf*prot_nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( prot_VXCz, prot_nbf*prot_nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( elec_EXC,  1,                 ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( prot_EXC,  1,                 ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_EL,     1,                 ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_PROT,   1,                 ReductionOp::Sum );

  });

}

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  neo_eval_exc_vxc_( int64_t elec_m, int64_t elec_n, int64_t prot_m, int64_t prot_n, 
                     const value_type* elec_Ps, int64_t elec_ldps,
                     const value_type* elec_Pz, int64_t elec_ldpz,
                     const value_type* prot_Ps, int64_t prot_ldps,
                     const value_type* prot_Pz, int64_t prot_ldpz,
                     value_type* elec_VXCs,     int64_t elec_ldvxcs,
                     value_type* elec_VXCz,     int64_t elec_ldvxcz,
                     value_type* prot_VXCs,     int64_t prot_ldvxcs,
                     value_type* prot_VXCz,     int64_t prot_ldvxcz,
                     value_type* elec_EXC,  value_type* prot_EXC, const IntegratorSettingsXC& settings ){
  
  const auto& elec_basis = this->load_balancer_->basis();
  const auto& prot_basis = this->load_balancer_->protonic_basis();

  // Check that P / VXC are sane
  const int64_t elec_nbf = elec_basis.nbf();
  const int64_t prot_nbf = prot_basis.nbf();

  if( elec_m != elec_n   | prot_m != prot_n)
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Be Square");
  if( elec_m != elec_nbf | prot_m != prot_nbf)
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Have Same Dimension as Basis");
  if( elec_ldps < elec_nbf | elec_ldpz < elec_nbf | prot_ldps < prot_nbf | prot_ldpz < prot_nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");
  if( elec_ldvxcs < elec_nbf | elec_ldvxcz < elec_nbf | prot_ldvxcs < prot_nbf | prot_ldvxcz < prot_nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXC");

  // Get Tasks
  this->load_balancer_->get_tasks();

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;
  // Temporary proton count to judge integrator accuracy
  value_type N_PROT;

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    neo_exc_vxc_local_work_( elec_Ps,   elec_ldps,
                             elec_Pz,   elec_ldpz,
                             prot_Ps,   prot_ldps,
                             prot_Pz,   prot_ldpz,
                             elec_VXCs, elec_ldvxcs,
                             elec_VXCz, elec_ldvxcz,
                             prot_VXCs, prot_ldvxcs,
                             prot_VXCz, prot_ldvxcz,
                             elec_EXC,  prot_EXC,
                             &N_EL, &N_PROT, settings  );
  });


  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");

    this->reduction_driver_->allreduce_inplace( elec_VXCs, elec_nbf*elec_nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( elec_VXCz, elec_nbf*elec_nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( prot_VXCs, prot_nbf*prot_nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( prot_VXCz, prot_nbf*prot_nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( elec_EXC,  1,                 ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( prot_EXC,  1,                 ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_EL,     1,                 ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_PROT,   1,                 ReductionOp::Sum );

  });

}


template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  neo_exc_vxc_local_work_( const value_type* elec_Ps, int64_t elec_ldps,
                           const value_type* elec_Pz, int64_t elec_ldpz,
                           const value_type* prot_Ps, int64_t prot_ldps,
                           const value_type* prot_Pz, int64_t prot_ldpz,
                           value_type* elec_VXCs,     int64_t elec_ldvxcs,
                           value_type* elec_VXCz,     int64_t elec_ldvxcz,
                           value_type* prot_VXCs,     int64_t prot_ldvxcs,
                           value_type* prot_VXCz,     int64_t prot_ldvxcz,
                           value_type* elec_EXC,  value_type* prot_EXC,  
                           value_type *N_EL,      value_type *N_PROT,
                           const IntegratorSettingsXC& settings ) {
  
  // Determine is electronic subsystem is RKS or UKS
  const bool is_uks = (elec_Pz != nullptr) and (elec_VXCz != nullptr);
  const bool is_rks = not is_uks; 
  // TODO: Integrate with GKS
  // TODO: Integrate with mGGA

  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  // Setup Aliases
  const auto& func   = *this->func_;
  const auto& epcfunc   = *this->epcfunc_;
  const auto& basis  = this->load_balancer_->basis();
  const auto& mol    = this->load_balancer_->molecule();

  // Get basis map
  BasisSetMap basis_map(basis,mol);
  
  const int32_t nbf = basis.nbf();

  // Get Protonic basis information
  const auto& protonic_basis = this->load_balancer_->protonic_basis();
  BasisSetMap protonic_basis_map(protonic_basis,mol);
  const int32_t protonic_nbf = protonic_basis.nbf();

  // Sort tasks on size (XXX: maybe doesnt matter?)
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.bfn_screening.nbe) > (b.points.size() * b.bfn_screening.nbe);
  };
  
  auto& tasks = this->load_balancer_->get_tasks();
  std::sort( tasks.begin(), tasks.end(), task_comparator );


  // Check that Partition Weights have been calculated
  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored ) {
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Beed Modified");
  }
  

  // Zero out integrands
  
  std::fill(elec_VXCs, elec_VXCs + nbf * elec_ldvxcs, 0.0);
  if(is_uks) std::fill(elec_VXCz, elec_VXCz + nbf * elec_ldvxcz, 0.0);

  std::fill(prot_VXCs, prot_VXCs + protonic_nbf * prot_ldvxcs, 0.0);
  std::fill(prot_VXCz, prot_VXCz + protonic_nbf * prot_ldvxcz, 0.0);

  *elec_EXC = 0.;
  *prot_EXC = 0.;

  double NEL_WORK = 0.0;
  double NPROT_WORK = 0.0;
  double EXC_WORK = 0.0;
  double EPC_WORK = 0.0;
    
  // Loop over tasks
  const size_t ntasks = tasks.size();

  #pragma omp parallel
  {

  XCHostData<value_type> host_data; // Thread local host data

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {

    //std::cout << iT << "/" << ntasks << std::endl;
    //printf("%lu / %lu\n", iT, ntasks);
    // Alias current task
    const auto& task = tasks[iT];

    // Get tasks constants
    const int32_t  npts     = task.points.size();
    const int32_t  nbe      = task.bfn_screening.nbe;
    const int32_t  nshells  = task.bfn_screening.shell_list.size();

    const auto* points        = task.points.data()->data();
    const auto* weights       = task.weights.data();
    const int32_t* shell_list = task.bfn_screening.shell_list.data();

    const int32_t  protonic_nbe        = task.protonic_bfn_screening.nbe;
    const int32_t  protonic_nshells    = task.protonic_bfn_screening.shell_list.size();
    const int32_t* protonic_shell_list = task.protonic_bfn_screening.shell_list.data();

    // Check if there's protonic shells to evaluate. If not, only electronic EXC/VXC will be calculated
    bool evalProtonic =  (protonic_nshells != 0);

    // Allocate enough memory for batch

    const size_t spin_dim_scal = is_rks ? 1 : 2;
    
    // Things that every calc needs
    

    // Use same scratch for both electronic and protonic
    const int32_t scr_dim = std::max(nbe, protonic_nbe);
    host_data.nbe_scr .resize(scr_dim  * scr_dim);
    
    //----------------------Start Electronic System Setup------------------------
    host_data.zmat    .resize(npts * nbe * spin_dim_scal); 
    host_data.eps     .resize(npts);
    host_data.vrho    .resize(npts * spin_dim_scal);

    // LDA data requirements
    if( func.is_lda() ){
      host_data.basis_eval .resize( npts * nbe );
      host_data.den_scr    .resize( npts * spin_dim_scal);
    }

    // GGA data requirements
    const size_t gga_dim_scal = is_rks ? 1 : 3;
    if( func.is_gga() ){
      host_data.basis_eval .resize( 4 * npts * nbe );
      host_data.den_scr    .resize( spin_dim_scal * 4 * npts );
      host_data.gamma      .resize( gga_dim_scal * npts );
      host_data.vgamma     .resize( gga_dim_scal * npts );
    }

    // Alias/Partition out scratch memory
    auto* basis_eval = host_data.basis_eval.data();
    auto* den_eval   = host_data.den_scr.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* zmat       = host_data.zmat.data();

    decltype(zmat) zmat_z = nullptr;
    if(!is_rks) {
      zmat_z = zmat + nbe * npts;
    }

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
      dden_x_eval   = den_eval    + spin_dim_scal * npts;
      dden_y_eval   = dden_x_eval + spin_dim_scal * npts;
      dden_z_eval   = dden_y_eval + spin_dim_scal * npts;
    }
    
    //----------------------End Electronic System Setup------------------------


    //----------------------Start Protonic System Setup------------------------
    // Set Up Memory (assuming UKS)
    host_data.epc                .resize( npts );
    host_data.protonic_vrho      .resize( npts * 2 );
    host_data.protonic_zmat      .resize( npts * protonic_nbe * 2 );
    // LDA
    host_data.protonic_basis_eval .resize( npts * protonic_nbe );
    host_data.protonic_den_scr    .resize( npts * 2);
    // Alias/Partition out scratch memory
    auto* protonic_basis_eval = host_data.protonic_basis_eval.data();
    auto* protonic_den_eval   = host_data.protonic_den_scr.data();
    auto* protonic_zmat       = host_data.protonic_zmat.data();
    decltype(protonic_zmat) protonic_zmat_z = protonic_zmat + protonic_nbe * npts;
    
    auto* epc                 = host_data.epc.data();
    auto* protonic_vrho       = host_data.protonic_vrho.data();
    // No GGA for NEO yet
    //----------------------End Protonic System Setup------------------------


    //std::cout << "Task: " << iT << "/" << ntasks << std::endl;
    //std::cout << "Electronic nbe: " << nbe << std::endl;
    //std::cout << "Protonic   nbe: " << protonic_nbe << std::endl;

    //----------------------Start Calculating Electronic Density & UV Variable------------------------
    // Get the submatrix map for batch
    std::vector< std::array<int32_t, 3> > submat_map;
    std::tie(submat_map, std::ignore) =
          gen_compressed_submat_map(basis_map, task.bfn_screening.shell_list, nbf, nbf);

    // Evaluate Collocation (+ Grad)
    if( func.is_gga() )
      lwd->eval_collocation_gradient( npts, nshells, nbe, points, basis, shell_list, 
        basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );
    else
      lwd->eval_collocation( npts, nshells, nbe, points, basis, shell_list, 
        basis_eval );
    
    // Evaluate X matrix (fac * P * B) -> store in Z
    const auto xmat_fac = is_rks ? 2.0 : 1.0; // TODO Fix for spinor RKS input
    lwd->eval_xmat( npts, nbf, nbe, submat_map, xmat_fac, elec_Ps, elec_ldps, basis_eval, nbe,
      zmat, nbe, nbe_scr );

    // X matrix for Pz
    if(not is_rks) {
      lwd->eval_xmat( npts, nbf, nbe, submat_map, 1.0, elec_Pz, elec_ldpz, basis_eval, nbe,
        zmat_z, nbe, nbe_scr );
    }

    // Evaluate U and V variables
    if( func.is_gga() ) {
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

    // Evaluate XC functional
    if( func.is_gga() )
      func.eval_exc_vxc( npts, den_eval, gamma, eps, vrho, vgamma );
    else
      func.eval_exc_vxc( npts, den_eval, eps, vrho );
    //----------------------End Calculating Electronic Density & UV Variable------------------------




    //----------------------Start Calculating Protonic Density & UV Variable------------------------
    std::vector< std::array<int32_t, 3> > protonic_submat_map;
    if(evalProtonic){
      std::tie(protonic_submat_map, std::ignore) =
        gen_compressed_submat_map(protonic_basis_map, task.protonic_bfn_screening.shell_list, protonic_nbf, protonic_nbf);

      // Evaluate Collocation 
      lwd->eval_collocation( npts, protonic_nshells, protonic_nbe, points, protonic_basis, protonic_shell_list,
        protonic_basis_eval );

      // Evaluate X matrix (P * B) -> store in Z
      lwd->eval_xmat( npts, protonic_nbf, protonic_nbe, protonic_submat_map, 1.0, prot_Ps, prot_ldps, protonic_basis_eval, protonic_nbe,
        protonic_zmat,   protonic_nbe, nbe_scr );
      lwd->eval_xmat( npts, protonic_nbf, protonic_nbe, protonic_submat_map, 1.0, prot_Pz, prot_ldpz, protonic_basis_eval, protonic_nbe,
        protonic_zmat_z, protonic_nbe, nbe_scr );

      // Evaluate U and V variables
      lwd->eval_uvvar_lda_uks( npts, protonic_nbe, protonic_basis_eval, protonic_zmat, protonic_nbe, protonic_zmat_z, 
        protonic_nbe, protonic_den_eval );

      // No protonic XC functional. Fill with eps and vrho to be 0.0
      std::fill_n(epc,  npts,   0.);
      std::fill_n(protonic_vrho, npts*2, 0.);
    }
    //----------------------End Calculating Protonic Density & UV Variable------------------------




    //----------------------Start EPC functional Evaluation---------------------------------------
    if(evalProtonic){
      // Prepare for kernal input
      for (int32_t iPt = 0; iPt < npts; iPt++ ){
        // Treat total erho as spin-up, treat total prho as spin down
        protonic_den_eval[2*iPt+1] = protonic_den_eval[2*iPt] + protonic_den_eval[2*iPt+1];
        protonic_den_eval[2*iPt]   = is_rks ? den_eval[iPt] : (den_eval[2*iPt] + den_eval[2*iPt+1]);
      }

      // EPC Functional Evaluation (Calling ExchCXX Builtin Function)
      epcfunc.eval_exc_vxc( npts, protonic_den_eval, epc, protonic_vrho );

      // Digest kernal output
      for (int32_t iPt = 0; iPt < npts; iPt++ ){
        // assign df/derho 
        vrho[spin_dim_scal*iPt]   += protonic_vrho[2*iPt];
        if(not is_rks) vrho[spin_dim_scal*iPt+1] += protonic_vrho[2*iPt];
        // assign df/dprho
        protonic_vrho[2*iPt] = protonic_vrho[2*iPt+1];
        protonic_vrho[2*iPt+1] = 0.0;
        // change back protonic density to original state
        protonic_den_eval[2*iPt] = protonic_den_eval[2*iPt+1];
        protonic_den_eval[2*iPt] = 0.0;
      }
    } // End if(evalProtonic)
    //----------------------End EPC functional Evaluation---------------------------------------
 


    //----------------------Begin Evaluating Electronic ZMat---------------------------- 
    // Factor weights into XC results
    for( int32_t i = 0; i < npts; ++i ) {
      eps[i]  *= weights[i];
      epc[i]  *= weights[i];
      vrho[spin_dim_scal*i] *= weights[i];
      if(not is_rks) vrho[spin_dim_scal*i+1] *= weights[i];
    }

    if( func.is_gga() ){
      for( int32_t i = 0; i < npts; ++i ) {
         vgamma[gga_dim_scal*i] *= weights[i];
         if(not is_rks) {
           vgamma[gga_dim_scal*i+1] *= weights[i];
           vgamma[gga_dim_scal*i+2] *= weights[i];
         }
      }
    }

    // Evaluate Z matrix for VXC
    if( func.is_gga() ) {
      if(is_rks) {
        lwd->eval_zmat_gga_vxc_rks( npts, nbe, vrho, vgamma, basis_eval, dbasis_x_eval,
                                dbasis_y_eval, dbasis_z_eval, dden_x_eval, dden_y_eval,
                                dden_z_eval, zmat, nbe);
      } else if(is_uks) {
        lwd->eval_zmat_gga_vxc_uks( npts, nbe, vrho, vgamma, basis_eval, dbasis_x_eval,
                                dbasis_y_eval, dbasis_z_eval, dden_x_eval, dden_y_eval,
                                dden_z_eval, zmat, nbe, zmat_z, nbe);
      }
    } else {
      if(is_rks) {
        lwd->eval_zmat_lda_vxc_rks( npts, nbe, vrho, basis_eval, zmat, nbe );
      } else if(is_uks) {
        lwd->eval_zmat_lda_vxc_uks( npts, nbe, vrho, basis_eval, zmat, nbe, zmat_z, nbe );
      }
    }
    //----------------------End Evaluating Electronic ZMat----------------------------




    //----------------------Begin Evaluating Protonic ZMat---------------------------- 
    if(evalProtonic){
      // Factor weights into XC results
      for( int32_t i = 0; i < npts; ++i ) {
        protonic_vrho[2*i]   *= weights[i];
        protonic_vrho[2*i+1] *= weights[i];
      }

      // Evaluate Z matrix for VXC
      lwd->eval_zmat_lda_vxc_uks( npts, protonic_nbe, protonic_vrho, protonic_basis_eval, protonic_zmat, protonic_nbe, 
        protonic_zmat_z, protonic_nbe );
    }
    //----------------------End Evaluating Protonic ZMat----------------------------

    // Scalar integrations
    double NEL_local   = 0.0;
    double NPROT_local = 0.0;
    double EXC_local   = 0.0;
    double EPC_local   = 0.0;
    for( int32_t i = 0; i < npts; ++i ) {
      const auto den  = is_rks ? den_eval[i] : (den_eval[2*i] + den_eval[2*i+1]);
      NEL_local    += weights[i] * den;
      EXC_local    += eps[i]     * den;
      EPC_local    += epc[i]     * den;;
    }
    // Protonic XC (EPC)
    if(evalProtonic){
      for( int32_t i = 0; i < npts; ++i ) {
        const auto protonic_den =  protonic_den_eval[2*i] + protonic_den_eval[2*i+1];
        NPROT_local            +=  weights[i]  * protonic_den;
        EPC_local              +=  epc[i]      * protonic_den;    
      }
    }

    // Atomic updates
    #pragma omp atomic
    NEL_WORK   += NEL_local;
    #pragma omp atomic
    NPROT_WORK += NPROT_local;
    #pragma omp atomic
    EXC_WORK   += EXC_local;
    #pragma omp atomic
    EPC_WORK   += EPC_local;

    
    // Incremeta LT of VXC
    {

      // Increment Electronic XC (VXC+EPC)
      lwd->inc_vxc( npts, nbf, nbe, basis_eval, submat_map, zmat, nbe, elec_VXCs, elec_ldvxcs, nbe_scr );
      if(not is_rks) 
        lwd->inc_vxc( npts, nbf, nbe, basis_eval, submat_map, zmat_z, nbe, elec_VXCz, elec_ldvxcz, nbe_scr );

      // Increment Protonic XC (EPC)
      // Scalar integrations
      if(evalProtonic){
        lwd->inc_vxc( npts, protonic_nbf, protonic_nbe, protonic_basis_eval, protonic_submat_map, 
          protonic_zmat,   protonic_nbe, prot_VXCs, prot_ldvxcs, nbe_scr );
        lwd->inc_vxc( npts, protonic_nbf, protonic_nbe, protonic_basis_eval, protonic_submat_map, 
          protonic_zmat_z, protonic_nbe, prot_VXCz, prot_ldvxcz, nbe_scr );
      }

    }

  } // Loop over tasks
  
  }  // End OpenMP region

  // Set scalar return values
  *N_EL      = NEL_WORK;
  *N_PROT    = NPROT_WORK;
  *elec_EXC  = EXC_WORK + EPC_WORK;
  *prot_EXC  = EPC_WORK;

  //std::cout << "N_EL = "     << std::setprecision(12) << std::scientific << *N_EL       << std::endl;
  //std::cout << "N_PROT = "   << std::setprecision(12) << std::scientific << *N_PROT     << std::endl;
  //std::cout << "elec_EXC = " << std::setprecision(12) << std::scientific << *elec_EXC   << std::endl;
  //std::cout << "prot_EXC = " << std::setprecision(12) << std::scientific << *prot_EXC   << std::endl;

  // Symmetrize Electronic VXC
  for( int32_t j = 0;   j < nbf; ++j )
  for( int32_t i = j+1; i < nbf; ++i )
    elec_VXCs[ j + i*elec_ldvxcs ] = elec_VXCs[ i + j*elec_ldvxcs ];
  if(not is_rks) {
  for( int32_t j = 0;   j < nbf; ++j ) {
  for( int32_t i = j+1; i < nbf; ++i ) {
    elec_VXCz[ j + i*elec_ldvxcz ] = elec_VXCz[ i + j*elec_ldvxcz ];
  }
  }
  }

  // Symmetrize Protonic VXC
  for( int32_t j = 0;   j < protonic_nbf; ++j ){
  for( int32_t i = j+1; i < protonic_nbf; ++i ){
    prot_VXCs[ j + i*prot_ldvxcs ] = prot_VXCs[ i + j*prot_ldvxcs ];
    prot_VXCz[ j + i*prot_ldvxcz ] = prot_VXCz[ i + j*prot_ldvxcz ];
  }
  }

}


}
}
