/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "incore_replicated_xc_device_integrator.hpp"
#include <gauxc/util/misc.hpp>
#include <gauxc/util/unused.hpp>
#include "device/scheme1_data_base.hpp"
#include "device/common/device_blas.hpp"
#include "integrator_util/onedft_util.hpp"
#include <cuda_runtime.h>
#include "device/cuda/cuda_backend.hpp"
#include <cstddef> // for size_t

namespace GauXC::detail {

FeatureDict prepare_onedft_features( const size_t natoms, const size_t total_npts, const size_t ndm,
  const at::TensorOptions options, const std::vector<std::string> feature_keys,
  double* den_eval, double* dden_eval, double* tau, double* grid_coords, 
  double* grid_weights, double* coords,
  const std::vector<int64_t>& atomic_grid_sizes = {},
  int64_t max_grid_size = 0,
  double* raw_grid_weights = nullptr );

size_t save_static_data_onedft_features (XCDeviceData* _data, const integrator_term_tracker enabled_terms, size_t offset);

void save_static_data_onedft_outputs(const at::Tensor exc, const FeatureDict& features_dict, XCDeviceData* _data);

size_t send_buffer_onedft_outputs (XCDeviceData* _data, const integrator_term_tracker enabled_terms, size_t offset);

void* my_malloc(size_t size, int device, cudaStream_t stream) {
  void *ptr;
  cudaMallocAsync(&ptr, size, stream);
  return ptr;
}
void my_free(void* ptr, ssize_t size, int device, cudaStream_t stream) {
  cudaFreeAsync(ptr, stream);
}
void init_custom_allocator() {
  setenv("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True", 1);
  // auto custom_allocator = torch::cuda::CUDAPluggableAllocator::createCustomAllocator(my_malloc, my_free);
  // torch::cuda::CUDAPluggableAllocator::changeCurrentAllocator(custom_allocator);
}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
eval_exc_vxc_onedft_( int64_t m, int64_t n, 
  const value_type* Ps, int64_t ldps,
  const value_type* Pz, int64_t ldpz,
  value_type* VXCs, int64_t ldvxcs,
  value_type* VXCz, int64_t ldvxcz,
  value_type* EXC, const IntegratorSettingsXC& settings ) {
  
  const bool is_uks = (Pz != nullptr);
  const bool is_rks = (Ps != nullptr) and (not is_uks);
  if (is_rks) { // TODO
    GAUXC_GENERIC_EXCEPTION("RKS Not Yet Implemented");
  }

  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
  const int64_t nbf = basis.nbf();
  if( m != n ) 
    GAUXC_GENERIC_EXCEPTION("P Must Be Square");
  if( m != nbf ) 
    GAUXC_GENERIC_EXCEPTION("P Must Have Same Dimension as Basis");

  if( ldps and ldps < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");
  if( ldvxcs < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCS");

  if( not is_rks ) {
    if( ldpz and ldpz < nbf )
      GAUXC_GENERIC_EXCEPTION("Invalid LDPZ");
    if( ldvxcz and ldvxcz < nbf )
      GAUXC_GENERIC_EXCEPTION("Invalid LDVXCZ");
  }

  // Get Tasks
  auto& tasks = this->load_balancer_->get_tasks();
  // Sort tasks by atom index so that grid points are grouped by atom.
  // build_atom_reorder_perm assumes this contiguous-by-atom layout.
  std::stable_sort(tasks.begin(), tasks.end(),
    [](const auto& a, const auto& b) { return a.iParent < b.iParent; });
  size_t total_npts = std::accumulate( tasks.begin(), tasks.end(), 0ul,
    [](const auto& a, const auto& b) { return a + b.npts; } );

  // Check that Partition Weights have been calculated
  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored ) {
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Been Modified"); 
  }

  // Allocate Device memory
  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
  auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
  auto device_data_ptr = lwd->create_device_data(rt);

  integrator_term_tracker enabled_terms;
  enabled_terms.exc_vxc = true;
  enabled_terms.onedft = true;

  if (is_rks) enabled_terms.ks_scheme = RKS;
  else if (is_uks) enabled_terms.ks_scheme = UKS;

  // load onedft model
  OneDFTSettings onedft_settings;
  if( auto* tmp = dynamic_cast<const OneDFTSettings*>(&settings) ) {
    onedft_settings = *tmp;
  }
  const auto model_path = onedft_settings.model;
  if (not torch::cuda::is_available()) {
    GAUXC_GENERIC_EXCEPTION("Torch CUDA Not Available");
  }
  torch::DeviceType torch_device = torch::kCUDA;
  init_custom_allocator();
  auto [exc_func, feature_keys] = load_model(model_path, torch_device);
  
  // determine what feature we need based on the keys
  if (feature_keys.size() == 0) {
    GAUXC_GENERIC_EXCEPTION("No feature keys found in model");
  }

  bool is_gga = false;
  bool is_mgga = false;
  bool is_lda = false;

  for (const auto& key : feature_keys) {
    if ( not valueExists(key) ) {
      GAUXC_GENERIC_EXCEPTION("Feature Key Required Not Implemented: " + key);
    }
    if (key == feat_map.at(ONEDFT_FEATURE::TAU)) {
      is_mgga = true;
    }
    if (key == feat_map.at(ONEDFT_FEATURE::DDEN)) {
      is_gga = true;
    }
  }

  if (is_mgga) {
    enabled_terms.xc_approx = integrator_xc_approx::MGGA_TAU;
    is_gga = false;
  } else if (is_gga)
    enabled_terms.xc_approx = integrator_xc_approx::GGA;
  else {
    is_lda = true;
    enabled_terms.xc_approx = integrator_xc_approx::LDA;
  }

  const auto& mol   = this->load_balancer_->molecule();
  const auto natoms = mol.natoms();
  const auto nshells = basis.nshells();
  // alocate onedft memory
  device_data_ptr->reset_allocations();
  device_data_ptr->allocate_static_data_onedft( nbf, nshells, natoms, total_npts, enabled_terms );
  device_data_ptr->send_static_data_onedft( mol, Ps, ldps, Pz, ldpz, nullptr, 0, nullptr, 0, basis );
  // Zero integrands
  device_data_ptr->zero_exc_vxc_integrands(enabled_terms);
     
  this->timer_.time_op("XCIntegrator.LocalWork_PreOneDFT", [&](){
    pre_onedft_local_work_( basis, Ps, ldps, Pz, ldpz, nullptr, 0, nullptr, 0,
      tasks.begin(), tasks.end(), *device_data_ptr, enabled_terms );
  });

  int32_t world_rank = rt.comm_rank();
  int32_t world_size = rt.comm_size();
  size_t ndm = enabled_terms.ks_scheme == UKS ? 2 : 1;

  std::vector<double> grid_weights, grid_coords, den_eval, dden_eval, tau;
  std::vector<int> displs(world_size), recvcounts(world_size);

  // Collect raw (pre-partition) quadrature weights from tasks
  std::vector<double> raw_grid_weights;
  raw_grid_weights.reserve(total_npts);
  for (const auto& task : tasks) {
    raw_grid_weights.insert(raw_grid_weights.end(), task.raw_weights.begin(), task.raw_weights.end());
  }

  // Compute per-atom grid sizes from tasks
  std::vector<int64_t> atomic_grid_sizes_vec(natoms, 0);
  for (const auto& task : tasks) {
    if (task.iParent >= 0 && task.iParent < (int)natoms) {
      atomic_grid_sizes_vec[task.iParent] += task.npts;
    }
  }

  // run onedft model on rank 0
  FeatureDict features_dict;
  std::vector<int64_t> atom_reorder_inv_perm;
  std::vector<int64_t> global_atomic_grid_sizes_vec = atomic_grid_sizes_vec;

  if ( world_size == 1 ) { // keep everything on device
    int64_t max_grid_size = *std::max_element(
      atomic_grid_sizes_vec.begin(), atomic_grid_sizes_vec.end());
    auto options = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCUDA);
    features_dict = prepare_onedft_features(
      natoms, total_npts, ndm, options, feature_keys, device_data_ptr->den_eval_device_data(),
      device_data_ptr->dden_eval_device_data(), device_data_ptr->tau_device_data(),
      device_data_ptr->grid_coords_device_data(), device_data_ptr->grid_weights_device_data(),
      device_data_ptr->coords_device_data(),
      atomic_grid_sizes_vec, max_grid_size,
      raw_grid_weights.data()
    );
  } else { // copy to host, gather, reorder, then back to device
    grid_weights.resize(total_npts);
    grid_coords.resize(total_npts * 3);
    den_eval.resize(total_npts * ndm);
    if (is_gga | is_mgga) {
      dden_eval.resize(total_npts * ndm * 3);
    }
    if (is_mgga) {
      tau.resize(total_npts * ndm);
    }
    device_data_ptr->retrieve_onedft_features( total_npts, 2, den_eval.data(),
      (is_gga || is_mgga) ? dden_eval.data() : nullptr,
      is_mgga ? tau.data() : nullptr,
      grid_coords.data(), grid_weights.data() );

    // Build host coords from molecule (avoids device-to-host copy of coords_device_data)
    std::vector<double> host_coords(natoms * 3);
    for (size_t i = 0; i < natoms; i++) {
      host_coords[3*i]   = mol[i].x;
      host_coords[3*i+1] = mol[i].y;
      host_coords[3*i+2] = mol[i].z;
    }

    auto reorder_result = mpi_gather_and_reorder_gpu(
      den_eval, dden_eval, tau, grid_coords, grid_weights,
      atomic_grid_sizes_vec, total_npts, natoms, rt, recvcounts, displs);
    int total_npts_sum = reorder_result.total_npts;
    atom_reorder_inv_perm = std::move(reorder_result.inv_perm);
    global_atomic_grid_sizes_vec = std::move(reorder_result.global_atomic_grid_sizes);

    // Gather and reorder raw_grid_weights using the same MPI layout
    GAUXC_MPI_CODE(
      if (world_size > 1) {
        int local_npts = (int)total_npts;
        std::vector<double> recv_raw(world_rank == 0 ? total_npts_sum : 0);
        MPI_Gatherv(raw_grid_weights.data(), local_npts, MPI_DOUBLE,
                    recv_raw.data(), recvcounts.data(), displs.data(),
                    MPI_DOUBLE, 0, rt.comm());
        if (world_rank == 0) {
          raw_grid_weights = std::move(recv_raw);
          std::vector<int64_t> perm(total_npts_sum);
          for (int64_t j = 0; j < total_npts_sum; j++) perm[atom_reorder_inv_perm[j]] = j;
          std::vector<double> tmp(total_npts_sum);
          for (int64_t i = 0; i < total_npts_sum; i++) tmp[perm[i]] = raw_grid_weights[i];
          raw_grid_weights = std::move(tmp);
        }
      }
    )

    if (world_rank == 0) {
      int64_t max_grid_size = *std::max_element(
        global_atomic_grid_sizes_vec.begin(), global_atomic_grid_sizes_vec.end());
      auto options = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);
      features_dict = prepare_onedft_features(
        natoms, total_npts_sum, ndm, options, feature_keys, den_eval.data(),
        dden_eval.data(), tau.data(), grid_coords.data(), grid_weights.data(),
        host_coords.data(),
        global_atomic_grid_sizes_vec, max_grid_size,
        raw_grid_weights.data()
      );
    }
  }
  if (world_rank == 0) {
    auto exc_on_grid = get_exc(exc_func, features_dict);
    auto exc = (exc_on_grid * features_dict.at(feat_map.at(ONEDFT_FEATURE::WEIGHTS))).sum();
    // if do_vxc
    exc.backward();
    EXC[0] = exc.item<double>();
  } else {
    EXC[0] = 0.0;
  }

  if ( world_size == 1 ) {
    double* den_grad, * dden_grad, * tau_grad;
    den_grad = features_dict.at(feat_map.at(ONEDFT_FEATURE::DEN)).grad().data_ptr<double>();
    if (features_dict.find(feat_map.at(ONEDFT_FEATURE::DDEN)) != features_dict.end()) {
      dden_grad = features_dict.at(feat_map.at(ONEDFT_FEATURE::DDEN)).grad().data_ptr<double>();
    } else {
      dden_grad = nullptr;
    }
    if (features_dict.find(feat_map.at(ONEDFT_FEATURE::TAU)) != features_dict.end()){
      tau_grad = features_dict.at(feat_map.at(ONEDFT_FEATURE::TAU)).grad().data_ptr<double>();
    } else {
      tau_grad = nullptr;
    }
    device_data_ptr->send_static_data_onedft_results( total_npts, ndm, EXC,
      den_grad, dden_grad, tau_grad );
  } else { 
    total_npts = mpi_scatter_onedft_outputs(features_dict, rt.comm_rank(), rt.comm_size(),
                                              recvcounts, displs, atom_reorder_inv_perm,
                                              den_eval, dden_eval, tau);
    device_data_ptr->send_static_data_onedft_results( total_npts, ndm, EXC,
      den_eval.data(), dden_eval.data(), tau.data());
  }

  this->timer_.time_op("XCIntegrator.LocalWork_PostOneDFT", [&](){
    post_onedft_local_work_( basis, Ps, ldps, Pz, ldpz, nullptr, 0, nullptr, 0,
      tasks.begin(), tasks.end(), *device_data_ptr, enabled_terms );
  });

  rt.device_backend()->master_queue_synchronize();

  value_type N_EL;

  this->timer_.time_op("XCIntegrator.DeviceToHostCopy_EXC_VXC",[&](){
    device_data_ptr->retrieve_exc_vxc_integrands( EXC, &N_EL, VXCs, ldvxcs, VXCz, ldvxcz, 
      nullptr, 0, nullptr, 0 );
  });

  
  GAUXC_MPI_CODE(
    this->timer_.time_op("XCIntegrator.ImbalanceWait_PostOneDFT",[&](){
      MPI_Barrier(rt.comm());
    });
  )

  this->timer_.time_op("XCIntegrator.Allreduce_OneDFT", [&](){
    this->reduction_driver_->allreduce_inplace( VXCs, nbf*nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( VXCz, nbf*nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( EXC, 1,       ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_EL, 1,       ReductionOp::Sum );
  });
  // std::cout << "exc: " << EXC[0] << std::endl;
} // eval_exc_vxc_onedft_

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
pre_onedft_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            const value_type* Pz, int64_t ldpz,
                            const value_type* Py, int64_t ldpy,
                            const value_type* Px, int64_t ldpx,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data, const integrator_term_tracker enabled_terms ) {

  bool is_mgga = enabled_terms.xc_approx == integrator_xc_approx::MGGA_TAU;
  bool is_gga = enabled_terms.xc_approx == integrator_xc_approx::GGA;
  bool is_lda = enabled_terms.xc_approx == integrator_xc_approx::LDA;
  bool is_rks = enabled_terms.ks_scheme == RKS;
  bool is_uks = enabled_terms.ks_scheme == UKS;
        
  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
  const auto& mol   = this->load_balancer_->molecule();
  const auto natoms = mol.natoms();
  BasisSetMap basis_map(basis,mol);
  device_data.populate_submat_maps( basis.nbf(), task_begin, task_end, basis_map );

  size_t total_npts = std::accumulate( task_begin, task_end, 0ul,
    [](const auto& a, const auto& b) { return a + b.npts; } );

  auto task_it = task_begin;
  size_t offset = 0;
  while( task_it != task_end ) {

    // Determine next task batch, send relevant data to device (EXC VXC only)
    task_it = 
      device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );


    const bool need_lapl = false;
    // Evaluate collocation
    if( is_mgga ) {
      if(need_lapl) lwd->eval_collocation_laplacian( &device_data );
      else          lwd->eval_collocation_gradient( &device_data );
    }
    else if( is_gga ) lwd->eval_collocation_gradient( &device_data );
    else                     lwd->eval_collocation( &device_data );
      
    const double xmat_fac = is_rks ? 2.0 : 1.0;
    const bool need_xmat_grad = is_mgga;

    // Evaluate X matrix and V vars
    auto do_xmat_vvar = [&](density_id den_id) {
      lwd->eval_xmat( xmat_fac, &device_data, need_xmat_grad, den_id );
      if(is_lda)      lwd->eval_vvars_lda( &device_data, den_id );
      else if(is_gga) lwd->eval_vvars_gga( &device_data, den_id ); 
      else                   lwd->eval_vvars_mgga( &device_data, den_id, need_lapl );
    };

    do_xmat_vvar(DEN_S);
    if (not is_rks) {
      do_xmat_vvar(DEN_Z);
      if (not is_uks) {
        do_xmat_vvar(DEN_Y);
        do_xmat_vvar(DEN_X);
      }
    }

    // Evaluate U variables
    if( is_mgga )      lwd->eval_uvars_mgga( &device_data, enabled_terms.ks_scheme, need_lapl );
    else if( is_gga )  lwd->eval_uvars_gga ( &device_data, enabled_terms.ks_scheme );
    else                      lwd->eval_uvars_lda ( &device_data, enabled_terms.ks_scheme );
    
    if (is_mgga or is_gga)
      lwd->sz_to_ab_onedft( &device_data, offset);
    offset = save_static_data_onedft_features( &device_data, enabled_terms, offset );
  }
  if (offset != total_npts) {
    GAUXC_GENERIC_EXCEPTION("eval_exc_vxc_onedft: Offset does not match total points");
  }
}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
post_onedft_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            const value_type* Pz, int64_t ldpz,
                            const value_type* Py, int64_t ldpy,
                            const value_type* Px, int64_t ldpx,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data, const integrator_term_tracker enabled_terms ) {

  bool is_mgga = enabled_terms.xc_approx == integrator_xc_approx::MGGA_TAU;
  bool is_gga = enabled_terms.xc_approx == integrator_xc_approx::GGA;
  bool is_lda = enabled_terms.xc_approx == integrator_xc_approx::LDA;
  bool is_rks = enabled_terms.ks_scheme == RKS;
  bool is_uks = enabled_terms.ks_scheme == UKS;

  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
  const auto& mol   = this->load_balancer_->molecule();
  BasisSetMap basis_map(basis,mol);
  device_data.populate_submat_maps( basis.nbf(), task_begin, task_end, basis_map );
  const auto nbf     = basis.nbf();
  const auto nshells = basis.nshells();
  size_t total_npts = std::accumulate( task_begin, task_end, 0ul,
    [](const auto& a, const auto& b) { return a + b.npts; } );

  // Check that Partition Weights have been calculated
  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored ) {
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Been Modified"); 
  }

  auto task_it = task_begin;
  size_t offset = 0;

  while( task_it != task_end ) {
    task_it = device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );
    // std::cout << offset << " offset: " << offset << std::endl;
    offset = send_buffer_onedft_outputs( &device_data, enabled_terms, offset );

    // Evaluate collocation
    if( is_mgga )     lwd->eval_collocation_gradient( &device_data );
    else if( is_gga ) lwd->eval_collocation_gradient( &device_data );
    else              lwd->eval_collocation( &device_data );
    auto do_zmat_vxc = [&](density_id den_id) {
      if( is_mgga ) {
        lwd->eval_zmat_onedft( &device_data, enabled_terms, den_id);
        lwd->eval_mmat_mgga_vxc( &device_data, enabled_terms.ks_scheme, false /*need_lapl*/, den_id);
      } else 
        lwd->eval_zmat_onedft( &device_data, enabled_terms, den_id);
      lwd->inc_vxc( &device_data, den_id, is_mgga );
   };
   do_zmat_vxc(DEN_S);
   if(not is_rks) {
     do_zmat_vxc(DEN_Z);
    }
   } // Loop over batches of batches 

  if (offset != total_npts) {
    GAUXC_GENERIC_EXCEPTION("eval_exc_vxc_onedft: Offset does not match total points");
  }

  // Symmetrize VXC in device memory
  lwd->symmetrize_vxc( &device_data, DEN_S );
  if (not is_rks) {
    lwd->symmetrize_vxc( &device_data, DEN_Z );
    if (not is_uks) {
      lwd->symmetrize_vxc( &device_data, DEN_Y );
      lwd->symmetrize_vxc( &device_data, DEN_X );
    }
  }
} // onedft_local_work_

void save_static_data_onedft_outputs(const at::Tensor EXC, const FeatureDict& features_dict, XCDeviceData* _data) {

  auto* data = dynamic_cast<Scheme1DataBase*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();
  auto backend = dynamic_cast<CUDABackend*>(data->device_backend_);
  auto static_stack = data->static_stack;
  size_t total_npts = data->global_dims.total_npts;

  // std::cout << "save_static_data_onedft_outputs EXC: " << *EXC << std::endl;
  backend->copy_async(1, EXC.data_ptr<double>(), static_stack.exc_device, "Copy OneDFT EXC");

  // std::cout << "den_grad: " << features_dict.at(feat_map.at(ONEDFT_FEATURE::DEN)).grad() << std::endl;
  // copy exc gradient to static stack
  backend->copy_async(2 * total_npts, features_dict.at(feat_map.at(ONEDFT_FEATURE::DEN)).grad().data_ptr<double>(), 
                      static_stack.den_grad_device, 
                      "Copy OneDFT den_grad_device");

  if (features_dict.find(feat_map.at(ONEDFT_FEATURE::DDEN)) != features_dict.end()){
    backend->copy_async(2 * 3 * total_npts, features_dict.at(feat_map.at(ONEDFT_FEATURE::DDEN)).grad().data_ptr<double>(), 
                        static_stack.dden_grad_device, 
                        "Copy OneDFT dden_grad_device");
    // std::cout << "dden_grad: " << features_dict.at(feat_map.at(ONEDFT_FEATURE::DDEN)).grad() << std::endl;
  }
  if (features_dict.find(feat_map.at(ONEDFT_FEATURE::TAU)) != features_dict.end()){
    backend->copy_async(2 * total_npts, features_dict.at(feat_map.at(ONEDFT_FEATURE::TAU)).grad().data_ptr<double>(), 
                        static_stack.tau_grad_device,
                        "Copy OneDFT tau_grad_device");
    // std::cout << "tau_grad: " << features_dict.at(feat_map.at(ONEDFT_FEATURE::TAU)).grad() << std::endl;
  }
}

size_t send_buffer_onedft_outputs(XCDeviceData* _data, const integrator_term_tracker enabled_terms, size_t offset) {
  auto* data = dynamic_cast<Scheme1DataBase*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();
  auto backend = dynamic_cast<CUDABackend*>(data->device_backend_);

  size_t npoints = data->total_npts_task_batch;
  size_t total_npts = data->global_dims.total_npts;
  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;

  const bool is_UKS  = data->allocated_terms.ks_scheme == UKS;
  size_t ndm = is_UKS ? 2 : 1;

  double* den_grad_a   = static_stack.den_grad_device + offset;
  double* den_grad_b   = static_stack.den_grad_device + total_npts + offset;

  double* dden_x_grad_a = static_stack.dden_grad_device + offset;
  double* dden_y_grad_a = static_stack.dden_grad_device + total_npts + offset;
  double* dden_z_grad_a = static_stack.dden_grad_device + total_npts*2 + offset;

  double* dden_x_grad_b = static_stack.dden_grad_device + total_npts*3 + offset;
  double* dden_y_grad_b = static_stack.dden_grad_device + total_npts*4 + offset;
  double* dden_z_grad_b = static_stack.dden_grad_device + total_npts*5 + offset;

  double* tau_a        = static_stack.tau_grad_device + offset;
  double* tau_b        = static_stack.tau_grad_device + total_npts + offset;

  backend->copy_async_2d(
    1, npoints, den_grad_a, 1, base_stack.vrho_pos_eval_device, 1, "Copy vrho_pos_eval_device");
  backend->copy_async_2d(
    1, npoints, den_grad_b, 1, base_stack.vrho_neg_eval_device, 1, "Copy vrho_neg_eval_device");
  
  if (dden_x_grad_a != nullptr && base_stack.gamma_pp_eval_device != nullptr) {
    backend->copy_async_2d(1, npoints, dden_x_grad_a, 1, base_stack.gamma_pp_eval_device, 1, "Copy dden_x_eval_a to gamma_pp_eval_device");
    backend->copy_async_2d(1, npoints, dden_x_grad_b, 1, base_stack.vgamma_pp_eval_device, 1, "Copy dden_x_eval_b to vgamma_pp_eval_device");
    backend->copy_async_2d(1, npoints, dden_y_grad_a, 1, base_stack.gamma_pm_eval_device, 1, "Copy dden_y_eval_a to gamma_pm_eval_device");
    backend->copy_async_2d(1, npoints, dden_y_grad_b, 1, base_stack.vgamma_pm_eval_device, 1, "Copy dden_y_eval_b to vgamma_pm_eval_device");
    backend->copy_async_2d(1, npoints, dden_z_grad_a, 1, base_stack.gamma_mm_eval_device, 1, "Copy dden_z_eval_a to gamma_mm_eval_device");
    backend->copy_async_2d(1, npoints, dden_z_grad_b, 1, base_stack.vgamma_mm_eval_device, 1, "Copy dden_z_eval_b to vgamma_mm_eval_device");
  }
  if ( tau_a != nullptr && base_stack.vtau_pos_eval_device != nullptr ) {
    backend->copy_async_2d(
      1, npoints, tau_a, 1, base_stack.vtau_pos_eval_device, 1, "Copy vtau_pos_eval_device");
    backend->copy_async_2d(
      1, npoints, tau_b, 1, base_stack.vtau_neg_eval_device, 1, "Copy vtau_neg_eval_device");
  }

  // concate den_eval_a and den_eval_b in device memory
  backend->master_queue_synchronize(); 
  return offset + npoints;
}

size_t save_static_data_onedft_features(XCDeviceData* _data, const integrator_term_tracker enabled_terms, size_t offset) {
  auto* data = dynamic_cast<Scheme1DataBase*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();
  auto backend = dynamic_cast<CUDABackend*>(data->device_backend_);

  size_t npoints = data->total_npts_task_batch;
  size_t total_npts = data->global_dims.total_npts;

  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;

  const bool is_UKS  = data->allocated_terms.ks_scheme == UKS;
  size_t ndm = is_UKS ? 2 : 1;

  double* grid_weights = static_stack.grid_weights_device + offset;
  double* grid_coords  = static_stack.grid_coords_device + offset * 3;

  double* den_eval_a   = static_stack.den_eval_device + offset;
  double* den_eval_b   = static_stack.den_eval_device + total_npts + offset;

  double* tau_a        = static_stack.tau_device + offset;
  double* tau_b        = static_stack.tau_device + total_npts + offset;
  
  backend->copy_async_2d(1, npoints, base_stack.weights_device, 1, grid_weights, 1, "Copy grid_weights");

  backend->copy_async_2d(
    1, npoints, base_stack.points_x_device, 1, grid_coords, 3, "Copy grid_coords x");
  backend->copy_async_2d(
    1, npoints, base_stack.points_y_device, 1, grid_coords + 1, 3, "Copy grid_coords y");
  backend->copy_async_2d(
    1, npoints, base_stack.points_z_device, 1, grid_coords + 2, 3, "Copy grid_coords z");

  backend->copy_async_2d(
    1, npoints, base_stack.den_s_eval_device, 1, den_eval_a, 1, "Copy den_eval_a");
  backend->copy_async_2d(
    1, npoints, base_stack.den_z_eval_device, 1, den_eval_b, 1, "Copy den_eval_b");

  if ( base_stack.tau_s_eval_device != nullptr ) {
    backend->copy_async_2d(
      1, npoints, base_stack.tau_s_eval_device, 1, tau_a, 1, "Copy tau_a");
    backend->copy_async_2d(
      1, npoints, base_stack.tau_z_eval_device, 1, tau_b, 1, "Copy tau_b");
  }

  // concate den_eval_a and den_eval_b in device memory
  backend->master_queue_synchronize(); 
  return offset + npoints;
}

FeatureDict prepare_onedft_features( const size_t natoms, const size_t total_npts, const size_t ndm,
  const at::TensorOptions options, const std::vector<std::string> feature_keys,
  double* den_eval, double* dden_eval, double* tau, double* grid_coords, 
  double* grid_weights, double* coords,
  const std::vector<int64_t>& atomic_grid_sizes,
  int64_t max_grid_size,
  double* raw_grid_weights ) {
  auto device = torch::Device(torch::kCUDA, 0);
  FeatureDict featmap;
  for (const auto& key : feature_keys) {
    auto enum_key = reverse_feat_map.at(key);
    switch (enum_key) {
    case ONEDFT_FEATURE::DEN: {
      auto flat_tensor = torch::from_blob(den_eval, {ndm * total_npts}, options);
      auto tensor = flat_tensor.view({ndm, total_npts}).to(device).requires_grad_(true);
      featmap.insert(key, tensor);
      break;
    }
    case ONEDFT_FEATURE::DDEN: {
      auto flat_tensor = torch::from_blob(dden_eval, {ndm * 3 * total_npts}, options);
      auto tensor = flat_tensor.view({ndm, 3, total_npts}).to(device).requires_grad_(true);
      featmap.insert(key, tensor);
      break;
    }
    case ONEDFT_FEATURE::TAU: {
      auto flat_tensor = torch::from_blob(tau, {ndm * total_npts}, options);
      auto tensor = flat_tensor.view({ndm, total_npts}).to(device).requires_grad_(true);
      featmap.insert(key, tensor);
      break;
    }
    case ONEDFT_FEATURE::POINTS: {
      auto flat_tensor = torch::from_blob(grid_coords, {total_npts * 3}, options);
      auto tensor = flat_tensor.view({total_npts, 3}).to(device);
      featmap.insert(key, tensor);
      break;
    }
    case ONEDFT_FEATURE::WEIGHTS: {
      auto flat_tensor = torch::from_blob(grid_weights, {total_npts}, options);
      auto tensor = flat_tensor.view({total_npts}).to(device);
      featmap.insert(key, tensor);
      break;
    }
    case ONEDFT_FEATURE::COORDS: {
      auto flat_tensor = torch::from_blob(coords, {natoms * 3}, options);
      auto tensor = flat_tensor.view({natoms, 3}).to(device);
      featmap.insert(key, tensor);
      break;
    }
    case ONEDFT_FEATURE::ATOMIC_GRID_WEIGHTS: {
      // Use raw (pre-partition) quadrature weights if available
      double* w_ptr = raw_grid_weights ? raw_grid_weights : grid_weights;
      auto w_opts = raw_grid_weights
        ? torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU)
        : options;
      auto flat_tensor = torch::from_blob(w_ptr, {total_npts}, w_opts);
      auto tensor = flat_tensor.view({total_npts}).to(device);
      featmap.insert(key, tensor);
      break;
    }
    case ONEDFT_FEATURE::ATOMIC_GRID_SIZES: {
      auto sizes_options = torch::TensorOptions().dtype(torch::kInt64).device(torch::kCPU);
      auto tensor = torch::from_blob(const_cast<int64_t*>(atomic_grid_sizes.data()),
                                     {static_cast<int64_t>(natoms)}, sizes_options).clone().to(device);
      featmap.insert(key, tensor);
      break;
    }
    case ONEDFT_FEATURE::ATOMIC_GRID_SIZE_BOUND_SHAPE: {
      auto sizes_options = torch::TensorOptions().dtype(torch::kInt64).device(torch::kCPU);
      auto tensor = torch::zeros({max_grid_size, 0}, sizes_options).to(device);
      featmap.insert(key, tensor);
      break;
    }
    default:
      GAUXC_GENERIC_EXCEPTION("Feature Key Not Implemented: " + key);
    }
  }
  return featmap;
}
} // namespace GauXC::detail (prepare_onedft_features)

// ============================================================================
// OneDFT EXC Gradient — Device Implementation
// ============================================================================

namespace GauXC::detail {

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
eval_exc_grad_onedft_( int64_t m, int64_t n, const value_type* Ps, int64_t ldps,
                       const value_type* Pz, int64_t ldpz, value_type* EXC_GRAD,
                       const IntegratorSettingsXC& settings ) {

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

  auto& tasks = this->load_balancer_->get_tasks();
  // Sort tasks by atom so grid points are contiguous per atom
  std::stable_sort(tasks.begin(), tasks.end(),
    [](const auto& a, const auto& b) { return a.iParent < b.iParent; });
  size_t total_npts = std::accumulate( tasks.begin(), tasks.end(), 0ul,
    [](const auto& a, const auto& b) { return a + b.npts; } );

  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get());
  auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
  auto device_data_ptr = lwd->create_device_data(rt);

  int32_t world_rank = rt.comm_rank();
  int32_t world_size = rt.comm_size();

  // Load model
  OneDFTSettings onedft_settings;
  if( auto* tmp = dynamic_cast<const OneDFTSettings*>(&settings) ) {
    onedft_settings = *tmp;
  }
  const auto model_path = onedft_settings.model;
  torch::DeviceType torch_device = torch::kCPU;
  auto [exc_func, feature_keys] = load_model(model_path, torch_device);

  // Determine feature requirements
  bool is_gga = false;
  bool is_mgga = false;
  for (const auto& key : feature_keys) {
    if ( not valueExists(key) ) GAUXC_GENERIC_EXCEPTION("Feature Key Required Not Implemented: " + key);
    if (key == feat_map.at(ONEDFT_FEATURE::TAU)) is_mgga = true;
    if (key == feat_map.at(ONEDFT_FEATURE::DDEN)) is_gga = true;
  }
  if (is_mgga) is_gga = false;

  const auto& mol   = this->load_balancer_->molecule();
  const auto natoms = mol.natoms();
  const auto nshells = basis.nshells();

  // Phase 1: Pre-work — compute density features on device
  integrator_term_tracker enabled_terms;
  enabled_terms.exc_vxc = true;
  enabled_terms.onedft = true;
  enabled_terms.exc_grad = true;  // Also allocate exc_grad_device
  if (is_uks) enabled_terms.ks_scheme = UKS;
  else        enabled_terms.ks_scheme = RKS;

  if (is_mgga)     enabled_terms.xc_approx = integrator_xc_approx::MGGA_TAU;
  else if (is_gga) enabled_terms.xc_approx = integrator_xc_approx::GGA;
  else             enabled_terms.xc_approx = integrator_xc_approx::LDA;

  device_data_ptr->reset_allocations();
  device_data_ptr->allocate_static_data_onedft( nbf, nshells, natoms, total_npts, enabled_terms );
  device_data_ptr->send_static_data_onedft( mol, Ps, ldps, Pz, ldpz, nullptr, 0, nullptr, 0, basis );
  device_data_ptr->zero_exc_grad_integrands();

  // Pre-OneDFT density computation on device
  integrator_term_tracker pre_terms = enabled_terms;
  pre_terms.exc_grad = false;  // pre-work doesn't need gradient buffers
  this->timer_.time_op("XCIntegrator.LocalWork_PreOneDFT", [&](){
    pre_onedft_local_work_( basis, Ps, ldps, Pz, ldpz, nullptr, 0, nullptr, 0,
      tasks.begin(), tasks.end(), *device_data_ptr, pre_terms );
  });

  // Collect raw quadrature weights and atomic grid sizes (for model features)
  std::vector<double> raw_grid_weights;
  raw_grid_weights.reserve(total_npts);
  for (const auto& task : tasks) {
    raw_grid_weights.insert(raw_grid_weights.end(), task.raw_weights.begin(), task.raw_weights.end());
  }
  std::vector<int64_t> atomic_grid_sizes_vec(natoms, 0);
  for (const auto& task : tasks) {
    if (task.iParent >= 0 && task.iParent < (int)natoms) {
      atomic_grid_sizes_vec[task.iParent] += task.npts;
    }
  }

  size_t ndm = 2; // UKS

  // Phase 2: Retrieve features and run model
  // For single rank: features stay on device → just build tensors from device pointers
  // For multi rank: retrieve to host, gather, reorder
  FeatureDict features_dict;
  std::vector<int64_t> atom_reorder_inv_perm;
  std::vector<double> den_eval, dden_eval, tau;
  std::vector<int> displs(world_size), recvcounts(world_size);

  if (world_size == 1) {
    // Keep on device for model inference
    int64_t max_grid_size = *std::max_element(
      atomic_grid_sizes_vec.begin(), atomic_grid_sizes_vec.end());
    // For gradient, we need CPU tensors since we need requires_grad on points/coords
    // Retrieve features to host
    std::vector<double> grid_weights(total_npts), grid_coords(total_npts * 3);
    den_eval.resize(total_npts * ndm);
    if (is_gga || is_mgga) dden_eval.resize(total_npts * ndm * 3);
    if (is_mgga) tau.resize(total_npts * ndm);

    device_data_ptr->retrieve_onedft_features( total_npts, 2, den_eval.data(),
      (is_gga || is_mgga) ? dden_eval.data() : nullptr,
      is_mgga ? tau.data() : nullptr,
      grid_coords.data(), grid_weights.data() );
    rt.device_backend()->master_queue_synchronize();

    std::vector<double> host_coords(natoms * 3);
    for (size_t i = 0; i < (size_t)natoms; i++) {
      host_coords[3*i]   = mol[i].x;
      host_coords[3*i+1] = mol[i].y;
      host_coords[3*i+2] = mol[i].z;
    }

    auto options = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);
    features_dict = prepare_onedft_features(
      natoms, total_npts, ndm, options, feature_keys, den_eval.data(),
      dden_eval.data(), tau.data(), grid_coords.data(), grid_weights.data(),
      host_coords.data(), atomic_grid_sizes_vec, max_grid_size,
      raw_grid_weights.data() );
  } else {
    // Multi-rank: retrieve features, gather, reorder
    std::vector<double> grid_weights(total_npts), grid_coords(total_npts * 3);
    den_eval.resize(total_npts * ndm);
    if (is_gga || is_mgga) dden_eval.resize(total_npts * ndm * 3);
    if (is_mgga) tau.resize(total_npts * ndm);

    device_data_ptr->retrieve_onedft_features( total_npts, 2, den_eval.data(),
      (is_gga || is_mgga) ? dden_eval.data() : nullptr,
      is_mgga ? tau.data() : nullptr,
      grid_coords.data(), grid_weights.data() );
    rt.device_backend()->master_queue_synchronize();

    std::vector<double> host_coords(natoms * 3);
    for (size_t i = 0; i < (size_t)natoms; i++) {
      host_coords[3*i]   = mol[i].x;
      host_coords[3*i+1] = mol[i].y;
      host_coords[3*i+2] = mol[i].z;
    }

    auto reorder_result = mpi_gather_and_reorder_gpu(
      den_eval, dden_eval, tau, grid_coords, grid_weights,
      atomic_grid_sizes_vec, total_npts, natoms, rt, recvcounts, displs);
    atom_reorder_inv_perm = std::move(reorder_result.inv_perm);

    if (world_rank == 0) {
      int64_t max_grid_size = *std::max_element(
        reorder_result.global_atomic_grid_sizes.begin(),
        reorder_result.global_atomic_grid_sizes.end());
      auto options = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);
      features_dict = prepare_onedft_features(
        natoms, reorder_result.total_npts, ndm, options, feature_keys,
        den_eval.data(), dden_eval.data(), tau.data(),
        grid_coords.data(), grid_weights.data(), host_coords.data(),
        reorder_result.global_atomic_grid_sizes, max_grid_size,
        raw_grid_weights.data() );
    }
  }

  // Phase 3: Forward + backward with requires_grad on POINTS and COORDS
  std::vector<double> eps_on_grid_global;
  std::vector<double> points_grad_global;
  std::vector<double> coords_grad_global;

  if (world_rank == 0) {
    // Enable requires_grad on points and coords
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

    // Extract eps_on_grid for weight derivative
    int total_npts_model = exc_on_grid.size(0);
    at::Tensor eps_cpu = exc_on_grid.detach().cpu().contiguous();
    eps_on_grid_global.resize(total_npts_model);
    std::memcpy(eps_on_grid_global.data(), eps_cpu.data_ptr<double>(),
                total_npts_model * sizeof(double));

    // Extract points.grad()
    if (features_dict.find(feat_map.at(ONEDFT_FEATURE::POINTS)) != features_dict.end()) {
      auto pg = features_dict.at(feat_map.at(ONEDFT_FEATURE::POINTS)).grad();
      if (pg.defined()) {
        at::Tensor pg_cpu = pg.cpu().contiguous();
        points_grad_global.resize(total_npts_model * 3);
        std::memcpy(points_grad_global.data(), pg_cpu.data_ptr<double>(),
                    total_npts_model * 3 * sizeof(double));
      }
    }

    // Extract coords.grad()
    if (features_dict.find(feat_map.at(ONEDFT_FEATURE::COORDS)) != features_dict.end()) {
      auto cg = features_dict.at(feat_map.at(ONEDFT_FEATURE::COORDS)).grad();
      if (cg.defined()) {
        at::Tensor cg_cpu = cg.cpu().contiguous();
        coords_grad_global.resize(natoms * 3);
        std::memcpy(coords_grad_global.data(), cg_cpu.data_ptr<double>(),
                    natoms * 3 * sizeof(double));
      }
    }

    // Reorder eps_on_grid and points_grad from atom-order back to rank-order
    if (!atom_reorder_inv_perm.empty()) {
      std::vector<double> tmp(total_npts_model);
      for (int64_t i = 0; i < total_npts_model; i++) {
        tmp[atom_reorder_inv_perm[i]] = eps_on_grid_global[i];
      }
      eps_on_grid_global = std::move(tmp);

      if (!points_grad_global.empty()) {
        std::vector<double> tmp3(total_npts_model * 3);
        for (int64_t i = 0; i < total_npts_model; i++) {
          int64_t j = atom_reorder_inv_perm[i];
          tmp3[j*3+0] = points_grad_global[i*3+0];
          tmp3[j*3+1] = points_grad_global[i*3+1];
          tmp3[j*3+2] = points_grad_global[i*3+2];
        }
        points_grad_global = std::move(tmp3);
      }
    }
  }

  // Phase 4: Send OneDFT Vxc outputs back to device
  if (world_size == 1) {
    double* den_grad  = features_dict.at(feat_map.at(ONEDFT_FEATURE::DEN)).grad().data_ptr<double>();
    double* dden_grad = nullptr;
    double* tau_grad  = nullptr;
    if (features_dict.find(feat_map.at(ONEDFT_FEATURE::DDEN)) != features_dict.end()) {
      dden_grad = features_dict.at(feat_map.at(ONEDFT_FEATURE::DDEN)).grad().data_ptr<double>();
    }
    if (features_dict.find(feat_map.at(ONEDFT_FEATURE::TAU)) != features_dict.end()) {
      tau_grad = features_dict.at(feat_map.at(ONEDFT_FEATURE::TAU)).grad().data_ptr<double>();
    }
    double exc_val = 0.0; // Not needed for gradient, but API requires it
    device_data_ptr->send_static_data_onedft_results( total_npts, ndm, &exc_val,
      den_grad, dden_grad, tau_grad );
  } else {
    total_npts = mpi_scatter_onedft_outputs(features_dict, rt.comm_rank(), rt.comm_size(),
                                              recvcounts, displs, atom_reorder_inv_perm,
                                              den_eval, dden_eval, tau);
    double exc_val = 0.0;
    device_data_ptr->send_static_data_onedft_results( total_npts, ndm, &exc_val,
      den_eval.data(), dden_eval.data(), tau.data());
  }

  // Scatter eps_on_grid to local tasks
  std::vector<double> eps_on_grid_local;
  if (world_size == 1) {
    eps_on_grid_local = std::move(eps_on_grid_global);
  } else {
    GAUXC_GENERIC_EXCEPTION("OneDFT gradient with MPI not yet implemented");
  }

  // Zero out EXC_GRAD on host
  for (int i = 0; i < 3*natoms; ++i) EXC_GRAD[i] = 0.0;

  // Phase 5: Pulay gradient + weight derivative on device
  // Use the standard gradient flow: for each batch, load OneDFT Vxc,
  // evaluate collocation (hessian), X-matrix, transform, and inc_exc_grad

  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored ) {
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Been Modified");
  }
  XCWeightAlg& weight_alg = lb_state.weight_alg;

  BasisSetMap basis_map(basis, mol);
  const auto& meta = this->load_balancer_->molmeta();
  device_data_ptr->populate_submat_maps( nbf, tasks.begin(), tasks.end(), basis_map );

  // Sort tasks by workload for load balancing
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.bfn_screening.nbe) > (b.points.size() * b.bfn_screening.nbe);
  };

  // First: distribute eps_on_grid to tasks (in iParent order, before sorting by workload)
  // Tasks are currently in iParent order (from stable_sort above)
  {
    size_t offset = 0;
    for (auto& task : tasks) {
      int64_t npts = task.points.size();
      task.feat.eps.resize(npts);
      std::copy(eps_on_grid_local.data() + offset,
                eps_on_grid_local.data() + offset + npts,
                task.feat.eps.begin());
      offset += npts;
    }
  }

  // Now sort by workload
  std::sort( tasks.begin(), tasks.end(), task_comparator );

  // Build concatenated eps_w array (eps_on_grid * weights) in task order after sorting
  // This will be used for the weight derivative
  std::vector<double> eps_w_all;
  eps_w_all.reserve(total_npts);
  for (const auto& task : tasks) {
    for (size_t ipt = 0; ipt < task.points.size(); ++ipt) {
      eps_w_all.push_back(task.feat.eps[ipt] * task.weights[ipt]);
    }
  }

  // Set up gradient computation terms
  integrator_term_tracker grad_terms;
  grad_terms.exc_grad = true;
  grad_terms.weights  = true;
  grad_terms.onedft   = true;
  grad_terms.ks_scheme = UKS;
  if (is_mgga)     grad_terms.xc_approx = integrator_xc_approx::MGGA_TAU;
  else if (is_gga) grad_terms.xc_approx = integrator_xc_approx::GGA;
  else             grad_terms.xc_approx = integrator_xc_approx::LDA;

  // Weight derivative needs: RAB, coords, dist scratch
  device_data_ptr->allocate_static_data_weights( natoms );
  device_data_ptr->send_static_data_weights( mol, meta );

  this->timer_.time_op("XCIntegrator.LocalWork_OneDFTGrad", [&](){
    auto task_it = tasks.begin();
    size_t vxc_offset = 0;
    size_t eps_offset = 0;
    while( task_it != tasks.end() ) {

      auto batch_begin = task_it;

      // Generate batch buffers (allocates collocation, xmat, hessian, etc.)
      task_it = device_data_ptr->generate_buffers( grad_terms, basis_map,
                                                    task_it, tasks.end() );

      // Load OneDFT Vxc for this batch
      vxc_offset = send_buffer_onedft_outputs( device_data_ptr.get(), grad_terms, vxc_offset );

      // Evaluate collocation (hessian for GGA/MGGA, gradient for LDA)
      if (is_gga || is_mgga) lwd->eval_collocation_hessian( device_data_ptr.get() );
      else                   lwd->eval_collocation_gradient( device_data_ptr.get() );

      // Evaluate X-matrices and save for gradient
      const double xmat_fac = 1.0; // UKS: factor 1
      const bool need_xmat_grad = is_gga || is_mgga;
      auto do_xmat = [&](density_id den_id) {
        lwd->eval_xmat( xmat_fac, device_data_ptr.get(), need_xmat_grad, den_id );
        lwd->save_xmat( device_data_ptr.get(), need_xmat_grad, den_id );
      };
      do_xmat(DEN_S);
      do_xmat(DEN_Z);

      // For GGA/MGGA: transform OneDFT per-direction Vxc to standard format
      if (is_gga || is_mgga) {
        lwd->transform_onedft_vxc_for_grad( device_data_ptr.get() );
      }

      // Increment EXC gradient using standard kernels
      const bool with_weight_derivatives = true;
      if (is_mgga)     lwd->inc_exc_grad_mgga( device_data_ptr.get(), UKS, false, with_weight_derivatives );
      else if (is_gga) lwd->inc_exc_grad_gga( device_data_ptr.get(), UKS, with_weight_derivatives );
      else             lwd->inc_exc_grad_lda( device_data_ptr.get(), UKS, with_weight_derivatives );

      // Weight derivative: load eps*w values for this batch into eps_eval,
      // set den_s=1.0 and den_z=0.0 so kernel computes eps*(1.0) = eps*w
      {
        size_t batch_npts = 0;
        for (auto it = batch_begin; it != task_it; ++it)
          batch_npts += it->points.size();

        auto* data = dynamic_cast<Scheme1DataBase*>(device_data_ptr.get());
        auto* backend = dynamic_cast<CUDABackend*>(data->device_backend_);
        auto base_stack = data->base_stack;

        // Copy eps*w for this batch to eps_eval_device
        backend->copy_async( batch_npts, eps_w_all.data() + eps_offset,
                            base_stack.eps_eval_device, "Copy OneDFT eps*w" );
        // Set den_s = 1.0 and den_z = 0.0 so kernel's (den_s += den_z; eps *= den_s) is identity
        std::vector<double> ones(batch_npts, 1.0);
        backend->copy_async( batch_npts, ones.data(),
                            base_stack.den_s_eval_device, "Set den_s=1 for weight deriv" );
        if (base_stack.den_z_eval_device) {
          backend->set_zero( batch_npts, base_stack.den_z_eval_device, "Zero den_z for weight deriv" );
        }
        backend->master_queue_synchronize();

        eps_offset += batch_npts;
      }

      lwd->eval_weight_1st_deriv_contracted( device_data_ptr.get(), weight_alg );

    } // batch loop
  });

  rt.device_backend()->master_queue_synchronize();

  // Retrieve gradient from device
  double N_EL;
  device_data_ptr->retrieve_exc_grad_integrands( EXC_GRAD, &N_EL );
  rt.device_backend()->master_queue_synchronize();

  // Phase 6: Add autograd forces (points -> parent atoms, coords -> direct)
  if (!points_grad_global.empty() && world_rank == 0) {
    size_t pg_offset = 0;
    // Iterate in iParent order — re-sort tasks
    std::stable_sort(tasks.begin(), tasks.end(),
      [](const auto& a, const auto& b) { return a.iParent < b.iParent; });
    for (const auto& task : tasks) {
      int iParent = task.iParent;
      for (size_t ipt = 0; ipt < task.points.size(); ++ipt) {
        EXC_GRAD[3*iParent + 0] += points_grad_global[(pg_offset + ipt)*3 + 0];
        EXC_GRAD[3*iParent + 1] += points_grad_global[(pg_offset + ipt)*3 + 1];
        EXC_GRAD[3*iParent + 2] += points_grad_global[(pg_offset + ipt)*3 + 2];
      }
      pg_offset += task.points.size();
    }
  }

  if (!coords_grad_global.empty() && world_rank == 0) {
    for (int a = 0; a < natoms; ++a) {
      EXC_GRAD[3*a + 0] += coords_grad_global[3*a + 0];
      EXC_GRAD[3*a + 1] += coords_grad_global[3*a + 1];
      EXC_GRAD[3*a + 2] += coords_grad_global[3*a + 2];
    }
  }

  // Phase 7: Allreduce
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){
    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");
    this->reduction_driver_->allreduce_inplace( EXC_GRAD, 3*natoms, ReductionOp::Sum );
  });
}

} // namespace GauXC::detail
