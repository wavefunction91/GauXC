/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "scheme1_base.hpp"
#include "device/common/zmat_vxc.hpp"
#include "device/common/collocation_device.hpp"
#include "device/common/device_blas.hpp"
#include "device/common/xc_functional_eval_wrapper.hpp"
#include "device/common/uvvars.hpp"
#include "device/common/pack_submat.hpp"
#include "device/common/inc_potential.hpp"
#include "device/common/symmetrize_mat.hpp"
#include "device/common/increment_exc_grad.hpp"
#include "device/common/exx_ek_screening.hpp"


#include "device/common/shell_pair_to_task.hpp"
#ifdef GAUXC_ENABLE_CUDA
#include "device_specific/cuda_util.hpp"
#include "gpu/integral_data_types.hpp"
#include "gpu/obara_saika_integrals.hpp"
#include "gpu/chebyshev_boys_computation.hpp"

#define GAUXC_ENABLE_EXX
#endif

#ifdef GAUXC_ENABLE_EXX
namespace XGPU {
  void integral_0_task_batched(
    size_t ntasks, size_t nsubtask,
    int max_primpairs, size_t max_nsp,
    GauXC::XCDeviceTask*                device_tasks,
    const GauXC::TaskToShellPairDevice* task2sp,
    const std::array<int32_t, 4>*  subtasks,
    const int32_t* nprim_pairs_device,
    shell_pair** sp_ptr_device,
    double* sp_X_AB_device,
    double* sp_Y_AB_device,
    double* sp_Z_AB_device,
    double *boys_table,
    cudaStream_t stream);

  void integral_1_task_batched(
    size_t ntasks, size_t nsubtask,
    int max_primpairs, size_t max_nsp,
    GauXC::XCDeviceTask*                device_tasks,
    const GauXC::TaskToShellPairDevice* task2sp,
    const std::array<int32_t, 4>*  subtasks,
    const int32_t* nprim_pairs_device,
    shell_pair** sp_ptr_device,
    double* sp_X_AB_device,
    double* sp_Y_AB_device,
    double* sp_Z_AB_device,
    double *boys_table,
    cudaStream_t stream);

  void integral_2_task_batched(
    size_t ntasks, size_t nsubtask,
    int max_primpairs, size_t max_nsp,
    GauXC::XCDeviceTask*                device_tasks,
    const GauXC::TaskToShellPairDevice* task2sp,
    const std::array<int32_t, 4>*  subtasks,
    const int32_t* nprim_pairs_device,
    shell_pair** sp_ptr_device,
    double* sp_X_AB_device,
    double* sp_Y_AB_device,
    double* sp_Z_AB_device,
    double *boys_table,
    cudaStream_t stream);


  void integral_0_0_task_batched(
        size_t ntasks,
        size_t nsubtasks,
        int max_primpairs, size_t max_nsp,
        GauXC::XCDeviceTask*                device_tasks,
        const GauXC::TaskToShellPairDevice* task2sp,
        const std::array<int32_t, 4>*  subtasks,
        const int32_t* nprim_pairs_device,
        shell_pair** sp_ptr_device,
        double* sp_X_AB_device,
        double* sp_Y_AB_device,
        double* sp_Z_AB_device,
        double *boys_table,
        cudaStream_t stream);

  void integral_0_0_shell_batched(
        size_t nsp,
        size_t max_ntask,
        const GauXC::ShellPairToTaskDevice* sp2task,
        GauXC::XCDeviceTask*                device_tasks,
		    double *boys_table,
        cudaStream_t stream); 

  void integral_1_1_task_batched(
        size_t ntasks,
        size_t nsubtasks,
        int max_primpairs, size_t max_nsp,
        GauXC::XCDeviceTask*                device_tasks,
        const GauXC::TaskToShellPairDevice* task2sp,
        const std::array<int32_t, 4>*  subtasks,
        const int32_t* nprim_pairs_device,
        shell_pair** sp_ptr_device,
        double* sp_X_AB_device,
        double* sp_Y_AB_device,
        double* sp_Z_AB_device,
        double *boys_table,
        cudaStream_t stream);

  void integral_1_1_shell_batched(
        size_t nsp,
        size_t max_ntask,
        const GauXC::ShellPairToTaskDevice* sp2task,
        GauXC::XCDeviceTask*                device_tasks,
		    double *boys_table,
        cudaStream_t stream); 

  void integral_2_2_task_batched(
        size_t ntasks,
        size_t nsubtasks,
        int max_primpairs, size_t max_nsp,
        GauXC::XCDeviceTask*                device_tasks,
        const GauXC::TaskToShellPairDevice* task2sp,
        const std::array<int32_t, 4>*  subtasks,
        const int32_t* nprim_pairs_device,
        shell_pair** sp_ptr_device,
        double* sp_X_AB_device,
        double* sp_Y_AB_device,
        double* sp_Z_AB_device,
        double *boys_table,
        cudaStream_t stream);

  void integral_2_2_shell_batched(
        size_t nsp,
        size_t max_ntask,
        const GauXC::ShellPairToTaskDevice* sp2task,
        GauXC::XCDeviceTask*                device_tasks,
		    double *boys_table,
        cudaStream_t stream); 
        
  void integral_1_0_task_batched(
        bool swap,
        size_t ntasks,
        size_t nsubtasks,
        int max_primpairs, size_t max_nsp,
        GauXC::XCDeviceTask*                device_tasks,
        const GauXC::TaskToShellPairDevice* task2sp,
        const std::array<int32_t, 4>*  subtasks,
        const int32_t* nprim_pairs_device,
        shell_pair** sp_ptr_device,
        double* sp_X_AB_device,
        double* sp_Y_AB_device,
        double* sp_Z_AB_device,
        double *boys_table,
        cudaStream_t stream);

  void integral_1_0_shell_batched(
        bool swap,
        size_t nsp,
        size_t max_ntask,
        const GauXC::ShellPairToTaskDevice* sp2task,
        GauXC::XCDeviceTask*                device_tasks,
		    double *boys_table,
        cudaStream_t stream); 

  void integral_2_0_task_batched(
        bool swap,
        size_t ntasks,
        size_t nsubtasks,
        int max_primpairs, size_t max_nsp,
        GauXC::XCDeviceTask*                device_tasks,
        const GauXC::TaskToShellPairDevice* task2sp,
        const std::array<int32_t, 4>*  subtasks,
        const int32_t* nprim_pairs_device,
        shell_pair** sp_ptr_device,
        double* sp_X_AB_device,
        double* sp_Y_AB_device,
        double* sp_Z_AB_device,
        double *boys_table,
        cudaStream_t stream);

  void integral_2_0_shell_batched(
        bool swap,
        size_t nsp,
        size_t max_ntask,
        const GauXC::ShellPairToTaskDevice* sp2task,
        GauXC::XCDeviceTask*                device_tasks,
		    double *boys_table,
        cudaStream_t stream); 

  void integral_2_1_task_batched(
        bool swap,
        size_t ntasks,
        size_t nsubtasks,
        int max_primpairs, size_t max_nsp,
        GauXC::XCDeviceTask*                device_tasks,
        const GauXC::TaskToShellPairDevice* task2sp,
        const std::array<int32_t, 4>*  subtasks,
        const int32_t* nprim_pairs_device,
        shell_pair** sp_ptr_device,
        double* sp_X_AB_device,
        double* sp_Y_AB_device,
        double* sp_Z_AB_device,
        double *boys_table,
        cudaStream_t stream);

  void integral_2_1_shell_batched(
        bool swap,
        size_t nsp,
        size_t max_ntask,
        const GauXC::ShellPairToTaskDevice* sp2task,
        GauXC::XCDeviceTask*                device_tasks,
		    double *boys_table,
        cudaStream_t stream); 
}
#endif


namespace GauXC {

AoSScheme1Base::AoSScheme1Base() {
#ifdef GAUXC_ENABLE_EXX
  dev_boys_table = XGPU::boys_init();
#endif
}

AoSScheme1Base::~AoSScheme1Base() noexcept {
#ifdef GAUXC_ENABLE_EXX
  XGPU::boys_finalize(dev_boys_table);
#endif
}

void AoSScheme1Base::eval_zmat_lda_vxc_rks( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.bfn_screening.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  auto aos_stack     = data->aos_stack;
  zmat_lda_vxc( ntasks, nbe_max, npts_max, aos_stack.device_tasks,
    data->device_backend_->queue() );

}

void AoSScheme1Base::eval_zmat_gga_vxc_rks( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.bfn_screening.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  auto aos_stack     = data->aos_stack;
  zmat_gga_vxc( ntasks, nbe_max, npts_max, aos_stack.device_tasks,
    data->device_backend_->queue() );

}

void AoSScheme1Base::eval_zmat_lda_vxc_uks( XCDeviceData* ){

  GAUXC_GENERIC_EXCEPTION("UKS NOT YET IMPLEMENTED FOR DEVICE");

}

void AoSScheme1Base::eval_zmat_gga_vxc_uks( XCDeviceData* ){

  GAUXC_GENERIC_EXCEPTION("UKS NOT YET IMPLEMENTED FOR DEVICE");

}

void AoSScheme1Base::eval_collocation( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  size_t npts_max = 0, nshells_max = 0;
  for( auto& task : tasks ) {
    npts_max    = std::max( npts_max, task.npts );
    nshells_max = std::max( nshells_max, task.bfn_screening.nshells );
  }

  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  if( ! static_stack.shells_device )
    GAUXC_GENERIC_EXCEPTION("Shells not Allocated");
  if( ! aos_stack.device_tasks )
    GAUXC_GENERIC_EXCEPTION("Device Tasks not Allocated");

  eval_collocation_masked_combined( ntasks, npts_max, nshells_max,
    static_stack.shells_device, aos_stack.device_tasks, 
    data->device_backend_->queue() );

}

void AoSScheme1Base::eval_collocation_gradient( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

#ifdef GAUXC_ENABLE_HIP
  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  size_t npts_max = 0, nshells_max = 0;
  for( auto& task : tasks ) {
    npts_max    = std::max( npts_max, task.npts );
    nshells_max = std::max( nshells_max, task.bfn_screening.nshells );
  }

  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  eval_collocation_masked_combined_deriv1( ntasks, npts_max, nshells_max,
    static_stack.shells_device, aos_stack.device_tasks, 
    data->device_backend_->queue() );
#else
  auto aos_stack     = data->aos_stack;

  auto max_l = data->l_batched_shell_to_task.size() - 1;
  eval_collocation_shell_to_task_gradient( max_l, 
    data->l_batched_shell_to_task.data(), aos_stack.device_tasks,
    data->device_backend_->queue() );
#endif
  
}

void AoSScheme1Base::eval_collocation_hessian( XCDeviceData* _data ) {
#ifdef GAUXC_ENABLE_HIP
  GAUXC_GENERIC_EXCEPTION("Hessian NYI for HIP Backends");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto aos_stack     = data->aos_stack;

  auto max_l = data->l_batched_shell_to_task.size() - 1;
  eval_collocation_shell_to_task_hessian( max_l, 
    data->l_batched_shell_to_task.data(), aos_stack.device_tasks,
    data->device_backend_->queue() );
#endif
}





void AoSScheme1Base::inc_exc( XCDeviceData* _data ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;
  gdot( data->device_backend_->master_blas_handle(), data->total_npts_task_batch,
    base_stack.eps_eval_device, 1, base_stack.den_eval_device, 1, 
    static_stack.acc_scr_device, static_stack.exc_device );

}
void AoSScheme1Base::inc_nel( XCDeviceData* _data ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;
  gdot( data->device_backend_->master_blas_handle(), data->total_npts_task_batch,
    base_stack.weights_device, 1, base_stack.den_eval_device, 1, 
    static_stack.acc_scr_device, static_stack.nel_device );

}















void AoSScheme1Base::eval_uvvar_lda_rks( XCDeviceData* _data ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.bfn_screening.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  // Zero density
  auto base_stack    = data->base_stack;
  data->device_backend_->set_zero_async_master_queue( data->total_npts_task_batch, base_stack.den_eval_device, "Den Zero" );
    

  // Evaluate U variables
  auto aos_stack     = data->aos_stack;
  eval_uvvars_lda( ntasks, nbe_max, npts_max, 
    aos_stack.device_tasks, data->device_backend_->queue() );

}




void AoSScheme1Base::eval_uvvar_gga_rks( XCDeviceData* _data ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  size_t nbe_max = 0, npts_max = 0;
  for( auto& task : tasks ) {
    nbe_max  = std::max( nbe_max, task.bfn_screening.nbe );
    npts_max = std::max( npts_max, task.npts );
  }

  // Zero density + gradient
  auto base_stack    = data->base_stack;
  data->device_backend_->set_zero_async_master_queue( data->total_npts_task_batch, base_stack.den_eval_device, "Den Zero" );
  data->device_backend_->set_zero_async_master_queue( data->total_npts_task_batch, base_stack.den_x_eval_device, "DenX Zero" );
  data->device_backend_->set_zero_async_master_queue( data->total_npts_task_batch, base_stack.den_y_eval_device, "DenY Zero" );
  data->device_backend_->set_zero_async_master_queue( data->total_npts_task_batch, base_stack.den_z_eval_device, "DenZ Zero" );

  // Evaluate U variables
  auto aos_stack     = data->aos_stack;
  eval_uvvars_gga( ntasks, data->total_npts_task_batch, nbe_max, npts_max, 
    aos_stack.device_tasks, base_stack.den_x_eval_device, base_stack.den_y_eval_device,
    base_stack.den_z_eval_device, base_stack.gamma_eval_device, 
    data->device_backend_->queue() );

}

void AoSScheme1Base::eval_uvvar_lda_uks( XCDeviceData* ){

  GAUXC_GENERIC_EXCEPTION("UKS NOT YET IMPLEMENTED FOR DEVICE"); 

}

void AoSScheme1Base::eval_uvvar_gga_uks( XCDeviceData* ){

  GAUXC_GENERIC_EXCEPTION("UKS NOT YET IMPLEMENTED FOR DEVICE");

}


void AoSScheme1Base::eval_kern_exc_vxc_lda( const functional_type& func, 
  XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  if( !func.is_lda() ) GAUXC_GENERIC_EXCEPTION("XC Kernel not LDA!");

  auto base_stack    = data->base_stack;
  GauXC::eval_kern_exc_vxc_lda( func, data->total_npts_task_batch, 
    base_stack.den_eval_device, base_stack.eps_eval_device, 
    base_stack.vrho_eval_device, data->device_backend_->queue() );

  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.eps_eval_device, 1 );
  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.vrho_eval_device, 1 );

}


void AoSScheme1Base::eval_kern_exc_vxc_gga( const functional_type& func, 
  XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  if( !func.is_gga() ) GAUXC_GENERIC_EXCEPTION("XC Kernel not GGA!");

  auto base_stack    = data->base_stack;
  GauXC::eval_kern_exc_vxc_gga( func, data->total_npts_task_batch, 
    base_stack.den_eval_device, base_stack.gamma_eval_device, 
    base_stack.eps_eval_device, base_stack.vrho_eval_device, 
    base_stack.vgamma_eval_device, data->device_backend_->queue() );

  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.eps_eval_device, 1 );
  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.vrho_eval_device, 1 );
  hadamard_product( data->device_backend_->master_blas_handle(), data->total_npts_task_batch, 1, 
                    base_stack.weights_device, 1, base_stack.vgamma_eval_device, 1 );

}










void AoSScheme1Base::eval_xmat( double fac, XCDeviceData* _data, bool do_grad ){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  // Pack density matrix 
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  sym_pack_submat( ntasks, aos_stack.device_tasks, static_stack.dmat_device, 
    nbf, submat_block_size, data->device_backend_->queue() );


  // Sync blas streams with master stream
  data->device_backend_->sync_blas_pool_with_master();

  auto do_gemm = [&]( auto& handle, size_t npts, size_t nbe, auto* bf_ptr, auto* den_ptr, int ldden, auto* x_ptr ) {
    gemm( handle, DeviceBlasOp::NoTrans, DeviceBlasOp::NoTrans, npts, nbe, nbe, fac, bf_ptr, npts,
      den_ptr, ldden, 0., x_ptr, npts ); 
  };

  // Launch GEMM in round-robin
  const auto n_blas_streams = data->device_backend_->blas_pool_size();
  //size_t nsingle = 0;
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];
      auto den_ptr = task.bfn_screening.ncut > 1 ? task.nbe_scr : static_stack.dmat_device + task.bfn_screening.ibf_begin*(nbf+1);
      int  ldden   = task.bfn_screening.ncut > 1 ? task.bfn_screening.nbe : nbf;
      auto handle = data->device_backend_->blas_pool_handle( iT % n_blas_streams );
      do_gemm( handle, task.npts, task.bfn_screening.nbe, task.bf, den_ptr, ldden, task.zmat );
      if( do_grad ) {
        do_gemm( handle, task.npts, task.bfn_screening.nbe, task.dbfx, den_ptr, ldden, task.xmat_x );
        do_gemm( handle, task.npts, task.bfn_screening.nbe, task.dbfy, den_ptr, ldden, task.xmat_y );
        do_gemm( handle, task.npts, task.bfn_screening.nbe, task.dbfz, den_ptr, ldden, task.xmat_z );
      }
  }

  // Record completion of BLAS ops on master stream
  data->device_backend_->sync_master_with_blas_pool();

}















void AoSScheme1Base::inc_vxc( XCDeviceData* _data){

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  // Sync blas streams with master stream
  data->device_backend_->sync_blas_pool_with_master();

  // Launch SYR2K in round robin
  const auto n_blas_streams = data->device_backend_->blas_pool_size();
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];
    syr2k( data->device_backend_->blas_pool_handle(iT % n_blas_streams), 
      DeviceBlasUplo::Lower, DeviceBlasOp::Trans, task.bfn_screening.nbe, task.npts, 1.,
      task.bf, task.npts, task.zmat, task.npts, 0., task.nbe_scr,
      task.bfn_screening.nbe );
  }

  // Record completion of BLAS ops on master stream
  data->device_backend_->sync_master_with_blas_pool();

  // Increment global VXC
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  sym_task_inc_potential( ntasks, aos_stack.device_tasks, 
    static_stack.vxc_device, nbf, submat_block_size, 
    data->device_backend_->queue() );
}













void AoSScheme1Base::symmetrize_vxc( XCDeviceData* _data) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  const auto nbf = data->global_dims.nbf;
  auto static_stack  = data->static_stack;
  symmetrize_matrix( nbf, static_stack.vxc_device, nbf, 
    data->device_backend_->queue() ); 

}




void AoSScheme1Base::inc_exc_grad_lda( XCDeviceData* _data ) {
#ifdef GAUXC_ENABLE_HIP
  GAUXC_GENERIC_EXCEPTION("LDA Grad NYI for HIP Backends");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  const auto nshell = data->global_dims.nshells;
  increment_exc_grad_lda( nshell, 
    data->shell_to_task_stack.shell_to_task_device, 
    data->aos_stack.device_tasks,
    data->static_stack.exc_grad_device,
    data->device_backend_->queue() ); 
#endif
}

void AoSScheme1Base::inc_exc_grad_gga( XCDeviceData* _data ) {
#ifdef GAUXC_ENABLE_HIP
  GAUXC_GENERIC_EXCEPTION("GGA Grad NYI for HIP Backends");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  const auto nshell = data->global_dims.nshells;
  increment_exc_grad_gga( nshell, 
    data->shell_to_task_stack.shell_to_task_device, 
    data->aos_stack.device_tasks,
    data->static_stack.exc_grad_device,
    data->device_backend_->queue() ); 
#endif
}


void AoSScheme1Base::eval_exx_fmat( XCDeviceData* _data ) {
#ifndef GAUXC_ENABLE_EXX
  GAUXC_GENERIC_EXCEPTION("EXX F-Matrix NYI for non-CUDA Backends");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();
  const auto nbf = data->global_dims.nbf;
  auto static_stack  = data->static_stack;

  // Pack the density matrix into (bfn, cou) shape
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  auto aos_stack     = data->aos_stack;
  asym_pack_submat( ntasks, aos_stack.device_tasks, static_stack.dmat_device,
    nbf, submat_block_size, data->device_backend_->queue() );

  // Sync blas streams with master stream
  data->device_backend_->sync_blas_pool_with_master();

  // Launch GEMM in round-robin
  const auto n_blas_streams = data->device_backend_->blas_pool_size();
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];
    auto handle = data->device_backend_->blas_pool_handle( iT % n_blas_streams );
    auto npts = task.npts;
    auto nbe_bfn = task.bfn_screening.nbe;
    auto nbe_cou = task.cou_screening.nbe;
    gemm( handle, DeviceBlasOp::NoTrans, DeviceBlasOp::NoTrans, 
      npts, nbe_cou, nbe_bfn, 1., task.bf, npts, task.nbe_scr, nbe_bfn, 
      0., task.fmat, npts );
  }

  // Record completion of BLAS ops on master stream
  data->device_backend_->sync_master_with_blas_pool();
#endif
}

void AoSScheme1Base::eval_exx_gmat( XCDeviceData* _data, 
  const BasisSetMap& basis_map ) {
#ifndef GAUXC_ENABLE_EXX
  GAUXC_GENERIC_EXCEPTION("EXX G-Matrix NYI for non-CUDA Backends");
#else

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  //const auto ntasks = tasks.size();
  const size_t nshells = data->global_dims.nshells;
  //auto static_stack  = data->static_stack;

  // XXX: Need to add screening capabilities, packing etc
  //const auto nbf = data->global_dims.nbf;

  // XXX: Need to add support for non-cartesian functions
  for( auto i = 0ul; i < nshells; ++i ) {
    if( basis_map.shell_pure(i) )
      GAUXC_GENERIC_EXCEPTION("GPU EXX + Spherical NYI");
  }

  if( basis_map.max_l() > 2 ) {
    GAUXC_GENERIC_EXCEPTION("GPU EXX + L>2 NYI");
  }

  // Zero out G
  for( auto& task : tasks ) {
    const size_t sz = task.npts*task.cou_screening.nbe;
    data->device_backend_->set_zero_async_master_queue( 
      sz, task.gmat, "Zero G" );
  }

  // Sync blas streams with master stream
  data->device_backend_->sync_blas_pool_with_master();

  // Launch Shell Pair Kernels in round-robin
  //const auto n_streams = data->device_backend_->blas_pool_size();

  auto& sp_to_task = data->shell_pair_to_task;
  #if 1
  bool do_batch = true;

  if( do_batch ) { // start batched code

    cudaStream_t stream = 
      data->device_backend_->queue().queue_as<util::cuda_stream>();

    XGPU::integral_0_task_batched(
      tasks.size(), data->subtask.size(),
      data->l_batch_diag_task_to_shell_pair_device[0].max_prim_pairs, 0,
      data->aos_stack.device_tasks,
      data->l_batch_diag_task_to_shell_pair_device[0].task_to_shell_pair_device,
      data->task_to_shell_pair_stack.subtask_device,
      data->task_to_shell_pair_stack.nprim_pairs_device,
      data->task_to_shell_pair_stack.sp_ptr_device,
      data->task_to_shell_pair_stack.sp_X_AB_device,
      data->task_to_shell_pair_stack.sp_Y_AB_device,
      data->task_to_shell_pair_stack.sp_Z_AB_device,
      dev_boys_table, stream
    );
    data->device_backend_->check_error("integral_0_task_batched" __FILE__ ": " + std::to_string(__LINE__));
    if(basis_map.max_l() > 0) {
    XGPU::integral_1_task_batched(
      tasks.size(), data->subtask.size(),
      data->l_batch_diag_task_to_shell_pair_device[1].max_prim_pairs, 0,
      data->aos_stack.device_tasks,
      data->l_batch_diag_task_to_shell_pair_device[1].task_to_shell_pair_device,
      data->task_to_shell_pair_stack.subtask_device,
      data->task_to_shell_pair_stack.nprim_pairs_device,
      data->task_to_shell_pair_stack.sp_ptr_device,
      data->task_to_shell_pair_stack.sp_X_AB_device,
      data->task_to_shell_pair_stack.sp_Y_AB_device,
      data->task_to_shell_pair_stack.sp_Z_AB_device,
      dev_boys_table, stream
    );
    data->device_backend_->check_error("integral_1_task_batched" __FILE__ ": " + std::to_string(__LINE__));
    }
    if(basis_map.max_l() > 1) {
    XGPU::integral_2_task_batched(
      tasks.size(), data->subtask.size(),
      data->l_batch_diag_task_to_shell_pair_device[2].max_prim_pairs, 0,
      data->aos_stack.device_tasks,
      data->l_batch_diag_task_to_shell_pair_device[2].task_to_shell_pair_device,
      data->task_to_shell_pair_stack.subtask_device,
      data->task_to_shell_pair_stack.nprim_pairs_device,
      data->task_to_shell_pair_stack.sp_ptr_device,
      data->task_to_shell_pair_stack.sp_X_AB_device,
      data->task_to_shell_pair_stack.sp_Y_AB_device,
      data->task_to_shell_pair_stack.sp_Z_AB_device,
      dev_boys_table, stream
    );
    data->device_backend_->check_error("integral_2_task_batched" __FILE__ ": " + std::to_string(__LINE__));
    }

  #define SP_LBATCH_IDX(I,J) (I*(basis_map.max_l()+1) + J)

    XGPU::integral_0_0_task_batched(
      tasks.size(), data->subtask.size(),
      data->l_batch_task_to_shell_pair_device[0].max_prim_pairs, 0,
      data->aos_stack.device_tasks,
      data->l_batch_task_to_shell_pair_device[0].task_to_shell_pair_device,
      data->task_to_shell_pair_stack.subtask_device,
      data->task_to_shell_pair_stack.nprim_pairs_device,
      data->task_to_shell_pair_stack.sp_ptr_device,
      data->task_to_shell_pair_stack.sp_X_AB_device,
      data->task_to_shell_pair_stack.sp_Y_AB_device,
      data->task_to_shell_pair_stack.sp_Z_AB_device,
      dev_boys_table, stream
    );
    data->device_backend_->check_error("integral_0_0_task_batched" __FILE__ ": " + std::to_string(__LINE__));

    if(basis_map.max_l() > 0) {
    XGPU::integral_1_1_task_batched(
      tasks.size(), data->subtask.size(),
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(1,1)].max_prim_pairs, 0,
      data->aos_stack.device_tasks,
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(1,1)].task_to_shell_pair_device,
      data->task_to_shell_pair_stack.subtask_device,
      data->task_to_shell_pair_stack.nprim_pairs_device,
      data->task_to_shell_pair_stack.sp_ptr_device,
      data->task_to_shell_pair_stack.sp_X_AB_device,
      data->task_to_shell_pair_stack.sp_Y_AB_device,
      data->task_to_shell_pair_stack.sp_Z_AB_device,
      dev_boys_table, stream
    );
    data->device_backend_->check_error("integral_1_1_task_batched" __FILE__ ": " + std::to_string(__LINE__));
    }

    if(basis_map.max_l() > 1) {
    XGPU::integral_2_2_task_batched(
      tasks.size(), data->subtask.size(),
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(2,2)].max_prim_pairs, 0,
      data->aos_stack.device_tasks,
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(2,2)].task_to_shell_pair_device,
      data->task_to_shell_pair_stack.subtask_device,
      data->task_to_shell_pair_stack.nprim_pairs_device,
      data->task_to_shell_pair_stack.sp_ptr_device,
      data->task_to_shell_pair_stack.sp_X_AB_device,
      data->task_to_shell_pair_stack.sp_Y_AB_device,
      data->task_to_shell_pair_stack.sp_Z_AB_device,
      dev_boys_table, stream
    );
    data->device_backend_->check_error("integral_2_2_task_batched" __FILE__ ": " + std::to_string(__LINE__));
    }

    if(basis_map.max_l() > 0) {
    XGPU::integral_1_0_task_batched( true,
      tasks.size(), data->subtask.size(),
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(0,1)].max_prim_pairs, 0,
      data->aos_stack.device_tasks,
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(0,1)].task_to_shell_pair_device,
      data->task_to_shell_pair_stack.subtask_device,
      data->task_to_shell_pair_stack.nprim_pairs_device,
      data->task_to_shell_pair_stack.sp_ptr_device,
      data->task_to_shell_pair_stack.sp_X_AB_device,
      data->task_to_shell_pair_stack.sp_Y_AB_device,
      data->task_to_shell_pair_stack.sp_Z_AB_device,
      dev_boys_table, stream
    );
    data->device_backend_->check_error("integral_1_0_task_batched" __FILE__ ": " + std::to_string(__LINE__));
    }

    if(basis_map.max_l() > 0) {
    XGPU::integral_1_0_task_batched( false,
      tasks.size(), data->subtask.size(),
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(1,0)].max_prim_pairs, 0,
      data->aos_stack.device_tasks,
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(1,0)].task_to_shell_pair_device,
      data->task_to_shell_pair_stack.subtask_device,
      data->task_to_shell_pair_stack.nprim_pairs_device,
      data->task_to_shell_pair_stack.sp_ptr_device,
      data->task_to_shell_pair_stack.sp_X_AB_device,
      data->task_to_shell_pair_stack.sp_Y_AB_device,
      data->task_to_shell_pair_stack.sp_Z_AB_device,
      dev_boys_table, stream
    );
    data->device_backend_->check_error("integral_1_0_task_batched" __FILE__ ": " + std::to_string(__LINE__));
    }

    if(basis_map.max_l() > 1) {
    XGPU::integral_2_0_task_batched( true,
      tasks.size(), data->subtask.size(),
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(0,2)].max_prim_pairs, 0,
      data->aos_stack.device_tasks,
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(0,2)].task_to_shell_pair_device,
      data->task_to_shell_pair_stack.subtask_device,
      data->task_to_shell_pair_stack.nprim_pairs_device,
      data->task_to_shell_pair_stack.sp_ptr_device,
      data->task_to_shell_pair_stack.sp_X_AB_device,
      data->task_to_shell_pair_stack.sp_Y_AB_device,
      data->task_to_shell_pair_stack.sp_Z_AB_device,
      dev_boys_table, stream
    );
    data->device_backend_->check_error("integral_2_0_task_batched" __FILE__ ": " + std::to_string(__LINE__));
    }

    if(basis_map.max_l() > 1) {
    XGPU::integral_2_0_task_batched( false,
      tasks.size(), data->subtask.size(),
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(2,0)].max_prim_pairs, 0,
      data->aos_stack.device_tasks,
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(2,0)].task_to_shell_pair_device,
      data->task_to_shell_pair_stack.subtask_device,
      data->task_to_shell_pair_stack.nprim_pairs_device,
      data->task_to_shell_pair_stack.sp_ptr_device,
      data->task_to_shell_pair_stack.sp_X_AB_device,
      data->task_to_shell_pair_stack.sp_Y_AB_device,
      data->task_to_shell_pair_stack.sp_Z_AB_device,
      dev_boys_table, stream
    );
    data->device_backend_->check_error("integral_2_0_task_batched" __FILE__ ": " + std::to_string(__LINE__));
    }

    if(basis_map.max_l() > 1) {
    XGPU::integral_2_1_task_batched( true,
      tasks.size(), data->subtask.size(),
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(1,2)].max_prim_pairs, 0,
      data->aos_stack.device_tasks,
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(1,2)].task_to_shell_pair_device,
      data->task_to_shell_pair_stack.subtask_device,
      data->task_to_shell_pair_stack.nprim_pairs_device,
      data->task_to_shell_pair_stack.sp_ptr_device,
      data->task_to_shell_pair_stack.sp_X_AB_device,
      data->task_to_shell_pair_stack.sp_Y_AB_device,
      data->task_to_shell_pair_stack.sp_Z_AB_device,
      dev_boys_table, stream
    );
    data->device_backend_->check_error("integral_2_1_task_batched" __FILE__ ": " + std::to_string(__LINE__));
    }

    if(basis_map.max_l() > 1) {
    XGPU::integral_2_1_task_batched( false,
      tasks.size(), data->subtask.size(),
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(2,1)].max_prim_pairs, 0,
      data->aos_stack.device_tasks,
      data->l_batch_task_to_shell_pair_device[SP_LBATCH_IDX(2,1)].task_to_shell_pair_device,
      data->task_to_shell_pair_stack.subtask_device,
      data->task_to_shell_pair_stack.nprim_pairs_device,
      data->task_to_shell_pair_stack.sp_ptr_device,
      data->task_to_shell_pair_stack.sp_X_AB_device,
      data->task_to_shell_pair_stack.sp_Y_AB_device,
      data->task_to_shell_pair_stack.sp_Z_AB_device,
      dev_boys_table, stream
    );
    data->device_backend_->check_error("integral_2_1_task_batched" __FILE__ ": " + std::to_string(__LINE__));
    }

  } else { // end batched start unbatched

    cudaStream_t stream = 
      data->device_backend_->queue().queue_as<util::cuda_stream>();
    for( auto& sptt : sp_to_task ) { 
      size_t ntask_sp = sptt.task_idx.size();
      auto ish = sptt.idx_bra;
      auto jsh = sptt.idx_ket;
      for( auto i = 0ul; i < ntask_sp; i++ ) {
        const auto iT = sptt.task_idx[i];
        const auto i_off = sptt.task_shell_off_row[i];
        const auto j_off = sptt.task_shell_off_col[i];

        const auto& task = tasks[iT];
        //cudaStream_t stream = 
          //data->device_backend_->blas_pool_queue(iT % n_streams)
          //  .queue_as<util::cuda_stream>();

        XGPU::compute_integral_shell_pair( ish == jsh,
          task.npts,
          task.points_x,
          task.points_y,
          task.points_z,
          sptt.lA, sptt.lB,
          sptt.rA, sptt.rB,
          sptt.shell_pair_device,
          task.fmat + i_off*task.npts,
          task.fmat + j_off*task.npts,
          task.npts,
          task.gmat + i_off*task.npts,
          task.gmat + j_off*task.npts,
          task.npts,
          task.weights,
          dev_boys_table, stream ); 
      } // Loop over tasks within a shell pair
    } // Loop over shell pair maps
  } // end unbatched
  #else
  size_t isptt = 0;
  for( auto& sptt : sp_to_task ) {
    size_t ntask_sp = sptt.task_idx.size();
    auto ish = sptt.idx_bra;
    auto jsh = sptt.idx_ket;
    //std::cout << "SH " << ish << " " << jsh << std::endl;
    if( true ) {

      cudaStream_t stream = 
        data->device_backend_->queue().queue_as<util::cuda_stream>();
      const auto X_AB = sptt.rA.x - sptt.rB.x;
      const auto Y_AB = sptt.rA.y - sptt.rB.y;
      const auto Z_AB = sptt.rA.z - sptt.rB.z;
      XGPU::compute_integral_shell_pair_batched( ish == jsh, ntask_sp, 
        sptt.lA, sptt.lB, X_AB, Y_AB, Z_AB,
        data->shell_pair_to_task_stack.shell_pair_to_task_device + isptt,
        data->aos_stack.device_tasks, dev_boys_table, stream );

    } else {

      for( auto i = 0ul; i < ntask_sp; i++ ) {
        const auto iT = sptt.task_idx[i];
        const auto i_off = sptt.task_shell_off_row[i];
        const auto j_off = sptt.task_shell_off_col[i];

        const auto& task = tasks[iT];
        cudaStream_t stream = 
          data->device_backend_->queue().queue_as<util::cuda_stream>();
          //data->device_backend_->blas_pool_queue(iT % n_streams)
          //  .queue_as<util::cuda_stream>();

        XGPU::compute_integral_shell_pair( ish == jsh,
          task.npts,
          task.points_x,
          task.points_y,
          task.points_z,
          sptt.lA, sptt.lB,
          sptt.rA, sptt.rB,
          sptt.shell_pair_device,
          task.fmat + i_off*task.npts,
          task.fmat + j_off*task.npts,
          task.npts,
          task.gmat + i_off*task.npts,
          task.gmat + j_off*task.npts,
          task.npts,
          task.weights,
          dev_boys_table, stream ); 
      
      }

    }
    isptt++;
  }
  #endif


  // Record completion of BLAS ops on master stream
  data->device_backend_->sync_master_with_blas_pool();
#endif
}



void AoSScheme1Base::inc_exx_k( XCDeviceData* _data ) {
#ifndef GAUXC_ENABLE_EXX
  GAUXC_GENERIC_EXCEPTION("EXX + non-CUDA NYI");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto& tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  // Sync blas streams with master stream
  data->device_backend_->sync_blas_pool_with_master();

  // Launch GEMM in round-robin
  const auto n_blas_streams = data->device_backend_->blas_pool_size();
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];
    auto handle = data->device_backend_->blas_pool_handle( iT % n_blas_streams );
    auto npts = task.npts;
    auto nbe_bfn = task.bfn_screening.nbe;
    auto nbe_cou = task.cou_screening.nbe;
    gemm( handle, DeviceBlasOp::Trans, DeviceBlasOp::NoTrans, 
      nbe_bfn, nbe_cou, npts, 1., task.bf, npts, task.gmat, npts, 0., 
      task.nbe_scr, nbe_bfn );
  }

  // Record completion of BLAS ops on master stream
  data->device_backend_->sync_master_with_blas_pool();

  // Increment EXX_K
  const auto nbf = data->global_dims.nbf;
  const auto submat_block_size = data->get_submat_chunk_size( nbf, 0 );
  auto static_stack  = data->static_stack;
  auto aos_stack     = data->aos_stack;
  asym_task_inc_potential( ntasks, aos_stack.device_tasks, 
    static_stack.exx_k_device, nbf, submat_block_size, 
    data->device_backend_->queue() );
#endif
}

void AoSScheme1Base::symmetrize_exx_k( XCDeviceData* _data ) {
#ifndef GAUXC_ENABLE_EXX
  GAUXC_GENERIC_EXCEPTION("EXX + non-CUDA NYI");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  const auto nbf = data->global_dims.nbf;
  auto static_stack  = data->static_stack;
  symmetrize_matrix_inc( nbf, static_stack.exx_k_device, nbf, 
    data->device_backend_->queue() ); 
#endif
}


void AoSScheme1Base::eval_exx_ek_screening_bfn_stats( XCDeviceData* _data ) {
#ifndef GAUXC_ENABLE_EXX
  GAUXC_GENERIC_EXCEPTION("EXX + non-CUDA NYI");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  auto tasks = data->host_device_tasks;
  const auto ntasks_ek = data->global_dims.ntask_ek;
  const auto ntasks = tasks.size();
  //const auto nbf = data->global_dims.nbf;
  auto aos_stack    = data->aos_stack;
  auto static_stack    = data->static_stack;
  GauXC::exx_ek_screening_bfn_stats( ntasks, aos_stack.device_tasks,
    static_stack.ek_max_bfn_sum_device, static_stack.ek_bfn_max_device, 
    ntasks_ek, data->device_backend_->queue() );
#endif
}


void AoSScheme1Base::exx_ek_shellpair_collision( double eps_E, double eps_K,
  XCDeviceData* _data, host_task_iterator tb, host_task_iterator te,
  const ShellPairCollection<double>& shpairs ) {
#ifndef GAUXC_ENABLE_EXX
  GAUXC_GENERIC_EXCEPTION("EXX + non-CUDA NYI");
#else
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  if( not data->device_backend_ ) GAUXC_UNINITIALIZED_DEVICE_BACKEND();

  const auto ntasks = std::distance(tb, te);
  if( ntasks > data->global_dims.ntask_ek ) 
    GAUXC_GENERIC_EXCEPTION("EK - Too Many Tasks");

  const auto nshells   = data->global_dims.nshells;
  const auto nbf   = data->global_dims.nbf;
  auto static_stack    = data->static_stack;

  GauXC::exx_ek_shellpair_collision( ntasks, nshells, nbf,
    static_stack.dmat_device, nbf,
    static_stack.vshell_max_sparse_device, 
    static_stack.shpair_row_ind_device,
    static_stack.shpair_col_ind_device,
    static_stack.ek_max_bfn_sum_device,
    static_stack.ek_bfn_max_device, data->global_dims.ntask_ek, 
    static_stack.shells_device, static_stack.shell_to_bf_device,
    static_stack.shell_sizes_device, eps_E, eps_K,
    data->dynmem_ptr, data->dynmem_sz,
    tb, te, shpairs,
    data->device_backend_->queue(),
    data->device_backend_->master_blas_handle()
   );
#endif
}


}
