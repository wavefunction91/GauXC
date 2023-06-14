#pragma once
#include <vector>
#include <cstdint>
#include <memory>
#include <gauxc/basisset.hpp>
#include <gauxc/xc_task.hpp>

#include <gauxc/gauxc_config.hpp>
#include <gauxc/util/cuda_util.hpp>
#include <gauxc/util/cublas_util.hpp>
#include <gauxc/util/magma_util.hpp>

#ifdef GAUXC_ENABLE_CUDA

namespace GauXC {

template <typename F>
struct XCCudaData {

  size_t nshells  = 0;
  size_t nbf      = 0;
  size_t n_deriv  = 0;
  size_t natoms   = 0;
  size_t LDatoms  = 0;

  bool batch_l3_blas = true;
  
  void* device_ptr = nullptr;
  void* dynmem_ptr = nullptr;
  size_t devmem_sz = 0;
  size_t dynmem_sz = 0;
   
  Shell<F>* shells_device             = nullptr;
  Shell<F>* important_shells_device   = nullptr;

  F*      vxc_device        = nullptr;
  F*      nbe_scr_device    = nullptr;
  F*      dmat_device       = nullptr;
  F*      zmat_device       = nullptr;
  F*      bf_eval_device    = nullptr;

  F*      dbf_x_eval_device = nullptr;
  F*      dbf_y_eval_device = nullptr;
  F*      dbf_z_eval_device = nullptr;

  F*      den_eval_device   = nullptr;
  F*      den_x_eval_device = nullptr;
  F*      den_y_eval_device = nullptr;
  F*      den_z_eval_device = nullptr;
  F*      eps_eval_device   = nullptr;
  F*      gamma_eval_device = nullptr;

  F*      vrho_eval_device    = nullptr;
  F*      vgamma_eval_device  = nullptr;


  F*     exc_device = nullptr;
  F*     nel_device = nullptr;
  F*     acc_scr_device = nullptr;

  F*     rab_device    = nullptr;
  F*     coords_device = nullptr;

  F**    dmat_array_device = nullptr;
  F**    zmat_array_device = nullptr;
  F**    bf_array_device   = nullptr;

  int*        m_array_device   = nullptr;
  int*        n_array_device   = nullptr;
  int*        k_array_device   = nullptr;
  int*        lda_array_device = nullptr;
  int*        ldb_array_device = nullptr;
  int*        ldc_array_device = nullptr;

  F*     dist_scratch_device = nullptr;

  // Buffer Vars
  F*           points_device_buffer     = nullptr;
  F*           weights_device_buffer    = nullptr;
  size_t*      shell_list_device_buffer = nullptr;
  size_t*      shell_offs_device_buffer = nullptr;
  int32_t*     submat_cut_device_buffer = nullptr;
  int32_t*     submat_block_device_buffer = nullptr;
  int32_t*     iparent_device_buffer    = nullptr;
  F*           dist_nearest_buffer      = nullptr;

  cuda::XCTaskDevice<F>* device_tasks  = nullptr;

  // Execution management
  std::unique_ptr<util::cuda_stream>   master_stream      = nullptr;
  std::unique_ptr<util::cublas_handle> master_handle      = nullptr;

#ifdef GAUXC_ENABLE_MAGMA
  std::unique_ptr<util::magma_queue>   master_magma_queue = nullptr;
#endif

  std::vector<util::cuda_stream>       blas_streams;
  std::vector<util::cublas_handle>     blas_handles;

  XCCudaData( bool _batch_l3_blas = true );

  ~XCCudaData() noexcept;
  XCCudaData( const XCCudaData& )          = delete;
  XCCudaData( XCCudaData&&      ) noexcept = delete;


  using task_iterator = std::vector< XCTask >::iterator;
  using device_task_container = std::vector< cuda::XCTaskDevice<F> >;


  void allocate_static_data( size_t _natoms,
                             size_t _n_deriv, 
                             size_t _nbf,
                             size_t _nshells );


  std::tuple< task_iterator, device_task_container >
    generate_buffers( const BasisSet<F>& basis,
                      task_iterator      task_begin,
                      task_iterator      task_end    );
 
};

}

#endif
