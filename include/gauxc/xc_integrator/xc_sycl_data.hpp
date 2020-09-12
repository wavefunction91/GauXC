#pragma once
#include <vector>
#include <cstdint>
#include <memory>
#include <gauxc/basisset.hpp>
#include <gauxc/xc_task.hpp>

#include <gauxc/gauxc_config.hpp>
#include <gauxc/util/sycl_util.hpp>


#ifdef GAUXC_ENABLE_SYCL

namespace GauXC {

template <typename F>
class XCSyclData {
public:
  size_t nshells  = 0;
  size_t nbf      = 0;
  size_t n_deriv  = 0;
  size_t natoms   = 0;

  bool denpack_host = false;
  bool vxcinc_host  = false;

  void* device_ptr = nullptr;
  void* dynmem_ptr = nullptr;
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
    F*        alpha_array_device = nullptr;
    F*        beta_array_device  = nullptr;
  oneapi::mkl::tranpose* transA_array_device = nullptr;
  oneapi::mkl::tranpose* transB_array_device = nullptr;

  F*     dist_scratch_device = nullptr;

  // Buffer Vars
  F*           points_device_buffer     = nullptr;
  F*           weights_device_buffer    = nullptr;
  size_t*      shell_list_device_buffer = nullptr;
  size_t*      shell_offs_device_buffer = nullptr;
  int64_t*     submat_cut_device_buffer = nullptr;
  int32_t*     iparent_device_buffer    = nullptr;
  F*           dist_nearest_buffer      = nullptr;

  sycl::XCTaskDevice<F>* device_tasks  = nullptr;

  // // Execution management
  std::unique_ptr<cl::sycl::queue> master_queue = nullptr;

  XCSyclData( size_t _natoms,
              size_t _n_deriv,
              size_t _nbf,
              size_t _nshells,
              bool _denpack_host = false,
              bool _vxcinc_host  = false );

  ~XCSyclData() noexcept;
  XCSyclData( const XCSyclData& )          = delete;
  XCSyclData( XCSyclData&&      ) noexcept = delete;


  using task_iterator = std::vector< XCTask >::iterator;
  using device_task_container = std::vector< sycl::XCTaskDevice<F> >;

  std::tuple< task_iterator, device_task_container >
    generate_buffers( const BasisSet<F>& basis,
                      task_iterator      task_begin,
                      task_iterator      task_end    );

};

}

#endif
