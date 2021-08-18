#pragma once
#include "device/local_device_work_driver_pimpl.hpp"
#include "device/cuda/xc_cuda_aos_data_base.hpp"

namespace GauXC {
namespace detail {

struct CudaAoSScheme1 : public LocalDeviceWorkDriverPIMPL {

  // Algorithmic constants
  static constexpr uint32_t weight_unroll = 4;
  static constexpr uint32_t weight_thread_block = 640;
  static constexpr uint32_t weight_thread_block_per_sm = 2;
  static constexpr uint32_t max_submat_blocks = 10;

  // API Overrides
  void partition_weights( XCDeviceData* ) override final;
  void eval_collocation( XCDeviceData* ) override final;
  void eval_collocation_gradient( XCDeviceData* ) override final;
  void eval_xmat( XCDeviceData* ) override final;
  void eval_uvvar_lda( XCDeviceData* ) override final;
  void eval_uvvar_gga( XCDeviceData* ) override final;
  void eval_zmat_lda_vxc( XCDeviceData* ) override final;
  void eval_zmat_gga_vxc( XCDeviceData* ) override final;
  void inc_vxc( XCDeviceData* ) override final;

  std::unique_ptr<XCDeviceData> create_device_data() override final;

  struct Data;
};


struct CudaAoSScheme1::Data : public XCCudaAoSDataBase {

  using base_type = XCCudaAoSDataBase;
  using base_type::host_task_type;
  using base_type::device_buffer_t;

  double*  dist_scratch_device = nullptr;
  double*  dist_nearest_device = nullptr;
  int32_t* iparent_device      = nullptr;

  std::vector<util::cuda_stream>   blas_streams;
  std::vector<util::cublas_handle> blas_handles;

  virtual ~Data() noexcept;
  Data(bool batch_l3_blas = true);

  size_t get_ldatoms();

  // Final overrides
  device_buffer_t add_extra_to_indirection(std::vector<XCDeviceTask>&, 
    device_buffer_t ) override final;

  void allocate_rab() override final;
  void send_rab(const MolMeta&) override final;
  size_t get_submat_chunk_size(int32_t,int32_t) override final;

  // Overrideable API's
  virtual size_t get_mem_req( const host_task_type&, const BasisSetMap&) override;
  virtual size_t get_static_mem_requirement() override; 


};


}
}
