#pragma once
#include "device/device_backend.hpp"
#include <memory>
#include <gauxc/util/cuda_util.hpp>
#include <gauxc/util/cublas_util.hpp>

namespace GauXC {

struct CUDABackend : public DeviceBackend {

  device_buffer_t   allocate_device_buffer(int64_t sz) override final;
  size_t            get_available_mem() override final;
  void              free_device_buffer( void* ptr ) override final;
  void              master_queue_synchronize() override final;
  void              create_blas_queue_pool(int32_t)   override final;
  type_erased_queue queue() override final;

  void copy_async_( size_t sz, const void* src, void* dest, 
                    std::string msg ) override final;
  void set_zero_( size_t sz, void* data, std::string msg) override final;

  void copy_async_2d_( size_t M, size_t N, const void* A, size_t LDA,
    void* B, size_t LDB, std::string msg ) override final;

  CUDABackend();
  ~CUDABackend() noexcept;

  // Execution management
  std::shared_ptr<util::cuda_stream>   master_stream = nullptr;
  std::shared_ptr<util::cublas_handle> master_handle = nullptr;

  std::vector<util::cuda_stream>   blas_streams;
  std::vector<util::cublas_handle> blas_handles;
};

}