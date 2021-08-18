#pragma once

#include "device/xc_device_aos_data.hpp"
#include <gauxc/util/cuda_util.hpp>
#include <gauxc/util/cublas_util.hpp>

namespace GauXC {


struct XCCudaAoSDataBase : public XCDeviceAoSData {

  // Final API overrides
  void allocate_device_buffer(int64_t sz = -1) override final;
  void master_queue_synchronize() override final;
  void copy_async_(size_t,const void*,void*,std::string) override final;
  void set_zero_( size_t, void*, std::string ) override final;


  XCCudaAoSDataBase();
  virtual ~XCCudaAoSDataBase() noexcept;

  // Execution management
  std::unique_ptr<util::cuda_stream>   master_stream = nullptr;
  std::unique_ptr<util::cublas_handle> master_handle = nullptr;

};

}
