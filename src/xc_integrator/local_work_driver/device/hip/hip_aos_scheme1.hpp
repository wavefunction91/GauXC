#pragma once
#include "device/scheme1_base.hpp"

namespace GauXC {

struct HipAoSScheme1 : public AoSScheme1Base {

  // Algorithmic constants
  //static constexpr uint32_t weight_unroll = 4;
  //static constexpr uint32_t weight_thread_block = 640;
  //static constexpr uint32_t weight_thread_block_per_sm = 2;
  static constexpr uint32_t max_submat_blocks = 10;

  // API Overrides
  void partition_weights( XCDeviceData* ) override final;

  std::unique_ptr<XCDeviceData> create_device_data() override final;

  struct Data;

};


struct HipAoSScheme1::Data : public Scheme1DataBase {

  virtual ~Data() noexcept;
  Data(bool batch_l3_blas = true);

  // Final overrides
  size_t get_submat_chunk_size(int32_t,int32_t) override final;
  size_t get_ldatoms() override final;
  size_t get_rab_align() override final;

};


}
