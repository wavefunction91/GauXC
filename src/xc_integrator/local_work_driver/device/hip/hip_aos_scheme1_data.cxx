#include "hip_aos_scheme1.hpp"
#include "buffer_adaptor.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device/hip/hip_backend.hpp"

namespace GauXC {

HipAoSScheme1::Data::~Data() noexcept = default;

HipAoSScheme1::Data::Data() :
  Scheme1DataBase( std::make_unique<HIPBackend>() ) { }

size_t HipAoSScheme1::Data::get_ldatoms() {
  //constexpr auto weight_unroll = HipAoSScheme1::weight_unroll;
  constexpr auto weight_unroll = 1;
  return util::div_ceil( global_dims.natoms, weight_unroll ) * weight_unroll;
}

size_t HipAoSScheme1::Data::get_rab_align() {
  return sizeof(double2);
}


size_t HipAoSScheme1::Data::get_submat_chunk_size(int32_t LDA, int32_t dev_id) {

  //constexpr auto max_submat_blocks = HipAoSScheme1::max_submat_blocks;

  int l2_cache_size;
  hipDeviceGetAttribute(&l2_cache_size, hipDeviceAttributeL2CacheSize, dev_id);

  int l2_block_size = (int) sqrt(0.75 * ((double) l2_cache_size / 8));
  int min_block_size = LDA / max_submat_blocks;

  int block_size = std::max(l2_block_size, min_block_size);
  block_size = std::min(block_size, LDA);

  return block_size;

}


}
