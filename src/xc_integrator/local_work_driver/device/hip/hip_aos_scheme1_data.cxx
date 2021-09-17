#include "hip_aos_scheme1.hpp"
#include "buffer_adaptor.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device/hip/hip_backend.hpp"

namespace GauXC {

template <typename Base>
HipAoSScheme1<Base>::Data::~Data() noexcept = default;

template <typename Base>
HipAoSScheme1<Base>::Data::Data() :
  Base::Data( std::make_unique<HIPBackend>() ) { }

template <typename Base>
size_t HipAoSScheme1<Base>::Data::get_ldatoms() {
  //constexpr auto weight_unroll = HipAoSScheme1<Base>::weight_unroll;
  constexpr auto weight_unroll = 1;
  return util::div_ceil( this->global_dims.natoms, weight_unroll ) * weight_unroll;
}

template <typename Base>
size_t HipAoSScheme1<Base>::Data::get_rab_align() {
  return sizeof(double2);
}


template <typename Base>
size_t HipAoSScheme1<Base>::Data::get_submat_chunk_size(int32_t LDA, int32_t dev_id) {

  constexpr int max_submat_blocks = 10;
  int l2_cache_size;
  hipDeviceGetAttribute(&l2_cache_size, hipDeviceAttributeL2CacheSize, dev_id);

  int l2_block_size = (int) sqrt(0.75 * ((double) l2_cache_size / 8));
  int min_block_size = LDA / max_submat_blocks;

  int block_size = std::max(l2_block_size, min_block_size);
  block_size = std::min(block_size, LDA);

  return block_size;

}

template struct HipAoSScheme1<AoSScheme1Base>::Data;
#ifdef GAUXC_ENABLE_MAGMA
template struct HipAoSScheme1<AoSScheme1MAGMABase>::Data;
#endif

}
