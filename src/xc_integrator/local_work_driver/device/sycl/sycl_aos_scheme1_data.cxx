/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "sycl_aos_scheme1.hpp"
#include "buffer_adaptor.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device/sycl/sycl_backend.hpp"
#include "device_specific/sycl_vector_types.hpp"

namespace GauXC {

template <typename Base>
SyclAoSScheme1<Base>::Data::~Data() noexcept = default;

template <typename Base>
SyclAoSScheme1<Base>::Data::Data(const DeviceRuntimeEnvironment& rt) :
  Base::Data( rt ) { }

template <typename Base>
size_t SyclAoSScheme1<Base>::Data::get_ldatoms() {
  constexpr auto weight_unroll = alg_constants::SyclAoSScheme1::weight_unroll;
  return util::div_ceil( this->global_dims.natoms, weight_unroll ) * weight_unroll;
}

template <typename Base>
size_t SyclAoSScheme1<Base>::Data::get_rab_align() {
  return sizeof(double2);
}

template <typename Base>
int SyclAoSScheme1<Base>::Data::get_points_per_subtask() {
  return alg_constants::SyclAoSScheme1::ObaraSaika::points_per_subtask;
}



template <typename Base>
size_t SyclAoSScheme1<Base>::Data::get_submat_chunk_size(int32_t LDA, 
  int32_t dev_id) {

  constexpr auto max_submat_blocks = 
    alg_constants::SyclAoSScheme1::max_submat_blocks;

  (void)(dev_id); // The chunk is sized against the device this Data is bound to

  auto* backend = dynamic_cast<SYCLBackend*>(this->device_backend_);
  if( !backend ) GAUXC_BAD_BACKEND_CAST();

  // SYCL exposes the last-level cache size rather than an explicit L2 query;
  // on the supported targets these are the same cache
  const size_t l2_cache_size =
    backend->device.get_info< ::sycl::info::device::global_mem_cache_size >();

  int l2_block_size = (int) sqrt(0.75 * ((double) l2_cache_size / 8));
  int min_block_size = LDA / max_submat_blocks;

  int block_size = std::max(l2_block_size, min_block_size);
  block_size = std::min(block_size, LDA);

  return block_size;

}

template struct SyclAoSScheme1<AoSScheme1Base>::Data;
template struct SyclAoSScheme1<AoSScheme1OneMKLBase>::Data;


}
