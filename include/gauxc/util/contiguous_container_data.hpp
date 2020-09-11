#pragma once
#include <gauxc/gauxc_config.hpp>
#include <array>

namespace GauXC  {
namespace detail {

template <typename T, size_t N>
inline HOST_DEVICE_ACCESSIBLE T* contiguous_data( const std::array<T,N>& arr ) {
  return reinterpret_cast<T*>( &const_cast<std::array<T,N>&>(arr) );
}

}
}
