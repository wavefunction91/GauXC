#pragma once

#include <utility>

namespace GauXC::util {

inline static void unused() { }

template <typename T, typename... Args>
inline static void unused( const T& t, Args&&... args ) {
  (void)(t);
  unused( std::forward<Args>(args)... );
}


}
