/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#pragma once
#include <type_traits>
#include <cstdlib>
#include <cstdint>
#include <cassert>

namespace GauXC  {
namespace util  {

namespace detail {

template <typename... Args>
struct are_integral;

template <typename Head, typename... Tail>
struct are_integral<Head, Tail...> {
  static constexpr bool value = std::is_integral<Head>::value and 
                         are_integral<Tail...>::value;
};

template <typename T>
struct are_integral<T> {
  static constexpr bool value = std::is_integral<T>::value;
};


template <typename... Args>
struct largest;

template <typename Head, typename... Tail>
struct largest< Head, Tail... > {
private:
  using tail_type = typename largest<Tail...>::type;
public:
  using type = std::conditional_t< 
                 (sizeof(Head) > sizeof(tail_type)),
                 Head, tail_type >;
};

template <typename T>
struct largest<T> {
  using type = T;
};

template <typename... Args>
using largest_t = typename largest<Args...>::type;

}

template <typename Integral1, typename Integral2>
intmax_t div_ceil( Integral1 i, Integral2 j ) {

  static_assert( detail::are_integral<Integral1, Integral2>::value );
  assert( i >= 0 );
  assert( j >  0 );

  intmax_t i_us = i;
  intmax_t j_us = j;

  auto d = std::div(i_us,j_us);
  return d.quot + !!d.rem;

}


}
}
