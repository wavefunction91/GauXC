#pragma once
#include <type_traits>

namespace GauXC {

namespace detail {
template <typename T, typename... Args>
using all_are_not = std::conjunction< std::negation<std::is_same<T,Args>>... >;

template <typename T, typename... Args>
using enable_if_all_are_not = std::enable_if< all_are_not<T,Args...>::value >;

template <typename T, typename... Args>
using enable_if_all_are_not_t = typename enable_if_all_are_not<T,Args...>::type;
}

}
