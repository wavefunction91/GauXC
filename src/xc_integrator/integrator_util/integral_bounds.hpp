#pragma once

#include <gauxc/shell.hpp>

namespace GauXC {
namespace util  {

template <typename T>
T max_coulomb( const Shell<T>& bra, const Shell<T>& ket);

extern template double max_coulomb( const Shell<double>&, const Shell<double>& );

}
}
