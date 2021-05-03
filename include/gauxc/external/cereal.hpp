#pragma once

#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/access.hpp>

namespace GauXC {
  template <typename F> class BasisSet; // fwd declare basis set
}

namespace cereal {
  template <typename Archive, typename F>
  struct specialize<Archive, GauXC::BasisSet<F>, specialization::member_serialize> {};
}
