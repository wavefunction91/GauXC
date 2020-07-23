#pragma once

#include "atom.hpp"
#include <vector>
#include <algorithm>

namespace GauXC {

class Molecule : public std::vector<Atom> {

public:

  template <typename... Args>
  Molecule( Args&&... args ) :
    std::vector<Atom>( std::forward<Args>(args)... ) { }

  Molecule( const Molecule& )          = default;
  Molecule( Molecule&&      ) noexcept = default;

  size_t natoms() const { return this->size(); }

  AtomicNumber maxZ() const {
    return std::max_element( this->cbegin(), this->cend(),
      []( const auto& a, const auto& b) {
        return a.Z.get() < b.Z.get();
      })->Z;
  }

};

}
