/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include <gauxc/atom.hpp>
#include <vector>
#include <algorithm>

#include <gauxc/gauxc_config.hpp>

namespace GauXC {

class Molecule : public std::vector<Atom> {

public:

  template <typename... Args>
  Molecule( Args&&... args ) :
    std::vector<Atom>( std::forward<Args>(args)... ) { }

  Molecule( const Molecule& )          = default;
  Molecule( Molecule&&      ) noexcept = default;

  Molecule& operator=( const Molecule& other ) = default;

  size_t natoms() const { return this->size(); }

  AtomicNumber maxZ() const {
    return std::max_element( this->cbegin(), this->cend(),
      []( const auto& a, const auto& b) {
        return a.Z.get() < b.Z.get();
      })->Z;
  }

  bool operator==(const Molecule& other) {
    if(other.size() != this->size()) return false;
    for( auto i = 0ul; i < this->size(); ++i )
      if( other[i] != operator[](i) ) return false;
    return true;
  }
};

}
