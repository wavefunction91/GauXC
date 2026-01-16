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
#pragma once

#include <gauxc/atom.hpp>
#include <vector>
#include <algorithm>

#include <gauxc/gauxc_config.hpp>

namespace GauXC {

/**
 *  @brief A class to represent a collection of atoms in a molecule.
 */
class Molecule : public std::vector<Atom> {
private:
  /// Tests if the base class can be constructed from @p Args.
  template <typename... Args>
  static constexpr auto can_construct_base_v =
    std::is_constructible_v<std::vector<Atom>, Args...>;

public:
  /**
   *  @brief Construct a molecule from any arguments accepted by
   *         `std::vector<Atom>`.
   *  @tparam Args Arguments forwarded to the base vector constructor.
   */
  template <typename... Args,
            typename = std::enable_if_t<can_construct_base_v<Args...>>>
  Molecule( Args&&... args ) :
    std::vector<Atom>( std::forward<Args>(args)... ) { }

  /**
   *  @brief Default copy constructor.
   */
  Molecule( const Molecule& )          = default;

  /**
   *  @brief Default move constructor.
   */
  Molecule( Molecule&&      ) noexcept = default;

  /**
   *  @brief Default copy assignment.
   */
  Molecule& operator=( const Molecule& other ) = default;

  /**
   *  @brief Get the number of atoms in the molecule.
   *  @return The number of atoms stored in the molecule.
   */
  size_t natoms() const { return this->size(); }

  /**
   *  @brief Get the maximum atomic number / nuclear charge in the molecule.
   *  @return The largest atomic number among all atoms in the molecule.
   *  @warning Undefined behavior if the molecule is empty.
   */
  AtomicNumber maxZ() const {
    return std::max_element( this->cbegin(), this->cend(),
      []( const auto& a, const auto& b) {
        return a.Z.get() < b.Z.get();
      })->Z;
  }

  /**
   *  @brief Compare two molecules for equality.
   *
   *  Checks if all atoms are equal in properties and order.
   *
   *  @param other Molecule to compare against.
   *  @return True if all atoms are equal and in the same order.
   */
  bool operator==(const Molecule& other) {
    if(other.size() != this->size()) return false;
    for( auto i = 0ul; i < this->size(); ++i )
      if( other[i] != operator[](i) ) return false;
    return true;
  }
};

}
