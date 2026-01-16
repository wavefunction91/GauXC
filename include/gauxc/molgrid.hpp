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

#include <gauxc/molecule.hpp>
#include <gauxc/grid.hpp>
#include <gauxc/grid_factory.hpp>

#include <unordered_map>

namespace GauXC {

/// Map from atomic number to concrete atom-centered Grid instances.
using atomic_grid_map = std::unordered_map< AtomicNumber, Grid >;
/// Map from atomic number to atom-centered grid specification variants.
using atomic_grid_spec_map = std::unordered_map< AtomicNumber, atomic_grid_variant>;

namespace detail {
  /// Forward declaration of the MolGrid implementation.
  class MolGridImpl;
}

/**
 *  @brief Molecular integration grid container.
 *
 *  MolGrid aggregates atom-centered quadrature grids for each atomic number
 *  and exposes access to per-element grids and batch sizing information for
 *  numerical integration.
 */
class MolGrid {

  /// Handle for internal implementation details.
  std::shared_ptr<detail::MolGridImpl> pimpl_;

public:

  /**
   *  @brief Construct from pre-built atom-centered grids.
   *  @param grids Map from atomic number to Grid instances.
   */
  MolGrid( const atomic_grid_map& );

  /**
   *  @brief Construct from grid specification variants.
   *  @param grid_specs Map from atomic number to grid specification variants.
   */
  MolGrid( const atomic_grid_spec_map& );

  /// Copy constructor.
  MolGrid( const MolGrid& );

  /// Move constructor.
  MolGrid( MolGrid&& ) noexcept;

  /// Destructor.
  ~MolGrid() noexcept;

  /**
   *  @brief Get the count of unique atomic numbers represented.
   *  @return Number of distinct atomic numbers in the grid map.
   */
  size_t natoms_uniq() const;

  /**
   *  @brief Retrieve the atom-centered grid for a given atomic number.
   *  @param Z Atomic number to look up.
   *  @return Const reference to the associated Grid.
   */
  const Grid& get_grid( AtomicNumber ) const;

  /**
   *  @brief Retrieve the atom-centered grid for a given atomic number (mutable).
   *  @param Z Atomic number to look up.
   *  @return Reference to the associated Grid.
   */
        Grid& get_grid( AtomicNumber );

  /**
   *  @brief Get the maximum batch count across all stored grids.
   *  @return Largest number of batches among all atom grids.
   */
  size_t max_nbatches() const;

};

}
