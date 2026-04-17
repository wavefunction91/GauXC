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

#include <gauxc/molgrid.hpp>

namespace GauXC {

  /// Return Slater radius (Bohr) for the given atomic number (64-element table).
  double slater_radius_64(AtomicNumber);

  /// Return Slater radius (Bohr) for the given atomic number (30-element table).
  double slater_radius_30(AtomicNumber);

  /// Return Clementi radius (Bohr) for the given atomic number (67-element table).
  double clementi_radius_67(AtomicNumber);

  /// Return UFF radius (Bohr) for the given atomic number (103-element table).
  double uff_radius_103(AtomicNumber);

  /// Return the default atomic radius (Bohr) for a given atomic number.
  double default_atomic_radius(AtomicNumber);

  /// Return the default Mura-Knowles radial scaling factor for an atomic number.
  RadialScale default_mk_radial_scaling_factor( AtomicNumber );

  /// Return the default Murray-Handy-Laming radial scaling factor for an atomic number.
  RadialScale default_mhl_radial_scaling_factor( AtomicNumber );

  /// Return the default Treutler-Ahlrichs radial scaling factor for an atomic number.
  RadialScale default_ta_radial_scaling_factor( AtomicNumber );

  /**
   *  @brief Return the default radial scaling factor for a radial quadrature and atomic number.
   *  @param rq Radial quadrature type.
   *  @param Z  Atomic number.
   *  @return Radial scaling factor.
   */
  RadialScale default_radial_scaling_factor( RadialQuad rq, AtomicNumber Z );

  /**
   *  @brief Return the default radial and angular grid sizes for a given element and quadrature.
   *  @param Z       Atomic number.
   *  @param rq      Radial quadrature type.
   *  @param gs_def  Atomic grid size preset.
   *  @return Tuple of (RadialSize, AngularSize).
   */
  std::tuple<RadialSize,AngularSize> 
    default_grid_size(AtomicNumber Z, RadialQuad rq, AtomicGridSizeDefault gs_def); 

  /**
   *  @brief Factory for constructing molecular grids with default settings.
   *
   *  Provides static methods to generate atomic grid specifications, grid maps,
   *  and complete MolGrid instances using sensible defaults for radial/angular
   *  quadratures and pruning schemes.
   */
  struct MolGridFactory {

    /**
     *  @brief Create an unpruned atomic grid specification with default scaling.
     *  @param Z      Atomic number.
     *  @param rq     Radial quadrature type.
     *  @param rsize  Number of radial points.
     *  @param asize  Number of angular points.
     *  @return Unpruned atomic grid specification.
     */
    static UnprunedAtomicGridSpecification create_default_unpruned_grid_spec(
      AtomicNumber Z, RadialQuad rq, RadialSize rsize, AngularSize asize
    );

    /**
     *  @brief Create an unpruned atomic grid specification using a preset size.
     *  @param Z       Atomic number.
     *  @param rq      Radial quadrature type.
     *  @param gs_def  Atomic grid size preset.
     *  @return Unpruned atomic grid specification.
     */
    static UnprunedAtomicGridSpecification create_default_unpruned_grid_spec(
      AtomicNumber Z, RadialQuad rq, AtomicGridSizeDefault gs_def
    );

    /**
     *  @brief Create a pruned atomic grid specification.
     *  @tparam Args Argument types forwarded to create_default_unpruned_grid_spec.
     *  @param scheme Pruning scheme to apply.
     *  @param args   Arguments for unpruned spec creation.
     *  @return Pruned or unpruned atomic grid specification variant.
     */
    template <typename... Args>
    inline static atomic_grid_variant 
      create_default_pruned_grid_spec( PruningScheme scheme, Args&&... args ) {
      return create_pruned_spec( scheme, 
        create_default_unpruned_grid_spec(std::forward<Args>(args)...)
      );
    }

    /**
     *  @brief Create a map of atomic grid specifications for all unique elements in a molecule.
     *  @tparam Args Argument types forwarded to create_default_pruned_grid_spec.
     *  @param mol    Molecule containing the atoms.
     *  @param scheme Pruning scheme to apply.
     *  @param args   Additional arguments for grid specification.
     *  @return Map from atomic number to grid specification variant.
     */
    template <typename... Args>
    inline static atomic_grid_spec_map create_default_grid_spec_map( 
      const Molecule& mol, PruningScheme scheme, Args&&... args ) {

      atomic_grid_spec_map molmap;
      for( const auto& atom : mol ) 
      if( !molmap.count(atom.Z) ) {
        molmap.emplace( atom.Z, 
          create_default_pruned_grid_spec(scheme, atom.Z, 
            std::forward<Args>(args)...)
        );
      }

      return molmap;
    }

    /**
     *  @brief Generate concrete grids from a grid specification map.
     *  @param gs_map Map from atomic number to grid specification.
     *  @param bsz    Batch size for the generated grids.
     *  @return Map from atomic number to Grid instances.
     */
    inline static atomic_grid_map generate_gridmap(
      const atomic_grid_spec_map& gs_map, BatchSize bsz ) {

      atomic_grid_map molmap;
      for( const auto& [key, val] : gs_map ) {
        molmap.emplace( key, AtomicGridFactory::generate_grid(val, bsz) );
      }
      return molmap;

    }

    /**
     *  @brief Create a default grid map for a molecule.
     *  @tparam Args Argument types forwarded to create_default_grid_spec_map.
     *  @param mol    Molecule containing the atoms.
     *  @param scheme Pruning scheme to apply.
     *  @param bsz    Batch size for the generated grids.
     *  @param args   Additional arguments for grid specification.
     *  @return Map from atomic number to Grid instances.
     */
    template <typename... Args>
    inline static atomic_grid_map create_default_gridmap( 
      const Molecule& mol, PruningScheme scheme, BatchSize bsz,
      Args&&... args ) {

      return generate_gridmap( create_default_grid_spec_map(mol, scheme, 
        std::forward<Args>(args)...), bsz );

    }

    /**
     *  @brief Create a default MolGrid for a molecule.
     *  @tparam Args Argument types forwarded to create_default_gridmap.
     *  @param args Arguments passed to create_default_gridmap.
     *  @return Constructed MolGrid instance.
     */
    template <typename... Args>
    inline static MolGrid create_default_molgrid( Args&&... args ) {
      return MolGrid( create_default_gridmap(std::forward<Args>(args)...) );
    }

  };

}

