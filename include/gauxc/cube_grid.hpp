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

#include <array>
#include <cstdint>
#include <vector>

#include <gauxc/molecule.hpp>

namespace GauXC {

/** @brief Specification of an axis-aligned 3D rectangular point grid.
 *
 *  The grid axes are parallel to the Cartesian frame, i.e. the axis vectors
 *  are (spacing[0], 0, 0), (0, spacing[1], 0) and (0, 0, spacing[2]). All
 *  quantities are in atomic units (Bohr).
 *
 *  Total number of points: nx*ny*nz. Storage cost per scalar field:
 *  8*nx*ny*nz bytes (double precision).
 */
struct CubeGrid {
  std::array<double, 3> origin{0.0, 0.0, 0.0};   ///< (0,0,0) corner, Bohr
  std::array<double, 3> spacing{0.2, 0.2, 0.2};  ///< Step on each axis, Bohr
  int64_t nx = 80;
  int64_t ny = 80;
  int64_t nz = 80;

  /// Total number of grid points.
  int64_t num_points() const noexcept { return nx * ny * nz; }

  /** @brief Build a default grid that tightly encloses a molecule.
   *
   *  The bounding box of the atomic centres is extended by `margin` Bohr on
   *  each side and discretised with the requested number of points along
   *  each axis. Spacing is chosen so that the first and last grid points
   *  coincide with the extended bounding-box corners (PySCF cubegen
   *  convention).
   *
   *  The margin is the same on all three axes, so a flat molecule gets a box
   *  that is thin in the direction it is flat in. An orbital that reaches out
   *  that way, such as the pi system of an aromatic ring, can still be large
   *  where the box ends and will look cut off when plotted. Use a bigger
   *  margin, or set origin, spacing and point counts yourself, if that matters.
   *
   *  @param mol     Molecule whose atomic centres define the bounding box.
   *  @param nx,ny,nz Number of grid points along each axis.
   *  @param margin  Margin (Bohr) added on each side. Default 3.0 matches
   *                 the PySCF cubegen default.
   */
  static CubeGrid from_molecule( const Molecule& mol,
                                 int64_t nx = 80, int64_t ny = 80,
                                 int64_t nz = 80, double margin = 3.0 );

  /** @brief Materialise the grid points as an AoS coordinate array.
   *
   *  Returns a vector of length 3*num_points() laid out in row-major
   *  (ix, iy, iz) order with iz varying fastest, matching the field-element
   *  ordering expected by `write_cube`. Suitable to be passed directly as
   *  the `points` argument to `OrbitalEvaluator::eval_orbital` /
   *  `eval_density`.
   */
  std::vector<double> points() const;

  /** @brief Write grid-point coordinates into a caller-supplied buffer.
   *
   *  Same layout as `points()` but writes into `out`, which must have room
   *  for at least 3*num_points() doubles. Use this to avoid a heap
   *  allocation when the caller already owns a suitably sized buffer.
   */
  void points_into( double* out ) const;
};

}  // namespace GauXC
