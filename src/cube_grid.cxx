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
#include <gauxc/cube_grid.hpp>

#include <algorithm>

#include <gauxc/exceptions.hpp>

namespace GauXC {

CubeGrid CubeGrid::from_molecule( const Molecule& mol, int64_t nx, int64_t ny,
                                  int64_t nz, double margin ) {
  if( mol.empty() ) {
    GAUXC_GENERIC_EXCEPTION("CubeGrid::from_molecule: molecule has no atoms.");
  }
  if( nx < 1 or ny < 1 or nz < 1 ) {
    GAUXC_GENERIC_EXCEPTION("CubeGrid::from_molecule: nx, ny, nz must be >= 1.");
  }

  double xmin = mol[0].x, xmax = mol[0].x;
  double ymin = mol[0].y, ymax = mol[0].y;
  double zmin = mol[0].z, zmax = mol[0].z;
  for( const auto& a : mol ) {
    xmin = std::min(xmin, a.x);
    xmax = std::max(xmax, a.x);
    ymin = std::min(ymin, a.y);
    ymax = std::max(ymax, a.y);
    zmin = std::min(zmin, a.z);
    zmax = std::max(zmax, a.z);
  }

  CubeGrid grid;
  grid.origin = {xmin - margin, ymin - margin, zmin - margin};
  grid.nx = nx;
  grid.ny = ny;
  grid.nz = nz;
  const double ex = (xmax - xmin) + 2.0 * margin;
  const double ey = (ymax - ymin) + 2.0 * margin;
  const double ez = (zmax - zmin) + 2.0 * margin;
  grid.spacing[0] = nx > 1 ? ex / static_cast<double>(nx - 1) : 0.0;
  grid.spacing[1] = ny > 1 ? ey / static_cast<double>(ny - 1) : 0.0;
  grid.spacing[2] = nz > 1 ? ez / static_cast<double>(nz - 1) : 0.0;
  return grid;
}

std::vector<double> CubeGrid::points() const {
  std::vector<double> pts(static_cast<size_t>(num_points()) * 3);
  points_into(pts.data());
  return pts;
}

void CubeGrid::points_into( double* out ) const {
  size_t k = 0;
  for( int64_t ix = 0; ix < nx; ++ix ) {
    const double x = origin[0] + spacing[0] * static_cast<double>(ix);
    for( int64_t iy = 0; iy < ny; ++iy ) {
      const double y = origin[1] + spacing[1] * static_cast<double>(iy);
      for( int64_t iz = 0; iz < nz; ++iz ) {
        out[3 * k + 0] = x;
        out[3 * k + 1] = y;
        out[3 * k + 2] = origin[2] + spacing[2] * static_cast<double>(iz);
        ++k;
      }
    }
  }
}

}  // namespace GauXC
