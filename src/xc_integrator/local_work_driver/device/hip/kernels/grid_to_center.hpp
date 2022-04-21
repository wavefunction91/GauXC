#pragma once

namespace GauXC {

void compute_grid_to_center_dist( int32_t npts, int32_t natoms,
  const double* coords, const double* points_x,  const double* points_y, 
  const double* points_z, double* dist, int32_t lddist, hipStream_t stream );

}
