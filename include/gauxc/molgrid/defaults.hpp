#pragma once

#include <gauxc/molgrid.hpp>

namespace GauXC {

  atomic_scal_factor_map get_default_scaling_factors( RadialQuad rq, 
    AtomicNumber maxZ );

  atomic_grid_size_map get_default_grid_sizes( AtomicGridSizeDefault G, 
    AtomicNumber maxZ );

  atomic_scal_factor_map slater_radii_64();
  atomic_scal_factor_map slater_radii_30();
  atomic_scal_factor_map clementi_radii_67();

}

