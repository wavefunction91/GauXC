#pragma once

#include "integrator_constants.hpp"
#include <gauxc/basisset.hpp>

namespace GauXC {

std::vector< std::pair<int32_t, int32_t> >
  gen_compressed_submat_map( const BasisSet<double>&       basis_set,
                             const std::vector< int32_t >& shell_mask ); 



}
