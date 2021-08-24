#pragma once

#include "integrator_constants.hpp"
#include <gauxc/basisset.hpp>

namespace GauXC      {
namespace integrator {

std::tuple< std::vector< std::array<int32_t, 3> > , std::vector< int32_t > >
  gen_compressed_submat_map( const BasisSet<double>&       basis_set,
                             const std::vector< int32_t >& shell_mask,
		             const int32_t LDA, const int32_t block_size ); 


}
}
