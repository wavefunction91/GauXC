#pragma once
#include <gauxc/gauxc_config.hpp>
#include <gauxc/shell.hpp>
#include <gauxc/enums.hpp>

namespace GauXC      {
namespace integrator {
namespace hip       {


void hip_reciprocal(size_t length, double* vec, hipStream_t stream); 

template <typename F>
void partition_weights_hip_SoA( XCWeightAlg    weight_alg,
                                 size_t         npts,
                                 size_t         LDatoms,
                                 size_t         natoms,
                                 const F*       points_device,
                                 const int32_t* iparent_device,
                                 const F*       dist_nearest_device,
                                 const F*       rab_device,
                                 const F*       atomic_coords_device,
                                       F*       weights_device,
                                       F*       dist_scratch_device,
                                 hipStream_t   stream );
                                 
                  
}
}
}
