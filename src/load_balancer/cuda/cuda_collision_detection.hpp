#pragma once

namespace GauXC      {
namespace integrator {
namespace cuda       {

void collision_detection( size_t ncubes,
                          size_t nspheres,
                          size_t LD_bit,
                          const double* low_points,
                          const double* high_points,
                          const double* centers,
                          const double* radii,
                                int32_t* counts,
                                int32_t** position_list);


}
}
}

