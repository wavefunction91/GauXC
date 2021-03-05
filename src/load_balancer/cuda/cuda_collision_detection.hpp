#pragma once

namespace GauXC      {
namespace load_balancer {
namespace cuda       {

struct CollisionDetectionCudaData {
    // Inputs
    double* low_points_device;
    double* high_points_device;
    double* centers_device;
    double* radii_device;
    size_t* shell_sizes_device;
    // Outputs
    int32_t position_list_length;
    int32_t* position_list_device;
    int32_t* counts_device;
    size_t*  nbe_list_device;
    // Intermediates
    int32_t* collisions_device;
    size_t temp_storage_bytes;
    void * temp_storage_device;
};

size_t compute_scratch( size_t ncubes, int32_t* counts_device );

void collision_detection( size_t ncubes,
                          size_t nspheres,
                          size_t LD_bit,
                          const double* low_points_device,
                          const double* high_points_device,
                          const double* centers_device,
                          const double* radii_device,
                                size_t  temp_storage_bytes,
                                 void * temp_storage_device,
                               int32_t* collisions_device, 
                               int32_t* counts_device,
                          cudaStream_t  stream);

void compute_position_list(size_t ncubes,
                           size_t nspheres,
                           size_t LD_bit,
                           const size_t* shell_sizes_device,
                           const int32_t* collisions_device,
                           const int32_t* counts_device,
                           int32_t* position_list_device,
                           size_t* nbe_list_device,
                          cudaStream_t  stream);

}
}
}

