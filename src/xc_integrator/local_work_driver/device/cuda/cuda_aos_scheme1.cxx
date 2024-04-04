/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "cuda_aos_scheme1.hpp"
#include "device/cuda/cuda_backend.hpp"
#include "cuda_aos_scheme1_weights.hpp"

namespace GauXC {

template <typename Base>
std::unique_ptr<XCDeviceData> CudaAoSScheme1<Base>::create_device_data(const DeviceRuntimeEnvironment& rt) {
  return std::make_unique<Data>(rt);
}

template <typename Base> 
void CudaAoSScheme1<Base>::partition_weights( XCDeviceData* _data ) {
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) GAUXC_BAD_LWD_DATA_CAST();

  auto device_backend = dynamic_cast<CUDABackend*>(data->device_backend_);
  if( !device_backend ) GAUXC_BAD_BACKEND_CAST();


  // Compute distances from grid to atomic centers
  const auto ldatoms = data->get_ldatoms();
  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;
  auto scheme1_stack = data->scheme1_stack;
  cuda_aos_scheme1_weights_wrapper( data->total_npts_task_batch, data->global_dims.natoms,
    base_stack.points_x_device, base_stack.points_y_device, base_stack.points_z_device,
    static_stack.rab_device, ldatoms, static_stack.coords_device, 
    scheme1_stack.dist_scratch_device, ldatoms, scheme1_stack.iparent_device, 
    scheme1_stack.dist_nearest_device, base_stack.weights_device, *device_backend->master_stream );
}


template struct CudaAoSScheme1<AoSScheme1Base>;
#ifdef GAUXC_HAS_MAGMA
template struct CudaAoSScheme1<AoSScheme1MAGMABase>;
#endif
#ifdef GAUXC_HAS_CUTLASS
template struct CudaAoSScheme1<AoSScheme1CUTLASSBase>;
#endif



}
