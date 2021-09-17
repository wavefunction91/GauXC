#include "hip_aos_scheme1.hpp"
#include "device/hip/hip_backend.hpp"
#include "kernels/grid_to_center.hpp"
#include "kernels/hip_ssf_1d.hpp"

namespace GauXC {

template <typename Base>
std::unique_ptr<XCDeviceData> HipAoSScheme1<Base>::create_device_data() {
  return std::make_unique<Data>();
}


template <typename Base>
void HipAoSScheme1<Base>::partition_weights( XCDeviceData* _data ) {
  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto device_backend = dynamic_cast<HIPBackend*>(data->device_backend_.get());
  if( !device_backend ) throw std::runtime_error("BAD BACKEND CAST");

  const auto ldatoms = data->get_ldatoms();
  auto base_stack    = data->base_stack;
  auto static_stack  = data->static_stack;
  auto scheme1_stack = data->scheme1_stack;

  // Compute distances from grid to atomic centers
  compute_grid_to_center_dist( data->total_npts_task_batch, data->global_dims.natoms,
    static_stack.coords_device, base_stack.points_device, 
    scheme1_stack.dist_scratch_device, ldatoms, *device_backend->master_stream );

  // Modify weights
  partition_weights_ssf_1d( data->total_npts_task_batch, data->global_dims.natoms,
    static_stack.rab_device, ldatoms, static_stack.coords_device, 
    scheme1_stack.dist_scratch_device, ldatoms, scheme1_stack.iparent_device, 
    scheme1_stack.dist_nearest_device, base_stack.weights_device,
    *device_backend->master_stream );

}


template struct HipAoSScheme1<AoSScheme1Base>;
#ifdef GAUXC_ENABLE_MAGMA
template struct HipAoSScheme1<AoSScheme1MAGMABase>;
#endif


}
