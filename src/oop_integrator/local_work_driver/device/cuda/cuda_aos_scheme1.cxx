#include "cuda_aos_scheme1.hpp"
#include "kernels/collocation_device.hpp"

namespace GauXC {
namespace detail {

std::unique_ptr<XCDeviceData> CudaAoSScheme1::create_device_data() {
  return std::make_unique<Data>();
}


void CudaAoSScheme1::eval_collocation( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  size_t npts_max, nshells_max;
  for( auto& task : tasks ) {
    npts_max    = std::max( npts_max, task.npts );
    nshells_max = std::max( nshells_max, task.nshells );
  }

  eval_collocation_masked_combined( ntasks, npts_max, nshells_max,
    data->shells_device, data->device_tasks, *data->master_stream );

}

void CudaAoSScheme1::eval_collocation_gradient( XCDeviceData* _data ) {

  auto* data = dynamic_cast<Data*>(_data);
  if( !data ) throw std::runtime_error("BAD DATA CAST");

  auto tasks = data->host_device_tasks;
  const auto ntasks = tasks.size();

  size_t npts_max, nshells_max;
  for( auto& task : tasks ) {
    npts_max    = std::max( npts_max, task.npts );
    nshells_max = std::max( nshells_max, task.nshells );
  }

  eval_collocation_masked_combined_deriv1( ntasks, npts_max, nshells_max,
    data->shells_device, data->device_tasks, *data->master_stream );

}

}
}
