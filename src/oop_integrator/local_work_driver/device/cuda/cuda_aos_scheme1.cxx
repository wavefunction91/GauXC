#include "cuda_aos_scheme1.hpp"

namespace GauXC {
namespace detail {

std::unique_ptr<XCDeviceData> CudaAoSScheme1::create_device_data() {
  return std::make_unique<Data>();
}

}
}
