#pragma once
#include <gauxc/gauxc_config.hpp>

namespace GauXC {

class RuntimeEnvironment;

#ifdef GAUXC_ENABLE_DEVICE
class DeviceRuntimeEnvironment;
class DeviceBackend;
#endif
}
