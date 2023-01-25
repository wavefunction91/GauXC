/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#pragma once
#include <gauxc/gauxc_config.hpp>

namespace GauXC {

class RuntimeEnvironment;

#ifdef GAUXC_ENABLE_DEVICE
class DeviceRuntimeEnvironment;
class DeviceBackend;
#endif
}
