/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#include "host_reduction_driver.hpp"

namespace GauXC {

HostReductionDriver::HostReductionDriver(const RuntimeEnvironment& rt) :
  detail::ReductionDriverImpl(rt) { }


HostReductionDriver::~HostReductionDriver() noexcept = default;



bool HostReductionDriver::takes_host_memory() const {return true;}; 
bool HostReductionDriver::takes_device_memory() const {return false;};

}
