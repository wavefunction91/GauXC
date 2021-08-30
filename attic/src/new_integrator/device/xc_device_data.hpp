#pragma once

#include <memory>
#include <iostream>

namespace GauXC {

template <typename F>
class XCDeviceData {

public:

  virtual void allocate_static_data( size_t _natoms,
                                     size_t _n_deriv,
                                     size_t _nbf,
                                     size_t _nshells ) = 0;

  virtual ~XCDeviceData() noexcept = default;

};

namespace integrator::device {

  template <typename T>
  std::shared_ptr<XCDeviceData<T>> make_device_data();

  extern template std::shared_ptr<XCDeviceData<double>> make_device_data();

}


}
