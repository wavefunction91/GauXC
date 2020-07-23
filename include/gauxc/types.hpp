#pragma once

#include <exchcxx/xc_functional.hpp>
#include <integratorxx/quadrature.hpp>
#include <integratorxx/batch/spherical_micro_batcher.hpp>

#include "named_type.hpp"

namespace GauXC {

using functional_type = ExchCXX::XCFunctional;
//using quadrature_type = IntegratorXX::QuadratureBase<
//  std::vector<std::array<double,3>>,
//  std::vector<double>
//>;
using quadrature_type = IntegratorXX::SphericalQuadratureBase<
  std::vector<std::array<double,3>>,
  std::vector<double>
>;

using batcher_type = IntegratorXX::SphericalMicroBatcher<
  typename quadrature_type::point_container,
  typename quadrature_type::weight_container
>;

}
