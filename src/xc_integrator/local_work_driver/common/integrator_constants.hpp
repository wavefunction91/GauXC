#pragma once

namespace GauXC      {
namespace integrator {

template <typename F = double>
constexpr F magic_ssf_factor = 0.64;

constexpr double ssf_weight_tol = 1e-10;

}
}
