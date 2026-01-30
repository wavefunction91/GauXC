/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include <gauxc/shell.h>
#include <gauxc/shell.hpp>
#include <gauxc/basisset.h>
#include <gauxc/basisset.hpp>

namespace GauXC::detail {
static inline BasisSet<double>* get_basisset_ptr(C::GauXCBasisSet basis) noexcept {
  return static_cast<BasisSet<double>*>(basis.ptr);
}
static inline Shell<double> convert_shell(C::GauXCShell shell, bool normalize) noexcept {
  Shell<double>::prim_array alpha{};
  Shell<double>::prim_array coeff{};
  Shell<double>::cart_array O{0.0, 0.0, 0.0};

  for( int32_t i = 0; i < shell.nprim; ++i ) {
    alpha[i] = shell.exponents[i];
    coeff[i] = shell.coefficients[i];
  }
  for( int32_t i = 0; i < 3; ++i ) {
    O[i] = shell.origin[i];
  }

  Shell<double> sh{ PrimSize{shell.nprim},
                    AngularMomentum{shell.l},
                    SphericalType{shell.pure},
                    alpha, coeff, O, normalize };

  if (shell.shell_tolerance >= 0.0) {
    sh.set_shell_tolerance( shell.shell_tolerance );
  }

  return sh;
}
} // namespace GauXC::detail