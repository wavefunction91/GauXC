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

#ifdef __cplusplus
#include <cstdint>
#include <cstdbool>
#else
#include <stdint.h>
#include <stdbool.h>
#endif

#ifdef __cplusplus
extern "C" {
namespace GauXC::C {
#endif

/**
 * @brief GauXC C API shell representation.
 */
typedef struct GauXCShell {
  int32_t l;        ///< Angular momentum.
  bool    pure;     ///< Spherical (true) or Cartesian (false) functions.
  int32_t nprim;    ///< Number of primitives.
  double  exponents[32];    ///< Pointer to array of primitive exponents.
  double  coefficients[32]; ///< Pointer to array of contraction coefficients.
  double  origin[3];        ///< Shell origin.
  double  shell_tolerance;  ///< Shell tolerance.
} GauXCShell;

#ifdef __cplusplus
} // namespace GauXC::C
} // extern "C"
#endif