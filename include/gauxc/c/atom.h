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
#include <cstddef>
#else
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#endif

#ifdef __cplusplus
extern "C" {
namespace GauXC::C {
#endif

/**
 * @brief GauXC C API Atom representation.
 */
typedef struct GauXCAtom {
  int64_t Z;   ///< Atomic number.
  double  x;   ///< X coordinate (Bohr).
  double  y;   ///< Y coordinate (Bohr).
  double  z;   ///< Z coordinate (Bohr).
} GauXCAtom;

#ifdef __cplusplus
} // namespace GauXC::C
} // extern "C"
#endif
