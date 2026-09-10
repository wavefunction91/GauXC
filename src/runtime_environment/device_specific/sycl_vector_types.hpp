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
#include <cstdint>

namespace GauXC {

/// Packed vector aggregates matching the layout / alignment of the CUDA and
/// HIP builtin vector types. sycl::vec is deliberately not used here: it pads
/// 3-component vectors out to 4 components, which would break the
/// reinterpret_cast reads of contiguous coordinate triples that the kernels
/// share with the CUDA backend.
struct alignas(16) double2 { double  x, y;    };
struct alignas(8)  double3 { double  x, y, z; };
struct alignas(4)  int3    { int32_t x, y, z; };

}
