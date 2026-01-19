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

#include <cstring>

namespace GauXC::detail {

static inline char* strdup(const char* str) {
  if(str == nullptr) return nullptr;
  size_t len = std::strlen(str);
  char* copy = new char[len + 1];
  std::strcpy(copy, str);
  return copy;
}

}