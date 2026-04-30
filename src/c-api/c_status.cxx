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
#include <gauxc/c/status.h>

#include <cstdlib>

namespace GauXC::C {
extern "C" {

void gauxc_status_delete(GauXCStatus* status) {
  if(status == nullptr) return;
  if(status->message != nullptr) {
    std::free(status->message);
    status->message = nullptr;
  }
  status->code = 0;
}

} // extern "C"
} // namespace GauXC::C