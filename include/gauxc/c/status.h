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
extern "C" {
namespace GauXC::C {
#endif

/**
 * @brief Status object for GauXC C API functions.
 *
 * Important, always initialize the status object before passing to any C API function.
 */
typedef struct GauXCStatus {
  int code;
  char* message;
} GauXCStatus;

/**
 * @brief Delete a GauXCStatus object.
 * @param status Pointer to the GauXCStatus object to delete.
 */
extern void gauxc_status_delete(GauXCStatus* status);

#ifdef __cplusplus
} // namespace GauXC::C
} // extern "C"
#endif