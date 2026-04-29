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
#include <cstdlib>

#include <gauxc/c/status.h>

#include <gauxc/exceptions.hpp>

namespace GauXC::detail {

static inline void gauxc_status_init(C::GauXCStatus* status) {
   if (status != nullptr) {
     status->code = 0;
     std::free(status->message);
     status->message = nullptr;
   }
}

static inline void gauxc_status_handle(C::GauXCStatus* status, int code, const char* message) {
   if (status != nullptr) {
     status->code = code;
     if (status->message != nullptr) {
       std::free(status->message);
     }
     status->message = (char*)std::malloc(std::strlen(message) + 1);
     std::strcpy(status->message, message);
   } else {
     GAUXC_GENERIC_EXCEPTION(message);
   }
}

}