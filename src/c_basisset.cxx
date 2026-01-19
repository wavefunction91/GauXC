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
#include <gauxc/shell.h>
#include <gauxc/shell.hpp>
#include <gauxc/basisset.h>
#include <gauxc/basisset.hpp>
#include <gauxc/util/c_basisset.hpp>
#include <gauxc/util/c_status.hpp>
#include <exception>

namespace GauXC::C {
extern "C" {

GauXCBasisSet gauxc_basisset_new(GauXCStatus* status) {
  status->code = 0;
  GauXCBasisSet basis{};
  basis.hdr = GauXCHeader{GauXC_Type_BasisSet};
  basis.ptr = nullptr;
  try {
    basis.ptr = new BasisSet<double>();
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return basis;
}

GauXCBasisSet gauxc_basisset_new_from_shells(GauXCStatus* status, GauXCShell* shells, size_t nshells, bool normalize) {
  status->code = 0;
  GauXCBasisSet basis{};
  basis.hdr = GauXCHeader{GauXC_Type_BasisSet};
  basis.ptr = nullptr;
  if (nshells > detail::shell_nprim_max) {
    status->code = 1;
    status->message = detail::strdup("Number of primitives in shell exceeds maximum allowed");
    return basis;
  }
  BasisSet<double>* basis_ptr = nullptr;
  try {
    basis_ptr = new BasisSet<double>();
    basis_ptr->reserve(nshells);
    for (size_t i = 0; i < nshells; ++i) {
      basis_ptr->push_back(detail::convert_shell(shells[i], normalize));
    }
    basis.ptr = basis_ptr;
  } catch (std::exception& e) {
    delete basis_ptr;
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return basis;
}

void gauxc_basisset_delete(GauXCStatus* status, GauXCBasisSet* basis) {
  status->code = 0;
  if (basis == nullptr) return;
  if (basis->ptr != nullptr)
    delete detail::get_basisset_ptr(*basis);
  basis->ptr = nullptr;
}

} // extern "C"
} // namespace GauXC::C