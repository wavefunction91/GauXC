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
#include <exception>

#include <gauxc/c/shell.h>
#include <gauxc/c/basisset.h>

#include <gauxc/shell.hpp>
#include <gauxc/basisset.hpp>

#include "c_basisset.hpp"
#include "c_status.hpp"

namespace GauXC::C {
extern "C" {

GauXCBasisSet gauxc_basisset_new(GauXCStatus* status) {
  detail::gauxc_status_init(status);
  GauXCBasisSet basis{};
  basis.hdr = GauXCHeader{GauXC_Type_BasisSet};
  basis.ptr = nullptr;
  try {
    basis.ptr = new BasisSet<double>();
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return basis;
}

GauXCBasisSet gauxc_basisset_new_from_shells(GauXCStatus* status, const GauXCShell* shells, size_t nshells, bool normalize) {
  detail::gauxc_status_init(status);
  GauXCBasisSet basis{};
  basis.hdr = GauXCHeader{GauXC_Type_BasisSet};
  basis.ptr = nullptr;
  if (shells == nullptr || nshells == 0) {
    detail::gauxc_status_handle(status, 1, "Shell list is null or empty");
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
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return basis;
}

void gauxc_basisset_delete(GauXCStatus* status, GauXCBasisSet* basis) {
  detail::gauxc_status_init(status);
  if (basis == nullptr) return;
  if (basis->ptr != nullptr)
    delete detail::get_basisset_ptr(*basis);
  basis->ptr = nullptr;
}

} // extern "C"
} // namespace GauXC::C