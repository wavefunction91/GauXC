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
#include <gauxc/c/enums.h>
#include <gauxc/c/molecule.h>
#include <gauxc/c/molgrid.h>

#include <gauxc/molgrid.hpp>
#include <gauxc/molgrid/defaults.hpp>

#include "c_molecule.hpp"
#include "c_molgrid.hpp"
#include "c_status.hpp"

namespace GauXC::C {
extern "C" {

GauXCMolGrid gauxc_molgrid_new_default(
  GauXCStatus* status,
  const GauXCMolecule mol,
  const enum GauXC_PruningScheme pruning_scheme,
  const int64_t batchsize,
  const enum GauXC_RadialQuad radial_quad,
  const enum GauXC_AtomicGridSizeDefault grid_size
) {
  detail::gauxc_status_init(status);
  GauXCMolGrid mg{};
  mg.hdr = GauXCHeader{GauXC_Type_MolGrid};
  mg.ptr = nullptr;

  try {
    auto grid_map = MolGridFactory::create_default_gridmap(
      *detail::get_molecule_ptr(mol),
      static_cast<PruningScheme>(pruning_scheme),
      static_cast<BatchSize>(batchsize),
      static_cast<RadialQuad>(radial_quad),
      static_cast<AtomicGridSizeDefault>(grid_size)
    );
    mg.ptr = new MolGrid(grid_map);
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return mg;
}

void gauxc_molgrid_delete(GauXCStatus* status, GauXCMolGrid* mg) {
  detail::gauxc_status_init(status);
  if (mg == nullptr) return;
  if (mg->ptr != nullptr)
    delete detail::get_molgrid_ptr(*mg);
  mg->ptr = nullptr;
}

} // extern "C"
} // namespace GauXC::C