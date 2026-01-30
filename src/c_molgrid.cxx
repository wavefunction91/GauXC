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
#include <gauxc/enums.h>
#include <gauxc/molecule.h>
#include <gauxc/molgrid.h>
#include <gauxc/molgrid.hpp>
#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/util/c_molecule.hpp>
#include <gauxc/util/c_molgrid.hpp>
#include <gauxc/util/c_status.hpp>

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
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return mg;
}

void gauxc_molgrid_delete(GauXCStatus* status, GauXCMolGrid* mg) {
  status->code = 0;
  if (mg == nullptr) return;
  if (mg->ptr != nullptr)
    delete detail::get_molgrid_ptr(*mg);
  mg->ptr = nullptr;
}

} // extern "C"
} // namespace GauXC::C