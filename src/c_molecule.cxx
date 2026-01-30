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
#include <gauxc/atom.h>
#include <gauxc/atom.hpp>
#include <gauxc/molecule.h>
#include <gauxc/molecule.hpp>
#include <gauxc/util/c_molecule.hpp>
#include <gauxc/util/c_status.hpp>

namespace GauXC::C {

extern "C" {

GauXCMolecule gauxc_molecule_new(GauXCStatus* status) {
  GauXCMolecule mol{};
  mol.hdr = GauXCHeader{GauXC_Type_Molecule};
  mol.ptr = nullptr;

  try {
    mol.ptr = new Molecule();
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return mol;
}

GauXCMolecule gauxc_molecule_new_from_atoms(GauXCStatus* status, GauXCAtom* atoms, size_t natoms) {
  GauXCMolecule mol{};
  mol.hdr = GauXCHeader{GauXC_Type_Molecule};
  mol.ptr = nullptr;
  Molecule* mol_ptr = nullptr;

  try {
    mol_ptr = new Molecule();
    mol_ptr->reserve(natoms);
    for (size_t i = 0; i < natoms; ++i) {
      mol_ptr->push_back(detail::convert_atom(atoms[i]));
    }
    mol.ptr = mol_ptr;
    status->code = 0;
  } catch (std::exception& e) {
    delete mol_ptr;
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return mol;
}

void gauxc_molecule_delete(GauXCStatus* status, GauXCMolecule* mol) {
  status->code = 0;
  if (mol == nullptr) return;
  if (mol->ptr != nullptr)
    delete detail::get_molecule_ptr(*mol);
  mol->ptr = nullptr;
}

size_t gauxc_molecule_natoms(GauXCStatus* status, const GauXCMolecule mol) {
  status->code = 0;
  if (mol.ptr == nullptr) return 0;
  return detail::get_molecule_ptr(mol)->natoms();
}

bool gauxc_molecule_equal(
  GauXCStatus* status,
  const GauXCMolecule mol_a,
  const GauXCMolecule mol_b
) {
  status->code = 0;
  if (mol_a.ptr == nullptr || mol_b.ptr == nullptr) return false;
  return *detail::get_molecule_ptr(mol_a) == *detail::get_molecule_ptr(mol_b);
}

} // extern "C"
} // namespace GauXC::C