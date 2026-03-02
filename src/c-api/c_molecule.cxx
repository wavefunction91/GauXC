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
#include <gauxc/c/atom.h>
#include <gauxc/c/molecule.h>

#include <gauxc/atom.hpp>
#include <gauxc/molecule.hpp>

#include "c_molecule.hpp"
#include "c_status.hpp"

namespace GauXC::C {

extern "C" {

GauXCMolecule gauxc_molecule_new(GauXCStatus* status) {
  detail::gauxc_status_init(status);
  GauXCMolecule mol{};
  mol.hdr = GauXCHeader{GauXC_Type_Molecule};
  mol.ptr = nullptr;

  try {
    mol.ptr = new Molecule();
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return mol;
}

GauXCMolecule gauxc_molecule_new_from_atoms(GauXCStatus* status, GauXCAtom* atoms, size_t natoms) {
  detail::gauxc_status_init(status);
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
  } catch (std::exception& e) {
    delete mol_ptr;
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return mol;
}

void gauxc_molecule_delete(GauXCStatus* status, GauXCMolecule* mol) {
  detail::gauxc_status_init(status);
  if (mol == nullptr) return;
  if (mol->ptr != nullptr)
    delete detail::get_molecule_ptr(*mol);
  mol->ptr = nullptr;
}

size_t gauxc_molecule_natoms(GauXCStatus* status, const GauXCMolecule mol) {
  detail::gauxc_status_init(status);
  if (mol.ptr == nullptr) return 0;
  return detail::get_molecule_ptr(mol)->natoms();
}

bool gauxc_molecule_equal(
  GauXCStatus* status,
  const GauXCMolecule mol_a,
  const GauXCMolecule mol_b
) {
  detail::gauxc_status_init(status);
  if (mol_a.ptr == nullptr || mol_b.ptr == nullptr) return false;
  return *detail::get_molecule_ptr(mol_a) == *detail::get_molecule_ptr(mol_b);
}

} // extern "C"
} // namespace GauXC::C