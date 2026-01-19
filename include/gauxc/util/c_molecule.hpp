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
#include <gauxc/molecule.h>
#include <gauxc/molecule.hpp>

namespace GauXC::detail {
static inline Molecule* get_molecule_ptr(C::GauXCMolecule mol) noexcept {
  return static_cast<Molecule*>(mol.ptr);
}
static inline Atom convert_atom(C::GauXCAtom atom) noexcept {
  return Atom{AtomicNumber(atom.Z), atom.x, atom.y, atom.z };
}
}