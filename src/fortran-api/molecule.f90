! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining molecular data structures for GauXC
module gauxc_molecule
  use iso_c_binding, only : c_ptr, c_size_t, c_bool, c_null_ptr
  use gauxc_atom, only : gauxc_atom_type
  use gauxc_status, only : gauxc_status_type
  use gauxc_types, only : gauxc_header_type, gauxc_type_molecule
  implicit none
  private

  public :: &
    & gauxc_molecule_new, &
    & gauxc_molecule_new_from_atoms, &
    & gauxc_molecule_delete, &
    & gauxc_molecule_natoms, &
    & gauxc_molecule_equal

  public :: &
    & gauxc_delete

  !> @brief C interoperable molecule type
  type, bind(c), public :: gauxc_molecule_type
    !> Header containing type information
    type(gauxc_header_type) :: hdr = gauxc_header_type(gauxc_type_molecule)
    !> Pointer to the internal molecule object
    type(c_ptr) :: ptr = c_null_ptr
  end type gauxc_molecule_type

  interface gauxc_molecule_new
    !> @brief Create a GauXC molecule object
    function gauxc_molecule_new(status) result(molecule) bind(c)
      import :: gauxc_status_type, gauxc_molecule_type
      implicit none
      !> @param status Status of the operation
      type(gauxc_status_type), intent(out) :: status
      !> @return Pointer to the newly created molecule object
      type(gauxc_molecule_type) :: molecule
    end function gauxc_molecule_new
  end interface gauxc_molecule_new

  interface gauxc_molecule_new
    !> @brief Create a new Molecule instance from an array of Atoms
    function gauxc_molecule_new_from_atoms(status, atoms, natoms) result(molecule) bind(c)
      import :: c_size_t, gauxc_status_type, gauxc_atom_type, gauxc_molecule_type
      implicit none
      !> @param status Status of the operation
      type(gauxc_status_type), intent(out) :: status
      !> @param atoms Pointer to an array of Atom objects
      type(gauxc_atom_type), intent(in) :: atoms(*)
      !> @param natoms Number of atoms in the array
      integer(c_size_t), value :: natoms
      !> @return Pointer to the newly created molecule object
      type(gauxc_molecule_type) :: molecule
    end function gauxc_molecule_new_from_atoms
  end interface gauxc_molecule_new

  interface gauxc_delete
    !> @brief Delete a GauXC molecule object
    subroutine gauxc_molecule_delete(status, molecule) bind(c)
      import :: c_ptr, gauxc_status_type, gauxc_molecule_type
      implicit none
      !> @param status Status of the operation
      type(gauxc_status_type), intent(out) :: status
      !> @param molecule Pointer to the molecule object to delete
      type(gauxc_molecule_type), intent(inout) :: molecule
    end subroutine gauxc_molecule_delete
  end interface gauxc_delete

  interface
    !> @brief Get the number of atoms in the molecule
    function gauxc_molecule_natoms(status, molecule) result(natoms) bind(c)
      import :: c_size_t, gauxc_status_type, gauxc_molecule_type
      implicit none
      !> @param status Status of the operation
      type(gauxc_status_type), intent(out) :: status
      !> @param molecule Pointer to the molecule object
      type(gauxc_molecule_type), value :: molecule
      !> @return Number of atoms in the molecule
      integer(c_size_t) :: natoms
    end function gauxc_molecule_natoms

    !> @brief Check if two Molecule instances are equal
    function gauxc_molecule_equal(status, mol_a, mol_b) result(is_equal) bind(c)
      import :: gauxc_status_type, gauxc_molecule_type, c_bool
      implicit none
      !> @param status Status of the operation
      type(gauxc_status_type), intent(out) :: status
      !> @param mol_a Pointer to the first molecule object
      type(gauxc_molecule_type), value :: mol_a
      !> @param mol_b Pointer to the second molecule object
      type(gauxc_molecule_type), value :: mol_b
      !> @return Logical indicating if the two molecules are equal
      logical(c_bool) :: is_equal
    end function gauxc_molecule_equal
  end interface

contains

end module gauxc_molecule
