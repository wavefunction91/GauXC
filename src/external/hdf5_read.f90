! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining HDF5 read functionality for GauXC
module gauxc_external_hdf5_read
  use iso_c_binding, only : c_ptr, c_char, c_null_char
  use gauxc_status, only : gauxc_status_type
  use gauxc_molecule, only : gauxc_molecule_type
  use gauxc_basisset, only : gauxc_basisset_type
  use gauxc_matrix, only : gauxc_matrix_type
  implicit none
  private

  public :: &
    & gauxc_molecule_read_hdf5_record, &
    & gauxc_basisset_read_hdf5_record, &
    & gauxc_matrix_read_hdf5_record

  public :: &
    & gauxc_read_hdf5_record

  interface
    !> @brief Read a GauXC molecule from an HDF5 file
    subroutine gauxc_molecule_read_hdf5_record_c(status, mol, fname, dset) &
      & bind(c, name="gauxc_molecule_read_hdf5_record")
      import :: gauxc_status_type, gauxc_molecule_type, c_char
      !> @param status GauXC status object
      type(gauxc_status_type), intent(out) :: status
      !> @param mol GauXC molecule object to populate
      type(gauxc_molecule_type), value :: mol
      !> @param fname Name of HDF5 file
      character(kind=c_char), intent(in) :: fname(*)
      !> @param dset Name of dataset within HDF5 file
      character(kind=c_char), intent(in) :: dset(*)
    end subroutine gauxc_molecule_read_hdf5_record_c

    !> @brief Read a GauXC basis set from an HDF5 file
    subroutine gauxc_basisset_read_hdf5_record_c(status, basis, fname, dset) &
      & bind(c, name="gauxc_basisset_read_hdf5_record")
      import :: gauxc_status_type, gauxc_basisset_type, c_char
      !> @param status GauXC status object
      type(gauxc_status_type), intent(out) :: status
      !> @param basis GauXC basis set object to populate
      type(gauxc_basisset_type), value :: basis
      !> @param fname Name of HDF5 file
      character(kind=c_char), intent(in) :: fname(*)
      !> @param dset Name of dataset within HDF5 file
      character(kind=c_char), intent(in) :: dset(*)
    end subroutine gauxc_basisset_read_hdf5_record_c

    !> @brief Read a GauXC matrix from an HDF5 file
    subroutine gauxc_matrix_read_hdf5_record_c(status, matrix, fname, dset) &
      & bind(c, name="gauxc_matrix_read_hdf5_record")
      import :: gauxc_status_type, gauxc_matrix_type, c_char
      !> @param status GauXC status object
      type(gauxc_status_type), intent(out) :: status
      !> @param matrix GauXC matrix object to populate
      type(gauxc_matrix_type), value :: matrix
      !> @param fname Name of HDF5 file
      character(kind=c_char), intent(in) :: fname(*)
      !> @param dset Name of dataset within HDF5 file
      character(kind=c_char), intent(in) :: dset(*)
    end subroutine gauxc_matrix_read_hdf5_record_c
  end interface

  interface gauxc_read_hdf5_record
    module procedure gauxc_molecule_read_hdf5_record
    module procedure gauxc_basisset_read_hdf5_record
    module procedure gauxc_matrix_read_hdf5_record
  end interface gauxc_read_hdf5_record

contains

  !> @brief Read a GauXC molecule from an HDF5 file
  subroutine gauxc_molecule_read_hdf5_record(status, mol, fname, dset)
    !> @param status GauXC status object
    type(gauxc_status_type), intent(out) :: status
    !> @param mol GauXC molecule object to populate
    type(gauxc_molecule_type), value :: mol
    !> @param fname Name of HDF5 file
    character(kind=c_char, len=*), intent(in) :: fname
    !> @param dset Name of dataset within HDF5 file
    character(kind=c_char, len=*), intent(in) :: dset

    character(kind=c_char), allocatable :: c_fname(:)
    character(kind=c_char), allocatable :: c_dset(:)

    c_fname = transfer(fname // c_null_char, [character(kind=c_char) ::], &
      & len(fname) + 1)
    c_dset = transfer(dset // c_null_char, [character(kind=c_char) ::], &
      & len(dset) + 1)

    call gauxc_molecule_read_hdf5_record_c(status, mol, c_fname, c_dset)
  end subroutine gauxc_molecule_read_hdf5_record

  !> @brief Read a GauXC basis set from an HDF5 file
  subroutine gauxc_basisset_read_hdf5_record(status, basis, fname, dset)
    !> @param status GauXC status object
    type(gauxc_status_type), intent(out) :: status
    !> @param basis GauXC basis set object to populate
    type(gauxc_basisset_type), value :: basis
    !> @param fname Name of HDF5 file
    character(kind=c_char, len=*), intent(in) :: fname
    !> @param dset Name of dataset within HDF5 file
    character(kind=c_char, len=*), intent(in) :: dset

    character(kind=c_char), allocatable :: c_fname(:)
    character(kind=c_char), allocatable :: c_dset(:)

    c_fname = transfer(fname // c_null_char, [character(kind=c_char) ::], &
      & len(fname) + 1)
    c_dset = transfer(dset // c_null_char, [character(kind=c_char) ::], &
      & len(dset) + 1)

    call gauxc_basisset_read_hdf5_record_c(status, basis, c_fname, c_dset)
  end subroutine gauxc_basisset_read_hdf5_record

  !> @brief Read a GauXC matrix from an HDF5 file
  subroutine gauxc_matrix_read_hdf5_record(status, matrix, fname, dset)
    !> @param status GauXC status object
    type(gauxc_status_type), intent(out) :: status
    !> @param matrix GauXC matrix object to populate
    type(gauxc_matrix_type), value :: matrix
    !> @param fname Name of HDF5 file
    character(kind=c_char, len=*), intent(in) :: fname
    !> @param dset Name of dataset within HDF5 file
    character(kind=c_char, len=*), intent(in) :: dset

    character(kind=c_char), allocatable :: c_fname(:)
    character(kind=c_char), allocatable :: c_dset(:)

    c_fname = transfer(fname // c_null_char, [character(kind=c_char) ::], &
      & len(fname) + 1)
    c_dset = transfer(dset // c_null_char, [character(kind=c_char) ::], &
      & len(dset) + 1)

    call gauxc_matrix_read_hdf5_record_c(status, matrix, c_fname, c_dset)
  end subroutine gauxc_matrix_read_hdf5_record
end module gauxc_external_hdf5_read
