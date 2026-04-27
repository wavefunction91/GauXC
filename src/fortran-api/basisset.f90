! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining basis set functionality for GauXC
module gauxc_basisset
  use iso_c_binding, only : c_ptr, c_null_ptr, c_size_t, c_bool
  use gauxc_status, only : gauxc_status_type
  use gauxc_types, only : gauxc_header_type, gauxc_type_basisset
  use gauxc_shell, only : gauxc_shell_type
  implicit none
  private

  public :: &
    & gauxc_basisset_new, &
    & gauxc_basisset_new_from_shells, &
    & gauxc_basisset_delete

  public :: &
    & gauxc_delete

  !> @brief C interoperable basis set type
  type, bind(c), public :: gauxc_basisset_type
    !> Header containing type information
    type(gauxc_header_type) :: hdr = gauxc_header_type(gauxc_type_basisset)
    !> Pointer to the internal basis set object
    type(c_ptr) :: ptr = c_null_ptr
  end type gauxc_basisset_type

  interface
    !> @brief Create new basis set object
    function gauxc_basisset_new(status) result(basis) bind(c)
      import :: gauxc_basisset_type, gauxc_status_type
      implicit none
      !> @param status Status of the operation
      type(gauxc_status_type), intent(inout) :: status
      !> @return Pointer to the newly created basis set object
      type(gauxc_basisset_type) :: basis
    end function gauxc_basisset_new
  end interface

  interface gauxc_basisset_new_from_shells
    !> @brief Create a new BasisSet instance from an array of Shells
    function gauxc_basisset_new_from_shells_c(&
      status, shells, nshells, normalize) result(basis) &

      & bind(c, name="gauxc_basisset_new_from_shells")
      import :: c_size_t, gauxc_status_type, gauxc_shell_type, &
                gauxc_basisset_type, c_bool
      implicit none
      !> @param status Status of the operation
      type(gauxc_status_type), intent(inout) :: status
      !> @param shells Pointer to an array of Shell objects
      type(gauxc_shell_type), intent(in) :: shells(*)
      !> @param nshells Number of shells in the array
      integer(c_size_t), value :: nshells
      !> @param normalize Whether to normalize the basis functions
      logical(c_bool), value :: normalize
      !> @return Pointer to the newly created basis set object
      type(gauxc_basisset_type) :: basis
    end function gauxc_basisset_new_from_shells_c
    module procedure gauxc_basisset_new_from_shells
  end interface gauxc_basisset_new_from_shells

  interface gauxc_delete
    !> @brief Delete a GauXC basis set object
    subroutine gauxc_basisset_delete(status, basis) bind(c)
      import :: gauxc_status_type, gauxc_basisset_type
      implicit none
      !> @param status Status of the operation
      type(gauxc_status_type), intent(inout) :: status
      !> @param basis Pointer to the basis set object to delete
      type(gauxc_basisset_type), intent(inout) :: basis
    end subroutine gauxc_basisset_delete
  end interface gauxc_delete

contains

  !> @brief Create a new BasisSet instance from an array of Shells
  function gauxc_basisset_new_from_shells(status, shells, normalize) &
    result(basis)
    !> @param status Status of the operation
    type(gauxc_status_type), intent(inout) :: status
    !> @param shells Pointer to an array of Shell objects
    type(gauxc_shell_type), intent(in) :: shells(:)
    !> @param normalize Whether to normalize the basis functions
    logical, intent(in), optional :: normalize
    !> @return Pointer to the newly created basis set object
    type(gauxc_basisset_type) :: basis

    logical :: normalize_resolved

    normalize_resolved = .true.
    if (present(normalize)) normalize_resolved = normalize

    basis = gauxc_basisset_new_from_shells_c( &
      status, shells, size(shells, kind=c_size_t), &
      logical(normalize, kind=c_bool))
  end function gauxc_basisset_new_from_shells

end module gauxc_basisset
