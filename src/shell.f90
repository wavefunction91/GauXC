! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module for shell-related data structures and routines
module gauxc_shell
  use iso_c_binding, only : c_ptr, c_int32_t, c_double, c_bool
  implicit none
  private

  !> @brief Number of maximum primitives per shell
  integer, parameter :: max_primitives_per_shell = 32

  !> @brief Default shell tolerance
  real(c_double), parameter :: default_shell_tolerance = 1.0e-10_c_double

  !> @brief Data structure representing a shell
  type, bind(c), public :: gauxc_shell_type
    !> Angular momentum quantum number
    integer(c_int32_t) :: l
    !> Is this a pure shell
    logical(c_bool) :: pure
    !> Number of primitives in the shell
    integer(c_int32_t) :: nprim
    !> Exponents of the primitives
    real(c_double) :: exponents(max_primitives_per_shell)
    !> Coefficients of the primitives
    real(c_double) :: coefficients(max_primitives_per_shell)
    !> Origin of the shell (x, y, z)
    real(c_double) :: origin(3)
    !> Shell tolerance
    real(c_double) :: shell_tolerance = default_shell_tolerance
  end type gauxc_shell_type

contains

end module gauxc_shell