! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining atomic data structures and routines for GauXC
module gauxc_atom
  use iso_c_binding, only : c_int64_t, c_double
  implicit none
  private

  !> @brief Data structure representing an atom
  type, bind(c), public :: gauxc_atom_type
    !> Atomic number
    integer(c_int64_t) :: atomic_number
    !> Cartesian x coordinate
    real(c_double) :: x
    !> Cartesian y coordinate
    real(c_double) :: y
    !> Cartesian z coordinate
    real(c_double) :: z 
  end type gauxc_atom_type

contains

end module gauxc_atom