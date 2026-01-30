! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining types used throughout GauXC
module gauxc_types
  use iso_c_binding, only : c_ptr, c_int
  implicit none
  private

  public :: gauxc_type_molecule, gauxc_type_basisset, gauxc_type_molgrid, &
    & gauxc_type_runtime_environment, gauxc_type_load_balancer, &
    & gauxc_type_load_balancer_factory, gauxc_type_molecular_weights, &
    & gauxc_type_molecular_weights_factory, gauxc_type_functional, &
    & gauxc_type_integrator, gauxc_type_integrator_factory, gauxc_type_matrix

  enum, bind(c)
    enumerator :: gauxc_type_molecule = 1_c_int
    enumerator :: gauxc_type_basisset = 2_c_int
    enumerator :: gauxc_type_molgrid  = 3_c_int
    enumerator :: gauxc_type_runtime_environment = 4_c_int
    enumerator :: gauxc_type_load_balancer = 5_c_int
    enumerator :: gauxc_type_load_balancer_factory = 6_c_int
    enumerator :: gauxc_type_molecular_weights = 7_c_int
    enumerator :: gauxc_type_molecular_weights_factory = 8_c_int
    enumerator :: gauxc_type_functional = 9_c_int
    enumerator :: gauxc_type_integrator = 10_c_int
    enumerator :: gauxc_type_integrator_factory = 11_c_int
    enumerator :: gauxc_type_matrix = 12_c_int
  end enum

  !> @brief C interoperable header type
  type, bind(c), public :: gauxc_header_type
    integer(c_int) :: type
  end type

contains

end module gauxc_types
