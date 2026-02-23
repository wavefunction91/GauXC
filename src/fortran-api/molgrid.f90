! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining molecular grid functionality for GauXC
module gauxc_molgrid
  use iso_c_binding, only : c_ptr, c_null_ptr, c_int64_t, c_int
  use gauxc_types, only : gauxc_header_type, gauxc_type_molgrid
  use gauxc_status, only : gauxc_status_type
  use gauxc_molecule, only : gauxc_molecule_type
  implicit none
  private

  public :: &
    & gauxc_molgrid_new_default, &
    & gauxc_molgrid_delete

  public :: &
    & gauxc_delete

  !> @brief GauXC C API Molecular Grid handle.
  type, bind(c), public :: gauxc_molgrid_type
    !> @brief Header for internal use.
    type(gauxc_header_type) :: hdr = gauxc_header_type(gauxc_type_molgrid)
    !> @brief Pointer to the Molecular Grid instance.
    type(c_ptr) :: ptr = c_null_ptr
  end type gauxc_molgrid_type

  interface
    !> @brief Create a new MolGrid instance from atomic grids
    function gauxc_molgrid_new_default(status, mol, pruning_scheme, &
      & batchsize, radial_quad, grid_size) result(molgrid) bind(c)
      import :: gauxc_status_type, gauxc_molgrid_type, gauxc_molecule_type, &
        & c_int, c_int64_t
      implicit none
      !> @param status GauXC status object
      type(gauxc_status_type), intent(out) :: status
      !> @param mol Molecule to generate grid for
      type(gauxc_molecule_type), value :: mol
      !> @param pruning_scheme Pruning scheme to use
      integer(c_int), value :: pruning_scheme
      !> @param batchsize Batch size for grid generation
      integer(c_int64_t), value :: batchsize
      !> @param radial_quad Radial quadrature to use
      integer(c_int), value :: radial_quad
      !> @param grid_size Estimated grid size
      integer(c_int), value :: grid_size
      !> @return MolGrid handle
      type(gauxc_molgrid_type) :: molgrid
    end function gauxc_molgrid_new_default
  end interface

  interface gauxc_delete
    !> @brief Delete a MolGrid instance
    subroutine gauxc_molgrid_delete(status, molgrid) bind(c)
      import :: gauxc_status_type, gauxc_molgrid_type
      implicit none
      !> @param status GauXC status object
      type(gauxc_status_type), intent(out) :: status
      !> @param molgrid MolGrid handle to delete
      type(gauxc_molgrid_type), intent(inout) :: molgrid
    end subroutine gauxc_molgrid_delete
  end interface gauxc_delete

contains

end module gauxc_molgrid
