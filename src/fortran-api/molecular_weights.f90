! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining molecular grid weights for GauXC
module gauxc_molecular_weights
  use iso_c_binding, only : c_ptr, c_null_ptr, c_bool, c_int, c_char, c_null_char
  use gauxc_enums, only : gauxc_xcweightalg
  use gauxc_status, only : gauxc_status_type
  use gauxc_types, only : gauxc_header_type, gauxc_type_molecular_weights, &
    & gauxc_type_molecular_weights_factory
  use gauxc_load_balancer, only : gauxc_load_balancer_type
  implicit none
  private

  public :: &
    & gauxc_molecular_weights_delete, &
    & gauxc_molecular_weights_modify_weights, &
    & gauxc_molecular_weights_factory_new, &
    & gauxc_molecular_weights_factory_delete, &
    & gauxc_molecular_weights_factory_get_instance

  public :: &
    & gauxc_delete, &
    & gauxc_get_instance

  !> @brief Settings for molecular grid weights
  type, bind(c), public :: gauxc_molecular_weights_settings
    integer(c_int) :: weight_alg = gauxc_xcweightalg%ssf
    logical(c_bool) :: becke_size_adjustment = .false.
  end type gauxc_molecular_weights_settings
   
  !> @brief GauXC MolecularWeights handle
  type, bind(c), public :: gauxc_molecular_weights_type
    type(gauxc_header_type) :: hdr = gauxc_header_type(gauxc_type_molecular_weights)
    type(c_ptr) :: ptr = c_null_ptr
  end type gauxc_molecular_weights_type

  !> @brief Create a MolecularWeightsFactory handle
  type, bind(c), public :: gauxc_molecular_weights_factory_type
    type(gauxc_header_type) :: hdr = gauxc_header_type(gauxc_type_molecular_weights_factory)
    type(c_ptr) :: ptr = c_null_ptr
  end type gauxc_molecular_weights_factory_type

  interface gauxc_delete
    !> @brief Delete a MolecularWeights instance.
    subroutine gauxc_molecular_weights_delete(status, mw) bind(c)
      import :: gauxc_status_type, gauxc_molecular_weights_type
      implicit none
      type(gauxc_status_type), intent(inout) :: status
      type(gauxc_molecular_weights_type), intent(inout) :: mw
    end subroutine gauxc_molecular_weights_delete
  end interface gauxc_delete

  interface
    !> @brief Apply molecular weights to a LoadBalancer's tasks.
    subroutine gauxc_molecular_weights_modify_weights(status, mw, lb) bind(c)
      import :: gauxc_status_type, gauxc_molecular_weights_type, gauxc_load_balancer_type
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param mw Handle to the MolecularWeights.
      type(gauxc_molecular_weights_type), value :: mw
      !> @param lb Handle to the LoadBalancer.
      type(gauxc_load_balancer_type), value :: lb
    end subroutine gauxc_molecular_weights_modify_weights

    !> @brief Create a new MolecularWeightsFactory instance.
    function gauxc_molecular_weights_factory_new_c(status, ex, local_work_kernel_name, settings) &
      & result(factory) bind(c, name="gauxc_molecular_weights_factory_new")
      import :: gauxc_status_type, gauxc_molecular_weights_factory_type, c_int, c_char, &
        & gauxc_molecular_weights_settings
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param ex Execution space.
      integer(c_int), value :: ex
      !> @param local_work_kernel_name Name of the LocalWorkDriver kernel to use.
      character(kind=c_char), intent(in) :: local_work_kernel_name(*)
      !> @param settings Settings for the MolecularWeights calculation.
      type(gauxc_molecular_weights_settings), value :: settings
      !> @return Handle to the created MolecularWeightsFactory.
      type(gauxc_molecular_weights_factory_type) :: factory
    end function gauxc_molecular_weights_factory_new_c
  end interface

  interface gauxc_delete
    !> @brief Delete a MolecularWeightsFactory instance.
    subroutine gauxc_molecular_weights_factory_delete(status, factory) bind(c)
      import :: gauxc_status_type, gauxc_molecular_weights_factory_type
      implicit none
      type(gauxc_status_type), intent(inout) :: status
      type(gauxc_molecular_weights_factory_type), intent(inout) :: factory
    end subroutine gauxc_molecular_weights_factory_delete
  end interface gauxc_delete

  interface gauxc_get_instance
    !> @brief Get MolecularWeights instance from a MolecularWeightsFactory.
    function gauxc_molecular_weights_factory_get_instance(status, factory) result(mw) bind(c)
      import :: gauxc_status_type, gauxc_molecular_weights_factory_type, gauxc_molecular_weights_type
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param factory Handle to the MolecularWeightsFactory.
      type(gauxc_molecular_weights_factory_type), value :: factory
      !> @return Handle to the created MolecularWeights.
      type(gauxc_molecular_weights_type) :: mw
    end function gauxc_molecular_weights_factory_get_instance
  end interface gauxc_get_instance

contains

  !> @brief Create a new MolecularWeightsFactory instance.
  function gauxc_molecular_weights_factory_new(status, ex, local_work_kernel_name, settings) &
    & result(factory)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param ex Execution space.
    integer(c_int), value :: ex
    !> @param local_work_kernel_name Name of the LocalWorkDriver kernel to use.
    character(kind=c_char, len=*), intent(in), optional :: local_work_kernel_name
    !> @param settings Settings for the MolecularWeights calculation.
    type(gauxc_molecular_weights_settings), intent(in), optional :: settings
    !> @return Handle to the created MolecularWeightsFactory.
    type(gauxc_molecular_weights_factory_type) :: factory

    type(gauxc_molecular_weights_settings) :: settings_
    character(kind=c_char, len=*), parameter :: default_kernel = "Default"
    character(kind=c_char), allocatable :: c_local_work_kernel_name(:)

    if (present(local_work_kernel_name)) then
      c_local_work_kernel_name = transfer(local_work_kernel_name // c_null_char, &
        & [character(kind=c_char) ::], len(local_work_kernel_name) + 1)
    else
      c_local_work_kernel_name = transfer(default_kernel // c_null_char, &
        & [character(kind=c_char) ::], len(default_kernel) + 1)
    end if

    if (present(settings)) then
      settings_ = settings
    else
      settings_ = gauxc_molecular_weights_settings()
    end if

    factory = gauxc_molecular_weights_factory_new_c(status, ex, &
      c_local_work_kernel_name, settings_)
  end function gauxc_molecular_weights_factory_new
end module gauxc_molecular_weights
