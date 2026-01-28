! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining load balancing functionality for GauXC
module gauxc_load_balancer
  use iso_c_binding, only : c_ptr, c_bool, c_int, c_char, c_null_char, c_null_ptr
  use gauxc_status, only : gauxc_status_type
  use gauxc_types, only : gauxc_header_type, gauxc_type_load_balancer, &
    & gauxc_type_load_balancer_factory
  use gauxc_runtime_environment, only : gauxc_runtime_environment_type
  use gauxc_molecule, only : gauxc_molecule_type
  use gauxc_basisset, only : gauxc_basisset_type
  use gauxc_molgrid, only : gauxc_molgrid_type
  implicit none
  private

  public :: &
    & gauxc_load_balancer_factory_new, &
    & gauxc_load_balancer_factory_delete, &
    & gauxc_load_balancer_factory_get_instance, &
    & gauxc_load_balancer_factory_get_shared_instance, &
    & gauxc_load_balancer_delete

  public :: &
    & gauxc_delete, &
    & gauxc_get_instance, &
    & gauxc_get_shared_instance

  type, bind(c), public :: gauxc_load_balancer_type
    type(gauxc_header_type) :: hdr = gauxc_header_type(gauxc_type_load_balancer)
    type(c_ptr) :: ptr = c_null_ptr
    logical(c_bool) :: owned = .true.
  end type gauxc_load_balancer_type

  type, bind(c), public :: gauxc_load_balancer_factory_type
    type(gauxc_header_type) :: hdr = gauxc_header_type(gauxc_type_load_balancer_factory)
    type(c_ptr) :: ptr = c_null_ptr
  end type gauxc_load_balancer_factory_type

  interface
    !> @brief Create a new LoadBalancerFactory instance.
    function gauxc_load_balancer_factory_new_c(status, ex, kernel_name) result(lbf) &
      & bind(c, name="gauxc_load_balancer_factory_new")
      import :: gauxc_status_type, gauxc_load_balancer_factory_type, c_int, c_char
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param ex Execution space.
      integer(c_int), value :: ex
      !> @param kernel_name Name of the load balancing kernel to use.
      character(kind=c_char), intent(in) :: kernel_name(*)
      !> @return Handle to the created LoadBalancerFactory.
      type(gauxc_load_balancer_factory_type) :: lbf
    end function gauxc_load_balancer_factory_new_c
  end interface

  interface gauxc_delete
    !> @brief Delete a LoadBalancerFactory instance.
    subroutine gauxc_load_balancer_factory_delete(status, factory) bind(c)
      import :: gauxc_status_type, gauxc_load_balancer_factory_type
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param factory Handle to the LoadBalancerFactory to delete.
      type(gauxc_load_balancer_factory_type), intent(inout) :: factory
    end subroutine gauxc_load_balancer_factory_delete
  end interface gauxc_delete

  interface gauxc_get_instance
    !> @brief Create a new LoadBalancer instance from a LoadBalancerFactory.
    function gauxc_load_balancer_factory_get_instance( &
      status, factory, env, mol, mg, basis) result(lb) bind(c)
      import :: gauxc_status_type, gauxc_load_balancer_factory_type, &
        & gauxc_load_balancer_type, gauxc_molecule_type, gauxc_molgrid_type, &
        & gauxc_basisset_type, gauxc_runtime_environment_type, c_ptr
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param factory Handle to the LoadBalancerFactory.
      type(gauxc_load_balancer_factory_type), value :: factory
      !> @param env Handle to the RuntimeEnvironment.
      type(gauxc_runtime_environment_type), value :: env
      !> @param mol Handle to the Molecule.
      type(gauxc_molecule_type), value :: mol
      !> @param mg Handle to the MolGrid.
      type(gauxc_molgrid_type), value :: mg
      !> @param basis Handle to the BasisSet.
      type(gauxc_basisset_type), value :: basis
      !> @return Handle to the created LoadBalancer.
      type(gauxc_load_balancer_type) :: lb
    end function gauxc_load_balancer_factory_get_instance
  end interface gauxc_get_instance

  interface gauxc_get_shared_instance
    !> @brief Create a shared LoadBalancer instance from a LoadBalancerFactory.
    function gauxc_load_balancer_factory_get_shared_instance( &
      status, factory, env, mol, mg, basis) result(lb) bind(c)
      import :: gauxc_status_type, gauxc_load_balancer_factory_type, &
        & gauxc_load_balancer_type, gauxc_molecule_type, gauxc_molgrid_type, &
        & gauxc_basisset_type, gauxc_runtime_environment_type, c_ptr
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param factory Handle to the LoadBalancerFactory.
      type(gauxc_load_balancer_factory_type), value :: factory
      !> @param env Handle to the RuntimeEnvironment.
      type(gauxc_runtime_environment_type), value :: env
      !> @param mol Handle to the Molecule.
      type(gauxc_molecule_type), value :: mol
      !> @param mg Handle to the MolGrid.
      type(gauxc_molgrid_type), value :: mg
      !> @param basis Handle to the BasisSet.
      type(gauxc_basisset_type), value :: basis
      !> @return Handle to the created LoadBalancer.
      type(gauxc_load_balancer_type) :: lb
    end function gauxc_load_balancer_factory_get_shared_instance
  end interface gauxc_get_shared_instance

  interface gauxc_delete
    !> @brief Delete a LoadBalancer instance.
    subroutine gauxc_load_balancer_delete(status, lb) bind(c)
      import :: gauxc_status_type, gauxc_load_balancer_type
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param lb Handle to the LoadBalancer to delete.
      type(gauxc_load_balancer_type), intent(inout) :: lb
    end subroutine gauxc_load_balancer_delete
  end interface gauxc_delete

contains

  !> @brief Create a new LoadBalancerFactory instance.
  function gauxc_load_balancer_factory_new(status, ex, kernel_name) result(lbf)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param ex Execution space.
    integer(c_int), value :: ex
    !> @param kernel_name Name of the load balancing kernel to use.
    character(kind=c_char, len=*), intent(in) :: kernel_name
    !> @return Handle to the created LoadBalancerFactory.
    type(gauxc_load_balancer_factory_type) :: lbf

    character(kind=c_char), allocatable :: c_kernel_name(:)

    c_kernel_name = transfer(kernel_name // c_null_char, &
      & [character(kind=c_char) ::], len(kernel_name) + 1)

    lbf = gauxc_load_balancer_factory_new_c(status, ex, c_kernel_name)
  end function gauxc_load_balancer_factory_new
end module gauxc_load_balancer
