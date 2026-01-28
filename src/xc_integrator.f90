! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining XC integrator functionality for GauXC
module gauxc_integrator
  use iso_c_binding, only : c_ptr, c_null_ptr, c_bool, c_int, c_char, c_double, c_null_char
  use gauxc_status, only : gauxc_status_type
  use gauxc_types, only : gauxc_header_type, gauxc_type_integrator, &
    & gauxc_type_integrator_factory
  use gauxc_xc_functional, only : gauxc_functional_type
  use gauxc_load_balancer, only : gauxc_load_balancer_type
  use gauxc_matrix, only : gauxc_matrix_type
  implicit none
  private

  public :: &
    & gauxc_integrator_factory_new, &
    & gauxc_integrator_factory_delete, &
    & gauxc_integrator_factory_get_instance, &
    & gauxc_integrator_factory_get_shared_instance, &
    & gauxc_integrator_delete, &
    & gauxc_integrator_integrate_den, &
    & gauxc_integrator_eval_exc_rks, &
    & gauxc_integrator_eval_exc_uks, &
    & gauxc_integrator_eval_exc_gks, &
    & gauxc_integrator_eval_exc_vxc_rks, &
    & gauxc_integrator_eval_exc_vxc_uks, &
    & gauxc_integrator_eval_exc_vxc_onedft_uks, &
    & gauxc_integrator_eval_exc_vxc_gks

  public :: &
    & gauxc_delete, &
    & gauxc_get_instance, &
    & gauxc_get_shared_instance, &
    & gauxc_eval_exc, &
    & gauxc_eval_exc_vxc

  !> @brief GauXC XCIntegrator handle
  type, bind(c), public :: gauxc_integrator_type
    type(gauxc_header_type) :: header = gauxc_header_type(gauxc_type_integrator)
    type(c_ptr) :: ptr = c_null_ptr
    logical(c_bool) :: owned = .true._c_bool
  end type gauxc_integrator_type

  !> @brief GauXC XCIntegratorFactory handle
  type, bind(c), public :: gauxc_integrator_factory_type
    type(gauxc_header_type) :: header = gauxc_header_type(gauxc_type_integrator_factory)
    type(c_ptr) :: ptr = c_null_ptr
  end type gauxc_integrator_factory_type

  interface gauxc_delete
    !> @brief Delete an XCIntegrator instance.
    subroutine gauxc_integrator_delete(status, integrator) bind(c)
      import :: gauxc_status_type, gauxc_integrator_type
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param integrator Handle to the XCIntegrator to delete.
      type(gauxc_integrator_type), intent(inout) :: integrator
    end subroutine gauxc_integrator_delete
  end interface gauxc_delete

  interface
    !> @brief Create a new XCIntegratorFactory instance.
    function gauxc_integrator_factory_new_c(status, execution_space, &
      integrator_input_type, integrator_kernel_name, local_work_kernel_name, &
      reduction_kernel_name) result(factory) bind(c, name="gauxc_integrator_factory_new")
      import :: gauxc_status_type, gauxc_integrator_factory_type, c_int, c_char
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param execution_space Execution space to use.
      integer(c_int), value :: execution_space
      !> @param integrator_input_type Type of integrator input.
      character(kind=c_char), intent(in) :: integrator_input_type(*)
      !> @param integrator_kernel_name Name of the integrator kernel.
      character(kind=c_char), intent(in) :: integrator_kernel_name(*)
      !> @param local_work_kernel_name Name of the local work kernel.
      character(kind=c_char), intent(in) :: local_work_kernel_name(*)
      !> @param reduction_kernel_name Name of the reduction kernel.
      character(kind=c_char), intent(in) :: reduction_kernel_name(*)
      !> @return Handle to the created XCIntegratorFactory.
      type(gauxc_integrator_factory_type) :: factory
    end function gauxc_integrator_factory_new_c
  end interface

  interface gauxc_delete
    !> @brief Delete an XCIntegratorFactory instance.
    subroutine gauxc_integrator_factory_delete(status, factory) bind(c)
      import :: gauxc_status_type, gauxc_integrator_factory_type
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param factory Handle to the XCIntegratorFactory to delete.
      type(gauxc_integrator_factory_type), intent(inout) :: factory
    end subroutine gauxc_integrator_factory_delete
  end interface

  interface gauxc_get_instance
    !> @brief Get an XCIntegrator instance from an XCIntegratorFactory.
    function gauxc_integrator_factory_get_instance(status, factory, func, lb) &
      result(integrator) bind(c)
      import :: gauxc_status_type, gauxc_integrator_factory_type, &
        & gauxc_integrator_type, gauxc_functional_type, gauxc_load_balancer_type
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param factory Handle to the XCIntegratorFactory.
      type(gauxc_integrator_factory_type), value :: factory
      !> @param func Handle to the XCFunctional.
      type(gauxc_functional_type), value :: func
      !> @param lb Handle to the LoadBalancer.
      type(gauxc_load_balancer_type), value :: lb
      !> @return Handle to the created XCIntegrator.
      type(gauxc_integrator_type) :: integrator
    end function gauxc_integrator_factory_get_instance
  end interface

  interface gauxc_get_shared_instance
    !> @brief Get a shared XCIntegrator instance from an XCIntegratorFactory.
    function gauxc_integrator_factory_get_shared_instance(status, factory, func, lb) &
      result(integrator) bind(c)
      import :: gauxc_status_type, gauxc_integrator_factory_type, &
        & gauxc_integrator_type, gauxc_functional_type, gauxc_load_balancer_type
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param factory Handle to the XCIntegratorFactory.
      type(gauxc_integrator_factory_type), value :: factory
      !> @param func Handle to the XCFunctional.
      type(gauxc_functional_type), value :: func
      !> @param lb Handle to the LoadBalancer.
      type(gauxc_load_balancer_type), value :: lb
      !> @return Handle to the created XCIntegrator.
      type(gauxc_integrator_type) :: integrator
    end function gauxc_integrator_factory_get_shared_instance
  end interface gauxc_get_shared_instance

  interface gauxc_integrate_den
    !> @brief Integrate the density matrix to get the number of electrons.
    subroutine gauxc_integrator_integrate_den(status, integrator, density_matrix, den) bind(c)
      import :: gauxc_status_type, gauxc_integrator_type, gauxc_matrix_type, c_double
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param density_matrix Density matrix container.
      type(gauxc_matrix_type), value :: density_matrix
      !> @param den Pointer to store the number of electrons.
      real(c_double), intent(out) :: den
    end subroutine gauxc_integrator_integrate_den
  end interface gauxc_integrate_den

  interface gauxc_eval_exc
    !> @brief Evaluate the exchange-correlation energy for RKS.
    subroutine gauxc_integrator_eval_exc_rks(status, integrator, density_matrix, exc) bind(c)
      import :: gauxc_status_type, gauxc_integrator_type, gauxc_matrix_type, c_double
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param density_matrix Density matrix container for RKS.
      type(gauxc_matrix_type), value :: density_matrix
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
    end subroutine gauxc_integrator_eval_exc_rks

    !> @brief Evaluate the exchange-correlation energy for UKS.
    subroutine gauxc_integrator_eval_exc_uks(status, integrator, &
      density_matrix_s, density_matrix_z, exc) bind(c)
      import :: gauxc_status_type, gauxc_integrator_type, gauxc_matrix_type, c_double
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param density_matrix_s Density matrix container for total density.
      type(gauxc_matrix_type), value :: density_matrix_s
      !> @param density_matrix_z Density matrix container for spin density.
      type(gauxc_matrix_type), value :: density_matrix_z
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
    end subroutine gauxc_integrator_eval_exc_uks

    !> @brief Evaluate the exchange-correlation energy for GKS.
    subroutine gauxc_integrator_eval_exc_gks(status, integrator, &
      density_matrix_s, density_matrix_z, density_matrix_x, density_matrix_y, exc) bind(c)
      import :: gauxc_status_type, gauxc_integrator_type, gauxc_matrix_type, c_double
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param density_matrix_s Density matrix container for total density.
      type(gauxc_matrix_type), value :: density_matrix_s
      !> @param density_matrix_z Density matrix container for spin z density.
      type(gauxc_matrix_type), value :: density_matrix_z
      !> @param density_matrix_x Density matrix container for spin x component.
      type(gauxc_matrix_type), value :: density_matrix_x
      !> @param density_matrix_y Density matrix container for spin y component.
      type(gauxc_matrix_type), value :: density_matrix_y
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
    end subroutine gauxc_integrator_eval_exc_gks
  end interface gauxc_eval_exc

  interface gauxc_eval_exc_vxc
    !> @brief Evaluate the exchange-correlation energy and potential for RKS.
    subroutine gauxc_integrator_eval_exc_vxc_rks(status, integrator, &
      density_matrix, exc, vxc_matrix) bind(c)
      import :: gauxc_status_type, gauxc_integrator_type, gauxc_matrix_type, c_double
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param density_matrix Density matrix container for RKS.
      type(gauxc_matrix_type), value :: density_matrix
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
      !> @param vxc_matrix Matrix container to store the exchange-correlation potential.
      type(gauxc_matrix_type), intent(out) :: vxc_matrix
    end subroutine gauxc_integrator_eval_exc_vxc_rks

    !> @brief Evaluate the exchange-correlation energy and potential for UKS.
    subroutine gauxc_integrator_eval_exc_vxc_uks(status, integrator, &
      density_matrix_s, density_matrix_z, exc, vxc_matrix_s, vxc_matrix_z) bind(c)
      import :: gauxc_status_type, gauxc_integrator_type, gauxc_matrix_type, c_double
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param density_matrix_s Density matrix container for total density.
      type(gauxc_matrix_type), value :: density_matrix_s
      !> @param density_matrix_z Density matrix container for spin density.
      type(gauxc_matrix_type), value :: density_matrix_z
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
      !> @param vxc_matrix_s Matrix container for total density potential.
      type(gauxc_matrix_type), intent(out) :: vxc_matrix_s
      !> @param vxc_matrix_z Matrix container for spin density potential.
      type(gauxc_matrix_type), intent(out) :: vxc_matrix_z
    end subroutine gauxc_integrator_eval_exc_vxc_uks

    !> @brief Evaluate the exchange-correlation energy and potential for GKS.
    subroutine gauxc_integrator_eval_exc_vxc_gks(status, integrator, &
      density_matrix_s, density_matrix_z, density_matrix_x, density_matrix_y, &
      exc, vxc_matrix_s, vxc_matrix_z, vxc_matrix_x, vxc_matrix_y) bind(c)
      import :: gauxc_status_type, gauxc_integrator_type, gauxc_matrix_type, c_double
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param density_matrix_s Density matrix container for total density.
      type(gauxc_matrix_type), value :: density_matrix_s
      !> @param density_matrix_z Density matrix container for spin z density.
      type(gauxc_matrix_type), value :: density_matrix_z
      !> @param density_matrix_x Density matrix container for spin x component.
      type(gauxc_matrix_type), value :: density_matrix_x
      !> @param density_matrix_y Density matrix container for spin y component.
      type(gauxc_matrix_type), value :: density_matrix_y
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
      !> @param vxc_matrix_s Matrix container for total density potential.
      type(gauxc_matrix_type), intent(out) :: vxc_matrix_s
      !> @param vxc_matrix_z Matrix container for spin z density potential.
      type(gauxc_matrix_type), intent(out) :: vxc_matrix_z
      !> @param vxc_matrix_x Matrix container for spin x component potential.
      type(gauxc_matrix_type), intent(out) :: vxc_matrix_x
      !> @param vxc_matrix_y Matrix container for spin y component potential.
      type(gauxc_matrix_type), intent(out) :: vxc_matrix_y
    end subroutine gauxc_integrator_eval_exc_vxc_gks
  end interface gauxc_eval_exc_vxc

  interface
    !> @brief Evaluate the exchange-correlation energy and potential for UKS.
    subroutine gauxc_integrator_eval_exc_vxc_onedft_uks_c(status, integrator, &
      density_matrix_s, density_matrix_z, model, exc, vxc_matrix_s, vxc_matrix_z) &
      bind(c, name="gauxc_integrator_eval_exc_vxc_onedft_uks")
      import :: gauxc_status_type, gauxc_integrator_type, gauxc_matrix_type, c_double, &
        & c_char
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(out) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param density_matrix_s Density matrix container for total density.
      type(gauxc_matrix_type), value :: density_matrix_s
      !> @param density_matrix_z Density matrix container for spin density.
      type(gauxc_matrix_type), value :: density_matrix_z
      !> @param String specifying the OneDFT model to use.
      character(kind=c_char), intent(in) :: model(*)
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
      !> @param vxc_matrix_s Matrix container for total density potential.
      type(gauxc_matrix_type), intent(out) :: vxc_matrix_s
      !> @param vxc_matrix_z Matrix container for spin density potential.
      type(gauxc_matrix_type), intent(out) :: vxc_matrix_z
    end subroutine gauxc_integrator_eval_exc_vxc_onedft_uks_c
  end interface

  interface gauxc_eval_exc_vxc
    module procedure gauxc_integrator_eval_exc_vxc_onedft_uks
  end interface gauxc_eval_exc_vxc

contains

  !> @brief Create a new XCIntegratorFactory instance.
  function gauxc_integrator_factory_new(status, execution_space, &
    integrator_input_type, integrator_kernel_name, local_work_kernel_name, &
    reduction_kernel_name) result(factory)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(out) :: status
    !> @param execution_space Execution space to use.
    integer(c_int), value :: execution_space
    !> @param integrator_input_type Type of integrator input.
    character(kind=c_char, len=*), intent(in) :: integrator_input_type
    !> @param integrator_kernel_name Name of the integrator kernel.
    character(kind=c_char, len=*), intent(in) :: integrator_kernel_name
    !> @param local_work_kernel_name Name of the local work kernel.
    character(kind=c_char, len=*), intent(in) :: local_work_kernel_name
    !> @param reduction_kernel_name Name of the reduction kernel.
    character(kind=c_char, len=*), intent(in) :: reduction_kernel_name
    !> @return Handle to the created XCIntegratorFactory.
    type(gauxc_integrator_factory_type) :: factory

    character(kind=c_char), allocatable :: c_integrator_input_type(:)
    character(kind=c_char), allocatable :: c_integrator_kernel_name(:)
    character(kind=c_char), allocatable :: c_local_work_kernel_name(:)
    character(kind=c_char), allocatable :: c_reduction_kernel_name(:)

    c_integrator_input_type = transfer(integrator_input_type//c_null_char, &
      & [character(kind=c_char) ::], len(integrator_input_type)+1)
    c_integrator_kernel_name = transfer(integrator_kernel_name//c_null_char, &
      & [character(kind=c_char) ::], len(integrator_kernel_name)+1)
    c_local_work_kernel_name = transfer(local_work_kernel_name//c_null_char, &
      & [character(kind=c_char) ::], len(local_work_kernel_name)+1)
    c_reduction_kernel_name = transfer(reduction_kernel_name//c_null_char, &
      & [character(kind=c_char) ::], len(reduction_kernel_name)+1)

    factory = gauxc_integrator_factory_new_c(status, execution_space, &
      c_integrator_input_type, c_integrator_kernel_name, &
      c_local_work_kernel_name, c_reduction_kernel_name)
  end function gauxc_integrator_factory_new

  !> @brief Evaluate the exchange-correlation energy and potential for UKS.
  subroutine gauxc_integrator_eval_exc_vxc_onedft_uks(status, integrator, &
    density_matrix_s, density_matrix_z, model, exc, vxc_matrix_s, vxc_matrix_z)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(out) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix_s Density matrix container for total density.
    type(gauxc_matrix_type), value :: density_matrix_s
    !> @param density_matrix_z Density matrix container for spin density.
    type(gauxc_matrix_type), value :: density_matrix_z
    !> @param String specifying the OneDFT model to use.
    character(kind=c_char, len=*), intent(in) :: model
    !> @param exc Pointer to store the exchange-correlation energy.
    real(c_double), intent(inout) :: exc
    !> @param vxc_matrix_s Matrix container for total density potential.
    type(gauxc_matrix_type), intent(inout) :: vxc_matrix_s
    !> @param vxc_matrix_z Matrix container for spin density potential.
    type(gauxc_matrix_type), intent(inout) :: vxc_matrix_z

    character(kind=c_char), allocatable :: c_model(:)

    c_model = transfer(model//c_null_char, [character(kind=c_char) ::], &
      & len(model)+1)

    call gauxc_integrator_eval_exc_vxc_onedft_uks_c(status, integrator, &
      density_matrix_s, density_matrix_z, c_model, exc, vxc_matrix_s, vxc_matrix_z)
  end subroutine gauxc_integrator_eval_exc_vxc_onedft_uks
end module gauxc_integrator