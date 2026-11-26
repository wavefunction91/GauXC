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
  use iso_c_binding, only : c_ptr, c_null_ptr, c_bool, c_int, c_int64_t, c_char, c_double, c_null_char
  use gauxc_status, only : gauxc_status_type
  use gauxc_types, only : gauxc_header_type, gauxc_type_integrator
  use gauxc_xc_functional, only : gauxc_functional_type
  use gauxc_load_balancer, only : gauxc_load_balancer_type
  implicit none
  private

  public :: &
    & gauxc_integrator_new, &
    & gauxc_integrator_delete, &
    & gauxc_integrator_integrate_den, &
    & gauxc_integrator_eval_exc_rks, &
    & gauxc_integrator_eval_exc_uks, &
    & gauxc_integrator_eval_exc_gks, &
    & gauxc_integrator_eval_exc_vxc_rks, &
    & gauxc_integrator_eval_exc_vxc_uks, &
    & gauxc_integrator_eval_exc_vxc_gks, &
    & gauxc_integrator_eval_exc_grad_rks, &
    & gauxc_integrator_eval_exc_grad_uks, &
    & gauxc_integrator_eval_exx_rks, &
    & gauxc_integrator_eval_fxc_contraction_rks, &
    & gauxc_integrator_eval_fxc_contraction_uks

  public :: &
    & gauxc_delete, &
    & gauxc_eval_exc, &
    & gauxc_eval_exc_vxc, &
    & gauxc_eval_exc_grad, &
    & gauxc_eval_exx, &
    & gauxc_eval_fxc_contraction

  !> @brief GauXC XCIntegrator handle
  type, bind(c), public :: gauxc_integrator_type
    type(gauxc_header_type) :: hdr = gauxc_header_type(gauxc_type_integrator)
    type(c_ptr) :: ptr = c_null_ptr
  end type gauxc_integrator_type

  interface gauxc_delete
    !> @brief Delete an XCIntegrator instance.
    subroutine gauxc_integrator_delete(status, integrator) bind(c)
      import :: gauxc_status_type, gauxc_integrator_type
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator to delete.
      type(gauxc_integrator_type), intent(inout) :: integrator
    end subroutine gauxc_integrator_delete
  end interface gauxc_delete

  interface
    !> @brief Create a new XCIntegrator instance.
    function gauxc_integrator_new_c(status, func, lb, execution_space, &
      integrator_input_type, integrator_kernel_name, local_work_kernel_name, &
      reduction_kernel_name) result(integrator) bind(c, name="gauxc_integrator_new")
      import :: gauxc_status_type, gauxc_integrator_type, gauxc_functional_type, &
        & gauxc_load_balancer_type, c_int, c_char
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param func Handle to the functional.
      type(gauxc_functional_type), value :: func
      !> @param lb Handle to the load balancer.
      type(gauxc_load_balancer_type), value :: lb
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
      !> @return Handle to the created XCIntegrator.
      type(gauxc_integrator_type) :: integrator
    end function gauxc_integrator_new_c
  end interface

  interface gauxc_integrate_den
    !> @brief Integrate the density matrix to get the number of electrons.
    subroutine gauxc_integrator_integrate_den_c(status, integrator, m, n, density_matrix, ldp, den) &
        & bind(c, name="gauxc_integrator_integrate_den")
      import :: gauxc_status_type, gauxc_integrator_type, c_double, c_int64_t
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param m Number of rows in the density matrix.
      integer(c_int64_t), value :: m
      !> @param n Number of columns in the density matrix.
      integer(c_int64_t), value :: n
      !> @param density_matrix Pointer to the density matrix data.
      real(c_double), intent(in) :: density_matrix(ldp, *)
      !> @param ldp Leading dimension of the density matrix.
      integer(c_int64_t), value :: ldp
      !> @param den Pointer to store the number of electrons.
      real(c_double), intent(out) :: den
    end subroutine gauxc_integrator_integrate_den_c
    module procedure gauxc_integrator_integrate_den
  end interface gauxc_integrate_den

  interface gauxc_eval_exc
    !> @brief Evaluate the exchange-correlation energy for RKS.
    subroutine gauxc_integrator_eval_exc_rks_c(status, integrator, m, n, density_matrix, ldp, exc) &
        & bind(c, name="gauxc_integrator_eval_exc_rks")
      import :: gauxc_status_type, gauxc_integrator_type, c_double, c_int64_t
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param m Number of rows in the density matrix.
      integer(c_int64_t), value :: m
      !> @param n Number of columns in the density matrix.
      integer(c_int64_t), value :: n
      !> @param density_matrix Pointer to the density matrix data.
      real(c_double), intent(in) :: density_matrix(ldp, *)
      !> @param ldp Leading dimension of the density matrix.
      integer(c_int64_t), value :: ldp
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
    end subroutine gauxc_integrator_eval_exc_rks_c
    module procedure gauxc_integrator_eval_exc_rks

    !> @brief Evaluate the exchange-correlation energy for UKS.
    subroutine gauxc_integrator_eval_exc_uks_c(status, integrator, m, n, &
        & density_matrix_s, ldp_s, density_matrix_z, ldp_z, exc) &
        & bind(c, name="gauxc_integrator_eval_exc_uks")
      import :: gauxc_status_type, gauxc_integrator_type, c_double, c_int64_t
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param m Number of rows in the density matrices.
      integer(c_int64_t), value :: m
      !> @param n Number of columns in the density matrices.
      integer(c_int64_t), value :: n
      !> @param density_matrix_s Pointer to the total density matrix data.
      real(c_double), intent(in) :: density_matrix_s(ldp_s, *)
      !> @param ldp_s Leading dimension of the total density matrix.
      integer(c_int64_t), value :: ldp_s
      !> @param density_matrix_z Pointer to the spin density matrix data.
      real(c_double), intent(in) :: density_matrix_z(ldp_z, *)
      !> @param ldp_z Leading dimension of the spin density matrix.
      integer(c_int64_t), value :: ldp_z
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
    end subroutine gauxc_integrator_eval_exc_uks_c
    module procedure gauxc_integrator_eval_exc_uks

    !> @brief Evaluate the exchange-correlation energy for GKS.
    subroutine gauxc_integrator_eval_exc_gks_c(status, integrator, m, n, &
        & density_matrix_s, ldp_s, density_matrix_z, ldp_z, &
        & density_matrix_y, ldp_y, density_matrix_x, ldp_x, exc) &
        & bind(c, name="gauxc_integrator_eval_exc_gks")
      import :: gauxc_status_type, gauxc_integrator_type, c_double, c_int64_t
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param m Number of rows in the density matrices.
      integer(c_int64_t), value :: m
      !> @param n Number of columns in the density matrices.
      integer(c_int64_t), value :: n
      !> @param density_matrix_s Pointer to the total density matrix data.
      real(c_double), intent(in) :: density_matrix_s(ldp_s, *)
      !> @param ldp_s Leading dimension of the total density matrix.
      integer(c_int64_t), value :: ldp_s
      !> @param density_matrix_z Pointer to the spin z density matrix data.
      real(c_double), intent(in) :: density_matrix_z(ldp_z, *)
      !> @param ldp_z Leading dimension of the spin z density matrix.
      integer(c_int64_t), value :: ldp_z
      !> @param density_matrix_y Pointer to the spin y density matrix data.
      real(c_double), intent(in) :: density_matrix_y(ldp_y, *)
      !> @param ldp_y Leading dimension of the spin y density matrix.
      integer(c_int64_t), value :: ldp_y
      !> @param density_matrix_x Pointer to the spin x density matrix data.
      real(c_double), intent(in) :: density_matrix_x(ldp_x, *)
      !> @param ldp_x Leading dimension of the spin x density matrix.
      integer(c_int64_t), value :: ldp_x
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
    end subroutine gauxc_integrator_eval_exc_gks_c
    module procedure gauxc_integrator_eval_exc_gks
  end interface gauxc_eval_exc

  interface gauxc_eval_exc_vxc
    !> @brief Evaluate the exchange-correlation energy and potential for RKS.
    subroutine gauxc_integrator_eval_exc_vxc_rks_c(status, integrator, m, n, &
        & density_matrix, ldp, exc, vxc_matrix, ldvxc) &
        & bind(c, name="gauxc_integrator_eval_exc_vxc_rks")
      import :: gauxc_status_type, gauxc_integrator_type, c_double, c_int64_t
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param m Number of rows in the density matrix.
      integer(c_int64_t), value :: m
      !> @param n Number of columns in the density matrix.
      integer(c_int64_t), value :: n
      !> @param density_matrix Pointer to the density matrix data.
      real(c_double), intent(in) :: density_matrix(ldp, *)
      !> @param ldp Leading dimension of the density matrix.
      integer(c_int64_t), value :: ldp
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
      !> @param vxc_matrix Pointer to the potential matrix data.
      real(c_double), intent(out) :: vxc_matrix(ldvxc, *)
      !> @param ldvxc Leading dimension of the potential matrix.
      integer(c_int64_t), value :: ldvxc
    end subroutine gauxc_integrator_eval_exc_vxc_rks_c
    module procedure gauxc_integrator_eval_exc_vxc_rks

    !> @brief Evaluate the exchange-correlation energy and potential for UKS.
    subroutine gauxc_integrator_eval_exc_vxc_uks_c(status, integrator, m, n, &
        & density_matrix_s, ldp_s, density_matrix_z, ldp_z, &
        & exc, vxc_matrix_s, ldvxc_s, vxc_matrix_z, ldvxc_z) &
        & bind(c, name="gauxc_integrator_eval_exc_vxc_uks")
      import :: gauxc_status_type, gauxc_integrator_type, c_double, c_int64_t
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param m Number of rows in the density matrices.
      integer(c_int64_t), value :: m
      !> @param n Number of columns in the density matrices.
      integer(c_int64_t), value :: n
      !> @param density_matrix_s Pointer to the total density matrix data.
      real(c_double), intent(in) :: density_matrix_s(ldp_s, *)
      !> @param ldp_s Leading dimension of the total density matrix.
      integer(c_int64_t), value :: ldp_s
      !> @param density_matrix_z Pointer to the spin density matrix data.
      real(c_double), intent(in) :: density_matrix_z(ldp_z, *)
      !> @param ldp_z Leading dimension of the spin density matrix.
      integer(c_int64_t), value :: ldp_z
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
      !> @param vxc_matrix_s Pointer to the total density potential matrix data.
      real(c_double), intent(out) :: vxc_matrix_s(ldvxc_s, *)
      !> @param ldvxc_s Leading dimension of the total density potential matrix.
      integer(c_int64_t), value :: ldvxc_s
      !> @param vxc_matrix_z Pointer to the spin density potential matrix data.
      real(c_double), intent(out) :: vxc_matrix_z(ldvxc_z, *)
      !> @param ldvxc_z Leading dimension of the spin density potential matrix.
      integer(c_int64_t), value :: ldvxc_z
    end subroutine gauxc_integrator_eval_exc_vxc_uks_c
    module procedure gauxc_integrator_eval_exc_vxc_uks

    !> @brief Evaluate the exchange-correlation energy and potential for GKS.
    subroutine gauxc_integrator_eval_exc_vxc_gks_c(status, integrator, m, n, &
        & density_matrix_s, ldp_s, density_matrix_z, ldp_z, density_matrix_y, ldp_y, density_matrix_x, ldp_x, &
        & exc, vxc_matrix_s, ldvxc_s, vxc_matrix_z, ldvxc_z, vxc_matrix_y, ldvxc_y, vxc_matrix_x, ldvxc_x) &
        & bind(c, name="gauxc_integrator_eval_exc_vxc_gks")
      import :: gauxc_status_type, gauxc_integrator_type, c_double, c_int64_t
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param m Number of rows in the density matrices.
      integer(c_int64_t), value :: m
      !> @param n Number of columns in the density matrices.
      integer(c_int64_t), value :: n
      !> @param density_matrix_s Pointer to the total density matrix data.
      real(c_double), intent(in) :: density_matrix_s(ldp_s, *)
      !> @param ldp_s Leading dimension of the total density matrix.
      integer(c_int64_t), value :: ldp_s
      !> @param density_matrix_z Pointer to the spin z density matrix data.
      real(c_double), intent(in) :: density_matrix_z(ldp_z, *)
      !> @param ldp_z Leading dimension of the spin z density matrix.
      integer(c_int64_t), value :: ldp_z
      !> @param density_matrix_y Pointer to the spin y density matrix data.
      real(c_double), intent(in) :: density_matrix_y(ldp_y, *)
      !> @param ldp_y Leading dimension of the spin y density matrix.
      integer(c_int64_t), value :: ldp_y
      !> @param density_matrix_x Pointer to the spin x density matrix data.
      real(c_double), intent(in) :: density_matrix_x(ldp_x, *)
      !> @param ldp_x Leading dimension of the spin x density matrix.
      integer(c_int64_t), value :: ldp_x
      !> @param exc Pointer to store the exchange-correlation energy.
      real(c_double), intent(out) :: exc
      !> @param vxc_matrix_s Pointer to the total density potential matrix data.
      real(c_double), intent(out) :: vxc_matrix_s(ldvxc_s, *)
      !> @param ldvxc_s Leading dimension of the total density potential matrix.
      integer(c_int64_t), value :: ldvxc_s
      !> @param vxc_matrix_z Pointer to the spin density potential matrix data.
      real(c_double), intent(out) :: vxc_matrix_z(ldvxc_z, *)
      !> @param ldvxc_z Leading dimension of the spin density potential matrix.
      integer(c_int64_t), value :: ldvxc_z
      !> @param vxc_matrix_y Pointer to the spin y density potential matrix data.
      real(c_double), intent(out) :: vxc_matrix_y(ldvxc_y, *)
      !> @param ldvxc_y Leading dimension of the spin y density potential matrix.
      integer(c_int64_t), value :: ldvxc_y
      !> @param vxc_matrix_x Pointer to the spin x density potential matrix data.
      real(c_double), intent(out) :: vxc_matrix_x(ldvxc_x, *)
      !> @param ldvxc_x Leading dimension of the spin x density potential matrix.
      integer(c_int64_t), value :: ldvxc_x
    end subroutine gauxc_integrator_eval_exc_vxc_gks_c
    module procedure gauxc_integrator_eval_exc_vxc_gks
  end interface gauxc_eval_exc_vxc

  interface gauxc_eval_exc_grad
    subroutine gauxc_integrator_eval_exc_grad_rks_c(status, integrator, m, n, &
        & density_matrix, ldp, exc_grad) &
        & bind(c, name="gauxc_integrator_eval_exc_grad_rks")
      import :: gauxc_status_type, gauxc_integrator_type, c_double, c_int64_t
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param m Number of rows in the density matrix.
      integer(c_int64_t), value :: m
      !> @param n Number of columns in the density matrix.
      integer(c_int64_t), value :: n
      !> @param density_matrix Pointer to the density matrix data.
      real(c_double), intent(in) :: density_matrix(ldp, *)
      !> @param ldp Leading dimension of the density matrix.
      integer(c_int64_t), value :: ldp
      !> @param exc_grad Pointer to the gradient.
      real(c_double), intent(out) :: exc_grad(*)
    end subroutine gauxc_integrator_eval_exc_grad_rks_c
    module procedure gauxc_integrator_eval_exc_grad_rks

    subroutine gauxc_integrator_eval_exc_grad_uks_c(status, integrator, m, n, &
        & density_matrix_s, ldp_s, density_matrix_z, ldp_z, exc_grad) &
        & bind(c, name="gauxc_integrator_eval_exc_grad_uks")
      import :: gauxc_status_type, gauxc_integrator_type, c_double, c_int64_t
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param m Number of rows in the density matrices.
      integer(c_int64_t), value :: m
      !> @param n Number of columns in the density matrices.
      integer(c_int64_t), value :: n
      !> @param density_matrix_s Pointer to the total density matrix data.
      real(c_double), intent(in) :: density_matrix_s(ldp_s, *)
      !> @param ldp_s Leading dimension of the total density matrix.
      integer(c_int64_t), value :: ldp_s
      !> @param density_matrix_z Pointer to the spin density matrix data.
      real(c_double), intent(in) :: density_matrix_z(ldp_z, *)
      !> @param ldp_z Leading dimension of the spin density matrix.
      integer(c_int64_t), value :: ldp_z
      !> @param exc_grad Pointer to the gradient.
      real(c_double), intent(out) :: exc_grad(*)
    end subroutine gauxc_integrator_eval_exc_grad_uks_c
    module procedure gauxc_integrator_eval_exc_grad_uks
  end interface gauxc_eval_exc_grad

  interface gauxc_eval_exx
    subroutine gauxc_integrator_eval_exx_rks_c(status, integrator, m, n, density_matrix, ldp, k_matrix, ldk) &
        & bind(c, name="gauxc_integrator_eval_exx_rks")
      import :: gauxc_status_type, gauxc_integrator_type, c_double, c_int64_t
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param m Number of rows in the density matrix.
      integer(c_int64_t), value :: m
      !> @param n Number of columns in the density matrix.
      integer(c_int64_t), value :: n
      !> @param density_matrix Pointer to the density matrix data.
      real(c_double), intent(in) :: density_matrix(ldp, *)
      !> @param ldp Leading dimension of the density matrix.
      integer(c_int64_t), value :: ldp
      !> @param k_matrix Pointer to the exchange matrix data.
      real(c_double), intent(out) :: k_matrix(ldk, *)
      !> @param ldk Leading dimension of the exchange matrix.
      integer(c_int64_t), value :: ldk
    end subroutine gauxc_integrator_eval_exx_rks_c
    module procedure gauxc_integrator_eval_exx_rks
  end interface gauxc_eval_exx

  interface gauxc_eval_fxc_contraction
    subroutine gauxc_integrator_eval_fxc_contraction_rks_c(status, integrator, m, n, &
        & density_matrix, ldp, t_density_matrix, ldtp, fxc_contraction_matrix, ldfxc) &
        & bind(c, name="gauxc_integrator_eval_fxc_contraction_rks")
      import :: gauxc_status_type, gauxc_integrator_type, c_double, c_int64_t
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param m Number of rows in the density matrix.
      integer(c_int64_t), value :: m
      !> @param n Number of columns in the density matrix.
      integer(c_int64_t), value :: n
      !> @param density_matrix Pointer to the density matrix data.
      real(c_double), intent(in) :: density_matrix(ldp, *)
      !> @param ldp Leading dimension of the density matrix.
      integer(c_int64_t), value :: ldp
      !> @param t_density_matrix Pointer to the density matrix data.
      real(c_double), intent(in) :: t_density_matrix(ldtp, *)
      !> @param ldtp Leading dimension of the t_density_matrix.
      integer(c_int64_t), value :: ldtp
      !> @param fxc_contraction_matrix Pointer to the fxc contraction matrix data.
      real(c_double), intent(out) :: fxc_contraction_matrix(ldfxc, *)
      !> @param ldfxc Leading dimension of the fxc contraction matrix.
      integer(c_int64_t), value :: ldfxc
    end subroutine gauxc_integrator_eval_fxc_contraction_rks_c
    module procedure gauxc_integrator_eval_fxc_contraction_rks

    subroutine gauxc_integrator_eval_fxc_contraction_uks_c(status, integrator, m, n, &
        & density_matrix_s, ldp_s, density_matrix_z, ldp_z, &
        & t_density_matrix_s, ldtp_s, t_density_matrix_z, ldtp_z, &
        & fxc_contraction_matrix_s, ldfxc_s, fxc_contraction_matrix_z, ldfxc_z) &
        & bind(c, name="gauxc_integrator_eval_fxc_contraction_uks")
      import :: gauxc_status_type, gauxc_integrator_type, c_double, c_int64_t
      implicit none
      !> @param status Status object to capture any errors.
      type(gauxc_status_type), intent(inout) :: status
      !> @param integrator Handle to the XCIntegrator.
      type(gauxc_integrator_type), value :: integrator
      !> @param m Number of rows in the density matrices.
      integer(c_int64_t), value :: m
      !> @param n Number of columns in the density matrices.
      integer(c_int64_t), value :: n
      !> @param density_matrix_s Pointer to the total density matrix data.
      real(c_double), intent(in) :: density_matrix_s(ldp_s, *)
      !> @param ldp_s Leading dimension of the total density matrix.
      integer(c_int64_t), value :: ldp_s
      !> @param density_matrix_z Pointer to the spin density matrix data.
      real(c_double), intent(in) :: density_matrix_z(ldp_z, *)
      !> @param ldp_z Leading dimension of the spin density matrix.
      integer(c_int64_t), value :: ldp_z
      !> @param t_density_matrix_s Pointer to the total density matrix data.
      real(c_double), intent(in) :: t_density_matrix_s(ldtp_s, *)
      !> @param ldtp_s Leading dimension of the total density matrix.
      integer(c_int64_t), value :: ldtp_s
      !> @param t_density_matrix_z Pointer to the spin density matrix data.
      real(c_double), intent(in) :: t_density_matrix_z(ldtp_z, *)
      !> @param ldtp_z Leading dimension of the spin density matrix.
      integer(c_int64_t), value :: ldtp_z
      !> @param fxc_contraction_matrix_s Pointer to the total density fxc contraction matrix data.
      real(c_double), intent(out) :: fxc_contraction_matrix_s(ldfxc_s, *)
      !> @param ldfxc_s Leading dimension of the total density fxc contraction matrix.
      integer(c_int64_t), value :: ldfxc_s
      !> @param fxc_contraction_matrix_z Pointer to the spin density fxc contraction matrix data.
      real(c_double), intent(out) :: fxc_contraction_matrix_z(ldfxc_z, *)
      !> @param ldfxc_z Leading dimension of the spin density fxc contraction matrix.
      integer(c_int64_t), value :: ldfxc_z
    end subroutine gauxc_integrator_eval_fxc_contraction_uks_c
    module procedure gauxc_integrator_eval_fxc_contraction_uks
  end interface gauxc_eval_fxc_contraction

contains

  !> @brief Create a new XCIntegrator instance.
  function gauxc_integrator_new(status, func, lb, execution_space, &
    integrator_input_type, integrator_kernel_name, local_work_kernel_name, &
    reduction_kernel_name) result(integrator)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param func Handle to the functional.
    type(gauxc_functional_type), value :: func
    !> @param lb Handle to the load balancer.
    type(gauxc_load_balancer_type), value :: lb
    !> @param execution_space Execution space to use.
    integer(c_int), value :: execution_space
    !> @param integrator_input_type Type of integrator input.
    character(kind=c_char, len=*), intent(in), optional :: integrator_input_type
    !> @param integrator_kernel_name Name of the integrator kernel.
    character(kind=c_char, len=*), intent(in), optional :: integrator_kernel_name
    !> @param local_work_kernel_name Name of the local work kernel.
    character(kind=c_char, len=*), intent(in), optional :: local_work_kernel_name
    !> @param reduction_kernel_name Name of the reduction kernel.
    character(kind=c_char, len=*), intent(in), optional :: reduction_kernel_name
    !> @return Handle to the created XCIntegrator.
    type(gauxc_integrator_type) :: integrator

    character(kind=c_char, len=*), parameter :: default_input_type = "Replicated"
    character(kind=c_char, len=*), parameter :: default_kernel = "Default"
    character(kind=c_char), allocatable :: c_integrator_input_type(:)
    character(kind=c_char), allocatable :: c_integrator_kernel_name(:)
    character(kind=c_char), allocatable :: c_local_work_kernel_name(:)
    character(kind=c_char), allocatable :: c_reduction_kernel_name(:)

    if (present(integrator_input_type)) then
      c_integrator_input_type = transfer(integrator_input_type//c_null_char, &
        & [character(kind=c_char) ::], len(integrator_input_type)+1)
    else
      c_integrator_input_type = transfer(default_input_type//c_null_char, &
        & [character(kind=c_char) ::], len(default_input_type)+1)
    end if

    if (present(integrator_kernel_name)) then
      c_integrator_kernel_name = transfer(integrator_kernel_name//c_null_char, &
        & [character(kind=c_char) ::], len(integrator_kernel_name)+1)
    else
      c_integrator_kernel_name = transfer(default_kernel//c_null_char, &
        & [character(kind=c_char) ::], len(default_kernel)+1)
    end if

    if (present(local_work_kernel_name)) then
      c_local_work_kernel_name = transfer(local_work_kernel_name//c_null_char, &
        & [character(kind=c_char) ::], len(local_work_kernel_name)+1)
    else
      c_local_work_kernel_name = transfer(default_kernel//c_null_char, &
        & [character(kind=c_char) ::], len(default_kernel)+1)
    end if

    if (present(reduction_kernel_name)) then
      c_reduction_kernel_name = transfer(reduction_kernel_name//c_null_char, &
        & [character(kind=c_char) ::], len(reduction_kernel_name)+1)
    else
      c_reduction_kernel_name = transfer(default_kernel//c_null_char, &
        & [character(kind=c_char) ::], len(default_kernel)+1)
    end if

    integrator = gauxc_integrator_new_c(status, func, lb, execution_space, &
      c_integrator_input_type, c_integrator_kernel_name, &
      c_local_work_kernel_name, c_reduction_kernel_name)
  end function gauxc_integrator_new

  !> @brief Integrate the density matrix to get the number of electrons.
  subroutine gauxc_integrator_integrate_den(status, integrator, density_matrix, den)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix Pointer to the density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix(:, :)
    !> @param den Pointer to store the number of electrons.
    real(c_double), intent(out) :: den

    integer(c_int64_t) :: m, n, ldp

    m = size(density_matrix, 1, kind=c_int64_t)
    n = size(density_matrix, 2, kind=c_int64_t)
    ldp = size(density_matrix, 1, kind=c_int64_t)
    call gauxc_integrator_integrate_den_c(status, integrator, m, n, density_matrix, ldp, den)
  end subroutine gauxc_integrator_integrate_den

  !> @brief Evaluate the exchange-correlation energy for RKS.
  subroutine gauxc_integrator_eval_exc_rks(status, integrator, density_matrix, exc)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix Pointer to the density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix(:, :)
    !> @param exc Pointer to store the exchange-correlation energy.
    real(c_double), intent(out) :: exc

    integer(c_int64_t) :: m, n, ldp

    m = size(density_matrix, 1, kind=c_int64_t)
    n = size(density_matrix, 2, kind=c_int64_t)
    ldp = size(density_matrix, 1, kind=c_int64_t)
    call gauxc_integrator_eval_exc_rks_c(status, integrator, m, n, density_matrix, ldp, exc)
  end subroutine gauxc_integrator_eval_exc_rks

  !> @brief Evaluate the exchange-correlation energy for UKS.
  subroutine gauxc_integrator_eval_exc_uks(status, integrator, &
      & density_matrix_s, density_matrix_z, exc)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix_s Pointer to the total density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_s(:, :)
    !> @param density_matrix_z Pointer to the spin density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_z(:, :)
    !> @param exc Pointer to store the exchange-correlation energy.
    real(c_double), intent(out) :: exc

    integer(c_int64_t) :: m, n, ldp_s, ldp_z

    m = size(density_matrix_s, 1, kind=c_int64_t)
    n = size(density_matrix_s, 2, kind=c_int64_t)
    ldp_s = size(density_matrix_s, 1, kind=c_int64_t)
    ldp_z = size(density_matrix_z, 1, kind=c_int64_t)
    call gauxc_integrator_eval_exc_uks_c(status, integrator, m, n, &
      & density_matrix_s, ldp_s, density_matrix_z, ldp_z, exc)
  end subroutine gauxc_integrator_eval_exc_uks

  !> @brief Evaluate the exchange-correlation energy for GKS.
  subroutine gauxc_integrator_eval_exc_gks(status, integrator, &
      & density_matrix_s, density_matrix_z, density_matrix_y, density_matrix_x, exc)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix_s Pointer to the total density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_s(:, :)
    !> @param density_matrix_z Pointer to the spin z density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_z(:, :)
    !> @param density_matrix_y Pointer to the spin y density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_y(:, :)
    !> @param density_matrix_x Pointer to the spin x density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_x(:, :)
    !> @param exc Pointer to store the exchange-correlation energy.
    real(c_double), intent(out) :: exc

    integer(c_int64_t) :: m, n, ldp_s, ldp_z, ldp_y, ldp_x

    m = size(density_matrix_s, 1, kind=c_int64_t)
    n = size(density_matrix_s, 2, kind=c_int64_t)
    ldp_s = size(density_matrix_s, 1, kind=c_int64_t)
    ldp_z = size(density_matrix_z, 1, kind=c_int64_t)
    ldp_y = size(density_matrix_y, 1, kind=c_int64_t)
    ldp_x = size(density_matrix_x, 1, kind=c_int64_t)
    call gauxc_integrator_eval_exc_gks_c(status, integrator, m, n, &
      & density_matrix_s, ldp_s, density_matrix_z, ldp_z, &
      & density_matrix_y, ldp_y, density_matrix_x, ldp_x, exc)
  end subroutine gauxc_integrator_eval_exc_gks

  !> @brief Evaluate the exchange-correlation energy and potential for RKS.
  subroutine gauxc_integrator_eval_exc_vxc_rks(status, integrator, &
      & density_matrix, exc, vxc_matrix)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix Pointer to the density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix(:, :)
    !> @param exc Pointer to store the exchange-correlation energy.
    real(c_double), intent(out) :: exc
    !> @param vxc_matrix Pointer to the potential matrix data.
    real(c_double), contiguous, intent(out) :: vxc_matrix(:, :)

    integer(c_int64_t) :: m, n, ldp, ldvxc

    m = size(density_matrix, 1, kind=c_int64_t)
    n = size(density_matrix, 2, kind=c_int64_t)
    ldp = size(density_matrix, 1, kind=c_int64_t)
    ldvxc = size(vxc_matrix, 1, kind=c_int64_t)

    call gauxc_integrator_eval_exc_vxc_rks_c(status, integrator, m, n, &
      & density_matrix, ldp, exc, vxc_matrix, ldvxc)
  end subroutine gauxc_integrator_eval_exc_vxc_rks

  !> @brief Evaluate the exchange-correlation energy and potential for UKS.
  subroutine gauxc_integrator_eval_exc_vxc_uks(status, integrator, &
      & density_matrix_s, density_matrix_z, &
      & exc, vxc_matrix_s, vxc_matrix_z)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix_s Pointer to the total density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_s(:, :)
    !> @param density_matrix_z Pointer to the spin density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_z(:, :)
    !> @param exc Pointer to store the exchange-correlation energy.
    real(c_double), intent(out) :: exc
    !> @param vxc_matrix_s Pointer to the total density potential matrix data.
    real(c_double), contiguous, intent(out) :: vxc_matrix_s(:, :)
    !> @param vxc_matrix_z Pointer to the spin density potential matrix data.
    real(c_double), contiguous, intent(out) :: vxc_matrix_z(:, :)

    integer(c_int64_t) :: m, n, ldp_s, ldp_z, ldvxc_s, ldvxc_z

    m = size(density_matrix_s, 1, kind=c_int64_t)
    n = size(density_matrix_s, 2, kind=c_int64_t)
    ldp_s = size(density_matrix_s, 1, kind=c_int64_t)
    ldp_z = size(density_matrix_z, 1, kind=c_int64_t)
    ldvxc_s = size(vxc_matrix_s, 1, kind=c_int64_t)
    ldvxc_z = size(vxc_matrix_z, 1, kind=c_int64_t)
    call gauxc_integrator_eval_exc_vxc_uks_c(status, integrator, m, n, &
      & density_matrix_s, ldp_s, density_matrix_z, ldp_z, &
      & exc, vxc_matrix_s, ldvxc_s, vxc_matrix_z, ldvxc_z)
  end subroutine gauxc_integrator_eval_exc_vxc_uks

  !> @brief Evaluate the exchange-correlation energy and potential for GKS.
  subroutine gauxc_integrator_eval_exc_vxc_gks(status, integrator, &
      & density_matrix_s, density_matrix_z, density_matrix_y, density_matrix_x, &
      & exc, vxc_matrix_s, vxc_matrix_z, vxc_matrix_y, vxc_matrix_x)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix_s Pointer to the total density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_s(:, :)
    !> @param density_matrix_z Pointer to the spin z density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_z(:, :)
    !> @param density_matrix_y Pointer to the spin y density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_y(:, :)
    !> @param density_matrix_x Pointer to the spin x density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_x(:, :)
    !> @param exc Pointer to store the exchange-correlation energy.
    real(c_double), intent(out) :: exc
    !> @param vxc_matrix_s Pointer to the total density potential matrix data.
    real(c_double), contiguous, intent(out) :: vxc_matrix_s(:, :)
    !> @param vxc_matrix_z Pointer to the spin density potential matrix data.
    real(c_double), contiguous, intent(out) :: vxc_matrix_z(:, :)
    !> @param vxc_matrix_y Pointer to the spin y density potential matrix data.
    real(c_double), contiguous, intent(out) :: vxc_matrix_y(:, :)
    !> @param vxc_matrix_x Pointer to the spin x density potential matrix data.
    real(c_double), contiguous, intent(out) :: vxc_matrix_x(:, :)

    integer(c_int64_t) :: m, n, ldp_s, ldp_z, ldp_y, ldp_x, &
      & ldvxc_s, ldvxc_z, ldvxc_y, ldvxc_x

    m = size(density_matrix_s, 1, kind=c_int64_t)
    n = size(density_matrix_s, 2, kind=c_int64_t)
    ldp_s = size(density_matrix_s, 1, kind=c_int64_t)
    ldp_z = size(density_matrix_z, 1, kind=c_int64_t)
    ldp_y = size(density_matrix_y, 1, kind=c_int64_t)
    ldp_x = size(density_matrix_x, 1, kind=c_int64_t)
    ldvxc_s = size(vxc_matrix_s, 1, kind=c_int64_t)
    ldvxc_z = size(vxc_matrix_z, 1, kind=c_int64_t)
    ldvxc_y = size(vxc_matrix_y, 1, kind=c_int64_t)
    ldvxc_x = size(vxc_matrix_x, 1, kind=c_int64_t)
    call gauxc_integrator_eval_exc_vxc_gks_c(status, integrator, m, n, &
      & density_matrix_s, ldp_s, density_matrix_z, ldp_z, density_matrix_y, ldp_y, density_matrix_x, ldp_x, &
      & exc, vxc_matrix_s, ldvxc_s, vxc_matrix_z, ldvxc_z, vxc_matrix_y, ldvxc_y, vxc_matrix_x, ldvxc_x)
  end subroutine gauxc_integrator_eval_exc_vxc_gks

  subroutine gauxc_integrator_eval_exc_grad_rks(status, integrator, &
      & density_matrix, exc_grad)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix Pointer to the density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix(:, :)
    !> @param exc_grad Pointer to the gradient.
    real(c_double), contiguous, intent(out) :: exc_grad(:)

    integer(c_int64_t) :: m, n, ldp

    m = size(density_matrix, 1, kind=c_int64_t)
    n = size(density_matrix, 2, kind=c_int64_t)
    ldp = size(density_matrix, 1, kind=c_int64_t)

    call gauxc_integrator_eval_exc_grad_rks_c(status, integrator, m, n, &
      & density_matrix, ldp, exc_grad)
  end subroutine gauxc_integrator_eval_exc_grad_rks

  subroutine gauxc_integrator_eval_exc_grad_uks(status, integrator, &
      & density_matrix_s, density_matrix_z, exc_grad)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix_s Pointer to the total density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_s(:, :)
    !> @param density_matrix_z Pointer to the spin density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_z(:, :)
    !> @param exc_grad Pointer to the gradient.
    real(c_double), contiguous, intent(out) :: exc_grad(:)

    integer(c_int64_t) :: m, n, ldp_s, ldp_z

    m = size(density_matrix_s, 1, kind=c_int64_t)
    n = size(density_matrix_s, 2, kind=c_int64_t)
    ldp_s = size(density_matrix_s, 1, kind=c_int64_t)
    ldp_z = size(density_matrix_z, 1, kind=c_int64_t)
    call gauxc_integrator_eval_exc_grad_uks_c(status, integrator, m, n, &
      & density_matrix_s, ldp_s, density_matrix_z, ldp_z, exc_grad)
  end subroutine gauxc_integrator_eval_exc_grad_uks

  subroutine gauxc_integrator_eval_exx_rks(status, integrator, &
      & density_matrix, k_matrix)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix Pointer to the density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix(:, :)
    !> @param k_matrix Pointer to the exchange matrix data.
    real(c_double), contiguous, intent(out) :: k_matrix(:, :)

    integer(c_int64_t) :: m, n, ldp, ldk

    m = size(density_matrix, 1, kind=c_int64_t)
    n = size(density_matrix, 2, kind=c_int64_t)
    ldp = size(density_matrix, 1, kind=c_int64_t)
    ldk = size(k_matrix, 1, kind=c_int64_t)

    call gauxc_integrator_eval_exx_rks_c(status, integrator, m, n, &
      & density_matrix, ldp, k_matrix, ldk)
  end subroutine gauxc_integrator_eval_exx_rks

  subroutine gauxc_integrator_eval_fxc_contraction_rks(status, integrator, &
      & density_matrix, t_density_matrix, fxc_contraction_matrix)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix Pointer to the density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix(:, :)
    !> @param t_density_matrix Pointer to the density matrix data.
    real(c_double), contiguous, intent(in) :: t_density_matrix(:, :)
    !> @param fxc_contraction_matrix Pointer to the fxc contraction matrix data.
    real(c_double), contiguous, intent(out) :: fxc_contraction_matrix(:, :)

    integer(c_int64_t) :: m, n, ldp, ldtp, ldfxc

    m = size(density_matrix, 1, kind=c_int64_t)
    n = size(density_matrix, 2, kind=c_int64_t)
    ldp = size(density_matrix, 1, kind=c_int64_t)
    ldtp = size(t_density_matrix, 1, kind=c_int64_t)
    ldfxc = size(fxc_contraction_matrix, 1, kind=c_int64_t)

    call gauxc_integrator_eval_fxc_contraction_rks_c(status, integrator, m, n, &
      & density_matrix, ldp, t_density_matrix, ldtp, fxc_contraction_matrix, ldfxc)
  end subroutine gauxc_integrator_eval_fxc_contraction_rks

  subroutine gauxc_integrator_eval_fxc_contraction_uks(status, integrator, &
      & density_matrix_s, density_matrix_z, &
      & t_density_matrix_s, t_density_matrix_z, &
      & fxc_contraction_matrix_s, fxc_contraction_matrix_z)
    !> @param status Status object to capture any errors.
    type(gauxc_status_type), intent(inout) :: status
    !> @param integrator Handle to the XCIntegrator.
    type(gauxc_integrator_type), value :: integrator
    !> @param density_matrix_s Pointer to the total density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_s(:, :)
    !> @param density_matrix_z Pointer to the spin density matrix data.
    real(c_double), contiguous, intent(in) :: density_matrix_z(:, :)
    !> @param t_density_matrix_s Pointer to the total density matrix data.
    real(c_double), contiguous, intent(in) :: t_density_matrix_s(:, :)
    !> @param t_density_matrix_z Pointer to the spin density matrix data.
    real(c_double), contiguous, intent(in) :: t_density_matrix_z(:, :)
    !> @param fxc_contraction_matrix_s Pointer to the total density fxc contraction matrix data.
    real(c_double), contiguous, intent(out) :: fxc_contraction_matrix_s(:, :)
    !> @param fxc_contraction_matrix_z Pointer to the spin density fxc contraction matrix data.
    real(c_double), contiguous, intent(out) :: fxc_contraction_matrix_z(:, :)

    integer(c_int64_t) :: m, n, ldp_s, ldp_z, ldtp_s, ldtp_z, &
      & ldfxc_s, ldfxc_z

    m = size(density_matrix_s, 1, kind=c_int64_t)
    n = size(density_matrix_s, 2, kind=c_int64_t)
    ldp_s = size(density_matrix_s, 1, kind=c_int64_t)
    ldp_z = size(density_matrix_z, 1, kind=c_int64_t)
    ldtp_s = size(t_density_matrix_s, 1, kind=c_int64_t)
    ldtp_z = size(t_density_matrix_z, 1, kind=c_int64_t)
    ldfxc_s = size(fxc_contraction_matrix_s, 1, kind=c_int64_t)
    ldfxc_z = size(fxc_contraction_matrix_z, 1, kind=c_int64_t)

    call gauxc_integrator_eval_fxc_contraction_uks_c(status, integrator, m, n, &
      & density_matrix_s, ldp_s, density_matrix_z, ldp_z, &
      & t_density_matrix_s, ldtp_s, t_density_matrix_z, ldtp_z, &
      & fxc_contraction_matrix_s, ldfxc_s, fxc_contraction_matrix_z, ldfxc_z)
  end subroutine gauxc_integrator_eval_fxc_contraction_uks

end module gauxc_integrator