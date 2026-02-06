! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

#include "gauxc/gauxc_config.f"

!> @brief Module defining the runtime environment for GauXC
module gauxc_runtime_environment
  use iso_c_binding, only : c_ptr, c_null_ptr, c_int, c_double, c_size_t
  use gauxc_status, only : gauxc_status_type
  use gauxc_types, only : gauxc_header_type, gauxc_type_runtime_environment
  implicit none
  private

  public :: &
    & gauxc_runtime_environment_new, &
    & gauxc_runtime_environment_delete, &
    & gauxc_runtime_environment_comm_rank, &
    & gauxc_runtime_environment_comm_size
#ifdef GAUXC_HAS_DEVICE
  public :: &
    & gauxc_device_runtime_environment_new, &
    & gauxc_device_runtime_environment_new_mem
#endif

  public :: &
    & gauxc_delete

  type, bind(c), public :: gauxc_runtime_environment_type
    type(gauxc_header_type) :: hdr = gauxc_header_type(gauxc_type_runtime_environment)
    type(c_ptr) :: ptr = c_null_ptr
#ifdef GAUXC_HAS_DEVICE
    type(c_ptr) :: device_ptr = c_null_ptr
#endif
  end type gauxc_runtime_environment_type

  interface
    function gauxc_runtime_environment_new &
#ifdef GAUXC_HAS_MPI
      & (status, comm) &
#else
      & (status) &
#endif
      & result(rt) bind(c)
      import :: gauxc_runtime_environment_type, gauxc_status_type, c_int
      implicit none
      type(gauxc_status_type), intent(out) :: status
#ifdef GAUXC_HAS_MPI
      integer(c_int), value :: comm
#endif
      type(gauxc_runtime_environment_type) :: rt
    end function gauxc_runtime_environment_new
  end interface

  interface gauxc_delete
    subroutine gauxc_runtime_environment_delete(status, rt) bind(c)
      import :: gauxc_runtime_environment_type, gauxc_status_type
      implicit none
      type(gauxc_status_type), intent(out) :: status
      type(gauxc_runtime_environment_type), intent(inout) :: rt
    end subroutine gauxc_runtime_environment_delete
  end interface gauxc_delete

  interface
    function gauxc_runtime_environment_comm_rank(status, rt) result(comm_size) bind(c)
      import :: gauxc_runtime_environment_type, gauxc_status_type, c_int
      implicit none
      type(gauxc_status_type), intent(out) :: status
      type(gauxc_runtime_environment_type), value :: rt
      integer(c_int) :: comm_size
    end function gauxc_runtime_environment_comm_rank
  end interface

  interface
    function gauxc_runtime_environment_comm_size(status, rt) result(comm_size) bind(c)
      import :: gauxc_runtime_environment_type, gauxc_status_type, c_int
      implicit none
      type(gauxc_status_type), intent(out) :: status
      type(gauxc_runtime_environment_type), value :: rt
      integer(c_int) :: comm_size
    end function gauxc_runtime_environment_comm_size
  end interface

#ifdef GAUXC_HAS_DEVICE
  interface gauxc_device_runtime_environment_new
    function gauxc_device_runtime_environment_new &
#ifdef GAUXC_HAS_MPI
      & (status, comm, fill_fraction) &
#else
      & (status, fill_fraction) &
#endif
      & result(rt) bind(c)
      import :: gauxc_runtime_environment_type, gauxc_status_type, c_int, c_double
      implicit none
      type(gauxc_status_type), intent(out) :: status
#ifdef GAUXC_HAS_MPI
      integer(c_int), value :: comm
#endif
      real(c_double), value :: fill_fraction
      type(gauxc_runtime_environment_type) :: rt
    end function gauxc_device_runtime_environment_new

    function gauxc_device_runtime_environment_new_mem &
#ifdef GAUXC_HAS_MPI
      & (status, comm, mem, mem_sz) &
#else
      & (status, mem, mem_sz) &
#endif
      & result(rt) bind(c)
      import :: gauxc_runtime_environment_type, gauxc_status_type, c_int, c_size_t
      implicit none
      type(gauxc_status_type), intent(out) :: status
#ifdef GAUXC_HAS_MPI
      integer(c_int), value :: comm
#endif
      type(c_ptr), value :: mem
      integer(c_size_t), value :: mem_sz
      type(gauxc_runtime_environment_type) :: rt
    end function gauxc_device_runtime_environment_new_mem
  end interface gauxc_device_runtime_environment_new
#endif

contains

end module gauxc_runtime_environment