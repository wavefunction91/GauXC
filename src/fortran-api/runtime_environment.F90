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
#ifdef GAUXC_HAS_MPI
#ifdef GAUXC_HAS_MPI_F08
  use mpi_f08, only : MPI_Comm, MPI_COMM_WORLD
#else
  use mpi, only : MPI_Comm, MPI_COMM_WORLD
#endif
#endif
  implicit none
  private

  public :: &
    & gauxc_runtime_environment_new, &
    & gauxc_runtime_environment_delete, &
    & gauxc_runtime_environment_comm_rank, &
    & gauxc_runtime_environment_comm_size
#ifdef GAUXC_HAS_DEVICE
  public :: &
    & gauxc_device_runtime_environment_new
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
    function gauxc_runtime_environment_new_c &
#ifdef GAUXC_HAS_MPI
      & (status, comm) &
#else
      & (status) &
#endif
      & result(rt) bind(c, name="gauxc_runtime_environment_new_f")
      import :: gauxc_runtime_environment_type, gauxc_status_type, c_int
      implicit none
      type(gauxc_status_type), intent(out) :: status
#ifdef GAUXC_HAS_MPI
      integer(c_int), value :: comm
#endif
      type(gauxc_runtime_environment_type) :: rt
    end function gauxc_runtime_environment_new_c
  end interface

  interface gauxc_runtime_environment_new
    module procedure gauxc_runtime_environment_new
#ifdef GAUXC_HAS_MPI
    module procedure gauxc_runtime_environment_new_mpi
#endif
#ifdef GAUXC_HAS_MPI_F08
    module procedure gauxc_runtime_environment_new_mpi_f08
#endif
  end interface gauxc_runtime_environment_new

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

  interface
    function gauxc_device_runtime_environment_new_c &
#ifdef GAUXC_HAS_MPI
      & (status, comm, fill_fraction) &
#else
      & (status, fill_fraction) &
#endif
      & result(rt) bind(c, name="gauxc_device_runtime_environment_new_f")
      import :: gauxc_runtime_environment_type, gauxc_status_type, c_int, c_double
      implicit none
      type(gauxc_status_type), intent(out) :: status
#ifdef GAUXC_HAS_MPI
      integer(c_int), value :: comm
#endif
      real(c_double), value :: fill_fraction
      type(gauxc_runtime_environment_type) :: rt
    end function gauxc_device_runtime_environment_new_c

    function gauxc_device_runtime_environment_new_mem_c &
#ifdef GAUXC_HAS_MPI
      & (status, comm, mem, mem_sz) &
#else
      & (status, mem, mem_sz) &
#endif
      & result(rt) bind(c, name="gauxc_device_runtime_environment_new_mem_f")
      import :: gauxc_runtime_environment_type, gauxc_status_type, c_int, c_size_t, c_ptr
      implicit none
      type(gauxc_status_type), intent(out) :: status
#ifdef GAUXC_HAS_MPI
      integer(c_int), value :: comm
#endif
      type(c_ptr), value :: mem
      integer(c_size_t), value :: mem_sz
      type(gauxc_runtime_environment_type) :: rt
    end function gauxc_device_runtime_environment_new_mem_c
  end interface

#ifdef GAUXC_HAS_DEVICE
  interface gauxc_device_runtime_environment_new
    module procedure gauxc_device_runtime_environment_new
#ifdef GAUXC_HAS_MPI
    module procedure gauxc_device_runtime_environment_new_mpi
#endif
#ifdef GAUXC_HAS_MPI_F08
    module procedure gauxc_device_runtime_environment_new_mpi_f08
#endif

    module procedure gauxc_device_runtime_environment_new_mem
#ifdef GAUXC_HAS_MPI
    module procedure gauxc_device_runtime_environment_new_mem_mpi
#endif
#ifdef GAUXC_HAS_MPI_F08
    module procedure gauxc_device_runtime_environment_new_mem_mpi_f08
#endif
  end interface gauxc_device_runtime_environment_new
#endif

contains

  function gauxc_runtime_environment_new(status) result(rt)
    type(gauxc_status_type), intent(out) :: status
    type(gauxc_runtime_environment_type) :: rt

#ifdef GAUXC_HAS_MPI
#ifdef GAUXC_HAS_MPI_F08
    rt = gauxc_runtime_environment_new_c(status, MPI_COMM_WORLD%mpi_val)
#else
    rt = gauxc_runtime_environment_new_c(status, MPI_COMM_WORLD)
#endif
#else
    rt = gauxc_runtime_environment_new_c(status)
#endif
  end function gauxc_runtime_environment_new

#ifdef GAUXC_HAS_MPI
  function gauxc_runtime_environment_new_mpi(status, comm) result(rt)
    type(gauxc_status_type), intent(out) :: status
    integer, intent(in) :: comm
    type(gauxc_runtime_environment_type) :: rt

    rt = gauxc_runtime_environment_new_c(status, comm)
  end function gauxc_runtime_environment_new_mpi
#endif

#ifdef GAUXC_HAS_MPI_F08
  function gauxc_runtime_environment_new_mpi_f08(status, comm) result(rt)
    type(gauxc_status_type), intent(out) :: status
    type(MPI_Comm), intent(in) :: comm
    type(gauxc_runtime_environment_type) :: rt

    rt = gauxc_runtime_environment_new_c(status, comm%mpi_val)
  end function gauxc_runtime_environment_new_mpi_f08
#endif

#ifdef GAUXC_HAS_DEVICE
  function gauxc_device_runtime_environment_new &
    & (status, fill_fraction) &
    & result(rt)
    type(gauxc_status_type), intent(out) :: status
    real(c_double), value :: fill_fraction
    type(gauxc_runtime_environment_type) :: rt

#ifdef GAUXC_HAS_MPI
#ifdef GAUXC_HAS_MPI_F08
    rt = gauxc_device_runtime_environment_new_c(status, MPI_COMM_WORLD%mpi_val, fill_fraction)
#else
    rt = gauxc_device_runtime_environment_new_c(status, MPI_COMM_WORLD, fill_fraction)
#endif
#else
    rt = gauxc_device_runtime_environment_new_c(status, fill_fraction)
#endif
  end function gauxc_device_runtime_environment_new

#ifdef GAUXC_HAS_MPI
  function gauxc_device_runtime_environment_new_mpi &
    & (status, comm, fill_fraction) &
    & result(rt)
    type(gauxc_status_type), intent(out) :: status
    integer, intent(in) :: comm
    real(c_double), value :: fill_fraction
    type(gauxc_runtime_environment_type) :: rt

    rt = gauxc_device_runtime_environment_new_c(status, comm, fill_fraction)
  end function gauxc_device_runtime_environment_new_mpi
#endif

#ifdef GAUXC_HAS_MPI_F08
  function gauxc_device_runtime_environment_new_mpi_f08 &
    & (status, comm, fill_fraction) &
    & result(rt)
    type(gauxc_status_type), intent(out) :: status
    type(MPI_Comm), intent(in) :: comm
    real(c_double), value :: fill_fraction
    type(gauxc_runtime_environment_type) :: rt

    rt = gauxc_device_runtime_environment_new_c(status, comm%mpi_val, fill_fraction)
  end function gauxc_device_runtime_environment_new_mpi_f08
#endif

  function gauxc_device_runtime_environment_new_mem &
    & (status, mem, mem_sz) &
    & result(rt)
    type(gauxc_status_type), intent(out) :: status
    type(c_ptr), value :: mem
    integer(c_size_t), value :: mem_sz
    type(gauxc_runtime_environment_type) :: rt

#ifdef GAUXC_HAS_MPI
#ifdef GAUXC_HAS_MPI_F08
    rt = gauxc_device_runtime_environment_new_mem_c(status, MPI_COMM_WORLD%mpi_val, mem, mem_sz)
#else
    rt = gauxc_device_runtime_environment_new_mem_c(status, MPI_COMM_WORLD, mem, mem_sz)
#endif
#else
    rt = gauxc_device_runtime_environment_new_mem_c(status, mem, mem_sz)
#endif
  end function gauxc_device_runtime_environment_new_mem

#ifdef GAUXC_HAS_MPI
  function gauxc_device_runtime_environment_new_mem_mpi &
    & (status, comm, mem, mem_sz) &
    & result(rt)
    type(gauxc_status_type), intent(out) :: status
    integer, intent(in) :: comm
    type(c_ptr), value :: mem
    integer(c_size_t), value :: mem_sz
    type(gauxc_runtime_environment_type) :: rt

    rt = gauxc_device_runtime_environment_new_mem_c(status, comm, mem, mem_sz)
  end function gauxc_device_runtime_environment_new_mem_mpi
#endif

#ifdef GAUXC_HAS_MPI_F08
  function gauxc_device_runtime_environment_new_mem_mpi_f08 &
    & (status, comm, mem, mem_sz) &
    & result(rt)
    type(gauxc_status_type), intent(out) :: status
    type(MPI_Comm), intent(in) :: comm
    type(c_ptr), value :: mem
    integer(c_size_t), value :: mem_sz
    type(gauxc_runtime_environment_type) :: rt

    rt = gauxc_device_runtime_environment_new_mem_c(status, comm%mpi_val, mem, mem_sz)
  end function gauxc_device_runtime_environment_new_mem_mpi_f08
#endif
#endif
end module gauxc_runtime_environment