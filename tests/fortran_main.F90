#include "gauxc/gauxc_config.f"

program gauxc_fortran_tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  use fortran_api_test, only : collect_fortran_api_suite
#ifdef GAUXC_HAS_MPI
  use mpi_f08, only : MPI_Initialized, MPI_Init, MPI_Finalized, MPI_Finalize, MPI_SUCCESS
#endif
  implicit none

  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'
#ifdef GAUXC_HAS_MPI
  integer :: mpi_ierr
#endif

#ifdef GAUXC_HAS_MPI
  call MPI_Init(mpi_ierr)
  if (mpi_ierr /= MPI_SUCCESS) error stop 1
#endif

  stat = 0

  testsuites = [ &
    new_testsuite("fortran-api", collect_fortran_api_suite) &
  ]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

#ifdef GAUXC_HAS_MPI
  call MPI_Finalize(mpi_ierr)
  if (mpi_ierr /= MPI_SUCCESS) error stop 1
#endif

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if

end program gauxc_fortran_tester
