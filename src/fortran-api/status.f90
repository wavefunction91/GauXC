! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining the status codes for GauXC
module gauxc_status
  use iso_c_binding, only : c_ptr, c_null_ptr, c_int, c_char, &
    & c_f_pointer, c_associated
  implicit none
  private

  public :: gauxc_status_message

  !> @brief C interoperable status type
  type, bind(c), public :: gauxc_status_type
    integer(c_int) :: code
    type(c_ptr) :: message = c_null_ptr
  end type gauxc_status_type

contains

  function gauxc_status_message(status) result(msg)
    type(gauxc_status_type), intent(in) :: status

    interface
      function strlen(str) result(len) bind(c)
        import :: c_ptr, c_int
        implicit none
        type(c_ptr), value :: str
        integer(c_int) :: len
      end function strlen
    end interface
    integer(c_int) :: slen
    character(kind=c_char, len=:), allocatable :: msg
    character(kind=c_char), pointer :: pmsg(:)

    if (c_associated(status%message)) then
      slen = strlen(status%message)
      call c_f_pointer(status%message, pmsg, [slen])
      msg = transfer(pmsg(1:slen), "")
    else
      msg = ""
    end if
  end function gauxc_status_message

end module gauxc_status