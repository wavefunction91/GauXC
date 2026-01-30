! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining matrix operations for GauXC
module gauxc_matrix
  use iso_c_binding, only : c_ptr, c_null_ptr, c_size_t, c_double, c_f_pointer
  use gauxc_types, only : gauxc_header_type, gauxc_type_matrix
  use gauxc_status, only : gauxc_status_type
  implicit none
  private

  public :: &
    & gauxc_matrix_empty, &
    & gauxc_matrix_new, &
    & gauxc_matrix_resize, &
    & gauxc_matrix_set_zero, &
    & gauxc_matrix_rows, &
    & gauxc_matrix_cols, &
    & gauxc_matrix_data, &
    & gauxc_matrix_delete

  public :: &
    & gauxc_delete

  !> @brief C interoperable matrix type
  type, bind(c), public :: gauxc_matrix_type
    !> Internal header
    type(gauxc_header_type) :: header = gauxc_header_type(gauxc_type_matrix)
    !> Pointer to internal matrix data
    type(c_ptr) :: ptr = c_null_ptr
  end type gauxc_matrix_type

  interface
    !> @brief Create an empty matrix object
    function gauxc_matrix_empty(status) result(mat) bind(c)
      import :: gauxc_matrix_type, gauxc_status_type
      implicit none
      !> @param status Output status object
      type(gauxc_status_type), intent(out) :: status
      !> @return Empty matrix object
      type(gauxc_matrix_type) :: mat
    end function gauxc_matrix_empty

    !> @brief Create a new matrix object
    function gauxc_matrix_new(status, rows, cols) result(mat) bind(c)
      import :: c_size_t, gauxc_matrix_type, gauxc_status_type
      implicit none
      !> @param status Output status object
      type(gauxc_status_type), intent(out) :: status
      !> @param rows Number of rows
      integer(c_size_t), value :: rows
      !> @param cols Number of columns
      integer(c_size_t), value :: cols
      !> @return New matrix object
      type(gauxc_matrix_type) :: mat
    end function gauxc_matrix_new

    !> @brief Resize an existing matrix object
    subroutine gauxc_matrix_resize(status, mat, rows, cols) bind(c)
      import :: c_size_t, gauxc_matrix_type, gauxc_status_type
      implicit none
      !> @param status Output status object
      type(gauxc_status_type), intent(out) :: status
      !> @param mat Matrix object to resize
      type(gauxc_matrix_type), value :: mat
      !> @param rows New number of rows
      integer(c_size_t), value :: rows
      !> @param cols New number of columns
      integer(c_size_t), value :: cols
    end subroutine gauxc_matrix_resize

    !> @brief Set all elements of the Matrix to zero
    subroutine gauxc_matrix_set_zero(status, mat) bind(c)
      import :: gauxc_matrix_type, gauxc_status_type
      implicit none
      !> @param status Output status object
      type(gauxc_status_type), intent(out) :: status
      !> @param mat Matrix object to set to zero
      type(gauxc_matrix_type), value :: mat
    end subroutine gauxc_matrix_set_zero

    !> @brief Get the number of rows in the matrix
    function gauxc_matrix_rows(status, mat) result(rows) bind(c)
      import :: c_size_t, gauxc_matrix_type, gauxc_status_type
      implicit none
      !> @param status Output status object
      type(gauxc_status_type), intent(out) :: status
      !> @param mat Matrix object
      type(gauxc_matrix_type), value :: mat
      !> @return Number of rows
      integer(c_size_t) :: rows
    end function gauxc_matrix_rows

    !> @brief Get the number of columns in the matrix
    function gauxc_matrix_cols(status, mat) result(cols) bind(c)
      import :: c_size_t, gauxc_matrix_type, gauxc_status_type
      implicit none
      !> @param status Output status object
      type(gauxc_status_type), intent(out) :: status
      !> @param mat Matrix object
      type(gauxc_matrix_type), value :: mat
      !> @return Number of columns
      integer(c_size_t) :: cols
    end function gauxc_matrix_cols

    !> @brief Get a pointer to the internal data of the matrix
    function gauxc_matrix_data_c(status, mat) result(data_ptr) &
      & bind(c, name="gauxc_matrix_data")
      import :: c_ptr, gauxc_matrix_type, gauxc_status_type
      implicit none
      !> @param status Output status object
      type(gauxc_status_type), intent(out) :: status
      !> @param mat Matrix object
      type(gauxc_matrix_type), value :: mat
      !> @return Pointer to internal data
      type(c_ptr) :: data_ptr
    end function gauxc_matrix_data_c
  end interface

  interface gauxc_delete
    !> @brief Delete a matrix object and free its resources
    subroutine gauxc_matrix_delete(status, mat) bind(c)
      import :: gauxc_matrix_type, gauxc_status_type
      implicit none
      !> @param status Output status object
      type(gauxc_status_type), intent(out) :: status
      !> @param mat Matrix object to delete
      type(gauxc_matrix_type), intent(inout) :: mat
    end subroutine gauxc_matrix_delete
  end interface gauxc_delete

contains

  !> @brief Get a pointer to the internal data of the matrix
  function gauxc_matrix_data(status, mat) result(data)
    !> @param status Output status object
    type(gauxc_status_type), intent(out) :: status
    !> @param mat Matrix object
    type(gauxc_matrix_type), value :: mat
    !> @return Pointer to internal data
    real(c_double), pointer :: data(:, :)

    type(c_ptr) :: data_ptr
    integer(c_size_t) :: rows, cols

    data => null()
    data_ptr = gauxc_matrix_data_c(status, mat)
    if (status%code == 0) &
      cols = gauxc_matrix_cols(status, mat)
    if (status%code == 0) &
      rows = gauxc_matrix_rows(status, mat)
    if (status%code == 0) &
      call c_f_pointer(data_ptr, data, [cols, rows])
  end function gauxc_matrix_data

end module gauxc_matrix