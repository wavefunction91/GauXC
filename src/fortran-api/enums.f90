! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining enumerations for GauXC
module gauxc_enums
  use iso_c_binding, only : c_int
  implicit none
  private

  public :: gauxc_radialquad_becke, gauxc_radialquad_mura_knowles, &
    & gauxc_radialquad_murray_handy_laming, gauxc_radialquad_treutler_ahlrichs
  enum, bind(c)
    !> @brief Becke radial quadrature
    enumerator :: gauxc_radialquad_becke
    !> @brief Mura-Knowles radial quadrature
    enumerator :: gauxc_radialquad_mura_knowles
    !> @brief Murray-Handy-Laming radial quadrature
    enumerator :: gauxc_radialquad_murray_handy_laming
    !> @brief Treutler-Ahlrichs radial quadrature
    enumerator :: gauxc_radialquad_treutler_ahlrichs
  end enum

  type :: gauxc_radialquad_enum
    integer(c_int) :: becke = gauxc_radialquad_becke
    integer(c_int) :: mura_knowles = gauxc_radialquad_mura_knowles
    integer(c_int) :: murray_handy_laming = gauxc_radialquad_murray_handy_laming
    integer(c_int) :: treutler_ahlrichs = gauxc_radialquad_treutler_ahlrichs
  end type gauxc_radialquad_enum
  type(gauxc_radialquad_enum), parameter, public :: gauxc_radialquad = gauxc_radialquad_enum()

  public :: gauxc_atomicgridsizedefault_finegrid, &
    & gauxc_atomicgridsizedefault_ultrafinegrid, &
    & gauxc_atomicgridsizedefault_superfinegrid, &
    & gauxc_atomicgridsizedefault_gm3, &
    & gauxc_atomicgridsizedefault_gm5
  enum, bind(c)
    !> @brief Fine grid size defaults for atomic grids
    enumerator :: gauxc_atomicgridsizedefault_finegrid
    !> @brief Ultrafine grid size defaults for atomic grids
    enumerator :: gauxc_atomicgridsizedefault_ultrafinegrid
    !> @brief Superfine grid size defaults for atomic grids
    enumerator :: gauxc_atomicgridsizedefault_superfinegrid
    !> @brief Treutler-Ahlrichs GM3 grid size defaults for atomic grids
    enumerator :: gauxc_atomicgridsizedefault_gm3
    !> @brief Treutler-Ahlrichs GM5 grid size defaults for atomic grids
    enumerator :: gauxc_atomicgridsizedefault_gm5
  end enum

  type :: gauxc_atomicgridsizedefault_enum
    integer(c_int) :: finegrid = gauxc_atomicgridsizedefault_finegrid
    integer(c_int) :: ultrafinegrid = gauxc_atomicgridsizedefault_ultrafinegrid
    integer(c_int) :: superfinegrid = gauxc_atomicgridsizedefault_superfinegrid
    integer(c_int) :: gm3 = gauxc_atomicgridsizedefault_gm3
    integer(c_int) :: gm5 = gauxc_atomicgridsizedefault_gm5
  end type gauxc_atomicgridsizedefault_enum
  type(gauxc_atomicgridsizedefault_enum), parameter, public :: &
    & gauxc_atomicgridsizedefault = gauxc_atomicgridsizedefault_enum()

  public :: gauxc_xcweightalg_notpartitioned, &
    & gauxc_xcweightalg_becke, &
    & gauxc_xcweightalg_ssf, &
    & gauxc_xcweightalg_lko
  enum, bind(c)
    !> @brief Exchange-correlation weight algorithm: Not partitioned
    enumerator :: gauxc_xcweightalg_notpartitioned
    !> @brief Exchange-correlation weight algorithm: Becke
    enumerator :: gauxc_xcweightalg_becke
    !> @brief Exchange-correlation weight algorithm: Stratmann-Scuseria-Frisch
    enumerator :: gauxc_xcweightalg_ssf
    !> @brief Exchange-correlation weight algorithm: Lauqua-Kuessman-Ochsenfeld
    enumerator :: gauxc_xcweightalg_lko
  end enum

  type :: gauxc_xcweightalg_enum
    integer(c_int) :: notpartitioned = gauxc_xcweightalg_notpartitioned
    integer(c_int) :: becke = gauxc_xcweightalg_becke
    integer(c_int) :: ssf = gauxc_xcweightalg_ssf
    integer(c_int) :: lko = gauxc_xcweightalg_lko
  end type gauxc_xcweightalg_enum
  type(gauxc_xcweightalg_enum), parameter, public :: &
    & gauxc_xcweightalg = gauxc_xcweightalg_enum()

  public :: gauxc_executionspace_host, gauxc_executionspace_device
  enum, bind(c)
    !> @brief Execute tasks on the host
    enumerator :: gauxc_executionspace_host
    !> @brief Execute tasks on the device
    enumerator :: gauxc_executionspace_device
  end enum

  type :: gauxc_executionspace_enum
    integer(c_int) :: host = gauxc_executionspace_host
    integer(c_int) :: device = gauxc_executionspace_device
  end type gauxc_executionspace_enum
  type(gauxc_executionspace_enum), parameter, public :: &
    & gauxc_executionspace = gauxc_executionspace_enum()

  public :: gauxc_supportedalg_xc, gauxc_supportedalg_den, gauxc_supportedalg_snlink
  enum, bind(c)
    !> @brief Supported algorithm: Exchange-correlation
    enumerator :: gauxc_supportedalg_xc
    !> @brief Supported algorithm: Density
    enumerator :: gauxc_supportedalg_den
    !> @brief Supported algorithm: SNLINK
    enumerator :: gauxc_supportedalg_snlink
  end enum

  type :: gauxc_supportedalg_enum
    integer(c_int) :: xc = gauxc_supportedalg_xc
    integer(c_int) :: den = gauxc_supportedalg_den
    integer(c_int) :: snlink = gauxc_supportedalg_snlink
  end type gauxc_supportedalg_enum
  type(gauxc_supportedalg_enum), parameter, public :: &
    & gauxc_supportedalg = gauxc_supportedalg_enum()

  public :: gauxc_pruningscheme_unpruned, &
    & gauxc_pruningscheme_robust, &
    & gauxc_pruningscheme_treutler
  enum, bind(c)
    !> @brief Pruning scheme: Unpruned atomic quadrature
    enumerator :: gauxc_pruningscheme_unpruned
    !> @brief Pruning scheme: The "Robust" scheme of Psi4
    enumerator :: gauxc_pruningscheme_robust
    !> @brief Pruning scheme: The Treutler-Ahlrichs scheme
    enumerator :: gauxc_pruningscheme_treutler
  end enum

  type :: gauxc_pruningscheme_enum
    integer(c_int) :: unpruned = gauxc_pruningscheme_unpruned
    integer(c_int) :: robust = gauxc_pruningscheme_robust
    integer(c_int) :: treutler = gauxc_pruningscheme_treutler
  end type gauxc_pruningscheme_enum
  type(gauxc_pruningscheme_enum), parameter, public :: &
    & gauxc_pruningscheme = gauxc_pruningscheme_enum()

end module gauxc_enums