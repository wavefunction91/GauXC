! GauXC Copyright (c) 2020-2024, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of
! any required approvals from the U.S. Dept. of Energy).
!
! (c) 2024-2025, Microsoft Corporation
!
! All rights reserved.
!
! See LICENSE.txt for details

!> @brief Module defining exchange-correlation functionals for GauXC
module gauxc_xc_functional
  use iso_c_binding, only : c_ptr, c_null_ptr, c_char, c_int, c_bool, c_null_char
  use gauxc_status, only : gauxc_status_type
  use gauxc_types, only : gauxc_header_type, gauxc_type_functional
  implicit none
  private

  public :: &
    & gauxc_functional_from_enum, &
    & gauxc_functional_from_string, &
    & gauxc_functional_delete

  public :: &
    & gauxc_delete

  enum, bind(c)
    !> @brief Slater exchange & Vosko, Wilk & Nusair correlation (VWN3)
    !> - libxc names: LDA_X (id=1) and LDA_C_VWN_3 (id=30)
    enumerator :: gauxc_functional_svwn3

    !> @brief Slater exchange & Vosko, Wilk & Nusair correlation (VWN5)
    !> - libxc names: LDA_X (id=1) and LDA_C_VWN (id=7)
    !> - xcfun names: SLATERX and VWN5C
    enumerator :: gauxc_functional_svwn5

    !> @brief Becke 88 exchange & Lee, Yang & Parr correlation
    !> - libxc names: GGA_X_B88 (id=106) and GGA_C_LYP (id=131)
    !> - xcfun names: BECKEX and LYPC
    enumerator :: gauxc_functional_blyp

    !> @brief Becke 88 exchange & Lee, Yang & Parr correlation, 3-parameter hybrid
    !> - libxc name: HYB_GGA_XC_B3LYP (id=402)
    enumerator :: gauxc_functional_b3lyp

    !> @brief Perdew-Burke-Ernzerhof exchange & correlation
    !> - libxc names: GGA_X_PBE (id=101) & GGA_C_PBE (id=130)
    !> - xcfun names: PBEEX and PBEC
    enumerator :: gauxc_functional_pbe

    !> @brief revised Perdew-Burke-Ernzerhof exchange & original PBE correlation
    !> - libxc names: GGA_X_PBE_R (id=102) & GGA_C_PBE (id=130)
    enumerator :: gauxc_functional_revpbe

    !> @brief Perdew-Burke-Ernzerhof exchange & correlation, 1-parameter hybrid
    !> - libxc name: HYB_GGA_XC_PBEH (id=406)
    enumerator :: gauxc_functional_pbe0

    !> @brief Strongly constrained and appropriately normed (SCAN) meta-GGA
    !> - libxc names: MGGA_X_SCAN (id=263) & MGGA_C_SCAN (id=267)
    enumerator :: gauxc_functional_scan

    !> @brief Regularized & restored strongly constrained and appropriately normed (R2SCAN) meta-GGA
    !> - libxc names: MGGA_X_R2SCAN (id=497) & MGGA_C_R2SCAN (id=498) 
    enumerator :: gauxc_functional_r2scan

    !> @brief Regularized & restored strongly constrained and appropriately normed (R2SCAN) meta-GGA, deorbitalized version
    !> - libxc names: MGGA_X_R2SCANL (id=718) & MGGA_C_R2SCANL (id=719)
    enumerator :: gauxc_functional_r2scanl

    !> @brief Minnesota 2006 hybrid functional
    !> - libxc names: HYB_MGGA_X_M06_2X (id=450) & MGGA_C_M06_2X (id=236)
    !> - xcfun names: M062X and M062C
    enumerator :: gauxc_functional_m062x

    !> @brief Perdew, Kurth, Zupan, and Blaha
    !> - libxc names: MGGA_X_PKZB (id=213) & MGGA_C_PKZB (id=239)
    enumerator :: gauxc_functional_pkzb

    !> @brief epc17(-1): electron-proton correlation 2017
    !> - libxc name: LDA_C_EPC17 (id=328)
    enumerator :: gauxc_functional_epc17_1

    !> @brief epc17-2: electron-proton correlation 2017 for proton affinities
    !> - libxc name: LDA_C_EPC17_2 (id=329)
    enumerator :: gauxc_functional_epc17_2

    !> @brief epc18-1: electron-proton correlation 2018
    !> - libxc name: LDA_C_EPC18_1 (id=330)
    enumerator :: gauxc_functional_epc18_1

    !> @brief epc18-2: electron-proton correlation 2018 for proton affinities
    !> - libxc name: LDA_C_EPC18_2 (id=331)
    enumerator :: gauxc_functional_epc18_2

    !> @brief Grimme's parametrization of the B97 functional, original D2 variant
    !> - libxc name: GGA_XC_B97_D (id=170)
    enumerator :: gauxc_functional_b97d

    !> @brief Grimme's parametrization of the B97 functional, D3(0) variant
    !> - libxc name: GGA_XC_B97_D (id=170)
    enumerator :: gauxc_functional_b97d3zero

    !> @brief Coulomb-attenuating method range-separated hybrid functional
    !> - libxc name: HYB_GGA_XC_CAM_B3LYP (id=433)
    enumerator :: gauxc_functional_camb3lyp

    !> @brief Slater exchange
    !> - libxc name: LDA_X (id=1)
    !> - xcfun name: SLATERX
    enumerator :: gauxc_functional_lda

    !> @brief Minnesota 2006 meta-GGA functional
    !> - libxc names: MGGA_X_M06_L (id=203) & MGGA_C_M06_L (id=233)
    enumerator :: gauxc_functional_m06l

    !> @brief Strongly constrained and appropriately normed (SCAN) meta-GGA, 1-parameter hybrid
    !> - libxc names: HYB_MGGA_X_SCAN0 (id=264) & MGGA_C_SCAN (id=267)
    enumerator :: gauxc_functional_scan0

    !> @brief Slater exchange & Perdew, Wang 92 correlation
    !> - libxc names: LDA_X (id=1) and LDA_C_PW (id=12)
    enumerator :: gauxc_functional_spw92

    !> @brief Tao, Perdew, Staroverov & Scuseria meta-GGA
    !> - libxc names: MGGA_X_TPSS (id=202) and MGGA_C_TPSS (id=231)
    enumerator :: gauxc_functional_tpss

    !> @brief Tao, Perdew, Staroverov & Scuseria meta-GGA, 1-parameter hybrid
    !> - libxc name: HYB_MGGA_XC_TPSSH (id=457)
    enumerator :: gauxc_functional_tpssh

    !> @brief Tao, Perdew, Staroverov & Scuseria meta-GGA, 1-parameter hybrid
    !> - libxc name: HYB_MGGA_XC_TPSS0 (id=396)
    enumerator :: gauxc_functional_tpss0

    !> @brief Vosko, Wilk & Nusair correlation (VWN3)
    !> - libxc name: LDA_C_VWN_3 (id=30)
    enumerator :: gauxc_functional_vwn3

    !> @brief Vosko, Wilk & Nusair correlation (VWN5)
    !> - libxc name: LDA_C_VWN (id=7)
    enumerator :: gauxc_functional_vwn5

    !> @brief HJS screened PBE exchange & original PBE correlation
    !> - libxc names: GGA_X_HJS_PBE (id=525) & GGA_C_PBE (id=130)
    enumerator :: gauxc_functional_lrcwpbe

    !> @brief HJS screened PBE exchange & original PBE correlation, hybrid version
    !> - libxc name: HYB_GGA_XC_HJS_PBE (id=429)
    enumerator :: gauxc_functional_lrcwpbeh

    !> @brief Becke 88 exchange and Perdew 86 correlation
    !> - libxc names: GGA_X_B88 (id=106) and GGA_C_P86 (id=132)
    enumerator :: gauxc_functional_bp86

    !> @brief Heyd-Scuseria-Ernzerhof screened hybrid functional (HSE03)
    !> - libxc name: HYB_GGA_XC_HSE03 (id=427)
    enumerator :: gauxc_functional_hse03

    !> @brief Heyd-Scuseria-Ernzerhof screened hybrid functional (HSE06)
    !> - libxc name: HYB_GGA_XC_HSE06 (id=428)
    enumerator :: gauxc_functional_hse06

    !> @brief Revised B3LYP
    !> - libxc name: HYB_GGA_XC_REVB3LYP (id=454)
    enumerator :: gauxc_functional_revb3lyp

    !> @brief revised Perdew-Burke-Ernzerhof exchange & original PBE correlation, hybrid version
    !> - libxc names: GGA_X_PBE_R (id=102) & GGA_C_PBE (id=130)
    enumerator :: gauxc_functional_revpbe0

    !> @brief revised Tao, Perdew, Staroverov & Scuseria
    !> - libxc names: MGGA_X_REVTPSS (id=212) & MGGA_C_REVTPSS (id=241)
    enumerator :: gauxc_functional_revtpss

    !> @brief revTPSSh
    !> - libxc name: HYB_MGGA_XC_REVTPSSH (id=458)
    enumerator :: gauxc_functional_revtpssh

    !> @brief Perdew-Wang 91 exchange and correlation
    !> - libxc name: GGA_X_PW91 (id=109) and GGA_C_PW91 (id=134)
    enumerator :: gauxc_functional_pw91

    !> @brief mBEEF exchange and Perdew, Burke & Ernzerhof SOL
    !> - libxc names: MGGA_X_MBEEF (id=249) and GGA_C_PBE_SOL (id=133)
    enumerator :: gauxc_functional_mbeef

    !> @brief The original (ACM, B3PW91) hybrid of Becke
    !> - libxc name: HYB_GGA_XC_B3PW91 (id=401)
    enumerator :: gauxc_functional_b3pw91

    !> @brief O3LYP
    !> - libxc name: HYB_GGA_XC_O3LYP (id=404)
    enumerator :: gauxc_functional_o3lyp

    !> @brief Handy & Cohen OPTX 01 exchange and Lee, Yang & Parr correlation
    !> - libxc names: GGA_X_OPTX (id=110) & GGA_C_LYP (id=131)
    enumerator :: gauxc_functional_olyp

    !> @brief Handy & Cohen OPTX 01 exchange and Perdew, Burke & Ernzerhof correlation
    !> - libxc names: GGA_X_OPTX (id=110) & GGA_C_PBE (id=130)
    enumerator :: gauxc_functional_opbe

    !> @brief mPW1K
    !> - libxc name: HYB_GGA_XC_MPW1K (id=405)
    enumerator :: gauxc_functional_mpw1k

    !> @brief Revised Perdew-Burke-Ernzerhof exchange by Hammer, Hansen, and Norskov
    !> - libxc name: GGA_X_RPBE (id=117)
    enumerator :: gauxc_functional_rpbe

    !> @brief Becke 88 exchange
    !> - libxc name: GGA_X_B88 (id=106)
    enumerator :: gauxc_functional_b88

    !> @brief modified Perdew-Wang 91 exchange by Adamo & Barone
    !> - libxc name: GGA_X_MPW91 (id=119)
    enumerator :: gauxc_functional_mpw91

    !> @brief Regularized strongly constrained and appropriately normed (RSCAN) meta-GGA by Bartok and Yates
    !> - libxc names: MGGA_X_RSCAN (id=493) and MGGA_C_RSCAN (id=494)
    enumerator :: gauxc_functional_rscan

    !> @brief CAM version of B3LYP, tuned for excitations and properties
    !> - libxc name: HYB_GGA_XC_TUNED_CAM_B3LYP (id=434)
    enumerator :: gauxc_functional_tunedcamb3lyp

    !> @brief wB97 range-separated functional
    !> - libxc name: HYB_GGA_XC_WB97 (id=463)
    enumerator :: gauxc_functional_wb97

    !> @brief wB97X range-separated functional
    !> - libxc name: HYB_GGA_XC_WB97X (id=464)
    enumerator :: gauxc_functional_wb97x

    !> @brief wB97X-D range-separated functional
    !> - libxc name: HYB_GGA_XC_WB97X_D (id=471)
    enumerator :: gauxc_functional_wb97xd

    !> @brief wB97X-D3 range-separated functional
    !> - libxc name: HYB_GGA_XC_WB97X_D3 (id=399)
    enumerator :: gauxc_functional_wb97xd3

    !> @brief Long-range corrected PBE (LC-wPBE) by Vydrov and Scuseria
    !> - libxc name: HYB_GGA_XC_LC_WPBE (id=478)
    enumerator :: gauxc_functional_lcwpbe

    !> @brief X3LYP
    !> - libxc name: HYB_GGA_XC_X3LYP (id=411)
    enumerator :: gauxc_functional_x3lyp

    !> @brief XLYP
    !> - libxc name: GGA_XC_XLYP (id=166)
    enumerator :: gauxc_functional_xlyp

    !> @brief BHandH i.e. BHLYP
    !> - libxc name: HYB_GGA_XC_BHANDH (id=435)
    enumerator :: gauxc_functional_bhandh

    !> @brief Boese-Martin for kinetics
    !> - libxc names: HYB_MGGA_X_BMK (id=279) & GGA_C_BMK (id=280)
    enumerator :: gauxc_functional_bmk

    !> @brief Becke 88 exchange and Perdew 86 based on VWN5 correlation, with more accurate value for ftilde
    !> - libxc names: GGA_X_B88 (id=106) & GGA_C_P86VWN_FT (id=253)
    enumerator :: gauxc_functional_bp86vwn

    !> @brief Mixture of PW86 with BC95
    !> - libxc name: HYB_MGGA_XC_PW86B95 (id=442)
    enumerator :: gauxc_functional_pw86b95

    !> @brief Perdew & Wang 86 exchange and PBE correlation
    !> - libxc names: GGA_X_PW86 (id=108) & GGA_C_PBE (id=130)
    enumerator :: gauxc_functional_pw86pbe

    !> @brief r2SCAN0: r2SCAN hybrid like PBE0 with 25% exact exchange
    !> - libxc name: HYB_MGGA_XC_R2SCAN0 (id=660)
    enumerator :: gauxc_functional_r2scan0

    !> @brief r2SCANh: r2SCAN hybrid like TPSSh with 10% exact exchange
    !> - libxc name: HYB_MGGA_XC_R2SCANH (id=659)
    enumerator :: gauxc_functional_r2scanh

    !> @brief r2SCAN50: r2SCAN hybrid like BHLYP with 50% exact exchange
    !> - libxc name: HYB_MGGA_XC_R2SCAN50 (id=661)
    enumerator :: gauxc_functional_r2scan50

    !> @brief Minnesota 2005 hybrid functional
    !> - libxc names: HYB_MGGA_X_M05 (id=438) & MGGA_C_M05 (id=237)
    enumerator :: gauxc_functional_m05

    !> @brief Minnesota 2008 hybrid functional
    !> - libxc names: HYB_MGGA_X_M06 (id=449) & MGGA_C_M06 (id=235)
    enumerator :: gauxc_functional_m06

    !> @brief Minnesota M08 hybrid functional
    !> - libxc names: HYB_MGGA_X_M08_HX (id=295) & MGGA_C_M08_HX (id=78)
    enumerator :: gauxc_functional_m08hx

    !> @brief Minnesota M08-SO hybrid exchange functional
    !> - libxc names: HYB_MGGA_X_M08_SO (id=296) & MGGA_C_M08_SO (id=77)
    enumerator :: gauxc_functional_m08so

    !> @brief Minnesota M05-2X hybrid exchange functional
    !> - libxc names: HYB_MGGA_X_M05_2X (id=439) & MGGA_C_M05_2X (id=238)
    enumerator :: gauxc_functional_m052x

    !> @brief Minnesota M06-SX short-range hybrid exchange functional
    !> - libxc names: HYB_MGGA_X_M06_SX (id=310) & MGGA_C_M06_SX (id=311)
    enumerator :: gauxc_functional_m06sx

    !> @brief Minnesota CF22D hybrid exchange functional
    !> - libxc names: HYB_MGGA_X_CF22D (id=340) & MGGA_C_CF22D (id=341)
    enumerator :: gauxc_functional_cf22d

    !> @brief Hybrid based on SOGGA11 form
    !> - libxc names: HYB_GGA_X_SOGGA11_X (id=426) & GGA_C_SOGGA11_X (id=159)
    enumerator :: gauxc_functional_sogga11x

    !> @brief Minnesota M06-HF hybrid exchange functional
    !> - libxc names: HYB_MGGA_X_M06_HF (id=444) & MGGA_C_M06_HF (id=234)
    enumerator :: gauxc_functional_m06hf

    !> @brief Minnesota M11 hybrid exchange functional
    !> - libxc names: HYB_MGGA_X_M11 (id=297) & MGGA_C_M11 (id=76)
    enumerator :: gauxc_functional_m11

    !> @brief Minnesota MN12-L exchange functional
    !> - libxc names: MGGA_X_MN12_L (id=227) & MGGA_C_MN12_L (id=74)
    enumerator :: gauxc_functional_mn12l

    !> @brief Minnesota MN12-SX hybrid exchange functional
    !> - libxc names: HYB_MGGA_X_MN12_SX (id=248) & MGGA_C_MN12_SX (id=73)
    enumerator :: gauxc_functional_mn12sx

    !> @brief Minnesota MN15 correlation functional
    !> - libxc names: HYB_MGGA_X_MN15 (id=268) & MGGA_C_MN15 (id=269)
    enumerator :: gauxc_functional_mn15

    !> @brief Minnesota MN15-L exchange functional
    !> - libxc names: MGGA_X_MN15_L (id=260) & MGGA_C_MN15_L (id=261)
    enumerator :: gauxc_functional_mn15l

    !> @brief Revised Minnesota 2006 meta-GGA functional
    !> - libxc names: MGGA_X_REVM06_L (id=293) & MGGA_C_REVM06_L (id=294)
    enumerator :: gauxc_functional_revm06l
  end enum

  type :: gauxc_functional_enum
    integer(c_int) :: svwn3 = gauxc_functional_vwn3
    integer(c_int) :: svwn5 = gauxc_functional_vwn5
    integer(c_int) :: blyp = gauxc_functional_blyp
    integer(c_int) :: b3lyp = gauxc_functional_b3lyp
    integer(c_int) :: pbe = gauxc_functional_pbe
    integer(c_int) :: revpbe = gauxc_functional_revpbe
    integer(c_int) :: pbe0 = gauxc_functional_pbe0
    integer(c_int) :: scan = gauxc_functional_scan
    integer(c_int) :: r2scan = gauxc_functional_r2scan
    integer(c_int) :: r2scanl = gauxc_functional_r2scanl
    integer(c_int) :: m062x = gauxc_functional_m062x
    integer(c_int) :: pkzb = gauxc_functional_pkzb
    integer(c_int) :: epc17_1 = gauxc_functional_epc17_1
    integer(c_int) :: epc17_2 = gauxc_functional_epc17_2
    integer(c_int) :: epc18_1 = gauxc_functional_epc18_1
    integer(c_int) :: epc18_2 = gauxc_functional_epc18_2
    integer(c_int) :: b97d = gauxc_functional_b97d
    integer(c_int) :: b97d3zero = gauxc_functional_b97d3zero
    integer(c_int) :: camb3lyp = gauxc_functional_camb3lyp
    integer(c_int) :: lda = gauxc_functional_lda
    integer(c_int) :: m06l = gauxc_functional_m06l
    integer(c_int) :: scan0 = gauxc_functional_scan0
    integer(c_int) :: spw92 = gauxc_functional_spw92
    integer(c_int) :: tpss = gauxc_functional_tpss
    integer(c_int) :: tpssh = gauxc_functional_tpssh
    integer(c_int) :: tpss0 = gauxc_functional_tpss0
    integer(c_int) :: vwn3 = gauxc_functional_vwn3
    integer(c_int) :: vwn5 = gauxc_functional_vwn5
    integer(c_int) :: lrcwpbe = gauxc_functional_lrcwpbe
    integer(c_int) :: lrcwpbeh = gauxc_functional_lrcwpbeh
    integer(c_int) :: bp86 = gauxc_functional_bp86
    integer(c_int) :: hse03 = gauxc_functional_hse03
    integer(c_int) :: hse06 = gauxc_functional_hse06
    integer(c_int) :: revb3lyp = gauxc_functional_revb3lyp
    integer(c_int) :: revpbe0 = gauxc_functional_revpbe0
    integer(c_int) :: revtpss = gauxc_functional_revtpss
    integer(c_int) :: revtpssh = gauxc_functional_revtpssh
    integer(c_int) :: pw91 = gauxc_functional_pw91
    integer(c_int) :: mbeef = gauxc_functional_mbeef
    integer(c_int) :: b3pw91 = gauxc_functional_b3pw91
    integer(c_int) :: o3lyp = gauxc_functional_o3lyp
    integer(c_int) :: olyp = gauxc_functional_olyp
    integer(c_int) :: opbe = gauxc_functional_opbe
    integer(c_int) :: mpw1k = gauxc_functional_mpw1k
    integer(c_int) :: rpbe = gauxc_functional_rpbe
    integer(c_int) :: b88 = gauxc_functional_b88
    integer(c_int) :: mpw91 = gauxc_functional_mpw91
    integer(c_int) :: rscan = gauxc_functional_rscan
    integer(c_int) :: tunedcamb3lyp = gauxc_functional_tunedcamb3lyp
    integer(c_int) :: wb97 = gauxc_functional_wb97
    integer(c_int) :: wb97x = gauxc_functional_wb97x
    integer(c_int) :: wb97xd = gauxc_functional_wb97xd
    integer(c_int) :: wb97xd3 = gauxc_functional_wb97xd3
    integer(c_int) :: lcwpbe = gauxc_functional_lcwpbe
    integer(c_int) :: x3lyp = gauxc_functional_x3lyp
    integer(c_int) :: xlyp = gauxc_functional_xlyp
    integer(c_int) :: bhandh = gauxc_functional_bhandh
    integer(c_int) :: bmk = gauxc_functional_bmk
    integer(c_int) :: bp86vwn = gauxc_functional_bp86vwn
    integer(c_int) :: pw86b95 = gauxc_functional_pw86b95
    integer(c_int) :: pw86pbe = gauxc_functional_pw86pbe
    integer(c_int) :: r2scan0 = gauxc_functional_r2scan0
    integer(c_int) :: r2scanh = gauxc_functional_r2scanh
    integer(c_int) :: r2scan50 = gauxc_functional_r2scan50
    integer(c_int) :: m05 = gauxc_functional_m05
    integer(c_int) :: m06 = gauxc_functional_m06
    integer(c_int) :: m08hx = gauxc_functional_m08hx
    integer(c_int) :: m08so = gauxc_functional_m08so
    integer(c_int) :: m052x = gauxc_functional_m052x
    integer(c_int) :: m06sx = gauxc_functional_m06sx
    integer(c_int) :: cf22d = gauxc_functional_cf22d
    integer(c_int) :: sogga11x = gauxc_functional_sogga11x
    integer(c_int) :: m06hf = gauxc_functional_m06hf
    integer(c_int) :: m11 = gauxc_functional_m11
    integer(c_int) :: mn12l = gauxc_functional_mn12l
    integer(c_int) :: mn12sx = gauxc_functional_mn12sx
    integer(c_int) :: mn15 = gauxc_functional_mn15
    integer(c_int) :: mn15l = gauxc_functional_mn15l
    integer(c_int) :: revm06l = gauxc_functional_revm06l
  end type gauxc_functional_enum
  type(gauxc_functional_enum), public, parameter :: gauxc_functional = gauxc_functional_enum()

  !> @brief GAUXC XC functional handle
  type, bind(c), public :: gauxc_functional_type
    type(gauxc_header_type) :: hdr = gauxc_header_type(gauxc_type_functional)
    type(c_ptr) :: ptr = c_null_ptr
  end type gauxc_functional_type

  interface
    !> @brief Create a GauXCFunctional from a string specification
    function gauxc_functional_from_string_c(status, functional_spec, polarized) &
      & result(func) bind(c, name="gauxc_functional_from_string")
      import :: gauxc_status_type, gauxc_functional_type, c_char, c_bool
      type(gauxc_status_type), intent(out) :: status
      character(kind=c_char), intent(in) :: functional_spec(*)
      logical(c_bool), value :: polarized
      type(gauxc_functional_type) :: func
    end function gauxc_functional_from_string_c

    !> @brief Create a GauXCFunctional from an enum specification
    function gauxc_functional_from_enum(status, functional_enum, polarized) &
      & result(func) bind(c)
      import :: gauxc_status_type, gauxc_functional_type, c_int, c_bool
      type(gauxc_status_type), intent(out) :: status
      integer(c_int), value :: functional_enum
      logical(c_bool), value :: polarized
      type(gauxc_functional_type) :: func
    end function gauxc_functional_from_enum
  end interface

  interface gauxc_delete
    !> @brief Destroy a GauXCFunctional
    subroutine gauxc_functional_delete(status, func) bind(c)
      import :: gauxc_status_type, gauxc_functional_type
      type(gauxc_status_type), intent(out) :: status
      type(gauxc_functional_type), intent(inout) :: func
    end subroutine gauxc_functional_delete
  end interface gauxc_delete

contains

  !> @brief Create a GauXCFunctional from a string specification
  function gauxc_functional_from_string(status, functional_spec, polarized) &
    & result(func)
    type(gauxc_status_type), intent(out) :: status
    character(kind=c_char, len=*), intent(in) :: functional_spec
    logical(c_bool), intent(in) :: polarized
    type(gauxc_functional_type) :: func

    character(kind=c_char), allocatable :: c_functional_spec(:)

    c_functional_spec = transfer(functional_spec // c_null_char, &
      [character(kind=c_char) ::], len(functional_spec)+1)

    func = gauxc_functional_from_string_c(status, c_functional_spec, polarized)
  end function gauxc_functional_from_string
end module gauxc_xc_functional