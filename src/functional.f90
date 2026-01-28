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
    !> - P. A. M. Dirac., Math. Proc. Cambridge Philos. Soc. 26, 376 (1930) (doi: 10.1017/S0305004100016108)
    !> - F. Bloch., Z. Phys. 57, 545 (1929) (doi: 10.1007/BF01340281)
    !> - S. H. Vosko, L. Wilk, and M. Nusair., Can. J. Phys. 58, 1200 (1980) (doi: 10.1139/p80-159)
    !> - libxc names: LDA_X (id=1) and LDA_C_VWN_3 (id=30)
    enumerator :: gauxc_functional_svwn3

    !> @brief Slater exchange & Vosko, Wilk & Nusair correlation (VWN5)
    !> - P. A. M. Dirac., Math. Proc. Cambridge Philos. Soc. 26, 376 (1930) (doi: 10.1017/S0305004100016108)
    !> - F. Bloch., Z. Phys. 57, 545 (1929) (doi: 10.1007/BF01340281)
    !> - S. H. Vosko, L. Wilk, and M. Nusair., Can. J. Phys. 58, 1200 (1980) (doi: 10.1139/p80-159)
    !> - libxc names: LDA_X (id=1) and LDA_C_VWN (id=7)
    !> - xcfun names: SLATERX and VWN5C
    enumerator :: gauxc_functional_svwn5

    !> @brief Becke 88 exchange & Lee, Yang & Parr correlation
    !> - A. D. Becke., Phys. Rev. A 38, 3098 (1988) (doi: 10.1103/PhysRevA.38.3098)
    !> - C. Lee, W. Yang, and R. G. Parr., Phys. Rev. B 37, 785 (1988) (doi: 10.1103/PhysRevB.37.785)
    !> - B. Miehlich, A. Savin, H. Stoll, and H. Preuss., Chem. Phys. Lett. 157, 200 (1989) (doi: 10.1016/0009-2614(89)87234-3)
    !> - libxc names: GGA_X_B88 (id=106) and GGA_C_LYP (id=131)
    !> - xcfun names: BECKEX and LYPC
    enumerator :: gauxc_functional_blyp

    !> @brief Becke 88 exchange & Lee, Yang & Parr correlation, 3-parameter hybrid
    !> - P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and M. J. Frisch., J. Phys. Chem. 98, 11623 (1994) (doi: 10.1021/j100096a001)
    !> - libxc name: HYB_GGA_XC_B3LYP (id=402)
    enumerator :: gauxc_functional_b3lyp

    !> @brief Perdew-Burke-Ernzerhof exchange & correlation
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
    !> - libxc names: GGA_X_PBE (id=101) & GGA_C_PBE (id=130)
    !> - xcfun names: PBEEX and PBEC
    enumerator :: gauxc_functional_pbe

    !> @brief revised Perdew-Burke-Ernzerhof exchange & original PBE correlation
    !> - Y. Zhang and W. Yang., Phys. Rev. Lett. 80, 890 (1998) (doi: 10.1103/PhysRevLett.80.890)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
    !> - libxc names: GGA_X_PBE_R (id=102) & GGA_C_PBE (id=130)
    enumerator :: gauxc_functional_revpbe

    !> @brief Perdew-Burke-Ernzerhof exchange & correlation, 1-parameter hybrid
    !> - C. Adamo and V. Barone., J. Chem. Phys. 110, 6158 (1999) (doi: 10.1063/1.478522)
    !> - M. Ernzerhof and G. E. Scuseria., J. Chem. Phys. 110, 5029 (1999) (doi: 10.1063/1.478401)
    !> - libxc name: HYB_GGA_XC_PBEH (id=406)
    enumerator :: gauxc_functional_pbe0

    !> @brief Strongly constrained and appropriately normed (SCAN) meta-GGA
    !> - J. Sun, A. Ruzsinszky, and J. P. Perdew., Phys. Rev. Lett. 115, 036402 (2015) (doi: 10.1103/PhysRevLett.115.036402)
    !> - libxc names: MGGA_X_SCAN (id=263) & MGGA_C_SCAN (id=267)
    enumerator :: gauxc_functional_scan

    !> @brief Regularized & restored strongly constrained and appropriately normed (R2SCAN) meta-GGA
    !> - J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, and J. Sun., J. Phys. Chem. Lett. 11, 8208-8215 (2020) (doi: 10.1021/acs.jpclett.0c02405)
    !> - J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, and J. Sun., J. Phys. Chem. Lett. 11, 9248-9248 (2020) (doi: 10.1021/acs.jpclett.0c03077)
    !> - libxc names: MGGA_X_R2SCAN (id=497) & MGGA_C_R2SCAN (id=498) 
    enumerator :: gauxc_functional_r2scan

    !> @brief Regularized & restored strongly constrained and appropriately normed (R2SCAN) meta-GGA, deorbitalized version
    !> - D. Mejía-Rodríguez and S. B. Trickey., Phys. Rev. B 102, 121109 (2020) (doi: 10.1103/PhysRevB.102.121109)
    !> - J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, and J. Sun., J. Phys. Chem. Lett. 11, 8208-8215 (2020) (doi: 10.1021/acs.jpclett.0c02405)
    !> - J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, and J. Sun., J. Phys. Chem. Lett. 11, 9248-9248 (2020) (doi: 10.1021/acs.jpclett.0c03077)
    !> - libxc names: MGGA_X_R2SCANL (id=718) & MGGA_C_R2SCANL (id=719)
    enumerator :: gauxc_functional_r2scanl

    !> @brief Minnesota 2006 hybrid functional
    !> - Y. Zhao and D. G. Truhlar., Theor. Chem. Acc. 120, 215 (2008) (doi: 10.1007/s00214-007-0310-x)
    !> - libxc names: HYB_MGGA_X_M06_2X (id=450) & MGGA_C_M06_2X (id=236)
    !> - xcfun names: M062X and M062C
    enumerator :: gauxc_functional_m062x

    !> @brief Perdew, Kurth, Zupan, and Blaha
    !> - J. P. Perdew, S. Kurth, A. Zupan, and P. Blaha., Phys. Rev. Lett. 82, 2544 (1999) (doi: 10.1103/PhysRevLett.82.2544)
    !> - libxc names: MGGA_X_PKZB (id=213) & MGGA_C_PKZB (id=239)
    enumerator :: gauxc_functional_pkzb

    !> @brief epc17(-1): electron-proton correlation 2017
    !> - Y. Yang, K. R. Brorsen, T. Culpitt, M. V. Pak, and S. Hammes-Schiffer., J. Chem. Phys. 147, 114113 (2017) (doi: 10.1063/1.4996038)
    !> - libxc name: LDA_C_EPC17 (id=328)
    enumerator :: gauxc_functional_epc17_1

    !> @brief epc17-2: electron-proton correlation 2017 for proton affinities
    !> - K. R. Brorsen, Y. Yang, and S. Hammes-Schiffer., J. Phys. Chem. Lett. 8, 3488-3493 (2017) (doi: 10.1021/acs.jpclett.7b01442)
    !> - libxc name: LDA_C_EPC17_2 (id=329)
    enumerator :: gauxc_functional_epc17_2

    !> @brief epc18-1: electron-proton correlation 2018
    !> - K. R. Brorsen, P. E. Schneider, and S. Hammes-Schiffer., J. Chem. Phys. 149, 044110 (2018) (doi: 10.1063/1.5037945)
    !> - libxc name: LDA_C_EPC18_1 (id=330)
    enumerator :: gauxc_functional_epc18_1

    !> @brief epc18-2: electron-proton correlation 2018 for proton affinities
    !> - K. R. Brorsen, P. E. Schneider, and S. Hammes-Schiffer., J. Chem. Phys. 149, 044110 (2018) (doi: 10.1063/1.5037945)
    !> - libxc name: LDA_C_EPC18_2 (id=331)
    enumerator :: gauxc_functional_epc18_2

    !> @brief Grimme's parametrization of the B97 functional, original D2 variant
    !> - S. Grimme., J. Comput. Chem. 27, 1787 (2006) (doi: 10.1002/jcc.20495)
    !> - libxc name: GGA_XC_B97_D (id=170)
    enumerator :: gauxc_functional_b97d

    !> @brief Grimme's parametrization of the B97 functional, D3(0) variant
    !> - S. Grimme., J. Comput. Chem. 27, 1787 (2006) (doi: 10.1002/jcc.20495)
    !> - libxc name: GGA_XC_B97_D (id=170)
    enumerator :: gauxc_functional_b97d3zero

    !> @brief Coulomb-attenuating method range-separated hybrid functional
    !> - T. Yanai, D. P. Tew, and N. C. Handy., Chem. Phys. Lett. 393, 51 (2004) (doi: 10.1016/j.cplett.2004.06.011)
    !> - libxc name: HYB_GGA_XC_CAM_B3LYP (id=433)
    enumerator :: gauxc_functional_camb3lyp

    !> @brief Slater exchange
    !> - P. A. M. Dirac., Math. Proc. Cambridge Philos. Soc. 26, 376 (1930) (doi: 10.1017/S0305004100016108)
    !> - F. Bloch., Z. Phys. 57, 545 (1929) (doi: 10.1007/BF01340281)
    !> - libxc name: LDA_X (id=1)
    !> - xcfun name: SLATERX
    enumerator :: gauxc_functional_lda

    !> @brief Minnesota 2006 meta-GGA functional
    !> - Y. Zhao and D. G. Truhlar., J. Chem. Phys. 125, 194101 (2006) (doi: 10.1063/1.2370993)
    !> - Y. Zhao and D. G. Truhlar., Theor. Chem. Acc. 120, 215 (2008) (doi: 10.1007/s00214-007-0310-x)
    !> - libxc names: MGGA_X_M06_L (id=203) & MGGA_C_M06_L (id=233)
    enumerator :: gauxc_functional_m06l

    !> @brief Strongly constrained and appropriately normed (SCAN) meta-GGA, 1-parameter hybrid
    !> - K. Hui and J.-D. Chai., J. Chem. Phys. 144, 044114 (2016) (doi: 10.1063/1.4940734)
    !> - J. Sun, A. Ruzsinszky, and J. P. Perdew., Phys. Rev. Lett. 115, 036402 (2015) (doi: 10.1103/PhysRevLett.115.036402)
    !> - libxc names: HYB_MGGA_X_SCAN0 (id=264) & MGGA_C_SCAN (id=267)
    enumerator :: gauxc_functional_scan0

    !> @brief Slater exchange & Perdew, Wang 92 correlation
    !> - P. A. M. Dirac., Math. Proc. Cambridge Philos. Soc. 26, 376 (1930) (doi: 10.1017/S0305004100016108)
    !> - F. Bloch., Z. Phys. 57, 545 (1929) (doi: 10.1007/BF01340281)
    !> - J. P. Perdew and Y. Wang., Phys. Rev. B 45, 13244 (1992) (doi: 10.1103/PhysRevB.45.13244)
    !> - libxc names: LDA_X (id=1) and LDA_C_PW (id=12)
    enumerator :: gauxc_functional_spw92

    !> @brief Tao, Perdew, Staroverov & Scuseria meta-GGA
    !> - J. Tao, J. P. Perdew, V. N. Staroverov, and G. E. Scuseria., Phys. Rev. Lett. 91, 146401 (2003) (doi: 10.1103/PhysRevLett.91.146401)
    !> - J. P. Perdew, J. Tao, V. N. Staroverov, and G. E. Scuseria., J. Chem. Phys. 120, 6898 (2004) (doi: 10.1063/1.1665298)
    !> - libxc names: MGGA_X_TPSS (id=202) and MGGA_C_TPSS (id=231)
    enumerator :: gauxc_functional_tpss

    !> @brief Tao, Perdew, Staroverov & Scuseria meta-GGA, 1-parameter hybrid
    !> - V. N. Staroverov, G. E. Scuseria, J. Tao, and J. P. Perdew., J. Chem. Phys. 119, 12129 (2003) (doi: 10.1063/1.1626543)
    !> - libxc name: HYB_MGGA_XC_TPSSH (id=457)
    enumerator :: gauxc_functional_tpssh

    !> @brief Tao, Perdew, Staroverov & Scuseria meta-GGA, 1-parameter hybrid
    !> - S. Grimme., J. Phys. Chem. A 109, 3067-3077 (2005) (doi: 10.1021/jp050036j)
    !> - libxc name: HYB_MGGA_XC_TPSS0 (id=396)
    enumerator :: gauxc_functional_tpss0

    !> @brief Vosko, Wilk & Nusair correlation (VWN3)
    !> - S. H. Vosko, L. Wilk, and M. Nusair., Can. J. Phys. 58, 1200 (1980) (doi: 10.1139/p80-159)
    !> - libxc name: LDA_C_VWN_3 (id=30)
    enumerator :: gauxc_functional_vwn3

    !> @brief Vosko, Wilk & Nusair correlation (VWN5)
    !> - S. H. Vosko, L. Wilk, and M. Nusair., Can. J. Phys. 58, 1200 (1980) (doi: 10.1139/p80-159)
    !> - libxc name: LDA_C_VWN (id=7)
    enumerator :: gauxc_functional_vwn5

    !> @brief HJS screened PBE exchange & original PBE correlation
    !> - T. M. Henderson, B. G. Janesko, and G. E. Scuseria., J. Chem. Phys. 128, 194105 (2008) (doi: 10.1063/1.2921797)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
    !> - libxc names: GGA_X_HJS_PBE (id=525) & GGA_C_PBE (id=130)
    enumerator :: gauxc_functional_lrcwpbe

    !> @brief HJS screened PBE exchange & original PBE correlation, hybrid version
    !> - T. M. Henderson, B. G. Janesko, and G. E. Scuseria., J. Chem. Phys. 128, 194105 (2008) (doi: 10.1063/1.2921797)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
    !> - libxc name: HYB_GGA_XC_HJS_PBE (id=429)
    enumerator :: gauxc_functional_lrcwpbeh

    !> @brief Becke 88 exchange and Perdew 86 correlation
    !> - A. D. Becke., Phys. Rev. A 38, 3098 (1988) (doi: 10.1103/PhysRevA.38.3098)
    !> - J. P. Perdew., Phys. Rev. B 33, 8822 (1986) (doi: 10.1103/PhysRevB.33.8822)
    !> - libxc names: GGA_X_B88 (id=106) and GGA_C_P86 (id=132)
    enumerator :: gauxc_functional_bp86

    !> @brief Heyd-Scuseria-Ernzerhof screened hybrid functional (HSE03)
    !> - J. Heyd, G. E. Scuseria, and M. Ernzerhof., J. Chem. Phys. 118, 8207 (2003) (doi: 10.1063/1.1564060)
    !> - J. Heyd, G. E. Scuseria, and M. Ernzerhof., J. Chem. Phys. 124, 219906 (2006) (doi: 10.1063/1.2204597)
    !> - libxc name: HYB_GGA_XC_HSE03 (id=427)
    enumerator :: gauxc_functional_hse03

    !> @brief Heyd-Scuseria-Ernzerhof screened hybrid functional (HSE06)
    !> - J. Heyd, G. E. Scuseria, and M. Ernzerhof., J. Chem. Phys. 118, 8207 (2003) (doi: 10.1063/1.1564060)
    !> - J. Heyd, G. E. Scuseria, and M. Ernzerhof., J. Chem. Phys. 124, 219906 (2006) (doi: 10.1063/1.2204597)
    !> - A. V. Krukau, O. A. Vydrov, A. F. Izmaylov, and G. E. Scuseria., J. Chem. Phys. 125, 224106 (2006) (doi: 10.1063/1.2404663)
    !> - libxc name: HYB_GGA_XC_HSE06 (id=428)
    enumerator :: gauxc_functional_hse06

    !> @brief Revised B3LYP
    !> - L. Lu, H. Hu, H. Hou, and B. Wang., Comput. Theor. Chem. 1015, 64 (2013) (doi: 10.1016/j.comptc.2013.04.009)
    !> - libxc name: HYB_GGA_XC_REVB3LYP (id=454)
    enumerator :: gauxc_functional_revb3lyp

    !> @brief revised Perdew-Burke-Ernzerhof exchange & original PBE correlation, hybrid version
    !> - Y. Zhang and W. Yang., Phys. Rev. Lett. 80, 890 (1998) (doi: 10.1103/PhysRevLett.80.890)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
    !> - libxc names: GGA_X_PBE_R (id=102) & GGA_C_PBE (id=130)
    enumerator :: gauxc_functional_revpbe0

    !> @brief revised Tao, Perdew, Staroverov & Scuseria
    !> - J. P. Perdew, A. Ruzsinszky, G. I. Csonka, L. A. Constantin, and J. Sun., Phys. Rev. Lett. 103, 026403 (2009) (doi: 10.1103/PhysRevLett.103.026403)
    !> - J. P. Perdew, A. Ruzsinszky, G. I. Csonka, L. A. Constantin, and J. Sun., Phys. Rev. Lett. 106, 179902 (2011) (doi: 10.1103/PhysRevLett.106.179902)
    !> - libxc names: MGGA_X_REVTPSS (id=212) & MGGA_C_REVTPSS (id=241)
    enumerator :: gauxc_functional_revtpss

    !> @brief revTPSSh
    !> - G. I. Csonka, J. P. Perdew, and A. Ruzsinszky., J. Chem. Theory Comput. 6, 3688 (2010) (doi: 10.1021/ct100488v)
    !> - libxc name: HYB_MGGA_XC_REVTPSSH (id=458)
    enumerator :: gauxc_functional_revtpssh

    !> @brief Perdew-Wang 91 exchange and correlation
    !> - J. P. Perdew. In P. Ziesche and H. Eschrig, editors, Proceedings of the 75. WE-Heraeus-Seminar and 21st Annual International Symposium on Electronic Structure of Solids, 11. Berlin, 1991. Akademie Verlag.
    !> - J. P. Perdew, J. A. Chevary, S. H. Vosko, K. A. Jackson, M. R. Pederson, D. J. Singh, and C. Fiolhais., Phys. Rev. B 46, 6671 (1992) (doi: 10.1103/PhysRevB.46.6671)
    !> - J. P. Perdew, J. A. Chevary, S. H. Vosko, K. A. Jackson, M. R. Pederson, D. J. Singh, and C. Fiolhais., Phys. Rev. B 48, 4978 (1993) (doi: 10.1103/PhysRevB.48.4978.2)
    !> - libxc name: GGA_X_PW91 (id=109) and GGA_C_PW91 (id=134)
    enumerator :: gauxc_functional_pw91

    !> @brief mBEEF exchange and Perdew, Burke & Ernzerhof SOL
    !> - J. Wellendorff, K. T. Lundgaard, K. W. Jacobsen, and T. Bligaard., J. Chem. Phys. 140, 144107 (2014) (doi: 10.1063/1.4870397)
    !> - J. P. Perdew, A. Ruzsinszky, G. I. Csonka, O. A. Vydrov, G. E. Scuseria, L. A. Constantin, X. Zhou, and K. Burke., Phys. Rev. Lett. 100, 136406 (2008) (doi: 10.1103/PhysRevLett.100.136406)
    !> - libxc names: MGGA_X_MBEEF (id=249) and GGA_C_PBE_SOL (id=133)
    enumerator :: gauxc_functional_mbeef

    !> @brief The original (ACM, B3PW91) hybrid of Becke
    !> - A. D. Becke., J. Chem. Phys. 98, 5648 (1993) (doi: 10.1063/1.464913)
    !> - libxc name: HYB_GGA_XC_B3PW91 (id=401)
    enumerator :: gauxc_functional_b3pw91

    !> @brief O3LYP
    !> - W.-M. Hoe, A. J. Cohen, and N. C. Handy., Chem. Phys. Lett. 341, 319–328 (2001) (doi: 10.1016/S0009-2614(01)00581-4)
    !> - A. J. Cohen and N. C. Handy., Mol. Phys. 99, 607 (2001) (doi: 10.1080/00268970010023435)
    !> - libxc name: HYB_GGA_XC_O3LYP (id=404)
    enumerator :: gauxc_functional_o3lyp

    !> @brief Handy & Cohen OPTX 01 exchange and Lee, Yang & Parr correlation
    !> - N. C. Handy and A. J. Cohen., Mol. Phys. 99, 403 (2001) (doi: 10.1080/00268970010018431)
    !> - C. Lee, W. Yang, and R. G. Parr., Phys. Rev. B 37, 785 (1988) (doi: 10.1103/PhysRevB.37.785)
    !> - B. Miehlich, A. Savin, H. Stoll, and H. Preuss., Chem. Phys. Lett. 157, 200 (1989) (doi: 10.1016/0009-2614(89)87234-3)
    !> - libxc names: GGA_X_OPTX (id=110) & GGA_C_LYP (id=131)
    enumerator :: gauxc_functional_olyp

    !> @brief Handy & Cohen OPTX 01 exchange and Perdew, Burke & Ernzerhof correlation
    !> - N. C. Handy and A. J. Cohen., Mol. Phys. 99, 403 (2001) (doi: 10.1080/00268970010018431)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
    !> - libxc names: GGA_X_OPTX (id=110) & GGA_C_PBE (id=130)
    enumerator :: gauxc_functional_opbe

    !> @brief mPW1K
    !> - B. J. Lynch, P. L. Fast, M. Harris, and D. G. Truhlar., J. Phys. Chem. A 104, 4811 (2000) (doi: 10.1021/jp000497z)
    !> - libxc name: HYB_GGA_XC_MPW1K (id=405)
    enumerator :: gauxc_functional_mpw1k

    !> @brief Revised Perdew-Burke-Ernzerhof exchange by Hammer, Hansen, and Norskov
    !> - B. Hammer, L. B. Hansen, and J. K. Nørskov., Phys. Rev. B 59, 7413 (1999) (doi: 10.1103/PhysRevB.59.7413)
    !> - libxc name: GGA_X_RPBE (id=117)
    enumerator :: gauxc_functional_rpbe

    !> @brief Becke 88 exchange
    !> - A. D. Becke., Phys. Rev. A 38, 3098 (1988) (doi: 10.1103/PhysRevA.38.3098)
    !> - libxc name: GGA_X_B88 (id=106)
    enumerator :: gauxc_functional_b88

    !> @brief modified Perdew-Wang 91 exchange by Adamo & Barone
    !> - C. Adamo and V. Barone., J. Chem. Phys. 108, 664 (1998) (doi: 10.1063/1.475428)
    !> - libxc name: GGA_X_MPW91 (id=119)
    enumerator :: gauxc_functional_mpw91

    !> @brief Regularized strongly constrained and appropriately normed (RSCAN) meta-GGA by Bartok and Yates
    !> - A. P. Bartók and J. R. Yates., J. Chem. Phys. 150, 161101 (2019) (doi: 10.1063/1.5094646)
    !> - libxc names: MGGA_X_RSCAN (id=493) and MGGA_C_RSCAN (id=494)
    enumerator :: gauxc_functional_rscan

    !> @brief CAM version of B3LYP, tuned for excitations and properties
    !> - K. Okuno, Y. Shigeta, R. Kishi, H. Miyasaka, and M. Nakano., J. Photochem. Photobiol., A 235, 29 (2012) (doi: 10.1016/j.jphotochem.2012.03.003)
    !> - libxc name: HYB_GGA_XC_TUNED_CAM_B3LYP (id=434)
    enumerator :: gauxc_functional_tunedcamb3lyp

    !> @brief wB97 range-separated functional
    !> - J.-D. Chai and M. Head-Gordon., J. Chem. Phys. 128, 084106 (2008) (doi: 10.1063/1.2834918)
    !> - libxc name: HYB_GGA_XC_WB97 (id=463)
    enumerator :: gauxc_functional_wb97

    !> @brief wB97X range-separated functional
    !> - J.-D. Chai and M. Head-Gordon., J. Chem. Phys. 128, 084106 (2008) (doi: 10.1063/1.2834918)
    !> - libxc name: HYB_GGA_XC_WB97X (id=464)
    enumerator :: gauxc_functional_wb97x

    !> @brief wB97X-D range-separated functional
    !> - J.-D. Chai and M. Head-Gordon., Phys. Chem. Chem. Phys. 10, 6615-6620 (2008) (doi: 10.1039/B810189B)
    !> - libxc name: HYB_GGA_XC_WB97X_D (id=471)
    enumerator :: gauxc_functional_wb97xd

    !> @brief wB97X-D3 range-separated functional
    !> - Y.-S. Lin, G.-D. Li, S.-P. Mao, and J.-D. Chai., J. Chem. Theory Comput. 9, 263-272 (2013) (doi: 10.1021/ct300715s)
    !> - libxc name: HYB_GGA_XC_WB97X_D3 (id=399)
    enumerator :: gauxc_functional_wb97xd3

    !> @brief Long-range corrected PBE (LC-wPBE) by Vydrov and Scuseria
    !> - O. A. Vydrov and G. E. Scuseria., J. Chem. Phys. 125, 234109 (2006) (doi: 10.1063/1.2409292)
    !> - libxc name: HYB_GGA_XC_LC_WPBE (id=478)
    enumerator :: gauxc_functional_lcwpbe

    !> @brief X3LYP
    !> - X. Xu and W. A. Goddard., Proc. Natl. Acad. Sci. U. S. A. 101, 2673 (2004) (doi: 10.1073/pnas.0308730100)
    !> - libxc name: HYB_GGA_XC_X3LYP (id=411)
    enumerator :: gauxc_functional_x3lyp

    !> @brief XLYP
    !> - X. Xu and W. A. Goddard., Proc. Natl. Acad. Sci. U. S. A. 101, 2673 (2004) (doi: 10.1073/pnas.0308730100)
    !> - libxc name: GGA_XC_XLYP (id=166)
    enumerator :: gauxc_functional_xlyp

    !> @brief BHandH i.e. BHLYP
    !> - A. D. Becke., J. Chem. Phys. 98, 1372 (1993) (doi: 10.1063/1.464304)
    !> - Defined through Gaussian implementation.
    !> - libxc name: HYB_GGA_XC_BHANDH (id=435)
    enumerator :: gauxc_functional_bhandh

    !> @brief Boese-Martin for kinetics
    !> - A. D. Boese and J. M. L. Martin., J. Chem. Phys. 121, 3405 (2004) (doi: 10.1063/1.1774975)
    !> - libxc names: HYB_MGGA_X_BMK (id=279) & GGA_C_BMK (id=280)
    enumerator :: gauxc_functional_bmk

    !> @brief Becke 88 exchange and Perdew 86 based on VWN5 correlation, with more accurate value for ftilde
    !> - A. D. Becke., Phys. Rev. A 38, 3098 (1988) (doi: 10.1103/PhysRevA.38.3098)
    !> - J. P. Perdew., Phys. Rev. B 33, 8822 (1986) (doi: 10.1103/PhysRevB.33.8822)
    !> - libxc names: GGA_X_B88 (id=106) & GGA_C_P86VWN_FT (id=253)
    enumerator :: gauxc_functional_bp86vwn

    !> @brief Mixture of PW86 with BC95
    !> - A. D. Becke., J. Chem. Phys. 104, 1040 (1996) (doi: 10.1063/1.470829)
    !> - libxc name: HYB_MGGA_XC_PW86B95 (id=442)
    enumerator :: gauxc_functional_pw86b95

    !> @brief Perdew & Wang 86 exchange and PBE correlation
    !> - J. P. Perdew and W. Yue., Phys. Rev. B 33, 8800 (1986) (doi: 10.1103/PhysRevB.33.8800)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
    !> - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
    !> - libxc names: GGA_X_PW86 (id=108) & GGA_C_PBE (id=130)
    enumerator :: gauxc_functional_pw86pbe

    !> @brief r2SCAN0: r2SCAN hybrid like PBE0 with 25% exact exchange
    !> - M. Bursch, H. Neugebauer, S. Ehlert, and S. Grimme., J. Chem. Phys. 156, 134105 (2022) (doi: 10.1063/5.0086040)
    !> - libxc name: HYB_MGGA_XC_R2SCAN0 (id=660)
    enumerator :: gauxc_functional_r2scan0

    !> @brief r2SCANh: r2SCAN hybrid like TPSSh with 10% exact exchange
    !> - M. Bursch, H. Neugebauer, S. Ehlert, and S. Grimme., J. Chem. Phys. 156, 134105 (2022) (doi: 10.1063/5.0086040)
    !> - libxc name: HYB_MGGA_XC_R2SCANH (id=659)
    enumerator :: gauxc_functional_r2scanh

    !> @brief r2SCAN50: r2SCAN hybrid like BHLYP with 50% exact exchange
    !> - M. Bursch, H. Neugebauer, S. Ehlert, and S. Grimme., J. Chem. Phys. 156, 134105 (2022) (doi: 10.1063/5.0086040)
    !> - libxc name: HYB_MGGA_XC_R2SCAN50 (id=661)
    enumerator :: gauxc_functional_r2scan50

    !> @brief Minnesota 2005 hybrid functional
    !> - Y. Zhao, N. E. Schultz, and D. G. Truhlar., J. Chem. Phys. 123, 161103 (2005) (doi: 10.1063/1.2126975)
    !> - libxc names: HYB_MGGA_X_M05 (id=438) & MGGA_C_M05 (id=237)
    enumerator :: gauxc_functional_m05

    !> @brief Minnesota 2008 hybrid functional
    !> - Y. Zhao and D. G. Truhlar., Theor. Chem. Acc. 120, 215 (2008) (doi: 10.1007/s00214-007-0310-x)
    !> - libxc names: HYB_MGGA_X_M06 (id=449) & MGGA_C_M06 (id=235)
    enumerator :: gauxc_functional_m06

    !> @brief Minnesota M08 hybrid functional
    !> - Y. Zhao and D. G. Truhlar., J. Chem. Theory Comput. 4, 1849 (2008) (doi: 10.1021/ct800246v)
    !> - libxc names: HYB_MGGA_X_M08_HX (id=295) & MGGA_C_M08_HX (id=78)
    enumerator :: gauxc_functional_m08hx

    !> @brief Minnesota M08-SO hybrid exchange functional
    !> - Y. Zhao and D. G. Truhlar., J. Chem. Theory Comput. 4, 1849 (2008) (doi: 10.1021/ct800246v)
    !> - libxc names: HYB_MGGA_X_M08_SO (id=296) & MGGA_C_M08_SO (id=77)
    enumerator :: gauxc_functional_m08so

    !> @brief Minnesota M05-2X hybrid exchange functional
    !> - Y. Zhao, N. E. Schultz, and D. G. Truhlar., J. Chem. Theory Comput. 2, 364 (2006) (doi: 10.1021/ct0502763)
    !> - libxc names: HYB_MGGA_X_M05_2X (id=439) & MGGA_C_M05_2X (id=238)
    enumerator :: gauxc_functional_m052x

    !> @brief Minnesota M06-SX short-range hybrid exchange functional
    !> - Y. Wang, P. Verma, L. Zhang, Y. Li, Z. Liu, D. G. Truhlar, and X. He., Proc. Natl. Acad. Sci. U. S. A. 117, 2294–2301 (2020) (doi: 10.1073/pnas.1913699117)
    !> - libxc names: HYB_MGGA_X_M06_SX (id=310) & MGGA_C_M06_SX (id=311)
    enumerator :: gauxc_functional_m06sx

    !> @brief Minnesota CF22D hybrid exchange functional
    !> - Y. Liu, C. Zhang, Z. Liu, D. G. Truhlar, Y. Wang, and X. He., Nature Computational Science 3, 48–58 (2022) (doi: 10.1038/s43588-022-00371-5)
    !> - libxc names: HYB_MGGA_X_CF22D (id=340) & MGGA_C_CF22D (id=341)
    enumerator :: gauxc_functional_cf22d

    !> @brief Hybrid based on SOGGA11 form
    !> - R. Peverati and D. G. Truhlar., J. Chem. Phys. 135, 191102 (2011) (doi: 10.1063/1.3663871)
    !> - libxc names: HYB_GGA_X_SOGGA11_X (id=426) & GGA_C_SOGGA11_X (id=159)
    enumerator :: gauxc_functional_sogga11x

    !> @brief Minnesota M06-HF hybrid exchange functional
    !> - Y. Zhao and D. G. Truhlar., J. Phys. Chem. A 110, 13126 (2006) (doi: 10.1021/jp066479k)
    !> - libxc names: HYB_MGGA_X_M06_HF (id=444) & MGGA_C_M06_HF (id=234)
    enumerator :: gauxc_functional_m06hf

    !> @brief Minnesota M11 hybrid exchange functional
    !> - R. Peverati and D. G. Truhlar., J. Phys. Chem. Lett. 2, 2810 (2011) (doi: 10.1021/jz201170d)
    !> - libxc names: HYB_MGGA_X_M11 (id=297) & MGGA_C_M11 (id=76)
    enumerator :: gauxc_functional_m11

    !> @brief Minnesota MN12-L exchange functional
    !> - R. Peverati and D. G. Truhlar., Phys. Chem. Chem. Phys. 14, 13171 (2012) (doi: 10.1039/C2CP42025B)
    !> - libxc names: MGGA_X_MN12_L (id=227) & MGGA_C_MN12_L (id=74)
    enumerator :: gauxc_functional_mn12l

    !> @brief Minnesota MN12-SX hybrid exchange functional
    !> - R. Peverati and D. G. Truhlar., Phys. Chem. Chem. Phys. 14, 16187 (2012) (doi: 10.1039/C2CP42576A)
    !> - libxc names: HYB_MGGA_X_MN12_SX (id=248) & MGGA_C_MN12_SX (id=73)
    enumerator :: gauxc_functional_mn12sx

    !> @brief Minnesota MN15 correlation functional
    !> - H. S. Yu, X. He, S. L. Li, and D. G. Truhlar., Chem. Sci. 7, 5032-5051 (2016) (doi: 10.1039/C6SC00705H)
    !> - libxc names: HYB_MGGA_X_MN15 (id=268) & MGGA_C_MN15 (id=269)
    enumerator :: gauxc_functional_mn15

    !> @brief Minnesota MN15-L exchange functional
    !> - H. S. Yu, X. He, and D. G. Truhlar., J. Chem. Theory Comput. 12, 1280-1293 (2016) (doi: 10.1021/acs.jctc.5b01082)
    !> - libxc names: MGGA_X_MN15_L (id=260) & MGGA_C_MN15_L (id=261)
    enumerator :: gauxc_functional_mn15l

    !> @brief Revised Minnesota 2006 meta-GGA functional
    !> - Y. Wang, X. Jin, H. S. Yu, D. G. Truhlar, and X. He., Proc. Natl. Acad. Sci. U. S. A. 114, 8487-8492 (2017) (doi: 10.1073/pnas.1705670114)
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
      logical(c_bool), intent(in) :: polarized
      type(gauxc_functional_type) :: func
    end function gauxc_functional_from_string_c

    !> @brief Create a GauXCFunctional from an enum specification
    function gauxc_functional_from_enum(status, functional_enum, polarized) &
      & result(func) bind(c)
      import :: gauxc_status_type, gauxc_functional_type, c_int, c_bool
      type(gauxc_status_type), intent(out) :: status
      integer(c_int), intent(in) :: functional_enum
      logical(c_bool), intent(in) :: polarized
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