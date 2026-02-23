/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include <gauxc/types.h>
#include <gauxc/status.h>

#ifdef __cplusplus
extern "C" {
namespace GauXC::C {
#endif

enum GauXC_Functional {
  /// @brief Slater exchange & Vosko, Wilk & Nusair correlation (VWN3)
  /// - libxc names: LDA_X (id=1) and LDA_C_VWN_3 (id=30)
  GauXC_Functional_SVWN3,

  /// @brief Slater exchange & Vosko, Wilk & Nusair correlation (VWN5)
  /// - libxc names: LDA_X (id=1) and LDA_C_VWN (id=7)
  /// - xcfun names: SLATERX and VWN5C
  GauXC_Functional_SVWN5,

  /// @brief Becke 88 exchange & Lee, Yang & Parr correlation
  /// - libxc names: GGA_X_B88 (id=106) and GGA_C_LYP (id=131)
  /// - xcfun names: BECKEX and LYPC
  GauXC_Functional_BLYP,

  /// @brief Becke 88 exchange & Lee, Yang & Parr correlation, 3-parameter hybrid
  /// - libxc name: HYB_GGA_XC_B3LYP (id=402)
  GauXC_Functional_B3LYP,

  /// @brief Perdew-Burke-Ernzerhof exchange & correlation
  /// - libxc names: GGA_X_PBE (id=101) & GGA_C_PBE (id=130)
  /// - xcfun names: PBEEX and PBEC
  GauXC_Functional_PBE,

  /// @brief revised Perdew-Burke-Ernzerhof exchange & original PBE correlation
  /// - libxc names: GGA_X_PBE_R (id=102) & GGA_C_PBE (id=130)
  GauXC_Functional_revPBE,

  /// @brief Perdew-Burke-Ernzerhof exchange & correlation, 1-parameter hybrid
  /// - libxc name: HYB_GGA_XC_PBEH (id=406)
  GauXC_Functional_PBE0,

  /// @brief Strongly constrained and appropriately normed (SCAN) meta-GGA
  /// - libxc names: MGGA_X_SCAN (id=263) & MGGA_C_SCAN (id=267)
  GauXC_Functional_SCAN,

  /// @brief Regularized & restored strongly constrained and appropriately normed (R2SCAN) meta-GGA
  /// - libxc names: MGGA_X_R2SCAN (id=497) & MGGA_C_R2SCAN (id=498) 
  GauXC_Functional_R2SCAN,

  /// @brief Regularized & restored strongly constrained and appropriately normed (R2SCAN) meta-GGA, deorbitalized version
  /// - libxc names: MGGA_X_R2SCANL (id=718) & MGGA_C_R2SCANL (id=719)
  GauXC_Functional_R2SCANL,

  /// @brief Minnesota 2006 hybrid functional
  /// - xcfun names: M062X and M062C
  GauXC_Functional_M062X,

  /// @brief Perdew, Kurth, Zupan, and Blaha
  /// - libxc names: MGGA_X_PKZB (id=213) & MGGA_C_PKZB (id=239)
  GauXC_Functional_PKZB,

  /// @brief epc17(-1): electron-proton correlation 2017
  /// - libxc name: LDA_C_EPC17 (id=328)
  GauXC_Functional_EPC17_1,

  /// @brief epc17-2: electron-proton correlation 2017 for proton affinities
  /// - libxc name: LDA_C_EPC17_2 (id=329)
  GauXC_Functional_EPC17_2,

  /// @brief epc18-1: electron-proton correlation 2018
  /// - libxc name: LDA_C_EPC18_1 (id=330)
  GauXC_Functional_EPC18_1,

  /// @brief epc18-2: electron-proton correlation 2018 for proton affinities
  /// - libxc name: LDA_C_EPC18_2 (id=331)
  GauXC_Functional_EPC18_2,

  /// @brief Grimme's parametrization of the B97 functional, original D2 variant
  /// - libxc name: GGA_XC_B97_D (id=170)
  GauXC_Functional_B97D,

  /// @brief Grimme's parametrization of the B97 functional, D3(0) variant
  /// - libxc name: GGA_XC_B97_D (id=170)
  GauXC_Functional_B97D3ZERO,

  /// @brief Coulomb-attenuating method range-separated hybrid functional
  /// - libxc name: HYB_GGA_XC_CAM_B3LYP (id=433)
  GauXC_Functional_CAMB3LYP,

  /// @brief Slater exchange
  /// - libxc name: LDA_X (id=1)
  /// - xcfun name: SLATERX
  GauXC_Functional_LDA,

  /// @brief Minnesota 2006 meta-GGA functional
  /// - libxc names: MGGA_X_M06_L (id=203) & MGGA_C_M06_L (id=233)
  GauXC_Functional_M06L,

  /// @brief Strongly constrained and appropriately normed (SCAN) meta-GGA, 1-parameter hybrid
  /// - libxc names: HYB_MGGA_X_SCAN0 (id=264) & MGGA_C_SCAN (id=267)
  GauXC_Functional_SCAN0,

  /// @brief Slater exchange & Perdew, Wang 92 correlation
  /// - libxc names: LDA_X (id=1) and LDA_C_PW (id=12)
  GauXC_Functional_SPW92,

  /// @brief Tao, Perdew, Staroverov & Scuseria meta-GGA
  /// - libxc names: MGGA_X_TPSS (id=202) and MGGA_C_TPSS (id=231)
  GauXC_Functional_TPSS,

  /// @brief Tao, Perdew, Staroverov & Scuseria meta-GGA, 1-parameter hybrid
  /// - libxc name: HYB_MGGA_XC_TPSSH (id=457)
  GauXC_Functional_TPSSh,

  /// @brief Tao, Perdew, Staroverov & Scuseria meta-GGA, 1-parameter hybrid
  /// - libxc name: HYB_MGGA_XC_TPSS0 (id=396)
  GauXC_Functional_TPSS0,

  /// @brief Vosko, Wilk & Nusair correlation (VWN3)
  /// - libxc name: LDA_C_VWN_3 (id=30)
  GauXC_Functional_VWN3,

  /// @brief Vosko, Wilk & Nusair correlation (VWN5)
  /// - libxc name: LDA_C_VWN (id=7)
  GauXC_Functional_VWN5,

  /// @brief HJS screened PBE exchange & original PBE correlation
  /// - libxc names: GGA_X_HJS_PBE (id=525) & GGA_C_PBE (id=130)
  GauXC_Functional_LRCwPBE,

  /// @brief HJS screened PBE exchange & original PBE correlation, hybrid version
  /// - libxc name: HYB_GGA_XC_HJS_PBE (id=429)
  GauXC_Functional_LRCwPBEh,

  /// @brief Becke 88 exchange and Perdew 86 correlation
  /// - libxc names: GGA_X_B88 (id=106) and GGA_C_P86 (id=132)
  GauXC_Functional_BP86,

  /// @brief Heyd-Scuseria-Ernzerhof screened hybrid functional (HSE03)
  /// - libxc name: HYB_GGA_XC_HSE03 (id=427)
  GauXC_Functional_HSE03,

  /// @brief Heyd-Scuseria-Ernzerhof screened hybrid functional (HSE06)
  /// - libxc name: HYB_GGA_XC_HSE06 (id=428)
  GauXC_Functional_HSE06,

  /// @brief Revised B3LYP
  /// - libxc name: HYB_GGA_XC_REVB3LYP (id=454)
  GauXC_Functional_revB3LYP,

  /// @brief revised Perdew-Burke-Ernzerhof exchange & original PBE correlation, hybrid version
  /// - libxc names: GGA_X_PBE_R (id=102) & GGA_C_PBE (id=130)
  GauXC_Functional_revPBE0,

  /// @brief revised Tao, Perdew, Staroverov & Scuseria
  /// - libxc names: MGGA_X_REVTPSS (id=212) & MGGA_C_REVTPSS (id=241)
  GauXC_Functional_revTPSS,

  /// @brief revTPSSh
  /// - libxc name: HYB_MGGA_XC_REVTPSSH (id=458)
  GauXC_Functional_revTPSSh,

  /// @brief Perdew-Wang 91 exchange and correlation
  /// - libxc name: GGA_X_PW91 (id=109) and GGA_C_PW91 (id=134)
  GauXC_Functional_PW91,

  /// @brief mBEEF exchange and Perdew, Burke & Ernzerhof SOL
  /// - libxc names: MGGA_X_MBEEF (id=249) and GGA_C_PBE_SOL (id=133)
  GauXC_Functional_mBEEF,

  /// @brief The original (ACM, B3PW91) hybrid of Becke
  /// - libxc name: HYB_GGA_XC_B3PW91 (id=401)
  GauXC_Functional_B3PW91,

  /// @brief O3LYP
  /// - libxc name: HYB_GGA_XC_O3LYP (id=404)
  GauXC_Functional_O3LYP,

  /// @brief Handy & Cohen OPTX 01 exchange and Lee, Yang & Parr correlation
  /// - libxc names: GGA_X_OPTX (id=110) & GGA_C_LYP (id=131)
  GauXC_Functional_OLYP,

  /// @brief Handy & Cohen OPTX 01 exchange and Perdew, Burke & Ernzerhof correlation
  /// - libxc names: GGA_X_OPTX (id=110) & GGA_C_PBE (id=130)
  GauXC_Functional_OPBE,

  /// @brief mPW1K
  /// - libxc name: HYB_GGA_XC_MPW1K (id=405)
  GauXC_Functional_MPW1K,

  /// @brief Revised Perdew-Burke-Ernzerhof exchange by Hammer, Hansen, and Norskov
  /// - libxc name: GGA_X_RPBE (id=117)
  GauXC_Functional_RPBE,

  /// @brief Becke 88 exchange
  /// - libxc name: GGA_X_B88 (id=106)
  GauXC_Functional_B88,

  /// @brief modified Perdew-Wang 91 exchange by Adamo & Barone
  /// - libxc name: GGA_X_MPW91 (id=119)
  GauXC_Functional_MPW91,

  /// @brief Regularized strongly constrained and appropriately normed (RSCAN) meta-GGA by Bartok and Yates
  /// - libxc names: MGGA_X_RSCAN (id=493) and MGGA_C_RSCAN (id=494)
  GauXC_Functional_RSCAN,

  /// @brief CAM version of B3LYP, tuned for excitations and properties
  /// - libxc name: HYB_GGA_XC_TUNED_CAM_B3LYP (id=434)
  GauXC_Functional_TUNEDCAMB3LYP,

  /// @brief wB97 range-separated functional
  /// - libxc name: HYB_GGA_XC_WB97 (id=463)
  GauXC_Functional_wB97,

  /// @brief wB97X range-separated functional
  /// - libxc name: HYB_GGA_XC_WB97X (id=464)
  GauXC_Functional_wB97X,

  /// @brief wB97X-D range-separated functional
  /// - libxc name: HYB_GGA_XC_WB97X_D (id=471)
  GauXC_Functional_wB97XD,

  /// @brief wB97X-D3 range-separated functional
  /// - libxc name: HYB_GGA_XC_WB97X_D3 (id=399)
  GauXC_Functional_wB97XD3,

  /// @brief Long-range corrected PBE (LC-wPBE) by Vydrov and Scuseria
  /// - libxc name: HYB_GGA_XC_LC_WPBE (id=478)
  GauXC_Functional_LCwPBE,

  /// @brief X3LYP
  /// - libxc name: HYB_GGA_XC_X3LYP (id=411)
  GauXC_Functional_X3LYP,

  /// @brief XLYP
  /// - libxc name: GGA_XC_XLYP (id=166)
  GauXC_Functional_XLYP,

  /// @brief BHandH i.e. BHLYP
  /// - libxc name: HYB_GGA_XC_BHANDH (id=435)
  GauXC_Functional_BHANDH,

  /// @brief Boese-Martin for kinetics
  /// - libxc names: HYB_MGGA_X_BMK (id=279) & GGA_C_BMK (id=280)
  GauXC_Functional_BMK,

  /// @brief Becke 88 exchange and Perdew 86 based on VWN5 correlation, with more accurate value for ftilde
  /// - libxc names: GGA_X_B88 (id=106) & GGA_C_P86VWN_FT (id=253)
  GauXC_Functional_BP86VWN,

  /// @brief Mixture of PW86 with BC95
  /// - libxc name: HYB_MGGA_XC_PW86B95 (id=442)
  GauXC_Functional_PW86B95,

  /// @brief Perdew & Wang 86 exchange and PBE correlation
  /// - libxc names: GGA_X_PW86 (id=108) & GGA_C_PBE (id=130)
  GauXC_Functional_PW86PBE,

  /// @brief r2SCAN0: r2SCAN hybrid like PBE0 with 25% exact exchange
  /// - libxc name: HYB_MGGA_XC_R2SCAN0 (id=660)
  GauXC_Functional_R2SCAN0,

  /// @brief r2SCANh: r2SCAN hybrid like TPSSh with 10% exact exchange
  /// - libxc name: HYB_MGGA_XC_R2SCANH (id=659)
  GauXC_Functional_R2SCANh,

  /// @brief r2SCAN50: r2SCAN hybrid like BHLYP with 50% exact exchange
  /// - libxc name: HYB_MGGA_XC_R2SCAN50 (id=661)
  GauXC_Functional_R2SCAN50,

  /// @brief Minnesota 2005 hybrid functional
  /// - libxc names: HYB_MGGA_X_M05 (id=438) & MGGA_C_M05 (id=237)
  GauXC_Functional_M05,

  /// @brief Minnesota 2008 hybrid functional
  /// - libxc names: HYB_MGGA_X_M06 (id=449) & MGGA_C_M06 (id=235)
  GauXC_Functional_M06,

  /// @brief Minnesota M08 hybrid functional
  /// - libxc names: HYB_MGGA_X_M08_HX (id=295) & MGGA_C_M08_HX (id=78)
  GauXC_Functional_M08HX,

  /// @brief Minnesota M08-SO hybrid exchange functional
  /// - libxc names: HYB_MGGA_X_M08_SO (id=296) & MGGA_C_M08_SO (id=77)
  GauXC_Functional_M08SO,

  /// @brief Minnesota M05-2X hybrid exchange functional
  /// - libxc names: HYB_MGGA_X_M05_2X (id=439) & MGGA_C_M05_2X (id=238)
  GauXC_Functional_M052X,

  /// @brief Minnesota M06-SX short-range hybrid exchange functional
  /// - libxc names: HYB_MGGA_X_M06_SX (id=310) & MGGA_C_M06_SX (id=311)
  GauXC_Functional_M06SX,

  /// @brief Minnesota CF22D hybrid exchange functional
  /// - libxc names: HYB_MGGA_X_CF22D (id=340) & MGGA_C_CF22D (id=341)
  GauXC_Functional_CF22D,

  /// @brief Hybrid based on SOGGA11 form
  /// - libxc names: HYB_GGA_X_SOGGA11_X (id=426) & GGA_C_SOGGA11_X (id=159)
  GauXC_Functional_SOGGA11X,

  /// @brief Minnesota M06-HF hybrid exchange functional
  /// - libxc names: HYB_MGGA_X_M06_HF (id=444) & MGGA_C_M06_HF (id=234)
  GauXC_Functional_M06HF,

  /// @brief Minnesota M11 hybrid exchange functional
  /// - libxc names: HYB_MGGA_X_M11 (id=297) & MGGA_C_M11 (id=76)
  GauXC_Functional_M11,

  /// @brief Minnesota MN12-L exchange functional
  /// - libxc names: MGGA_X_MN12_L (id=227) & MGGA_C_MN12_L (id=74)
  GauXC_Functional_MN12L,

  /// @brief Minnesota MN12-SX hybrid exchange functional
  /// - libxc names: HYB_MGGA_X_MN12_SX (id=248) & MGGA_C_MN12_SX (id=73)
  GauXC_Functional_MN12SX,

  /// @brief Minnesota MN15 correlation functional
  /// - libxc names: HYB_MGGA_X_MN15 (id=268) & MGGA_C_MN15 (id=269)
  GauXC_Functional_MN15,

  /// @brief Minnesota MN15-L exchange functional
  /// - libxc names: MGGA_X_MN15_L (id=260) & MGGA_C_MN15_L (id=261)
  GauXC_Functional_MN15L,

  /// @brief Revised Minnesota 2006 meta-GGA functional
  /// - libxc names: MGGA_X_REVM06_L (id=293) & MGGA_C_REVM06_L (id=294)
  GauXC_Functional_revM06L,
};

/**
 * @brief GauXC C API Functional handle.
 */
typedef struct GauXCFunctional {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr;  ///< Pointer to the Functional instance.
} GauXCFunctional;

/**
 * @brief Create a GauXCFunctional from a string specification.
 * @param status Pointer to GauXCStatus for error handling.
 * @param functional_spec String specification of the functional.
 * @param polarized Whether the functional is spin-polarized.
 * @return A handle to the created GauXCFunctional.
 */
extern GauXCFunctional gauxc_functional_from_string(
  GauXCStatus* status,
  const char* functional_spec,
  bool polarized
);

/**
 * @brief Create a GauXCFunctional from a GauXC_Functional enum.
 * @param status Pointer to GauXCStatus for error handling.
 * @param functional_type The type of functional to create.
 * @param polarized Whether the functional is spin-polarized.
 * @return A handle to the created GauXCFunctional.
 */
extern GauXCFunctional gauxc_functional_from_enum(
  GauXCStatus* status,
  enum GauXC_Functional functional_type,
  bool polarized
);

/**
 * @brief Delete a GauXCFunctional handle.
 * @param status Pointer to GauXCStatus for error handling.
 * @param functional The GauXCFunctional handle to delete.
 */
extern void gauxc_functional_delete(
  GauXCStatus* status,
  GauXCFunctional* functional
);

#ifdef __cplusplus
} // namespace GauXC::C
} // extern "C"
#endif