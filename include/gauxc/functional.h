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
namespace GauXC::C {
extern "C" {
#endif

enum GauXC_Functional {
  /// @brief Slater exchange & Vosko, Wilk & Nusair correlation (VWN3)
  /// - P. A. M. Dirac., Math. Proc. Cambridge Philos. Soc. 26, 376 (1930) (doi: 10.1017/S0305004100016108)
  /// - F. Bloch., Z. Phys. 57, 545 (1929) (doi: 10.1007/BF01340281)
  /// - S. H. Vosko, L. Wilk, and M. Nusair., Can. J. Phys. 58, 1200 (1980) (doi: 10.1139/p80-159)
  /// - libxc names: LDA_X (id=1) and LDA_C_VWN_3 (id=30)
  GauXC_Functional_SVWN3,

  /// @brief Slater exchange & Vosko, Wilk & Nusair correlation (VWN5)
  /// - P. A. M. Dirac., Math. Proc. Cambridge Philos. Soc. 26, 376 (1930) (doi: 10.1017/S0305004100016108)
  /// - F. Bloch., Z. Phys. 57, 545 (1929) (doi: 10.1007/BF01340281)
  /// - S. H. Vosko, L. Wilk, and M. Nusair., Can. J. Phys. 58, 1200 (1980) (doi: 10.1139/p80-159)
  /// - libxc names: LDA_X (id=1) and LDA_C_VWN (id=7)
  /// - xcfun names: SLATERX and VWN5C
  GauXC_Functional_SVWN5,

  /// @brief Becke 88 exchange & Lee, Yang & Parr correlation
  /// - A. D. Becke., Phys. Rev. A 38, 3098 (1988) (doi: 10.1103/PhysRevA.38.3098)
  /// - C. Lee, W. Yang, and R. G. Parr., Phys. Rev. B 37, 785 (1988) (doi: 10.1103/PhysRevB.37.785)
  /// - B. Miehlich, A. Savin, H. Stoll, and H. Preuss., Chem. Phys. Lett. 157, 200 (1989) (doi: 10.1016/0009-2614(89)87234-3)
  /// - libxc names: GGA_X_B88 (id=106) and GGA_C_LYP (id=131)
  /// - xcfun names: BECKEX and LYPC
  GauXC_Functional_BLYP,

  /// @brief Becke 88 exchange & Lee, Yang & Parr correlation, 3-parameter hybrid
  /// - P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and M. J. Frisch., J. Phys. Chem. 98, 11623 (1994) (doi: 10.1021/j100096a001)
  /// - libxc name: HYB_GGA_XC_B3LYP (id=402)
  GauXC_Functional_B3LYP,

  /// @brief Perdew-Burke-Ernzerhof exchange & correlation
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
  /// - libxc names: GGA_X_PBE (id=101) & GGA_C_PBE (id=130)
  /// - xcfun names: PBEEX and PBEC
  GauXC_Functional_PBE,

  /// @brief revised Perdew-Burke-Ernzerhof exchange & original PBE correlation
  /// - Y. Zhang and W. Yang., Phys. Rev. Lett. 80, 890 (1998) (doi: 10.1103/PhysRevLett.80.890)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
  /// - libxc names: GGA_X_PBE_R (id=102) & GGA_C_PBE (id=130)
  GauXC_Functional_revPBE,

  /// @brief Perdew-Burke-Ernzerhof exchange & correlation, 1-parameter hybrid
  /// - C. Adamo and V. Barone., J. Chem. Phys. 110, 6158 (1999) (doi: 10.1063/1.478522)
  /// - M. Ernzerhof and G. E. Scuseria., J. Chem. Phys. 110, 5029 (1999) (doi: 10.1063/1.478401)
  /// - libxc name: HYB_GGA_XC_PBEH (id=406)
  GauXC_Functional_PBE0,

  /// @brief Strongly constrained and appropriately normed (SCAN) meta-GGA
  /// - J. Sun, A. Ruzsinszky, and J. P. Perdew., Phys. Rev. Lett. 115, 036402 (2015) (doi: 10.1103/PhysRevLett.115.036402)
  /// - libxc names: MGGA_X_SCAN (id=263) & MGGA_C_SCAN (id=267)
  GauXC_Functional_SCAN,

  /// @brief Regularized & restored strongly constrained and appropriately normed (R2SCAN) meta-GGA
  /// - J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, and J. Sun., J. Phys. Chem. Lett. 11, 8208-8215 (2020) (doi: 10.1021/acs.jpclett.0c02405)
  /// - J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, and J. Sun., J. Phys. Chem. Lett. 11, 9248-9248 (2020) (doi: 10.1021/acs.jpclett.0c03077)
  /// - libxc names: MGGA_X_R2SCAN (id=497) & MGGA_C_R2SCAN (id=498) 
  GauXC_Functional_R2SCAN,

  /// @brief Regularized & restored strongly constrained and appropriately normed (R2SCAN) meta-GGA, deorbitalized version
  /// - D. Mejía-Rodríguez and S. B. Trickey., Phys. Rev. B 102, 121109 (2020) (doi: 10.1103/PhysRevB.102.121109)
  /// - J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, and J. Sun., J. Phys. Chem. Lett. 11, 8208-8215 (2020) (doi: 10.1021/acs.jpclett.0c02405)
  /// - J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, and J. Sun., J. Phys. Chem. Lett. 11, 9248-9248 (2020) (doi: 10.1021/acs.jpclett.0c03077)
  /// - libxc names: MGGA_X_R2SCANL (id=718) & MGGA_C_R2SCANL (id=719)
  GauXC_Functional_R2SCANL,

  /// @brief Minnesota 2006 hybrid functional
  /// - Y. Zhao and D. G. Truhlar., Theor. Chem. Acc. 120, 215 (2008) (doi: 10.1007/s00214-007-0310-x)
  /// - libxc names: HYB_MGGA_X_M06_2X (id=450) & MGGA_C_M06_2X (id=236)
  /// - xcfun names: M062X and M062C
  GauXC_Functional_M062X,

  /// @brief Perdew, Kurth, Zupan, and Blaha
  /// - J. P. Perdew, S. Kurth, A. Zupan, and P. Blaha., Phys. Rev. Lett. 82, 2544 (1999) (doi: 10.1103/PhysRevLett.82.2544)
  /// - libxc names: MGGA_X_PKZB (id=213) & MGGA_C_PKZB (id=239)
  GauXC_Functional_PKZB,

  /// @brief epc17(-1): electron-proton correlation 2017
  /// - Y. Yang, K. R. Brorsen, T. Culpitt, M. V. Pak, and S. Hammes-Schiffer., J. Chem. Phys. 147, 114113 (2017) (doi: 10.1063/1.4996038)
  /// - libxc name: LDA_C_EPC17 (id=328)
  GauXC_Functional_EPC17_1,

  /// @brief epc17-2: electron-proton correlation 2017 for proton affinities
  /// - K. R. Brorsen, Y. Yang, and S. Hammes-Schiffer., J. Phys. Chem. Lett. 8, 3488-3493 (2017) (doi: 10.1021/acs.jpclett.7b01442)
  /// - libxc name: LDA_C_EPC17_2 (id=329)
  GauXC_Functional_EPC17_2,

  /// @brief epc18-1: electron-proton correlation 2018
  /// - K. R. Brorsen, P. E. Schneider, and S. Hammes-Schiffer., J. Chem. Phys. 149, 044110 (2018) (doi: 10.1063/1.5037945)
  /// - libxc name: LDA_C_EPC18_1 (id=330)
  GauXC_Functional_EPC18_1,

  /// @brief epc18-2: electron-proton correlation 2018 for proton affinities
  /// - K. R. Brorsen, P. E. Schneider, and S. Hammes-Schiffer., J. Chem. Phys. 149, 044110 (2018) (doi: 10.1063/1.5037945)
  /// - libxc name: LDA_C_EPC18_2 (id=331)
  GauXC_Functional_EPC18_2,

  /// @brief Grimme's parametrization of the B97 functional, original D2 variant
  /// - S. Grimme., J. Comput. Chem. 27, 1787 (2006) (doi: 10.1002/jcc.20495)
  /// - libxc name: GGA_XC_B97_D (id=170)
  GauXC_Functional_B97D,

  /// @brief Grimme's parametrization of the B97 functional, D3(0) variant
  /// - S. Grimme., J. Comput. Chem. 27, 1787 (2006) (doi: 10.1002/jcc.20495)
  /// - libxc name: GGA_XC_B97_D (id=170)
  GauXC_Functional_B97D3ZERO,

  /// @brief Coulomb-attenuating method range-separated hybrid functional
  /// - T. Yanai, D. P. Tew, and N. C. Handy., Chem. Phys. Lett. 393, 51 (2004) (doi: 10.1016/j.cplett.2004.06.011)
  /// - libxc name: HYB_GGA_XC_CAM_B3LYP (id=433)
  GauXC_Functional_CAMB3LYP,

  /// @brief Slater exchange
  /// - P. A. M. Dirac., Math. Proc. Cambridge Philos. Soc. 26, 376 (1930) (doi: 10.1017/S0305004100016108)
  /// - F. Bloch., Z. Phys. 57, 545 (1929) (doi: 10.1007/BF01340281)
  /// - libxc name: LDA_X (id=1)
  /// - xcfun name: SLATERX
  GauXC_Functional_LDA,

  /// @brief Minnesota 2006 meta-GGA functional
  /// - Y. Zhao and D. G. Truhlar., J. Chem. Phys. 125, 194101 (2006) (doi: 10.1063/1.2370993)
  /// - Y. Zhao and D. G. Truhlar., Theor. Chem. Acc. 120, 215 (2008) (doi: 10.1007/s00214-007-0310-x)
  /// - libxc names: MGGA_X_M06_L (id=203) & MGGA_C_M06_L (id=233)
  GauXC_Functional_M06L,

  /// @brief Strongly constrained and appropriately normed (SCAN) meta-GGA, 1-parameter hybrid
  /// - K. Hui and J.-D. Chai., J. Chem. Phys. 144, 044114 (2016) (doi: 10.1063/1.4940734)
  /// - J. Sun, A. Ruzsinszky, and J. P. Perdew., Phys. Rev. Lett. 115, 036402 (2015) (doi: 10.1103/PhysRevLett.115.036402)
  /// - libxc names: HYB_MGGA_X_SCAN0 (id=264) & MGGA_C_SCAN (id=267)
  GauXC_Functional_SCAN0,

  /// @brief Slater exchange & Perdew, Wang 92 correlation
  /// - P. A. M. Dirac., Math. Proc. Cambridge Philos. Soc. 26, 376 (1930) (doi: 10.1017/S0305004100016108)
  /// - F. Bloch., Z. Phys. 57, 545 (1929) (doi: 10.1007/BF01340281)
  /// - J. P. Perdew and Y. Wang., Phys. Rev. B 45, 13244 (1992) (doi: 10.1103/PhysRevB.45.13244)
  /// - libxc names: LDA_X (id=1) and LDA_C_PW (id=12)
  GauXC_Functional_SPW92,

  /// @brief Tao, Perdew, Staroverov & Scuseria meta-GGA
  /// - J. Tao, J. P. Perdew, V. N. Staroverov, and G. E. Scuseria., Phys. Rev. Lett. 91, 146401 (2003) (doi: 10.1103/PhysRevLett.91.146401)
  /// - J. P. Perdew, J. Tao, V. N. Staroverov, and G. E. Scuseria., J. Chem. Phys. 120, 6898 (2004) (doi: 10.1063/1.1665298)
  /// - libxc names: MGGA_X_TPSS (id=202) and MGGA_C_TPSS (id=231)
  GauXC_Functional_TPSS,

  /// @brief Tao, Perdew, Staroverov & Scuseria meta-GGA, 1-parameter hybrid
  /// - V. N. Staroverov, G. E. Scuseria, J. Tao, and J. P. Perdew., J. Chem. Phys. 119, 12129 (2003) (doi: 10.1063/1.1626543)
  /// - libxc name: HYB_MGGA_XC_TPSSH (id=457)
  GauXC_Functional_TPSSh,

  /// @brief Tao, Perdew, Staroverov & Scuseria meta-GGA, 1-parameter hybrid
  /// - S. Grimme., J. Phys. Chem. A 109, 3067-3077 (2005) (doi: 10.1021/jp050036j)
  /// - libxc name: HYB_MGGA_XC_TPSS0 (id=396)
  GauXC_Functional_TPSS0,

  /// @brief Vosko, Wilk & Nusair correlation (VWN3)
  /// - S. H. Vosko, L. Wilk, and M. Nusair., Can. J. Phys. 58, 1200 (1980) (doi: 10.1139/p80-159)
  /// - libxc name: LDA_C_VWN_3 (id=30)
  GauXC_Functional_VWN3,

  /// @brief Vosko, Wilk & Nusair correlation (VWN5)
  /// - S. H. Vosko, L. Wilk, and M. Nusair., Can. J. Phys. 58, 1200 (1980) (doi: 10.1139/p80-159)
  /// - libxc name: LDA_C_VWN (id=7)
  GauXC_Functional_VWN5,

  /// @brief HJS screened PBE exchange & original PBE correlation
  /// - T. M. Henderson, B. G. Janesko, and G. E. Scuseria., J. Chem. Phys. 128, 194105 (2008) (doi: 10.1063/1.2921797)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
  /// - libxc names: GGA_X_HJS_PBE (id=525) & GGA_C_PBE (id=130)
  GauXC_Functional_LRCwPBE,

  /// @brief HJS screened PBE exchange & original PBE correlation, hybrid version
  /// - T. M. Henderson, B. G. Janesko, and G. E. Scuseria., J. Chem. Phys. 128, 194105 (2008) (doi: 10.1063/1.2921797)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
  /// - libxc name: HYB_GGA_XC_HJS_PBE (id=429)
  GauXC_Functional_LRCwPBEh,

  /// @brief Becke 88 exchange and Perdew 86 correlation
  /// - A. D. Becke., Phys. Rev. A 38, 3098 (1988) (doi: 10.1103/PhysRevA.38.3098)
  /// - J. P. Perdew., Phys. Rev. B 33, 8822 (1986) (doi: 10.1103/PhysRevB.33.8822)
  /// - libxc names: GGA_X_B88 (id=106) and GGA_C_P86 (id=132)
  GauXC_Functional_BP86,

  /// @brief Heyd-Scuseria-Ernzerhof screened hybrid functional (HSE03)
  /// - J. Heyd, G. E. Scuseria, and M. Ernzerhof., J. Chem. Phys. 118, 8207 (2003) (doi: 10.1063/1.1564060)
  /// - J. Heyd, G. E. Scuseria, and M. Ernzerhof., J. Chem. Phys. 124, 219906 (2006) (doi: 10.1063/1.2204597)
  /// - libxc name: HYB_GGA_XC_HSE03 (id=427)
  GauXC_Functional_HSE03,

  /// @brief Heyd-Scuseria-Ernzerhof screened hybrid functional (HSE06)
  /// - J. Heyd, G. E. Scuseria, and M. Ernzerhof., J. Chem. Phys. 118, 8207 (2003) (doi: 10.1063/1.1564060)
  /// - J. Heyd, G. E. Scuseria, and M. Ernzerhof., J. Chem. Phys. 124, 219906 (2006) (doi: 10.1063/1.2204597)
  /// - A. V. Krukau, O. A. Vydrov, A. F. Izmaylov, and G. E. Scuseria., J. Chem. Phys. 125, 224106 (2006) (doi: 10.1063/1.2404663)
  /// - libxc name: HYB_GGA_XC_HSE06 (id=428)
  GauXC_Functional_HSE06,

  /// @brief Revised B3LYP
  /// - L. Lu, H. Hu, H. Hou, and B. Wang., Comput. Theor. Chem. 1015, 64 (2013) (doi: 10.1016/j.comptc.2013.04.009)
  /// - libxc name: HYB_GGA_XC_REVB3LYP (id=454)
  GauXC_Functional_revB3LYP,

  /// @brief revised Perdew-Burke-Ernzerhof exchange & original PBE correlation, hybrid version
  /// - Y. Zhang and W. Yang., Phys. Rev. Lett. 80, 890 (1998) (doi: 10.1103/PhysRevLett.80.890)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
  /// - libxc names: GGA_X_PBE_R (id=102) & GGA_C_PBE (id=130)
  GauXC_Functional_revPBE0,

  /// @brief revised Tao, Perdew, Staroverov & Scuseria
  /// - J. P. Perdew, A. Ruzsinszky, G. I. Csonka, L. A. Constantin, and J. Sun., Phys. Rev. Lett. 103, 026403 (2009) (doi: 10.1103/PhysRevLett.103.026403)
  /// - J. P. Perdew, A. Ruzsinszky, G. I. Csonka, L. A. Constantin, and J. Sun., Phys. Rev. Lett. 106, 179902 (2011) (doi: 10.1103/PhysRevLett.106.179902)
  /// - libxc names: MGGA_X_REVTPSS (id=212) & MGGA_C_REVTPSS (id=241)
  GauXC_Functional_revTPSS,

  /// @brief revTPSSh
  /// - G. I. Csonka, J. P. Perdew, and A. Ruzsinszky., J. Chem. Theory Comput. 6, 3688 (2010) (doi: 10.1021/ct100488v)
  /// - libxc name: HYB_MGGA_XC_REVTPSSH (id=458)
  GauXC_Functional_revTPSSh,

  /// @brief Perdew-Wang 91 exchange and correlation
  /// - J. P. Perdew. In P. Ziesche and H. Eschrig, editors, Proceedings of the 75. WE-Heraeus-Seminar and 21st Annual International Symposium on Electronic Structure of Solids, 11. Berlin, 1991. Akademie Verlag.
  /// - J. P. Perdew, J. A. Chevary, S. H. Vosko, K. A. Jackson, M. R. Pederson, D. J. Singh, and C. Fiolhais., Phys. Rev. B 46, 6671 (1992) (doi: 10.1103/PhysRevB.46.6671)
  /// - J. P. Perdew, J. A. Chevary, S. H. Vosko, K. A. Jackson, M. R. Pederson, D. J. Singh, and C. Fiolhais., Phys. Rev. B 48, 4978 (1993) (doi: 10.1103/PhysRevB.48.4978.2)
  /// - libxc name: GGA_X_PW91 (id=109) and GGA_C_PW91 (id=134)
  GauXC_Functional_PW91,

  /// @brief mBEEF exchange and Perdew, Burke & Ernzerhof SOL
  /// - J. Wellendorff, K. T. Lundgaard, K. W. Jacobsen, and T. Bligaard., J. Chem. Phys. 140, 144107 (2014) (doi: 10.1063/1.4870397)
  /// - J. P. Perdew, A. Ruzsinszky, G. I. Csonka, O. A. Vydrov, G. E. Scuseria, L. A. Constantin, X. Zhou, and K. Burke., Phys. Rev. Lett. 100, 136406 (2008) (doi: 10.1103/PhysRevLett.100.136406)
  /// - libxc names: MGGA_X_MBEEF (id=249) and GGA_C_PBE_SOL (id=133)
  GauXC_Functional_mBEEF,

  /// @brief The original (ACM, B3PW91) hybrid of Becke
  /// - A. D. Becke., J. Chem. Phys. 98, 5648 (1993) (doi: 10.1063/1.464913)
  /// - libxc name: HYB_GGA_XC_B3PW91 (id=401)
  GauXC_Functional_B3PW91,

  /// @brief O3LYP
  /// - W.-M. Hoe, A. J. Cohen, and N. C. Handy., Chem. Phys. Lett. 341, 319–328 (2001) (doi: 10.1016/S0009-2614(01)00581-4)
  /// - A. J. Cohen and N. C. Handy., Mol. Phys. 99, 607 (2001) (doi: 10.1080/00268970010023435)
  /// - libxc name: HYB_GGA_XC_O3LYP (id=404)
  GauXC_Functional_O3LYP,

  /// @brief Handy & Cohen OPTX 01 exchange and Lee, Yang & Parr correlation
  /// - N. C. Handy and A. J. Cohen., Mol. Phys. 99, 403 (2001) (doi: 10.1080/00268970010018431)
  /// - C. Lee, W. Yang, and R. G. Parr., Phys. Rev. B 37, 785 (1988) (doi: 10.1103/PhysRevB.37.785)
  /// - B. Miehlich, A. Savin, H. Stoll, and H. Preuss., Chem. Phys. Lett. 157, 200 (1989) (doi: 10.1016/0009-2614(89)87234-3)
  /// - libxc names: GGA_X_OPTX (id=110) & GGA_C_LYP (id=131)
  GauXC_Functional_OLYP,

  /// @brief Handy & Cohen OPTX 01 exchange and Perdew, Burke & Ernzerhof correlation
  /// - N. C. Handy and A. J. Cohen., Mol. Phys. 99, 403 (2001) (doi: 10.1080/00268970010018431)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
  /// - libxc names: GGA_X_OPTX (id=110) & GGA_C_PBE (id=130)
  GauXC_Functional_OPBE,

  /// @brief mPW1K
  /// - B. J. Lynch, P. L. Fast, M. Harris, and D. G. Truhlar., J. Phys. Chem. A 104, 4811 (2000) (doi: 10.1021/jp000497z)
  /// - libxc name: HYB_GGA_XC_MPW1K (id=405)
  GauXC_Functional_MPW1K,

  /// @brief Revised Perdew-Burke-Ernzerhof exchange by Hammer, Hansen, and Norskov
  /// - B. Hammer, L. B. Hansen, and J. K. Nørskov., Phys. Rev. B 59, 7413 (1999) (doi: 10.1103/PhysRevB.59.7413)
  /// - libxc name: GGA_X_RPBE (id=117)
  GauXC_Functional_RPBE,

  /// @brief Becke 88 exchange
  /// - A. D. Becke., Phys. Rev. A 38, 3098 (1988) (doi: 10.1103/PhysRevA.38.3098)
  /// - libxc name: GGA_X_B88 (id=106)
  GauXC_Functional_B88,

  /// @brief modified Perdew-Wang 91 exchange by Adamo & Barone
  /// - C. Adamo and V. Barone., J. Chem. Phys. 108, 664 (1998) (doi: 10.1063/1.475428)
  /// - libxc name: GGA_X_MPW91 (id=119)
  GauXC_Functional_MPW91,

  /// @brief Regularized strongly constrained and appropriately normed (RSCAN) meta-GGA by Bartok and Yates
  /// - A. P. Bartók and J. R. Yates., J. Chem. Phys. 150, 161101 (2019) (doi: 10.1063/1.5094646)
  /// - libxc names: MGGA_X_RSCAN (id=493) and MGGA_C_RSCAN (id=494)
  GauXC_Functional_RSCAN,

  /// @brief CAM version of B3LYP, tuned for excitations and properties
  /// - K. Okuno, Y. Shigeta, R. Kishi, H. Miyasaka, and M. Nakano., J. Photochem. Photobiol., A 235, 29 (2012) (doi: 10.1016/j.jphotochem.2012.03.003)
  /// - libxc name: HYB_GGA_XC_TUNED_CAM_B3LYP (id=434)
  GauXC_Functional_TUNEDCAMB3LYP,

  /// @brief wB97 range-separated functional
  /// - J.-D. Chai and M. Head-Gordon., J. Chem. Phys. 128, 084106 (2008) (doi: 10.1063/1.2834918)
  /// - libxc name: HYB_GGA_XC_WB97 (id=463)
  GauXC_Functional_wB97,

  /// @brief wB97X range-separated functional
  /// - J.-D. Chai and M. Head-Gordon., J. Chem. Phys. 128, 084106 (2008) (doi: 10.1063/1.2834918)
  /// - libxc name: HYB_GGA_XC_WB97X (id=464)
  GauXC_Functional_wB97X,

  /// @brief wB97X-D range-separated functional
  /// - J.-D. Chai and M. Head-Gordon., Phys. Chem. Chem. Phys. 10, 6615-6620 (2008) (doi: 10.1039/B810189B)
  /// - libxc name: HYB_GGA_XC_WB97X_D (id=471)
  GauXC_Functional_wB97XD,

  /// @brief wB97X-D3 range-separated functional
  /// - Y.-S. Lin, G.-D. Li, S.-P. Mao, and J.-D. Chai., J. Chem. Theory Comput. 9, 263-272 (2013) (doi: 10.1021/ct300715s)
  /// - libxc name: HYB_GGA_XC_WB97X_D3 (id=399)
  GauXC_Functional_wB97XD3,

  /// @brief Long-range corrected PBE (LC-wPBE) by Vydrov and Scuseria
  /// - O. A. Vydrov and G. E. Scuseria., J. Chem. Phys. 125, 234109 (2006) (doi: 10.1063/1.2409292)
  /// - libxc name: HYB_GGA_XC_LC_WPBE (id=478)
  GauXC_Functional_LCwPBE,

  /// @brief X3LYP
  /// - X. Xu and W. A. Goddard., Proc. Natl. Acad. Sci. U. S. A. 101, 2673 (2004) (doi: 10.1073/pnas.0308730100)
  /// - libxc name: HYB_GGA_XC_X3LYP (id=411)
  GauXC_Functional_X3LYP,

  /// @brief XLYP
  /// - X. Xu and W. A. Goddard., Proc. Natl. Acad. Sci. U. S. A. 101, 2673 (2004) (doi: 10.1073/pnas.0308730100)
  /// - libxc name: GGA_XC_XLYP (id=166)
  GauXC_Functional_XLYP,

  /// @brief BHandH i.e. BHLYP
  /// - A. D. Becke., J. Chem. Phys. 98, 1372 (1993) (doi: 10.1063/1.464304)
  /// - Defined through Gaussian implementation.
  /// - libxc name: HYB_GGA_XC_BHANDH (id=435)
  GauXC_Functional_BHANDH,

  /// @brief Boese-Martin for kinetics
  /// - A. D. Boese and J. M. L. Martin., J. Chem. Phys. 121, 3405 (2004) (doi: 10.1063/1.1774975)
  /// - libxc names: HYB_MGGA_X_BMK (id=279) & GGA_C_BMK (id=280)
  GauXC_Functional_BMK,

  /// @brief Becke 88 exchange and Perdew 86 based on VWN5 correlation, with more accurate value for ftilde
  /// - A. D. Becke., Phys. Rev. A 38, 3098 (1988) (doi: 10.1103/PhysRevA.38.3098)
  /// - J. P. Perdew., Phys. Rev. B 33, 8822 (1986) (doi: 10.1103/PhysRevB.33.8822)
  /// - libxc names: GGA_X_B88 (id=106) & GGA_C_P86VWN_FT (id=253)
  GauXC_Functional_BP86VWN,

  /// @brief Mixture of PW86 with BC95
  /// - A. D. Becke., J. Chem. Phys. 104, 1040 (1996) (doi: 10.1063/1.470829)
  /// - libxc name: HYB_MGGA_XC_PW86B95 (id=442)
  GauXC_Functional_PW86B95,

  /// @brief Perdew & Wang 86 exchange and PBE correlation
  /// - J. P. Perdew and W. Yue., Phys. Rev. B 33, 8800 (1986) (doi: 10.1103/PhysRevB.33.8800)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 77, 3865 (1996) (doi: 10.1103/PhysRevLett.77.3865)
  /// - J. P. Perdew, K. Burke, and M. Ernzerhof., Phys. Rev. Lett. 78, 1396 (1997) (doi: 10.1103/PhysRevLett.78.1396)
  /// - libxc names: GGA_X_PW86 (id=108) & GGA_C_PBE (id=130)
  GauXC_Functional_PW86PBE,

  /// @brief r2SCAN0: r2SCAN hybrid like PBE0 with 25% exact exchange
  /// - M. Bursch, H. Neugebauer, S. Ehlert, and S. Grimme., J. Chem. Phys. 156, 134105 (2022) (doi: 10.1063/5.0086040)
  /// - libxc name: HYB_MGGA_XC_R2SCAN0 (id=660)
  GauXC_Functional_R2SCAN0,

  /// @brief r2SCANh: r2SCAN hybrid like TPSSh with 10% exact exchange
  /// - M. Bursch, H. Neugebauer, S. Ehlert, and S. Grimme., J. Chem. Phys. 156, 134105 (2022) (doi: 10.1063/5.0086040)
  /// - libxc name: HYB_MGGA_XC_R2SCANH (id=659)
  GauXC_Functional_R2SCANh,

  /// @brief r2SCAN50: r2SCAN hybrid like BHLYP with 50% exact exchange
  /// - M. Bursch, H. Neugebauer, S. Ehlert, and S. Grimme., J. Chem. Phys. 156, 134105 (2022) (doi: 10.1063/5.0086040)
  /// - libxc name: HYB_MGGA_XC_R2SCAN50 (id=661)
  GauXC_Functional_R2SCAN50,

  /// @brief Minnesota 2005 hybrid functional
  /// - Y. Zhao, N. E. Schultz, and D. G. Truhlar., J. Chem. Phys. 123, 161103 (2005) (doi: 10.1063/1.2126975)
  /// - libxc names: HYB_MGGA_X_M05 (id=438) & MGGA_C_M05 (id=237)
  GauXC_Functional_M05,

  /// @brief Minnesota 2008 hybrid functional
  /// - Y. Zhao and D. G. Truhlar., Theor. Chem. Acc. 120, 215 (2008) (doi: 10.1007/s00214-007-0310-x)
  /// - libxc names: HYB_MGGA_X_M06 (id=449) & MGGA_C_M06 (id=235)
  GauXC_Functional_M06,

  /// @brief Minnesota M08 hybrid functional
  /// - Y. Zhao and D. G. Truhlar., J. Chem. Theory Comput. 4, 1849 (2008) (doi: 10.1021/ct800246v)
  /// - libxc names: HYB_MGGA_X_M08_HX (id=295) & MGGA_C_M08_HX (id=78)
  GauXC_Functional_M08HX,

  /// @brief Minnesota M08-SO hybrid exchange functional
  /// - Y. Zhao and D. G. Truhlar., J. Chem. Theory Comput. 4, 1849 (2008) (doi: 10.1021/ct800246v)
  /// - libxc names: HYB_MGGA_X_M08_SO (id=296) & MGGA_C_M08_SO (id=77)
  GauXC_Functional_M08SO,

  /// @brief Minnesota M05-2X hybrid exchange functional
  /// - Y. Zhao, N. E. Schultz, and D. G. Truhlar., J. Chem. Theory Comput. 2, 364 (2006) (doi: 10.1021/ct0502763)
  /// - libxc names: HYB_MGGA_X_M05_2X (id=439) & MGGA_C_M05_2X (id=238)
  GauXC_Functional_M052X,

  /// @brief Minnesota M06-SX short-range hybrid exchange functional
  /// - Y. Wang, P. Verma, L. Zhang, Y. Li, Z. Liu, D. G. Truhlar, and X. He., Proc. Natl. Acad. Sci. U. S. A. 117, 2294–2301 (2020) (doi: 10.1073/pnas.1913699117)
  /// - libxc names: HYB_MGGA_X_M06_SX (id=310) & MGGA_C_M06_SX (id=311)
  GauXC_Functional_M06SX,

  /// @brief Minnesota CF22D hybrid exchange functional
  /// - Y. Liu, C. Zhang, Z. Liu, D. G. Truhlar, Y. Wang, and X. He., Nature Computational Science 3, 48–58 (2022) (doi: 10.1038/s43588-022-00371-5)
  /// - libxc names: HYB_MGGA_X_CF22D (id=340) & MGGA_C_CF22D (id=341)
  GauXC_Functional_CF22D,

  /// @brief Hybrid based on SOGGA11 form
  /// - R. Peverati and D. G. Truhlar., J. Chem. Phys. 135, 191102 (2011) (doi: 10.1063/1.3663871)
  /// - libxc names: HYB_GGA_X_SOGGA11_X (id=426) & GGA_C_SOGGA11_X (id=159)
  GauXC_Functional_SOGGA11X,

  /// @brief Minnesota M06-HF hybrid exchange functional
  /// - Y. Zhao and D. G. Truhlar., J. Phys. Chem. A 110, 13126 (2006) (doi: 10.1021/jp066479k)
  /// - libxc names: HYB_MGGA_X_M06_HF (id=444) & MGGA_C_M06_HF (id=234)
  GauXC_Functional_M06HF,

  /// @brief Minnesota M11 hybrid exchange functional
  /// - R. Peverati and D. G. Truhlar., J. Phys. Chem. Lett. 2, 2810 (2011) (doi: 10.1021/jz201170d)
  /// - libxc names: HYB_MGGA_X_M11 (id=297) & MGGA_C_M11 (id=76)
  GauXC_Functional_M11,

  /// @brief Minnesota MN12-L exchange functional
  /// - R. Peverati and D. G. Truhlar., Phys. Chem. Chem. Phys. 14, 13171 (2012) (doi: 10.1039/C2CP42025B)
  /// - libxc names: MGGA_X_MN12_L (id=227) & MGGA_C_MN12_L (id=74)
  GauXC_Functional_MN12L,

  /// @brief Minnesota MN12-SX hybrid exchange functional
  /// - R. Peverati and D. G. Truhlar., Phys. Chem. Chem. Phys. 14, 16187 (2012) (doi: 10.1039/C2CP42576A)
  /// - libxc names: HYB_MGGA_X_MN12_SX (id=248) & MGGA_C_MN12_SX (id=73)
  GauXC_Functional_MN12SX,

  /// @brief Minnesota MN15 correlation functional
  /// - H. S. Yu, X. He, S. L. Li, and D. G. Truhlar., Chem. Sci. 7, 5032-5051 (2016) (doi: 10.1039/C6SC00705H)
  /// - libxc names: HYB_MGGA_X_MN15 (id=268) & MGGA_C_MN15 (id=269)
  GauXC_Functional_MN15,

  /// @brief Minnesota MN15-L exchange functional
  /// - H. S. Yu, X. He, and D. G. Truhlar., J. Chem. Theory Comput. 12, 1280-1293 (2016) (doi: 10.1021/acs.jctc.5b01082)
  /// - libxc names: MGGA_X_MN15_L (id=260) & MGGA_C_MN15_L (id=261)
  GauXC_Functional_MN15L,

  /// @brief Revised Minnesota 2006 meta-GGA functional
  /// - Y. Wang, X. Jin, H. S. Yu, D. G. Truhlar, and X. He., Proc. Natl. Acad. Sci. U. S. A. 114, 8487-8492 (2017) (doi: 10.1073/pnas.1705670114)
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
GauXCFunctional gauxc_functional_from_string(
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
GauXCFunctional gauxc_functional_from_enum(
  GauXCStatus* status,
  enum GauXC_Functional functional_type,
  bool polarized
);

/**
 * @brief Delete a GauXCFunctional handle.
 * @param status Pointer to GauXCStatus for error handling.
 * @param functional The GauXCFunctional handle to delete.
 */
void gauxc_functional_delete(
  GauXCStatus* status,
  GauXCFunctional* functional
);

#ifdef __cplusplus
} // extern "C"
} // namespace GauXC::C
#endif