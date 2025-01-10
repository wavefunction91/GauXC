/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

namespace GauXC {

struct XCDeviceTask {

  size_t npts = 0;

  struct screening_quantities {
    size_t nbe             = 0;
    size_t ncut            = 0;
    size_t nblock          = 0;
    size_t nshells         = 0;
    size_t ibf_begin       = 0;
    size_t* shell_list     = nullptr;
    size_t* shell_offs     = nullptr;
    int32_t* submat_cut    = nullptr;
    int32_t* submat_block  = nullptr;
  };

  screening_quantities bfn_screening;
  screening_quantities cou_screening;

  double* points_x       = nullptr;
  double* points_y       = nullptr;
  double* points_z       = nullptr;
  double* weights        = nullptr;

  double*   nbe_scr = nullptr;
  double*   zmat    = nullptr;
  double*   fmat    = nullptr;
  double*   gmat    = nullptr;
  double*   xmat_x  = nullptr;
  double*   xmat_y  = nullptr;
  double*   xmat_z  = nullptr;
  double*   bf      = nullptr;
  double*   dbfx    = nullptr;
  double*   dbfy    = nullptr;
  double*   dbfz    = nullptr;
  double*   d2bfxx    = nullptr;
  double*   d2bfxy    = nullptr;
  double*   d2bfxz    = nullptr;
  double*   d2bfyy    = nullptr;
  double*   d2bfyz    = nullptr;
  double*   d2bfzz    = nullptr;
  double*   eps     = nullptr;

  double* den       = nullptr;
  double* gamma     = nullptr;
  double*   vrho    = nullptr;
  double*   vgamma  = nullptr;
    
  // (S,Z,Y,X) densities
  double* den_s     = nullptr;
  double* den_z     = nullptr;
  double* den_y     = nullptr;
  double* den_x     = nullptr;
  // Del(S,Z,Y,X) Gradients
  double* dden_sx   = nullptr;
  double* dden_sy   = nullptr;
  double* dden_sz   = nullptr;
  double* dden_zx   = nullptr;
  double* dden_zy   = nullptr;
  double* dden_zz   = nullptr;
  double* dden_yx   = nullptr;
  double* dden_yy   = nullptr;
  double* dden_yz   = nullptr;
  double* dden_xx   = nullptr;
  double* dden_xy   = nullptr;
  double* dden_xz   = nullptr;
  
  // 2C U vars
  double* vrho_pos  = nullptr;
  double* vrho_neg  = nullptr;
  double* gamma_pp  = nullptr;
  double* gamma_pm  = nullptr;
  double* gamma_mm  = nullptr;
  double* vgamma_pp  = nullptr;
  double* vgamma_pm  = nullptr;
  double* vgamma_mm  = nullptr;

  // GKS K,H matrices
  double* K_z        = nullptr;
  double* K_y        = nullptr;
  double* K_x        = nullptr;
  double* H_z        = nullptr;
  double* H_y        = nullptr;
  double* H_x        = nullptr;

  // MGGA
  double*   d2bflapl    = nullptr;
  double*   d3bflapl_x    = nullptr;
  double*   d3bflapl_y    = nullptr;
  double*   d3bflapl_z    = nullptr;
  double*   denlapl     = nullptr;
  double*   tau   = nullptr;
  double*   vtau  = nullptr;
  double*   vlapl  = nullptr;

  int32_t iParent       = -1;
  double dist_nearest   = 0.;
  double * dist_scratch = nullptr;

  int32_t* bfn_shell_indirection = nullptr;

};

}
