/*
 * BSD 3-Clause License
 * 
 * Copyright (c) 2017, Daniel Smith
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * * Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * This is a Gau2Grid automatically generated C file.
 *
 * More details can found at the following repo:
 *   https://github.com/dgasmith/gau2grid
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef GAU2GRID_GUARD_H
#define GAU2GRID_GUARD_H

#include "gau2grid/gau2grid_pragma.h"

// Order definitions
#define GG_SPHERICAL_CCA 300
#define GG_SPHERICAL_GAUSSIAN 301
#define GG_CARTESIAN_CCA 400
#define GG_CARTESIAN_MOLDEN 401
// Information helpers
int gg_max_L();

int gg_ncomponents(const int L, const int spherical);

// Fast transposers
void gg_naive_transpose(unsigned long n, unsigned long m, const double* PRAGMA_RESTRICT input, double* PRAGMA_RESTRICT output);
void gg_fast_transpose(unsigned long n, unsigned long m, const double* PRAGMA_RESTRICT input, double* PRAGMA_RESTRICT output);

// Fast segment copiers
void block_copy(unsigned long n, unsigned long m, const double* PRAGMA_RESTRICT input, unsigned long is, double* PRAGMA_RESTRICT output, unsigned long os, const int trans);


// Orbitals on a grid
void gg_orbitals(int L, const double* PRAGMA_RESTRICT C, const unsigned long norbitals, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT orbital_out);

// Collocation matrix functions
void gg_collocation(int L, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out);

void gg_collocation_deriv1(int L, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out);

void gg_collocation_deriv2(int L, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out);

void gg_collocation_deriv3(int L, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out, double* PRAGMA_RESTRICT phi_xxx_out, double* PRAGMA_RESTRICT phi_xxy_out, double* PRAGMA_RESTRICT phi_xxz_out, double* PRAGMA_RESTRICT phi_xyy_out, double* PRAGMA_RESTRICT phi_xyz_out, double* PRAGMA_RESTRICT phi_xzz_out, double* PRAGMA_RESTRICT phi_yyy_out, double* PRAGMA_RESTRICT phi_yyz_out, double* PRAGMA_RESTRICT phi_yzz_out, double* PRAGMA_RESTRICT phi_zzz_out);

#ifdef __cplusplus
}
#endif
#endif /* GAU2GRID_GUARD_H */