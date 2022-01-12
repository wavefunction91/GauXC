/*
 * This is a Gau2Grid automatically generated C file.
 *
 * More details can found at the following repo:
 *   https://github.com/dgasmith/gau2grid
 */

#include <math.h>
#if defined(__clang__) && defined(_MSC_VER)
#include <malloc.h>
#elif defined __clang__
#include <mm_malloc.h>
#elif defined _MSC_VER
#include <malloc.h>
#else
#include <stdlib.h>
#endif

#include "gau2grid/gau2grid.h"
#include "gau2grid/gau2grid_utility.h"
#include "gau2grid/gau2grid_pragma.h"

// Information helpers
int gg_max_L() { return 6; }

int gg_ncomponents(const int L, const int spherical) {
    if (spherical) {
    return 2 * L + 1;
    } else {
    return (L + 2) * (L + 1) / 2;
    }
}

// Collocation selector functions
void gg_orbitals(int L, const double* PRAGMA_RESTRICT C, const unsigned long norbitals, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT orbital_out) {
    // Chooses the correct function for a given L
    if (L == 0) {
        gg_orbitals_L0(C, norbitals, npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, orbital_out);
    } else if (L == 1) {
        gg_orbitals_L1(C, norbitals, npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, orbital_out);
    } else if (L == 2) {
        gg_orbitals_L2(C, norbitals, npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, orbital_out);
    } else if (L == 3) {
        gg_orbitals_L3(C, norbitals, npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, orbital_out);
    } else if (L == 4) {
        gg_orbitals_L4(C, norbitals, npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, orbital_out);
    } else if (L == 5) {
        gg_orbitals_L5(C, norbitals, npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, orbital_out);
    } else if (L == 6) {
        gg_orbitals_L6(C, norbitals, npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, orbital_out);
    } else {
        exit(0);
    }
}
void gg_collocation(int L, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out) {
    // Chooses the correct function for a given L
    if (L == 0) {
        gg_collocation_L0(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out);
    } else if (L == 1) {
        gg_collocation_L1(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out);
    } else if (L == 2) {
        gg_collocation_L2(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out);
    } else if (L == 3) {
        gg_collocation_L3(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out);
    } else if (L == 4) {
        gg_collocation_L4(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out);
    } else if (L == 5) {
        gg_collocation_L5(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out);
    } else if (L == 6) {
        gg_collocation_L6(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out);
    } else {
        exit(0);
    }
}
void gg_collocation_deriv1(int L, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out) {
    // Chooses the correct function for a given L
    if (L == 0) {
        gg_collocation_L0_deriv1(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out);
    } else if (L == 1) {
        gg_collocation_L1_deriv1(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out);
    } else if (L == 2) {
        gg_collocation_L2_deriv1(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out);
    } else if (L == 3) {
        gg_collocation_L3_deriv1(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out);
    } else if (L == 4) {
        gg_collocation_L4_deriv1(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out);
    } else if (L == 5) {
        gg_collocation_L5_deriv1(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out);
    } else if (L == 6) {
        gg_collocation_L6_deriv1(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out);
    } else {
        exit(0);
    }
}
void gg_collocation_deriv2(int L, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out) {
    // Chooses the correct function for a given L
    if (L == 0) {
        gg_collocation_L0_deriv2(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out);
    } else if (L == 1) {
        gg_collocation_L1_deriv2(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out);
    } else if (L == 2) {
        gg_collocation_L2_deriv2(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out);
    } else if (L == 3) {
        gg_collocation_L3_deriv2(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out);
    } else if (L == 4) {
        gg_collocation_L4_deriv2(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out);
    } else if (L == 5) {
        gg_collocation_L5_deriv2(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out);
    } else if (L == 6) {
        gg_collocation_L6_deriv2(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out);
    } else {
        exit(0);
    }
}
void gg_collocation_deriv3(int L, const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out, double* PRAGMA_RESTRICT phi_xxx_out, double* PRAGMA_RESTRICT phi_xxy_out, double* PRAGMA_RESTRICT phi_xxz_out, double* PRAGMA_RESTRICT phi_xyy_out, double* PRAGMA_RESTRICT phi_xyz_out, double* PRAGMA_RESTRICT phi_xzz_out, double* PRAGMA_RESTRICT phi_yyy_out, double* PRAGMA_RESTRICT phi_yyz_out, double* PRAGMA_RESTRICT phi_yzz_out, double* PRAGMA_RESTRICT phi_zzz_out) {
    // Chooses the correct function for a given L
    if (L == 0) {
        gg_collocation_L0_deriv3(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out, phi_xxx_out, phi_xxy_out, phi_xxz_out, phi_xyy_out, phi_xyz_out, phi_xzz_out, phi_yyy_out, phi_yyz_out, phi_yzz_out, phi_zzz_out);
    } else if (L == 1) {
        gg_collocation_L1_deriv3(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out, phi_xxx_out, phi_xxy_out, phi_xxz_out, phi_xyy_out, phi_xyz_out, phi_xzz_out, phi_yyy_out, phi_yyz_out, phi_yzz_out, phi_zzz_out);
    } else if (L == 2) {
        gg_collocation_L2_deriv3(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out, phi_xxx_out, phi_xxy_out, phi_xxz_out, phi_xyy_out, phi_xyz_out, phi_xzz_out, phi_yyy_out, phi_yyz_out, phi_yzz_out, phi_zzz_out);
    } else if (L == 3) {
        gg_collocation_L3_deriv3(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out, phi_xxx_out, phi_xxy_out, phi_xxz_out, phi_xyy_out, phi_xyz_out, phi_xzz_out, phi_yyy_out, phi_yyz_out, phi_yzz_out, phi_zzz_out);
    } else if (L == 4) {
        gg_collocation_L4_deriv3(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out, phi_xxx_out, phi_xxy_out, phi_xxz_out, phi_xyy_out, phi_xyz_out, phi_xzz_out, phi_yyy_out, phi_yyz_out, phi_yzz_out, phi_zzz_out);
    } else if (L == 5) {
        gg_collocation_L5_deriv3(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out, phi_xxx_out, phi_xxy_out, phi_xxz_out, phi_xyy_out, phi_xyz_out, phi_xzz_out, phi_yyy_out, phi_yyz_out, phi_yzz_out, phi_zzz_out);
    } else if (L == 6) {
        gg_collocation_L6_deriv3(npoints, xyz, xyz_stride, nprim, coeffs, exponents, center, order, phi_out, phi_x_out, phi_y_out, phi_z_out, phi_xx_out, phi_xy_out, phi_xz_out, phi_yy_out, phi_yz_out, phi_zz_out, phi_xxx_out, phi_xxy_out, phi_xxz_out, phi_xyy_out, phi_xyz_out, phi_xzz_out, phi_yyy_out, phi_yyz_out, phi_yzz_out, phi_zzz_out);
    } else {
        exit(0);
    }
}