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



void gg_collocation_L0_deriv3(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out, double* PRAGMA_RESTRICT phi_xxx_out, double* PRAGMA_RESTRICT phi_xxy_out, double* PRAGMA_RESTRICT phi_xxz_out, double* PRAGMA_RESTRICT phi_xyy_out, double* PRAGMA_RESTRICT phi_xyz_out, double* PRAGMA_RESTRICT phi_xzz_out, double* PRAGMA_RESTRICT phi_yyy_out, double* PRAGMA_RESTRICT phi_yyz_out, double* PRAGMA_RESTRICT phi_yzz_out, double* PRAGMA_RESTRICT phi_zzz_out) {

    // Sizing
    unsigned long nblocks = npoints / 32;
    nblocks += (npoints % 32) ? 1 : 0;
    const unsigned long ncart = 1;
    const unsigned long nspherical = 1;
    unsigned long nout;

    if ((order == GG_SPHERICAL_CCA) || (order == GG_SPHERICAL_GAUSSIAN)) {
        nout = nspherical;
        } else {
        nout = ncart;
    }

    // Allocate S temporaries, single block to stay on cache
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 288 * sizeof(double));
    double* PRAGMA_RESTRICT xc = cache_data + 0;
    ASSUME_ALIGNED(xc, 64);
    double* PRAGMA_RESTRICT yc = cache_data + 32;
    ASSUME_ALIGNED(yc, 64);
    double* PRAGMA_RESTRICT zc = cache_data + 64;
    ASSUME_ALIGNED(zc, 64);
    double* PRAGMA_RESTRICT R2 = cache_data + 96;
    ASSUME_ALIGNED(R2, 64);
    double* PRAGMA_RESTRICT S0 = cache_data + 128;
    ASSUME_ALIGNED(S0, 64);
    double* PRAGMA_RESTRICT tmp1 = cache_data + 160;
    ASSUME_ALIGNED(tmp1, 64);
    double* PRAGMA_RESTRICT S1 = cache_data + 192;
    ASSUME_ALIGNED(S1, 64);
    double* PRAGMA_RESTRICT S2 = cache_data + 224;
    ASSUME_ALIGNED(S2, 64);
    double* PRAGMA_RESTRICT S3 = cache_data + 256;
    ASSUME_ALIGNED(S3, 64);

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));
    double* PRAGMA_RESTRICT expn2 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_tmp, 64);
    double* PRAGMA_RESTRICT phi_x_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_x_tmp, 64);
    double* PRAGMA_RESTRICT phi_y_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_y_tmp, 64);
    double* PRAGMA_RESTRICT phi_z_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_z_tmp, 64);
    double* PRAGMA_RESTRICT phi_xx_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_xx_tmp, 64);
    double* PRAGMA_RESTRICT phi_xy_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_xy_tmp, 64);
    double* PRAGMA_RESTRICT phi_xz_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_xz_tmp, 64);
    double* PRAGMA_RESTRICT phi_yy_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_yy_tmp, 64);
    double* PRAGMA_RESTRICT phi_yz_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_yz_tmp, 64);
    double* PRAGMA_RESTRICT phi_zz_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_zz_tmp, 64);
    double* PRAGMA_RESTRICT phi_xxx_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_xxx_tmp, 64);
    double* PRAGMA_RESTRICT phi_xxy_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_xxy_tmp, 64);
    double* PRAGMA_RESTRICT phi_xxz_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_xxz_tmp, 64);
    double* PRAGMA_RESTRICT phi_xyy_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_xyy_tmp, 64);
    double* PRAGMA_RESTRICT phi_xyz_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_xyz_tmp, 64);
    double* PRAGMA_RESTRICT phi_xzz_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_xzz_tmp, 64);
    double* PRAGMA_RESTRICT phi_yyy_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_yyy_tmp, 64);
    double* PRAGMA_RESTRICT phi_yyz_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_yyz_tmp, 64);
    double* PRAGMA_RESTRICT phi_yzz_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_yzz_tmp, 64);
    double* PRAGMA_RESTRICT phi_zzz_tmp = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(phi_zzz_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;
    double AXX, AXY, AXZ, AYY, AYZ, AZZ;
    double AXXX, XXY, XXZ, XYY, XYZ, XZZ, YYY, YYZ, YZZ, ZZZ;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
        expn2[i] = -2.0 * exponents[i];
    }

    // Start outer block loop
    for (unsigned long block = 0; block < nblocks; block++) {


        // Copy data into inner temps
        const unsigned long start = block * 32;
        const unsigned long remain = ((start + 32) > npoints) ? (npoints - start) : 32;

        // Handle non-AM dependant temps
        if (xyz_stride == 1) {
            const double* PRAGMA_RESTRICT x = xyz + start;
            const double* PRAGMA_RESTRICT y = xyz + npoints + start;
            const double* PRAGMA_RESTRICT z = xyz + 2 * npoints + start;
            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = x[i] - center_x;
                yc[i] = y[i] - center_y;
                zc[i] = z[i] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
            } else {
            unsigned int start_shift = start * xyz_stride;

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = xyz[start_shift + i * xyz_stride] - center_x;
                yc[i] = xyz[start_shift + i * xyz_stride + 1] - center_y;
                zc[i] = xyz[start_shift + i * xyz_stride + 2] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];
            const double alpha_n2 = expn2[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
                const double T2 = alpha_n2 * T1;
                S1[i] += T2;
                const double T3 = alpha_n2 * T2;
                S2[i] += T3;
                const double T4 = alpha_n2 * T3;
                S3[i] += T4;
            }

        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            // Gaussian derivs (gradients)
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];

            // Gaussians derivs (Hessians)
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];

            // Gaussians 3rd derivs)
            const double SXXX = S3[i] * xc[i] * xc[i] * xc[i] + 3 * xc[i] * S2[i];
            const double SXXY = S3[i] * xc[i] * xc[i] * yc[i] + yc[i] * S2[i];
            const double SXXZ = S3[i] * xc[i] * xc[i] * zc[i] + zc[i] * S2[i];
            const double SXYY = S3[i] * xc[i] * yc[i] * yc[i] + xc[i] * S2[i];
            const double SXYZ = S3[i] * xc[i] * yc[i] * zc[i];
            const double SXZZ = S3[i] * xc[i] * zc[i] * zc[i] + xc[i] * S2[i];
            const double SYYY = S3[i] * yc[i] * yc[i] * yc[i] + 3 * yc[i] * S2[i];
            const double SYYZ = S3[i] * yc[i] * yc[i] * zc[i] + zc[i] * S2[i];
            const double SYZZ = S3[i] * yc[i] * zc[i] * zc[i] + yc[i] * S2[i];
            const double SZZZ = S3[i] * zc[i] * zc[i] * zc[i] + 3 * zc[i] * S2[i];
            phi_out[start + i] = S0[i];

            // Gradient AM=0 Component=0
            phi_x_out[start + i] = SX;
            phi_y_out[start + i] = SY;
            phi_z_out[start + i] = SZ;

            // Hessian AM=0 Component=0
            phi_xx_out[start + i] = SXX;
            phi_yy_out[start + i] = SYY;
            phi_zz_out[start + i] = SZZ;
            phi_xy_out[start + i] = SXY;
            phi_xz_out[start + i] = SXZ;
            phi_yz_out[start + i] = SYZ;

            // Der3 AM=0 Component=0
            phi_xxx_out[start + i] = SXXX;
            phi_xxy_out[start + i] = SXXY;
            phi_xxz_out[start + i] = SXXZ;
            phi_xyy_out[start + i] = SXYY;
            phi_xyz_out[start + i] = SXYZ;
            phi_xzz_out[start + i] = SXZZ;
            phi_yyy_out[start + i] = SYYY;
            phi_yyz_out[start + i] = SYYZ;
            phi_yzz_out[start + i] = SYZZ;
            phi_zzz_out[start + i] = SZZZ;
        }
    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);
    ALIGNED_FREE(expn2);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);
    ALIGNED_FREE(phi_x_tmp);
    ALIGNED_FREE(phi_y_tmp);
    ALIGNED_FREE(phi_z_tmp);
    ALIGNED_FREE(phi_xx_tmp);
    ALIGNED_FREE(phi_xy_tmp);
    ALIGNED_FREE(phi_xz_tmp);
    ALIGNED_FREE(phi_yy_tmp);
    ALIGNED_FREE(phi_yz_tmp);
    ALIGNED_FREE(phi_zz_tmp);
    ALIGNED_FREE(phi_xxx_tmp);
    ALIGNED_FREE(phi_xxy_tmp);
    ALIGNED_FREE(phi_xxz_tmp);
    ALIGNED_FREE(phi_xyy_tmp);
    ALIGNED_FREE(phi_xyz_tmp);
    ALIGNED_FREE(phi_xzz_tmp);
    ALIGNED_FREE(phi_yyy_tmp);
    ALIGNED_FREE(phi_yyz_tmp);
    ALIGNED_FREE(phi_yzz_tmp);
    ALIGNED_FREE(phi_zzz_tmp);

}

void gg_collocation_L1_deriv3(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out, double* PRAGMA_RESTRICT phi_xxx_out, double* PRAGMA_RESTRICT phi_xxy_out, double* PRAGMA_RESTRICT phi_xxz_out, double* PRAGMA_RESTRICT phi_xyy_out, double* PRAGMA_RESTRICT phi_xyz_out, double* PRAGMA_RESTRICT phi_xzz_out, double* PRAGMA_RESTRICT phi_yyy_out, double* PRAGMA_RESTRICT phi_yyz_out, double* PRAGMA_RESTRICT phi_yzz_out, double* PRAGMA_RESTRICT phi_zzz_out) {

    // Sizing
    unsigned long nblocks = npoints / 32;
    nblocks += (npoints % 32) ? 1 : 0;
    const unsigned long ncart = 3;
    const unsigned long nspherical = 3;
    unsigned long nout;

    if ((order == GG_SPHERICAL_CCA) || (order == GG_SPHERICAL_GAUSSIAN)) {
        nout = nspherical;
        } else {
        nout = ncart;
    }

    // Allocate S temporaries, single block to stay on cache
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 288 * sizeof(double));
    double* PRAGMA_RESTRICT xc = cache_data + 0;
    ASSUME_ALIGNED(xc, 64);
    double* PRAGMA_RESTRICT yc = cache_data + 32;
    ASSUME_ALIGNED(yc, 64);
    double* PRAGMA_RESTRICT zc = cache_data + 64;
    ASSUME_ALIGNED(zc, 64);
    double* PRAGMA_RESTRICT R2 = cache_data + 96;
    ASSUME_ALIGNED(R2, 64);
    double* PRAGMA_RESTRICT S0 = cache_data + 128;
    ASSUME_ALIGNED(S0, 64);
    double* PRAGMA_RESTRICT tmp1 = cache_data + 160;
    ASSUME_ALIGNED(tmp1, 64);
    double* PRAGMA_RESTRICT S1 = cache_data + 192;
    ASSUME_ALIGNED(S1, 64);
    double* PRAGMA_RESTRICT S2 = cache_data + 224;
    ASSUME_ALIGNED(S2, 64);
    double* PRAGMA_RESTRICT S3 = cache_data + 256;
    ASSUME_ALIGNED(S3, 64);

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));
    double* PRAGMA_RESTRICT expn2 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_tmp, 64);
    double* PRAGMA_RESTRICT phi_x_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_x_tmp, 64);
    double* PRAGMA_RESTRICT phi_y_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_y_tmp, 64);
    double* PRAGMA_RESTRICT phi_z_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_z_tmp, 64);
    double* PRAGMA_RESTRICT phi_xx_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_xx_tmp, 64);
    double* PRAGMA_RESTRICT phi_xy_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_xy_tmp, 64);
    double* PRAGMA_RESTRICT phi_xz_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_xz_tmp, 64);
    double* PRAGMA_RESTRICT phi_yy_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_yy_tmp, 64);
    double* PRAGMA_RESTRICT phi_yz_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_yz_tmp, 64);
    double* PRAGMA_RESTRICT phi_zz_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_zz_tmp, 64);
    double* PRAGMA_RESTRICT phi_xxx_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_xxx_tmp, 64);
    double* PRAGMA_RESTRICT phi_xxy_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_xxy_tmp, 64);
    double* PRAGMA_RESTRICT phi_xxz_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_xxz_tmp, 64);
    double* PRAGMA_RESTRICT phi_xyy_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_xyy_tmp, 64);
    double* PRAGMA_RESTRICT phi_xyz_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_xyz_tmp, 64);
    double* PRAGMA_RESTRICT phi_xzz_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_xzz_tmp, 64);
    double* PRAGMA_RESTRICT phi_yyy_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_yyy_tmp, 64);
    double* PRAGMA_RESTRICT phi_yyz_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_yyz_tmp, 64);
    double* PRAGMA_RESTRICT phi_yzz_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_yzz_tmp, 64);
    double* PRAGMA_RESTRICT phi_zzz_tmp = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(phi_zzz_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;
    double AXX, AXY, AXZ, AYY, AYZ, AZZ;
    double AXXX, XXY, XXZ, XYY, XYZ, XZZ, YYY, YYZ, YZZ, ZZZ;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
        expn2[i] = -2.0 * exponents[i];
    }

    // Start outer block loop
    for (unsigned long block = 0; block < nblocks; block++) {


        // Copy data into inner temps
        const unsigned long start = block * 32;
        const unsigned long remain = ((start + 32) > npoints) ? (npoints - start) : 32;

        // Handle non-AM dependant temps
        if (xyz_stride == 1) {
            const double* PRAGMA_RESTRICT x = xyz + start;
            const double* PRAGMA_RESTRICT y = xyz + npoints + start;
            const double* PRAGMA_RESTRICT z = xyz + 2 * npoints + start;
            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = x[i] - center_x;
                yc[i] = y[i] - center_y;
                zc[i] = z[i] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
            } else {
            unsigned int start_shift = start * xyz_stride;

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = xyz[start_shift + i * xyz_stride] - center_x;
                yc[i] = xyz[start_shift + i * xyz_stride + 1] - center_y;
                zc[i] = xyz[start_shift + i * xyz_stride + 2] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];
            const double alpha_n2 = expn2[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
                const double T2 = alpha_n2 * T1;
                S1[i] += T2;
                const double T3 = alpha_n2 * T2;
                S2[i] += T3;
                const double T4 = alpha_n2 * T3;
                S3[i] += T4;
            }

        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            // Gaussian derivs (gradients)
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];

            // Gaussians derivs (Hessians)
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];

            // Gaussians 3rd derivs)
            const double SXXX = S3[i] * xc[i] * xc[i] * xc[i] + 3 * xc[i] * S2[i];
            const double SXXY = S3[i] * xc[i] * xc[i] * yc[i] + yc[i] * S2[i];
            const double SXXZ = S3[i] * xc[i] * xc[i] * zc[i] + zc[i] * S2[i];
            const double SXYY = S3[i] * xc[i] * yc[i] * yc[i] + xc[i] * S2[i];
            const double SXYZ = S3[i] * xc[i] * yc[i] * zc[i];
            const double SXZZ = S3[i] * xc[i] * zc[i] * zc[i] + xc[i] * S2[i];
            const double SYYY = S3[i] * yc[i] * yc[i] * yc[i] + 3 * yc[i] * S2[i];
            const double SYYZ = S3[i] * yc[i] * yc[i] * zc[i] + zc[i] * S2[i];
            const double SYZZ = S3[i] * yc[i] * zc[i] * zc[i] + yc[i] * S2[i];
            const double SZZZ = S3[i] * zc[i] * zc[i] * zc[i] + 3 * zc[i] * S2[i];

            // Density AM=1 Component=X
            phi_tmp[i] = S0[i] * xc[i];

            // Gradient AM=1 Component=X
            phi_x_tmp[i] = SX * xc[i];
            phi_y_tmp[i] = SY * xc[i];
            phi_z_tmp[i] = SZ * xc[i];
            phi_x_tmp[i] += S0[i];

            // Hessian AM=1 Component=X
            phi_xx_tmp[i] = SXX * xc[i];
            phi_xx_tmp[i] += SX;
            phi_xx_tmp[i] += SX;
            phi_yy_tmp[i] = SYY * xc[i];
            phi_zz_tmp[i] = SZZ * xc[i];
            phi_xy_tmp[i] = SXY * xc[i];
            phi_xy_tmp[i] += SY;
            phi_xz_tmp[i] = SXZ * xc[i];
            phi_xz_tmp[i] += SZ;
            phi_yz_tmp[i] = SYZ * xc[i];
            phi_xyz_tmp[i] = SXYZ * xc[i];
            phi_xyz_tmp[i] += SYZ;
            phi_xxy_tmp[i] = SXXY * xc[i];
            phi_xxy_tmp[i] += 2.0 * SXY;
            phi_xxz_tmp[i] = SXXZ * xc[i];
            phi_xxz_tmp[i] += 2.0 * SXZ;
            phi_xyy_tmp[i] = SXYY * xc[i];
            phi_xyy_tmp[i] += SYY;
            phi_xzz_tmp[i] = SXZZ * xc[i];
            phi_xzz_tmp[i] += SZZ;
            phi_yyz_tmp[i] = SYYZ * xc[i];
            phi_yzz_tmp[i] = SYZZ * xc[i];
            phi_xxx_tmp[i] = SXXX * xc[i];
            phi_xxx_tmp[i] += 3.0 * SXX;
            phi_yyy_tmp[i] = SYYY * xc[i];
            phi_zzz_tmp[i] = SZZZ * xc[i];


            // Density AM=1 Component=Y
            phi_tmp[32 + i] = S0[i] * yc[i];

            // Gradient AM=1 Component=Y
            phi_x_tmp[32 + i] = SX * yc[i];
            phi_y_tmp[32 + i] = SY * yc[i];
            phi_z_tmp[32 + i] = SZ * yc[i];
            phi_y_tmp[32 + i] += S0[i];

            // Hessian AM=1 Component=Y
            phi_xx_tmp[32 + i] = SXX * yc[i];
            phi_yy_tmp[32 + i] = SYY * yc[i];
            phi_yy_tmp[32 + i] += SY;
            phi_yy_tmp[32 + i] += SY;
            phi_zz_tmp[32 + i] = SZZ * yc[i];
            phi_xy_tmp[32 + i] = SXY * yc[i];
            phi_xy_tmp[32 + i] += SX;
            phi_xz_tmp[32 + i] = SXZ * yc[i];
            phi_yz_tmp[32 + i] = SYZ * yc[i];
            phi_yz_tmp[32 + i] += SZ;
            phi_xyz_tmp[32 + i] = SXYZ * yc[i];
            phi_xyz_tmp[32 + i] += SXZ;
            phi_xxy_tmp[32 + i] = SXXY * yc[i];
            phi_xxy_tmp[32 + i] += SXX;
            phi_xxz_tmp[32 + i] = SXXZ * yc[i];
            phi_xyy_tmp[32 + i] = SXYY * yc[i];
            phi_xyy_tmp[32 + i] += 2.0 * SXY;
            phi_xzz_tmp[32 + i] = SXZZ * yc[i];
            phi_yyz_tmp[32 + i] = SYYZ * yc[i];
            phi_yyz_tmp[32 + i] += 2.0 * SYZ;
            phi_yzz_tmp[32 + i] = SYZZ * yc[i];
            phi_yzz_tmp[32 + i] += SZZ;
            phi_xxx_tmp[32 + i] = SXXX * yc[i];
            phi_yyy_tmp[32 + i] = SYYY * yc[i];
            phi_yyy_tmp[32 + i] += 3.0 * SYY;
            phi_zzz_tmp[32 + i] = SZZZ * yc[i];


            // Density AM=1 Component=Z
            phi_tmp[64 + i] = S0[i] * zc[i];

            // Gradient AM=1 Component=Z
            phi_x_tmp[64 + i] = SX * zc[i];
            phi_y_tmp[64 + i] = SY * zc[i];
            phi_z_tmp[64 + i] = SZ * zc[i];
            phi_z_tmp[64 + i] += S0[i];

            // Hessian AM=1 Component=Z
            phi_xx_tmp[64 + i] = SXX * zc[i];
            phi_yy_tmp[64 + i] = SYY * zc[i];
            phi_zz_tmp[64 + i] = SZZ * zc[i];
            phi_zz_tmp[64 + i] += SZ;
            phi_zz_tmp[64 + i] += SZ;
            phi_xy_tmp[64 + i] = SXY * zc[i];
            phi_xz_tmp[64 + i] = SXZ * zc[i];
            phi_xz_tmp[64 + i] += SX;
            phi_yz_tmp[64 + i] = SYZ * zc[i];
            phi_yz_tmp[64 + i] += SY;
            phi_xyz_tmp[64 + i] = SXYZ * zc[i];
            phi_xyz_tmp[64 + i] += SXY;
            phi_xxy_tmp[64 + i] = SXXY * zc[i];
            phi_xxz_tmp[64 + i] = SXXZ * zc[i];
            phi_xxz_tmp[64 + i] += SXX;
            phi_xyy_tmp[64 + i] = SXYY * zc[i];
            phi_xzz_tmp[64 + i] = SXZZ * zc[i];
            phi_xzz_tmp[64 + i] += 2.0 * SXZ;
            phi_yyz_tmp[64 + i] = SYYZ * zc[i];
            phi_yyz_tmp[64 + i] += SYY;
            phi_yzz_tmp[64 + i] = SYZZ * zc[i];
            phi_yzz_tmp[64 + i] += 2.0 * SYZ;
            phi_xxx_tmp[64 + i] = SXXX * zc[i];
            phi_yyy_tmp[64 + i] = SYYY * zc[i];
            phi_zzz_tmp[64 + i] = SZZZ * zc[i];
            phi_zzz_tmp[64 + i] += 3.0 * SZZ;


        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L1(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_to_spherical_L1(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);

            // Hessian, transform data to outer temps
            gg_cca_cart_to_spherical_L1(remain, phi_xx_tmp, 32, (phi_xx_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_xy_tmp, 32, (phi_xy_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_xz_tmp, 32, (phi_xz_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_yy_tmp, 32, (phi_yy_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_yz_tmp, 32, (phi_yz_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_zz_tmp, 32, (phi_zz_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_xxx_tmp, 32, (phi_xxx_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_xxy_tmp, 32, (phi_xxy_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_xxz_tmp, 32, (phi_xxz_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_xyy_tmp, 32, (phi_xyy_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_xyz_tmp, 32, (phi_xyz_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_xzz_tmp, 32, (phi_xzz_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_yyy_tmp, 32, (phi_yyy_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_yyz_tmp, 32, (phi_yyz_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_yzz_tmp, 32, (phi_yzz_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_zzz_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L1(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_gaussian_cart_to_spherical_L1(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);

            // Hessian, transform data to outer temps
            gg_gaussian_cart_to_spherical_L1(remain, phi_xx_tmp, 32, (phi_xx_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_xy_tmp, 32, (phi_xy_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_xz_tmp, 32, (phi_xz_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_yy_tmp, 32, (phi_yy_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_yz_tmp, 32, (phi_yz_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_zz_tmp, 32, (phi_zz_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_xxx_tmp, 32, (phi_xxx_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_xxy_tmp, 32, (phi_xxy_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_xxz_tmp, 32, (phi_xxz_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_xyy_tmp, 32, (phi_xyy_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_xyz_tmp, 32, (phi_xyz_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_xzz_tmp, 32, (phi_xzz_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_yyy_tmp, 32, (phi_yyy_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_yyz_tmp, 32, (phi_yyz_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_yzz_tmp, 32, (phi_yzz_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_zzz_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L1(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_copy_L1(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);

            // Hessian, transform data to outer temps
            gg_cca_cart_copy_L1(remain, phi_xx_tmp, 32, (phi_xx_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_xy_tmp, 32, (phi_xy_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_xz_tmp, 32, (phi_xz_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_yy_tmp, 32, (phi_yy_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_yz_tmp, 32, (phi_yz_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_zz_tmp, 32, (phi_zz_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_xxx_tmp, 32, (phi_xxx_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_xxy_tmp, 32, (phi_xxy_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_xxz_tmp, 32, (phi_xxz_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_xyy_tmp, 32, (phi_xyy_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_xyz_tmp, 32, (phi_xyz_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_xzz_tmp, 32, (phi_xzz_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_yyy_tmp, 32, (phi_yyy_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_yyz_tmp, 32, (phi_yyz_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_yzz_tmp, 32, (phi_yzz_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_zzz_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L1(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_molden_cart_copy_L1(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);

            // Hessian, transform data to outer temps
            gg_molden_cart_copy_L1(remain, phi_xx_tmp, 32, (phi_xx_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_xy_tmp, 32, (phi_xy_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_xz_tmp, 32, (phi_xz_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_yy_tmp, 32, (phi_yy_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_yz_tmp, 32, (phi_yz_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_zz_tmp, 32, (phi_zz_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_xxx_tmp, 32, (phi_xxx_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_xxy_tmp, 32, (phi_xxy_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_xxz_tmp, 32, (phi_xxz_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_xyy_tmp, 32, (phi_xyy_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_xyz_tmp, 32, (phi_xyz_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_xzz_tmp, 32, (phi_xzz_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_yyy_tmp, 32, (phi_yyy_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_yyz_tmp, 32, (phi_yyz_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_yzz_tmp, 32, (phi_yzz_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_zzz_tmp, 32, (phi_zzz_out + start), npoints);
        }

    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);
    ALIGNED_FREE(expn2);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);
    ALIGNED_FREE(phi_x_tmp);
    ALIGNED_FREE(phi_y_tmp);
    ALIGNED_FREE(phi_z_tmp);
    ALIGNED_FREE(phi_xx_tmp);
    ALIGNED_FREE(phi_xy_tmp);
    ALIGNED_FREE(phi_xz_tmp);
    ALIGNED_FREE(phi_yy_tmp);
    ALIGNED_FREE(phi_yz_tmp);
    ALIGNED_FREE(phi_zz_tmp);
    ALIGNED_FREE(phi_xxx_tmp);
    ALIGNED_FREE(phi_xxy_tmp);
    ALIGNED_FREE(phi_xxz_tmp);
    ALIGNED_FREE(phi_xyy_tmp);
    ALIGNED_FREE(phi_xyz_tmp);
    ALIGNED_FREE(phi_xzz_tmp);
    ALIGNED_FREE(phi_yyy_tmp);
    ALIGNED_FREE(phi_yyz_tmp);
    ALIGNED_FREE(phi_yzz_tmp);
    ALIGNED_FREE(phi_zzz_tmp);

}

void gg_collocation_L2_deriv3(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out, double* PRAGMA_RESTRICT phi_xxx_out, double* PRAGMA_RESTRICT phi_xxy_out, double* PRAGMA_RESTRICT phi_xxz_out, double* PRAGMA_RESTRICT phi_xyy_out, double* PRAGMA_RESTRICT phi_xyz_out, double* PRAGMA_RESTRICT phi_xzz_out, double* PRAGMA_RESTRICT phi_yyy_out, double* PRAGMA_RESTRICT phi_yyz_out, double* PRAGMA_RESTRICT phi_yzz_out, double* PRAGMA_RESTRICT phi_zzz_out) {

    // Sizing
    unsigned long nblocks = npoints / 32;
    nblocks += (npoints % 32) ? 1 : 0;
    const unsigned long ncart = 6;
    const unsigned long nspherical = 5;
    unsigned long nout;

    if ((order == GG_SPHERICAL_CCA) || (order == GG_SPHERICAL_GAUSSIAN)) {
        nout = nspherical;
        } else {
        nout = ncart;
    }

    // Allocate S temporaries, single block to stay on cache
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 288 * sizeof(double));
    double* PRAGMA_RESTRICT xc = cache_data + 0;
    ASSUME_ALIGNED(xc, 64);
    double* PRAGMA_RESTRICT yc = cache_data + 32;
    ASSUME_ALIGNED(yc, 64);
    double* PRAGMA_RESTRICT zc = cache_data + 64;
    ASSUME_ALIGNED(zc, 64);
    double* PRAGMA_RESTRICT R2 = cache_data + 96;
    ASSUME_ALIGNED(R2, 64);
    double* PRAGMA_RESTRICT S0 = cache_data + 128;
    ASSUME_ALIGNED(S0, 64);
    double* PRAGMA_RESTRICT tmp1 = cache_data + 160;
    ASSUME_ALIGNED(tmp1, 64);
    double* PRAGMA_RESTRICT S1 = cache_data + 192;
    ASSUME_ALIGNED(S1, 64);
    double* PRAGMA_RESTRICT S2 = cache_data + 224;
    ASSUME_ALIGNED(S2, 64);
    double* PRAGMA_RESTRICT S3 = cache_data + 256;
    ASSUME_ALIGNED(S3, 64);

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));
    double* PRAGMA_RESTRICT expn2 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));

    // Allocate power temporaries
    double* PRAGMA_RESTRICT xc_pow = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(xc_pow, 64);
    double* PRAGMA_RESTRICT yc_pow = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(yc_pow, 64);
    double* PRAGMA_RESTRICT zc_pow = (double*)ALIGNED_MALLOC(64, 32 * sizeof(double));
    ASSUME_ALIGNED(zc_pow, 64);

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, 192 * sizeof(double));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;
    double AXX, AXY, AXZ, AYY, AYZ, AZZ;
    double AXXX, XXY, XXZ, XYY, XYZ, XZZ, YYY, YYZ, YZZ, ZZZ;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
        expn2[i] = -2.0 * exponents[i];
    }

    // Start outer block loop
    for (unsigned long block = 0; block < nblocks; block++) {


        // Copy data into inner temps
        const unsigned long start = block * 32;
        const unsigned long remain = ((start + 32) > npoints) ? (npoints - start) : 32;

        // Handle non-AM dependant temps
        if (xyz_stride == 1) {
            const double* PRAGMA_RESTRICT x = xyz + start;
            const double* PRAGMA_RESTRICT y = xyz + npoints + start;
            const double* PRAGMA_RESTRICT z = xyz + 2 * npoints + start;
            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = x[i] - center_x;
                yc[i] = y[i] - center_y;
                zc[i] = z[i] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
            } else {
            unsigned int start_shift = start * xyz_stride;

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = xyz[start_shift + i * xyz_stride] - center_x;
                yc[i] = xyz[start_shift + i * xyz_stride + 1] - center_y;
                zc[i] = xyz[start_shift + i * xyz_stride + 2] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];
            const double alpha_n2 = expn2[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
                const double T2 = alpha_n2 * T1;
                S1[i] += T2;
                const double T3 = alpha_n2 * T2;
                S2[i] += T3;
                const double T4 = alpha_n2 * T3;
                S3[i] += T4;
            }

        }

        // Build powers
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {

            // Cartesian derivs
            xc_pow[i] = xc[i] * xc[i];
            yc_pow[i] = yc[i] * yc[i];
            zc_pow[i] = zc[i] * zc[i];

        }
        // Combine A blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {

            phi_tmp[i] = xc_pow[i] * S0[i];
            phi_tmp[32 + i] = xc[i] * yc[i] * S0[i];
            phi_tmp[64 + i] = xc[i] * zc[i] * S0[i];
            phi_tmp[96 + i] = yc_pow[i] * S0[i];
            phi_tmp[128 + i] = yc[i] * zc[i] * S0[i];
            phi_tmp[160 + i] = zc_pow[i] * S0[i];
        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_out + start), npoints);
        }

        // Combine X blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];

            phi_tmp[i] = xc_pow[i] * SX;
            phi_tmp[i] += 2.0 * xc[i] * S0[i];

            phi_tmp[32 + i] = xc[i] * yc[i] * SX;
            phi_tmp[32 + i] += yc[i] * S0[i];

            phi_tmp[64 + i] = xc[i] * zc[i] * SX;
            phi_tmp[64 + i] += zc[i] * S0[i];

            phi_tmp[96 + i] = yc_pow[i] * SX;

            phi_tmp[128 + i] = yc[i] * zc[i] * SX;

            phi_tmp[160 + i] = zc_pow[i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_x_out + start), npoints);
        }

        // Combine Y blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];

            phi_tmp[i] = xc_pow[i] * SY;

            phi_tmp[32 + i] = xc[i] * yc[i] * SY;
            phi_tmp[32 + i] += xc[i] * S0[i];

            phi_tmp[64 + i] = xc[i] * zc[i] * SY;

            phi_tmp[96 + i] = yc_pow[i] * SY;
            phi_tmp[96 + i] += 2.0 * yc[i] * S0[i];

            phi_tmp[128 + i] = yc[i] * zc[i] * SY;
            phi_tmp[128 + i] += zc[i] * S0[i];

            phi_tmp[160 + i] = zc_pow[i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_y_out + start), npoints);
        }

        // Combine Z blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];

            phi_tmp[i] = xc_pow[i] * SZ;

            phi_tmp[32 + i] = xc[i] * yc[i] * SZ;

            phi_tmp[64 + i] = xc[i] * zc[i] * SZ;
            phi_tmp[64 + i] += xc[i] * S0[i];

            phi_tmp[96 + i] = yc_pow[i] * SZ;

            phi_tmp[128 + i] = yc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += yc[i] * S0[i];

            phi_tmp[160 + i] = zc_pow[i] * SZ;
            phi_tmp[160 + i] += 2.0 * zc[i] * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_z_out + start), npoints);
        }

        // Combine XX blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];

            phi_tmp[i] = xc_pow[i] * SXX;
            phi_tmp[i] += 4.0 * xc[i] * SX;
            phi_tmp[i] += 2.0 * S0[i];

            phi_tmp[32 + i] = xc[i] * yc[i] * SXX;
            phi_tmp[32 + i] += 2.0 * yc[i] * SX;

            phi_tmp[64 + i] = xc[i] * zc[i] * SXX;
            phi_tmp[64 + i] += 2.0 * zc[i] * SX;

            phi_tmp[96 + i] = yc_pow[i] * SXX;

            phi_tmp[128 + i] = yc[i] * zc[i] * SXX;

            phi_tmp[160 + i] = zc_pow[i] * SXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
        }

        // Combine XY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];

            phi_tmp[i] = xc_pow[i] * SXY;
            phi_tmp[i] += 2.0 * xc[i] * SY;

            phi_tmp[32 + i] = xc[i] * yc[i] * SXY;
            phi_tmp[32 + i] += xc[i] * SX;
            phi_tmp[32 + i] += yc[i] * SY;
            phi_tmp[32 + i] +=  1 * S0[i];

            phi_tmp[64 + i] = xc[i] * zc[i] * SXY;
            phi_tmp[64 + i] += zc[i] * SY;

            phi_tmp[96 + i] = yc_pow[i] * SXY;
            phi_tmp[96 + i] += 2.0 * yc[i] * SX;

            phi_tmp[128 + i] = yc[i] * zc[i] * SXY;
            phi_tmp[128 + i] += zc[i] * SX;

            phi_tmp[160 + i] = zc_pow[i] * SXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
        }

        // Combine XZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];

            phi_tmp[i] = xc_pow[i] * SXZ;
            phi_tmp[i] += 2.0 * xc[i] * SZ;

            phi_tmp[32 + i] = xc[i] * yc[i] * SXZ;
            phi_tmp[32 + i] += yc[i] * SZ;

            phi_tmp[64 + i] = xc[i] * zc[i] * SXZ;
            phi_tmp[64 + i] += xc[i] * SX;
            phi_tmp[64 + i] += zc[i] * SZ;
            phi_tmp[64 + i] +=  1 * S0[i];

            phi_tmp[96 + i] = yc_pow[i] * SXZ;

            phi_tmp[128 + i] = yc[i] * zc[i] * SXZ;
            phi_tmp[128 + i] += yc[i] * SX;

            phi_tmp[160 + i] = zc_pow[i] * SXZ;
            phi_tmp[160 + i] += 2.0 * zc[i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
        }

        // Combine YY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];

            phi_tmp[i] = xc_pow[i] * SYY;

            phi_tmp[32 + i] = xc[i] * yc[i] * SYY;
            phi_tmp[32 + i] += 2.0 * xc[i] * SY;

            phi_tmp[64 + i] = xc[i] * zc[i] * SYY;

            phi_tmp[96 + i] = yc_pow[i] * SYY;
            phi_tmp[96 + i] += 4.0 * yc[i] * SY;
            phi_tmp[96 + i] += 2.0 * S0[i];

            phi_tmp[128 + i] = yc[i] * zc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * zc[i] * SY;

            phi_tmp[160 + i] = zc_pow[i] * SYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
        }

        // Combine YZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];

            phi_tmp[i] = xc_pow[i] * SYZ;

            phi_tmp[32 + i] = xc[i] * yc[i] * SYZ;
            phi_tmp[32 + i] += xc[i] * SZ;

            phi_tmp[64 + i] = xc[i] * zc[i] * SYZ;
            phi_tmp[64 + i] += xc[i] * SY;

            phi_tmp[96 + i] = yc_pow[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * yc[i] * SZ;

            phi_tmp[128 + i] = yc[i] * zc[i] * SYZ;
            phi_tmp[128 + i] += yc[i] * SY;
            phi_tmp[128 + i] += zc[i] * SZ;
            phi_tmp[128 + i] +=  1 * S0[i];

            phi_tmp[160 + i] = zc_pow[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * zc[i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
        }

        // Combine ZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];

            phi_tmp[i] = xc_pow[i] * SZZ;

            phi_tmp[32 + i] = xc[i] * yc[i] * SZZ;

            phi_tmp[64 + i] = xc[i] * zc[i] * SZZ;
            phi_tmp[64 + i] += 2.0 * xc[i] * SZ;

            phi_tmp[96 + i] = yc_pow[i] * SZZ;

            phi_tmp[128 + i] = yc[i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * yc[i] * SZ;

            phi_tmp[160 + i] = zc_pow[i] * SZZ;
            phi_tmp[160 + i] += 4.0 * zc[i] * SZ;
            phi_tmp[160 + i] += 2.0 * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
        }

        // Combine XXX blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXX = S3[i] * xc[i] * xc[i] * xc[i] + 3 * xc[i] * S2[i];

            phi_tmp[i] = xc_pow[i] * SXXX;
            phi_tmp[i] += 3.0 * 2.0 * xc[i] * SXX;
            phi_tmp[i] += 3.0 * 2.0 * SX;

            phi_tmp[32 + i] = xc[i] * yc[i] * SXXX;
            phi_tmp[32 + i] += 3.0 * yc[i] * SXX;

            phi_tmp[64 + i] = xc[i] * zc[i] * SXXX;
            phi_tmp[64 + i] += 3.0 * zc[i] * SXX;

            phi_tmp[96 + i] = yc_pow[i] * SXXX;

            phi_tmp[128 + i] = yc[i] * zc[i] * SXXX;

            phi_tmp[160 + i] = zc_pow[i] * SXXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
        }

        // Combine XXY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXY = S3[i] * xc[i] * xc[i] * yc[i] + S2[i] * yc[i];

            phi_tmp[i] = xc_pow[i] * SXXY;
            phi_tmp[i] += 2.0 * 2.0 * xc[i] * SXY;
            phi_tmp[i] += 2.0 * SY;

            phi_tmp[32 + i] = xc[i] * yc[i] * SXXY;
            phi_tmp[32 + i] += 2.0 * yc[i] * SXY;
            phi_tmp[32 + i] += xc[i] * SXX;
            phi_tmp[32 + i] += 2.0 *  1 * SX;

            phi_tmp[64 + i] = xc[i] * zc[i] * SXXY;
            phi_tmp[64 + i] += 2.0 * zc[i] * SXY;

            phi_tmp[96 + i] = yc_pow[i] * SXXY;
            phi_tmp[96 + i] += 2.0 * yc[i] * SXX;

            phi_tmp[128 + i] = yc[i] * zc[i] * SXXY;
            phi_tmp[128 + i] += zc[i] * SXX;

            phi_tmp[160 + i] = zc_pow[i] * SXXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
        }

        // Combine XXZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXZ = S3[i] * xc[i] * xc[i] * zc[i] + S2[i] * zc[i];

            phi_tmp[i] = xc_pow[i] * SXXZ;
            phi_tmp[i] += 2.0 * 2.0 * xc[i] * SXZ;
            phi_tmp[i] += 2.0 * SZ;

            phi_tmp[32 + i] = xc[i] * yc[i] * SXXZ;
            phi_tmp[32 + i] += 2.0 * yc[i] * SXZ;

            phi_tmp[64 + i] = xc[i] * zc[i] * SXXZ;
            phi_tmp[64 + i] += 2.0 * zc[i] * SXZ;
            phi_tmp[64 + i] += xc[i] * SXX;
            phi_tmp[64 + i] += 2.0 *  1 * SX;

            phi_tmp[96 + i] = yc_pow[i] * SXXZ;

            phi_tmp[128 + i] = yc[i] * zc[i] * SXXZ;
            phi_tmp[128 + i] += yc[i] * SXX;

            phi_tmp[160 + i] = zc_pow[i] * SXXZ;
            phi_tmp[160 + i] += 2.0 * zc[i] * SXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
        }

        // Combine XYY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SXYY = S3[i] * xc[i] * yc[i] * yc[i] + S2[i] * xc[i];

            phi_tmp[i] = xc_pow[i] * SXYY;
            phi_tmp[i] += 2.0 * xc[i] * SYY;

            phi_tmp[32 + i] = xc[i] * yc[i] * SXYY;
            phi_tmp[32 + i] += 2.0 * xc[i] * SXY;
            phi_tmp[32 + i] += yc[i] * SYY;
            phi_tmp[32 + i] += 2.0 *  1 * SY;

            phi_tmp[64 + i] = xc[i] * zc[i] * SXYY;
            phi_tmp[64 + i] += zc[i] * SYY;

            phi_tmp[96 + i] = yc_pow[i] * SXYY;
            phi_tmp[96 + i] += 2.0 * 2.0 * yc[i] * SXY;
            phi_tmp[96 + i] += 2.0 * SX;

            phi_tmp[128 + i] = yc[i] * zc[i] * SXYY;
            phi_tmp[128 + i] += 2.0 * zc[i] * SXY;

            phi_tmp[160 + i] = zc_pow[i] * SXYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
        }

        // Combine XYZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SXYZ = S3[i] * xc[i] * yc[i] * zc[i];

            phi_tmp[i] = xc_pow[i] * SXYZ;
            phi_tmp[i] += 2.0 * xc[i] * SYZ;

            phi_tmp[32 + i] = xc[i] * yc[i] * SXYZ;
            phi_tmp[32 + i] += yc[i] * SYZ;
            phi_tmp[32 + i] += xc[i] * SXZ;
            phi_tmp[32 + i] +=  1 * SZ;

            phi_tmp[64 + i] = xc[i] * zc[i] * SXYZ;
            phi_tmp[64 + i] += zc[i] * SYZ;
            phi_tmp[64 + i] += xc[i] * SXY;
            phi_tmp[64 + i] +=  1 * SY;

            phi_tmp[96 + i] = yc_pow[i] * SXYZ;
            phi_tmp[96 + i] += 2.0 * yc[i] * SXZ;

            phi_tmp[128 + i] = yc[i] * zc[i] * SXYZ;
            phi_tmp[128 + i] += zc[i] * SXZ;
            phi_tmp[128 + i] += yc[i] * SXY;
            phi_tmp[128 + i] +=  1 * SX;

            phi_tmp[160 + i] = zc_pow[i] * SXYZ;
            phi_tmp[160 + i] += 2.0 * zc[i] * SXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
        }

        // Combine XZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SXZZ = S3[i] * xc[i] * zc[i] * zc[i] + S2[i] * xc[i];

            phi_tmp[i] = xc_pow[i] * SXZZ;
            phi_tmp[i] += 2.0 * xc[i] * SZZ;

            phi_tmp[32 + i] = xc[i] * yc[i] * SXZZ;
            phi_tmp[32 + i] += yc[i] * SZZ;

            phi_tmp[64 + i] = xc[i] * zc[i] * SXZZ;
            phi_tmp[64 + i] += 2.0 * xc[i] * SXZ;
            phi_tmp[64 + i] += zc[i] * SZZ;
            phi_tmp[64 + i] += 2.0 *  1 * SZ;

            phi_tmp[96 + i] = yc_pow[i] * SXZZ;

            phi_tmp[128 + i] = yc[i] * zc[i] * SXZZ;
            phi_tmp[128 + i] += 2.0 * yc[i] * SXZ;

            phi_tmp[160 + i] = zc_pow[i] * SXZZ;
            phi_tmp[160 + i] += 2.0 * 2.0 * zc[i] * SXZ;
            phi_tmp[160 + i] += 2.0 * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
        }

        // Combine YYY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SYYY = S3[i] * yc[i] * yc[i] * yc[i] + 3 * yc[i] * S2[i];

            phi_tmp[i] = xc_pow[i] * SYYY;

            phi_tmp[32 + i] = xc[i] * yc[i] * SYYY;
            phi_tmp[32 + i] += 3.0 * xc[i] * SYY;

            phi_tmp[64 + i] = xc[i] * zc[i] * SYYY;

            phi_tmp[96 + i] = yc_pow[i] * SYYY;
            phi_tmp[96 + i] += 3.0 * 2.0 * yc[i] * SYY;
            phi_tmp[96 + i] += 3.0 * 2.0 * SY;

            phi_tmp[128 + i] = yc[i] * zc[i] * SYYY;
            phi_tmp[128 + i] += 3.0 * zc[i] * SYY;

            phi_tmp[160 + i] = zc_pow[i] * SYYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
        }

        // Combine YYZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SYYZ = S3[i] * yc[i] * yc[i] * zc[i] + S2[i] * zc[i];

            phi_tmp[i] = xc_pow[i] * SYYZ;

            phi_tmp[32 + i] = xc[i] * yc[i] * SYYZ;
            phi_tmp[32 + i] += 2.0 * xc[i] * SYZ;

            phi_tmp[64 + i] = xc[i] * zc[i] * SYYZ;
            phi_tmp[64 + i] += xc[i] * SYY;

            phi_tmp[96 + i] = yc_pow[i] * SYYZ;
            phi_tmp[96 + i] += 2.0 * 2.0 * yc[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * SZ;

            phi_tmp[128 + i] = yc[i] * zc[i] * SYYZ;
            phi_tmp[128 + i] += 2.0 * zc[i] * SYZ;
            phi_tmp[128 + i] += yc[i] * SYY;
            phi_tmp[128 + i] += 2.0 *  1 * SY;

            phi_tmp[160 + i] = zc_pow[i] * SYYZ;
            phi_tmp[160 + i] += 2.0 * zc[i] * SYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
        }

        // Combine YZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SYZZ = S3[i] * yc[i] * zc[i] * zc[i] + S2[i] * yc[i];

            phi_tmp[i] = xc_pow[i] * SYZZ;

            phi_tmp[32 + i] = xc[i] * yc[i] * SYZZ;
            phi_tmp[32 + i] += xc[i] * SZZ;

            phi_tmp[64 + i] = xc[i] * zc[i] * SYZZ;
            phi_tmp[64 + i] += 2.0 * xc[i] * SYZ;

            phi_tmp[96 + i] = yc_pow[i] * SYZZ;
            phi_tmp[96 + i] += 2.0 * yc[i] * SZZ;

            phi_tmp[128 + i] = yc[i] * zc[i] * SYZZ;
            phi_tmp[128 + i] += 2.0 * yc[i] * SYZ;
            phi_tmp[128 + i] += zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 *  1 * SZ;

            phi_tmp[160 + i] = zc_pow[i] * SYZZ;
            phi_tmp[160 + i] += 2.0 * 2.0 * zc[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
        }

        // Combine ZZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SZZZ = S3[i] * zc[i] * zc[i] * zc[i] + 3 * zc[i] * S2[i];

            phi_tmp[i] = xc_pow[i] * SZZZ;

            phi_tmp[32 + i] = xc[i] * yc[i] * SZZZ;

            phi_tmp[64 + i] = xc[i] * zc[i] * SZZZ;
            phi_tmp[64 + i] += 3.0 * xc[i] * SZZ;

            phi_tmp[96 + i] = yc_pow[i] * SZZZ;

            phi_tmp[128 + i] = yc[i] * zc[i] * SZZZ;
            phi_tmp[128 + i] += 3.0 * yc[i] * SZZ;

            phi_tmp[160 + i] = zc_pow[i] * SZZZ;
            phi_tmp[160 + i] += 3.0 * 2.0 * zc[i] * SZZ;
            phi_tmp[160 + i] += 3.0 * 2.0 * SZ;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
        }

    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);
    ALIGNED_FREE(expn2);

    // Free Power temporaries
    ALIGNED_FREE(xc_pow);
    ALIGNED_FREE(yc_pow);
    ALIGNED_FREE(zc_pow);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);

}

void gg_collocation_L3_deriv3(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out, double* PRAGMA_RESTRICT phi_xxx_out, double* PRAGMA_RESTRICT phi_xxy_out, double* PRAGMA_RESTRICT phi_xxz_out, double* PRAGMA_RESTRICT phi_xyy_out, double* PRAGMA_RESTRICT phi_xyz_out, double* PRAGMA_RESTRICT phi_xzz_out, double* PRAGMA_RESTRICT phi_yyy_out, double* PRAGMA_RESTRICT phi_yyz_out, double* PRAGMA_RESTRICT phi_yzz_out, double* PRAGMA_RESTRICT phi_zzz_out) {

    // Sizing
    unsigned long nblocks = npoints / 32;
    nblocks += (npoints % 32) ? 1 : 0;
    const unsigned long ncart = 10;
    const unsigned long nspherical = 7;
    unsigned long nout;

    if ((order == GG_SPHERICAL_CCA) || (order == GG_SPHERICAL_GAUSSIAN)) {
        nout = nspherical;
        } else {
        nout = ncart;
    }

    // Allocate S temporaries, single block to stay on cache
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 288 * sizeof(double));
    double* PRAGMA_RESTRICT xc = cache_data + 0;
    ASSUME_ALIGNED(xc, 64);
    double* PRAGMA_RESTRICT yc = cache_data + 32;
    ASSUME_ALIGNED(yc, 64);
    double* PRAGMA_RESTRICT zc = cache_data + 64;
    ASSUME_ALIGNED(zc, 64);
    double* PRAGMA_RESTRICT R2 = cache_data + 96;
    ASSUME_ALIGNED(R2, 64);
    double* PRAGMA_RESTRICT S0 = cache_data + 128;
    ASSUME_ALIGNED(S0, 64);
    double* PRAGMA_RESTRICT tmp1 = cache_data + 160;
    ASSUME_ALIGNED(tmp1, 64);
    double* PRAGMA_RESTRICT S1 = cache_data + 192;
    ASSUME_ALIGNED(S1, 64);
    double* PRAGMA_RESTRICT S2 = cache_data + 224;
    ASSUME_ALIGNED(S2, 64);
    double* PRAGMA_RESTRICT S3 = cache_data + 256;
    ASSUME_ALIGNED(S3, 64);

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));
    double* PRAGMA_RESTRICT expn2 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));

    // Allocate power temporaries
    double* PRAGMA_RESTRICT xc_pow = (double*)ALIGNED_MALLOC(64, 64 * sizeof(double));
    ASSUME_ALIGNED(xc_pow, 64);
    double* PRAGMA_RESTRICT yc_pow = (double*)ALIGNED_MALLOC(64, 64 * sizeof(double));
    ASSUME_ALIGNED(yc_pow, 64);
    double* PRAGMA_RESTRICT zc_pow = (double*)ALIGNED_MALLOC(64, 64 * sizeof(double));
    ASSUME_ALIGNED(zc_pow, 64);

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, 320 * sizeof(double));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;
    double AXX, AXY, AXZ, AYY, AYZ, AZZ;
    double AXXX, XXY, XXZ, XYY, XYZ, XZZ, YYY, YYZ, YZZ, ZZZ;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
        expn2[i] = -2.0 * exponents[i];
    }

    // Start outer block loop
    for (unsigned long block = 0; block < nblocks; block++) {


        // Copy data into inner temps
        const unsigned long start = block * 32;
        const unsigned long remain = ((start + 32) > npoints) ? (npoints - start) : 32;

        // Handle non-AM dependant temps
        if (xyz_stride == 1) {
            const double* PRAGMA_RESTRICT x = xyz + start;
            const double* PRAGMA_RESTRICT y = xyz + npoints + start;
            const double* PRAGMA_RESTRICT z = xyz + 2 * npoints + start;
            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = x[i] - center_x;
                yc[i] = y[i] - center_y;
                zc[i] = z[i] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
            } else {
            unsigned int start_shift = start * xyz_stride;

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = xyz[start_shift + i * xyz_stride] - center_x;
                yc[i] = xyz[start_shift + i * xyz_stride + 1] - center_y;
                zc[i] = xyz[start_shift + i * xyz_stride + 2] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];
            const double alpha_n2 = expn2[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
                const double T2 = alpha_n2 * T1;
                S1[i] += T2;
                const double T3 = alpha_n2 * T2;
                S2[i] += T3;
                const double T4 = alpha_n2 * T3;
                S3[i] += T4;
            }

        }

        // Build powers
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {

            // Cartesian derivs
            xc_pow[i] = xc[i] * xc[i];
            yc_pow[i] = yc[i] * yc[i];
            zc_pow[i] = zc[i] * zc[i];
            xc_pow[32 + i] = xc_pow[i] * xc[i];
            yc_pow[32 + i] = yc_pow[i] * yc[i];
            zc_pow[32 + i] = zc_pow[i] * zc[i];
        }
        // Combine A blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {

            phi_tmp[i] = xc_pow[32 + i] * S0[i];
            phi_tmp[32 + i] = xc_pow[i] * yc[i] * S0[i];
            phi_tmp[64 + i] = xc_pow[i] * zc[i] * S0[i];
            phi_tmp[96 + i] = xc[i] * yc_pow[i] * S0[i];
            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * S0[i];
            phi_tmp[160 + i] = xc[i] * zc_pow[i] * S0[i];
            phi_tmp[192 + i] = yc_pow[32 + i] * S0[i];
            phi_tmp[224 + i] = yc_pow[i] * zc[i] * S0[i];
            phi_tmp[256 + i] = yc[i] * zc_pow[i] * S0[i];
            phi_tmp[288 + i] = zc_pow[32 + i] * S0[i];
        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_out + start), npoints);
        }

        // Combine X blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];

            phi_tmp[i] = xc_pow[32 + i] * SX;
            phi_tmp[i] += 3.0 * xc_pow[i] * S0[i];

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SX;
            phi_tmp[32 + i] += 2.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SX;
            phi_tmp[64 + i] += 2.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SX;
            phi_tmp[96 + i] += yc_pow[i] * S0[i];

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SX;
            phi_tmp[128 + i] += yc[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SX;
            phi_tmp[160 + i] += zc_pow[i] * S0[i];

            phi_tmp[192 + i] = yc_pow[32 + i] * SX;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SX;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SX;

            phi_tmp[288 + i] = zc_pow[32 + i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_x_out + start), npoints);
        }

        // Combine Y blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];

            phi_tmp[i] = xc_pow[32 + i] * SY;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SY;
            phi_tmp[32 + i] += xc_pow[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SY;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SY;
            phi_tmp[96 + i] += 2.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SY;
            phi_tmp[128 + i] += xc[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SY;

            phi_tmp[192 + i] = yc_pow[32 + i] * SY;
            phi_tmp[192 + i] += 3.0 * yc_pow[i] * S0[i];

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SY;
            phi_tmp[224 + i] += 2.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SY;
            phi_tmp[256 + i] += zc_pow[i] * S0[i];

            phi_tmp[288 + i] = zc_pow[32 + i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_y_out + start), npoints);
        }

        // Combine Z blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];

            phi_tmp[i] = xc_pow[32 + i] * SZ;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SZ;
            phi_tmp[64 + i] += xc_pow[i] * S0[i];

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SZ;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += xc[i] * yc[i] * S0[i];

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SZ;
            phi_tmp[160 + i] += 2.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[192 + i] = yc_pow[32 + i] * SZ;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SZ;
            phi_tmp[224 + i] += yc_pow[i] * S0[i];

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = zc_pow[32 + i] * SZ;
            phi_tmp[288 + i] += 3.0 * zc_pow[i] * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_z_out + start), npoints);
        }

        // Combine XX blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];

            phi_tmp[i] = xc_pow[32 + i] * SXX;
            phi_tmp[i] += 6.0 * xc_pow[i] * SX;
            phi_tmp[i] += 6.0 * xc[i] * S0[i];

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SXX;
            phi_tmp[32 + i] += 4.0 * xc[i] * yc[i] * SX;
            phi_tmp[32 + i] += 2.0 * yc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SXX;
            phi_tmp[64 + i] += 4.0 * xc[i] * zc[i] * SX;
            phi_tmp[64 + i] += 2.0 * zc[i] * S0[i];

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SXX;
            phi_tmp[96 + i] += 2.0 * yc_pow[i] * SX;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SXX;
            phi_tmp[128 + i] += 2.0 * yc[i] * zc[i] * SX;

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SXX;
            phi_tmp[160 + i] += 2.0 * zc_pow[i] * SX;

            phi_tmp[192 + i] = yc_pow[32 + i] * SXX;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SXX;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SXX;

            phi_tmp[288 + i] = zc_pow[32 + i] * SXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
        }

        // Combine XY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];

            phi_tmp[i] = xc_pow[32 + i] * SXY;
            phi_tmp[i] += 3.0 * xc_pow[i] * SY;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SXY;
            phi_tmp[32 + i] += xc_pow[i] * SX;
            phi_tmp[32 + i] += 2.0 * xc[i] * yc[i] * SY;
            phi_tmp[32 + i] += 2.0 * xc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SXY;
            phi_tmp[64 + i] += 2.0 * xc[i] * zc[i] * SY;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SXY;
            phi_tmp[96 + i] += 2.0 * xc[i] * yc[i] * SX;
            phi_tmp[96 + i] += yc_pow[i] * SY;
            phi_tmp[96 + i] += 2.0 * yc[i] * S0[i];

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SXY;
            phi_tmp[128 + i] += xc[i] * zc[i] * SX;
            phi_tmp[128 + i] += yc[i] * zc[i] * SY;
            phi_tmp[128 + i] += zc[i] * S0[i];

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SXY;
            phi_tmp[160 + i] += zc_pow[i] * SY;

            phi_tmp[192 + i] = yc_pow[32 + i] * SXY;
            phi_tmp[192 + i] += 3.0 * yc_pow[i] * SX;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SXY;
            phi_tmp[224 + i] += 2.0 * yc[i] * zc[i] * SX;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SXY;
            phi_tmp[256 + i] += zc_pow[i] * SX;

            phi_tmp[288 + i] = zc_pow[32 + i] * SXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
        }

        // Combine XZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];

            phi_tmp[i] = xc_pow[32 + i] * SXZ;
            phi_tmp[i] += 3.0 * xc_pow[i] * SZ;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SXZ;
            phi_tmp[32 + i] += 2.0 * xc[i] * yc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SXZ;
            phi_tmp[64 + i] += xc_pow[i] * SX;
            phi_tmp[64 + i] += 2.0 * xc[i] * zc[i] * SZ;
            phi_tmp[64 + i] += 2.0 * xc[i] * S0[i];

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SXZ;
            phi_tmp[96 + i] += yc_pow[i] * SZ;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SXZ;
            phi_tmp[128 + i] += xc[i] * yc[i] * SX;
            phi_tmp[128 + i] += yc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += yc[i] * S0[i];

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SXZ;
            phi_tmp[160 + i] += 2.0 * xc[i] * zc[i] * SX;
            phi_tmp[160 + i] += zc_pow[i] * SZ;
            phi_tmp[160 + i] += 2.0 * zc[i] * S0[i];

            phi_tmp[192 + i] = yc_pow[32 + i] * SXZ;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SXZ;
            phi_tmp[224 + i] += yc_pow[i] * SX;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc[i] * SX;

            phi_tmp[288 + i] = zc_pow[32 + i] * SXZ;
            phi_tmp[288 + i] += 3.0 * zc_pow[i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
        }

        // Combine YY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];

            phi_tmp[i] = xc_pow[32 + i] * SYY;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SYY;
            phi_tmp[32 + i] += 2.0 * xc_pow[i] * SY;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SYY;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SYY;
            phi_tmp[96 + i] += 4.0 * xc[i] * yc[i] * SY;
            phi_tmp[96 + i] += 2.0 * xc[i] * S0[i];

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * xc[i] * zc[i] * SY;

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SYY;

            phi_tmp[192 + i] = yc_pow[32 + i] * SYY;
            phi_tmp[192 + i] += 6.0 * yc_pow[i] * SY;
            phi_tmp[192 + i] += 6.0 * yc[i] * S0[i];

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SYY;
            phi_tmp[224 + i] += 4.0 * yc[i] * zc[i] * SY;
            phi_tmp[224 + i] += 2.0 * zc[i] * S0[i];

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SYY;
            phi_tmp[256 + i] += 2.0 * zc_pow[i] * SY;

            phi_tmp[288 + i] = zc_pow[32 + i] * SYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
        }

        // Combine YZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];

            phi_tmp[i] = xc_pow[32 + i] * SYZ;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SYZ;
            phi_tmp[32 + i] += xc_pow[i] * SZ;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SYZ;
            phi_tmp[64 + i] += xc_pow[i] * SY;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * xc[i] * yc[i] * SZ;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SYZ;
            phi_tmp[128 + i] += xc[i] * yc[i] * SY;
            phi_tmp[128 + i] += xc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += xc[i] * S0[i];

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * xc[i] * zc[i] * SY;

            phi_tmp[192 + i] = yc_pow[32 + i] * SYZ;
            phi_tmp[192 + i] += 3.0 * yc_pow[i] * SZ;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SYZ;
            phi_tmp[224 + i] += yc_pow[i] * SY;
            phi_tmp[224 + i] += 2.0 * yc[i] * zc[i] * SZ;
            phi_tmp[224 + i] += 2.0 * yc[i] * S0[i];

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SYZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc[i] * SY;
            phi_tmp[256 + i] += zc_pow[i] * SZ;
            phi_tmp[256 + i] += 2.0 * zc[i] * S0[i];

            phi_tmp[288 + i] = zc_pow[32 + i] * SYZ;
            phi_tmp[288 + i] += 3.0 * zc_pow[i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
        }

        // Combine ZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];

            phi_tmp[i] = xc_pow[32 + i] * SZZ;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SZZ;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SZZ;
            phi_tmp[64 + i] += 2.0 * xc_pow[i] * SZ;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SZZ;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * xc[i] * yc[i] * SZ;

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SZZ;
            phi_tmp[160 + i] += 4.0 * xc[i] * zc[i] * SZ;
            phi_tmp[160 + i] += 2.0 * xc[i] * S0[i];

            phi_tmp[192 + i] = yc_pow[32 + i] * SZZ;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SZZ;
            phi_tmp[224 + i] += 2.0 * yc_pow[i] * SZ;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SZZ;
            phi_tmp[256 + i] += 4.0 * yc[i] * zc[i] * SZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * S0[i];

            phi_tmp[288 + i] = zc_pow[32 + i] * SZZ;
            phi_tmp[288 + i] += 6.0 * zc_pow[i] * SZ;
            phi_tmp[288 + i] += 6.0 * zc[i] * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
        }

        // Combine XXX blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXX = S3[i] * xc[i] * xc[i] * xc[i] + 3 * xc[i] * S2[i];

            phi_tmp[i] = xc_pow[32 + i] * SXXX;
            phi_tmp[i] += 3.0 * 3.0 * xc_pow[i] * SXX;
            phi_tmp[i] += 3.0 * 6.0 * xc[i] * SX;
            phi_tmp[i] += 6.0 * S0[i];

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SXXX;
            phi_tmp[32 + i] += 3.0 * 2.0 * xc[i] * yc[i] * SXX;
            phi_tmp[32 + i] += 3.0 * 2.0 * yc[i] * SX;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SXXX;
            phi_tmp[64 + i] += 3.0 * 2.0 * xc[i] * zc[i] * SXX;
            phi_tmp[64 + i] += 3.0 * 2.0 * zc[i] * SX;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SXXX;
            phi_tmp[96 + i] += 3.0 * yc_pow[i] * SXX;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SXXX;
            phi_tmp[128 + i] += 3.0 * yc[i] * zc[i] * SXX;

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SXXX;
            phi_tmp[160 + i] += 3.0 * zc_pow[i] * SXX;

            phi_tmp[192 + i] = yc_pow[32 + i] * SXXX;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SXXX;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SXXX;

            phi_tmp[288 + i] = zc_pow[32 + i] * SXXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
        }

        // Combine XXY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXY = S3[i] * xc[i] * xc[i] * yc[i] + S2[i] * yc[i];

            phi_tmp[i] = xc_pow[32 + i] * SXXY;
            phi_tmp[i] += 2.0 * 3.0 * xc_pow[i] * SXY;
            phi_tmp[i] += 6.0 * xc[i] * SY;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SXXY;
            phi_tmp[32 + i] += 2.0 * 2.0 * xc[i] * yc[i] * SXY;
            phi_tmp[32 + i] += xc_pow[i] * SXX;
            phi_tmp[32 + i] += 2.0 * yc[i] * SY;
            phi_tmp[32 + i] += 2.0 * 2.0 * xc[i] * SX;
            phi_tmp[32 + i] += 2.0 * S0[i];

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SXXY;
            phi_tmp[64 + i] += 2.0 * 2.0 * xc[i] * zc[i] * SXY;
            phi_tmp[64 + i] += 2.0 * zc[i] * SY;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SXXY;
            phi_tmp[96 + i] += 2.0 * yc_pow[i] * SXY;
            phi_tmp[96 + i] += 2.0 * xc[i] * yc[i] * SXX;
            phi_tmp[96 + i] += 2.0 * 2.0 * yc[i] * SX;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SXXY;
            phi_tmp[128 + i] += 2.0 * yc[i] * zc[i] * SXY;
            phi_tmp[128 + i] += xc[i] * zc[i] * SXX;
            phi_tmp[128 + i] += 2.0 * zc[i] * SX;

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SXXY;
            phi_tmp[160 + i] += 2.0 * zc_pow[i] * SXY;

            phi_tmp[192 + i] = yc_pow[32 + i] * SXXY;
            phi_tmp[192 + i] += 3.0 * yc_pow[i] * SXX;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SXXY;
            phi_tmp[224 + i] += 2.0 * yc[i] * zc[i] * SXX;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SXXY;
            phi_tmp[256 + i] += zc_pow[i] * SXX;

            phi_tmp[288 + i] = zc_pow[32 + i] * SXXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
        }

        // Combine XXZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXZ = S3[i] * xc[i] * xc[i] * zc[i] + S2[i] * zc[i];

            phi_tmp[i] = xc_pow[32 + i] * SXXZ;
            phi_tmp[i] += 2.0 * 3.0 * xc_pow[i] * SXZ;
            phi_tmp[i] += 6.0 * xc[i] * SZ;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SXXZ;
            phi_tmp[32 + i] += 2.0 * 2.0 * xc[i] * yc[i] * SXZ;
            phi_tmp[32 + i] += 2.0 * yc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SXXZ;
            phi_tmp[64 + i] += 2.0 * 2.0 * xc[i] * zc[i] * SXZ;
            phi_tmp[64 + i] += xc_pow[i] * SXX;
            phi_tmp[64 + i] += 2.0 * zc[i] * SZ;
            phi_tmp[64 + i] += 2.0 * 2.0 * xc[i] * SX;
            phi_tmp[64 + i] += 2.0 * S0[i];

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SXXZ;
            phi_tmp[96 + i] += 2.0 * yc_pow[i] * SXZ;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SXXZ;
            phi_tmp[128 + i] += 2.0 * yc[i] * zc[i] * SXZ;
            phi_tmp[128 + i] += xc[i] * yc[i] * SXX;
            phi_tmp[128 + i] += 2.0 * yc[i] * SX;

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SXXZ;
            phi_tmp[160 + i] += 2.0 * zc_pow[i] * SXZ;
            phi_tmp[160 + i] += 2.0 * xc[i] * zc[i] * SXX;
            phi_tmp[160 + i] += 2.0 * 2.0 * zc[i] * SX;

            phi_tmp[192 + i] = yc_pow[32 + i] * SXXZ;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SXXZ;
            phi_tmp[224 + i] += yc_pow[i] * SXX;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SXXZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc[i] * SXX;

            phi_tmp[288 + i] = zc_pow[32 + i] * SXXZ;
            phi_tmp[288 + i] += 3.0 * zc_pow[i] * SXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
        }

        // Combine XYY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SXYY = S3[i] * xc[i] * yc[i] * yc[i] + S2[i] * xc[i];

            phi_tmp[i] = xc_pow[32 + i] * SXYY;
            phi_tmp[i] += 3.0 * xc_pow[i] * SYY;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SXYY;
            phi_tmp[32 + i] += 2.0 * xc_pow[i] * SXY;
            phi_tmp[32 + i] += 2.0 * xc[i] * yc[i] * SYY;
            phi_tmp[32 + i] += 2.0 * 2.0 * xc[i] * SY;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SXYY;
            phi_tmp[64 + i] += 2.0 * xc[i] * zc[i] * SYY;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SXYY;
            phi_tmp[96 + i] += 2.0 * 2.0 * xc[i] * yc[i] * SXY;
            phi_tmp[96 + i] += yc_pow[i] * SYY;
            phi_tmp[96 + i] += 2.0 * xc[i] * SX;
            phi_tmp[96 + i] += 2.0 * 2.0 * yc[i] * SY;
            phi_tmp[96 + i] += 2.0 * S0[i];

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SXYY;
            phi_tmp[128 + i] += 2.0 * xc[i] * zc[i] * SXY;
            phi_tmp[128 + i] += yc[i] * zc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * zc[i] * SY;

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SXYY;
            phi_tmp[160 + i] += zc_pow[i] * SYY;

            phi_tmp[192 + i] = yc_pow[32 + i] * SXYY;
            phi_tmp[192 + i] += 2.0 * 3.0 * yc_pow[i] * SXY;
            phi_tmp[192 + i] += 6.0 * yc[i] * SX;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SXYY;
            phi_tmp[224 + i] += 2.0 * 2.0 * yc[i] * zc[i] * SXY;
            phi_tmp[224 + i] += 2.0 * zc[i] * SX;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SXYY;
            phi_tmp[256 + i] += 2.0 * zc_pow[i] * SXY;

            phi_tmp[288 + i] = zc_pow[32 + i] * SXYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
        }

        // Combine XYZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SXYZ = S3[i] * xc[i] * yc[i] * zc[i];

            phi_tmp[i] = xc_pow[32 + i] * SXYZ;
            phi_tmp[i] += 3.0 * xc_pow[i] * SYZ;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SXYZ;
            phi_tmp[32 + i] += 2.0 * xc[i] * yc[i] * SYZ;
            phi_tmp[32 + i] += xc_pow[i] * SXZ;
            phi_tmp[32 + i] += 2.0 * xc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SXYZ;
            phi_tmp[64 + i] += 2.0 * xc[i] * zc[i] * SYZ;
            phi_tmp[64 + i] += xc_pow[i] * SXY;
            phi_tmp[64 + i] += 2.0 * xc[i] * SY;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SXYZ;
            phi_tmp[96 + i] += yc_pow[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * xc[i] * yc[i] * SXZ;
            phi_tmp[96 + i] += 2.0 * yc[i] * SZ;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SXYZ;
            phi_tmp[128 + i] += yc[i] * zc[i] * SYZ;
            phi_tmp[128 + i] += xc[i] * zc[i] * SXZ;
            phi_tmp[128 + i] += xc[i] * yc[i] * SXY;
            phi_tmp[128 + i] += zc[i] * SZ;
            phi_tmp[128 + i] += yc[i] * SY;
            phi_tmp[128 + i] += xc[i] * SX;
            phi_tmp[128 + i] +=  1 * S0[i];

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SXYZ;
            phi_tmp[160 + i] += zc_pow[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * xc[i] * zc[i] * SXY;
            phi_tmp[160 + i] += 2.0 * zc[i] * SY;

            phi_tmp[192 + i] = yc_pow[32 + i] * SXYZ;
            phi_tmp[192 + i] += 3.0 * yc_pow[i] * SXZ;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SXYZ;
            phi_tmp[224 + i] += 2.0 * yc[i] * zc[i] * SXZ;
            phi_tmp[224 + i] += yc_pow[i] * SXY;
            phi_tmp[224 + i] += 2.0 * yc[i] * SX;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SXYZ;
            phi_tmp[256 + i] += zc_pow[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc[i] * SXY;
            phi_tmp[256 + i] += 2.0 * zc[i] * SX;

            phi_tmp[288 + i] = zc_pow[32 + i] * SXYZ;
            phi_tmp[288 + i] += 3.0 * zc_pow[i] * SXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
        }

        // Combine XZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SXZZ = S3[i] * xc[i] * zc[i] * zc[i] + S2[i] * xc[i];

            phi_tmp[i] = xc_pow[32 + i] * SXZZ;
            phi_tmp[i] += 3.0 * xc_pow[i] * SZZ;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SXZZ;
            phi_tmp[32 + i] += 2.0 * xc[i] * yc[i] * SZZ;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SXZZ;
            phi_tmp[64 + i] += 2.0 * xc_pow[i] * SXZ;
            phi_tmp[64 + i] += 2.0 * xc[i] * zc[i] * SZZ;
            phi_tmp[64 + i] += 2.0 * 2.0 * xc[i] * SZ;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SXZZ;
            phi_tmp[96 + i] += yc_pow[i] * SZZ;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SXZZ;
            phi_tmp[128 + i] += 2.0 * xc[i] * yc[i] * SXZ;
            phi_tmp[128 + i] += yc[i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * yc[i] * SZ;

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SXZZ;
            phi_tmp[160 + i] += 2.0 * 2.0 * xc[i] * zc[i] * SXZ;
            phi_tmp[160 + i] += zc_pow[i] * SZZ;
            phi_tmp[160 + i] += 2.0 * xc[i] * SX;
            phi_tmp[160 + i] += 2.0 * 2.0 * zc[i] * SZ;
            phi_tmp[160 + i] += 2.0 * S0[i];

            phi_tmp[192 + i] = yc_pow[32 + i] * SXZZ;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SXZZ;
            phi_tmp[224 + i] += 2.0 * yc_pow[i] * SXZ;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SXZZ;
            phi_tmp[256 + i] += 2.0 * 2.0 * yc[i] * zc[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * SX;

            phi_tmp[288 + i] = zc_pow[32 + i] * SXZZ;
            phi_tmp[288 + i] += 2.0 * 3.0 * zc_pow[i] * SXZ;
            phi_tmp[288 + i] += 6.0 * zc[i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
        }

        // Combine YYY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SYYY = S3[i] * yc[i] * yc[i] * yc[i] + 3 * yc[i] * S2[i];

            phi_tmp[i] = xc_pow[32 + i] * SYYY;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SYYY;
            phi_tmp[32 + i] += 3.0 * xc_pow[i] * SYY;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SYYY;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SYYY;
            phi_tmp[96 + i] += 3.0 * 2.0 * xc[i] * yc[i] * SYY;
            phi_tmp[96 + i] += 3.0 * 2.0 * xc[i] * SY;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SYYY;
            phi_tmp[128 + i] += 3.0 * xc[i] * zc[i] * SYY;

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SYYY;

            phi_tmp[192 + i] = yc_pow[32 + i] * SYYY;
            phi_tmp[192 + i] += 3.0 * 3.0 * yc_pow[i] * SYY;
            phi_tmp[192 + i] += 3.0 * 6.0 * yc[i] * SY;
            phi_tmp[192 + i] += 6.0 * S0[i];

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SYYY;
            phi_tmp[224 + i] += 3.0 * 2.0 * yc[i] * zc[i] * SYY;
            phi_tmp[224 + i] += 3.0 * 2.0 * zc[i] * SY;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SYYY;
            phi_tmp[256 + i] += 3.0 * zc_pow[i] * SYY;

            phi_tmp[288 + i] = zc_pow[32 + i] * SYYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
        }

        // Combine YYZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SYYZ = S3[i] * yc[i] * yc[i] * zc[i] + S2[i] * zc[i];

            phi_tmp[i] = xc_pow[32 + i] * SYYZ;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SYYZ;
            phi_tmp[32 + i] += 2.0 * xc_pow[i] * SYZ;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SYYZ;
            phi_tmp[64 + i] += xc_pow[i] * SYY;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SYYZ;
            phi_tmp[96 + i] += 2.0 * 2.0 * xc[i] * yc[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * xc[i] * SZ;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SYYZ;
            phi_tmp[128 + i] += 2.0 * xc[i] * zc[i] * SYZ;
            phi_tmp[128 + i] += xc[i] * yc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * xc[i] * SY;

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SYYZ;
            phi_tmp[160 + i] += 2.0 * xc[i] * zc[i] * SYY;

            phi_tmp[192 + i] = yc_pow[32 + i] * SYYZ;
            phi_tmp[192 + i] += 2.0 * 3.0 * yc_pow[i] * SYZ;
            phi_tmp[192 + i] += 6.0 * yc[i] * SZ;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SYYZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * yc[i] * zc[i] * SYZ;
            phi_tmp[224 + i] += yc_pow[i] * SYY;
            phi_tmp[224 + i] += 2.0 * zc[i] * SZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * yc[i] * SY;
            phi_tmp[224 + i] += 2.0 * S0[i];

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SYYZ;
            phi_tmp[256 + i] += 2.0 * zc_pow[i] * SYZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc[i] * SYY;
            phi_tmp[256 + i] += 2.0 * 2.0 * zc[i] * SY;

            phi_tmp[288 + i] = zc_pow[32 + i] * SYYZ;
            phi_tmp[288 + i] += 3.0 * zc_pow[i] * SYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
        }

        // Combine YZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SYZZ = S3[i] * yc[i] * zc[i] * zc[i] + S2[i] * yc[i];

            phi_tmp[i] = xc_pow[32 + i] * SYZZ;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SYZZ;
            phi_tmp[32 + i] += xc_pow[i] * SZZ;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SYZZ;
            phi_tmp[64 + i] += 2.0 * xc_pow[i] * SYZ;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SYZZ;
            phi_tmp[96 + i] += 2.0 * xc[i] * yc[i] * SZZ;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SYZZ;
            phi_tmp[128 + i] += 2.0 * xc[i] * yc[i] * SYZ;
            phi_tmp[128 + i] += xc[i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * xc[i] * SZ;

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SYZZ;
            phi_tmp[160 + i] += 2.0 * 2.0 * xc[i] * zc[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * xc[i] * SY;

            phi_tmp[192 + i] = yc_pow[32 + i] * SYZZ;
            phi_tmp[192 + i] += 3.0 * yc_pow[i] * SZZ;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SYZZ;
            phi_tmp[224 + i] += 2.0 * yc_pow[i] * SYZ;
            phi_tmp[224 + i] += 2.0 * yc[i] * zc[i] * SZZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * yc[i] * SZ;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SYZZ;
            phi_tmp[256 + i] += 2.0 * 2.0 * yc[i] * zc[i] * SYZ;
            phi_tmp[256 + i] += zc_pow[i] * SZZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * SY;
            phi_tmp[256 + i] += 2.0 * 2.0 * zc[i] * SZ;
            phi_tmp[256 + i] += 2.0 * S0[i];

            phi_tmp[288 + i] = zc_pow[32 + i] * SYZZ;
            phi_tmp[288 + i] += 2.0 * 3.0 * zc_pow[i] * SYZ;
            phi_tmp[288 + i] += 6.0 * zc[i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
        }

        // Combine ZZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SZZZ = S3[i] * zc[i] * zc[i] * zc[i] + 3 * zc[i] * S2[i];

            phi_tmp[i] = xc_pow[32 + i] * SZZZ;

            phi_tmp[32 + i] = xc_pow[i] * yc[i] * SZZZ;

            phi_tmp[64 + i] = xc_pow[i] * zc[i] * SZZZ;
            phi_tmp[64 + i] += 3.0 * xc_pow[i] * SZZ;

            phi_tmp[96 + i] = xc[i] * yc_pow[i] * SZZZ;

            phi_tmp[128 + i] = xc[i] * yc[i] * zc[i] * SZZZ;
            phi_tmp[128 + i] += 3.0 * xc[i] * yc[i] * SZZ;

            phi_tmp[160 + i] = xc[i] * zc_pow[i] * SZZZ;
            phi_tmp[160 + i] += 3.0 * 2.0 * xc[i] * zc[i] * SZZ;
            phi_tmp[160 + i] += 3.0 * 2.0 * xc[i] * SZ;

            phi_tmp[192 + i] = yc_pow[32 + i] * SZZZ;

            phi_tmp[224 + i] = yc_pow[i] * zc[i] * SZZZ;
            phi_tmp[224 + i] += 3.0 * yc_pow[i] * SZZ;

            phi_tmp[256 + i] = yc[i] * zc_pow[i] * SZZZ;
            phi_tmp[256 + i] += 3.0 * 2.0 * yc[i] * zc[i] * SZZ;
            phi_tmp[256 + i] += 3.0 * 2.0 * yc[i] * SZ;

            phi_tmp[288 + i] = zc_pow[32 + i] * SZZZ;
            phi_tmp[288 + i] += 3.0 * 3.0 * zc_pow[i] * SZZ;
            phi_tmp[288 + i] += 3.0 * 6.0 * zc[i] * SZ;
            phi_tmp[288 + i] += 6.0 * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
        }

    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);
    ALIGNED_FREE(expn2);

    // Free Power temporaries
    ALIGNED_FREE(xc_pow);
    ALIGNED_FREE(yc_pow);
    ALIGNED_FREE(zc_pow);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);

}

void gg_collocation_L4_deriv3(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out, double* PRAGMA_RESTRICT phi_xxx_out, double* PRAGMA_RESTRICT phi_xxy_out, double* PRAGMA_RESTRICT phi_xxz_out, double* PRAGMA_RESTRICT phi_xyy_out, double* PRAGMA_RESTRICT phi_xyz_out, double* PRAGMA_RESTRICT phi_xzz_out, double* PRAGMA_RESTRICT phi_yyy_out, double* PRAGMA_RESTRICT phi_yyz_out, double* PRAGMA_RESTRICT phi_yzz_out, double* PRAGMA_RESTRICT phi_zzz_out) {

    // Sizing
    unsigned long nblocks = npoints / 32;
    nblocks += (npoints % 32) ? 1 : 0;
    const unsigned long ncart = 15;
    const unsigned long nspherical = 9;
    unsigned long nout;

    if ((order == GG_SPHERICAL_CCA) || (order == GG_SPHERICAL_GAUSSIAN)) {
        nout = nspherical;
        } else {
        nout = ncart;
    }

    // Allocate S temporaries, single block to stay on cache
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 288 * sizeof(double));
    double* PRAGMA_RESTRICT xc = cache_data + 0;
    ASSUME_ALIGNED(xc, 64);
    double* PRAGMA_RESTRICT yc = cache_data + 32;
    ASSUME_ALIGNED(yc, 64);
    double* PRAGMA_RESTRICT zc = cache_data + 64;
    ASSUME_ALIGNED(zc, 64);
    double* PRAGMA_RESTRICT R2 = cache_data + 96;
    ASSUME_ALIGNED(R2, 64);
    double* PRAGMA_RESTRICT S0 = cache_data + 128;
    ASSUME_ALIGNED(S0, 64);
    double* PRAGMA_RESTRICT tmp1 = cache_data + 160;
    ASSUME_ALIGNED(tmp1, 64);
    double* PRAGMA_RESTRICT S1 = cache_data + 192;
    ASSUME_ALIGNED(S1, 64);
    double* PRAGMA_RESTRICT S2 = cache_data + 224;
    ASSUME_ALIGNED(S2, 64);
    double* PRAGMA_RESTRICT S3 = cache_data + 256;
    ASSUME_ALIGNED(S3, 64);

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));
    double* PRAGMA_RESTRICT expn2 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));

    // Allocate power temporaries
    double* PRAGMA_RESTRICT xc_pow = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(xc_pow, 64);
    double* PRAGMA_RESTRICT yc_pow = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(yc_pow, 64);
    double* PRAGMA_RESTRICT zc_pow = (double*)ALIGNED_MALLOC(64, 96 * sizeof(double));
    ASSUME_ALIGNED(zc_pow, 64);

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, 480 * sizeof(double));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;
    double AXX, AXY, AXZ, AYY, AYZ, AZZ;
    double AXXX, XXY, XXZ, XYY, XYZ, XZZ, YYY, YYZ, YZZ, ZZZ;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
        expn2[i] = -2.0 * exponents[i];
    }

    // Start outer block loop
    for (unsigned long block = 0; block < nblocks; block++) {


        // Copy data into inner temps
        const unsigned long start = block * 32;
        const unsigned long remain = ((start + 32) > npoints) ? (npoints - start) : 32;

        // Handle non-AM dependant temps
        if (xyz_stride == 1) {
            const double* PRAGMA_RESTRICT x = xyz + start;
            const double* PRAGMA_RESTRICT y = xyz + npoints + start;
            const double* PRAGMA_RESTRICT z = xyz + 2 * npoints + start;
            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = x[i] - center_x;
                yc[i] = y[i] - center_y;
                zc[i] = z[i] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
            } else {
            unsigned int start_shift = start * xyz_stride;

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = xyz[start_shift + i * xyz_stride] - center_x;
                yc[i] = xyz[start_shift + i * xyz_stride + 1] - center_y;
                zc[i] = xyz[start_shift + i * xyz_stride + 2] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];
            const double alpha_n2 = expn2[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
                const double T2 = alpha_n2 * T1;
                S1[i] += T2;
                const double T3 = alpha_n2 * T2;
                S2[i] += T3;
                const double T4 = alpha_n2 * T3;
                S3[i] += T4;
            }

        }

        // Build powers
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {

            // Cartesian derivs
            xc_pow[i] = xc[i] * xc[i];
            yc_pow[i] = yc[i] * yc[i];
            zc_pow[i] = zc[i] * zc[i];
            xc_pow[32 + i] = xc_pow[i] * xc[i];
            yc_pow[32 + i] = yc_pow[i] * yc[i];
            zc_pow[32 + i] = zc_pow[i] * zc[i];
            xc_pow[64 + i] = xc_pow[32 + i] * xc[i];
            yc_pow[64 + i] = yc_pow[32 + i] * yc[i];
            zc_pow[64 + i] = zc_pow[32 + i] * zc[i];
        }
        // Combine A blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {

            phi_tmp[i] = xc_pow[64 + i] * S0[i];
            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * S0[i];
            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * S0[i];
            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * S0[i];
            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * S0[i];
            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * S0[i];
            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * S0[i];
            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * S0[i];
            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * S0[i];
            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * S0[i];
            phi_tmp[320 + i] = yc_pow[64 + i] * S0[i];
            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * S0[i];
            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * S0[i];
            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * S0[i];
            phi_tmp[448 + i] = zc_pow[64 + i] * S0[i];
        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_out + start), npoints);
        }

        // Combine X blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];

            phi_tmp[i] = xc_pow[64 + i] * SX;
            phi_tmp[i] += 4.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SX;
            phi_tmp[32 + i] += 3.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SX;
            phi_tmp[64 + i] += 3.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SX;
            phi_tmp[96 + i] += 2.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SX;
            phi_tmp[128 + i] += 2.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SX;
            phi_tmp[160 + i] += 2.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SX;
            phi_tmp[192 + i] += yc_pow[32 + i] * S0[i];

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[224 + i] += yc_pow[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[256 + i] += yc[i] * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SX;
            phi_tmp[288 + i] += zc_pow[32 + i] * S0[i];

            phi_tmp[320 + i] = yc_pow[64 + i] * SX;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SX;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SX;

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SX;

            phi_tmp[448 + i] = zc_pow[64 + i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_x_out + start), npoints);
        }

        // Combine Y blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];

            phi_tmp[i] = xc_pow[64 + i] * SY;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SY;
            phi_tmp[32 + i] += xc_pow[32 + i] * S0[i];

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SY;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SY;
            phi_tmp[96 + i] += 2.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SY;
            phi_tmp[128 + i] += xc_pow[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SY;

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SY;
            phi_tmp[192 + i] += 3.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[256 + i] += xc[i] * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SY;

            phi_tmp[320 + i] = yc_pow[64 + i] * SY;
            phi_tmp[320 + i] += 4.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[352 + i] += 3.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[384 + i] += 2.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[416 + i] += zc_pow[32 + i] * S0[i];

            phi_tmp[448 + i] = zc_pow[64 + i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_y_out + start), npoints);
        }

        // Combine Z blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];

            phi_tmp[i] = xc_pow[64 + i] * SZ;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[64 + i] += xc_pow[32 + i] * S0[i];

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SZ;

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += xc_pow[i] * yc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SZ;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[224 + i] += xc[i] * yc_pow[i] * S0[i];

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[288 + i] += 3.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[320 + i] = yc_pow[64 + i] * SZ;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[352 + i] += yc_pow[32 + i] * S0[i];

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[416 + i] += 3.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[448 + i] = zc_pow[64 + i] * SZ;
            phi_tmp[448 + i] += 4.0 * zc_pow[32 + i] * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_z_out + start), npoints);
        }

        // Combine XX blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];

            phi_tmp[i] = xc_pow[64 + i] * SXX;
            phi_tmp[i] += 8.0 * xc_pow[32 + i] * SX;
            phi_tmp[i] += 12.0 * xc_pow[i] * S0[i];

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SXX;
            phi_tmp[32 + i] += 6.0 * xc_pow[i] * yc[i] * SX;
            phi_tmp[32 + i] += 6.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SXX;
            phi_tmp[64 + i] += 6.0 * xc_pow[i] * zc[i] * SX;
            phi_tmp[64 + i] += 6.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SXX;
            phi_tmp[96 + i] += 4.0 * xc[i] * yc_pow[i] * SX;
            phi_tmp[96 + i] += 2.0 * yc_pow[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SXX;
            phi_tmp[128 + i] += 4.0 * xc[i] * yc[i] * zc[i] * SX;
            phi_tmp[128 + i] += 2.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SXX;
            phi_tmp[160 + i] += 4.0 * xc[i] * zc_pow[i] * SX;
            phi_tmp[160 + i] += 2.0 * zc_pow[i] * S0[i];

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SXX;
            phi_tmp[192 + i] += 2.0 * yc_pow[32 + i] * SX;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SXX;
            phi_tmp[224 + i] += 2.0 * yc_pow[i] * zc[i] * SX;

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SXX;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc_pow[i] * SX;

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SXX;
            phi_tmp[288 + i] += 2.0 * zc_pow[32 + i] * SX;

            phi_tmp[320 + i] = yc_pow[64 + i] * SXX;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SXX;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SXX;

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SXX;

            phi_tmp[448 + i] = zc_pow[64 + i] * SXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
        }

        // Combine XY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];

            phi_tmp[i] = xc_pow[64 + i] * SXY;
            phi_tmp[i] += 4.0 * xc_pow[32 + i] * SY;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SXY;
            phi_tmp[32 + i] += xc_pow[32 + i] * SX;
            phi_tmp[32 + i] += 3.0 * xc_pow[i] * yc[i] * SY;
            phi_tmp[32 + i] += 3.0 * xc_pow[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[64 + i] += 3.0 * xc_pow[i] * zc[i] * SY;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SXY;
            phi_tmp[96 + i] += 2.0 * xc_pow[i] * yc[i] * SX;
            phi_tmp[96 + i] += 2.0 * xc[i] * yc_pow[i] * SY;
            phi_tmp[96 + i] += 4.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SXY;
            phi_tmp[128 + i] += xc_pow[i] * zc[i] * SX;
            phi_tmp[128 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SY;
            phi_tmp[128 + i] += 2.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[160 + i] += 2.0 * xc[i] * zc_pow[i] * SY;

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SXY;
            phi_tmp[192 + i] += 3.0 * xc[i] * yc_pow[i] * SX;
            phi_tmp[192 + i] += yc_pow[32 + i] * SY;
            phi_tmp[192 + i] += 3.0 * yc_pow[i] * S0[i];

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SXY;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SX;
            phi_tmp[224 + i] += yc_pow[i] * zc[i] * SY;
            phi_tmp[224 + i] += 2.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SXY;
            phi_tmp[256 + i] += xc[i] * zc_pow[i] * SX;
            phi_tmp[256 + i] += yc[i] * zc_pow[i] * SY;
            phi_tmp[256 + i] += zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[288 + i] += zc_pow[32 + i] * SY;

            phi_tmp[320 + i] = yc_pow[64 + i] * SXY;
            phi_tmp[320 + i] += 4.0 * yc_pow[32 + i] * SX;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[352 + i] += 3.0 * yc_pow[i] * zc[i] * SX;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[384 + i] += 2.0 * yc[i] * zc_pow[i] * SX;

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[416 + i] += zc_pow[32 + i] * SX;

            phi_tmp[448 + i] = zc_pow[64 + i] * SXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
        }

        // Combine XZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];

            phi_tmp[i] = xc_pow[64 + i] * SXZ;
            phi_tmp[i] += 4.0 * xc_pow[32 + i] * SZ;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SXZ;
            phi_tmp[32 + i] += 3.0 * xc_pow[i] * yc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[64 + i] += xc_pow[32 + i] * SX;
            phi_tmp[64 + i] += 3.0 * xc_pow[i] * zc[i] * SZ;
            phi_tmp[64 + i] += 3.0 * xc_pow[i] * S0[i];

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SXZ;
            phi_tmp[96 + i] += 2.0 * xc[i] * yc_pow[i] * SZ;

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SXZ;
            phi_tmp[128 + i] += xc_pow[i] * yc[i] * SX;
            phi_tmp[128 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += 2.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[i] * zc[i] * SX;
            phi_tmp[160 + i] += 2.0 * xc[i] * zc_pow[i] * SZ;
            phi_tmp[160 + i] += 4.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SXZ;
            phi_tmp[192 + i] += yc_pow[32 + i] * SZ;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SXZ;
            phi_tmp[224 + i] += xc[i] * yc_pow[i] * SX;
            phi_tmp[224 + i] += yc_pow[i] * zc[i] * SZ;
            phi_tmp[224 + i] += yc_pow[i] * S0[i];

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SX;
            phi_tmp[256 + i] += yc[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[288 + i] += 3.0 * xc[i] * zc_pow[i] * SX;
            phi_tmp[288 + i] += zc_pow[32 + i] * SZ;
            phi_tmp[288 + i] += 3.0 * zc_pow[i] * S0[i];

            phi_tmp[320 + i] = yc_pow[64 + i] * SXZ;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[352 + i] += yc_pow[32 + i] * SX;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc[i] * SX;

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[416 + i] += 3.0 * yc[i] * zc_pow[i] * SX;

            phi_tmp[448 + i] = zc_pow[64 + i] * SXZ;
            phi_tmp[448 + i] += 4.0 * zc_pow[32 + i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
        }

        // Combine YY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];

            phi_tmp[i] = xc_pow[64 + i] * SYY;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SYY;
            phi_tmp[32 + i] += 2.0 * xc_pow[32 + i] * SY;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SYY;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SYY;
            phi_tmp[96 + i] += 4.0 * xc_pow[i] * yc[i] * SY;
            phi_tmp[96 + i] += 2.0 * xc_pow[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * xc_pow[i] * zc[i] * SY;

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SYY;

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SYY;
            phi_tmp[192 + i] += 6.0 * xc[i] * yc_pow[i] * SY;
            phi_tmp[192 + i] += 6.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SYY;
            phi_tmp[224 + i] += 4.0 * xc[i] * yc[i] * zc[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SYY;
            phi_tmp[256 + i] += 2.0 * xc[i] * zc_pow[i] * SY;

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SYY;

            phi_tmp[320 + i] = yc_pow[64 + i] * SYY;
            phi_tmp[320 + i] += 8.0 * yc_pow[32 + i] * SY;
            phi_tmp[320 + i] += 12.0 * yc_pow[i] * S0[i];

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SYY;
            phi_tmp[352 + i] += 6.0 * yc_pow[i] * zc[i] * SY;
            phi_tmp[352 + i] += 6.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SYY;
            phi_tmp[384 + i] += 4.0 * yc[i] * zc_pow[i] * SY;
            phi_tmp[384 + i] += 2.0 * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SYY;
            phi_tmp[416 + i] += 2.0 * zc_pow[32 + i] * SY;

            phi_tmp[448 + i] = zc_pow[64 + i] * SYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
        }

        // Combine YZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];

            phi_tmp[i] = xc_pow[64 + i] * SYZ;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SYZ;
            phi_tmp[32 + i] += xc_pow[32 + i] * SZ;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[64 + i] += xc_pow[32 + i] * SY;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * xc_pow[i] * yc[i] * SZ;

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SYZ;
            phi_tmp[128 + i] += xc_pow[i] * yc[i] * SY;
            phi_tmp[128 + i] += xc_pow[i] * zc[i] * SZ;
            phi_tmp[128 + i] += xc_pow[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[i] * zc[i] * SY;

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SYZ;
            phi_tmp[192 + i] += 3.0 * xc[i] * yc_pow[i] * SZ;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SYZ;
            phi_tmp[224 + i] += xc[i] * yc_pow[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SZ;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SYZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SY;
            phi_tmp[256 + i] += xc[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[288 + i] += 3.0 * xc[i] * zc_pow[i] * SY;

            phi_tmp[320 + i] = yc_pow[64 + i] * SYZ;
            phi_tmp[320 + i] += 4.0 * yc_pow[32 + i] * SZ;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[352 + i] += yc_pow[32 + i] * SY;
            phi_tmp[352 + i] += 3.0 * yc_pow[i] * zc[i] * SZ;
            phi_tmp[352 + i] += 3.0 * yc_pow[i] * S0[i];

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc[i] * SY;
            phi_tmp[384 + i] += 2.0 * yc[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 4.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[416 + i] += 3.0 * yc[i] * zc_pow[i] * SY;
            phi_tmp[416 + i] += zc_pow[32 + i] * SZ;
            phi_tmp[416 + i] += 3.0 * zc_pow[i] * S0[i];

            phi_tmp[448 + i] = zc_pow[64 + i] * SYZ;
            phi_tmp[448 + i] += 4.0 * zc_pow[32 + i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
        }

        // Combine ZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];

            phi_tmp[i] = xc_pow[64 + i] * SZZ;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SZZ;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[64 + i] += 2.0 * xc_pow[32 + i] * SZ;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SZZ;

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[i] * yc[i] * SZ;

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[160 + i] += 4.0 * xc_pow[i] * zc[i] * SZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[i] * S0[i];

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SZZ;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SZZ;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc_pow[i] * SZ;

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SZZ;
            phi_tmp[256 + i] += 4.0 * xc[i] * yc[i] * zc[i] * SZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[288 + i] += 6.0 * xc[i] * zc_pow[i] * SZ;
            phi_tmp[288 + i] += 6.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[320 + i] = yc_pow[64 + i] * SZZ;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[352 + i] += 2.0 * yc_pow[32 + i] * SZ;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[384 + i] += 4.0 * yc_pow[i] * zc[i] * SZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * S0[i];

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[416 + i] += 6.0 * yc[i] * zc_pow[i] * SZ;
            phi_tmp[416 + i] += 6.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[448 + i] = zc_pow[64 + i] * SZZ;
            phi_tmp[448 + i] += 8.0 * zc_pow[32 + i] * SZ;
            phi_tmp[448 + i] += 12.0 * zc_pow[i] * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
        }

        // Combine XXX blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXX = S3[i] * xc[i] * xc[i] * xc[i] + 3 * xc[i] * S2[i];

            phi_tmp[i] = xc_pow[64 + i] * SXXX;
            phi_tmp[i] += 3.0 * 4.0 * xc_pow[32 + i] * SXX;
            phi_tmp[i] += 3.0 * 12.0 * xc_pow[i] * SX;
            phi_tmp[i] += 24.0 * xc[i] * S0[i];

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SXXX;
            phi_tmp[32 + i] += 3.0 * 3.0 * xc_pow[i] * yc[i] * SXX;
            phi_tmp[32 + i] += 3.0 * 6.0 * xc[i] * yc[i] * SX;
            phi_tmp[32 + i] += 6.0 * yc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SXXX;
            phi_tmp[64 + i] += 3.0 * 3.0 * xc_pow[i] * zc[i] * SXX;
            phi_tmp[64 + i] += 3.0 * 6.0 * xc[i] * zc[i] * SX;
            phi_tmp[64 + i] += 6.0 * zc[i] * S0[i];

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SXXX;
            phi_tmp[96 + i] += 3.0 * 2.0 * xc[i] * yc_pow[i] * SXX;
            phi_tmp[96 + i] += 3.0 * 2.0 * yc_pow[i] * SX;

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SXXX;
            phi_tmp[128 + i] += 3.0 * 2.0 * xc[i] * yc[i] * zc[i] * SXX;
            phi_tmp[128 + i] += 3.0 * 2.0 * yc[i] * zc[i] * SX;

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SXXX;
            phi_tmp[160 + i] += 3.0 * 2.0 * xc[i] * zc_pow[i] * SXX;
            phi_tmp[160 + i] += 3.0 * 2.0 * zc_pow[i] * SX;

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SXXX;
            phi_tmp[192 + i] += 3.0 * yc_pow[32 + i] * SXX;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SXXX;
            phi_tmp[224 + i] += 3.0 * yc_pow[i] * zc[i] * SXX;

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SXXX;
            phi_tmp[256 + i] += 3.0 * yc[i] * zc_pow[i] * SXX;

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SXXX;
            phi_tmp[288 + i] += 3.0 * zc_pow[32 + i] * SXX;

            phi_tmp[320 + i] = yc_pow[64 + i] * SXXX;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SXXX;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SXXX;

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SXXX;

            phi_tmp[448 + i] = zc_pow[64 + i] * SXXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
        }

        // Combine XXY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXY = S3[i] * xc[i] * xc[i] * yc[i] + S2[i] * yc[i];

            phi_tmp[i] = xc_pow[64 + i] * SXXY;
            phi_tmp[i] += 2.0 * 4.0 * xc_pow[32 + i] * SXY;
            phi_tmp[i] += 12.0 * xc_pow[i] * SY;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SXXY;
            phi_tmp[32 + i] += 2.0 * 3.0 * xc_pow[i] * yc[i] * SXY;
            phi_tmp[32 + i] += xc_pow[32 + i] * SXX;
            phi_tmp[32 + i] += 6.0 * xc[i] * yc[i] * SY;
            phi_tmp[32 + i] += 2.0 * 3.0 * xc_pow[i] * SX;
            phi_tmp[32 + i] += 6.0 * xc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SXXY;
            phi_tmp[64 + i] += 2.0 * 3.0 * xc_pow[i] * zc[i] * SXY;
            phi_tmp[64 + i] += 6.0 * xc[i] * zc[i] * SY;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SXXY;
            phi_tmp[96 + i] += 2.0 * 2.0 * xc[i] * yc_pow[i] * SXY;
            phi_tmp[96 + i] += 2.0 * xc_pow[i] * yc[i] * SXX;
            phi_tmp[96 + i] += 2.0 * yc_pow[i] * SY;
            phi_tmp[96 + i] += 2.0 * 4.0 * xc[i] * yc[i] * SX;
            phi_tmp[96 + i] += 4.0 * yc[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SXXY;
            phi_tmp[128 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc[i] * SXY;
            phi_tmp[128 + i] += xc_pow[i] * zc[i] * SXX;
            phi_tmp[128 + i] += 2.0 * yc[i] * zc[i] * SY;
            phi_tmp[128 + i] += 2.0 * 2.0 * xc[i] * zc[i] * SX;
            phi_tmp[128 + i] += 2.0 * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SXXY;
            phi_tmp[160 + i] += 2.0 * 2.0 * xc[i] * zc_pow[i] * SXY;
            phi_tmp[160 + i] += 2.0 * zc_pow[i] * SY;

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SXXY;
            phi_tmp[192 + i] += 2.0 * yc_pow[32 + i] * SXY;
            phi_tmp[192 + i] += 3.0 * xc[i] * yc_pow[i] * SXX;
            phi_tmp[192 + i] += 2.0 * 3.0 * yc_pow[i] * SX;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SXXY;
            phi_tmp[224 + i] += 2.0 * yc_pow[i] * zc[i] * SXY;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SXX;
            phi_tmp[224 + i] += 2.0 * 2.0 * yc[i] * zc[i] * SX;

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SXXY;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc_pow[i] * SXY;
            phi_tmp[256 + i] += xc[i] * zc_pow[i] * SXX;
            phi_tmp[256 + i] += 2.0 * zc_pow[i] * SX;

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SXXY;
            phi_tmp[288 + i] += 2.0 * zc_pow[32 + i] * SXY;

            phi_tmp[320 + i] = yc_pow[64 + i] * SXXY;
            phi_tmp[320 + i] += 4.0 * yc_pow[32 + i] * SXX;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SXXY;
            phi_tmp[352 + i] += 3.0 * yc_pow[i] * zc[i] * SXX;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SXXY;
            phi_tmp[384 + i] += 2.0 * yc[i] * zc_pow[i] * SXX;

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SXXY;
            phi_tmp[416 + i] += zc_pow[32 + i] * SXX;

            phi_tmp[448 + i] = zc_pow[64 + i] * SXXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
        }

        // Combine XXZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXZ = S3[i] * xc[i] * xc[i] * zc[i] + S2[i] * zc[i];

            phi_tmp[i] = xc_pow[64 + i] * SXXZ;
            phi_tmp[i] += 2.0 * 4.0 * xc_pow[32 + i] * SXZ;
            phi_tmp[i] += 12.0 * xc_pow[i] * SZ;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SXXZ;
            phi_tmp[32 + i] += 2.0 * 3.0 * xc_pow[i] * yc[i] * SXZ;
            phi_tmp[32 + i] += 6.0 * xc[i] * yc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SXXZ;
            phi_tmp[64 + i] += 2.0 * 3.0 * xc_pow[i] * zc[i] * SXZ;
            phi_tmp[64 + i] += xc_pow[32 + i] * SXX;
            phi_tmp[64 + i] += 6.0 * xc[i] * zc[i] * SZ;
            phi_tmp[64 + i] += 2.0 * 3.0 * xc_pow[i] * SX;
            phi_tmp[64 + i] += 6.0 * xc[i] * S0[i];

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SXXZ;
            phi_tmp[96 + i] += 2.0 * 2.0 * xc[i] * yc_pow[i] * SXZ;
            phi_tmp[96 + i] += 2.0 * yc_pow[i] * SZ;

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SXXZ;
            phi_tmp[128 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc[i] * SXZ;
            phi_tmp[128 + i] += xc_pow[i] * yc[i] * SXX;
            phi_tmp[128 + i] += 2.0 * yc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += 2.0 * 2.0 * xc[i] * yc[i] * SX;
            phi_tmp[128 + i] += 2.0 * yc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SXXZ;
            phi_tmp[160 + i] += 2.0 * 2.0 * xc[i] * zc_pow[i] * SXZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[i] * zc[i] * SXX;
            phi_tmp[160 + i] += 2.0 * zc_pow[i] * SZ;
            phi_tmp[160 + i] += 2.0 * 4.0 * xc[i] * zc[i] * SX;
            phi_tmp[160 + i] += 4.0 * zc[i] * S0[i];

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SXXZ;
            phi_tmp[192 + i] += 2.0 * yc_pow[32 + i] * SXZ;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SXXZ;
            phi_tmp[224 + i] += 2.0 * yc_pow[i] * zc[i] * SXZ;
            phi_tmp[224 + i] += xc[i] * yc_pow[i] * SXX;
            phi_tmp[224 + i] += 2.0 * yc_pow[i] * SX;

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SXXZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc_pow[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SXX;
            phi_tmp[256 + i] += 2.0 * 2.0 * yc[i] * zc[i] * SX;

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SXXZ;
            phi_tmp[288 + i] += 2.0 * zc_pow[32 + i] * SXZ;
            phi_tmp[288 + i] += 3.0 * xc[i] * zc_pow[i] * SXX;
            phi_tmp[288 + i] += 2.0 * 3.0 * zc_pow[i] * SX;

            phi_tmp[320 + i] = yc_pow[64 + i] * SXXZ;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SXXZ;
            phi_tmp[352 + i] += yc_pow[32 + i] * SXX;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SXXZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc[i] * SXX;

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SXXZ;
            phi_tmp[416 + i] += 3.0 * yc[i] * zc_pow[i] * SXX;

            phi_tmp[448 + i] = zc_pow[64 + i] * SXXZ;
            phi_tmp[448 + i] += 4.0 * zc_pow[32 + i] * SXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
        }

        // Combine XYY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SXYY = S3[i] * xc[i] * yc[i] * yc[i] + S2[i] * xc[i];

            phi_tmp[i] = xc_pow[64 + i] * SXYY;
            phi_tmp[i] += 4.0 * xc_pow[32 + i] * SYY;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SXYY;
            phi_tmp[32 + i] += 2.0 * xc_pow[32 + i] * SXY;
            phi_tmp[32 + i] += 3.0 * xc_pow[i] * yc[i] * SYY;
            phi_tmp[32 + i] += 2.0 * 3.0 * xc_pow[i] * SY;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SXYY;
            phi_tmp[64 + i] += 3.0 * xc_pow[i] * zc[i] * SYY;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SXYY;
            phi_tmp[96 + i] += 2.0 * 2.0 * xc_pow[i] * yc[i] * SXY;
            phi_tmp[96 + i] += 2.0 * xc[i] * yc_pow[i] * SYY;
            phi_tmp[96 + i] += 2.0 * xc_pow[i] * SX;
            phi_tmp[96 + i] += 2.0 * 4.0 * xc[i] * yc[i] * SY;
            phi_tmp[96 + i] += 4.0 * xc[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SXYY;
            phi_tmp[128 + i] += 2.0 * xc_pow[i] * zc[i] * SXY;
            phi_tmp[128 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * 2.0 * xc[i] * zc[i] * SY;

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SXYY;
            phi_tmp[160 + i] += 2.0 * xc[i] * zc_pow[i] * SYY;

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SXYY;
            phi_tmp[192 + i] += 2.0 * 3.0 * xc[i] * yc_pow[i] * SXY;
            phi_tmp[192 + i] += yc_pow[32 + i] * SYY;
            phi_tmp[192 + i] += 6.0 * xc[i] * yc[i] * SX;
            phi_tmp[192 + i] += 2.0 * 3.0 * yc_pow[i] * SY;
            phi_tmp[192 + i] += 6.0 * yc[i] * S0[i];

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SXYY;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc[i] * SXY;
            phi_tmp[224 + i] += yc_pow[i] * zc[i] * SYY;
            phi_tmp[224 + i] += 2.0 * xc[i] * zc[i] * SX;
            phi_tmp[224 + i] += 2.0 * 2.0 * yc[i] * zc[i] * SY;
            phi_tmp[224 + i] += 2.0 * zc[i] * S0[i];

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SXYY;
            phi_tmp[256 + i] += 2.0 * xc[i] * zc_pow[i] * SXY;
            phi_tmp[256 + i] += yc[i] * zc_pow[i] * SYY;
            phi_tmp[256 + i] += 2.0 * zc_pow[i] * SY;

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SXYY;
            phi_tmp[288 + i] += zc_pow[32 + i] * SYY;

            phi_tmp[320 + i] = yc_pow[64 + i] * SXYY;
            phi_tmp[320 + i] += 2.0 * 4.0 * yc_pow[32 + i] * SXY;
            phi_tmp[320 + i] += 12.0 * yc_pow[i] * SX;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SXYY;
            phi_tmp[352 + i] += 2.0 * 3.0 * yc_pow[i] * zc[i] * SXY;
            phi_tmp[352 + i] += 6.0 * yc[i] * zc[i] * SX;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SXYY;
            phi_tmp[384 + i] += 2.0 * 2.0 * yc[i] * zc_pow[i] * SXY;
            phi_tmp[384 + i] += 2.0 * zc_pow[i] * SX;

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SXYY;
            phi_tmp[416 + i] += 2.0 * zc_pow[32 + i] * SXY;

            phi_tmp[448 + i] = zc_pow[64 + i] * SXYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
        }

        // Combine XYZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SXYZ = S3[i] * xc[i] * yc[i] * zc[i];

            phi_tmp[i] = xc_pow[64 + i] * SXYZ;
            phi_tmp[i] += 4.0 * xc_pow[32 + i] * SYZ;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SXYZ;
            phi_tmp[32 + i] += 3.0 * xc_pow[i] * yc[i] * SYZ;
            phi_tmp[32 + i] += xc_pow[32 + i] * SXZ;
            phi_tmp[32 + i] += 3.0 * xc_pow[i] * SZ;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SXYZ;
            phi_tmp[64 + i] += 3.0 * xc_pow[i] * zc[i] * SYZ;
            phi_tmp[64 + i] += xc_pow[32 + i] * SXY;
            phi_tmp[64 + i] += 3.0 * xc_pow[i] * SY;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SXYZ;
            phi_tmp[96 + i] += 2.0 * xc[i] * yc_pow[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * xc_pow[i] * yc[i] * SXZ;
            phi_tmp[96 + i] += 4.0 * xc[i] * yc[i] * SZ;

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SXYZ;
            phi_tmp[128 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SYZ;
            phi_tmp[128 + i] += xc_pow[i] * zc[i] * SXZ;
            phi_tmp[128 + i] += xc_pow[i] * yc[i] * SXY;
            phi_tmp[128 + i] += 2.0 * xc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += 2.0 * xc[i] * yc[i] * SY;
            phi_tmp[128 + i] += xc_pow[i] * SX;
            phi_tmp[128 + i] += 2.0 * xc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SXYZ;
            phi_tmp[160 + i] += 2.0 * xc[i] * zc_pow[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[i] * zc[i] * SXY;
            phi_tmp[160 + i] += 4.0 * xc[i] * zc[i] * SY;

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SXYZ;
            phi_tmp[192 + i] += yc_pow[32 + i] * SYZ;
            phi_tmp[192 + i] += 3.0 * xc[i] * yc_pow[i] * SXZ;
            phi_tmp[192 + i] += 3.0 * yc_pow[i] * SZ;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SXYZ;
            phi_tmp[224 + i] += yc_pow[i] * zc[i] * SYZ;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SXZ;
            phi_tmp[224 + i] += xc[i] * yc_pow[i] * SXY;
            phi_tmp[224 + i] += 2.0 * yc[i] * zc[i] * SZ;
            phi_tmp[224 + i] += yc_pow[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc[i] * SX;
            phi_tmp[224 + i] += 4.0 * yc[i] * S0[i];

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SXYZ;
            phi_tmp[256 + i] += yc[i] * zc_pow[i] * SYZ;
            phi_tmp[256 + i] += xc[i] * zc_pow[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SXY;
            phi_tmp[256 + i] += zc_pow[i] * SZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc[i] * SY;
            phi_tmp[256 + i] += 2.0 * xc[i] * zc[i] * SX;
            phi_tmp[256 + i] += zc[i] * S0[i];

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SXYZ;
            phi_tmp[288 + i] += zc_pow[32 + i] * SYZ;
            phi_tmp[288 + i] += 3.0 * xc[i] * zc_pow[i] * SXY;
            phi_tmp[288 + i] += 3.0 * zc_pow[i] * SY;

            phi_tmp[320 + i] = yc_pow[64 + i] * SXYZ;
            phi_tmp[320 + i] += 4.0 * yc_pow[32 + i] * SXZ;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SXYZ;
            phi_tmp[352 + i] += 3.0 * yc_pow[i] * zc[i] * SXZ;
            phi_tmp[352 + i] += yc_pow[32 + i] * SXY;
            phi_tmp[352 + i] += 3.0 * yc_pow[i] * SX;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SXYZ;
            phi_tmp[384 + i] += 2.0 * yc[i] * zc_pow[i] * SXZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc[i] * SXY;
            phi_tmp[384 + i] += 4.0 * yc[i] * zc[i] * SX;

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SXYZ;
            phi_tmp[416 + i] += zc_pow[32 + i] * SXZ;
            phi_tmp[416 + i] += 3.0 * yc[i] * zc_pow[i] * SXY;
            phi_tmp[416 + i] += 3.0 * zc_pow[i] * SX;

            phi_tmp[448 + i] = zc_pow[64 + i] * SXYZ;
            phi_tmp[448 + i] += 4.0 * zc_pow[32 + i] * SXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
        }

        // Combine XZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SXZZ = S3[i] * xc[i] * zc[i] * zc[i] + S2[i] * xc[i];

            phi_tmp[i] = xc_pow[64 + i] * SXZZ;
            phi_tmp[i] += 4.0 * xc_pow[32 + i] * SZZ;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SXZZ;
            phi_tmp[32 + i] += 3.0 * xc_pow[i] * yc[i] * SZZ;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SXZZ;
            phi_tmp[64 + i] += 2.0 * xc_pow[32 + i] * SXZ;
            phi_tmp[64 + i] += 3.0 * xc_pow[i] * zc[i] * SZZ;
            phi_tmp[64 + i] += 2.0 * 3.0 * xc_pow[i] * SZ;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SXZZ;
            phi_tmp[96 + i] += 2.0 * xc[i] * yc_pow[i] * SZZ;

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SXZZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[i] * yc[i] * SXZ;
            phi_tmp[128 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * 2.0 * xc[i] * yc[i] * SZ;

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SXZZ;
            phi_tmp[160 + i] += 2.0 * 2.0 * xc_pow[i] * zc[i] * SXZ;
            phi_tmp[160 + i] += 2.0 * xc[i] * zc_pow[i] * SZZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[i] * SX;
            phi_tmp[160 + i] += 2.0 * 4.0 * xc[i] * zc[i] * SZ;
            phi_tmp[160 + i] += 4.0 * xc[i] * S0[i];

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SXZZ;
            phi_tmp[192 + i] += yc_pow[32 + i] * SZZ;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SXZZ;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc_pow[i] * SXZ;
            phi_tmp[224 + i] += yc_pow[i] * zc[i] * SZZ;
            phi_tmp[224 + i] += 2.0 * yc_pow[i] * SZ;

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SXZZ;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc[i] * SXZ;
            phi_tmp[256 + i] += yc[i] * zc_pow[i] * SZZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * SX;
            phi_tmp[256 + i] += 2.0 * 2.0 * yc[i] * zc[i] * SZ;
            phi_tmp[256 + i] += 2.0 * yc[i] * S0[i];

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SXZZ;
            phi_tmp[288 + i] += 2.0 * 3.0 * xc[i] * zc_pow[i] * SXZ;
            phi_tmp[288 + i] += zc_pow[32 + i] * SZZ;
            phi_tmp[288 + i] += 6.0 * xc[i] * zc[i] * SX;
            phi_tmp[288 + i] += 2.0 * 3.0 * zc_pow[i] * SZ;
            phi_tmp[288 + i] += 6.0 * zc[i] * S0[i];

            phi_tmp[320 + i] = yc_pow[64 + i] * SXZZ;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SXZZ;
            phi_tmp[352 + i] += 2.0 * yc_pow[32 + i] * SXZ;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SXZZ;
            phi_tmp[384 + i] += 2.0 * 2.0 * yc_pow[i] * zc[i] * SXZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * SX;

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SXZZ;
            phi_tmp[416 + i] += 2.0 * 3.0 * yc[i] * zc_pow[i] * SXZ;
            phi_tmp[416 + i] += 6.0 * yc[i] * zc[i] * SX;

            phi_tmp[448 + i] = zc_pow[64 + i] * SXZZ;
            phi_tmp[448 + i] += 2.0 * 4.0 * zc_pow[32 + i] * SXZ;
            phi_tmp[448 + i] += 12.0 * zc_pow[i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
        }

        // Combine YYY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SYYY = S3[i] * yc[i] * yc[i] * yc[i] + 3 * yc[i] * S2[i];

            phi_tmp[i] = xc_pow[64 + i] * SYYY;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SYYY;
            phi_tmp[32 + i] += 3.0 * xc_pow[32 + i] * SYY;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SYYY;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SYYY;
            phi_tmp[96 + i] += 3.0 * 2.0 * xc_pow[i] * yc[i] * SYY;
            phi_tmp[96 + i] += 3.0 * 2.0 * xc_pow[i] * SY;

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SYYY;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * zc[i] * SYY;

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SYYY;

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SYYY;
            phi_tmp[192 + i] += 3.0 * 3.0 * xc[i] * yc_pow[i] * SYY;
            phi_tmp[192 + i] += 3.0 * 6.0 * xc[i] * yc[i] * SY;
            phi_tmp[192 + i] += 6.0 * xc[i] * S0[i];

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SYYY;
            phi_tmp[224 + i] += 3.0 * 2.0 * xc[i] * yc[i] * zc[i] * SYY;
            phi_tmp[224 + i] += 3.0 * 2.0 * xc[i] * zc[i] * SY;

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SYYY;
            phi_tmp[256 + i] += 3.0 * xc[i] * zc_pow[i] * SYY;

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SYYY;

            phi_tmp[320 + i] = yc_pow[64 + i] * SYYY;
            phi_tmp[320 + i] += 3.0 * 4.0 * yc_pow[32 + i] * SYY;
            phi_tmp[320 + i] += 3.0 * 12.0 * yc_pow[i] * SY;
            phi_tmp[320 + i] += 24.0 * yc[i] * S0[i];

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SYYY;
            phi_tmp[352 + i] += 3.0 * 3.0 * yc_pow[i] * zc[i] * SYY;
            phi_tmp[352 + i] += 3.0 * 6.0 * yc[i] * zc[i] * SY;
            phi_tmp[352 + i] += 6.0 * zc[i] * S0[i];

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SYYY;
            phi_tmp[384 + i] += 3.0 * 2.0 * yc[i] * zc_pow[i] * SYY;
            phi_tmp[384 + i] += 3.0 * 2.0 * zc_pow[i] * SY;

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SYYY;
            phi_tmp[416 + i] += 3.0 * zc_pow[32 + i] * SYY;

            phi_tmp[448 + i] = zc_pow[64 + i] * SYYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
        }

        // Combine YYZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SYYZ = S3[i] * yc[i] * yc[i] * zc[i] + S2[i] * zc[i];

            phi_tmp[i] = xc_pow[64 + i] * SYYZ;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SYYZ;
            phi_tmp[32 + i] += 2.0 * xc_pow[32 + i] * SYZ;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SYYZ;
            phi_tmp[64 + i] += xc_pow[32 + i] * SYY;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SYYZ;
            phi_tmp[96 + i] += 2.0 * 2.0 * xc_pow[i] * yc[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * xc_pow[i] * SZ;

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SYYZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[i] * zc[i] * SYZ;
            phi_tmp[128 + i] += xc_pow[i] * yc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * xc_pow[i] * SY;

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SYYZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[i] * zc[i] * SYY;

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SYYZ;
            phi_tmp[192 + i] += 2.0 * 3.0 * xc[i] * yc_pow[i] * SYZ;
            phi_tmp[192 + i] += 6.0 * xc[i] * yc[i] * SZ;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SYYZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc[i] * SYZ;
            phi_tmp[224 + i] += xc[i] * yc_pow[i] * SYY;
            phi_tmp[224 + i] += 2.0 * xc[i] * zc[i] * SZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc[i] * yc[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc[i] * S0[i];

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SYYZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * zc_pow[i] * SYZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SYY;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc[i] * zc[i] * SY;

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SYYZ;
            phi_tmp[288 + i] += 3.0 * xc[i] * zc_pow[i] * SYY;

            phi_tmp[320 + i] = yc_pow[64 + i] * SYYZ;
            phi_tmp[320 + i] += 2.0 * 4.0 * yc_pow[32 + i] * SYZ;
            phi_tmp[320 + i] += 12.0 * yc_pow[i] * SZ;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SYYZ;
            phi_tmp[352 + i] += 2.0 * 3.0 * yc_pow[i] * zc[i] * SYZ;
            phi_tmp[352 + i] += yc_pow[32 + i] * SYY;
            phi_tmp[352 + i] += 6.0 * yc[i] * zc[i] * SZ;
            phi_tmp[352 + i] += 2.0 * 3.0 * yc_pow[i] * SY;
            phi_tmp[352 + i] += 6.0 * yc[i] * S0[i];

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SYYZ;
            phi_tmp[384 + i] += 2.0 * 2.0 * yc[i] * zc_pow[i] * SYZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc[i] * SYY;
            phi_tmp[384 + i] += 2.0 * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 2.0 * 4.0 * yc[i] * zc[i] * SY;
            phi_tmp[384 + i] += 4.0 * zc[i] * S0[i];

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SYYZ;
            phi_tmp[416 + i] += 2.0 * zc_pow[32 + i] * SYZ;
            phi_tmp[416 + i] += 3.0 * yc[i] * zc_pow[i] * SYY;
            phi_tmp[416 + i] += 2.0 * 3.0 * zc_pow[i] * SY;

            phi_tmp[448 + i] = zc_pow[64 + i] * SYYZ;
            phi_tmp[448 + i] += 4.0 * zc_pow[32 + i] * SYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
        }

        // Combine YZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SYZZ = S3[i] * yc[i] * zc[i] * zc[i] + S2[i] * yc[i];

            phi_tmp[i] = xc_pow[64 + i] * SYZZ;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SYZZ;
            phi_tmp[32 + i] += xc_pow[32 + i] * SZZ;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SYZZ;
            phi_tmp[64 + i] += 2.0 * xc_pow[32 + i] * SYZ;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SYZZ;
            phi_tmp[96 + i] += 2.0 * xc_pow[i] * yc[i] * SZZ;

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SYZZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[i] * yc[i] * SYZ;
            phi_tmp[128 + i] += xc_pow[i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[i] * SZ;

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SYZZ;
            phi_tmp[160 + i] += 2.0 * 2.0 * xc_pow[i] * zc[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[i] * SY;

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SYZZ;
            phi_tmp[192 + i] += 3.0 * xc[i] * yc_pow[i] * SZZ;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SYZZ;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc_pow[i] * SYZ;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc[i] * zc[i] * SZZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc[i] * yc[i] * SZ;

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SYZZ;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc[i] * SYZ;
            phi_tmp[256 + i] += xc[i] * zc_pow[i] * SZZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * SY;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc[i] * zc[i] * SZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * S0[i];

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SYZZ;
            phi_tmp[288 + i] += 2.0 * 3.0 * xc[i] * zc_pow[i] * SYZ;
            phi_tmp[288 + i] += 6.0 * xc[i] * zc[i] * SY;

            phi_tmp[320 + i] = yc_pow[64 + i] * SYZZ;
            phi_tmp[320 + i] += 4.0 * yc_pow[32 + i] * SZZ;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SYZZ;
            phi_tmp[352 + i] += 2.0 * yc_pow[32 + i] * SYZ;
            phi_tmp[352 + i] += 3.0 * yc_pow[i] * zc[i] * SZZ;
            phi_tmp[352 + i] += 2.0 * 3.0 * yc_pow[i] * SZ;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SYZZ;
            phi_tmp[384 + i] += 2.0 * 2.0 * yc_pow[i] * zc[i] * SYZ;
            phi_tmp[384 + i] += 2.0 * yc[i] * zc_pow[i] * SZZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * SY;
            phi_tmp[384 + i] += 2.0 * 4.0 * yc[i] * zc[i] * SZ;
            phi_tmp[384 + i] += 4.0 * yc[i] * S0[i];

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SYZZ;
            phi_tmp[416 + i] += 2.0 * 3.0 * yc[i] * zc_pow[i] * SYZ;
            phi_tmp[416 + i] += zc_pow[32 + i] * SZZ;
            phi_tmp[416 + i] += 6.0 * yc[i] * zc[i] * SY;
            phi_tmp[416 + i] += 2.0 * 3.0 * zc_pow[i] * SZ;
            phi_tmp[416 + i] += 6.0 * zc[i] * S0[i];

            phi_tmp[448 + i] = zc_pow[64 + i] * SYZZ;
            phi_tmp[448 + i] += 2.0 * 4.0 * zc_pow[32 + i] * SYZ;
            phi_tmp[448 + i] += 12.0 * zc_pow[i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
        }

        // Combine ZZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SZZZ = S3[i] * zc[i] * zc[i] * zc[i] + 3 * zc[i] * S2[i];

            phi_tmp[i] = xc_pow[64 + i] * SZZZ;

            phi_tmp[32 + i] = xc_pow[32 + i] * yc[i] * SZZZ;

            phi_tmp[64 + i] = xc_pow[32 + i] * zc[i] * SZZZ;
            phi_tmp[64 + i] += 3.0 * xc_pow[32 + i] * SZZ;

            phi_tmp[96 + i] = xc_pow[i] * yc_pow[i] * SZZZ;

            phi_tmp[128 + i] = xc_pow[i] * yc[i] * zc[i] * SZZZ;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * yc[i] * SZZ;

            phi_tmp[160 + i] = xc_pow[i] * zc_pow[i] * SZZZ;
            phi_tmp[160 + i] += 3.0 * 2.0 * xc_pow[i] * zc[i] * SZZ;
            phi_tmp[160 + i] += 3.0 * 2.0 * xc_pow[i] * SZ;

            phi_tmp[192 + i] = xc[i] * yc_pow[32 + i] * SZZZ;

            phi_tmp[224 + i] = xc[i] * yc_pow[i] * zc[i] * SZZZ;
            phi_tmp[224 + i] += 3.0 * xc[i] * yc_pow[i] * SZZ;

            phi_tmp[256 + i] = xc[i] * yc[i] * zc_pow[i] * SZZZ;
            phi_tmp[256 + i] += 3.0 * 2.0 * xc[i] * yc[i] * zc[i] * SZZ;
            phi_tmp[256 + i] += 3.0 * 2.0 * xc[i] * yc[i] * SZ;

            phi_tmp[288 + i] = xc[i] * zc_pow[32 + i] * SZZZ;
            phi_tmp[288 + i] += 3.0 * 3.0 * xc[i] * zc_pow[i] * SZZ;
            phi_tmp[288 + i] += 3.0 * 6.0 * xc[i] * zc[i] * SZ;
            phi_tmp[288 + i] += 6.0 * xc[i] * S0[i];

            phi_tmp[320 + i] = yc_pow[64 + i] * SZZZ;

            phi_tmp[352 + i] = yc_pow[32 + i] * zc[i] * SZZZ;
            phi_tmp[352 + i] += 3.0 * yc_pow[32 + i] * SZZ;

            phi_tmp[384 + i] = yc_pow[i] * zc_pow[i] * SZZZ;
            phi_tmp[384 + i] += 3.0 * 2.0 * yc_pow[i] * zc[i] * SZZ;
            phi_tmp[384 + i] += 3.0 * 2.0 * yc_pow[i] * SZ;

            phi_tmp[416 + i] = yc[i] * zc_pow[32 + i] * SZZZ;
            phi_tmp[416 + i] += 3.0 * 3.0 * yc[i] * zc_pow[i] * SZZ;
            phi_tmp[416 + i] += 3.0 * 6.0 * yc[i] * zc[i] * SZ;
            phi_tmp[416 + i] += 6.0 * yc[i] * S0[i];

            phi_tmp[448 + i] = zc_pow[64 + i] * SZZZ;
            phi_tmp[448 + i] += 3.0 * 4.0 * zc_pow[32 + i] * SZZ;
            phi_tmp[448 + i] += 3.0 * 12.0 * zc_pow[i] * SZ;
            phi_tmp[448 + i] += 24.0 * zc[i] * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
        }

    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);
    ALIGNED_FREE(expn2);

    // Free Power temporaries
    ALIGNED_FREE(xc_pow);
    ALIGNED_FREE(yc_pow);
    ALIGNED_FREE(zc_pow);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);

}

void gg_collocation_L5_deriv3(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out, double* PRAGMA_RESTRICT phi_xxx_out, double* PRAGMA_RESTRICT phi_xxy_out, double* PRAGMA_RESTRICT phi_xxz_out, double* PRAGMA_RESTRICT phi_xyy_out, double* PRAGMA_RESTRICT phi_xyz_out, double* PRAGMA_RESTRICT phi_xzz_out, double* PRAGMA_RESTRICT phi_yyy_out, double* PRAGMA_RESTRICT phi_yyz_out, double* PRAGMA_RESTRICT phi_yzz_out, double* PRAGMA_RESTRICT phi_zzz_out) {

    // Sizing
    unsigned long nblocks = npoints / 32;
    nblocks += (npoints % 32) ? 1 : 0;
    const unsigned long ncart = 21;
    const unsigned long nspherical = 11;
    unsigned long nout;

    if ((order == GG_SPHERICAL_CCA) || (order == GG_SPHERICAL_GAUSSIAN)) {
        nout = nspherical;
        } else {
        nout = ncart;
    }

    // Allocate S temporaries, single block to stay on cache
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 288 * sizeof(double));
    double* PRAGMA_RESTRICT xc = cache_data + 0;
    ASSUME_ALIGNED(xc, 64);
    double* PRAGMA_RESTRICT yc = cache_data + 32;
    ASSUME_ALIGNED(yc, 64);
    double* PRAGMA_RESTRICT zc = cache_data + 64;
    ASSUME_ALIGNED(zc, 64);
    double* PRAGMA_RESTRICT R2 = cache_data + 96;
    ASSUME_ALIGNED(R2, 64);
    double* PRAGMA_RESTRICT S0 = cache_data + 128;
    ASSUME_ALIGNED(S0, 64);
    double* PRAGMA_RESTRICT tmp1 = cache_data + 160;
    ASSUME_ALIGNED(tmp1, 64);
    double* PRAGMA_RESTRICT S1 = cache_data + 192;
    ASSUME_ALIGNED(S1, 64);
    double* PRAGMA_RESTRICT S2 = cache_data + 224;
    ASSUME_ALIGNED(S2, 64);
    double* PRAGMA_RESTRICT S3 = cache_data + 256;
    ASSUME_ALIGNED(S3, 64);

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));
    double* PRAGMA_RESTRICT expn2 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));

    // Allocate power temporaries
    double* PRAGMA_RESTRICT xc_pow = (double*)ALIGNED_MALLOC(64, 128 * sizeof(double));
    ASSUME_ALIGNED(xc_pow, 64);
    double* PRAGMA_RESTRICT yc_pow = (double*)ALIGNED_MALLOC(64, 128 * sizeof(double));
    ASSUME_ALIGNED(yc_pow, 64);
    double* PRAGMA_RESTRICT zc_pow = (double*)ALIGNED_MALLOC(64, 128 * sizeof(double));
    ASSUME_ALIGNED(zc_pow, 64);

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, 672 * sizeof(double));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;
    double AXX, AXY, AXZ, AYY, AYZ, AZZ;
    double AXXX, XXY, XXZ, XYY, XYZ, XZZ, YYY, YYZ, YZZ, ZZZ;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
        expn2[i] = -2.0 * exponents[i];
    }

    // Start outer block loop
    for (unsigned long block = 0; block < nblocks; block++) {


        // Copy data into inner temps
        const unsigned long start = block * 32;
        const unsigned long remain = ((start + 32) > npoints) ? (npoints - start) : 32;

        // Handle non-AM dependant temps
        if (xyz_stride == 1) {
            const double* PRAGMA_RESTRICT x = xyz + start;
            const double* PRAGMA_RESTRICT y = xyz + npoints + start;
            const double* PRAGMA_RESTRICT z = xyz + 2 * npoints + start;
            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = x[i] - center_x;
                yc[i] = y[i] - center_y;
                zc[i] = z[i] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
            } else {
            unsigned int start_shift = start * xyz_stride;

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = xyz[start_shift + i * xyz_stride] - center_x;
                yc[i] = xyz[start_shift + i * xyz_stride + 1] - center_y;
                zc[i] = xyz[start_shift + i * xyz_stride + 2] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];
            const double alpha_n2 = expn2[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
                const double T2 = alpha_n2 * T1;
                S1[i] += T2;
                const double T3 = alpha_n2 * T2;
                S2[i] += T3;
                const double T4 = alpha_n2 * T3;
                S3[i] += T4;
            }

        }

        // Build powers
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {

            // Cartesian derivs
            xc_pow[i] = xc[i] * xc[i];
            yc_pow[i] = yc[i] * yc[i];
            zc_pow[i] = zc[i] * zc[i];
            xc_pow[32 + i] = xc_pow[i] * xc[i];
            yc_pow[32 + i] = yc_pow[i] * yc[i];
            zc_pow[32 + i] = zc_pow[i] * zc[i];
            xc_pow[64 + i] = xc_pow[32 + i] * xc[i];
            yc_pow[64 + i] = yc_pow[32 + i] * yc[i];
            zc_pow[64 + i] = zc_pow[32 + i] * zc[i];
            xc_pow[96 + i] = xc_pow[64 + i] * xc[i];
            yc_pow[96 + i] = yc_pow[64 + i] * yc[i];
            zc_pow[96 + i] = zc_pow[64 + i] * zc[i];
        }
        // Combine A blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {

            phi_tmp[i] = xc_pow[96 + i] * S0[i];
            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * S0[i];
            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * S0[i];
            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * S0[i];
            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * S0[i];
            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * S0[i];
            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * S0[i];
            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * S0[i];
            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * S0[i];
            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * S0[i];
            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * S0[i];
            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * S0[i];
            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * S0[i];
            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * S0[i];
            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * S0[i];
            phi_tmp[480 + i] = yc_pow[96 + i] * S0[i];
            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * S0[i];
            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * S0[i];
            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * S0[i];
            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * S0[i];
            phi_tmp[640 + i] = zc_pow[96 + i] * S0[i];
        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_out + start), npoints);
        }

        // Combine X blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];

            phi_tmp[i] = xc_pow[96 + i] * SX;
            phi_tmp[i] += 5.0 * xc_pow[64 + i] * S0[i];

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SX;
            phi_tmp[32 + i] += 4.0 * xc_pow[32 + i] * yc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SX;
            phi_tmp[64 + i] += 4.0 * xc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SX;
            phi_tmp[96 + i] += 3.0 * xc_pow[i] * yc_pow[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SX;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SX;
            phi_tmp[160 + i] += 3.0 * xc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SX;
            phi_tmp[192 + i] += 2.0 * xc[i] * yc_pow[32 + i] * S0[i];

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SX;
            phi_tmp[288 + i] += 2.0 * xc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SX;
            phi_tmp[320 + i] += yc_pow[64 + i] * S0[i];

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SX;
            phi_tmp[352 + i] += yc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SX;
            phi_tmp[384 + i] += yc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SX;
            phi_tmp[416 + i] += yc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SX;
            phi_tmp[448 + i] += zc_pow[64 + i] * S0[i];

            phi_tmp[480 + i] = yc_pow[96 + i] * SX;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SX;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SX;

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SX;

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SX;

            phi_tmp[640 + i] = zc_pow[96 + i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_x_out + start), npoints);
        }

        // Combine Y blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];

            phi_tmp[i] = xc_pow[96 + i] * SY;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SY;
            phi_tmp[32 + i] += xc_pow[64 + i] * S0[i];

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SY;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SY;
            phi_tmp[96 + i] += 2.0 * xc_pow[32 + i] * yc[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SY;
            phi_tmp[128 + i] += xc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SY;

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SY;
            phi_tmp[192 + i] += 3.0 * xc_pow[i] * yc_pow[i] * S0[i];

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[256 + i] += xc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SY;

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SY;
            phi_tmp[320 + i] += 4.0 * xc[i] * yc_pow[32 + i] * S0[i];

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[352 + i] += 3.0 * xc[i] * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[416 + i] += xc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SY;

            phi_tmp[480 + i] = yc_pow[96 + i] * SY;
            phi_tmp[480 + i] += 5.0 * yc_pow[64 + i] * S0[i];

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SY;
            phi_tmp[512 + i] += 4.0 * yc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SY;
            phi_tmp[544 + i] += 3.0 * yc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SY;
            phi_tmp[576 + i] += 2.0 * yc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SY;
            phi_tmp[608 + i] += zc_pow[64 + i] * S0[i];

            phi_tmp[640 + i] = zc_pow[96 + i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_y_out + start), npoints);
        }

        // Combine Z blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];

            phi_tmp[i] = xc_pow[96 + i] * SZ;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SZ;
            phi_tmp[64 + i] += xc_pow[64 + i] * S0[i];

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SZ;

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += xc_pow[32 + i] * yc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SZ;

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[224 + i] += xc_pow[i] * yc_pow[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SZ;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[352 + i] += xc[i] * yc_pow[32 + i] * S0[i];

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[416 + i] += 3.0 * xc[i] * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SZ;
            phi_tmp[448 + i] += 4.0 * xc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[480 + i] = yc_pow[96 + i] * SZ;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SZ;
            phi_tmp[512 + i] += yc_pow[64 + i] * S0[i];

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SZ;
            phi_tmp[576 + i] += 3.0 * yc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SZ;
            phi_tmp[608 + i] += 4.0 * yc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[640 + i] = zc_pow[96 + i] * SZ;
            phi_tmp[640 + i] += 5.0 * zc_pow[64 + i] * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_z_out + start), npoints);
        }

        // Combine XX blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];

            phi_tmp[i] = xc_pow[96 + i] * SXX;
            phi_tmp[i] += 10.0 * xc_pow[64 + i] * SX;
            phi_tmp[i] += 20.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SXX;
            phi_tmp[32 + i] += 8.0 * xc_pow[32 + i] * yc[i] * SX;
            phi_tmp[32 + i] += 12.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SXX;
            phi_tmp[64 + i] += 8.0 * xc_pow[32 + i] * zc[i] * SX;
            phi_tmp[64 + i] += 12.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SXX;
            phi_tmp[96 + i] += 6.0 * xc_pow[i] * yc_pow[i] * SX;
            phi_tmp[96 + i] += 6.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SXX;
            phi_tmp[128 + i] += 6.0 * xc_pow[i] * yc[i] * zc[i] * SX;
            phi_tmp[128 + i] += 6.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SXX;
            phi_tmp[160 + i] += 6.0 * xc_pow[i] * zc_pow[i] * SX;
            phi_tmp[160 + i] += 6.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SXX;
            phi_tmp[192 + i] += 4.0 * xc[i] * yc_pow[32 + i] * SX;
            phi_tmp[192 + i] += 2.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SXX;
            phi_tmp[224 + i] += 4.0 * xc[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[224 + i] += 2.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SXX;
            phi_tmp[256 + i] += 4.0 * xc[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SXX;
            phi_tmp[288 + i] += 4.0 * xc[i] * zc_pow[32 + i] * SX;
            phi_tmp[288 + i] += 2.0 * zc_pow[32 + i] * S0[i];

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SXX;
            phi_tmp[320 + i] += 2.0 * yc_pow[64 + i] * SX;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SXX;
            phi_tmp[352 + i] += 2.0 * yc_pow[32 + i] * zc[i] * SX;

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SXX;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc_pow[i] * SX;

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SXX;
            phi_tmp[416 + i] += 2.0 * yc[i] * zc_pow[32 + i] * SX;

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SXX;
            phi_tmp[448 + i] += 2.0 * zc_pow[64 + i] * SX;

            phi_tmp[480 + i] = yc_pow[96 + i] * SXX;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SXX;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SXX;

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SXX;

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SXX;

            phi_tmp[640 + i] = zc_pow[96 + i] * SXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
        }

        // Combine XY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];

            phi_tmp[i] = xc_pow[96 + i] * SXY;
            phi_tmp[i] += 5.0 * xc_pow[64 + i] * SY;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SXY;
            phi_tmp[32 + i] += xc_pow[64 + i] * SX;
            phi_tmp[32 + i] += 4.0 * xc_pow[32 + i] * yc[i] * SY;
            phi_tmp[32 + i] += 4.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SXY;
            phi_tmp[64 + i] += 4.0 * xc_pow[32 + i] * zc[i] * SY;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SXY;
            phi_tmp[96 + i] += 2.0 * xc_pow[32 + i] * yc[i] * SX;
            phi_tmp[96 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SY;
            phi_tmp[96 + i] += 6.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SXY;
            phi_tmp[128 + i] += xc_pow[32 + i] * zc[i] * SX;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * yc[i] * zc[i] * SY;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SXY;
            phi_tmp[160 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SY;

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SXY;
            phi_tmp[192 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SX;
            phi_tmp[192 + i] += 2.0 * xc[i] * yc_pow[32 + i] * SY;
            phi_tmp[192 + i] += 6.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SXY;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * yc[i] * zc[i] * SX;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[224 + i] += 4.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SXY;
            phi_tmp[256 + i] += xc_pow[i] * zc_pow[i] * SX;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[256 + i] += 2.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SXY;
            phi_tmp[288 + i] += 2.0 * xc[i] * zc_pow[32 + i] * SY;

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SXY;
            phi_tmp[320 + i] += 4.0 * xc[i] * yc_pow[32 + i] * SX;
            phi_tmp[320 + i] += yc_pow[64 + i] * SY;
            phi_tmp[320 + i] += 4.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[352 + i] += 3.0 * xc[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[352 + i] += yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[352 + i] += 3.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[384 + i] += yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[384 + i] += 2.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[416 + i] += xc[i] * zc_pow[32 + i] * SX;
            phi_tmp[416 + i] += yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[416 + i] += zc_pow[32 + i] * S0[i];

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SXY;
            phi_tmp[448 + i] += zc_pow[64 + i] * SY;

            phi_tmp[480 + i] = yc_pow[96 + i] * SXY;
            phi_tmp[480 + i] += 5.0 * yc_pow[64 + i] * SX;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SXY;
            phi_tmp[512 + i] += 4.0 * yc_pow[32 + i] * zc[i] * SX;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SXY;
            phi_tmp[544 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SX;

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SXY;
            phi_tmp[576 + i] += 2.0 * yc[i] * zc_pow[32 + i] * SX;

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SXY;
            phi_tmp[608 + i] += zc_pow[64 + i] * SX;

            phi_tmp[640 + i] = zc_pow[96 + i] * SXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
        }

        // Combine XZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];

            phi_tmp[i] = xc_pow[96 + i] * SXZ;
            phi_tmp[i] += 5.0 * xc_pow[64 + i] * SZ;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SXZ;
            phi_tmp[32 + i] += 4.0 * xc_pow[32 + i] * yc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SXZ;
            phi_tmp[64 + i] += xc_pow[64 + i] * SX;
            phi_tmp[64 + i] += 4.0 * xc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[64 + i] += 4.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SXZ;
            phi_tmp[96 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SZ;

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SXZ;
            phi_tmp[128 + i] += xc_pow[32 + i] * yc[i] * SX;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * yc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SXZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[32 + i] * zc[i] * SX;
            phi_tmp[160 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[160 + i] += 6.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SXZ;
            phi_tmp[192 + i] += 2.0 * xc[i] * yc_pow[32 + i] * SZ;

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SXZ;
            phi_tmp[224 + i] += xc_pow[i] * yc_pow[i] * SX;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * yc[i] * zc[i] * SX;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 4.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SX;
            phi_tmp[288 + i] += 2.0 * xc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[288 + i] += 6.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SXZ;
            phi_tmp[320 + i] += yc_pow[64 + i] * SZ;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[352 + i] += xc[i] * yc_pow[32 + i] * SX;
            phi_tmp[352 + i] += yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[352 + i] += yc_pow[32 + i] * S0[i];

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[384 + i] += yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[416 + i] += 3.0 * xc[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[416 + i] += yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[416 + i] += 3.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SXZ;
            phi_tmp[448 + i] += 4.0 * xc[i] * zc_pow[32 + i] * SX;
            phi_tmp[448 + i] += zc_pow[64 + i] * SZ;
            phi_tmp[448 + i] += 4.0 * zc_pow[32 + i] * S0[i];

            phi_tmp[480 + i] = yc_pow[96 + i] * SXZ;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SXZ;
            phi_tmp[512 + i] += yc_pow[64 + i] * SX;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SXZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * zc[i] * SX;

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[576 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SX;

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SXZ;
            phi_tmp[608 + i] += 4.0 * yc[i] * zc_pow[32 + i] * SX;

            phi_tmp[640 + i] = zc_pow[96 + i] * SXZ;
            phi_tmp[640 + i] += 5.0 * zc_pow[64 + i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
        }

        // Combine YY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];

            phi_tmp[i] = xc_pow[96 + i] * SYY;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SYY;
            phi_tmp[32 + i] += 2.0 * xc_pow[64 + i] * SY;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SYY;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SYY;
            phi_tmp[96 + i] += 4.0 * xc_pow[32 + i] * yc[i] * SY;
            phi_tmp[96 + i] += 2.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * xc_pow[32 + i] * zc[i] * SY;

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SYY;

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SYY;
            phi_tmp[192 + i] += 6.0 * xc_pow[i] * yc_pow[i] * SY;
            phi_tmp[192 + i] += 6.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SYY;
            phi_tmp[224 + i] += 4.0 * xc_pow[i] * yc[i] * zc[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SYY;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * zc_pow[i] * SY;

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SYY;

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SYY;
            phi_tmp[320 + i] += 8.0 * xc[i] * yc_pow[32 + i] * SY;
            phi_tmp[320 + i] += 12.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SYY;
            phi_tmp[352 + i] += 6.0 * xc[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[352 + i] += 6.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SYY;
            phi_tmp[384 + i] += 4.0 * xc[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[384 + i] += 2.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SYY;
            phi_tmp[416 + i] += 2.0 * xc[i] * zc_pow[32 + i] * SY;

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SYY;

            phi_tmp[480 + i] = yc_pow[96 + i] * SYY;
            phi_tmp[480 + i] += 10.0 * yc_pow[64 + i] * SY;
            phi_tmp[480 + i] += 20.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SYY;
            phi_tmp[512 + i] += 8.0 * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[512 + i] += 12.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SYY;
            phi_tmp[544 + i] += 6.0 * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[544 + i] += 6.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SYY;
            phi_tmp[576 + i] += 4.0 * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[576 + i] += 2.0 * zc_pow[32 + i] * S0[i];

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SYY;
            phi_tmp[608 + i] += 2.0 * zc_pow[64 + i] * SY;

            phi_tmp[640 + i] = zc_pow[96 + i] * SYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
        }

        // Combine YZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];

            phi_tmp[i] = xc_pow[96 + i] * SYZ;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SYZ;
            phi_tmp[32 + i] += xc_pow[64 + i] * SZ;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SYZ;
            phi_tmp[64 + i] += xc_pow[64 + i] * SY;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * xc_pow[32 + i] * yc[i] * SZ;

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SYZ;
            phi_tmp[128 + i] += xc_pow[32 + i] * yc[i] * SY;
            phi_tmp[128 + i] += xc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[128 + i] += xc_pow[32 + i] * S0[i];

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[32 + i] * zc[i] * SY;

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SYZ;
            phi_tmp[192 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SZ;

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SYZ;
            phi_tmp[224 + i] += xc_pow[i] * yc_pow[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * yc[i] * zc[i] * SZ;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SYZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * yc[i] * zc[i] * SY;
            phi_tmp[256 + i] += xc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SY;

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SYZ;
            phi_tmp[320 + i] += 4.0 * xc[i] * yc_pow[32 + i] * SZ;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[352 + i] += xc[i] * yc_pow[32 + i] * SY;
            phi_tmp[352 + i] += 3.0 * xc[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[352 + i] += 3.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 4.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[416 + i] += 3.0 * xc[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[416 + i] += xc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[416 + i] += 3.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SYZ;
            phi_tmp[448 + i] += 4.0 * xc[i] * zc_pow[32 + i] * SY;

            phi_tmp[480 + i] = yc_pow[96 + i] * SYZ;
            phi_tmp[480 + i] += 5.0 * yc_pow[64 + i] * SZ;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SYZ;
            phi_tmp[512 + i] += yc_pow[64 + i] * SY;
            phi_tmp[512 + i] += 4.0 * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[512 + i] += 4.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SYZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[544 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[544 + i] += 6.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[576 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[576 + i] += 2.0 * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[576 + i] += 6.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SYZ;
            phi_tmp[608 + i] += 4.0 * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[608 + i] += zc_pow[64 + i] * SZ;
            phi_tmp[608 + i] += 4.0 * zc_pow[32 + i] * S0[i];

            phi_tmp[640 + i] = zc_pow[96 + i] * SYZ;
            phi_tmp[640 + i] += 5.0 * zc_pow[64 + i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
        }

        // Combine ZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];

            phi_tmp[i] = xc_pow[96 + i] * SZZ;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SZZ;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SZZ;
            phi_tmp[64 + i] += 2.0 * xc_pow[64 + i] * SZ;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SZZ;

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[32 + i] * yc[i] * SZ;

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SZZ;
            phi_tmp[160 + i] += 4.0 * xc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SZZ;

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SZZ;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * yc_pow[i] * SZ;

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SZZ;
            phi_tmp[256 + i] += 4.0 * xc_pow[i] * yc[i] * zc[i] * SZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[288 + i] += 6.0 * xc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[288 + i] += 6.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SZZ;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[352 + i] += 2.0 * xc[i] * yc_pow[32 + i] * SZ;

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[384 + i] += 4.0 * xc[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[416 + i] += 6.0 * xc[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[416 + i] += 6.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SZZ;
            phi_tmp[448 + i] += 8.0 * xc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[448 + i] += 12.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[480 + i] = yc_pow[96 + i] * SZZ;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SZZ;
            phi_tmp[512 + i] += 2.0 * yc_pow[64 + i] * SZ;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SZZ;
            phi_tmp[544 + i] += 4.0 * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[576 + i] += 6.0 * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[576 + i] += 6.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SZZ;
            phi_tmp[608 + i] += 8.0 * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[608 + i] += 12.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[640 + i] = zc_pow[96 + i] * SZZ;
            phi_tmp[640 + i] += 10.0 * zc_pow[64 + i] * SZ;
            phi_tmp[640 + i] += 20.0 * zc_pow[32 + i] * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
        }

        // Combine XXX blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXX = S3[i] * xc[i] * xc[i] * xc[i] + 3 * xc[i] * S2[i];

            phi_tmp[i] = xc_pow[96 + i] * SXXX;
            phi_tmp[i] += 3.0 * 5.0 * xc_pow[64 + i] * SXX;
            phi_tmp[i] += 3.0 * 20.0 * xc_pow[32 + i] * SX;
            phi_tmp[i] += 60.0 * xc_pow[i] * S0[i];

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SXXX;
            phi_tmp[32 + i] += 3.0 * 4.0 * xc_pow[32 + i] * yc[i] * SXX;
            phi_tmp[32 + i] += 3.0 * 12.0 * xc_pow[i] * yc[i] * SX;
            phi_tmp[32 + i] += 24.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SXXX;
            phi_tmp[64 + i] += 3.0 * 4.0 * xc_pow[32 + i] * zc[i] * SXX;
            phi_tmp[64 + i] += 3.0 * 12.0 * xc_pow[i] * zc[i] * SX;
            phi_tmp[64 + i] += 24.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SXXX;
            phi_tmp[96 + i] += 3.0 * 3.0 * xc_pow[i] * yc_pow[i] * SXX;
            phi_tmp[96 + i] += 3.0 * 6.0 * xc[i] * yc_pow[i] * SX;
            phi_tmp[96 + i] += 6.0 * yc_pow[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SXXX;
            phi_tmp[128 + i] += 3.0 * 3.0 * xc_pow[i] * yc[i] * zc[i] * SXX;
            phi_tmp[128 + i] += 3.0 * 6.0 * xc[i] * yc[i] * zc[i] * SX;
            phi_tmp[128 + i] += 6.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SXXX;
            phi_tmp[160 + i] += 3.0 * 3.0 * xc_pow[i] * zc_pow[i] * SXX;
            phi_tmp[160 + i] += 3.0 * 6.0 * xc[i] * zc_pow[i] * SX;
            phi_tmp[160 + i] += 6.0 * zc_pow[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SXXX;
            phi_tmp[192 + i] += 3.0 * 2.0 * xc[i] * yc_pow[32 + i] * SXX;
            phi_tmp[192 + i] += 3.0 * 2.0 * yc_pow[32 + i] * SX;

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SXXX;
            phi_tmp[224 + i] += 3.0 * 2.0 * xc[i] * yc_pow[i] * zc[i] * SXX;
            phi_tmp[224 + i] += 3.0 * 2.0 * yc_pow[i] * zc[i] * SX;

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SXXX;
            phi_tmp[256 + i] += 3.0 * 2.0 * xc[i] * yc[i] * zc_pow[i] * SXX;
            phi_tmp[256 + i] += 3.0 * 2.0 * yc[i] * zc_pow[i] * SX;

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SXXX;
            phi_tmp[288 + i] += 3.0 * 2.0 * xc[i] * zc_pow[32 + i] * SXX;
            phi_tmp[288 + i] += 3.0 * 2.0 * zc_pow[32 + i] * SX;

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SXXX;
            phi_tmp[320 + i] += 3.0 * yc_pow[64 + i] * SXX;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SXXX;
            phi_tmp[352 + i] += 3.0 * yc_pow[32 + i] * zc[i] * SXX;

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SXXX;
            phi_tmp[384 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SXX;

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SXXX;
            phi_tmp[416 + i] += 3.0 * yc[i] * zc_pow[32 + i] * SXX;

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SXXX;
            phi_tmp[448 + i] += 3.0 * zc_pow[64 + i] * SXX;

            phi_tmp[480 + i] = yc_pow[96 + i] * SXXX;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SXXX;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SXXX;

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SXXX;

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SXXX;

            phi_tmp[640 + i] = zc_pow[96 + i] * SXXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
        }

        // Combine XXY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXY = S3[i] * xc[i] * xc[i] * yc[i] + S2[i] * yc[i];

            phi_tmp[i] = xc_pow[96 + i] * SXXY;
            phi_tmp[i] += 2.0 * 5.0 * xc_pow[64 + i] * SXY;
            phi_tmp[i] += 20.0 * xc_pow[32 + i] * SY;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SXXY;
            phi_tmp[32 + i] += 2.0 * 4.0 * xc_pow[32 + i] * yc[i] * SXY;
            phi_tmp[32 + i] += xc_pow[64 + i] * SXX;
            phi_tmp[32 + i] += 12.0 * xc_pow[i] * yc[i] * SY;
            phi_tmp[32 + i] += 2.0 * 4.0 * xc_pow[32 + i] * SX;
            phi_tmp[32 + i] += 12.0 * xc_pow[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SXXY;
            phi_tmp[64 + i] += 2.0 * 4.0 * xc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[64 + i] += 12.0 * xc_pow[i] * zc[i] * SY;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SXXY;
            phi_tmp[96 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[i] * SXY;
            phi_tmp[96 + i] += 2.0 * xc_pow[32 + i] * yc[i] * SXX;
            phi_tmp[96 + i] += 6.0 * xc[i] * yc_pow[i] * SY;
            phi_tmp[96 + i] += 2.0 * 6.0 * xc_pow[i] * yc[i] * SX;
            phi_tmp[96 + i] += 12.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SXXY;
            phi_tmp[128 + i] += 2.0 * 3.0 * xc_pow[i] * yc[i] * zc[i] * SXY;
            phi_tmp[128 + i] += xc_pow[32 + i] * zc[i] * SXX;
            phi_tmp[128 + i] += 6.0 * xc[i] * yc[i] * zc[i] * SY;
            phi_tmp[128 + i] += 2.0 * 3.0 * xc_pow[i] * zc[i] * SX;
            phi_tmp[128 + i] += 6.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SXXY;
            phi_tmp[160 + i] += 2.0 * 3.0 * xc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[160 + i] += 6.0 * xc[i] * zc_pow[i] * SY;

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SXXY;
            phi_tmp[192 + i] += 2.0 * 2.0 * xc[i] * yc_pow[32 + i] * SXY;
            phi_tmp[192 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SXX;
            phi_tmp[192 + i] += 2.0 * yc_pow[32 + i] * SY;
            phi_tmp[192 + i] += 2.0 * 6.0 * xc[i] * yc_pow[i] * SX;
            phi_tmp[192 + i] += 6.0 * yc_pow[i] * S0[i];

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SXXY;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc[i] * yc_pow[i] * zc[i] * SXY;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * yc[i] * zc[i] * SXX;
            phi_tmp[224 + i] += 2.0 * yc_pow[i] * zc[i] * SY;
            phi_tmp[224 + i] += 2.0 * 4.0 * xc[i] * yc[i] * zc[i] * SX;
            phi_tmp[224 + i] += 4.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SXXY;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc_pow[i] * SXY;
            phi_tmp[256 + i] += xc_pow[i] * zc_pow[i] * SXX;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc_pow[i] * SY;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc[i] * zc_pow[i] * SX;
            phi_tmp[256 + i] += 2.0 * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SXXY;
            phi_tmp[288 + i] += 2.0 * 2.0 * xc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[288 + i] += 2.0 * zc_pow[32 + i] * SY;

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SXXY;
            phi_tmp[320 + i] += 2.0 * yc_pow[64 + i] * SXY;
            phi_tmp[320 + i] += 4.0 * xc[i] * yc_pow[32 + i] * SXX;
            phi_tmp[320 + i] += 2.0 * 4.0 * yc_pow[32 + i] * SX;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SXXY;
            phi_tmp[352 + i] += 2.0 * yc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[352 + i] += 3.0 * xc[i] * yc_pow[i] * zc[i] * SXX;
            phi_tmp[352 + i] += 2.0 * 3.0 * yc_pow[i] * zc[i] * SX;

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SXXY;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc[i] * zc_pow[i] * SXX;
            phi_tmp[384 + i] += 2.0 * 2.0 * yc[i] * zc_pow[i] * SX;

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SXXY;
            phi_tmp[416 + i] += 2.0 * yc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[416 + i] += xc[i] * zc_pow[32 + i] * SXX;
            phi_tmp[416 + i] += 2.0 * zc_pow[32 + i] * SX;

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SXXY;
            phi_tmp[448 + i] += 2.0 * zc_pow[64 + i] * SXY;

            phi_tmp[480 + i] = yc_pow[96 + i] * SXXY;
            phi_tmp[480 + i] += 5.0 * yc_pow[64 + i] * SXX;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SXXY;
            phi_tmp[512 + i] += 4.0 * yc_pow[32 + i] * zc[i] * SXX;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SXXY;
            phi_tmp[544 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SXX;

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SXXY;
            phi_tmp[576 + i] += 2.0 * yc[i] * zc_pow[32 + i] * SXX;

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SXXY;
            phi_tmp[608 + i] += zc_pow[64 + i] * SXX;

            phi_tmp[640 + i] = zc_pow[96 + i] * SXXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
        }

        // Combine XXZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXZ = S3[i] * xc[i] * xc[i] * zc[i] + S2[i] * zc[i];

            phi_tmp[i] = xc_pow[96 + i] * SXXZ;
            phi_tmp[i] += 2.0 * 5.0 * xc_pow[64 + i] * SXZ;
            phi_tmp[i] += 20.0 * xc_pow[32 + i] * SZ;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SXXZ;
            phi_tmp[32 + i] += 2.0 * 4.0 * xc_pow[32 + i] * yc[i] * SXZ;
            phi_tmp[32 + i] += 12.0 * xc_pow[i] * yc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SXXZ;
            phi_tmp[64 + i] += 2.0 * 4.0 * xc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[64 + i] += xc_pow[64 + i] * SXX;
            phi_tmp[64 + i] += 12.0 * xc_pow[i] * zc[i] * SZ;
            phi_tmp[64 + i] += 2.0 * 4.0 * xc_pow[32 + i] * SX;
            phi_tmp[64 + i] += 12.0 * xc_pow[i] * S0[i];

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SXXZ;
            phi_tmp[96 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[i] * SXZ;
            phi_tmp[96 + i] += 6.0 * xc[i] * yc_pow[i] * SZ;

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SXXZ;
            phi_tmp[128 + i] += 2.0 * 3.0 * xc_pow[i] * yc[i] * zc[i] * SXZ;
            phi_tmp[128 + i] += xc_pow[32 + i] * yc[i] * SXX;
            phi_tmp[128 + i] += 6.0 * xc[i] * yc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += 2.0 * 3.0 * xc_pow[i] * yc[i] * SX;
            phi_tmp[128 + i] += 6.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SXXZ;
            phi_tmp[160 + i] += 2.0 * 3.0 * xc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[32 + i] * zc[i] * SXX;
            phi_tmp[160 + i] += 6.0 * xc[i] * zc_pow[i] * SZ;
            phi_tmp[160 + i] += 2.0 * 6.0 * xc_pow[i] * zc[i] * SX;
            phi_tmp[160 + i] += 12.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SXXZ;
            phi_tmp[192 + i] += 2.0 * 2.0 * xc[i] * yc_pow[32 + i] * SXZ;
            phi_tmp[192 + i] += 2.0 * yc_pow[32 + i] * SZ;

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SXXZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc[i] * yc_pow[i] * zc[i] * SXZ;
            phi_tmp[224 + i] += xc_pow[i] * yc_pow[i] * SXX;
            phi_tmp[224 + i] += 2.0 * yc_pow[i] * zc[i] * SZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc[i] * yc_pow[i] * SX;
            phi_tmp[224 + i] += 2.0 * yc_pow[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SXXZ;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc_pow[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * yc[i] * zc[i] * SXX;
            phi_tmp[256 + i] += 2.0 * yc[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 2.0 * 4.0 * xc[i] * yc[i] * zc[i] * SX;
            phi_tmp[256 + i] += 4.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SXXZ;
            phi_tmp[288 + i] += 2.0 * 2.0 * xc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SXX;
            phi_tmp[288 + i] += 2.0 * zc_pow[32 + i] * SZ;
            phi_tmp[288 + i] += 2.0 * 6.0 * xc[i] * zc_pow[i] * SX;
            phi_tmp[288 + i] += 6.0 * zc_pow[i] * S0[i];

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SXXZ;
            phi_tmp[320 + i] += 2.0 * yc_pow[64 + i] * SXZ;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SXXZ;
            phi_tmp[352 + i] += 2.0 * yc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[352 + i] += xc[i] * yc_pow[32 + i] * SXX;
            phi_tmp[352 + i] += 2.0 * yc_pow[32 + i] * SX;

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SXXZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * zc[i] * SXX;
            phi_tmp[384 + i] += 2.0 * 2.0 * yc_pow[i] * zc[i] * SX;

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SXXZ;
            phi_tmp[416 + i] += 2.0 * yc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[416 + i] += 3.0 * xc[i] * yc[i] * zc_pow[i] * SXX;
            phi_tmp[416 + i] += 2.0 * 3.0 * yc[i] * zc_pow[i] * SX;

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SXXZ;
            phi_tmp[448 + i] += 2.0 * zc_pow[64 + i] * SXZ;
            phi_tmp[448 + i] += 4.0 * xc[i] * zc_pow[32 + i] * SXX;
            phi_tmp[448 + i] += 2.0 * 4.0 * zc_pow[32 + i] * SX;

            phi_tmp[480 + i] = yc_pow[96 + i] * SXXZ;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SXXZ;
            phi_tmp[512 + i] += yc_pow[64 + i] * SXX;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SXXZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * zc[i] * SXX;

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SXXZ;
            phi_tmp[576 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SXX;

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SXXZ;
            phi_tmp[608 + i] += 4.0 * yc[i] * zc_pow[32 + i] * SXX;

            phi_tmp[640 + i] = zc_pow[96 + i] * SXXZ;
            phi_tmp[640 + i] += 5.0 * zc_pow[64 + i] * SXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
        }

        // Combine XYY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SXYY = S3[i] * xc[i] * yc[i] * yc[i] + S2[i] * xc[i];

            phi_tmp[i] = xc_pow[96 + i] * SXYY;
            phi_tmp[i] += 5.0 * xc_pow[64 + i] * SYY;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SXYY;
            phi_tmp[32 + i] += 2.0 * xc_pow[64 + i] * SXY;
            phi_tmp[32 + i] += 4.0 * xc_pow[32 + i] * yc[i] * SYY;
            phi_tmp[32 + i] += 2.0 * 4.0 * xc_pow[32 + i] * SY;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SXYY;
            phi_tmp[64 + i] += 4.0 * xc_pow[32 + i] * zc[i] * SYY;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SXYY;
            phi_tmp[96 + i] += 2.0 * 2.0 * xc_pow[32 + i] * yc[i] * SXY;
            phi_tmp[96 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SYY;
            phi_tmp[96 + i] += 2.0 * xc_pow[32 + i] * SX;
            phi_tmp[96 + i] += 2.0 * 6.0 * xc_pow[i] * yc[i] * SY;
            phi_tmp[96 + i] += 6.0 * xc_pow[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SXYY;
            phi_tmp[128 + i] += 2.0 * xc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * yc[i] * zc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * 3.0 * xc_pow[i] * zc[i] * SY;

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SXYY;
            phi_tmp[160 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SYY;

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SXYY;
            phi_tmp[192 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[i] * SXY;
            phi_tmp[192 + i] += 2.0 * xc[i] * yc_pow[32 + i] * SYY;
            phi_tmp[192 + i] += 6.0 * xc_pow[i] * yc[i] * SX;
            phi_tmp[192 + i] += 2.0 * 6.0 * xc[i] * yc_pow[i] * SY;
            phi_tmp[192 + i] += 12.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SXYY;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc_pow[i] * yc[i] * zc[i] * SXY;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc_pow[i] * zc[i] * SYY;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * zc[i] * SX;
            phi_tmp[224 + i] += 2.0 * 4.0 * xc[i] * yc[i] * zc[i] * SY;
            phi_tmp[224 + i] += 4.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SXYY;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * zc_pow[i] * SYY;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc[i] * zc_pow[i] * SY;

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SXYY;
            phi_tmp[288 + i] += 2.0 * xc[i] * zc_pow[32 + i] * SYY;

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SXYY;
            phi_tmp[320 + i] += 2.0 * 4.0 * xc[i] * yc_pow[32 + i] * SXY;
            phi_tmp[320 + i] += yc_pow[64 + i] * SYY;
            phi_tmp[320 + i] += 12.0 * xc[i] * yc_pow[i] * SX;
            phi_tmp[320 + i] += 2.0 * 4.0 * yc_pow[32 + i] * SY;
            phi_tmp[320 + i] += 12.0 * yc_pow[i] * S0[i];

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SXYY;
            phi_tmp[352 + i] += 2.0 * 3.0 * xc[i] * yc_pow[i] * zc[i] * SXY;
            phi_tmp[352 + i] += yc_pow[32 + i] * zc[i] * SYY;
            phi_tmp[352 + i] += 6.0 * xc[i] * yc[i] * zc[i] * SX;
            phi_tmp[352 + i] += 2.0 * 3.0 * yc_pow[i] * zc[i] * SY;
            phi_tmp[352 + i] += 6.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SXYY;
            phi_tmp[384 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc_pow[i] * SXY;
            phi_tmp[384 + i] += yc_pow[i] * zc_pow[i] * SYY;
            phi_tmp[384 + i] += 2.0 * xc[i] * zc_pow[i] * SX;
            phi_tmp[384 + i] += 2.0 * 2.0 * yc[i] * zc_pow[i] * SY;
            phi_tmp[384 + i] += 2.0 * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SXYY;
            phi_tmp[416 + i] += 2.0 * xc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[416 + i] += yc[i] * zc_pow[32 + i] * SYY;
            phi_tmp[416 + i] += 2.0 * zc_pow[32 + i] * SY;

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SXYY;
            phi_tmp[448 + i] += zc_pow[64 + i] * SYY;

            phi_tmp[480 + i] = yc_pow[96 + i] * SXYY;
            phi_tmp[480 + i] += 2.0 * 5.0 * yc_pow[64 + i] * SXY;
            phi_tmp[480 + i] += 20.0 * yc_pow[32 + i] * SX;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SXYY;
            phi_tmp[512 + i] += 2.0 * 4.0 * yc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[512 + i] += 12.0 * yc_pow[i] * zc[i] * SX;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SXYY;
            phi_tmp[544 + i] += 2.0 * 3.0 * yc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[544 + i] += 6.0 * yc[i] * zc_pow[i] * SX;

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SXYY;
            phi_tmp[576 + i] += 2.0 * 2.0 * yc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[576 + i] += 2.0 * zc_pow[32 + i] * SX;

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SXYY;
            phi_tmp[608 + i] += 2.0 * zc_pow[64 + i] * SXY;

            phi_tmp[640 + i] = zc_pow[96 + i] * SXYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
        }

        // Combine XYZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SXYZ = S3[i] * xc[i] * yc[i] * zc[i];

            phi_tmp[i] = xc_pow[96 + i] * SXYZ;
            phi_tmp[i] += 5.0 * xc_pow[64 + i] * SYZ;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SXYZ;
            phi_tmp[32 + i] += 4.0 * xc_pow[32 + i] * yc[i] * SYZ;
            phi_tmp[32 + i] += xc_pow[64 + i] * SXZ;
            phi_tmp[32 + i] += 4.0 * xc_pow[32 + i] * SZ;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SXYZ;
            phi_tmp[64 + i] += 4.0 * xc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[64 + i] += xc_pow[64 + i] * SXY;
            phi_tmp[64 + i] += 4.0 * xc_pow[32 + i] * SY;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SXYZ;
            phi_tmp[96 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * xc_pow[32 + i] * yc[i] * SXZ;
            phi_tmp[96 + i] += 6.0 * xc_pow[i] * yc[i] * SZ;

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SXYZ;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * yc[i] * zc[i] * SYZ;
            phi_tmp[128 + i] += xc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[128 + i] += xc_pow[32 + i] * yc[i] * SXY;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * zc[i] * SZ;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * yc[i] * SY;
            phi_tmp[128 + i] += xc_pow[32 + i] * SX;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SXYZ;
            phi_tmp[160 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[160 + i] += 6.0 * xc_pow[i] * zc[i] * SY;

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SXYZ;
            phi_tmp[192 + i] += 2.0 * xc[i] * yc_pow[32 + i] * SYZ;
            phi_tmp[192 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SXZ;
            phi_tmp[192 + i] += 6.0 * xc[i] * yc_pow[i] * SZ;

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SXYZ;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc_pow[i] * zc[i] * SYZ;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * yc[i] * zc[i] * SXZ;
            phi_tmp[224 + i] += xc_pow[i] * yc_pow[i] * SXY;
            phi_tmp[224 + i] += 4.0 * xc[i] * yc[i] * zc[i] * SZ;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc_pow[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * yc[i] * SX;
            phi_tmp[224 + i] += 8.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SXYZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * zc_pow[i] * SYZ;
            phi_tmp[256 + i] += xc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * yc[i] * zc[i] * SXY;
            phi_tmp[256 + i] += 2.0 * xc[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 4.0 * xc[i] * yc[i] * zc[i] * SY;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * zc[i] * SX;
            phi_tmp[256 + i] += 2.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SXYZ;
            phi_tmp[288 + i] += 2.0 * xc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[288 + i] += 6.0 * xc[i] * zc_pow[i] * SY;

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SXYZ;
            phi_tmp[320 + i] += yc_pow[64 + i] * SYZ;
            phi_tmp[320 + i] += 4.0 * xc[i] * yc_pow[32 + i] * SXZ;
            phi_tmp[320 + i] += 4.0 * yc_pow[32 + i] * SZ;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SXYZ;
            phi_tmp[352 + i] += yc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[352 + i] += 3.0 * xc[i] * yc_pow[i] * zc[i] * SXZ;
            phi_tmp[352 + i] += xc[i] * yc_pow[32 + i] * SXY;
            phi_tmp[352 + i] += 3.0 * yc_pow[i] * zc[i] * SZ;
            phi_tmp[352 + i] += yc_pow[32 + i] * SY;
            phi_tmp[352 + i] += 3.0 * xc[i] * yc_pow[i] * SX;
            phi_tmp[352 + i] += 9.0 * yc_pow[i] * S0[i];

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SXYZ;
            phi_tmp[384 + i] += yc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc[i] * zc_pow[i] * SXZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * zc[i] * SXY;
            phi_tmp[384 + i] += 2.0 * yc[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc[i] * SY;
            phi_tmp[384 + i] += 4.0 * xc[i] * yc[i] * zc[i] * SX;
            phi_tmp[384 + i] += 4.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SXYZ;
            phi_tmp[416 + i] += yc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[416 + i] += xc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[416 + i] += 3.0 * xc[i] * yc[i] * zc_pow[i] * SXY;
            phi_tmp[416 + i] += zc_pow[32 + i] * SZ;
            phi_tmp[416 + i] += 3.0 * yc[i] * zc_pow[i] * SY;
            phi_tmp[416 + i] += 3.0 * xc[i] * zc_pow[i] * SX;
            phi_tmp[416 + i] += zc_pow[i] * S0[i];

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SXYZ;
            phi_tmp[448 + i] += zc_pow[64 + i] * SYZ;
            phi_tmp[448 + i] += 4.0 * xc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[448 + i] += 4.0 * zc_pow[32 + i] * SY;

            phi_tmp[480 + i] = yc_pow[96 + i] * SXYZ;
            phi_tmp[480 + i] += 5.0 * yc_pow[64 + i] * SXZ;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SXYZ;
            phi_tmp[512 + i] += 4.0 * yc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[512 + i] += yc_pow[64 + i] * SXY;
            phi_tmp[512 + i] += 4.0 * yc_pow[32 + i] * SX;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SXYZ;
            phi_tmp[544 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[544 + i] += 6.0 * yc_pow[i] * zc[i] * SX;

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SXYZ;
            phi_tmp[576 + i] += 2.0 * yc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[576 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[576 + i] += 6.0 * yc[i] * zc_pow[i] * SX;

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SXYZ;
            phi_tmp[608 + i] += zc_pow[64 + i] * SXZ;
            phi_tmp[608 + i] += 4.0 * yc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[608 + i] += 4.0 * zc_pow[32 + i] * SX;

            phi_tmp[640 + i] = zc_pow[96 + i] * SXYZ;
            phi_tmp[640 + i] += 5.0 * zc_pow[64 + i] * SXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
        }

        // Combine XZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SXZZ = S3[i] * xc[i] * zc[i] * zc[i] + S2[i] * xc[i];

            phi_tmp[i] = xc_pow[96 + i] * SXZZ;
            phi_tmp[i] += 5.0 * xc_pow[64 + i] * SZZ;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SXZZ;
            phi_tmp[32 + i] += 4.0 * xc_pow[32 + i] * yc[i] * SZZ;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SXZZ;
            phi_tmp[64 + i] += 2.0 * xc_pow[64 + i] * SXZ;
            phi_tmp[64 + i] += 4.0 * xc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[64 + i] += 2.0 * 4.0 * xc_pow[32 + i] * SZ;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SXZZ;
            phi_tmp[96 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SZZ;

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SXZZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[32 + i] * yc[i] * SXZ;
            phi_tmp[128 + i] += 3.0 * xc_pow[i] * yc[i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * 3.0 * xc_pow[i] * yc[i] * SZ;

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SXZZ;
            phi_tmp[160 + i] += 2.0 * 2.0 * xc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[160 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[32 + i] * SX;
            phi_tmp[160 + i] += 2.0 * 6.0 * xc_pow[i] * zc[i] * SZ;
            phi_tmp[160 + i] += 6.0 * xc_pow[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SXZZ;
            phi_tmp[192 + i] += 2.0 * xc[i] * yc_pow[32 + i] * SZZ;

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SXZZ;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * yc_pow[i] * SXZ;
            phi_tmp[224 + i] += 2.0 * xc[i] * yc_pow[i] * zc[i] * SZZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc[i] * yc_pow[i] * SZ;

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SXZZ;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc_pow[i] * yc[i] * zc[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * xc[i] * yc[i] * zc_pow[i] * SZZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * yc[i] * SX;
            phi_tmp[256 + i] += 2.0 * 4.0 * xc[i] * yc[i] * zc[i] * SZ;
            phi_tmp[256 + i] += 4.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SXZZ;
            phi_tmp[288 + i] += 2.0 * 3.0 * xc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[288 + i] += 2.0 * xc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[288 + i] += 6.0 * xc_pow[i] * zc[i] * SX;
            phi_tmp[288 + i] += 2.0 * 6.0 * xc[i] * zc_pow[i] * SZ;
            phi_tmp[288 + i] += 12.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SXZZ;
            phi_tmp[320 + i] += yc_pow[64 + i] * SZZ;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SXZZ;
            phi_tmp[352 + i] += 2.0 * xc[i] * yc_pow[32 + i] * SXZ;
            phi_tmp[352 + i] += yc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[352 + i] += 2.0 * yc_pow[32 + i] * SZ;

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SXZZ;
            phi_tmp[384 + i] += 2.0 * 2.0 * xc[i] * yc_pow[i] * zc[i] * SXZ;
            phi_tmp[384 + i] += yc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * SX;
            phi_tmp[384 + i] += 2.0 * 2.0 * yc_pow[i] * zc[i] * SZ;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SXZZ;
            phi_tmp[416 + i] += 2.0 * 3.0 * xc[i] * yc[i] * zc_pow[i] * SXZ;
            phi_tmp[416 + i] += yc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[416 + i] += 6.0 * xc[i] * yc[i] * zc[i] * SX;
            phi_tmp[416 + i] += 2.0 * 3.0 * yc[i] * zc_pow[i] * SZ;
            phi_tmp[416 + i] += 6.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SXZZ;
            phi_tmp[448 + i] += 2.0 * 4.0 * xc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[448 + i] += zc_pow[64 + i] * SZZ;
            phi_tmp[448 + i] += 12.0 * xc[i] * zc_pow[i] * SX;
            phi_tmp[448 + i] += 2.0 * 4.0 * zc_pow[32 + i] * SZ;
            phi_tmp[448 + i] += 12.0 * zc_pow[i] * S0[i];

            phi_tmp[480 + i] = yc_pow[96 + i] * SXZZ;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SXZZ;
            phi_tmp[512 + i] += 2.0 * yc_pow[64 + i] * SXZ;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SXZZ;
            phi_tmp[544 + i] += 2.0 * 2.0 * yc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * SX;

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SXZZ;
            phi_tmp[576 + i] += 2.0 * 3.0 * yc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[576 + i] += 6.0 * yc_pow[i] * zc[i] * SX;

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SXZZ;
            phi_tmp[608 + i] += 2.0 * 4.0 * yc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[608 + i] += 12.0 * yc[i] * zc_pow[i] * SX;

            phi_tmp[640 + i] = zc_pow[96 + i] * SXZZ;
            phi_tmp[640 + i] += 2.0 * 5.0 * zc_pow[64 + i] * SXZ;
            phi_tmp[640 + i] += 20.0 * zc_pow[32 + i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
        }

        // Combine YYY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SYYY = S3[i] * yc[i] * yc[i] * yc[i] + 3 * yc[i] * S2[i];

            phi_tmp[i] = xc_pow[96 + i] * SYYY;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SYYY;
            phi_tmp[32 + i] += 3.0 * xc_pow[64 + i] * SYY;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SYYY;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SYYY;
            phi_tmp[96 + i] += 3.0 * 2.0 * xc_pow[32 + i] * yc[i] * SYY;
            phi_tmp[96 + i] += 3.0 * 2.0 * xc_pow[32 + i] * SY;

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SYYY;
            phi_tmp[128 + i] += 3.0 * xc_pow[32 + i] * zc[i] * SYY;

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SYYY;

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SYYY;
            phi_tmp[192 + i] += 3.0 * 3.0 * xc_pow[i] * yc_pow[i] * SYY;
            phi_tmp[192 + i] += 3.0 * 6.0 * xc_pow[i] * yc[i] * SY;
            phi_tmp[192 + i] += 6.0 * xc_pow[i] * S0[i];

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SYYY;
            phi_tmp[224 + i] += 3.0 * 2.0 * xc_pow[i] * yc[i] * zc[i] * SYY;
            phi_tmp[224 + i] += 3.0 * 2.0 * xc_pow[i] * zc[i] * SY;

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SYYY;
            phi_tmp[256 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SYY;

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SYYY;

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SYYY;
            phi_tmp[320 + i] += 3.0 * 4.0 * xc[i] * yc_pow[32 + i] * SYY;
            phi_tmp[320 + i] += 3.0 * 12.0 * xc[i] * yc_pow[i] * SY;
            phi_tmp[320 + i] += 24.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SYYY;
            phi_tmp[352 + i] += 3.0 * 3.0 * xc[i] * yc_pow[i] * zc[i] * SYY;
            phi_tmp[352 + i] += 3.0 * 6.0 * xc[i] * yc[i] * zc[i] * SY;
            phi_tmp[352 + i] += 6.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SYYY;
            phi_tmp[384 + i] += 3.0 * 2.0 * xc[i] * yc[i] * zc_pow[i] * SYY;
            phi_tmp[384 + i] += 3.0 * 2.0 * xc[i] * zc_pow[i] * SY;

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SYYY;
            phi_tmp[416 + i] += 3.0 * xc[i] * zc_pow[32 + i] * SYY;

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SYYY;

            phi_tmp[480 + i] = yc_pow[96 + i] * SYYY;
            phi_tmp[480 + i] += 3.0 * 5.0 * yc_pow[64 + i] * SYY;
            phi_tmp[480 + i] += 3.0 * 20.0 * yc_pow[32 + i] * SY;
            phi_tmp[480 + i] += 60.0 * yc_pow[i] * S0[i];

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SYYY;
            phi_tmp[512 + i] += 3.0 * 4.0 * yc_pow[32 + i] * zc[i] * SYY;
            phi_tmp[512 + i] += 3.0 * 12.0 * yc_pow[i] * zc[i] * SY;
            phi_tmp[512 + i] += 24.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SYYY;
            phi_tmp[544 + i] += 3.0 * 3.0 * yc_pow[i] * zc_pow[i] * SYY;
            phi_tmp[544 + i] += 3.0 * 6.0 * yc[i] * zc_pow[i] * SY;
            phi_tmp[544 + i] += 6.0 * zc_pow[i] * S0[i];

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SYYY;
            phi_tmp[576 + i] += 3.0 * 2.0 * yc[i] * zc_pow[32 + i] * SYY;
            phi_tmp[576 + i] += 3.0 * 2.0 * zc_pow[32 + i] * SY;

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SYYY;
            phi_tmp[608 + i] += 3.0 * zc_pow[64 + i] * SYY;

            phi_tmp[640 + i] = zc_pow[96 + i] * SYYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
        }

        // Combine YYZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SYYZ = S3[i] * yc[i] * yc[i] * zc[i] + S2[i] * zc[i];

            phi_tmp[i] = xc_pow[96 + i] * SYYZ;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SYYZ;
            phi_tmp[32 + i] += 2.0 * xc_pow[64 + i] * SYZ;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SYYZ;
            phi_tmp[64 + i] += xc_pow[64 + i] * SYY;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SYYZ;
            phi_tmp[96 + i] += 2.0 * 2.0 * xc_pow[32 + i] * yc[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * xc_pow[32 + i] * SZ;

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SYYZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[128 + i] += xc_pow[32 + i] * yc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * xc_pow[32 + i] * SY;

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SYYZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[32 + i] * zc[i] * SYY;

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SYYZ;
            phi_tmp[192 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[i] * SYZ;
            phi_tmp[192 + i] += 6.0 * xc_pow[i] * yc[i] * SZ;

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SYYZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc_pow[i] * yc[i] * zc[i] * SYZ;
            phi_tmp[224 + i] += xc_pow[i] * yc_pow[i] * SYY;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * zc[i] * SZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc_pow[i] * yc[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SYYZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * yc[i] * zc[i] * SYY;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc_pow[i] * zc[i] * SY;

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SYYZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SYY;

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SYYZ;
            phi_tmp[320 + i] += 2.0 * 4.0 * xc[i] * yc_pow[32 + i] * SYZ;
            phi_tmp[320 + i] += 12.0 * xc[i] * yc_pow[i] * SZ;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SYYZ;
            phi_tmp[352 + i] += 2.0 * 3.0 * xc[i] * yc_pow[i] * zc[i] * SYZ;
            phi_tmp[352 + i] += xc[i] * yc_pow[32 + i] * SYY;
            phi_tmp[352 + i] += 6.0 * xc[i] * yc[i] * zc[i] * SZ;
            phi_tmp[352 + i] += 2.0 * 3.0 * xc[i] * yc_pow[i] * SY;
            phi_tmp[352 + i] += 6.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SYYZ;
            phi_tmp[384 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc_pow[i] * SYZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * zc[i] * SYY;
            phi_tmp[384 + i] += 2.0 * xc[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 2.0 * 4.0 * xc[i] * yc[i] * zc[i] * SY;
            phi_tmp[384 + i] += 4.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SYYZ;
            phi_tmp[416 + i] += 2.0 * xc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[416 + i] += 3.0 * xc[i] * yc[i] * zc_pow[i] * SYY;
            phi_tmp[416 + i] += 2.0 * 3.0 * xc[i] * zc_pow[i] * SY;

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SYYZ;
            phi_tmp[448 + i] += 4.0 * xc[i] * zc_pow[32 + i] * SYY;

            phi_tmp[480 + i] = yc_pow[96 + i] * SYYZ;
            phi_tmp[480 + i] += 2.0 * 5.0 * yc_pow[64 + i] * SYZ;
            phi_tmp[480 + i] += 20.0 * yc_pow[32 + i] * SZ;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SYYZ;
            phi_tmp[512 + i] += 2.0 * 4.0 * yc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[512 + i] += yc_pow[64 + i] * SYY;
            phi_tmp[512 + i] += 12.0 * yc_pow[i] * zc[i] * SZ;
            phi_tmp[512 + i] += 2.0 * 4.0 * yc_pow[32 + i] * SY;
            phi_tmp[512 + i] += 12.0 * yc_pow[i] * S0[i];

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SYYZ;
            phi_tmp[544 + i] += 2.0 * 3.0 * yc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * zc[i] * SYY;
            phi_tmp[544 + i] += 6.0 * yc[i] * zc_pow[i] * SZ;
            phi_tmp[544 + i] += 2.0 * 6.0 * yc_pow[i] * zc[i] * SY;
            phi_tmp[544 + i] += 12.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SYYZ;
            phi_tmp[576 + i] += 2.0 * 2.0 * yc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[576 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SYY;
            phi_tmp[576 + i] += 2.0 * zc_pow[32 + i] * SZ;
            phi_tmp[576 + i] += 2.0 * 6.0 * yc[i] * zc_pow[i] * SY;
            phi_tmp[576 + i] += 6.0 * zc_pow[i] * S0[i];

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SYYZ;
            phi_tmp[608 + i] += 2.0 * zc_pow[64 + i] * SYZ;
            phi_tmp[608 + i] += 4.0 * yc[i] * zc_pow[32 + i] * SYY;
            phi_tmp[608 + i] += 2.0 * 4.0 * zc_pow[32 + i] * SY;

            phi_tmp[640 + i] = zc_pow[96 + i] * SYYZ;
            phi_tmp[640 + i] += 5.0 * zc_pow[64 + i] * SYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
        }

        // Combine YZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SYZZ = S3[i] * yc[i] * zc[i] * zc[i] + S2[i] * yc[i];

            phi_tmp[i] = xc_pow[96 + i] * SYZZ;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SYZZ;
            phi_tmp[32 + i] += xc_pow[64 + i] * SZZ;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SYZZ;
            phi_tmp[64 + i] += 2.0 * xc_pow[64 + i] * SYZ;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SYZZ;
            phi_tmp[96 + i] += 2.0 * xc_pow[32 + i] * yc[i] * SZZ;

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SYZZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[32 + i] * yc[i] * SYZ;
            phi_tmp[128 + i] += xc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[32 + i] * SZ;

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SYZZ;
            phi_tmp[160 + i] += 2.0 * 2.0 * xc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[32 + i] * SY;

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SYZZ;
            phi_tmp[192 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SZZ;

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SYZZ;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * yc_pow[i] * SYZ;
            phi_tmp[224 + i] += 2.0 * xc_pow[i] * yc[i] * zc[i] * SZZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc_pow[i] * yc[i] * SZ;

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SYZZ;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc_pow[i] * yc[i] * zc[i] * SYZ;
            phi_tmp[256 + i] += xc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * yc[i] * SY;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc_pow[i] * zc[i] * SZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SYZZ;
            phi_tmp[288 + i] += 2.0 * 3.0 * xc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[288 + i] += 6.0 * xc_pow[i] * zc[i] * SY;

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SYZZ;
            phi_tmp[320 + i] += 4.0 * xc[i] * yc_pow[32 + i] * SZZ;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SYZZ;
            phi_tmp[352 + i] += 2.0 * xc[i] * yc_pow[32 + i] * SYZ;
            phi_tmp[352 + i] += 3.0 * xc[i] * yc_pow[i] * zc[i] * SZZ;
            phi_tmp[352 + i] += 2.0 * 3.0 * xc[i] * yc_pow[i] * SZ;

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SYZZ;
            phi_tmp[384 + i] += 2.0 * 2.0 * xc[i] * yc_pow[i] * zc[i] * SYZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc[i] * zc_pow[i] * SZZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * SY;
            phi_tmp[384 + i] += 2.0 * 4.0 * xc[i] * yc[i] * zc[i] * SZ;
            phi_tmp[384 + i] += 4.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SYZZ;
            phi_tmp[416 + i] += 2.0 * 3.0 * xc[i] * yc[i] * zc_pow[i] * SYZ;
            phi_tmp[416 + i] += xc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[416 + i] += 6.0 * xc[i] * yc[i] * zc[i] * SY;
            phi_tmp[416 + i] += 2.0 * 3.0 * xc[i] * zc_pow[i] * SZ;
            phi_tmp[416 + i] += 6.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SYZZ;
            phi_tmp[448 + i] += 2.0 * 4.0 * xc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[448 + i] += 12.0 * xc[i] * zc_pow[i] * SY;

            phi_tmp[480 + i] = yc_pow[96 + i] * SYZZ;
            phi_tmp[480 + i] += 5.0 * yc_pow[64 + i] * SZZ;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SYZZ;
            phi_tmp[512 + i] += 2.0 * yc_pow[64 + i] * SYZ;
            phi_tmp[512 + i] += 4.0 * yc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[512 + i] += 2.0 * 4.0 * yc_pow[32 + i] * SZ;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SYZZ;
            phi_tmp[544 + i] += 2.0 * 2.0 * yc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[544 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * SY;
            phi_tmp[544 + i] += 2.0 * 6.0 * yc_pow[i] * zc[i] * SZ;
            phi_tmp[544 + i] += 6.0 * yc_pow[i] * S0[i];

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SYZZ;
            phi_tmp[576 + i] += 2.0 * 3.0 * yc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[576 + i] += 2.0 * yc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[576 + i] += 6.0 * yc_pow[i] * zc[i] * SY;
            phi_tmp[576 + i] += 2.0 * 6.0 * yc[i] * zc_pow[i] * SZ;
            phi_tmp[576 + i] += 12.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SYZZ;
            phi_tmp[608 + i] += 2.0 * 4.0 * yc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[608 + i] += zc_pow[64 + i] * SZZ;
            phi_tmp[608 + i] += 12.0 * yc[i] * zc_pow[i] * SY;
            phi_tmp[608 + i] += 2.0 * 4.0 * zc_pow[32 + i] * SZ;
            phi_tmp[608 + i] += 12.0 * zc_pow[i] * S0[i];

            phi_tmp[640 + i] = zc_pow[96 + i] * SYZZ;
            phi_tmp[640 + i] += 2.0 * 5.0 * zc_pow[64 + i] * SYZ;
            phi_tmp[640 + i] += 20.0 * zc_pow[32 + i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
        }

        // Combine ZZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SZZZ = S3[i] * zc[i] * zc[i] * zc[i] + 3 * zc[i] * S2[i];

            phi_tmp[i] = xc_pow[96 + i] * SZZZ;

            phi_tmp[32 + i] = xc_pow[64 + i] * yc[i] * SZZZ;

            phi_tmp[64 + i] = xc_pow[64 + i] * zc[i] * SZZZ;
            phi_tmp[64 + i] += 3.0 * xc_pow[64 + i] * SZZ;

            phi_tmp[96 + i] = xc_pow[32 + i] * yc_pow[i] * SZZZ;

            phi_tmp[128 + i] = xc_pow[32 + i] * yc[i] * zc[i] * SZZZ;
            phi_tmp[128 + i] += 3.0 * xc_pow[32 + i] * yc[i] * SZZ;

            phi_tmp[160 + i] = xc_pow[32 + i] * zc_pow[i] * SZZZ;
            phi_tmp[160 + i] += 3.0 * 2.0 * xc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[160 + i] += 3.0 * 2.0 * xc_pow[32 + i] * SZ;

            phi_tmp[192 + i] = xc_pow[i] * yc_pow[32 + i] * SZZZ;

            phi_tmp[224 + i] = xc_pow[i] * yc_pow[i] * zc[i] * SZZZ;
            phi_tmp[224 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SZZ;

            phi_tmp[256 + i] = xc_pow[i] * yc[i] * zc_pow[i] * SZZZ;
            phi_tmp[256 + i] += 3.0 * 2.0 * xc_pow[i] * yc[i] * zc[i] * SZZ;
            phi_tmp[256 + i] += 3.0 * 2.0 * xc_pow[i] * yc[i] * SZ;

            phi_tmp[288 + i] = xc_pow[i] * zc_pow[32 + i] * SZZZ;
            phi_tmp[288 + i] += 3.0 * 3.0 * xc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[288 + i] += 3.0 * 6.0 * xc_pow[i] * zc[i] * SZ;
            phi_tmp[288 + i] += 6.0 * xc_pow[i] * S0[i];

            phi_tmp[320 + i] = xc[i] * yc_pow[64 + i] * SZZZ;

            phi_tmp[352 + i] = xc[i] * yc_pow[32 + i] * zc[i] * SZZZ;
            phi_tmp[352 + i] += 3.0 * xc[i] * yc_pow[32 + i] * SZZ;

            phi_tmp[384 + i] = xc[i] * yc_pow[i] * zc_pow[i] * SZZZ;
            phi_tmp[384 + i] += 3.0 * 2.0 * xc[i] * yc_pow[i] * zc[i] * SZZ;
            phi_tmp[384 + i] += 3.0 * 2.0 * xc[i] * yc_pow[i] * SZ;

            phi_tmp[416 + i] = xc[i] * yc[i] * zc_pow[32 + i] * SZZZ;
            phi_tmp[416 + i] += 3.0 * 3.0 * xc[i] * yc[i] * zc_pow[i] * SZZ;
            phi_tmp[416 + i] += 3.0 * 6.0 * xc[i] * yc[i] * zc[i] * SZ;
            phi_tmp[416 + i] += 6.0 * xc[i] * yc[i] * S0[i];

            phi_tmp[448 + i] = xc[i] * zc_pow[64 + i] * SZZZ;
            phi_tmp[448 + i] += 3.0 * 4.0 * xc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[448 + i] += 3.0 * 12.0 * xc[i] * zc_pow[i] * SZ;
            phi_tmp[448 + i] += 24.0 * xc[i] * zc[i] * S0[i];

            phi_tmp[480 + i] = yc_pow[96 + i] * SZZZ;

            phi_tmp[512 + i] = yc_pow[64 + i] * zc[i] * SZZZ;
            phi_tmp[512 + i] += 3.0 * yc_pow[64 + i] * SZZ;

            phi_tmp[544 + i] = yc_pow[32 + i] * zc_pow[i] * SZZZ;
            phi_tmp[544 + i] += 3.0 * 2.0 * yc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[544 + i] += 3.0 * 2.0 * yc_pow[32 + i] * SZ;

            phi_tmp[576 + i] = yc_pow[i] * zc_pow[32 + i] * SZZZ;
            phi_tmp[576 + i] += 3.0 * 3.0 * yc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[576 + i] += 3.0 * 6.0 * yc_pow[i] * zc[i] * SZ;
            phi_tmp[576 + i] += 6.0 * yc_pow[i] * S0[i];

            phi_tmp[608 + i] = yc[i] * zc_pow[64 + i] * SZZZ;
            phi_tmp[608 + i] += 3.0 * 4.0 * yc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[608 + i] += 3.0 * 12.0 * yc[i] * zc_pow[i] * SZ;
            phi_tmp[608 + i] += 24.0 * yc[i] * zc[i] * S0[i];

            phi_tmp[640 + i] = zc_pow[96 + i] * SZZZ;
            phi_tmp[640 + i] += 3.0 * 5.0 * zc_pow[64 + i] * SZZ;
            phi_tmp[640 + i] += 3.0 * 20.0 * zc_pow[32 + i] * SZ;
            phi_tmp[640 + i] += 60.0 * zc_pow[i] * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
        }

    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);
    ALIGNED_FREE(expn2);

    // Free Power temporaries
    ALIGNED_FREE(xc_pow);
    ALIGNED_FREE(yc_pow);
    ALIGNED_FREE(zc_pow);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);

}

void gg_collocation_L6_deriv3(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out, double* PRAGMA_RESTRICT phi_xx_out, double* PRAGMA_RESTRICT phi_xy_out, double* PRAGMA_RESTRICT phi_xz_out, double* PRAGMA_RESTRICT phi_yy_out, double* PRAGMA_RESTRICT phi_yz_out, double* PRAGMA_RESTRICT phi_zz_out, double* PRAGMA_RESTRICT phi_xxx_out, double* PRAGMA_RESTRICT phi_xxy_out, double* PRAGMA_RESTRICT phi_xxz_out, double* PRAGMA_RESTRICT phi_xyy_out, double* PRAGMA_RESTRICT phi_xyz_out, double* PRAGMA_RESTRICT phi_xzz_out, double* PRAGMA_RESTRICT phi_yyy_out, double* PRAGMA_RESTRICT phi_yyz_out, double* PRAGMA_RESTRICT phi_yzz_out, double* PRAGMA_RESTRICT phi_zzz_out) {

    // Sizing
    unsigned long nblocks = npoints / 32;
    nblocks += (npoints % 32) ? 1 : 0;
    const unsigned long ncart = 28;
    const unsigned long nspherical = 13;
    unsigned long nout;

    if ((order == GG_SPHERICAL_CCA) || (order == GG_SPHERICAL_GAUSSIAN)) {
        nout = nspherical;
        } else {
        nout = ncart;
    }

    // Allocate S temporaries, single block to stay on cache
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 288 * sizeof(double));
    double* PRAGMA_RESTRICT xc = cache_data + 0;
    ASSUME_ALIGNED(xc, 64);
    double* PRAGMA_RESTRICT yc = cache_data + 32;
    ASSUME_ALIGNED(yc, 64);
    double* PRAGMA_RESTRICT zc = cache_data + 64;
    ASSUME_ALIGNED(zc, 64);
    double* PRAGMA_RESTRICT R2 = cache_data + 96;
    ASSUME_ALIGNED(R2, 64);
    double* PRAGMA_RESTRICT S0 = cache_data + 128;
    ASSUME_ALIGNED(S0, 64);
    double* PRAGMA_RESTRICT tmp1 = cache_data + 160;
    ASSUME_ALIGNED(tmp1, 64);
    double* PRAGMA_RESTRICT S1 = cache_data + 192;
    ASSUME_ALIGNED(S1, 64);
    double* PRAGMA_RESTRICT S2 = cache_data + 224;
    ASSUME_ALIGNED(S2, 64);
    double* PRAGMA_RESTRICT S3 = cache_data + 256;
    ASSUME_ALIGNED(S3, 64);

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));
    double* PRAGMA_RESTRICT expn2 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));

    // Allocate power temporaries
    double* PRAGMA_RESTRICT xc_pow = (double*)ALIGNED_MALLOC(64, 160 * sizeof(double));
    ASSUME_ALIGNED(xc_pow, 64);
    double* PRAGMA_RESTRICT yc_pow = (double*)ALIGNED_MALLOC(64, 160 * sizeof(double));
    ASSUME_ALIGNED(yc_pow, 64);
    double* PRAGMA_RESTRICT zc_pow = (double*)ALIGNED_MALLOC(64, 160 * sizeof(double));
    ASSUME_ALIGNED(zc_pow, 64);

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, 896 * sizeof(double));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;
    double AXX, AXY, AXZ, AYY, AYZ, AZZ;
    double AXXX, XXY, XXZ, XYY, XYZ, XZZ, YYY, YYZ, YZZ, ZZZ;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
        expn2[i] = -2.0 * exponents[i];
    }

    // Start outer block loop
    for (unsigned long block = 0; block < nblocks; block++) {


        // Copy data into inner temps
        const unsigned long start = block * 32;
        const unsigned long remain = ((start + 32) > npoints) ? (npoints - start) : 32;

        // Handle non-AM dependant temps
        if (xyz_stride == 1) {
            const double* PRAGMA_RESTRICT x = xyz + start;
            const double* PRAGMA_RESTRICT y = xyz + npoints + start;
            const double* PRAGMA_RESTRICT z = xyz + 2 * npoints + start;
            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = x[i] - center_x;
                yc[i] = y[i] - center_y;
                zc[i] = z[i] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
            } else {
            unsigned int start_shift = start * xyz_stride;

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                xc[i] = xyz[start_shift + i * xyz_stride] - center_x;
                yc[i] = xyz[start_shift + i * xyz_stride + 1] - center_y;
                zc[i] = xyz[start_shift + i * xyz_stride + 2] - center_z;

                // Distance
                R2[i] = xc[i] * xc[i];
                R2[i] += yc[i] * yc[i];
                R2[i] += zc[i] * zc[i];

                // Zero out S tmps
                S0[i] = 0.0;
                S1[i] = 0.0;
                S2[i] = 0.0;
                S3[i] = 0.0;
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];
            const double alpha_n2 = expn2[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
                const double T2 = alpha_n2 * T1;
                S1[i] += T2;
                const double T3 = alpha_n2 * T2;
                S2[i] += T3;
                const double T4 = alpha_n2 * T3;
                S3[i] += T4;
            }

        }

        // Build powers
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {

            // Cartesian derivs
            xc_pow[i] = xc[i] * xc[i];
            yc_pow[i] = yc[i] * yc[i];
            zc_pow[i] = zc[i] * zc[i];
            xc_pow[32 + i] = xc_pow[i] * xc[i];
            yc_pow[32 + i] = yc_pow[i] * yc[i];
            zc_pow[32 + i] = zc_pow[i] * zc[i];
            xc_pow[64 + i] = xc_pow[32 + i] * xc[i];
            yc_pow[64 + i] = yc_pow[32 + i] * yc[i];
            zc_pow[64 + i] = zc_pow[32 + i] * zc[i];
            xc_pow[96 + i] = xc_pow[64 + i] * xc[i];
            yc_pow[96 + i] = yc_pow[64 + i] * yc[i];
            zc_pow[96 + i] = zc_pow[64 + i] * zc[i];
            xc_pow[128 + i] = xc_pow[96 + i] * xc[i];
            yc_pow[128 + i] = yc_pow[96 + i] * yc[i];
            zc_pow[128 + i] = zc_pow[96 + i] * zc[i];
        }
        // Combine A blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {

            phi_tmp[i] = xc_pow[128 + i] * S0[i];
            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * S0[i];
            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * S0[i];
            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * S0[i];
            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * S0[i];
            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * S0[i];
            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * S0[i];
            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * S0[i];
            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * S0[i];
            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * S0[i];
            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * S0[i];
            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * S0[i];
            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * S0[i];
            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * S0[i];
            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * S0[i];
            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * S0[i];
            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * S0[i];
            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * S0[i];
            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * S0[i];
            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * S0[i];
            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * S0[i];
            phi_tmp[672 + i] = yc_pow[128 + i] * S0[i];
            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * S0[i];
            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * S0[i];
            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * S0[i];
            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * S0[i];
            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * S0[i];
            phi_tmp[864 + i] = zc_pow[128 + i] * S0[i];
        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_out + start), npoints);
        }

        // Combine X blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];

            phi_tmp[i] = xc_pow[128 + i] * SX;
            phi_tmp[i] += 6.0 * xc_pow[96 + i] * S0[i];

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SX;
            phi_tmp[32 + i] += 5.0 * xc_pow[64 + i] * yc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SX;
            phi_tmp[64 + i] += 5.0 * xc_pow[64 + i] * zc[i] * S0[i];

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SX;
            phi_tmp[96 + i] += 4.0 * xc_pow[32 + i] * yc_pow[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SX;
            phi_tmp[128 + i] += 4.0 * xc_pow[32 + i] * yc[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SX;
            phi_tmp[160 + i] += 4.0 * xc_pow[32 + i] * zc_pow[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SX;
            phi_tmp[192 + i] += 3.0 * xc_pow[i] * yc_pow[32 + i] * S0[i];

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[224 + i] += 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[256 + i] += 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SX;
            phi_tmp[288 + i] += 3.0 * xc_pow[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SX;
            phi_tmp[320 + i] += 2.0 * xc[i] * yc_pow[64 + i] * S0[i];

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SX;
            phi_tmp[352 + i] += 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SX;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SX;
            phi_tmp[416 + i] += 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SX;
            phi_tmp[448 + i] += 2.0 * xc[i] * zc_pow[64 + i] * S0[i];

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SX;
            phi_tmp[480 + i] += yc_pow[96 + i] * S0[i];

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SX;
            phi_tmp[512 + i] += yc_pow[64 + i] * zc[i] * S0[i];

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SX;
            phi_tmp[544 + i] += yc_pow[32 + i] * zc_pow[i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SX;
            phi_tmp[576 + i] += yc_pow[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SX;
            phi_tmp[608 + i] += yc[i] * zc_pow[64 + i] * S0[i];

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SX;
            phi_tmp[640 + i] += zc_pow[96 + i] * S0[i];

            phi_tmp[672 + i] = yc_pow[128 + i] * SX;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SX;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SX;

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SX;

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SX;

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SX;

            phi_tmp[864 + i] = zc_pow[128 + i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_x_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_x_out + start), npoints);
        }

        // Combine Y blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];

            phi_tmp[i] = xc_pow[128 + i] * SY;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SY;
            phi_tmp[32 + i] += xc_pow[96 + i] * S0[i];

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SY;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SY;
            phi_tmp[96 + i] += 2.0 * xc_pow[64 + i] * yc[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SY;
            phi_tmp[128 + i] += xc_pow[64 + i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SY;

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SY;
            phi_tmp[192 + i] += 3.0 * xc_pow[32 + i] * yc_pow[i] * S0[i];

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[256 + i] += xc_pow[32 + i] * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SY;

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SY;
            phi_tmp[320 + i] += 4.0 * xc_pow[i] * yc_pow[32 + i] * S0[i];

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[352 + i] += 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[416 + i] += xc_pow[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SY;

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SY;
            phi_tmp[480 + i] += 5.0 * xc[i] * yc_pow[64 + i] * S0[i];

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SY;
            phi_tmp[512 + i] += 4.0 * xc[i] * yc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SY;
            phi_tmp[544 + i] += 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SY;
            phi_tmp[576 + i] += 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SY;
            phi_tmp[608 + i] += xc[i] * zc_pow[64 + i] * S0[i];

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SY;

            phi_tmp[672 + i] = yc_pow[128 + i] * SY;
            phi_tmp[672 + i] += 6.0 * yc_pow[96 + i] * S0[i];

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SY;
            phi_tmp[704 + i] += 5.0 * yc_pow[64 + i] * zc[i] * S0[i];

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SY;
            phi_tmp[736 + i] += 4.0 * yc_pow[32 + i] * zc_pow[i] * S0[i];

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SY;
            phi_tmp[768 + i] += 3.0 * yc_pow[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SY;
            phi_tmp[800 + i] += 2.0 * yc[i] * zc_pow[64 + i] * S0[i];

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SY;
            phi_tmp[832 + i] += zc_pow[96 + i] * S0[i];

            phi_tmp[864 + i] = zc_pow[128 + i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_y_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_y_out + start), npoints);
        }

        // Combine Z blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];

            phi_tmp[i] = xc_pow[128 + i] * SZ;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SZ;
            phi_tmp[64 + i] += xc_pow[96 + i] * S0[i];

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SZ;

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += xc_pow[64 + i] * yc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[64 + i] * zc[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SZ;

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[224 + i] += xc_pow[32 + i] * yc_pow[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[32 + i] * zc_pow[i] * S0[i];

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SZ;

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[352 + i] += xc_pow[i] * yc_pow[32 + i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[416 + i] += 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SZ;
            phi_tmp[448 + i] += 4.0 * xc_pow[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SZ;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SZ;
            phi_tmp[512 + i] += xc[i] * yc_pow[64 + i] * S0[i];

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SZ;
            phi_tmp[544 + i] += 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SZ;
            phi_tmp[576 + i] += 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SZ;
            phi_tmp[608 + i] += 4.0 * xc[i] * yc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SZ;
            phi_tmp[640 + i] += 5.0 * xc[i] * zc_pow[64 + i] * S0[i];

            phi_tmp[672 + i] = yc_pow[128 + i] * SZ;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SZ;
            phi_tmp[704 + i] += yc_pow[96 + i] * S0[i];

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SZ;
            phi_tmp[736 + i] += 2.0 * yc_pow[64 + i] * zc[i] * S0[i];

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SZ;
            phi_tmp[768 + i] += 3.0 * yc_pow[32 + i] * zc_pow[i] * S0[i];

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SZ;
            phi_tmp[800 + i] += 4.0 * yc_pow[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SZ;
            phi_tmp[832 + i] += 5.0 * yc[i] * zc_pow[64 + i] * S0[i];

            phi_tmp[864 + i] = zc_pow[128 + i] * SZ;
            phi_tmp[864 + i] += 6.0 * zc_pow[96 + i] * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_z_out + start), npoints);
        }

        // Combine XX blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];

            phi_tmp[i] = xc_pow[128 + i] * SXX;
            phi_tmp[i] += 12.0 * xc_pow[96 + i] * SX;
            phi_tmp[i] += 30.0 * xc_pow[64 + i] * S0[i];

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SXX;
            phi_tmp[32 + i] += 10.0 * xc_pow[64 + i] * yc[i] * SX;
            phi_tmp[32 + i] += 20.0 * xc_pow[32 + i] * yc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SXX;
            phi_tmp[64 + i] += 10.0 * xc_pow[64 + i] * zc[i] * SX;
            phi_tmp[64 + i] += 20.0 * xc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SXX;
            phi_tmp[96 + i] += 8.0 * xc_pow[32 + i] * yc_pow[i] * SX;
            phi_tmp[96 + i] += 12.0 * xc_pow[i] * yc_pow[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SXX;
            phi_tmp[128 + i] += 8.0 * xc_pow[32 + i] * yc[i] * zc[i] * SX;
            phi_tmp[128 + i] += 12.0 * xc_pow[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SXX;
            phi_tmp[160 + i] += 8.0 * xc_pow[32 + i] * zc_pow[i] * SX;
            phi_tmp[160 + i] += 12.0 * xc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SXX;
            phi_tmp[192 + i] += 6.0 * xc_pow[i] * yc_pow[32 + i] * SX;
            phi_tmp[192 + i] += 6.0 * xc[i] * yc_pow[32 + i] * S0[i];

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SXX;
            phi_tmp[224 + i] += 6.0 * xc_pow[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[224 + i] += 6.0 * xc[i] * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SXX;
            phi_tmp[256 + i] += 6.0 * xc_pow[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[256 + i] += 6.0 * xc[i] * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SXX;
            phi_tmp[288 + i] += 6.0 * xc_pow[i] * zc_pow[32 + i] * SX;
            phi_tmp[288 + i] += 6.0 * xc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SXX;
            phi_tmp[320 + i] += 4.0 * xc[i] * yc_pow[64 + i] * SX;
            phi_tmp[320 + i] += 2.0 * yc_pow[64 + i] * S0[i];

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SXX;
            phi_tmp[352 + i] += 4.0 * xc[i] * yc_pow[32 + i] * zc[i] * SX;
            phi_tmp[352 + i] += 2.0 * yc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SXX;
            phi_tmp[384 + i] += 4.0 * xc[i] * yc_pow[i] * zc_pow[i] * SX;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SXX;
            phi_tmp[416 + i] += 4.0 * xc[i] * yc[i] * zc_pow[32 + i] * SX;
            phi_tmp[416 + i] += 2.0 * yc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SXX;
            phi_tmp[448 + i] += 4.0 * xc[i] * zc_pow[64 + i] * SX;
            phi_tmp[448 + i] += 2.0 * zc_pow[64 + i] * S0[i];

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SXX;
            phi_tmp[480 + i] += 2.0 * yc_pow[96 + i] * SX;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SXX;
            phi_tmp[512 + i] += 2.0 * yc_pow[64 + i] * zc[i] * SX;

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SXX;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * zc_pow[i] * SX;

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SXX;
            phi_tmp[576 + i] += 2.0 * yc_pow[i] * zc_pow[32 + i] * SX;

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SXX;
            phi_tmp[608 + i] += 2.0 * yc[i] * zc_pow[64 + i] * SX;

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SXX;
            phi_tmp[640 + i] += 2.0 * zc_pow[96 + i] * SX;

            phi_tmp[672 + i] = yc_pow[128 + i] * SXX;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SXX;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SXX;

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SXX;

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SXX;

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SXX;

            phi_tmp[864 + i] = zc_pow[128 + i] * SXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_xx_out + start), npoints);
        }

        // Combine XY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];

            phi_tmp[i] = xc_pow[128 + i] * SXY;
            phi_tmp[i] += 6.0 * xc_pow[96 + i] * SY;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SXY;
            phi_tmp[32 + i] += xc_pow[96 + i] * SX;
            phi_tmp[32 + i] += 5.0 * xc_pow[64 + i] * yc[i] * SY;
            phi_tmp[32 + i] += 5.0 * xc_pow[64 + i] * S0[i];

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SXY;
            phi_tmp[64 + i] += 5.0 * xc_pow[64 + i] * zc[i] * SY;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SXY;
            phi_tmp[96 + i] += 2.0 * xc_pow[64 + i] * yc[i] * SX;
            phi_tmp[96 + i] += 4.0 * xc_pow[32 + i] * yc_pow[i] * SY;
            phi_tmp[96 + i] += 8.0 * xc_pow[32 + i] * yc[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SXY;
            phi_tmp[128 + i] += xc_pow[64 + i] * zc[i] * SX;
            phi_tmp[128 + i] += 4.0 * xc_pow[32 + i] * yc[i] * zc[i] * SY;
            phi_tmp[128 + i] += 4.0 * xc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SXY;
            phi_tmp[160 + i] += 4.0 * xc_pow[32 + i] * zc_pow[i] * SY;

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SXY;
            phi_tmp[192 + i] += 3.0 * xc_pow[32 + i] * yc_pow[i] * SX;
            phi_tmp[192 + i] += 3.0 * xc_pow[i] * yc_pow[32 + i] * SY;
            phi_tmp[192 + i] += 9.0 * xc_pow[i] * yc_pow[i] * S0[i];

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SXY;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SX;
            phi_tmp[224 + i] += 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[224 + i] += 6.0 * xc_pow[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SXY;
            phi_tmp[256 + i] += xc_pow[32 + i] * zc_pow[i] * SX;
            phi_tmp[256 + i] += 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[256 + i] += 3.0 * xc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SXY;
            phi_tmp[288 + i] += 3.0 * xc_pow[i] * zc_pow[32 + i] * SY;

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SXY;
            phi_tmp[320 + i] += 4.0 * xc_pow[i] * yc_pow[32 + i] * SX;
            phi_tmp[320 + i] += 2.0 * xc[i] * yc_pow[64 + i] * SY;
            phi_tmp[320 + i] += 8.0 * xc[i] * yc_pow[32 + i] * S0[i];

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[352 + i] += 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[352 + i] += 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[352 + i] += 6.0 * xc[i] * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[384 + i] += 4.0 * xc[i] * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[416 + i] += xc_pow[i] * zc_pow[32 + i] * SX;
            phi_tmp[416 + i] += 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[416 + i] += 2.0 * xc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SXY;
            phi_tmp[448 + i] += 2.0 * xc[i] * zc_pow[64 + i] * SY;

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SXY;
            phi_tmp[480 + i] += 5.0 * xc[i] * yc_pow[64 + i] * SX;
            phi_tmp[480 + i] += yc_pow[96 + i] * SY;
            phi_tmp[480 + i] += 5.0 * yc_pow[64 + i] * S0[i];

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SXY;
            phi_tmp[512 + i] += 4.0 * xc[i] * yc_pow[32 + i] * zc[i] * SX;
            phi_tmp[512 + i] += yc_pow[64 + i] * zc[i] * SY;
            phi_tmp[512 + i] += 4.0 * yc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SXY;
            phi_tmp[544 + i] += 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SX;
            phi_tmp[544 + i] += yc_pow[32 + i] * zc_pow[i] * SY;
            phi_tmp[544 + i] += 3.0 * yc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SXY;
            phi_tmp[576 + i] += 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SX;
            phi_tmp[576 + i] += yc_pow[i] * zc_pow[32 + i] * SY;
            phi_tmp[576 + i] += 2.0 * yc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SXY;
            phi_tmp[608 + i] += xc[i] * zc_pow[64 + i] * SX;
            phi_tmp[608 + i] += yc[i] * zc_pow[64 + i] * SY;
            phi_tmp[608 + i] += zc_pow[64 + i] * S0[i];

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SXY;
            phi_tmp[640 + i] += zc_pow[96 + i] * SY;

            phi_tmp[672 + i] = yc_pow[128 + i] * SXY;
            phi_tmp[672 + i] += 6.0 * yc_pow[96 + i] * SX;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SXY;
            phi_tmp[704 + i] += 5.0 * yc_pow[64 + i] * zc[i] * SX;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SXY;
            phi_tmp[736 + i] += 4.0 * yc_pow[32 + i] * zc_pow[i] * SX;

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SXY;
            phi_tmp[768 + i] += 3.0 * yc_pow[i] * zc_pow[32 + i] * SX;

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SXY;
            phi_tmp[800 + i] += 2.0 * yc[i] * zc_pow[64 + i] * SX;

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SXY;
            phi_tmp[832 + i] += zc_pow[96 + i] * SX;

            phi_tmp[864 + i] = zc_pow[128 + i] * SXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_xy_out + start), npoints);
        }

        // Combine XZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];

            phi_tmp[i] = xc_pow[128 + i] * SXZ;
            phi_tmp[i] += 6.0 * xc_pow[96 + i] * SZ;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SXZ;
            phi_tmp[32 + i] += 5.0 * xc_pow[64 + i] * yc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SXZ;
            phi_tmp[64 + i] += xc_pow[96 + i] * SX;
            phi_tmp[64 + i] += 5.0 * xc_pow[64 + i] * zc[i] * SZ;
            phi_tmp[64 + i] += 5.0 * xc_pow[64 + i] * S0[i];

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SXZ;
            phi_tmp[96 + i] += 4.0 * xc_pow[32 + i] * yc_pow[i] * SZ;

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SXZ;
            phi_tmp[128 + i] += xc_pow[64 + i] * yc[i] * SX;
            phi_tmp[128 + i] += 4.0 * xc_pow[32 + i] * yc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += 4.0 * xc_pow[32 + i] * yc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SXZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[64 + i] * zc[i] * SX;
            phi_tmp[160 + i] += 4.0 * xc_pow[32 + i] * zc_pow[i] * SZ;
            phi_tmp[160 + i] += 8.0 * xc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SXZ;
            phi_tmp[192 + i] += 3.0 * xc_pow[i] * yc_pow[32 + i] * SZ;

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SXZ;
            phi_tmp[224 + i] += xc_pow[32 + i] * yc_pow[i] * SX;
            phi_tmp[224 + i] += 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[224 + i] += 3.0 * xc_pow[i] * yc_pow[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SX;
            phi_tmp[256 + i] += 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 6.0 * xc_pow[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SXZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[32 + i] * zc_pow[i] * SX;
            phi_tmp[288 + i] += 3.0 * xc_pow[i] * zc_pow[32 + i] * SZ;
            phi_tmp[288 + i] += 9.0 * xc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SXZ;
            phi_tmp[320 + i] += 2.0 * xc[i] * yc_pow[64 + i] * SZ;

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[352 + i] += xc_pow[i] * yc_pow[32 + i] * SX;
            phi_tmp[352 + i] += 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[352 + i] += 2.0 * xc[i] * yc_pow[32 + i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 4.0 * xc[i] * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[416 + i] += 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[416 + i] += 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[416 + i] += 6.0 * xc[i] * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SXZ;
            phi_tmp[448 + i] += 4.0 * xc_pow[i] * zc_pow[32 + i] * SX;
            phi_tmp[448 + i] += 2.0 * xc[i] * zc_pow[64 + i] * SZ;
            phi_tmp[448 + i] += 8.0 * xc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SXZ;
            phi_tmp[480 + i] += yc_pow[96 + i] * SZ;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SXZ;
            phi_tmp[512 + i] += xc[i] * yc_pow[64 + i] * SX;
            phi_tmp[512 + i] += yc_pow[64 + i] * zc[i] * SZ;
            phi_tmp[512 + i] += yc_pow[64 + i] * S0[i];

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SXZ;
            phi_tmp[544 + i] += 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SX;
            phi_tmp[544 + i] += yc_pow[32 + i] * zc_pow[i] * SZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[576 + i] += 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SX;
            phi_tmp[576 + i] += yc_pow[i] * zc_pow[32 + i] * SZ;
            phi_tmp[576 + i] += 3.0 * yc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SXZ;
            phi_tmp[608 + i] += 4.0 * xc[i] * yc[i] * zc_pow[32 + i] * SX;
            phi_tmp[608 + i] += yc[i] * zc_pow[64 + i] * SZ;
            phi_tmp[608 + i] += 4.0 * yc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SXZ;
            phi_tmp[640 + i] += 5.0 * xc[i] * zc_pow[64 + i] * SX;
            phi_tmp[640 + i] += zc_pow[96 + i] * SZ;
            phi_tmp[640 + i] += 5.0 * zc_pow[64 + i] * S0[i];

            phi_tmp[672 + i] = yc_pow[128 + i] * SXZ;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SXZ;
            phi_tmp[704 + i] += yc_pow[96 + i] * SX;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SXZ;
            phi_tmp[736 + i] += 2.0 * yc_pow[64 + i] * zc[i] * SX;

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SXZ;
            phi_tmp[768 + i] += 3.0 * yc_pow[32 + i] * zc_pow[i] * SX;

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SXZ;
            phi_tmp[800 + i] += 4.0 * yc_pow[i] * zc_pow[32 + i] * SX;

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SXZ;
            phi_tmp[832 + i] += 5.0 * yc[i] * zc_pow[64 + i] * SX;

            phi_tmp[864 + i] = zc_pow[128 + i] * SXZ;
            phi_tmp[864 + i] += 6.0 * zc_pow[96 + i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_xz_out + start), npoints);
        }

        // Combine YY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];

            phi_tmp[i] = xc_pow[128 + i] * SYY;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SYY;
            phi_tmp[32 + i] += 2.0 * xc_pow[96 + i] * SY;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SYY;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SYY;
            phi_tmp[96 + i] += 4.0 * xc_pow[64 + i] * yc[i] * SY;
            phi_tmp[96 + i] += 2.0 * xc_pow[64 + i] * S0[i];

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * xc_pow[64 + i] * zc[i] * SY;

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SYY;

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SYY;
            phi_tmp[192 + i] += 6.0 * xc_pow[32 + i] * yc_pow[i] * SY;
            phi_tmp[192 + i] += 6.0 * xc_pow[32 + i] * yc[i] * S0[i];

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SYY;
            phi_tmp[224 + i] += 4.0 * xc_pow[32 + i] * yc[i] * zc[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SYY;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * zc_pow[i] * SY;

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SYY;

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SYY;
            phi_tmp[320 + i] += 8.0 * xc_pow[i] * yc_pow[32 + i] * SY;
            phi_tmp[320 + i] += 12.0 * xc_pow[i] * yc_pow[i] * S0[i];

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SYY;
            phi_tmp[352 + i] += 6.0 * xc_pow[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[352 + i] += 6.0 * xc_pow[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SYY;
            phi_tmp[384 + i] += 4.0 * xc_pow[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SYY;
            phi_tmp[416 + i] += 2.0 * xc_pow[i] * zc_pow[32 + i] * SY;

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SYY;

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SYY;
            phi_tmp[480 + i] += 10.0 * xc[i] * yc_pow[64 + i] * SY;
            phi_tmp[480 + i] += 20.0 * xc[i] * yc_pow[32 + i] * S0[i];

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SYY;
            phi_tmp[512 + i] += 8.0 * xc[i] * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[512 + i] += 12.0 * xc[i] * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SYY;
            phi_tmp[544 + i] += 6.0 * xc[i] * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[544 + i] += 6.0 * xc[i] * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SYY;
            phi_tmp[576 + i] += 4.0 * xc[i] * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[576 + i] += 2.0 * xc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SYY;
            phi_tmp[608 + i] += 2.0 * xc[i] * zc_pow[64 + i] * SY;

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SYY;

            phi_tmp[672 + i] = yc_pow[128 + i] * SYY;
            phi_tmp[672 + i] += 12.0 * yc_pow[96 + i] * SY;
            phi_tmp[672 + i] += 30.0 * yc_pow[64 + i] * S0[i];

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SYY;
            phi_tmp[704 + i] += 10.0 * yc_pow[64 + i] * zc[i] * SY;
            phi_tmp[704 + i] += 20.0 * yc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SYY;
            phi_tmp[736 + i] += 8.0 * yc_pow[32 + i] * zc_pow[i] * SY;
            phi_tmp[736 + i] += 12.0 * yc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SYY;
            phi_tmp[768 + i] += 6.0 * yc_pow[i] * zc_pow[32 + i] * SY;
            phi_tmp[768 + i] += 6.0 * yc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SYY;
            phi_tmp[800 + i] += 4.0 * yc[i] * zc_pow[64 + i] * SY;
            phi_tmp[800 + i] += 2.0 * zc_pow[64 + i] * S0[i];

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SYY;
            phi_tmp[832 + i] += 2.0 * zc_pow[96 + i] * SY;

            phi_tmp[864 + i] = zc_pow[128 + i] * SYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_yy_out + start), npoints);
        }

        // Combine YZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];

            phi_tmp[i] = xc_pow[128 + i] * SYZ;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SYZ;
            phi_tmp[32 + i] += xc_pow[96 + i] * SZ;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SYZ;
            phi_tmp[64 + i] += xc_pow[96 + i] * SY;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * xc_pow[64 + i] * yc[i] * SZ;

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SYZ;
            phi_tmp[128 + i] += xc_pow[64 + i] * yc[i] * SY;
            phi_tmp[128 + i] += xc_pow[64 + i] * zc[i] * SZ;
            phi_tmp[128 + i] += xc_pow[64 + i] * S0[i];

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[64 + i] * zc[i] * SY;

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SYZ;
            phi_tmp[192 + i] += 3.0 * xc_pow[32 + i] * yc_pow[i] * SZ;

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SYZ;
            phi_tmp[224 + i] += xc_pow[32 + i] * yc_pow[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SZ;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * yc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SYZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SY;
            phi_tmp[256 + i] += xc_pow[32 + i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SYZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[32 + i] * zc_pow[i] * SY;

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SYZ;
            phi_tmp[320 + i] += 4.0 * xc_pow[i] * yc_pow[32 + i] * SZ;

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[352 + i] += xc_pow[i] * yc_pow[32 + i] * SY;
            phi_tmp[352 + i] += 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[352 + i] += 3.0 * xc_pow[i] * yc_pow[i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 4.0 * xc_pow[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[416 + i] += 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[416 + i] += xc_pow[i] * zc_pow[32 + i] * SZ;
            phi_tmp[416 + i] += 3.0 * xc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SYZ;
            phi_tmp[448 + i] += 4.0 * xc_pow[i] * zc_pow[32 + i] * SY;

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SYZ;
            phi_tmp[480 + i] += 5.0 * xc[i] * yc_pow[64 + i] * SZ;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SYZ;
            phi_tmp[512 + i] += xc[i] * yc_pow[64 + i] * SY;
            phi_tmp[512 + i] += 4.0 * xc[i] * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[512 + i] += 4.0 * xc[i] * yc_pow[32 + i] * S0[i];

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SYZ;
            phi_tmp[544 + i] += 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[544 + i] += 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[544 + i] += 6.0 * xc[i] * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[576 + i] += 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[576 + i] += 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[576 + i] += 6.0 * xc[i] * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SYZ;
            phi_tmp[608 + i] += 4.0 * xc[i] * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[608 + i] += xc[i] * zc_pow[64 + i] * SZ;
            phi_tmp[608 + i] += 4.0 * xc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SYZ;
            phi_tmp[640 + i] += 5.0 * xc[i] * zc_pow[64 + i] * SY;

            phi_tmp[672 + i] = yc_pow[128 + i] * SYZ;
            phi_tmp[672 + i] += 6.0 * yc_pow[96 + i] * SZ;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SYZ;
            phi_tmp[704 + i] += yc_pow[96 + i] * SY;
            phi_tmp[704 + i] += 5.0 * yc_pow[64 + i] * zc[i] * SZ;
            phi_tmp[704 + i] += 5.0 * yc_pow[64 + i] * S0[i];

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SYZ;
            phi_tmp[736 + i] += 2.0 * yc_pow[64 + i] * zc[i] * SY;
            phi_tmp[736 + i] += 4.0 * yc_pow[32 + i] * zc_pow[i] * SZ;
            phi_tmp[736 + i] += 8.0 * yc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SYZ;
            phi_tmp[768 + i] += 3.0 * yc_pow[32 + i] * zc_pow[i] * SY;
            phi_tmp[768 + i] += 3.0 * yc_pow[i] * zc_pow[32 + i] * SZ;
            phi_tmp[768 + i] += 9.0 * yc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SYZ;
            phi_tmp[800 + i] += 4.0 * yc_pow[i] * zc_pow[32 + i] * SY;
            phi_tmp[800 + i] += 2.0 * yc[i] * zc_pow[64 + i] * SZ;
            phi_tmp[800 + i] += 8.0 * yc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SYZ;
            phi_tmp[832 + i] += 5.0 * yc[i] * zc_pow[64 + i] * SY;
            phi_tmp[832 + i] += zc_pow[96 + i] * SZ;
            phi_tmp[832 + i] += 5.0 * zc_pow[64 + i] * S0[i];

            phi_tmp[864 + i] = zc_pow[128 + i] * SYZ;
            phi_tmp[864 + i] += 6.0 * zc_pow[96 + i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_yz_out + start), npoints);
        }

        // Combine ZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];

            phi_tmp[i] = xc_pow[128 + i] * SZZ;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SZZ;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SZZ;
            phi_tmp[64 + i] += 2.0 * xc_pow[96 + i] * SZ;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SZZ;

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[64 + i] * yc[i] * SZ;

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SZZ;
            phi_tmp[160 + i] += 4.0 * xc_pow[64 + i] * zc[i] * SZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[64 + i] * S0[i];

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SZZ;

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SZZ;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * yc_pow[i] * SZ;

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SZZ;
            phi_tmp[256 + i] += 4.0 * xc_pow[32 + i] * yc[i] * zc[i] * SZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * yc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SZZ;
            phi_tmp[288 + i] += 6.0 * xc_pow[32 + i] * zc_pow[i] * SZ;
            phi_tmp[288 + i] += 6.0 * xc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SZZ;

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[352 + i] += 2.0 * xc_pow[i] * yc_pow[32 + i] * SZ;

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[384 + i] += 4.0 * xc_pow[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[416 + i] += 6.0 * xc_pow[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[416 + i] += 6.0 * xc_pow[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SZZ;
            phi_tmp[448 + i] += 8.0 * xc_pow[i] * zc_pow[32 + i] * SZ;
            phi_tmp[448 + i] += 12.0 * xc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SZZ;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SZZ;
            phi_tmp[512 + i] += 2.0 * xc[i] * yc_pow[64 + i] * SZ;

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SZZ;
            phi_tmp[544 + i] += 4.0 * xc[i] * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[544 + i] += 2.0 * xc[i] * yc_pow[32 + i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[576 + i] += 6.0 * xc[i] * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[576 + i] += 6.0 * xc[i] * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SZZ;
            phi_tmp[608 + i] += 8.0 * xc[i] * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[608 + i] += 12.0 * xc[i] * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SZZ;
            phi_tmp[640 + i] += 10.0 * xc[i] * zc_pow[64 + i] * SZ;
            phi_tmp[640 + i] += 20.0 * xc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[672 + i] = yc_pow[128 + i] * SZZ;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SZZ;
            phi_tmp[704 + i] += 2.0 * yc_pow[96 + i] * SZ;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SZZ;
            phi_tmp[736 + i] += 4.0 * yc_pow[64 + i] * zc[i] * SZ;
            phi_tmp[736 + i] += 2.0 * yc_pow[64 + i] * S0[i];

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SZZ;
            phi_tmp[768 + i] += 6.0 * yc_pow[32 + i] * zc_pow[i] * SZ;
            phi_tmp[768 + i] += 6.0 * yc_pow[32 + i] * zc[i] * S0[i];

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SZZ;
            phi_tmp[800 + i] += 8.0 * yc_pow[i] * zc_pow[32 + i] * SZ;
            phi_tmp[800 + i] += 12.0 * yc_pow[i] * zc_pow[i] * S0[i];

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SZZ;
            phi_tmp[832 + i] += 10.0 * yc[i] * zc_pow[64 + i] * SZ;
            phi_tmp[832 + i] += 20.0 * yc[i] * zc_pow[32 + i] * S0[i];

            phi_tmp[864 + i] = zc_pow[128 + i] * SZZ;
            phi_tmp[864 + i] += 12.0 * zc_pow[96 + i] * SZ;
            phi_tmp[864 + i] += 30.0 * zc_pow[64 + i] * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_zz_out + start), npoints);
        }

        // Combine XXX blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXX = S3[i] * xc[i] * xc[i] * xc[i] + 3 * xc[i] * S2[i];

            phi_tmp[i] = xc_pow[128 + i] * SXXX;
            phi_tmp[i] += 3.0 * 6.0 * xc_pow[96 + i] * SXX;
            phi_tmp[i] += 3.0 * 30.0 * xc_pow[64 + i] * SX;
            phi_tmp[i] += 120.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SXXX;
            phi_tmp[32 + i] += 3.0 * 5.0 * xc_pow[64 + i] * yc[i] * SXX;
            phi_tmp[32 + i] += 3.0 * 20.0 * xc_pow[32 + i] * yc[i] * SX;
            phi_tmp[32 + i] += 60.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SXXX;
            phi_tmp[64 + i] += 3.0 * 5.0 * xc_pow[64 + i] * zc[i] * SXX;
            phi_tmp[64 + i] += 3.0 * 20.0 * xc_pow[32 + i] * zc[i] * SX;
            phi_tmp[64 + i] += 60.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SXXX;
            phi_tmp[96 + i] += 3.0 * 4.0 * xc_pow[32 + i] * yc_pow[i] * SXX;
            phi_tmp[96 + i] += 3.0 * 12.0 * xc_pow[i] * yc_pow[i] * SX;
            phi_tmp[96 + i] += 24.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SXXX;
            phi_tmp[128 + i] += 3.0 * 4.0 * xc_pow[32 + i] * yc[i] * zc[i] * SXX;
            phi_tmp[128 + i] += 3.0 * 12.0 * xc_pow[i] * yc[i] * zc[i] * SX;
            phi_tmp[128 + i] += 24.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SXXX;
            phi_tmp[160 + i] += 3.0 * 4.0 * xc_pow[32 + i] * zc_pow[i] * SXX;
            phi_tmp[160 + i] += 3.0 * 12.0 * xc_pow[i] * zc_pow[i] * SX;
            phi_tmp[160 + i] += 24.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SXXX;
            phi_tmp[192 + i] += 3.0 * 3.0 * xc_pow[i] * yc_pow[32 + i] * SXX;
            phi_tmp[192 + i] += 3.0 * 6.0 * xc[i] * yc_pow[32 + i] * SX;
            phi_tmp[192 + i] += 6.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SXXX;
            phi_tmp[224 + i] += 3.0 * 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SXX;
            phi_tmp[224 + i] += 3.0 * 6.0 * xc[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[224 + i] += 6.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SXXX;
            phi_tmp[256 + i] += 3.0 * 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SXX;
            phi_tmp[256 + i] += 3.0 * 6.0 * xc[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[256 + i] += 6.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SXXX;
            phi_tmp[288 + i] += 3.0 * 3.0 * xc_pow[i] * zc_pow[32 + i] * SXX;
            phi_tmp[288 + i] += 3.0 * 6.0 * xc[i] * zc_pow[32 + i] * SX;
            phi_tmp[288 + i] += 6.0 * zc_pow[32 + i] * S0[i];

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SXXX;
            phi_tmp[320 + i] += 3.0 * 2.0 * xc[i] * yc_pow[64 + i] * SXX;
            phi_tmp[320 + i] += 3.0 * 2.0 * yc_pow[64 + i] * SX;

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SXXX;
            phi_tmp[352 + i] += 3.0 * 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SXX;
            phi_tmp[352 + i] += 3.0 * 2.0 * yc_pow[32 + i] * zc[i] * SX;

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SXXX;
            phi_tmp[384 + i] += 3.0 * 2.0 * xc[i] * yc_pow[i] * zc_pow[i] * SXX;
            phi_tmp[384 + i] += 3.0 * 2.0 * yc_pow[i] * zc_pow[i] * SX;

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SXXX;
            phi_tmp[416 + i] += 3.0 * 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SXX;
            phi_tmp[416 + i] += 3.0 * 2.0 * yc[i] * zc_pow[32 + i] * SX;

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SXXX;
            phi_tmp[448 + i] += 3.0 * 2.0 * xc[i] * zc_pow[64 + i] * SXX;
            phi_tmp[448 + i] += 3.0 * 2.0 * zc_pow[64 + i] * SX;

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SXXX;
            phi_tmp[480 + i] += 3.0 * yc_pow[96 + i] * SXX;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SXXX;
            phi_tmp[512 + i] += 3.0 * yc_pow[64 + i] * zc[i] * SXX;

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SXXX;
            phi_tmp[544 + i] += 3.0 * yc_pow[32 + i] * zc_pow[i] * SXX;

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SXXX;
            phi_tmp[576 + i] += 3.0 * yc_pow[i] * zc_pow[32 + i] * SXX;

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SXXX;
            phi_tmp[608 + i] += 3.0 * yc[i] * zc_pow[64 + i] * SXX;

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SXXX;
            phi_tmp[640 + i] += 3.0 * zc_pow[96 + i] * SXX;

            phi_tmp[672 + i] = yc_pow[128 + i] * SXXX;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SXXX;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SXXX;

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SXXX;

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SXXX;

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SXXX;

            phi_tmp[864 + i] = zc_pow[128 + i] * SXXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_xxx_out + start), npoints);
        }

        // Combine XXY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXY = S3[i] * xc[i] * xc[i] * yc[i] + S2[i] * yc[i];

            phi_tmp[i] = xc_pow[128 + i] * SXXY;
            phi_tmp[i] += 2.0 * 6.0 * xc_pow[96 + i] * SXY;
            phi_tmp[i] += 30.0 * xc_pow[64 + i] * SY;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SXXY;
            phi_tmp[32 + i] += 2.0 * 5.0 * xc_pow[64 + i] * yc[i] * SXY;
            phi_tmp[32 + i] += xc_pow[96 + i] * SXX;
            phi_tmp[32 + i] += 20.0 * xc_pow[32 + i] * yc[i] * SY;
            phi_tmp[32 + i] += 2.0 * 5.0 * xc_pow[64 + i] * SX;
            phi_tmp[32 + i] += 20.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SXXY;
            phi_tmp[64 + i] += 2.0 * 5.0 * xc_pow[64 + i] * zc[i] * SXY;
            phi_tmp[64 + i] += 20.0 * xc_pow[32 + i] * zc[i] * SY;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SXXY;
            phi_tmp[96 + i] += 2.0 * 4.0 * xc_pow[32 + i] * yc_pow[i] * SXY;
            phi_tmp[96 + i] += 2.0 * xc_pow[64 + i] * yc[i] * SXX;
            phi_tmp[96 + i] += 12.0 * xc_pow[i] * yc_pow[i] * SY;
            phi_tmp[96 + i] += 2.0 * 8.0 * xc_pow[32 + i] * yc[i] * SX;
            phi_tmp[96 + i] += 24.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SXXY;
            phi_tmp[128 + i] += 2.0 * 4.0 * xc_pow[32 + i] * yc[i] * zc[i] * SXY;
            phi_tmp[128 + i] += xc_pow[64 + i] * zc[i] * SXX;
            phi_tmp[128 + i] += 12.0 * xc_pow[i] * yc[i] * zc[i] * SY;
            phi_tmp[128 + i] += 2.0 * 4.0 * xc_pow[32 + i] * zc[i] * SX;
            phi_tmp[128 + i] += 12.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SXXY;
            phi_tmp[160 + i] += 2.0 * 4.0 * xc_pow[32 + i] * zc_pow[i] * SXY;
            phi_tmp[160 + i] += 12.0 * xc_pow[i] * zc_pow[i] * SY;

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SXXY;
            phi_tmp[192 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[32 + i] * SXY;
            phi_tmp[192 + i] += 3.0 * xc_pow[32 + i] * yc_pow[i] * SXX;
            phi_tmp[192 + i] += 6.0 * xc[i] * yc_pow[32 + i] * SY;
            phi_tmp[192 + i] += 2.0 * 9.0 * xc_pow[i] * yc_pow[i] * SX;
            phi_tmp[192 + i] += 18.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SXXY;
            phi_tmp[224 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SXY;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SXX;
            phi_tmp[224 + i] += 6.0 * xc[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[224 + i] += 2.0 * 6.0 * xc_pow[i] * yc[i] * zc[i] * SX;
            phi_tmp[224 + i] += 12.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SXXY;
            phi_tmp[256 + i] += 2.0 * 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SXY;
            phi_tmp[256 + i] += xc_pow[32 + i] * zc_pow[i] * SXX;
            phi_tmp[256 + i] += 6.0 * xc[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[256 + i] += 2.0 * 3.0 * xc_pow[i] * zc_pow[i] * SX;
            phi_tmp[256 + i] += 6.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SXXY;
            phi_tmp[288 + i] += 2.0 * 3.0 * xc_pow[i] * zc_pow[32 + i] * SXY;
            phi_tmp[288 + i] += 6.0 * xc[i] * zc_pow[32 + i] * SY;

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SXXY;
            phi_tmp[320 + i] += 2.0 * 2.0 * xc[i] * yc_pow[64 + i] * SXY;
            phi_tmp[320 + i] += 4.0 * xc_pow[i] * yc_pow[32 + i] * SXX;
            phi_tmp[320 + i] += 2.0 * yc_pow[64 + i] * SY;
            phi_tmp[320 + i] += 2.0 * 8.0 * xc[i] * yc_pow[32 + i] * SX;
            phi_tmp[320 + i] += 8.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SXXY;
            phi_tmp[352 + i] += 2.0 * 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[352 + i] += 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SXX;
            phi_tmp[352 + i] += 2.0 * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[352 + i] += 2.0 * 6.0 * xc[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[352 + i] += 6.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SXXY;
            phi_tmp[384 + i] += 2.0 * 2.0 * xc[i] * yc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc[i] * zc_pow[i] * SXX;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[384 + i] += 2.0 * 4.0 * xc[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[384 + i] += 4.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SXXY;
            phi_tmp[416 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[416 + i] += xc_pow[i] * zc_pow[32 + i] * SXX;
            phi_tmp[416 + i] += 2.0 * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[416 + i] += 2.0 * 2.0 * xc[i] * zc_pow[32 + i] * SX;
            phi_tmp[416 + i] += 2.0 * zc_pow[32 + i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SXXY;
            phi_tmp[448 + i] += 2.0 * 2.0 * xc[i] * zc_pow[64 + i] * SXY;
            phi_tmp[448 + i] += 2.0 * zc_pow[64 + i] * SY;

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SXXY;
            phi_tmp[480 + i] += 2.0 * yc_pow[96 + i] * SXY;
            phi_tmp[480 + i] += 5.0 * xc[i] * yc_pow[64 + i] * SXX;
            phi_tmp[480 + i] += 2.0 * 5.0 * yc_pow[64 + i] * SX;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SXXY;
            phi_tmp[512 + i] += 2.0 * yc_pow[64 + i] * zc[i] * SXY;
            phi_tmp[512 + i] += 4.0 * xc[i] * yc_pow[32 + i] * zc[i] * SXX;
            phi_tmp[512 + i] += 2.0 * 4.0 * yc_pow[32 + i] * zc[i] * SX;

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SXXY;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * zc_pow[i] * SXY;
            phi_tmp[544 + i] += 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SXX;
            phi_tmp[544 + i] += 2.0 * 3.0 * yc_pow[i] * zc_pow[i] * SX;

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SXXY;
            phi_tmp[576 + i] += 2.0 * yc_pow[i] * zc_pow[32 + i] * SXY;
            phi_tmp[576 + i] += 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SXX;
            phi_tmp[576 + i] += 2.0 * 2.0 * yc[i] * zc_pow[32 + i] * SX;

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SXXY;
            phi_tmp[608 + i] += 2.0 * yc[i] * zc_pow[64 + i] * SXY;
            phi_tmp[608 + i] += xc[i] * zc_pow[64 + i] * SXX;
            phi_tmp[608 + i] += 2.0 * zc_pow[64 + i] * SX;

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SXXY;
            phi_tmp[640 + i] += 2.0 * zc_pow[96 + i] * SXY;

            phi_tmp[672 + i] = yc_pow[128 + i] * SXXY;
            phi_tmp[672 + i] += 6.0 * yc_pow[96 + i] * SXX;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SXXY;
            phi_tmp[704 + i] += 5.0 * yc_pow[64 + i] * zc[i] * SXX;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SXXY;
            phi_tmp[736 + i] += 4.0 * yc_pow[32 + i] * zc_pow[i] * SXX;

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SXXY;
            phi_tmp[768 + i] += 3.0 * yc_pow[i] * zc_pow[32 + i] * SXX;

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SXXY;
            phi_tmp[800 + i] += 2.0 * yc[i] * zc_pow[64 + i] * SXX;

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SXXY;
            phi_tmp[832 + i] += zc_pow[96 + i] * SXX;

            phi_tmp[864 + i] = zc_pow[128 + i] * SXXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_xxy_out + start), npoints);
        }

        // Combine XXZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SXX = S2[i] * xc[i] * xc[i] + S1[i];
            const double SXXZ = S3[i] * xc[i] * xc[i] * zc[i] + S2[i] * zc[i];

            phi_tmp[i] = xc_pow[128 + i] * SXXZ;
            phi_tmp[i] += 2.0 * 6.0 * xc_pow[96 + i] * SXZ;
            phi_tmp[i] += 30.0 * xc_pow[64 + i] * SZ;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SXXZ;
            phi_tmp[32 + i] += 2.0 * 5.0 * xc_pow[64 + i] * yc[i] * SXZ;
            phi_tmp[32 + i] += 20.0 * xc_pow[32 + i] * yc[i] * SZ;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SXXZ;
            phi_tmp[64 + i] += 2.0 * 5.0 * xc_pow[64 + i] * zc[i] * SXZ;
            phi_tmp[64 + i] += xc_pow[96 + i] * SXX;
            phi_tmp[64 + i] += 20.0 * xc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[64 + i] += 2.0 * 5.0 * xc_pow[64 + i] * SX;
            phi_tmp[64 + i] += 20.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SXXZ;
            phi_tmp[96 + i] += 2.0 * 4.0 * xc_pow[32 + i] * yc_pow[i] * SXZ;
            phi_tmp[96 + i] += 12.0 * xc_pow[i] * yc_pow[i] * SZ;

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SXXZ;
            phi_tmp[128 + i] += 2.0 * 4.0 * xc_pow[32 + i] * yc[i] * zc[i] * SXZ;
            phi_tmp[128 + i] += xc_pow[64 + i] * yc[i] * SXX;
            phi_tmp[128 + i] += 12.0 * xc_pow[i] * yc[i] * zc[i] * SZ;
            phi_tmp[128 + i] += 2.0 * 4.0 * xc_pow[32 + i] * yc[i] * SX;
            phi_tmp[128 + i] += 12.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SXXZ;
            phi_tmp[160 + i] += 2.0 * 4.0 * xc_pow[32 + i] * zc_pow[i] * SXZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[64 + i] * zc[i] * SXX;
            phi_tmp[160 + i] += 12.0 * xc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[160 + i] += 2.0 * 8.0 * xc_pow[32 + i] * zc[i] * SX;
            phi_tmp[160 + i] += 24.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SXXZ;
            phi_tmp[192 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[32 + i] * SXZ;
            phi_tmp[192 + i] += 6.0 * xc[i] * yc_pow[32 + i] * SZ;

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SXXZ;
            phi_tmp[224 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SXZ;
            phi_tmp[224 + i] += xc_pow[32 + i] * yc_pow[i] * SXX;
            phi_tmp[224 + i] += 6.0 * xc[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[224 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[i] * SX;
            phi_tmp[224 + i] += 6.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SXXZ;
            phi_tmp[256 + i] += 2.0 * 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SXX;
            phi_tmp[256 + i] += 6.0 * xc[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 2.0 * 6.0 * xc_pow[i] * yc[i] * zc[i] * SX;
            phi_tmp[256 + i] += 12.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SXXZ;
            phi_tmp[288 + i] += 2.0 * 3.0 * xc_pow[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[32 + i] * zc_pow[i] * SXX;
            phi_tmp[288 + i] += 6.0 * xc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[288 + i] += 2.0 * 9.0 * xc_pow[i] * zc_pow[i] * SX;
            phi_tmp[288 + i] += 18.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SXXZ;
            phi_tmp[320 + i] += 2.0 * 2.0 * xc[i] * yc_pow[64 + i] * SXZ;
            phi_tmp[320 + i] += 2.0 * yc_pow[64 + i] * SZ;

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SXXZ;
            phi_tmp[352 + i] += 2.0 * 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[352 + i] += xc_pow[i] * yc_pow[32 + i] * SXX;
            phi_tmp[352 + i] += 2.0 * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[352 + i] += 2.0 * 2.0 * xc[i] * yc_pow[32 + i] * SX;
            phi_tmp[352 + i] += 2.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SXXZ;
            phi_tmp[384 + i] += 2.0 * 2.0 * xc[i] * yc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc_pow[i] * zc[i] * SXX;
            phi_tmp[384 + i] += 2.0 * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 2.0 * 4.0 * xc[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[384 + i] += 4.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SXXZ;
            phi_tmp[416 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[416 + i] += 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SXX;
            phi_tmp[416 + i] += 2.0 * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[416 + i] += 2.0 * 6.0 * xc[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[416 + i] += 6.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SXXZ;
            phi_tmp[448 + i] += 2.0 * 2.0 * xc[i] * zc_pow[64 + i] * SXZ;
            phi_tmp[448 + i] += 4.0 * xc_pow[i] * zc_pow[32 + i] * SXX;
            phi_tmp[448 + i] += 2.0 * zc_pow[64 + i] * SZ;
            phi_tmp[448 + i] += 2.0 * 8.0 * xc[i] * zc_pow[32 + i] * SX;
            phi_tmp[448 + i] += 8.0 * zc_pow[32 + i] * S0[i];

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SXXZ;
            phi_tmp[480 + i] += 2.0 * yc_pow[96 + i] * SXZ;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SXXZ;
            phi_tmp[512 + i] += 2.0 * yc_pow[64 + i] * zc[i] * SXZ;
            phi_tmp[512 + i] += xc[i] * yc_pow[64 + i] * SXX;
            phi_tmp[512 + i] += 2.0 * yc_pow[64 + i] * SX;

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SXXZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * zc_pow[i] * SXZ;
            phi_tmp[544 + i] += 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SXX;
            phi_tmp[544 + i] += 2.0 * 2.0 * yc_pow[32 + i] * zc[i] * SX;

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SXXZ;
            phi_tmp[576 + i] += 2.0 * yc_pow[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[576 + i] += 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SXX;
            phi_tmp[576 + i] += 2.0 * 3.0 * yc_pow[i] * zc_pow[i] * SX;

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SXXZ;
            phi_tmp[608 + i] += 2.0 * yc[i] * zc_pow[64 + i] * SXZ;
            phi_tmp[608 + i] += 4.0 * xc[i] * yc[i] * zc_pow[32 + i] * SXX;
            phi_tmp[608 + i] += 2.0 * 4.0 * yc[i] * zc_pow[32 + i] * SX;

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SXXZ;
            phi_tmp[640 + i] += 2.0 * zc_pow[96 + i] * SXZ;
            phi_tmp[640 + i] += 5.0 * xc[i] * zc_pow[64 + i] * SXX;
            phi_tmp[640 + i] += 2.0 * 5.0 * zc_pow[64 + i] * SX;

            phi_tmp[672 + i] = yc_pow[128 + i] * SXXZ;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SXXZ;
            phi_tmp[704 + i] += yc_pow[96 + i] * SXX;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SXXZ;
            phi_tmp[736 + i] += 2.0 * yc_pow[64 + i] * zc[i] * SXX;

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SXXZ;
            phi_tmp[768 + i] += 3.0 * yc_pow[32 + i] * zc_pow[i] * SXX;

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SXXZ;
            phi_tmp[800 + i] += 4.0 * yc_pow[i] * zc_pow[32 + i] * SXX;

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SXXZ;
            phi_tmp[832 + i] += 5.0 * yc[i] * zc_pow[64 + i] * SXX;

            phi_tmp[864 + i] = zc_pow[128 + i] * SXXZ;
            phi_tmp[864 + i] += 6.0 * zc_pow[96 + i] * SXX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_xxz_out + start), npoints);
        }

        // Combine XYY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SXYY = S3[i] * xc[i] * yc[i] * yc[i] + S2[i] * xc[i];

            phi_tmp[i] = xc_pow[128 + i] * SXYY;
            phi_tmp[i] += 6.0 * xc_pow[96 + i] * SYY;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SXYY;
            phi_tmp[32 + i] += 2.0 * xc_pow[96 + i] * SXY;
            phi_tmp[32 + i] += 5.0 * xc_pow[64 + i] * yc[i] * SYY;
            phi_tmp[32 + i] += 2.0 * 5.0 * xc_pow[64 + i] * SY;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SXYY;
            phi_tmp[64 + i] += 5.0 * xc_pow[64 + i] * zc[i] * SYY;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SXYY;
            phi_tmp[96 + i] += 2.0 * 2.0 * xc_pow[64 + i] * yc[i] * SXY;
            phi_tmp[96 + i] += 4.0 * xc_pow[32 + i] * yc_pow[i] * SYY;
            phi_tmp[96 + i] += 2.0 * xc_pow[64 + i] * SX;
            phi_tmp[96 + i] += 2.0 * 8.0 * xc_pow[32 + i] * yc[i] * SY;
            phi_tmp[96 + i] += 8.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SXYY;
            phi_tmp[128 + i] += 2.0 * xc_pow[64 + i] * zc[i] * SXY;
            phi_tmp[128 + i] += 4.0 * xc_pow[32 + i] * yc[i] * zc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * 4.0 * xc_pow[32 + i] * zc[i] * SY;

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SXYY;
            phi_tmp[160 + i] += 4.0 * xc_pow[32 + i] * zc_pow[i] * SYY;

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SXYY;
            phi_tmp[192 + i] += 2.0 * 3.0 * xc_pow[32 + i] * yc_pow[i] * SXY;
            phi_tmp[192 + i] += 3.0 * xc_pow[i] * yc_pow[32 + i] * SYY;
            phi_tmp[192 + i] += 6.0 * xc_pow[32 + i] * yc[i] * SX;
            phi_tmp[192 + i] += 2.0 * 9.0 * xc_pow[i] * yc_pow[i] * SY;
            phi_tmp[192 + i] += 18.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SXYY;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SXY;
            phi_tmp[224 + i] += 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SYY;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * zc[i] * SX;
            phi_tmp[224 + i] += 2.0 * 6.0 * xc_pow[i] * yc[i] * zc[i] * SY;
            phi_tmp[224 + i] += 6.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SXYY;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * zc_pow[i] * SXY;
            phi_tmp[256 + i] += 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SYY;
            phi_tmp[256 + i] += 2.0 * 3.0 * xc_pow[i] * zc_pow[i] * SY;

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SXYY;
            phi_tmp[288 + i] += 3.0 * xc_pow[i] * zc_pow[32 + i] * SYY;

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SXYY;
            phi_tmp[320 + i] += 2.0 * 4.0 * xc_pow[i] * yc_pow[32 + i] * SXY;
            phi_tmp[320 + i] += 2.0 * xc[i] * yc_pow[64 + i] * SYY;
            phi_tmp[320 + i] += 12.0 * xc_pow[i] * yc_pow[i] * SX;
            phi_tmp[320 + i] += 2.0 * 8.0 * xc[i] * yc_pow[32 + i] * SY;
            phi_tmp[320 + i] += 24.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SXYY;
            phi_tmp[352 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SXY;
            phi_tmp[352 + i] += 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SYY;
            phi_tmp[352 + i] += 6.0 * xc_pow[i] * yc[i] * zc[i] * SX;
            phi_tmp[352 + i] += 2.0 * 6.0 * xc[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[352 + i] += 12.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SXYY;
            phi_tmp[384 + i] += 2.0 * 2.0 * xc_pow[i] * yc[i] * zc_pow[i] * SXY;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * zc_pow[i] * SYY;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * zc_pow[i] * SX;
            phi_tmp[384 + i] += 2.0 * 4.0 * xc[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[384 + i] += 4.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SXYY;
            phi_tmp[416 + i] += 2.0 * xc_pow[i] * zc_pow[32 + i] * SXY;
            phi_tmp[416 + i] += 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SYY;
            phi_tmp[416 + i] += 2.0 * 2.0 * xc[i] * zc_pow[32 + i] * SY;

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SXYY;
            phi_tmp[448 + i] += 2.0 * xc[i] * zc_pow[64 + i] * SYY;

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SXYY;
            phi_tmp[480 + i] += 2.0 * 5.0 * xc[i] * yc_pow[64 + i] * SXY;
            phi_tmp[480 + i] += yc_pow[96 + i] * SYY;
            phi_tmp[480 + i] += 20.0 * xc[i] * yc_pow[32 + i] * SX;
            phi_tmp[480 + i] += 2.0 * 5.0 * yc_pow[64 + i] * SY;
            phi_tmp[480 + i] += 20.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SXYY;
            phi_tmp[512 + i] += 2.0 * 4.0 * xc[i] * yc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[512 + i] += yc_pow[64 + i] * zc[i] * SYY;
            phi_tmp[512 + i] += 12.0 * xc[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[512 + i] += 2.0 * 4.0 * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[512 + i] += 12.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SXYY;
            phi_tmp[544 + i] += 2.0 * 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[544 + i] += yc_pow[32 + i] * zc_pow[i] * SYY;
            phi_tmp[544 + i] += 6.0 * xc[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[544 + i] += 2.0 * 3.0 * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[544 + i] += 6.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SXYY;
            phi_tmp[576 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[576 + i] += yc_pow[i] * zc_pow[32 + i] * SYY;
            phi_tmp[576 + i] += 2.0 * xc[i] * zc_pow[32 + i] * SX;
            phi_tmp[576 + i] += 2.0 * 2.0 * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[576 + i] += 2.0 * zc_pow[32 + i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SXYY;
            phi_tmp[608 + i] += 2.0 * xc[i] * zc_pow[64 + i] * SXY;
            phi_tmp[608 + i] += yc[i] * zc_pow[64 + i] * SYY;
            phi_tmp[608 + i] += 2.0 * zc_pow[64 + i] * SY;

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SXYY;
            phi_tmp[640 + i] += zc_pow[96 + i] * SYY;

            phi_tmp[672 + i] = yc_pow[128 + i] * SXYY;
            phi_tmp[672 + i] += 2.0 * 6.0 * yc_pow[96 + i] * SXY;
            phi_tmp[672 + i] += 30.0 * yc_pow[64 + i] * SX;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SXYY;
            phi_tmp[704 + i] += 2.0 * 5.0 * yc_pow[64 + i] * zc[i] * SXY;
            phi_tmp[704 + i] += 20.0 * yc_pow[32 + i] * zc[i] * SX;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SXYY;
            phi_tmp[736 + i] += 2.0 * 4.0 * yc_pow[32 + i] * zc_pow[i] * SXY;
            phi_tmp[736 + i] += 12.0 * yc_pow[i] * zc_pow[i] * SX;

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SXYY;
            phi_tmp[768 + i] += 2.0 * 3.0 * yc_pow[i] * zc_pow[32 + i] * SXY;
            phi_tmp[768 + i] += 6.0 * yc[i] * zc_pow[32 + i] * SX;

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SXYY;
            phi_tmp[800 + i] += 2.0 * 2.0 * yc[i] * zc_pow[64 + i] * SXY;
            phi_tmp[800 + i] += 2.0 * zc_pow[64 + i] * SX;

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SXYY;
            phi_tmp[832 + i] += 2.0 * zc_pow[96 + i] * SXY;

            phi_tmp[864 + i] = zc_pow[128 + i] * SXYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_xyy_out + start), npoints);
        }

        // Combine XYZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SXY = S2[i] * xc[i] * yc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SXYZ = S3[i] * xc[i] * yc[i] * zc[i];

            phi_tmp[i] = xc_pow[128 + i] * SXYZ;
            phi_tmp[i] += 6.0 * xc_pow[96 + i] * SYZ;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SXYZ;
            phi_tmp[32 + i] += 5.0 * xc_pow[64 + i] * yc[i] * SYZ;
            phi_tmp[32 + i] += xc_pow[96 + i] * SXZ;
            phi_tmp[32 + i] += 5.0 * xc_pow[64 + i] * SZ;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SXYZ;
            phi_tmp[64 + i] += 5.0 * xc_pow[64 + i] * zc[i] * SYZ;
            phi_tmp[64 + i] += xc_pow[96 + i] * SXY;
            phi_tmp[64 + i] += 5.0 * xc_pow[64 + i] * SY;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SXYZ;
            phi_tmp[96 + i] += 4.0 * xc_pow[32 + i] * yc_pow[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * xc_pow[64 + i] * yc[i] * SXZ;
            phi_tmp[96 + i] += 8.0 * xc_pow[32 + i] * yc[i] * SZ;

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SXYZ;
            phi_tmp[128 + i] += 4.0 * xc_pow[32 + i] * yc[i] * zc[i] * SYZ;
            phi_tmp[128 + i] += xc_pow[64 + i] * zc[i] * SXZ;
            phi_tmp[128 + i] += xc_pow[64 + i] * yc[i] * SXY;
            phi_tmp[128 + i] += 4.0 * xc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[128 + i] += 4.0 * xc_pow[32 + i] * yc[i] * SY;
            phi_tmp[128 + i] += xc_pow[64 + i] * SX;
            phi_tmp[128 + i] += 4.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SXYZ;
            phi_tmp[160 + i] += 4.0 * xc_pow[32 + i] * zc_pow[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[64 + i] * zc[i] * SXY;
            phi_tmp[160 + i] += 8.0 * xc_pow[32 + i] * zc[i] * SY;

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SXYZ;
            phi_tmp[192 + i] += 3.0 * xc_pow[i] * yc_pow[32 + i] * SYZ;
            phi_tmp[192 + i] += 3.0 * xc_pow[32 + i] * yc_pow[i] * SXZ;
            phi_tmp[192 + i] += 9.0 * xc_pow[i] * yc_pow[i] * SZ;

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SXYZ;
            phi_tmp[224 + i] += 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SYZ;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SXZ;
            phi_tmp[224 + i] += xc_pow[32 + i] * yc_pow[i] * SXY;
            phi_tmp[224 + i] += 6.0 * xc_pow[i] * yc[i] * zc[i] * SZ;
            phi_tmp[224 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * yc[i] * SX;
            phi_tmp[224 + i] += 12.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SXYZ;
            phi_tmp[256 + i] += 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SYZ;
            phi_tmp[256 + i] += xc_pow[32 + i] * zc_pow[i] * SXZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SXY;
            phi_tmp[256 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[256 + i] += 6.0 * xc_pow[i] * yc[i] * zc[i] * SY;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * zc[i] * SX;
            phi_tmp[256 + i] += 3.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SXYZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[32 + i] * zc_pow[i] * SXY;
            phi_tmp[288 + i] += 9.0 * xc_pow[i] * zc_pow[i] * SY;

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SXYZ;
            phi_tmp[320 + i] += 2.0 * xc[i] * yc_pow[64 + i] * SYZ;
            phi_tmp[320 + i] += 4.0 * xc_pow[i] * yc_pow[32 + i] * SXZ;
            phi_tmp[320 + i] += 8.0 * xc[i] * yc_pow[32 + i] * SZ;

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SXYZ;
            phi_tmp[352 + i] += 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[352 + i] += 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SXZ;
            phi_tmp[352 + i] += xc_pow[i] * yc_pow[32 + i] * SXY;
            phi_tmp[352 + i] += 6.0 * xc[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[352 + i] += 2.0 * xc[i] * yc_pow[32 + i] * SY;
            phi_tmp[352 + i] += 3.0 * xc_pow[i] * yc_pow[i] * SX;
            phi_tmp[352 + i] += 18.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SXYZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc[i] * zc_pow[i] * SXZ;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc_pow[i] * zc[i] * SXY;
            phi_tmp[384 + i] += 4.0 * xc[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 4.0 * xc[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[384 + i] += 4.0 * xc_pow[i] * yc[i] * zc[i] * SX;
            phi_tmp[384 + i] += 8.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SXYZ;
            phi_tmp[416 + i] += 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[416 + i] += xc_pow[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[416 + i] += 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SXY;
            phi_tmp[416 + i] += 2.0 * xc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[416 + i] += 6.0 * xc[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[416 + i] += 3.0 * xc_pow[i] * zc_pow[i] * SX;
            phi_tmp[416 + i] += 2.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SXYZ;
            phi_tmp[448 + i] += 2.0 * xc[i] * zc_pow[64 + i] * SYZ;
            phi_tmp[448 + i] += 4.0 * xc_pow[i] * zc_pow[32 + i] * SXY;
            phi_tmp[448 + i] += 8.0 * xc[i] * zc_pow[32 + i] * SY;

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SXYZ;
            phi_tmp[480 + i] += yc_pow[96 + i] * SYZ;
            phi_tmp[480 + i] += 5.0 * xc[i] * yc_pow[64 + i] * SXZ;
            phi_tmp[480 + i] += 5.0 * yc_pow[64 + i] * SZ;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SXYZ;
            phi_tmp[512 + i] += yc_pow[64 + i] * zc[i] * SYZ;
            phi_tmp[512 + i] += 4.0 * xc[i] * yc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[512 + i] += xc[i] * yc_pow[64 + i] * SXY;
            phi_tmp[512 + i] += 4.0 * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[512 + i] += yc_pow[64 + i] * SY;
            phi_tmp[512 + i] += 4.0 * xc[i] * yc_pow[32 + i] * SX;
            phi_tmp[512 + i] += 16.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SXYZ;
            phi_tmp[544 + i] += yc_pow[32 + i] * zc_pow[i] * SYZ;
            phi_tmp[544 + i] += 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[544 + i] += 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SXY;
            phi_tmp[544 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[544 + i] += 6.0 * xc[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[544 + i] += 9.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SXYZ;
            phi_tmp[576 + i] += yc_pow[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[576 + i] += 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[576 + i] += 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SXY;
            phi_tmp[576 + i] += 2.0 * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[576 + i] += 3.0 * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[576 + i] += 6.0 * xc[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[576 + i] += 4.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SXYZ;
            phi_tmp[608 + i] += yc[i] * zc_pow[64 + i] * SYZ;
            phi_tmp[608 + i] += xc[i] * zc_pow[64 + i] * SXZ;
            phi_tmp[608 + i] += 4.0 * xc[i] * yc[i] * zc_pow[32 + i] * SXY;
            phi_tmp[608 + i] += zc_pow[64 + i] * SZ;
            phi_tmp[608 + i] += 4.0 * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[608 + i] += 4.0 * xc[i] * zc_pow[32 + i] * SX;
            phi_tmp[608 + i] += zc_pow[32 + i] * S0[i];

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SXYZ;
            phi_tmp[640 + i] += zc_pow[96 + i] * SYZ;
            phi_tmp[640 + i] += 5.0 * xc[i] * zc_pow[64 + i] * SXY;
            phi_tmp[640 + i] += 5.0 * zc_pow[64 + i] * SY;

            phi_tmp[672 + i] = yc_pow[128 + i] * SXYZ;
            phi_tmp[672 + i] += 6.0 * yc_pow[96 + i] * SXZ;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SXYZ;
            phi_tmp[704 + i] += 5.0 * yc_pow[64 + i] * zc[i] * SXZ;
            phi_tmp[704 + i] += yc_pow[96 + i] * SXY;
            phi_tmp[704 + i] += 5.0 * yc_pow[64 + i] * SX;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SXYZ;
            phi_tmp[736 + i] += 4.0 * yc_pow[32 + i] * zc_pow[i] * SXZ;
            phi_tmp[736 + i] += 2.0 * yc_pow[64 + i] * zc[i] * SXY;
            phi_tmp[736 + i] += 8.0 * yc_pow[32 + i] * zc[i] * SX;

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SXYZ;
            phi_tmp[768 + i] += 3.0 * yc_pow[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[768 + i] += 3.0 * yc_pow[32 + i] * zc_pow[i] * SXY;
            phi_tmp[768 + i] += 9.0 * yc_pow[i] * zc_pow[i] * SX;

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SXYZ;
            phi_tmp[800 + i] += 2.0 * yc[i] * zc_pow[64 + i] * SXZ;
            phi_tmp[800 + i] += 4.0 * yc_pow[i] * zc_pow[32 + i] * SXY;
            phi_tmp[800 + i] += 8.0 * yc[i] * zc_pow[32 + i] * SX;

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SXYZ;
            phi_tmp[832 + i] += zc_pow[96 + i] * SXZ;
            phi_tmp[832 + i] += 5.0 * yc[i] * zc_pow[64 + i] * SXY;
            phi_tmp[832 + i] += 5.0 * zc_pow[64 + i] * SX;

            phi_tmp[864 + i] = zc_pow[128 + i] * SXYZ;
            phi_tmp[864 + i] += 6.0 * zc_pow[96 + i] * SXY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_xyz_out + start), npoints);
        }

        // Combine XZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SX = S1[i] * xc[i];
            const double SZ = S1[i] * zc[i];
            const double SXZ = S2[i] * xc[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SXZZ = S3[i] * xc[i] * zc[i] * zc[i] + S2[i] * xc[i];

            phi_tmp[i] = xc_pow[128 + i] * SXZZ;
            phi_tmp[i] += 6.0 * xc_pow[96 + i] * SZZ;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SXZZ;
            phi_tmp[32 + i] += 5.0 * xc_pow[64 + i] * yc[i] * SZZ;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SXZZ;
            phi_tmp[64 + i] += 2.0 * xc_pow[96 + i] * SXZ;
            phi_tmp[64 + i] += 5.0 * xc_pow[64 + i] * zc[i] * SZZ;
            phi_tmp[64 + i] += 2.0 * 5.0 * xc_pow[64 + i] * SZ;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SXZZ;
            phi_tmp[96 + i] += 4.0 * xc_pow[32 + i] * yc_pow[i] * SZZ;

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SXZZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[64 + i] * yc[i] * SXZ;
            phi_tmp[128 + i] += 4.0 * xc_pow[32 + i] * yc[i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * 4.0 * xc_pow[32 + i] * yc[i] * SZ;

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SXZZ;
            phi_tmp[160 + i] += 2.0 * 2.0 * xc_pow[64 + i] * zc[i] * SXZ;
            phi_tmp[160 + i] += 4.0 * xc_pow[32 + i] * zc_pow[i] * SZZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[64 + i] * SX;
            phi_tmp[160 + i] += 2.0 * 8.0 * xc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[160 + i] += 8.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SXZZ;
            phi_tmp[192 + i] += 3.0 * xc_pow[i] * yc_pow[32 + i] * SZZ;

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SXZZ;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * yc_pow[i] * SXZ;
            phi_tmp[224 + i] += 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SZZ;
            phi_tmp[224 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[i] * SZ;

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SXZZ;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SXZ;
            phi_tmp[256 + i] += 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SZZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * yc[i] * SX;
            phi_tmp[256 + i] += 2.0 * 6.0 * xc_pow[i] * yc[i] * zc[i] * SZ;
            phi_tmp[256 + i] += 6.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SXZZ;
            phi_tmp[288 + i] += 2.0 * 3.0 * xc_pow[32 + i] * zc_pow[i] * SXZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[288 + i] += 6.0 * xc_pow[32 + i] * zc[i] * SX;
            phi_tmp[288 + i] += 2.0 * 9.0 * xc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[288 + i] += 18.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SXZZ;
            phi_tmp[320 + i] += 2.0 * xc[i] * yc_pow[64 + i] * SZZ;

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SXZZ;
            phi_tmp[352 + i] += 2.0 * xc_pow[i] * yc_pow[32 + i] * SXZ;
            phi_tmp[352 + i] += 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[352 + i] += 2.0 * 2.0 * xc[i] * yc_pow[32 + i] * SZ;

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SXZZ;
            phi_tmp[384 + i] += 2.0 * 2.0 * xc_pow[i] * yc_pow[i] * zc[i] * SXZ;
            phi_tmp[384 + i] += 2.0 * xc[i] * yc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc_pow[i] * SX;
            phi_tmp[384 + i] += 2.0 * 4.0 * xc[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[384 + i] += 4.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SXZZ;
            phi_tmp[416 + i] += 2.0 * 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SXZ;
            phi_tmp[416 + i] += 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[416 + i] += 6.0 * xc_pow[i] * yc[i] * zc[i] * SX;
            phi_tmp[416 + i] += 2.0 * 6.0 * xc[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[416 + i] += 12.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SXZZ;
            phi_tmp[448 + i] += 2.0 * 4.0 * xc_pow[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[448 + i] += 2.0 * xc[i] * zc_pow[64 + i] * SZZ;
            phi_tmp[448 + i] += 12.0 * xc_pow[i] * zc_pow[i] * SX;
            phi_tmp[448 + i] += 2.0 * 8.0 * xc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[448 + i] += 24.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SXZZ;
            phi_tmp[480 + i] += yc_pow[96 + i] * SZZ;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SXZZ;
            phi_tmp[512 + i] += 2.0 * xc[i] * yc_pow[64 + i] * SXZ;
            phi_tmp[512 + i] += yc_pow[64 + i] * zc[i] * SZZ;
            phi_tmp[512 + i] += 2.0 * yc_pow[64 + i] * SZ;

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SXZZ;
            phi_tmp[544 + i] += 2.0 * 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SXZ;
            phi_tmp[544 + i] += yc_pow[32 + i] * zc_pow[i] * SZZ;
            phi_tmp[544 + i] += 2.0 * xc[i] * yc_pow[32 + i] * SX;
            phi_tmp[544 + i] += 2.0 * 2.0 * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[544 + i] += 2.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SXZZ;
            phi_tmp[576 + i] += 2.0 * 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SXZ;
            phi_tmp[576 + i] += yc_pow[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[576 + i] += 6.0 * xc[i] * yc_pow[i] * zc[i] * SX;
            phi_tmp[576 + i] += 2.0 * 3.0 * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[576 + i] += 6.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SXZZ;
            phi_tmp[608 + i] += 2.0 * 4.0 * xc[i] * yc[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[608 + i] += yc[i] * zc_pow[64 + i] * SZZ;
            phi_tmp[608 + i] += 12.0 * xc[i] * yc[i] * zc_pow[i] * SX;
            phi_tmp[608 + i] += 2.0 * 4.0 * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[608 + i] += 12.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SXZZ;
            phi_tmp[640 + i] += 2.0 * 5.0 * xc[i] * zc_pow[64 + i] * SXZ;
            phi_tmp[640 + i] += zc_pow[96 + i] * SZZ;
            phi_tmp[640 + i] += 20.0 * xc[i] * zc_pow[32 + i] * SX;
            phi_tmp[640 + i] += 2.0 * 5.0 * zc_pow[64 + i] * SZ;
            phi_tmp[640 + i] += 20.0 * zc_pow[32 + i] * S0[i];

            phi_tmp[672 + i] = yc_pow[128 + i] * SXZZ;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SXZZ;
            phi_tmp[704 + i] += 2.0 * yc_pow[96 + i] * SXZ;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SXZZ;
            phi_tmp[736 + i] += 2.0 * 2.0 * yc_pow[64 + i] * zc[i] * SXZ;
            phi_tmp[736 + i] += 2.0 * yc_pow[64 + i] * SX;

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SXZZ;
            phi_tmp[768 + i] += 2.0 * 3.0 * yc_pow[32 + i] * zc_pow[i] * SXZ;
            phi_tmp[768 + i] += 6.0 * yc_pow[32 + i] * zc[i] * SX;

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SXZZ;
            phi_tmp[800 + i] += 2.0 * 4.0 * yc_pow[i] * zc_pow[32 + i] * SXZ;
            phi_tmp[800 + i] += 12.0 * yc_pow[i] * zc_pow[i] * SX;

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SXZZ;
            phi_tmp[832 + i] += 2.0 * 5.0 * yc[i] * zc_pow[64 + i] * SXZ;
            phi_tmp[832 + i] += 20.0 * yc[i] * zc_pow[32 + i] * SX;

            phi_tmp[864 + i] = zc_pow[128 + i] * SXZZ;
            phi_tmp[864 + i] += 2.0 * 6.0 * zc_pow[96 + i] * SXZ;
            phi_tmp[864 + i] += 30.0 * zc_pow[64 + i] * SX;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_xzz_out + start), npoints);
        }

        // Combine YYY blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SYYY = S3[i] * yc[i] * yc[i] * yc[i] + 3 * yc[i] * S2[i];

            phi_tmp[i] = xc_pow[128 + i] * SYYY;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SYYY;
            phi_tmp[32 + i] += 3.0 * xc_pow[96 + i] * SYY;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SYYY;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SYYY;
            phi_tmp[96 + i] += 3.0 * 2.0 * xc_pow[64 + i] * yc[i] * SYY;
            phi_tmp[96 + i] += 3.0 * 2.0 * xc_pow[64 + i] * SY;

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SYYY;
            phi_tmp[128 + i] += 3.0 * xc_pow[64 + i] * zc[i] * SYY;

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SYYY;

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SYYY;
            phi_tmp[192 + i] += 3.0 * 3.0 * xc_pow[32 + i] * yc_pow[i] * SYY;
            phi_tmp[192 + i] += 3.0 * 6.0 * xc_pow[32 + i] * yc[i] * SY;
            phi_tmp[192 + i] += 6.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SYYY;
            phi_tmp[224 + i] += 3.0 * 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SYY;
            phi_tmp[224 + i] += 3.0 * 2.0 * xc_pow[32 + i] * zc[i] * SY;

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SYYY;
            phi_tmp[256 + i] += 3.0 * xc_pow[32 + i] * zc_pow[i] * SYY;

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SYYY;

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SYYY;
            phi_tmp[320 + i] += 3.0 * 4.0 * xc_pow[i] * yc_pow[32 + i] * SYY;
            phi_tmp[320 + i] += 3.0 * 12.0 * xc_pow[i] * yc_pow[i] * SY;
            phi_tmp[320 + i] += 24.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SYYY;
            phi_tmp[352 + i] += 3.0 * 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SYY;
            phi_tmp[352 + i] += 3.0 * 6.0 * xc_pow[i] * yc[i] * zc[i] * SY;
            phi_tmp[352 + i] += 6.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SYYY;
            phi_tmp[384 + i] += 3.0 * 2.0 * xc_pow[i] * yc[i] * zc_pow[i] * SYY;
            phi_tmp[384 + i] += 3.0 * 2.0 * xc_pow[i] * zc_pow[i] * SY;

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SYYY;
            phi_tmp[416 + i] += 3.0 * xc_pow[i] * zc_pow[32 + i] * SYY;

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SYYY;

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SYYY;
            phi_tmp[480 + i] += 3.0 * 5.0 * xc[i] * yc_pow[64 + i] * SYY;
            phi_tmp[480 + i] += 3.0 * 20.0 * xc[i] * yc_pow[32 + i] * SY;
            phi_tmp[480 + i] += 60.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SYYY;
            phi_tmp[512 + i] += 3.0 * 4.0 * xc[i] * yc_pow[32 + i] * zc[i] * SYY;
            phi_tmp[512 + i] += 3.0 * 12.0 * xc[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[512 + i] += 24.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SYYY;
            phi_tmp[544 + i] += 3.0 * 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SYY;
            phi_tmp[544 + i] += 3.0 * 6.0 * xc[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[544 + i] += 6.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SYYY;
            phi_tmp[576 + i] += 3.0 * 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SYY;
            phi_tmp[576 + i] += 3.0 * 2.0 * xc[i] * zc_pow[32 + i] * SY;

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SYYY;
            phi_tmp[608 + i] += 3.0 * xc[i] * zc_pow[64 + i] * SYY;

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SYYY;

            phi_tmp[672 + i] = yc_pow[128 + i] * SYYY;
            phi_tmp[672 + i] += 3.0 * 6.0 * yc_pow[96 + i] * SYY;
            phi_tmp[672 + i] += 3.0 * 30.0 * yc_pow[64 + i] * SY;
            phi_tmp[672 + i] += 120.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SYYY;
            phi_tmp[704 + i] += 3.0 * 5.0 * yc_pow[64 + i] * zc[i] * SYY;
            phi_tmp[704 + i] += 3.0 * 20.0 * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[704 + i] += 60.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SYYY;
            phi_tmp[736 + i] += 3.0 * 4.0 * yc_pow[32 + i] * zc_pow[i] * SYY;
            phi_tmp[736 + i] += 3.0 * 12.0 * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[736 + i] += 24.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SYYY;
            phi_tmp[768 + i] += 3.0 * 3.0 * yc_pow[i] * zc_pow[32 + i] * SYY;
            phi_tmp[768 + i] += 3.0 * 6.0 * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[768 + i] += 6.0 * zc_pow[32 + i] * S0[i];

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SYYY;
            phi_tmp[800 + i] += 3.0 * 2.0 * yc[i] * zc_pow[64 + i] * SYY;
            phi_tmp[800 + i] += 3.0 * 2.0 * zc_pow[64 + i] * SY;

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SYYY;
            phi_tmp[832 + i] += 3.0 * zc_pow[96 + i] * SYY;

            phi_tmp[864 + i] = zc_pow[128 + i] * SYYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_yyy_out + start), npoints);
        }

        // Combine YYZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SYY = S2[i] * yc[i] * yc[i] + S1[i];
            const double SYYZ = S3[i] * yc[i] * yc[i] * zc[i] + S2[i] * zc[i];

            phi_tmp[i] = xc_pow[128 + i] * SYYZ;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SYYZ;
            phi_tmp[32 + i] += 2.0 * xc_pow[96 + i] * SYZ;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SYYZ;
            phi_tmp[64 + i] += xc_pow[96 + i] * SYY;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SYYZ;
            phi_tmp[96 + i] += 2.0 * 2.0 * xc_pow[64 + i] * yc[i] * SYZ;
            phi_tmp[96 + i] += 2.0 * xc_pow[64 + i] * SZ;

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SYYZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[64 + i] * zc[i] * SYZ;
            phi_tmp[128 + i] += xc_pow[64 + i] * yc[i] * SYY;
            phi_tmp[128 + i] += 2.0 * xc_pow[64 + i] * SY;

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SYYZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[64 + i] * zc[i] * SYY;

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SYYZ;
            phi_tmp[192 + i] += 2.0 * 3.0 * xc_pow[32 + i] * yc_pow[i] * SYZ;
            phi_tmp[192 + i] += 6.0 * xc_pow[32 + i] * yc[i] * SZ;

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SYYZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SYZ;
            phi_tmp[224 + i] += xc_pow[32 + i] * yc_pow[i] * SYY;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc_pow[32 + i] * yc[i] * SY;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SYYZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * zc_pow[i] * SYZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SYY;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc_pow[32 + i] * zc[i] * SY;

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SYYZ;
            phi_tmp[288 + i] += 3.0 * xc_pow[32 + i] * zc_pow[i] * SYY;

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SYYZ;
            phi_tmp[320 + i] += 2.0 * 4.0 * xc_pow[i] * yc_pow[32 + i] * SYZ;
            phi_tmp[320 + i] += 12.0 * xc_pow[i] * yc_pow[i] * SZ;

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SYYZ;
            phi_tmp[352 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SYZ;
            phi_tmp[352 + i] += xc_pow[i] * yc_pow[32 + i] * SYY;
            phi_tmp[352 + i] += 6.0 * xc_pow[i] * yc[i] * zc[i] * SZ;
            phi_tmp[352 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[i] * SY;
            phi_tmp[352 + i] += 6.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SYYZ;
            phi_tmp[384 + i] += 2.0 * 2.0 * xc_pow[i] * yc[i] * zc_pow[i] * SYZ;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc_pow[i] * zc[i] * SYY;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[384 + i] += 2.0 * 4.0 * xc_pow[i] * yc[i] * zc[i] * SY;
            phi_tmp[384 + i] += 4.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SYYZ;
            phi_tmp[416 + i] += 2.0 * xc_pow[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[416 + i] += 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SYY;
            phi_tmp[416 + i] += 2.0 * 3.0 * xc_pow[i] * zc_pow[i] * SY;

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SYYZ;
            phi_tmp[448 + i] += 4.0 * xc_pow[i] * zc_pow[32 + i] * SYY;

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SYYZ;
            phi_tmp[480 + i] += 2.0 * 5.0 * xc[i] * yc_pow[64 + i] * SYZ;
            phi_tmp[480 + i] += 20.0 * xc[i] * yc_pow[32 + i] * SZ;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SYYZ;
            phi_tmp[512 + i] += 2.0 * 4.0 * xc[i] * yc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[512 + i] += xc[i] * yc_pow[64 + i] * SYY;
            phi_tmp[512 + i] += 12.0 * xc[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[512 + i] += 2.0 * 4.0 * xc[i] * yc_pow[32 + i] * SY;
            phi_tmp[512 + i] += 12.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SYYZ;
            phi_tmp[544 + i] += 2.0 * 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[544 + i] += 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SYY;
            phi_tmp[544 + i] += 6.0 * xc[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[544 + i] += 2.0 * 6.0 * xc[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[544 + i] += 12.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SYYZ;
            phi_tmp[576 + i] += 2.0 * 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[576 + i] += 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SYY;
            phi_tmp[576 + i] += 2.0 * xc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[576 + i] += 2.0 * 6.0 * xc[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[576 + i] += 6.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SYYZ;
            phi_tmp[608 + i] += 2.0 * xc[i] * zc_pow[64 + i] * SYZ;
            phi_tmp[608 + i] += 4.0 * xc[i] * yc[i] * zc_pow[32 + i] * SYY;
            phi_tmp[608 + i] += 2.0 * 4.0 * xc[i] * zc_pow[32 + i] * SY;

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SYYZ;
            phi_tmp[640 + i] += 5.0 * xc[i] * zc_pow[64 + i] * SYY;

            phi_tmp[672 + i] = yc_pow[128 + i] * SYYZ;
            phi_tmp[672 + i] += 2.0 * 6.0 * yc_pow[96 + i] * SYZ;
            phi_tmp[672 + i] += 30.0 * yc_pow[64 + i] * SZ;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SYYZ;
            phi_tmp[704 + i] += 2.0 * 5.0 * yc_pow[64 + i] * zc[i] * SYZ;
            phi_tmp[704 + i] += yc_pow[96 + i] * SYY;
            phi_tmp[704 + i] += 20.0 * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[704 + i] += 2.0 * 5.0 * yc_pow[64 + i] * SY;
            phi_tmp[704 + i] += 20.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SYYZ;
            phi_tmp[736 + i] += 2.0 * 4.0 * yc_pow[32 + i] * zc_pow[i] * SYZ;
            phi_tmp[736 + i] += 2.0 * yc_pow[64 + i] * zc[i] * SYY;
            phi_tmp[736 + i] += 12.0 * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[736 + i] += 2.0 * 8.0 * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[736 + i] += 24.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SYYZ;
            phi_tmp[768 + i] += 2.0 * 3.0 * yc_pow[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[768 + i] += 3.0 * yc_pow[32 + i] * zc_pow[i] * SYY;
            phi_tmp[768 + i] += 6.0 * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[768 + i] += 2.0 * 9.0 * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[768 + i] += 18.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SYYZ;
            phi_tmp[800 + i] += 2.0 * 2.0 * yc[i] * zc_pow[64 + i] * SYZ;
            phi_tmp[800 + i] += 4.0 * yc_pow[i] * zc_pow[32 + i] * SYY;
            phi_tmp[800 + i] += 2.0 * zc_pow[64 + i] * SZ;
            phi_tmp[800 + i] += 2.0 * 8.0 * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[800 + i] += 8.0 * zc_pow[32 + i] * S0[i];

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SYYZ;
            phi_tmp[832 + i] += 2.0 * zc_pow[96 + i] * SYZ;
            phi_tmp[832 + i] += 5.0 * yc[i] * zc_pow[64 + i] * SYY;
            phi_tmp[832 + i] += 2.0 * 5.0 * zc_pow[64 + i] * SY;

            phi_tmp[864 + i] = zc_pow[128 + i] * SYYZ;
            phi_tmp[864 + i] += 6.0 * zc_pow[96 + i] * SYY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_yyz_out + start), npoints);
        }

        // Combine YZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            const double SYZ = S2[i] * yc[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SYZZ = S3[i] * yc[i] * zc[i] * zc[i] + S2[i] * yc[i];

            phi_tmp[i] = xc_pow[128 + i] * SYZZ;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SYZZ;
            phi_tmp[32 + i] += xc_pow[96 + i] * SZZ;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SYZZ;
            phi_tmp[64 + i] += 2.0 * xc_pow[96 + i] * SYZ;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SYZZ;
            phi_tmp[96 + i] += 2.0 * xc_pow[64 + i] * yc[i] * SZZ;

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SYZZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[64 + i] * yc[i] * SYZ;
            phi_tmp[128 + i] += xc_pow[64 + i] * zc[i] * SZZ;
            phi_tmp[128 + i] += 2.0 * xc_pow[64 + i] * SZ;

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SYZZ;
            phi_tmp[160 + i] += 2.0 * 2.0 * xc_pow[64 + i] * zc[i] * SYZ;
            phi_tmp[160 + i] += 2.0 * xc_pow[64 + i] * SY;

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SYZZ;
            phi_tmp[192 + i] += 3.0 * xc_pow[32 + i] * yc_pow[i] * SZZ;

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SYZZ;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * yc_pow[i] * SYZ;
            phi_tmp[224 + i] += 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SZZ;
            phi_tmp[224 + i] += 2.0 * 2.0 * xc_pow[32 + i] * yc[i] * SZ;

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SYZZ;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SYZ;
            phi_tmp[256 + i] += xc_pow[32 + i] * zc_pow[i] * SZZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * yc[i] * SY;
            phi_tmp[256 + i] += 2.0 * 2.0 * xc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[256 + i] += 2.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SYZZ;
            phi_tmp[288 + i] += 2.0 * 3.0 * xc_pow[32 + i] * zc_pow[i] * SYZ;
            phi_tmp[288 + i] += 6.0 * xc_pow[32 + i] * zc[i] * SY;

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SYZZ;
            phi_tmp[320 + i] += 4.0 * xc_pow[i] * yc_pow[32 + i] * SZZ;

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SYZZ;
            phi_tmp[352 + i] += 2.0 * xc_pow[i] * yc_pow[32 + i] * SYZ;
            phi_tmp[352 + i] += 3.0 * xc_pow[i] * yc_pow[i] * zc[i] * SZZ;
            phi_tmp[352 + i] += 2.0 * 3.0 * xc_pow[i] * yc_pow[i] * SZ;

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SYZZ;
            phi_tmp[384 + i] += 2.0 * 2.0 * xc_pow[i] * yc_pow[i] * zc[i] * SYZ;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc[i] * zc_pow[i] * SZZ;
            phi_tmp[384 + i] += 2.0 * xc_pow[i] * yc_pow[i] * SY;
            phi_tmp[384 + i] += 2.0 * 4.0 * xc_pow[i] * yc[i] * zc[i] * SZ;
            phi_tmp[384 + i] += 4.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SYZZ;
            phi_tmp[416 + i] += 2.0 * 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SYZ;
            phi_tmp[416 + i] += xc_pow[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[416 + i] += 6.0 * xc_pow[i] * yc[i] * zc[i] * SY;
            phi_tmp[416 + i] += 2.0 * 3.0 * xc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[416 + i] += 6.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SYZZ;
            phi_tmp[448 + i] += 2.0 * 4.0 * xc_pow[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[448 + i] += 12.0 * xc_pow[i] * zc_pow[i] * SY;

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SYZZ;
            phi_tmp[480 + i] += 5.0 * xc[i] * yc_pow[64 + i] * SZZ;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SYZZ;
            phi_tmp[512 + i] += 2.0 * xc[i] * yc_pow[64 + i] * SYZ;
            phi_tmp[512 + i] += 4.0 * xc[i] * yc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[512 + i] += 2.0 * 4.0 * xc[i] * yc_pow[32 + i] * SZ;

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SYZZ;
            phi_tmp[544 + i] += 2.0 * 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SYZ;
            phi_tmp[544 + i] += 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[544 + i] += 2.0 * xc[i] * yc_pow[32 + i] * SY;
            phi_tmp[544 + i] += 2.0 * 6.0 * xc[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[544 + i] += 6.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SYZZ;
            phi_tmp[576 + i] += 2.0 * 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SYZ;
            phi_tmp[576 + i] += 2.0 * xc[i] * yc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[576 + i] += 6.0 * xc[i] * yc_pow[i] * zc[i] * SY;
            phi_tmp[576 + i] += 2.0 * 6.0 * xc[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[576 + i] += 12.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SYZZ;
            phi_tmp[608 + i] += 2.0 * 4.0 * xc[i] * yc[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[608 + i] += xc[i] * zc_pow[64 + i] * SZZ;
            phi_tmp[608 + i] += 12.0 * xc[i] * yc[i] * zc_pow[i] * SY;
            phi_tmp[608 + i] += 2.0 * 4.0 * xc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[608 + i] += 12.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SYZZ;
            phi_tmp[640 + i] += 2.0 * 5.0 * xc[i] * zc_pow[64 + i] * SYZ;
            phi_tmp[640 + i] += 20.0 * xc[i] * zc_pow[32 + i] * SY;

            phi_tmp[672 + i] = yc_pow[128 + i] * SYZZ;
            phi_tmp[672 + i] += 6.0 * yc_pow[96 + i] * SZZ;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SYZZ;
            phi_tmp[704 + i] += 2.0 * yc_pow[96 + i] * SYZ;
            phi_tmp[704 + i] += 5.0 * yc_pow[64 + i] * zc[i] * SZZ;
            phi_tmp[704 + i] += 2.0 * 5.0 * yc_pow[64 + i] * SZ;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SYZZ;
            phi_tmp[736 + i] += 2.0 * 2.0 * yc_pow[64 + i] * zc[i] * SYZ;
            phi_tmp[736 + i] += 4.0 * yc_pow[32 + i] * zc_pow[i] * SZZ;
            phi_tmp[736 + i] += 2.0 * yc_pow[64 + i] * SY;
            phi_tmp[736 + i] += 2.0 * 8.0 * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[736 + i] += 8.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SYZZ;
            phi_tmp[768 + i] += 2.0 * 3.0 * yc_pow[32 + i] * zc_pow[i] * SYZ;
            phi_tmp[768 + i] += 3.0 * yc_pow[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[768 + i] += 6.0 * yc_pow[32 + i] * zc[i] * SY;
            phi_tmp[768 + i] += 2.0 * 9.0 * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[768 + i] += 18.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SYZZ;
            phi_tmp[800 + i] += 2.0 * 4.0 * yc_pow[i] * zc_pow[32 + i] * SYZ;
            phi_tmp[800 + i] += 2.0 * yc[i] * zc_pow[64 + i] * SZZ;
            phi_tmp[800 + i] += 12.0 * yc_pow[i] * zc_pow[i] * SY;
            phi_tmp[800 + i] += 2.0 * 8.0 * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[800 + i] += 24.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SYZZ;
            phi_tmp[832 + i] += 2.0 * 5.0 * yc[i] * zc_pow[64 + i] * SYZ;
            phi_tmp[832 + i] += zc_pow[96 + i] * SZZ;
            phi_tmp[832 + i] += 20.0 * yc[i] * zc_pow[32 + i] * SY;
            phi_tmp[832 + i] += 2.0 * 5.0 * zc_pow[64 + i] * SZ;
            phi_tmp[832 + i] += 20.0 * zc_pow[32 + i] * S0[i];

            phi_tmp[864 + i] = zc_pow[128 + i] * SYZZ;
            phi_tmp[864 + i] += 2.0 * 6.0 * zc_pow[96 + i] * SYZ;
            phi_tmp[864 + i] += 30.0 * zc_pow[64 + i] * SY;

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_yzz_out + start), npoints);
        }

        // Combine ZZZ blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            const double SZ = S1[i] * zc[i];
            const double SZZ = S2[i] * zc[i] * zc[i] + S1[i];
            const double SZZZ = S3[i] * zc[i] * zc[i] * zc[i] + 3 * zc[i] * S2[i];

            phi_tmp[i] = xc_pow[128 + i] * SZZZ;

            phi_tmp[32 + i] = xc_pow[96 + i] * yc[i] * SZZZ;

            phi_tmp[64 + i] = xc_pow[96 + i] * zc[i] * SZZZ;
            phi_tmp[64 + i] += 3.0 * xc_pow[96 + i] * SZZ;

            phi_tmp[96 + i] = xc_pow[64 + i] * yc_pow[i] * SZZZ;

            phi_tmp[128 + i] = xc_pow[64 + i] * yc[i] * zc[i] * SZZZ;
            phi_tmp[128 + i] += 3.0 * xc_pow[64 + i] * yc[i] * SZZ;

            phi_tmp[160 + i] = xc_pow[64 + i] * zc_pow[i] * SZZZ;
            phi_tmp[160 + i] += 3.0 * 2.0 * xc_pow[64 + i] * zc[i] * SZZ;
            phi_tmp[160 + i] += 3.0 * 2.0 * xc_pow[64 + i] * SZ;

            phi_tmp[192 + i] = xc_pow[32 + i] * yc_pow[32 + i] * SZZZ;

            phi_tmp[224 + i] = xc_pow[32 + i] * yc_pow[i] * zc[i] * SZZZ;
            phi_tmp[224 + i] += 3.0 * xc_pow[32 + i] * yc_pow[i] * SZZ;

            phi_tmp[256 + i] = xc_pow[32 + i] * yc[i] * zc_pow[i] * SZZZ;
            phi_tmp[256 + i] += 3.0 * 2.0 * xc_pow[32 + i] * yc[i] * zc[i] * SZZ;
            phi_tmp[256 + i] += 3.0 * 2.0 * xc_pow[32 + i] * yc[i] * SZ;

            phi_tmp[288 + i] = xc_pow[32 + i] * zc_pow[32 + i] * SZZZ;
            phi_tmp[288 + i] += 3.0 * 3.0 * xc_pow[32 + i] * zc_pow[i] * SZZ;
            phi_tmp[288 + i] += 3.0 * 6.0 * xc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[288 + i] += 6.0 * xc_pow[32 + i] * S0[i];

            phi_tmp[320 + i] = xc_pow[i] * yc_pow[64 + i] * SZZZ;

            phi_tmp[352 + i] = xc_pow[i] * yc_pow[32 + i] * zc[i] * SZZZ;
            phi_tmp[352 + i] += 3.0 * xc_pow[i] * yc_pow[32 + i] * SZZ;

            phi_tmp[384 + i] = xc_pow[i] * yc_pow[i] * zc_pow[i] * SZZZ;
            phi_tmp[384 + i] += 3.0 * 2.0 * xc_pow[i] * yc_pow[i] * zc[i] * SZZ;
            phi_tmp[384 + i] += 3.0 * 2.0 * xc_pow[i] * yc_pow[i] * SZ;

            phi_tmp[416 + i] = xc_pow[i] * yc[i] * zc_pow[32 + i] * SZZZ;
            phi_tmp[416 + i] += 3.0 * 3.0 * xc_pow[i] * yc[i] * zc_pow[i] * SZZ;
            phi_tmp[416 + i] += 3.0 * 6.0 * xc_pow[i] * yc[i] * zc[i] * SZ;
            phi_tmp[416 + i] += 6.0 * xc_pow[i] * yc[i] * S0[i];

            phi_tmp[448 + i] = xc_pow[i] * zc_pow[64 + i] * SZZZ;
            phi_tmp[448 + i] += 3.0 * 4.0 * xc_pow[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[448 + i] += 3.0 * 12.0 * xc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[448 + i] += 24.0 * xc_pow[i] * zc[i] * S0[i];

            phi_tmp[480 + i] = xc[i] * yc_pow[96 + i] * SZZZ;

            phi_tmp[512 + i] = xc[i] * yc_pow[64 + i] * zc[i] * SZZZ;
            phi_tmp[512 + i] += 3.0 * xc[i] * yc_pow[64 + i] * SZZ;

            phi_tmp[544 + i] = xc[i] * yc_pow[32 + i] * zc_pow[i] * SZZZ;
            phi_tmp[544 + i] += 3.0 * 2.0 * xc[i] * yc_pow[32 + i] * zc[i] * SZZ;
            phi_tmp[544 + i] += 3.0 * 2.0 * xc[i] * yc_pow[32 + i] * SZ;

            phi_tmp[576 + i] = xc[i] * yc_pow[i] * zc_pow[32 + i] * SZZZ;
            phi_tmp[576 + i] += 3.0 * 3.0 * xc[i] * yc_pow[i] * zc_pow[i] * SZZ;
            phi_tmp[576 + i] += 3.0 * 6.0 * xc[i] * yc_pow[i] * zc[i] * SZ;
            phi_tmp[576 + i] += 6.0 * xc[i] * yc_pow[i] * S0[i];

            phi_tmp[608 + i] = xc[i] * yc[i] * zc_pow[64 + i] * SZZZ;
            phi_tmp[608 + i] += 3.0 * 4.0 * xc[i] * yc[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[608 + i] += 3.0 * 12.0 * xc[i] * yc[i] * zc_pow[i] * SZ;
            phi_tmp[608 + i] += 24.0 * xc[i] * yc[i] * zc[i] * S0[i];

            phi_tmp[640 + i] = xc[i] * zc_pow[96 + i] * SZZZ;
            phi_tmp[640 + i] += 3.0 * 5.0 * xc[i] * zc_pow[64 + i] * SZZ;
            phi_tmp[640 + i] += 3.0 * 20.0 * xc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[640 + i] += 60.0 * xc[i] * zc_pow[i] * S0[i];

            phi_tmp[672 + i] = yc_pow[128 + i] * SZZZ;

            phi_tmp[704 + i] = yc_pow[96 + i] * zc[i] * SZZZ;
            phi_tmp[704 + i] += 3.0 * yc_pow[96 + i] * SZZ;

            phi_tmp[736 + i] = yc_pow[64 + i] * zc_pow[i] * SZZZ;
            phi_tmp[736 + i] += 3.0 * 2.0 * yc_pow[64 + i] * zc[i] * SZZ;
            phi_tmp[736 + i] += 3.0 * 2.0 * yc_pow[64 + i] * SZ;

            phi_tmp[768 + i] = yc_pow[32 + i] * zc_pow[32 + i] * SZZZ;
            phi_tmp[768 + i] += 3.0 * 3.0 * yc_pow[32 + i] * zc_pow[i] * SZZ;
            phi_tmp[768 + i] += 3.0 * 6.0 * yc_pow[32 + i] * zc[i] * SZ;
            phi_tmp[768 + i] += 6.0 * yc_pow[32 + i] * S0[i];

            phi_tmp[800 + i] = yc_pow[i] * zc_pow[64 + i] * SZZZ;
            phi_tmp[800 + i] += 3.0 * 4.0 * yc_pow[i] * zc_pow[32 + i] * SZZ;
            phi_tmp[800 + i] += 3.0 * 12.0 * yc_pow[i] * zc_pow[i] * SZ;
            phi_tmp[800 + i] += 24.0 * yc_pow[i] * zc[i] * S0[i];

            phi_tmp[832 + i] = yc[i] * zc_pow[96 + i] * SZZZ;
            phi_tmp[832 + i] += 3.0 * 5.0 * yc[i] * zc_pow[64 + i] * SZZ;
            phi_tmp[832 + i] += 3.0 * 20.0 * yc[i] * zc_pow[32 + i] * SZ;
            phi_tmp[832 + i] += 60.0 * yc[i] * zc_pow[i] * S0[i];

            phi_tmp[864 + i] = zc_pow[128 + i] * SZZZ;
            phi_tmp[864 + i] += 3.0 * 6.0 * zc_pow[96 + i] * SZZ;
            phi_tmp[864 + i] += 3.0 * 30.0 * zc_pow[64 + i] * SZ;
            phi_tmp[864 + i] += 120.0 * zc_pow[32 + i] * S0[i];

        }

        if (order == GG_SPHERICAL_CCA) {
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_zzz_out + start), npoints);
        }

    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);
    ALIGNED_FREE(expn2);

    // Free Power temporaries
    ALIGNED_FREE(xc_pow);
    ALIGNED_FREE(yc_pow);
    ALIGNED_FREE(zc_pow);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);

}
