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



void gg_collocation_L0_deriv1(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out) {

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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 224 * sizeof(double));
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

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;

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
            }

        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            // Gaussian derivs (gradients)
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];
            phi_out[start + i] = S0[i];

            // Gradient AM=0 Component=0
            phi_x_out[start + i] = SX;
            phi_y_out[start + i] = SY;
            phi_z_out[start + i] = SZ;
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

}

void gg_collocation_L1_deriv1(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out) {

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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 224 * sizeof(double));
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

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;

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
            }

        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            // Gaussian derivs (gradients)
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];

            // Density AM=1 Component=X
            phi_tmp[i] = S0[i] * xc[i];

            // Gradient AM=1 Component=X
            phi_x_tmp[i] = SX * xc[i];
            phi_y_tmp[i] = SY * xc[i];
            phi_z_tmp[i] = SZ * xc[i];
            phi_x_tmp[i] += S0[i];

            // Density AM=1 Component=Y
            phi_tmp[32 + i] = S0[i] * yc[i];

            // Gradient AM=1 Component=Y
            phi_x_tmp[32 + i] = SX * yc[i];
            phi_y_tmp[32 + i] = SY * yc[i];
            phi_z_tmp[32 + i] = SZ * yc[i];
            phi_y_tmp[32 + i] += S0[i];

            // Density AM=1 Component=Z
            phi_tmp[64 + i] = S0[i] * zc[i];

            // Gradient AM=1 Component=Z
            phi_x_tmp[64 + i] = SX * zc[i];
            phi_y_tmp[64 + i] = SY * zc[i];
            phi_z_tmp[64 + i] = SZ * zc[i];
            phi_z_tmp[64 + i] += S0[i];

        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L1(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_to_spherical_L1(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_to_spherical_L1(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L1(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_gaussian_cart_to_spherical_L1(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_gaussian_cart_to_spherical_L1(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L1(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_copy_L1(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_copy_L1(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L1(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_molden_cart_copy_L1(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_molden_cart_copy_L1(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
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

}

void gg_collocation_L2_deriv1(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out) {

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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 224 * sizeof(double));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));
    double* PRAGMA_RESTRICT expn2 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, 192 * sizeof(double));
    ASSUME_ALIGNED(phi_tmp, 64);
    double* PRAGMA_RESTRICT phi_x_tmp = (double*)ALIGNED_MALLOC(64, 192 * sizeof(double));
    ASSUME_ALIGNED(phi_x_tmp, 64);
    double* PRAGMA_RESTRICT phi_y_tmp = (double*)ALIGNED_MALLOC(64, 192 * sizeof(double));
    ASSUME_ALIGNED(phi_y_tmp, 64);
    double* PRAGMA_RESTRICT phi_z_tmp = (double*)ALIGNED_MALLOC(64, 192 * sizeof(double));
    ASSUME_ALIGNED(phi_z_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;

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
            }

        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            // Gaussian derivs (gradients)
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];

            // Cartesian derivs
            const double xc_pow2 = xc[i] * xc[i];
            const double yc_pow2 = yc[i] * yc[i];
            const double zc_pow2 = zc[i] * zc[i];


            // Density AM=2 Component=XX
            phi_tmp[i] = S0[i] * xc_pow2;

            // Gradient AM=2 Component=XX
            phi_x_tmp[i] = SX * xc_pow2;
            phi_y_tmp[i] = SY * xc_pow2;
            phi_z_tmp[i] = SZ * xc_pow2;
            AX = 2.0 * xc[i];
            phi_x_tmp[i] += S0[i] * AX;

            // Density AM=2 Component=XY
            A = xc[i] * yc[i];
            phi_tmp[32 + i] = S0[i] * A;

            // Gradient AM=2 Component=XY
            phi_x_tmp[32 + i] = SX * A;
            phi_y_tmp[32 + i] = SY * A;
            phi_z_tmp[32 + i] = SZ * A;
            phi_x_tmp[32 + i] += S0[i] * yc[i];
            phi_y_tmp[32 + i] += S0[i] * xc[i];

            // Density AM=2 Component=XZ
            A = xc[i] * zc[i];
            phi_tmp[64 + i] = S0[i] * A;

            // Gradient AM=2 Component=XZ
            phi_x_tmp[64 + i] = SX * A;
            phi_y_tmp[64 + i] = SY * A;
            phi_z_tmp[64 + i] = SZ * A;
            phi_x_tmp[64 + i] += S0[i] * zc[i];
            phi_z_tmp[64 + i] += S0[i] * xc[i];

            // Density AM=2 Component=YY
            phi_tmp[96 + i] = S0[i] * yc_pow2;

            // Gradient AM=2 Component=YY
            phi_x_tmp[96 + i] = SX * yc_pow2;
            phi_y_tmp[96 + i] = SY * yc_pow2;
            phi_z_tmp[96 + i] = SZ * yc_pow2;
            AY = 2.0 * yc[i];
            phi_y_tmp[96 + i] += S0[i] * AY;

            // Density AM=2 Component=YZ
            A = yc[i] * zc[i];
            phi_tmp[128 + i] = S0[i] * A;

            // Gradient AM=2 Component=YZ
            phi_x_tmp[128 + i] = SX * A;
            phi_y_tmp[128 + i] = SY * A;
            phi_z_tmp[128 + i] = SZ * A;
            phi_y_tmp[128 + i] += S0[i] * zc[i];
            phi_z_tmp[128 + i] += S0[i] * yc[i];

            // Density AM=2 Component=ZZ
            phi_tmp[160 + i] = S0[i] * zc_pow2;

            // Gradient AM=2 Component=ZZ
            phi_x_tmp[160 + i] = SX * zc_pow2;
            phi_y_tmp[160 + i] = SY * zc_pow2;
            phi_z_tmp[160 + i] = SZ * zc_pow2;
            AZ = 2.0 * zc[i];
            phi_z_tmp[160 + i] += S0[i] * AZ;

        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_to_spherical_L2(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_to_spherical_L2(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_to_spherical_L2(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_gaussian_cart_to_spherical_L2(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_gaussian_cart_to_spherical_L2(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_gaussian_cart_to_spherical_L2(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_copy_L2(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_copy_L2(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_copy_L2(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_molden_cart_copy_L2(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_molden_cart_copy_L2(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_molden_cart_copy_L2(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
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

}

void gg_collocation_L3_deriv1(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out) {

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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 224 * sizeof(double));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));
    double* PRAGMA_RESTRICT expn2 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, 320 * sizeof(double));
    ASSUME_ALIGNED(phi_tmp, 64);
    double* PRAGMA_RESTRICT phi_x_tmp = (double*)ALIGNED_MALLOC(64, 320 * sizeof(double));
    ASSUME_ALIGNED(phi_x_tmp, 64);
    double* PRAGMA_RESTRICT phi_y_tmp = (double*)ALIGNED_MALLOC(64, 320 * sizeof(double));
    ASSUME_ALIGNED(phi_y_tmp, 64);
    double* PRAGMA_RESTRICT phi_z_tmp = (double*)ALIGNED_MALLOC(64, 320 * sizeof(double));
    ASSUME_ALIGNED(phi_z_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;

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
            }

        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            // Gaussian derivs (gradients)
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];

            // Cartesian derivs
            const double xc_pow2 = xc[i] * xc[i];
            const double yc_pow2 = yc[i] * yc[i];
            const double zc_pow2 = zc[i] * zc[i];

            const double xc_pow3 = xc_pow2 * xc[i];
            const double yc_pow3 = yc_pow2 * yc[i];
            const double zc_pow3 = zc_pow2 * zc[i];


            // Density AM=3 Component=XXX
            phi_tmp[i] = S0[i] * xc_pow3;

            // Gradient AM=3 Component=XXX
            phi_x_tmp[i] = SX * xc_pow3;
            phi_y_tmp[i] = SY * xc_pow3;
            phi_z_tmp[i] = SZ * xc_pow3;
            AX = 3.0 * xc_pow2;
            phi_x_tmp[i] += S0[i] * AX;

            // Density AM=3 Component=XXY
            A = xc_pow2 * yc[i];
            phi_tmp[32 + i] = S0[i] * A;

            // Gradient AM=3 Component=XXY
            phi_x_tmp[32 + i] = SX * A;
            phi_y_tmp[32 + i] = SY * A;
            phi_z_tmp[32 + i] = SZ * A;
            AX = 2.0 * xc[i] * yc[i];
            phi_x_tmp[32 + i] += S0[i] * AX;
            phi_y_tmp[32 + i] += S0[i] * xc_pow2;

            // Density AM=3 Component=XXZ
            A = xc_pow2 * zc[i];
            phi_tmp[64 + i] = S0[i] * A;

            // Gradient AM=3 Component=XXZ
            phi_x_tmp[64 + i] = SX * A;
            phi_y_tmp[64 + i] = SY * A;
            phi_z_tmp[64 + i] = SZ * A;
            AX = 2.0 * xc[i] * zc[i];
            phi_x_tmp[64 + i] += S0[i] * AX;
            phi_z_tmp[64 + i] += S0[i] * xc_pow2;

            // Density AM=3 Component=XYY
            A = xc[i] * yc_pow2;
            phi_tmp[96 + i] = S0[i] * A;

            // Gradient AM=3 Component=XYY
            phi_x_tmp[96 + i] = SX * A;
            phi_y_tmp[96 + i] = SY * A;
            phi_z_tmp[96 + i] = SZ * A;
            phi_x_tmp[96 + i] += S0[i] * yc_pow2;
            AY = 2.0 * xc[i] * yc[i];
            phi_y_tmp[96 + i] += S0[i] * AY;

            // Density AM=3 Component=XYZ
            A = xc[i] * yc[i] * zc[i];
            phi_tmp[128 + i] = S0[i] * A;

            // Gradient AM=3 Component=XYZ
            phi_x_tmp[128 + i] = SX * A;
            phi_y_tmp[128 + i] = SY * A;
            phi_z_tmp[128 + i] = SZ * A;
            AX = yc[i] * zc[i];
            phi_x_tmp[128 + i] += S0[i] * AX;
            AY = xc[i] * zc[i];
            phi_y_tmp[128 + i] += S0[i] * AY;
            AZ = xc[i] * yc[i];
            phi_z_tmp[128 + i] += S0[i] * AZ;

            // Density AM=3 Component=XZZ
            A = xc[i] * zc_pow2;
            phi_tmp[160 + i] = S0[i] * A;

            // Gradient AM=3 Component=XZZ
            phi_x_tmp[160 + i] = SX * A;
            phi_y_tmp[160 + i] = SY * A;
            phi_z_tmp[160 + i] = SZ * A;
            phi_x_tmp[160 + i] += S0[i] * zc_pow2;
            AZ = 2.0 * xc[i] * zc[i];
            phi_z_tmp[160 + i] += S0[i] * AZ;

            // Density AM=3 Component=YYY
            phi_tmp[192 + i] = S0[i] * yc_pow3;

            // Gradient AM=3 Component=YYY
            phi_x_tmp[192 + i] = SX * yc_pow3;
            phi_y_tmp[192 + i] = SY * yc_pow3;
            phi_z_tmp[192 + i] = SZ * yc_pow3;
            AY = 3.0 * yc_pow2;
            phi_y_tmp[192 + i] += S0[i] * AY;

            // Density AM=3 Component=YYZ
            A = yc_pow2 * zc[i];
            phi_tmp[224 + i] = S0[i] * A;

            // Gradient AM=3 Component=YYZ
            phi_x_tmp[224 + i] = SX * A;
            phi_y_tmp[224 + i] = SY * A;
            phi_z_tmp[224 + i] = SZ * A;
            AY = 2.0 * yc[i] * zc[i];
            phi_y_tmp[224 + i] += S0[i] * AY;
            phi_z_tmp[224 + i] += S0[i] * yc_pow2;

            // Density AM=3 Component=YZZ
            A = yc[i] * zc_pow2;
            phi_tmp[256 + i] = S0[i] * A;

            // Gradient AM=3 Component=YZZ
            phi_x_tmp[256 + i] = SX * A;
            phi_y_tmp[256 + i] = SY * A;
            phi_z_tmp[256 + i] = SZ * A;
            phi_y_tmp[256 + i] += S0[i] * zc_pow2;
            AZ = 2.0 * yc[i] * zc[i];
            phi_z_tmp[256 + i] += S0[i] * AZ;

            // Density AM=3 Component=ZZZ
            phi_tmp[288 + i] = S0[i] * zc_pow3;

            // Gradient AM=3 Component=ZZZ
            phi_x_tmp[288 + i] = SX * zc_pow3;
            phi_y_tmp[288 + i] = SY * zc_pow3;
            phi_z_tmp[288 + i] = SZ * zc_pow3;
            AZ = 3.0 * zc_pow2;
            phi_z_tmp[288 + i] += S0[i] * AZ;

        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_to_spherical_L3(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_to_spherical_L3(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_to_spherical_L3(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_gaussian_cart_to_spherical_L3(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_gaussian_cart_to_spherical_L3(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_gaussian_cart_to_spherical_L3(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_copy_L3(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_copy_L3(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_copy_L3(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_molden_cart_copy_L3(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_molden_cart_copy_L3(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_molden_cart_copy_L3(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
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

}

void gg_collocation_L4_deriv1(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out) {

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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 224 * sizeof(double));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));
    double* PRAGMA_RESTRICT expn2 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, 480 * sizeof(double));
    ASSUME_ALIGNED(phi_tmp, 64);
    double* PRAGMA_RESTRICT phi_x_tmp = (double*)ALIGNED_MALLOC(64, 480 * sizeof(double));
    ASSUME_ALIGNED(phi_x_tmp, 64);
    double* PRAGMA_RESTRICT phi_y_tmp = (double*)ALIGNED_MALLOC(64, 480 * sizeof(double));
    ASSUME_ALIGNED(phi_y_tmp, 64);
    double* PRAGMA_RESTRICT phi_z_tmp = (double*)ALIGNED_MALLOC(64, 480 * sizeof(double));
    ASSUME_ALIGNED(phi_z_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;

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
            }

        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            // Gaussian derivs (gradients)
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];

            // Cartesian derivs
            const double xc_pow2 = xc[i] * xc[i];
            const double yc_pow2 = yc[i] * yc[i];
            const double zc_pow2 = zc[i] * zc[i];

            const double xc_pow3 = xc_pow2 * xc[i];
            const double yc_pow3 = yc_pow2 * yc[i];
            const double zc_pow3 = zc_pow2 * zc[i];

            const double xc_pow4 = xc_pow3 * xc[i];
            const double yc_pow4 = yc_pow3 * yc[i];
            const double zc_pow4 = zc_pow3 * zc[i];


            // Density AM=4 Component=XXXX
            phi_tmp[i] = S0[i] * xc_pow4;

            // Gradient AM=4 Component=XXXX
            phi_x_tmp[i] = SX * xc_pow4;
            phi_y_tmp[i] = SY * xc_pow4;
            phi_z_tmp[i] = SZ * xc_pow4;
            AX = 4.0 * xc_pow3;
            phi_x_tmp[i] += S0[i] * AX;

            // Density AM=4 Component=XXXY
            A = xc_pow3 * yc[i];
            phi_tmp[32 + i] = S0[i] * A;

            // Gradient AM=4 Component=XXXY
            phi_x_tmp[32 + i] = SX * A;
            phi_y_tmp[32 + i] = SY * A;
            phi_z_tmp[32 + i] = SZ * A;
            AX = 3.0 * xc_pow2 * yc[i];
            phi_x_tmp[32 + i] += S0[i] * AX;
            phi_y_tmp[32 + i] += S0[i] * xc_pow3;

            // Density AM=4 Component=XXXZ
            A = xc_pow3 * zc[i];
            phi_tmp[64 + i] = S0[i] * A;

            // Gradient AM=4 Component=XXXZ
            phi_x_tmp[64 + i] = SX * A;
            phi_y_tmp[64 + i] = SY * A;
            phi_z_tmp[64 + i] = SZ * A;
            AX = 3.0 * xc_pow2 * zc[i];
            phi_x_tmp[64 + i] += S0[i] * AX;
            phi_z_tmp[64 + i] += S0[i] * xc_pow3;

            // Density AM=4 Component=XXYY
            A = xc_pow2 * yc_pow2;
            phi_tmp[96 + i] = S0[i] * A;

            // Gradient AM=4 Component=XXYY
            phi_x_tmp[96 + i] = SX * A;
            phi_y_tmp[96 + i] = SY * A;
            phi_z_tmp[96 + i] = SZ * A;
            AX = 2.0 * xc[i] * yc_pow2;
            phi_x_tmp[96 + i] += S0[i] * AX;
            AY = 2.0 * xc_pow2 * yc[i];
            phi_y_tmp[96 + i] += S0[i] * AY;

            // Density AM=4 Component=XXYZ
            A = xc_pow2 * yc[i] * zc[i];
            phi_tmp[128 + i] = S0[i] * A;

            // Gradient AM=4 Component=XXYZ
            phi_x_tmp[128 + i] = SX * A;
            phi_y_tmp[128 + i] = SY * A;
            phi_z_tmp[128 + i] = SZ * A;
            AX = 2.0 * xc[i] * yc[i] * zc[i];
            phi_x_tmp[128 + i] += S0[i] * AX;
            AY = xc_pow2 * zc[i];
            phi_y_tmp[128 + i] += S0[i] * AY;
            AZ = xc_pow2 * yc[i];
            phi_z_tmp[128 + i] += S0[i] * AZ;

            // Density AM=4 Component=XXZZ
            A = xc_pow2 * zc_pow2;
            phi_tmp[160 + i] = S0[i] * A;

            // Gradient AM=4 Component=XXZZ
            phi_x_tmp[160 + i] = SX * A;
            phi_y_tmp[160 + i] = SY * A;
            phi_z_tmp[160 + i] = SZ * A;
            AX = 2.0 * xc[i] * zc_pow2;
            phi_x_tmp[160 + i] += S0[i] * AX;
            AZ = 2.0 * xc_pow2 * zc[i];
            phi_z_tmp[160 + i] += S0[i] * AZ;

            // Density AM=4 Component=XYYY
            A = xc[i] * yc_pow3;
            phi_tmp[192 + i] = S0[i] * A;

            // Gradient AM=4 Component=XYYY
            phi_x_tmp[192 + i] = SX * A;
            phi_y_tmp[192 + i] = SY * A;
            phi_z_tmp[192 + i] = SZ * A;
            phi_x_tmp[192 + i] += S0[i] * yc_pow3;
            AY = 3.0 * xc[i] * yc_pow2;
            phi_y_tmp[192 + i] += S0[i] * AY;

            // Density AM=4 Component=XYYZ
            A = xc[i] * yc_pow2 * zc[i];
            phi_tmp[224 + i] = S0[i] * A;

            // Gradient AM=4 Component=XYYZ
            phi_x_tmp[224 + i] = SX * A;
            phi_y_tmp[224 + i] = SY * A;
            phi_z_tmp[224 + i] = SZ * A;
            AX = yc_pow2 * zc[i];
            phi_x_tmp[224 + i] += S0[i] * AX;
            AY = 2.0 * xc[i] * yc[i] * zc[i];
            phi_y_tmp[224 + i] += S0[i] * AY;
            AZ = xc[i] * yc_pow2;
            phi_z_tmp[224 + i] += S0[i] * AZ;

            // Density AM=4 Component=XYZZ
            A = xc[i] * yc[i] * zc_pow2;
            phi_tmp[256 + i] = S0[i] * A;

            // Gradient AM=4 Component=XYZZ
            phi_x_tmp[256 + i] = SX * A;
            phi_y_tmp[256 + i] = SY * A;
            phi_z_tmp[256 + i] = SZ * A;
            AX = yc[i] * zc_pow2;
            phi_x_tmp[256 + i] += S0[i] * AX;
            AY = xc[i] * zc_pow2;
            phi_y_tmp[256 + i] += S0[i] * AY;
            AZ = 2.0 * xc[i] * yc[i] * zc[i];
            phi_z_tmp[256 + i] += S0[i] * AZ;

            // Density AM=4 Component=XZZZ
            A = xc[i] * zc_pow3;
            phi_tmp[288 + i] = S0[i] * A;

            // Gradient AM=4 Component=XZZZ
            phi_x_tmp[288 + i] = SX * A;
            phi_y_tmp[288 + i] = SY * A;
            phi_z_tmp[288 + i] = SZ * A;
            phi_x_tmp[288 + i] += S0[i] * zc_pow3;
            AZ = 3.0 * xc[i] * zc_pow2;
            phi_z_tmp[288 + i] += S0[i] * AZ;

            // Density AM=4 Component=YYYY
            phi_tmp[320 + i] = S0[i] * yc_pow4;

            // Gradient AM=4 Component=YYYY
            phi_x_tmp[320 + i] = SX * yc_pow4;
            phi_y_tmp[320 + i] = SY * yc_pow4;
            phi_z_tmp[320 + i] = SZ * yc_pow4;
            AY = 4.0 * yc_pow3;
            phi_y_tmp[320 + i] += S0[i] * AY;

            // Density AM=4 Component=YYYZ
            A = yc_pow3 * zc[i];
            phi_tmp[352 + i] = S0[i] * A;

            // Gradient AM=4 Component=YYYZ
            phi_x_tmp[352 + i] = SX * A;
            phi_y_tmp[352 + i] = SY * A;
            phi_z_tmp[352 + i] = SZ * A;
            AY = 3.0 * yc_pow2 * zc[i];
            phi_y_tmp[352 + i] += S0[i] * AY;
            phi_z_tmp[352 + i] += S0[i] * yc_pow3;

            // Density AM=4 Component=YYZZ
            A = yc_pow2 * zc_pow2;
            phi_tmp[384 + i] = S0[i] * A;

            // Gradient AM=4 Component=YYZZ
            phi_x_tmp[384 + i] = SX * A;
            phi_y_tmp[384 + i] = SY * A;
            phi_z_tmp[384 + i] = SZ * A;
            AY = 2.0 * yc[i] * zc_pow2;
            phi_y_tmp[384 + i] += S0[i] * AY;
            AZ = 2.0 * yc_pow2 * zc[i];
            phi_z_tmp[384 + i] += S0[i] * AZ;

            // Density AM=4 Component=YZZZ
            A = yc[i] * zc_pow3;
            phi_tmp[416 + i] = S0[i] * A;

            // Gradient AM=4 Component=YZZZ
            phi_x_tmp[416 + i] = SX * A;
            phi_y_tmp[416 + i] = SY * A;
            phi_z_tmp[416 + i] = SZ * A;
            phi_y_tmp[416 + i] += S0[i] * zc_pow3;
            AZ = 3.0 * yc[i] * zc_pow2;
            phi_z_tmp[416 + i] += S0[i] * AZ;

            // Density AM=4 Component=ZZZZ
            phi_tmp[448 + i] = S0[i] * zc_pow4;

            // Gradient AM=4 Component=ZZZZ
            phi_x_tmp[448 + i] = SX * zc_pow4;
            phi_y_tmp[448 + i] = SY * zc_pow4;
            phi_z_tmp[448 + i] = SZ * zc_pow4;
            AZ = 4.0 * zc_pow3;
            phi_z_tmp[448 + i] += S0[i] * AZ;

        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_to_spherical_L4(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_to_spherical_L4(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_to_spherical_L4(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_gaussian_cart_to_spherical_L4(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_gaussian_cart_to_spherical_L4(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_gaussian_cart_to_spherical_L4(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_copy_L4(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_copy_L4(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_copy_L4(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_molden_cart_copy_L4(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_molden_cart_copy_L4(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_molden_cart_copy_L4(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
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

}

void gg_collocation_L5_deriv1(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out) {

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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 224 * sizeof(double));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));
    double* PRAGMA_RESTRICT expn2 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, 672 * sizeof(double));
    ASSUME_ALIGNED(phi_tmp, 64);
    double* PRAGMA_RESTRICT phi_x_tmp = (double*)ALIGNED_MALLOC(64, 672 * sizeof(double));
    ASSUME_ALIGNED(phi_x_tmp, 64);
    double* PRAGMA_RESTRICT phi_y_tmp = (double*)ALIGNED_MALLOC(64, 672 * sizeof(double));
    ASSUME_ALIGNED(phi_y_tmp, 64);
    double* PRAGMA_RESTRICT phi_z_tmp = (double*)ALIGNED_MALLOC(64, 672 * sizeof(double));
    ASSUME_ALIGNED(phi_z_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;

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
            }

        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            // Gaussian derivs (gradients)
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];

            // Cartesian derivs
            const double xc_pow2 = xc[i] * xc[i];
            const double yc_pow2 = yc[i] * yc[i];
            const double zc_pow2 = zc[i] * zc[i];

            const double xc_pow3 = xc_pow2 * xc[i];
            const double yc_pow3 = yc_pow2 * yc[i];
            const double zc_pow3 = zc_pow2 * zc[i];

            const double xc_pow4 = xc_pow3 * xc[i];
            const double yc_pow4 = yc_pow3 * yc[i];
            const double zc_pow4 = zc_pow3 * zc[i];

            const double xc_pow5 = xc_pow4 * xc[i];
            const double yc_pow5 = yc_pow4 * yc[i];
            const double zc_pow5 = zc_pow4 * zc[i];


            // Density AM=5 Component=XXXXX
            phi_tmp[i] = S0[i] * xc_pow5;

            // Gradient AM=5 Component=XXXXX
            phi_x_tmp[i] = SX * xc_pow5;
            phi_y_tmp[i] = SY * xc_pow5;
            phi_z_tmp[i] = SZ * xc_pow5;
            AX = 5.0 * xc_pow4;
            phi_x_tmp[i] += S0[i] * AX;

            // Density AM=5 Component=XXXXY
            A = xc_pow4 * yc[i];
            phi_tmp[32 + i] = S0[i] * A;

            // Gradient AM=5 Component=XXXXY
            phi_x_tmp[32 + i] = SX * A;
            phi_y_tmp[32 + i] = SY * A;
            phi_z_tmp[32 + i] = SZ * A;
            AX = 4.0 * xc_pow3 * yc[i];
            phi_x_tmp[32 + i] += S0[i] * AX;
            phi_y_tmp[32 + i] += S0[i] * xc_pow4;

            // Density AM=5 Component=XXXXZ
            A = xc_pow4 * zc[i];
            phi_tmp[64 + i] = S0[i] * A;

            // Gradient AM=5 Component=XXXXZ
            phi_x_tmp[64 + i] = SX * A;
            phi_y_tmp[64 + i] = SY * A;
            phi_z_tmp[64 + i] = SZ * A;
            AX = 4.0 * xc_pow3 * zc[i];
            phi_x_tmp[64 + i] += S0[i] * AX;
            phi_z_tmp[64 + i] += S0[i] * xc_pow4;

            // Density AM=5 Component=XXXYY
            A = xc_pow3 * yc_pow2;
            phi_tmp[96 + i] = S0[i] * A;

            // Gradient AM=5 Component=XXXYY
            phi_x_tmp[96 + i] = SX * A;
            phi_y_tmp[96 + i] = SY * A;
            phi_z_tmp[96 + i] = SZ * A;
            AX = 3.0 * xc_pow2 * yc_pow2;
            phi_x_tmp[96 + i] += S0[i] * AX;
            AY = 2.0 * xc_pow3 * yc[i];
            phi_y_tmp[96 + i] += S0[i] * AY;

            // Density AM=5 Component=XXXYZ
            A = xc_pow3 * yc[i] * zc[i];
            phi_tmp[128 + i] = S0[i] * A;

            // Gradient AM=5 Component=XXXYZ
            phi_x_tmp[128 + i] = SX * A;
            phi_y_tmp[128 + i] = SY * A;
            phi_z_tmp[128 + i] = SZ * A;
            AX = 3.0 * xc_pow2 * yc[i] * zc[i];
            phi_x_tmp[128 + i] += S0[i] * AX;
            AY = xc_pow3 * zc[i];
            phi_y_tmp[128 + i] += S0[i] * AY;
            AZ = xc_pow3 * yc[i];
            phi_z_tmp[128 + i] += S0[i] * AZ;

            // Density AM=5 Component=XXXZZ
            A = xc_pow3 * zc_pow2;
            phi_tmp[160 + i] = S0[i] * A;

            // Gradient AM=5 Component=XXXZZ
            phi_x_tmp[160 + i] = SX * A;
            phi_y_tmp[160 + i] = SY * A;
            phi_z_tmp[160 + i] = SZ * A;
            AX = 3.0 * xc_pow2 * zc_pow2;
            phi_x_tmp[160 + i] += S0[i] * AX;
            AZ = 2.0 * xc_pow3 * zc[i];
            phi_z_tmp[160 + i] += S0[i] * AZ;

            // Density AM=5 Component=XXYYY
            A = xc_pow2 * yc_pow3;
            phi_tmp[192 + i] = S0[i] * A;

            // Gradient AM=5 Component=XXYYY
            phi_x_tmp[192 + i] = SX * A;
            phi_y_tmp[192 + i] = SY * A;
            phi_z_tmp[192 + i] = SZ * A;
            AX = 2.0 * xc[i] * yc_pow3;
            phi_x_tmp[192 + i] += S0[i] * AX;
            AY = 3.0 * xc_pow2 * yc_pow2;
            phi_y_tmp[192 + i] += S0[i] * AY;

            // Density AM=5 Component=XXYYZ
            A = xc_pow2 * yc_pow2 * zc[i];
            phi_tmp[224 + i] = S0[i] * A;

            // Gradient AM=5 Component=XXYYZ
            phi_x_tmp[224 + i] = SX * A;
            phi_y_tmp[224 + i] = SY * A;
            phi_z_tmp[224 + i] = SZ * A;
            AX = 2.0 * xc[i] * yc_pow2 * zc[i];
            phi_x_tmp[224 + i] += S0[i] * AX;
            AY = 2.0 * xc_pow2 * yc[i] * zc[i];
            phi_y_tmp[224 + i] += S0[i] * AY;
            AZ = xc_pow2 * yc_pow2;
            phi_z_tmp[224 + i] += S0[i] * AZ;

            // Density AM=5 Component=XXYZZ
            A = xc_pow2 * yc[i] * zc_pow2;
            phi_tmp[256 + i] = S0[i] * A;

            // Gradient AM=5 Component=XXYZZ
            phi_x_tmp[256 + i] = SX * A;
            phi_y_tmp[256 + i] = SY * A;
            phi_z_tmp[256 + i] = SZ * A;
            AX = 2.0 * xc[i] * yc[i] * zc_pow2;
            phi_x_tmp[256 + i] += S0[i] * AX;
            AY = xc_pow2 * zc_pow2;
            phi_y_tmp[256 + i] += S0[i] * AY;
            AZ = 2.0 * xc_pow2 * yc[i] * zc[i];
            phi_z_tmp[256 + i] += S0[i] * AZ;

            // Density AM=5 Component=XXZZZ
            A = xc_pow2 * zc_pow3;
            phi_tmp[288 + i] = S0[i] * A;

            // Gradient AM=5 Component=XXZZZ
            phi_x_tmp[288 + i] = SX * A;
            phi_y_tmp[288 + i] = SY * A;
            phi_z_tmp[288 + i] = SZ * A;
            AX = 2.0 * xc[i] * zc_pow3;
            phi_x_tmp[288 + i] += S0[i] * AX;
            AZ = 3.0 * xc_pow2 * zc_pow2;
            phi_z_tmp[288 + i] += S0[i] * AZ;

            // Density AM=5 Component=XYYYY
            A = xc[i] * yc_pow4;
            phi_tmp[320 + i] = S0[i] * A;

            // Gradient AM=5 Component=XYYYY
            phi_x_tmp[320 + i] = SX * A;
            phi_y_tmp[320 + i] = SY * A;
            phi_z_tmp[320 + i] = SZ * A;
            phi_x_tmp[320 + i] += S0[i] * yc_pow4;
            AY = 4.0 * xc[i] * yc_pow3;
            phi_y_tmp[320 + i] += S0[i] * AY;

            // Density AM=5 Component=XYYYZ
            A = xc[i] * yc_pow3 * zc[i];
            phi_tmp[352 + i] = S0[i] * A;

            // Gradient AM=5 Component=XYYYZ
            phi_x_tmp[352 + i] = SX * A;
            phi_y_tmp[352 + i] = SY * A;
            phi_z_tmp[352 + i] = SZ * A;
            AX = yc_pow3 * zc[i];
            phi_x_tmp[352 + i] += S0[i] * AX;
            AY = 3.0 * xc[i] * yc_pow2 * zc[i];
            phi_y_tmp[352 + i] += S0[i] * AY;
            AZ = xc[i] * yc_pow3;
            phi_z_tmp[352 + i] += S0[i] * AZ;

            // Density AM=5 Component=XYYZZ
            A = xc[i] * yc_pow2 * zc_pow2;
            phi_tmp[384 + i] = S0[i] * A;

            // Gradient AM=5 Component=XYYZZ
            phi_x_tmp[384 + i] = SX * A;
            phi_y_tmp[384 + i] = SY * A;
            phi_z_tmp[384 + i] = SZ * A;
            AX = yc_pow2 * zc_pow2;
            phi_x_tmp[384 + i] += S0[i] * AX;
            AY = 2.0 * xc[i] * yc[i] * zc_pow2;
            phi_y_tmp[384 + i] += S0[i] * AY;
            AZ = 2.0 * xc[i] * yc_pow2 * zc[i];
            phi_z_tmp[384 + i] += S0[i] * AZ;

            // Density AM=5 Component=XYZZZ
            A = xc[i] * yc[i] * zc_pow3;
            phi_tmp[416 + i] = S0[i] * A;

            // Gradient AM=5 Component=XYZZZ
            phi_x_tmp[416 + i] = SX * A;
            phi_y_tmp[416 + i] = SY * A;
            phi_z_tmp[416 + i] = SZ * A;
            AX = yc[i] * zc_pow3;
            phi_x_tmp[416 + i] += S0[i] * AX;
            AY = xc[i] * zc_pow3;
            phi_y_tmp[416 + i] += S0[i] * AY;
            AZ = 3.0 * xc[i] * yc[i] * zc_pow2;
            phi_z_tmp[416 + i] += S0[i] * AZ;

            // Density AM=5 Component=XZZZZ
            A = xc[i] * zc_pow4;
            phi_tmp[448 + i] = S0[i] * A;

            // Gradient AM=5 Component=XZZZZ
            phi_x_tmp[448 + i] = SX * A;
            phi_y_tmp[448 + i] = SY * A;
            phi_z_tmp[448 + i] = SZ * A;
            phi_x_tmp[448 + i] += S0[i] * zc_pow4;
            AZ = 4.0 * xc[i] * zc_pow3;
            phi_z_tmp[448 + i] += S0[i] * AZ;

            // Density AM=5 Component=YYYYY
            phi_tmp[480 + i] = S0[i] * yc_pow5;

            // Gradient AM=5 Component=YYYYY
            phi_x_tmp[480 + i] = SX * yc_pow5;
            phi_y_tmp[480 + i] = SY * yc_pow5;
            phi_z_tmp[480 + i] = SZ * yc_pow5;
            AY = 5.0 * yc_pow4;
            phi_y_tmp[480 + i] += S0[i] * AY;

            // Density AM=5 Component=YYYYZ
            A = yc_pow4 * zc[i];
            phi_tmp[512 + i] = S0[i] * A;

            // Gradient AM=5 Component=YYYYZ
            phi_x_tmp[512 + i] = SX * A;
            phi_y_tmp[512 + i] = SY * A;
            phi_z_tmp[512 + i] = SZ * A;
            AY = 4.0 * yc_pow3 * zc[i];
            phi_y_tmp[512 + i] += S0[i] * AY;
            phi_z_tmp[512 + i] += S0[i] * yc_pow4;

            // Density AM=5 Component=YYYZZ
            A = yc_pow3 * zc_pow2;
            phi_tmp[544 + i] = S0[i] * A;

            // Gradient AM=5 Component=YYYZZ
            phi_x_tmp[544 + i] = SX * A;
            phi_y_tmp[544 + i] = SY * A;
            phi_z_tmp[544 + i] = SZ * A;
            AY = 3.0 * yc_pow2 * zc_pow2;
            phi_y_tmp[544 + i] += S0[i] * AY;
            AZ = 2.0 * yc_pow3 * zc[i];
            phi_z_tmp[544 + i] += S0[i] * AZ;

            // Density AM=5 Component=YYZZZ
            A = yc_pow2 * zc_pow3;
            phi_tmp[576 + i] = S0[i] * A;

            // Gradient AM=5 Component=YYZZZ
            phi_x_tmp[576 + i] = SX * A;
            phi_y_tmp[576 + i] = SY * A;
            phi_z_tmp[576 + i] = SZ * A;
            AY = 2.0 * yc[i] * zc_pow3;
            phi_y_tmp[576 + i] += S0[i] * AY;
            AZ = 3.0 * yc_pow2 * zc_pow2;
            phi_z_tmp[576 + i] += S0[i] * AZ;

            // Density AM=5 Component=YZZZZ
            A = yc[i] * zc_pow4;
            phi_tmp[608 + i] = S0[i] * A;

            // Gradient AM=5 Component=YZZZZ
            phi_x_tmp[608 + i] = SX * A;
            phi_y_tmp[608 + i] = SY * A;
            phi_z_tmp[608 + i] = SZ * A;
            phi_y_tmp[608 + i] += S0[i] * zc_pow4;
            AZ = 4.0 * yc[i] * zc_pow3;
            phi_z_tmp[608 + i] += S0[i] * AZ;

            // Density AM=5 Component=ZZZZZ
            phi_tmp[640 + i] = S0[i] * zc_pow5;

            // Gradient AM=5 Component=ZZZZZ
            phi_x_tmp[640 + i] = SX * zc_pow5;
            phi_y_tmp[640 + i] = SY * zc_pow5;
            phi_z_tmp[640 + i] = SZ * zc_pow5;
            AZ = 5.0 * zc_pow4;
            phi_z_tmp[640 + i] += S0[i] * AZ;

        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_to_spherical_L5(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_to_spherical_L5(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_to_spherical_L5(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_gaussian_cart_to_spherical_L5(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_gaussian_cart_to_spherical_L5(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_gaussian_cart_to_spherical_L5(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_copy_L5(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_copy_L5(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_copy_L5(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_molden_cart_copy_L5(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_molden_cart_copy_L5(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_molden_cart_copy_L5(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
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

}

void gg_collocation_L6_deriv1(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out, double* PRAGMA_RESTRICT phi_x_out, double* PRAGMA_RESTRICT phi_y_out, double* PRAGMA_RESTRICT phi_z_out) {

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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, 224 * sizeof(double));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));
    double* PRAGMA_RESTRICT expn2 = (double*)ALIGNED_MALLOC(64, nprim * sizeof(double));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, 896 * sizeof(double));
    ASSUME_ALIGNED(phi_tmp, 64);
    double* PRAGMA_RESTRICT phi_x_tmp = (double*)ALIGNED_MALLOC(64, 896 * sizeof(double));
    ASSUME_ALIGNED(phi_x_tmp, 64);
    double* PRAGMA_RESTRICT phi_y_tmp = (double*)ALIGNED_MALLOC(64, 896 * sizeof(double));
    ASSUME_ALIGNED(phi_y_tmp, 64);
    double* PRAGMA_RESTRICT phi_z_tmp = (double*)ALIGNED_MALLOC(64, 896 * sizeof(double));
    ASSUME_ALIGNED(phi_z_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;
    double AX, AY, AZ;

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
            }

        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            // Gaussian derivs (gradients)
            const double SX = S1[i] * xc[i];
            const double SY = S1[i] * yc[i];
            const double SZ = S1[i] * zc[i];

            // Cartesian derivs
            const double xc_pow2 = xc[i] * xc[i];
            const double yc_pow2 = yc[i] * yc[i];
            const double zc_pow2 = zc[i] * zc[i];

            const double xc_pow3 = xc_pow2 * xc[i];
            const double yc_pow3 = yc_pow2 * yc[i];
            const double zc_pow3 = zc_pow2 * zc[i];

            const double xc_pow4 = xc_pow3 * xc[i];
            const double yc_pow4 = yc_pow3 * yc[i];
            const double zc_pow4 = zc_pow3 * zc[i];

            const double xc_pow5 = xc_pow4 * xc[i];
            const double yc_pow5 = yc_pow4 * yc[i];
            const double zc_pow5 = zc_pow4 * zc[i];

            const double xc_pow6 = xc_pow5 * xc[i];
            const double yc_pow6 = yc_pow5 * yc[i];
            const double zc_pow6 = zc_pow5 * zc[i];


            // Density AM=6 Component=XXXXXX
            phi_tmp[i] = S0[i] * xc_pow6;

            // Gradient AM=6 Component=XXXXXX
            phi_x_tmp[i] = SX * xc_pow6;
            phi_y_tmp[i] = SY * xc_pow6;
            phi_z_tmp[i] = SZ * xc_pow6;
            AX = 6.0 * xc_pow5;
            phi_x_tmp[i] += S0[i] * AX;

            // Density AM=6 Component=XXXXXY
            A = xc_pow5 * yc[i];
            phi_tmp[32 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXXXXY
            phi_x_tmp[32 + i] = SX * A;
            phi_y_tmp[32 + i] = SY * A;
            phi_z_tmp[32 + i] = SZ * A;
            AX = 5.0 * xc_pow4 * yc[i];
            phi_x_tmp[32 + i] += S0[i] * AX;
            phi_y_tmp[32 + i] += S0[i] * xc_pow5;

            // Density AM=6 Component=XXXXXZ
            A = xc_pow5 * zc[i];
            phi_tmp[64 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXXXXZ
            phi_x_tmp[64 + i] = SX * A;
            phi_y_tmp[64 + i] = SY * A;
            phi_z_tmp[64 + i] = SZ * A;
            AX = 5.0 * xc_pow4 * zc[i];
            phi_x_tmp[64 + i] += S0[i] * AX;
            phi_z_tmp[64 + i] += S0[i] * xc_pow5;

            // Density AM=6 Component=XXXXYY
            A = xc_pow4 * yc_pow2;
            phi_tmp[96 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXXXYY
            phi_x_tmp[96 + i] = SX * A;
            phi_y_tmp[96 + i] = SY * A;
            phi_z_tmp[96 + i] = SZ * A;
            AX = 4.0 * xc_pow3 * yc_pow2;
            phi_x_tmp[96 + i] += S0[i] * AX;
            AY = 2.0 * xc_pow4 * yc[i];
            phi_y_tmp[96 + i] += S0[i] * AY;

            // Density AM=6 Component=XXXXYZ
            A = xc_pow4 * yc[i] * zc[i];
            phi_tmp[128 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXXXYZ
            phi_x_tmp[128 + i] = SX * A;
            phi_y_tmp[128 + i] = SY * A;
            phi_z_tmp[128 + i] = SZ * A;
            AX = 4.0 * xc_pow3 * yc[i] * zc[i];
            phi_x_tmp[128 + i] += S0[i] * AX;
            AY = xc_pow4 * zc[i];
            phi_y_tmp[128 + i] += S0[i] * AY;
            AZ = xc_pow4 * yc[i];
            phi_z_tmp[128 + i] += S0[i] * AZ;

            // Density AM=6 Component=XXXXZZ
            A = xc_pow4 * zc_pow2;
            phi_tmp[160 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXXXZZ
            phi_x_tmp[160 + i] = SX * A;
            phi_y_tmp[160 + i] = SY * A;
            phi_z_tmp[160 + i] = SZ * A;
            AX = 4.0 * xc_pow3 * zc_pow2;
            phi_x_tmp[160 + i] += S0[i] * AX;
            AZ = 2.0 * xc_pow4 * zc[i];
            phi_z_tmp[160 + i] += S0[i] * AZ;

            // Density AM=6 Component=XXXYYY
            A = xc_pow3 * yc_pow3;
            phi_tmp[192 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXXYYY
            phi_x_tmp[192 + i] = SX * A;
            phi_y_tmp[192 + i] = SY * A;
            phi_z_tmp[192 + i] = SZ * A;
            AX = 3.0 * xc_pow2 * yc_pow3;
            phi_x_tmp[192 + i] += S0[i] * AX;
            AY = 3.0 * xc_pow3 * yc_pow2;
            phi_y_tmp[192 + i] += S0[i] * AY;

            // Density AM=6 Component=XXXYYZ
            A = xc_pow3 * yc_pow2 * zc[i];
            phi_tmp[224 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXXYYZ
            phi_x_tmp[224 + i] = SX * A;
            phi_y_tmp[224 + i] = SY * A;
            phi_z_tmp[224 + i] = SZ * A;
            AX = 3.0 * xc_pow2 * yc_pow2 * zc[i];
            phi_x_tmp[224 + i] += S0[i] * AX;
            AY = 2.0 * xc_pow3 * yc[i] * zc[i];
            phi_y_tmp[224 + i] += S0[i] * AY;
            AZ = xc_pow3 * yc_pow2;
            phi_z_tmp[224 + i] += S0[i] * AZ;

            // Density AM=6 Component=XXXYZZ
            A = xc_pow3 * yc[i] * zc_pow2;
            phi_tmp[256 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXXYZZ
            phi_x_tmp[256 + i] = SX * A;
            phi_y_tmp[256 + i] = SY * A;
            phi_z_tmp[256 + i] = SZ * A;
            AX = 3.0 * xc_pow2 * yc[i] * zc_pow2;
            phi_x_tmp[256 + i] += S0[i] * AX;
            AY = xc_pow3 * zc_pow2;
            phi_y_tmp[256 + i] += S0[i] * AY;
            AZ = 2.0 * xc_pow3 * yc[i] * zc[i];
            phi_z_tmp[256 + i] += S0[i] * AZ;

            // Density AM=6 Component=XXXZZZ
            A = xc_pow3 * zc_pow3;
            phi_tmp[288 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXXZZZ
            phi_x_tmp[288 + i] = SX * A;
            phi_y_tmp[288 + i] = SY * A;
            phi_z_tmp[288 + i] = SZ * A;
            AX = 3.0 * xc_pow2 * zc_pow3;
            phi_x_tmp[288 + i] += S0[i] * AX;
            AZ = 3.0 * xc_pow3 * zc_pow2;
            phi_z_tmp[288 + i] += S0[i] * AZ;

            // Density AM=6 Component=XXYYYY
            A = xc_pow2 * yc_pow4;
            phi_tmp[320 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXYYYY
            phi_x_tmp[320 + i] = SX * A;
            phi_y_tmp[320 + i] = SY * A;
            phi_z_tmp[320 + i] = SZ * A;
            AX = 2.0 * xc[i] * yc_pow4;
            phi_x_tmp[320 + i] += S0[i] * AX;
            AY = 4.0 * xc_pow2 * yc_pow3;
            phi_y_tmp[320 + i] += S0[i] * AY;

            // Density AM=6 Component=XXYYYZ
            A = xc_pow2 * yc_pow3 * zc[i];
            phi_tmp[352 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXYYYZ
            phi_x_tmp[352 + i] = SX * A;
            phi_y_tmp[352 + i] = SY * A;
            phi_z_tmp[352 + i] = SZ * A;
            AX = 2.0 * xc[i] * yc_pow3 * zc[i];
            phi_x_tmp[352 + i] += S0[i] * AX;
            AY = 3.0 * xc_pow2 * yc_pow2 * zc[i];
            phi_y_tmp[352 + i] += S0[i] * AY;
            AZ = xc_pow2 * yc_pow3;
            phi_z_tmp[352 + i] += S0[i] * AZ;

            // Density AM=6 Component=XXYYZZ
            A = xc_pow2 * yc_pow2 * zc_pow2;
            phi_tmp[384 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXYYZZ
            phi_x_tmp[384 + i] = SX * A;
            phi_y_tmp[384 + i] = SY * A;
            phi_z_tmp[384 + i] = SZ * A;
            AX = 2.0 * xc[i] * yc_pow2 * zc_pow2;
            phi_x_tmp[384 + i] += S0[i] * AX;
            AY = 2.0 * xc_pow2 * yc[i] * zc_pow2;
            phi_y_tmp[384 + i] += S0[i] * AY;
            AZ = 2.0 * xc_pow2 * yc_pow2 * zc[i];
            phi_z_tmp[384 + i] += S0[i] * AZ;

            // Density AM=6 Component=XXYZZZ
            A = xc_pow2 * yc[i] * zc_pow3;
            phi_tmp[416 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXYZZZ
            phi_x_tmp[416 + i] = SX * A;
            phi_y_tmp[416 + i] = SY * A;
            phi_z_tmp[416 + i] = SZ * A;
            AX = 2.0 * xc[i] * yc[i] * zc_pow3;
            phi_x_tmp[416 + i] += S0[i] * AX;
            AY = xc_pow2 * zc_pow3;
            phi_y_tmp[416 + i] += S0[i] * AY;
            AZ = 3.0 * xc_pow2 * yc[i] * zc_pow2;
            phi_z_tmp[416 + i] += S0[i] * AZ;

            // Density AM=6 Component=XXZZZZ
            A = xc_pow2 * zc_pow4;
            phi_tmp[448 + i] = S0[i] * A;

            // Gradient AM=6 Component=XXZZZZ
            phi_x_tmp[448 + i] = SX * A;
            phi_y_tmp[448 + i] = SY * A;
            phi_z_tmp[448 + i] = SZ * A;
            AX = 2.0 * xc[i] * zc_pow4;
            phi_x_tmp[448 + i] += S0[i] * AX;
            AZ = 4.0 * xc_pow2 * zc_pow3;
            phi_z_tmp[448 + i] += S0[i] * AZ;

            // Density AM=6 Component=XYYYYY
            A = xc[i] * yc_pow5;
            phi_tmp[480 + i] = S0[i] * A;

            // Gradient AM=6 Component=XYYYYY
            phi_x_tmp[480 + i] = SX * A;
            phi_y_tmp[480 + i] = SY * A;
            phi_z_tmp[480 + i] = SZ * A;
            phi_x_tmp[480 + i] += S0[i] * yc_pow5;
            AY = 5.0 * xc[i] * yc_pow4;
            phi_y_tmp[480 + i] += S0[i] * AY;

            // Density AM=6 Component=XYYYYZ
            A = xc[i] * yc_pow4 * zc[i];
            phi_tmp[512 + i] = S0[i] * A;

            // Gradient AM=6 Component=XYYYYZ
            phi_x_tmp[512 + i] = SX * A;
            phi_y_tmp[512 + i] = SY * A;
            phi_z_tmp[512 + i] = SZ * A;
            AX = yc_pow4 * zc[i];
            phi_x_tmp[512 + i] += S0[i] * AX;
            AY = 4.0 * xc[i] * yc_pow3 * zc[i];
            phi_y_tmp[512 + i] += S0[i] * AY;
            AZ = xc[i] * yc_pow4;
            phi_z_tmp[512 + i] += S0[i] * AZ;

            // Density AM=6 Component=XYYYZZ
            A = xc[i] * yc_pow3 * zc_pow2;
            phi_tmp[544 + i] = S0[i] * A;

            // Gradient AM=6 Component=XYYYZZ
            phi_x_tmp[544 + i] = SX * A;
            phi_y_tmp[544 + i] = SY * A;
            phi_z_tmp[544 + i] = SZ * A;
            AX = yc_pow3 * zc_pow2;
            phi_x_tmp[544 + i] += S0[i] * AX;
            AY = 3.0 * xc[i] * yc_pow2 * zc_pow2;
            phi_y_tmp[544 + i] += S0[i] * AY;
            AZ = 2.0 * xc[i] * yc_pow3 * zc[i];
            phi_z_tmp[544 + i] += S0[i] * AZ;

            // Density AM=6 Component=XYYZZZ
            A = xc[i] * yc_pow2 * zc_pow3;
            phi_tmp[576 + i] = S0[i] * A;

            // Gradient AM=6 Component=XYYZZZ
            phi_x_tmp[576 + i] = SX * A;
            phi_y_tmp[576 + i] = SY * A;
            phi_z_tmp[576 + i] = SZ * A;
            AX = yc_pow2 * zc_pow3;
            phi_x_tmp[576 + i] += S0[i] * AX;
            AY = 2.0 * xc[i] * yc[i] * zc_pow3;
            phi_y_tmp[576 + i] += S0[i] * AY;
            AZ = 3.0 * xc[i] * yc_pow2 * zc_pow2;
            phi_z_tmp[576 + i] += S0[i] * AZ;

            // Density AM=6 Component=XYZZZZ
            A = xc[i] * yc[i] * zc_pow4;
            phi_tmp[608 + i] = S0[i] * A;

            // Gradient AM=6 Component=XYZZZZ
            phi_x_tmp[608 + i] = SX * A;
            phi_y_tmp[608 + i] = SY * A;
            phi_z_tmp[608 + i] = SZ * A;
            AX = yc[i] * zc_pow4;
            phi_x_tmp[608 + i] += S0[i] * AX;
            AY = xc[i] * zc_pow4;
            phi_y_tmp[608 + i] += S0[i] * AY;
            AZ = 4.0 * xc[i] * yc[i] * zc_pow3;
            phi_z_tmp[608 + i] += S0[i] * AZ;

            // Density AM=6 Component=XZZZZZ
            A = xc[i] * zc_pow5;
            phi_tmp[640 + i] = S0[i] * A;

            // Gradient AM=6 Component=XZZZZZ
            phi_x_tmp[640 + i] = SX * A;
            phi_y_tmp[640 + i] = SY * A;
            phi_z_tmp[640 + i] = SZ * A;
            phi_x_tmp[640 + i] += S0[i] * zc_pow5;
            AZ = 5.0 * xc[i] * zc_pow4;
            phi_z_tmp[640 + i] += S0[i] * AZ;

            // Density AM=6 Component=YYYYYY
            phi_tmp[672 + i] = S0[i] * yc_pow6;

            // Gradient AM=6 Component=YYYYYY
            phi_x_tmp[672 + i] = SX * yc_pow6;
            phi_y_tmp[672 + i] = SY * yc_pow6;
            phi_z_tmp[672 + i] = SZ * yc_pow6;
            AY = 6.0 * yc_pow5;
            phi_y_tmp[672 + i] += S0[i] * AY;

            // Density AM=6 Component=YYYYYZ
            A = yc_pow5 * zc[i];
            phi_tmp[704 + i] = S0[i] * A;

            // Gradient AM=6 Component=YYYYYZ
            phi_x_tmp[704 + i] = SX * A;
            phi_y_tmp[704 + i] = SY * A;
            phi_z_tmp[704 + i] = SZ * A;
            AY = 5.0 * yc_pow4 * zc[i];
            phi_y_tmp[704 + i] += S0[i] * AY;
            phi_z_tmp[704 + i] += S0[i] * yc_pow5;

            // Density AM=6 Component=YYYYZZ
            A = yc_pow4 * zc_pow2;
            phi_tmp[736 + i] = S0[i] * A;

            // Gradient AM=6 Component=YYYYZZ
            phi_x_tmp[736 + i] = SX * A;
            phi_y_tmp[736 + i] = SY * A;
            phi_z_tmp[736 + i] = SZ * A;
            AY = 4.0 * yc_pow3 * zc_pow2;
            phi_y_tmp[736 + i] += S0[i] * AY;
            AZ = 2.0 * yc_pow4 * zc[i];
            phi_z_tmp[736 + i] += S0[i] * AZ;

            // Density AM=6 Component=YYYZZZ
            A = yc_pow3 * zc_pow3;
            phi_tmp[768 + i] = S0[i] * A;

            // Gradient AM=6 Component=YYYZZZ
            phi_x_tmp[768 + i] = SX * A;
            phi_y_tmp[768 + i] = SY * A;
            phi_z_tmp[768 + i] = SZ * A;
            AY = 3.0 * yc_pow2 * zc_pow3;
            phi_y_tmp[768 + i] += S0[i] * AY;
            AZ = 3.0 * yc_pow3 * zc_pow2;
            phi_z_tmp[768 + i] += S0[i] * AZ;

            // Density AM=6 Component=YYZZZZ
            A = yc_pow2 * zc_pow4;
            phi_tmp[800 + i] = S0[i] * A;

            // Gradient AM=6 Component=YYZZZZ
            phi_x_tmp[800 + i] = SX * A;
            phi_y_tmp[800 + i] = SY * A;
            phi_z_tmp[800 + i] = SZ * A;
            AY = 2.0 * yc[i] * zc_pow4;
            phi_y_tmp[800 + i] += S0[i] * AY;
            AZ = 4.0 * yc_pow2 * zc_pow3;
            phi_z_tmp[800 + i] += S0[i] * AZ;

            // Density AM=6 Component=YZZZZZ
            A = yc[i] * zc_pow5;
            phi_tmp[832 + i] = S0[i] * A;

            // Gradient AM=6 Component=YZZZZZ
            phi_x_tmp[832 + i] = SX * A;
            phi_y_tmp[832 + i] = SY * A;
            phi_z_tmp[832 + i] = SZ * A;
            phi_y_tmp[832 + i] += S0[i] * zc_pow5;
            AZ = 5.0 * yc[i] * zc_pow4;
            phi_z_tmp[832 + i] += S0[i] * AZ;

            // Density AM=6 Component=ZZZZZZ
            phi_tmp[864 + i] = S0[i] * zc_pow6;

            // Gradient AM=6 Component=ZZZZZZ
            phi_x_tmp[864 + i] = SX * zc_pow6;
            phi_y_tmp[864 + i] = SY * zc_pow6;
            phi_z_tmp[864 + i] = SZ * zc_pow6;
            AZ = 6.0 * zc_pow5;
            phi_z_tmp[864 + i] += S0[i] * AZ;

        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_to_spherical_L6(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_to_spherical_L6(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_to_spherical_L6(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_gaussian_cart_to_spherical_L6(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_gaussian_cart_to_spherical_L6(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_gaussian_cart_to_spherical_L6(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_cca_cart_copy_L6(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_cca_cart_copy_L6(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_cca_cart_copy_L6(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
            } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_out + start), npoints);

            // Gradient, transform data to outer temps
            gg_molden_cart_copy_L6(remain, phi_x_tmp, 32, (phi_x_out + start), npoints);
            gg_molden_cart_copy_L6(remain, phi_y_tmp, 32, (phi_y_out + start), npoints);
            gg_molden_cart_copy_L6(remain, phi_z_tmp, 32, (phi_z_out + start), npoints);
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

}
