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

void gg_collocation_L0(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride,
                       const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents,
                       const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out) {
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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, ((((192 * sizeof(double)) + 64 - 1) / 64) * 64));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, ((((nprim * sizeof(double)) + 64 - 1) / 64) * 64));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, ((((32 * sizeof(double)) + 64 - 1) / 64) * 64));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
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
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
            }
        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            phi_out[start + i] = S0[i];
        }
    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);
}

void gg_collocation_L1(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride,
                       const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents,
                       const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out) {
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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, ((((192 * sizeof(double)) + 64 - 1) / 64) * 64));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, ((((nprim * sizeof(double)) + 64 - 1) / 64) * 64));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, ((((96 * sizeof(double)) + 64 - 1) / 64) * 64));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
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
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
            }
        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            // Density AM=1 Component=X
            phi_tmp[i] = S0[i] * xc[i];

            // Density AM=1 Component=Y
            phi_tmp[32 + i] = S0[i] * yc[i];

            // Density AM=1 Component=Z
            phi_tmp[64 + i] = S0[i] * zc[i];
        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L1(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L1(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L1(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L1(remain, phi_tmp, 32, (phi_out + start), npoints);
        }
    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);
}

void gg_collocation_L2(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride,
                       const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents,
                       const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out) {
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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, ((((192 * sizeof(double)) + 64 - 1) / 64) * 64));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, ((((nprim * sizeof(double)) + 64 - 1) / 64) * 64));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, ((((192 * sizeof(double)) + 64 - 1) / 64) * 64));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
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
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
            }
        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            // Cartesian derivs
            const double xc_pow2 = xc[i] * xc[i];
            const double yc_pow2 = yc[i] * yc[i];
            const double zc_pow2 = zc[i] * zc[i];

            // Density AM=2 Component=XX
            phi_tmp[i] = S0[i] * xc_pow2;

            // Density AM=2 Component=XY
            A = xc[i] * yc[i];
            phi_tmp[32 + i] = S0[i] * A;

            // Density AM=2 Component=XZ
            A = xc[i] * zc[i];
            phi_tmp[64 + i] = S0[i] * A;

            // Density AM=2 Component=YY
            phi_tmp[96 + i] = S0[i] * yc_pow2;

            // Density AM=2 Component=YZ
            A = yc[i] * zc[i];
            phi_tmp[128 + i] = S0[i] * A;

            // Density AM=2 Component=ZZ
            phi_tmp[160 + i] = S0[i] * zc_pow2;
        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L2(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L2(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L2(remain, phi_tmp, 32, (phi_out + start), npoints);
        }
    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);
}

void gg_collocation_L3(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride,
                       const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents,
                       const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out) {
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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, ((((192 * sizeof(double)) + 64 - 1) / 64) * 64));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, ((((nprim * sizeof(double)) + 64 - 1) / 64) * 64));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, ((((320 * sizeof(double)) + 64 - 1) / 64) * 64));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
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
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
            }
        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
            // Cartesian derivs
            const double xc_pow2 = xc[i] * xc[i];
            const double yc_pow2 = yc[i] * yc[i];
            const double zc_pow2 = zc[i] * zc[i];

            const double xc_pow3 = xc_pow2 * xc[i];
            const double yc_pow3 = yc_pow2 * yc[i];
            const double zc_pow3 = zc_pow2 * zc[i];

            // Density AM=3 Component=XXX
            phi_tmp[i] = S0[i] * xc_pow3;

            // Density AM=3 Component=XXY
            A = xc_pow2 * yc[i];
            phi_tmp[32 + i] = S0[i] * A;

            // Density AM=3 Component=XXZ
            A = xc_pow2 * zc[i];
            phi_tmp[64 + i] = S0[i] * A;

            // Density AM=3 Component=XYY
            A = xc[i] * yc_pow2;
            phi_tmp[96 + i] = S0[i] * A;

            // Density AM=3 Component=XYZ
            A = xc[i] * yc[i] * zc[i];
            phi_tmp[128 + i] = S0[i] * A;

            // Density AM=3 Component=XZZ
            A = xc[i] * zc_pow2;
            phi_tmp[160 + i] = S0[i] * A;

            // Density AM=3 Component=YYY
            phi_tmp[192 + i] = S0[i] * yc_pow3;

            // Density AM=3 Component=YYZ
            A = yc_pow2 * zc[i];
            phi_tmp[224 + i] = S0[i] * A;

            // Density AM=3 Component=YZZ
            A = yc[i] * zc_pow2;
            phi_tmp[256 + i] = S0[i] * A;

            // Density AM=3 Component=ZZZ
            phi_tmp[288 + i] = S0[i] * zc_pow3;
        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L3(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L3(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L3(remain, phi_tmp, 32, (phi_out + start), npoints);
        }
    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);
}

void gg_collocation_L4(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride,
                       const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents,
                       const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out) {
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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, ((((192 * sizeof(double)) + 64 - 1) / 64) * 64));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, ((((nprim * sizeof(double)) + 64 - 1) / 64) * 64));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, ((((480 * sizeof(double)) + 64 - 1) / 64) * 64));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
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
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
            }
        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
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

            // Density AM=4 Component=XXXY
            A = xc_pow3 * yc[i];
            phi_tmp[32 + i] = S0[i] * A;

            // Density AM=4 Component=XXXZ
            A = xc_pow3 * zc[i];
            phi_tmp[64 + i] = S0[i] * A;

            // Density AM=4 Component=XXYY
            A = xc_pow2 * yc_pow2;
            phi_tmp[96 + i] = S0[i] * A;

            // Density AM=4 Component=XXYZ
            A = xc_pow2 * yc[i] * zc[i];
            phi_tmp[128 + i] = S0[i] * A;

            // Density AM=4 Component=XXZZ
            A = xc_pow2 * zc_pow2;
            phi_tmp[160 + i] = S0[i] * A;

            // Density AM=4 Component=XYYY
            A = xc[i] * yc_pow3;
            phi_tmp[192 + i] = S0[i] * A;

            // Density AM=4 Component=XYYZ
            A = xc[i] * yc_pow2 * zc[i];
            phi_tmp[224 + i] = S0[i] * A;

            // Density AM=4 Component=XYZZ
            A = xc[i] * yc[i] * zc_pow2;
            phi_tmp[256 + i] = S0[i] * A;

            // Density AM=4 Component=XZZZ
            A = xc[i] * zc_pow3;
            phi_tmp[288 + i] = S0[i] * A;

            // Density AM=4 Component=YYYY
            phi_tmp[320 + i] = S0[i] * yc_pow4;

            // Density AM=4 Component=YYYZ
            A = yc_pow3 * zc[i];
            phi_tmp[352 + i] = S0[i] * A;

            // Density AM=4 Component=YYZZ
            A = yc_pow2 * zc_pow2;
            phi_tmp[384 + i] = S0[i] * A;

            // Density AM=4 Component=YZZZ
            A = yc[i] * zc_pow3;
            phi_tmp[416 + i] = S0[i] * A;

            // Density AM=4 Component=ZZZZ
            phi_tmp[448 + i] = S0[i] * zc_pow4;
        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L4(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L4(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L4(remain, phi_tmp, 32, (phi_out + start), npoints);
        }
    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);
}

void gg_collocation_L5(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride,
                       const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents,
                       const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out) {
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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, ((((192 * sizeof(double)) + 64 - 1) / 64) * 64));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, ((((nprim * sizeof(double)) + 64 - 1) / 64) * 64));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, ((((672 * sizeof(double)) + 64 - 1) / 64) * 64));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
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
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
            }
        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
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

            // Density AM=5 Component=XXXXY
            A = xc_pow4 * yc[i];
            phi_tmp[32 + i] = S0[i] * A;

            // Density AM=5 Component=XXXXZ
            A = xc_pow4 * zc[i];
            phi_tmp[64 + i] = S0[i] * A;

            // Density AM=5 Component=XXXYY
            A = xc_pow3 * yc_pow2;
            phi_tmp[96 + i] = S0[i] * A;

            // Density AM=5 Component=XXXYZ
            A = xc_pow3 * yc[i] * zc[i];
            phi_tmp[128 + i] = S0[i] * A;

            // Density AM=5 Component=XXXZZ
            A = xc_pow3 * zc_pow2;
            phi_tmp[160 + i] = S0[i] * A;

            // Density AM=5 Component=XXYYY
            A = xc_pow2 * yc_pow3;
            phi_tmp[192 + i] = S0[i] * A;

            // Density AM=5 Component=XXYYZ
            A = xc_pow2 * yc_pow2 * zc[i];
            phi_tmp[224 + i] = S0[i] * A;

            // Density AM=5 Component=XXYZZ
            A = xc_pow2 * yc[i] * zc_pow2;
            phi_tmp[256 + i] = S0[i] * A;

            // Density AM=5 Component=XXZZZ
            A = xc_pow2 * zc_pow3;
            phi_tmp[288 + i] = S0[i] * A;

            // Density AM=5 Component=XYYYY
            A = xc[i] * yc_pow4;
            phi_tmp[320 + i] = S0[i] * A;

            // Density AM=5 Component=XYYYZ
            A = xc[i] * yc_pow3 * zc[i];
            phi_tmp[352 + i] = S0[i] * A;

            // Density AM=5 Component=XYYZZ
            A = xc[i] * yc_pow2 * zc_pow2;
            phi_tmp[384 + i] = S0[i] * A;

            // Density AM=5 Component=XYZZZ
            A = xc[i] * yc[i] * zc_pow3;
            phi_tmp[416 + i] = S0[i] * A;

            // Density AM=5 Component=XZZZZ
            A = xc[i] * zc_pow4;
            phi_tmp[448 + i] = S0[i] * A;

            // Density AM=5 Component=YYYYY
            phi_tmp[480 + i] = S0[i] * yc_pow5;

            // Density AM=5 Component=YYYYZ
            A = yc_pow4 * zc[i];
            phi_tmp[512 + i] = S0[i] * A;

            // Density AM=5 Component=YYYZZ
            A = yc_pow3 * zc_pow2;
            phi_tmp[544 + i] = S0[i] * A;

            // Density AM=5 Component=YYZZZ
            A = yc_pow2 * zc_pow3;
            phi_tmp[576 + i] = S0[i] * A;

            // Density AM=5 Component=YZZZZ
            A = yc[i] * zc_pow4;
            phi_tmp[608 + i] = S0[i] * A;

            // Density AM=5 Component=ZZZZZ
            phi_tmp[640 + i] = S0[i] * zc_pow5;
        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L5(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L5(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L5(remain, phi_tmp, 32, (phi_out + start), npoints);
        }
    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);
}

void gg_collocation_L6(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride,
                       const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents,
                       const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out) {
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
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, ((((192 * sizeof(double)) + 64 - 1) / 64) * 64));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, ((((nprim * sizeof(double)) + 64 - 1) / 64) * 64));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, ((((896 * sizeof(double)) + 64 - 1) / 64) * 64));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
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
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
            }
        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
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

            // Density AM=6 Component=XXXXXY
            A = xc_pow5 * yc[i];
            phi_tmp[32 + i] = S0[i] * A;

            // Density AM=6 Component=XXXXXZ
            A = xc_pow5 * zc[i];
            phi_tmp[64 + i] = S0[i] * A;

            // Density AM=6 Component=XXXXYY
            A = xc_pow4 * yc_pow2;
            phi_tmp[96 + i] = S0[i] * A;

            // Density AM=6 Component=XXXXYZ
            A = xc_pow4 * yc[i] * zc[i];
            phi_tmp[128 + i] = S0[i] * A;

            // Density AM=6 Component=XXXXZZ
            A = xc_pow4 * zc_pow2;
            phi_tmp[160 + i] = S0[i] * A;

            // Density AM=6 Component=XXXYYY
            A = xc_pow3 * yc_pow3;
            phi_tmp[192 + i] = S0[i] * A;

            // Density AM=6 Component=XXXYYZ
            A = xc_pow3 * yc_pow2 * zc[i];
            phi_tmp[224 + i] = S0[i] * A;

            // Density AM=6 Component=XXXYZZ
            A = xc_pow3 * yc[i] * zc_pow2;
            phi_tmp[256 + i] = S0[i] * A;

            // Density AM=6 Component=XXXZZZ
            A = xc_pow3 * zc_pow3;
            phi_tmp[288 + i] = S0[i] * A;

            // Density AM=6 Component=XXYYYY
            A = xc_pow2 * yc_pow4;
            phi_tmp[320 + i] = S0[i] * A;

            // Density AM=6 Component=XXYYYZ
            A = xc_pow2 * yc_pow3 * zc[i];
            phi_tmp[352 + i] = S0[i] * A;

            // Density AM=6 Component=XXYYZZ
            A = xc_pow2 * yc_pow2 * zc_pow2;
            phi_tmp[384 + i] = S0[i] * A;

            // Density AM=6 Component=XXYZZZ
            A = xc_pow2 * yc[i] * zc_pow3;
            phi_tmp[416 + i] = S0[i] * A;

            // Density AM=6 Component=XXZZZZ
            A = xc_pow2 * zc_pow4;
            phi_tmp[448 + i] = S0[i] * A;

            // Density AM=6 Component=XYYYYY
            A = xc[i] * yc_pow5;
            phi_tmp[480 + i] = S0[i] * A;

            // Density AM=6 Component=XYYYYZ
            A = xc[i] * yc_pow4 * zc[i];
            phi_tmp[512 + i] = S0[i] * A;

            // Density AM=6 Component=XYYYZZ
            A = xc[i] * yc_pow3 * zc_pow2;
            phi_tmp[544 + i] = S0[i] * A;

            // Density AM=6 Component=XYYZZZ
            A = xc[i] * yc_pow2 * zc_pow3;
            phi_tmp[576 + i] = S0[i] * A;

            // Density AM=6 Component=XYZZZZ
            A = xc[i] * yc[i] * zc_pow4;
            phi_tmp[608 + i] = S0[i] * A;

            // Density AM=6 Component=XZZZZZ
            A = xc[i] * zc_pow5;
            phi_tmp[640 + i] = S0[i] * A;

            // Density AM=6 Component=YYYYYY
            phi_tmp[672 + i] = S0[i] * yc_pow6;

            // Density AM=6 Component=YYYYYZ
            A = yc_pow5 * zc[i];
            phi_tmp[704 + i] = S0[i] * A;

            // Density AM=6 Component=YYYYZZ
            A = yc_pow4 * zc_pow2;
            phi_tmp[736 + i] = S0[i] * A;

            // Density AM=6 Component=YYYZZZ
            A = yc_pow3 * zc_pow3;
            phi_tmp[768 + i] = S0[i] * A;

            // Density AM=6 Component=YYZZZZ
            A = yc_pow2 * zc_pow4;
            phi_tmp[800 + i] = S0[i] * A;

            // Density AM=6 Component=YZZZZZ
            A = yc[i] * zc_pow5;
            phi_tmp[832 + i] = S0[i] * A;

            // Density AM=6 Component=ZZZZZZ
            phi_tmp[864 + i] = S0[i] * zc_pow6;
        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L6(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L6(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L6(remain, phi_tmp, 32, (phi_out + start), npoints);
        }
    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);
}

void gg_collocation_L7(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride,
                       const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents,
                       const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out) {
    // Sizing
    unsigned long nblocks = npoints / 32;
    nblocks += (npoints % 32) ? 1 : 0;
    const unsigned long ncart = 36;
    const unsigned long nspherical = 15;
    unsigned long nout;

    if ((order == GG_SPHERICAL_CCA) || (order == GG_SPHERICAL_GAUSSIAN)) {
        nout = nspherical;
    } else {
        nout = ncart;
    }

    // Allocate S temporaries, single block to stay on cache
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, ((((192 * sizeof(double)) + 64 - 1) / 64) * 64));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, ((((nprim * sizeof(double)) + 64 - 1) / 64) * 64));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, ((((1152 * sizeof(double)) + 64 - 1) / 64) * 64));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
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
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
            }
        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
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

            const double xc_pow7 = xc_pow6 * xc[i];
            const double yc_pow7 = yc_pow6 * yc[i];
            const double zc_pow7 = zc_pow6 * zc[i];

            // Density AM=7 Component=XXXXXXX
            phi_tmp[i] = S0[i] * xc_pow7;

            // Density AM=7 Component=XXXXXXY
            A = xc_pow6 * yc[i];
            phi_tmp[32 + i] = S0[i] * A;

            // Density AM=7 Component=XXXXXXZ
            A = xc_pow6 * zc[i];
            phi_tmp[64 + i] = S0[i] * A;

            // Density AM=7 Component=XXXXXYY
            A = xc_pow5 * yc_pow2;
            phi_tmp[96 + i] = S0[i] * A;

            // Density AM=7 Component=XXXXXYZ
            A = xc_pow5 * yc[i] * zc[i];
            phi_tmp[128 + i] = S0[i] * A;

            // Density AM=7 Component=XXXXXZZ
            A = xc_pow5 * zc_pow2;
            phi_tmp[160 + i] = S0[i] * A;

            // Density AM=7 Component=XXXXYYY
            A = xc_pow4 * yc_pow3;
            phi_tmp[192 + i] = S0[i] * A;

            // Density AM=7 Component=XXXXYYZ
            A = xc_pow4 * yc_pow2 * zc[i];
            phi_tmp[224 + i] = S0[i] * A;

            // Density AM=7 Component=XXXXYZZ
            A = xc_pow4 * yc[i] * zc_pow2;
            phi_tmp[256 + i] = S0[i] * A;

            // Density AM=7 Component=XXXXZZZ
            A = xc_pow4 * zc_pow3;
            phi_tmp[288 + i] = S0[i] * A;

            // Density AM=7 Component=XXXYYYY
            A = xc_pow3 * yc_pow4;
            phi_tmp[320 + i] = S0[i] * A;

            // Density AM=7 Component=XXXYYYZ
            A = xc_pow3 * yc_pow3 * zc[i];
            phi_tmp[352 + i] = S0[i] * A;

            // Density AM=7 Component=XXXYYZZ
            A = xc_pow3 * yc_pow2 * zc_pow2;
            phi_tmp[384 + i] = S0[i] * A;

            // Density AM=7 Component=XXXYZZZ
            A = xc_pow3 * yc[i] * zc_pow3;
            phi_tmp[416 + i] = S0[i] * A;

            // Density AM=7 Component=XXXZZZZ
            A = xc_pow3 * zc_pow4;
            phi_tmp[448 + i] = S0[i] * A;

            // Density AM=7 Component=XXYYYYY
            A = xc_pow2 * yc_pow5;
            phi_tmp[480 + i] = S0[i] * A;

            // Density AM=7 Component=XXYYYYZ
            A = xc_pow2 * yc_pow4 * zc[i];
            phi_tmp[512 + i] = S0[i] * A;

            // Density AM=7 Component=XXYYYZZ
            A = xc_pow2 * yc_pow3 * zc_pow2;
            phi_tmp[544 + i] = S0[i] * A;

            // Density AM=7 Component=XXYYZZZ
            A = xc_pow2 * yc_pow2 * zc_pow3;
            phi_tmp[576 + i] = S0[i] * A;

            // Density AM=7 Component=XXYZZZZ
            A = xc_pow2 * yc[i] * zc_pow4;
            phi_tmp[608 + i] = S0[i] * A;

            // Density AM=7 Component=XXZZZZZ
            A = xc_pow2 * zc_pow5;
            phi_tmp[640 + i] = S0[i] * A;

            // Density AM=7 Component=XYYYYYY
            A = xc[i] * yc_pow6;
            phi_tmp[672 + i] = S0[i] * A;

            // Density AM=7 Component=XYYYYYZ
            A = xc[i] * yc_pow5 * zc[i];
            phi_tmp[704 + i] = S0[i] * A;

            // Density AM=7 Component=XYYYYZZ
            A = xc[i] * yc_pow4 * zc_pow2;
            phi_tmp[736 + i] = S0[i] * A;

            // Density AM=7 Component=XYYYZZZ
            A = xc[i] * yc_pow3 * zc_pow3;
            phi_tmp[768 + i] = S0[i] * A;

            // Density AM=7 Component=XYYZZZZ
            A = xc[i] * yc_pow2 * zc_pow4;
            phi_tmp[800 + i] = S0[i] * A;

            // Density AM=7 Component=XYZZZZZ
            A = xc[i] * yc[i] * zc_pow5;
            phi_tmp[832 + i] = S0[i] * A;

            // Density AM=7 Component=XZZZZZZ
            A = xc[i] * zc_pow6;
            phi_tmp[864 + i] = S0[i] * A;

            // Density AM=7 Component=YYYYYYY
            phi_tmp[896 + i] = S0[i] * yc_pow7;

            // Density AM=7 Component=YYYYYYZ
            A = yc_pow6 * zc[i];
            phi_tmp[928 + i] = S0[i] * A;

            // Density AM=7 Component=YYYYYZZ
            A = yc_pow5 * zc_pow2;
            phi_tmp[960 + i] = S0[i] * A;

            // Density AM=7 Component=YYYYZZZ
            A = yc_pow4 * zc_pow3;
            phi_tmp[992 + i] = S0[i] * A;

            // Density AM=7 Component=YYYZZZZ
            A = yc_pow3 * zc_pow4;
            phi_tmp[1024 + i] = S0[i] * A;

            // Density AM=7 Component=YYZZZZZ
            A = yc_pow2 * zc_pow5;
            phi_tmp[1056 + i] = S0[i] * A;

            // Density AM=7 Component=YZZZZZZ
            A = yc[i] * zc_pow6;
            phi_tmp[1088 + i] = S0[i] * A;

            // Density AM=7 Component=ZZZZZZZ
            phi_tmp[1120 + i] = S0[i] * zc_pow7;
        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L7(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L7(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L7(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L7(remain, phi_tmp, 32, (phi_out + start), npoints);
        }
    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);
}

void gg_collocation_L8(const unsigned long npoints, const double* PRAGMA_RESTRICT xyz, const unsigned long xyz_stride,
                       const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents,
                       const double* PRAGMA_RESTRICT center, const int order, double* PRAGMA_RESTRICT phi_out) {
    // Sizing
    unsigned long nblocks = npoints / 32;
    nblocks += (npoints % 32) ? 1 : 0;
    const unsigned long ncart = 45;
    const unsigned long nspherical = 17;
    unsigned long nout;

    if ((order == GG_SPHERICAL_CCA) || (order == GG_SPHERICAL_GAUSSIAN)) {
        nout = nspherical;
    } else {
        nout = ncart;
    }

    // Allocate S temporaries, single block to stay on cache
    double* PRAGMA_RESTRICT cache_data = (double*)ALIGNED_MALLOC(64, ((((192 * sizeof(double)) + 64 - 1) / 64) * 64));
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

    // Allocate exponential temporaries
    double* PRAGMA_RESTRICT expn1 = (double*)ALIGNED_MALLOC(64, ((((nprim * sizeof(double)) + 64 - 1) / 64) * 64));

    // Allocate output temporaries
    double* PRAGMA_RESTRICT phi_tmp = (double*)ALIGNED_MALLOC(64, ((((1440 * sizeof(double)) + 64 - 1) / 64) * 64));
    ASSUME_ALIGNED(phi_tmp, 64);

    // Declare doubles
    const double center_x = center[0];
    const double center_y = center[1];
    const double center_z = center[2];
    double A;

    // Build negative exponents
    for (unsigned long i = 0; i < nprim; i++) {
        expn1[i] = -1.0 * exponents[i];
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
            }
        }

        // Start exponential block loop
        for (unsigned long n = 0; n < nprim; n++) {
            const double coef = coeffs[n];
            const double alpha_n1 = expn1[n];

            PRAGMA_VECTORIZE
            for (unsigned long i = 0; i < remain; i++) {
                const double width = alpha_n1 * R2[i];
                const double T1 = coef * exp(width);
                S0[i] += T1;
            }
        }

        // Combine blocks
        PRAGMA_VECTORIZE
        for (unsigned long i = 0; i < remain; i++) {
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

            const double xc_pow7 = xc_pow6 * xc[i];
            const double yc_pow7 = yc_pow6 * yc[i];
            const double zc_pow7 = zc_pow6 * zc[i];

            const double xc_pow8 = xc_pow7 * xc[i];
            const double yc_pow8 = yc_pow7 * yc[i];
            const double zc_pow8 = zc_pow7 * zc[i];

            // Density AM=8 Component=XXXXXXXX
            phi_tmp[i] = S0[i] * xc_pow8;

            // Density AM=8 Component=XXXXXXXY
            A = xc_pow7 * yc[i];
            phi_tmp[32 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXXXXZ
            A = xc_pow7 * zc[i];
            phi_tmp[64 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXXXYY
            A = xc_pow6 * yc_pow2;
            phi_tmp[96 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXXXYZ
            A = xc_pow6 * yc[i] * zc[i];
            phi_tmp[128 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXXXZZ
            A = xc_pow6 * zc_pow2;
            phi_tmp[160 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXXYYY
            A = xc_pow5 * yc_pow3;
            phi_tmp[192 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXXYYZ
            A = xc_pow5 * yc_pow2 * zc[i];
            phi_tmp[224 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXXYZZ
            A = xc_pow5 * yc[i] * zc_pow2;
            phi_tmp[256 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXXZZZ
            A = xc_pow5 * zc_pow3;
            phi_tmp[288 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXYYYY
            A = xc_pow4 * yc_pow4;
            phi_tmp[320 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXYYYZ
            A = xc_pow4 * yc_pow3 * zc[i];
            phi_tmp[352 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXYYZZ
            A = xc_pow4 * yc_pow2 * zc_pow2;
            phi_tmp[384 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXYZZZ
            A = xc_pow4 * yc[i] * zc_pow3;
            phi_tmp[416 + i] = S0[i] * A;

            // Density AM=8 Component=XXXXZZZZ
            A = xc_pow4 * zc_pow4;
            phi_tmp[448 + i] = S0[i] * A;

            // Density AM=8 Component=XXXYYYYY
            A = xc_pow3 * yc_pow5;
            phi_tmp[480 + i] = S0[i] * A;

            // Density AM=8 Component=XXXYYYYZ
            A = xc_pow3 * yc_pow4 * zc[i];
            phi_tmp[512 + i] = S0[i] * A;

            // Density AM=8 Component=XXXYYYZZ
            A = xc_pow3 * yc_pow3 * zc_pow2;
            phi_tmp[544 + i] = S0[i] * A;

            // Density AM=8 Component=XXXYYZZZ
            A = xc_pow3 * yc_pow2 * zc_pow3;
            phi_tmp[576 + i] = S0[i] * A;

            // Density AM=8 Component=XXXYZZZZ
            A = xc_pow3 * yc[i] * zc_pow4;
            phi_tmp[608 + i] = S0[i] * A;

            // Density AM=8 Component=XXXZZZZZ
            A = xc_pow3 * zc_pow5;
            phi_tmp[640 + i] = S0[i] * A;

            // Density AM=8 Component=XXYYYYYY
            A = xc_pow2 * yc_pow6;
            phi_tmp[672 + i] = S0[i] * A;

            // Density AM=8 Component=XXYYYYYZ
            A = xc_pow2 * yc_pow5 * zc[i];
            phi_tmp[704 + i] = S0[i] * A;

            // Density AM=8 Component=XXYYYYZZ
            A = xc_pow2 * yc_pow4 * zc_pow2;
            phi_tmp[736 + i] = S0[i] * A;

            // Density AM=8 Component=XXYYYZZZ
            A = xc_pow2 * yc_pow3 * zc_pow3;
            phi_tmp[768 + i] = S0[i] * A;

            // Density AM=8 Component=XXYYZZZZ
            A = xc_pow2 * yc_pow2 * zc_pow4;
            phi_tmp[800 + i] = S0[i] * A;

            // Density AM=8 Component=XXYZZZZZ
            A = xc_pow2 * yc[i] * zc_pow5;
            phi_tmp[832 + i] = S0[i] * A;

            // Density AM=8 Component=XXZZZZZZ
            A = xc_pow2 * zc_pow6;
            phi_tmp[864 + i] = S0[i] * A;

            // Density AM=8 Component=XYYYYYYY
            A = xc[i] * yc_pow7;
            phi_tmp[896 + i] = S0[i] * A;

            // Density AM=8 Component=XYYYYYYZ
            A = xc[i] * yc_pow6 * zc[i];
            phi_tmp[928 + i] = S0[i] * A;

            // Density AM=8 Component=XYYYYYZZ
            A = xc[i] * yc_pow5 * zc_pow2;
            phi_tmp[960 + i] = S0[i] * A;

            // Density AM=8 Component=XYYYYZZZ
            A = xc[i] * yc_pow4 * zc_pow3;
            phi_tmp[992 + i] = S0[i] * A;

            // Density AM=8 Component=XYYYZZZZ
            A = xc[i] * yc_pow3 * zc_pow4;
            phi_tmp[1024 + i] = S0[i] * A;

            // Density AM=8 Component=XYYZZZZZ
            A = xc[i] * yc_pow2 * zc_pow5;
            phi_tmp[1056 + i] = S0[i] * A;

            // Density AM=8 Component=XYZZZZZZ
            A = xc[i] * yc[i] * zc_pow6;
            phi_tmp[1088 + i] = S0[i] * A;

            // Density AM=8 Component=XZZZZZZZ
            A = xc[i] * zc_pow7;
            phi_tmp[1120 + i] = S0[i] * A;

            // Density AM=8 Component=YYYYYYYY
            phi_tmp[1152 + i] = S0[i] * yc_pow8;

            // Density AM=8 Component=YYYYYYYZ
            A = yc_pow7 * zc[i];
            phi_tmp[1184 + i] = S0[i] * A;

            // Density AM=8 Component=YYYYYYZZ
            A = yc_pow6 * zc_pow2;
            phi_tmp[1216 + i] = S0[i] * A;

            // Density AM=8 Component=YYYYYZZZ
            A = yc_pow5 * zc_pow3;
            phi_tmp[1248 + i] = S0[i] * A;

            // Density AM=8 Component=YYYYZZZZ
            A = yc_pow4 * zc_pow4;
            phi_tmp[1280 + i] = S0[i] * A;

            // Density AM=8 Component=YYYZZZZZ
            A = yc_pow3 * zc_pow5;
            phi_tmp[1312 + i] = S0[i] * A;

            // Density AM=8 Component=YYZZZZZZ
            A = yc_pow2 * zc_pow6;
            phi_tmp[1344 + i] = S0[i] * A;

            // Density AM=8 Component=YZZZZZZZ
            A = yc[i] * zc_pow7;
            phi_tmp[1376 + i] = S0[i] * A;

            // Density AM=8 Component=ZZZZZZZZ
            phi_tmp[1408 + i] = S0[i] * zc_pow8;
        }

        // Copy data back into outer temps
        if (order == GG_SPHERICAL_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_to_spherical_L8(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_SPHERICAL_GAUSSIAN) {
            // Phi, transform data to outer temps
            gg_gaussian_cart_to_spherical_L8(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_CCA) {
            // Phi, transform data to outer temps
            gg_cca_cart_copy_L8(remain, phi_tmp, 32, (phi_out + start), npoints);
        } else if (order == GG_CARTESIAN_MOLDEN) {
            // Phi, transform data to outer temps
            gg_molden_cart_copy_L8(remain, phi_tmp, 32, (phi_out + start), npoints);
        }
    }

    // Free S temporaries
    ALIGNED_FREE(cache_data);
    ALIGNED_FREE(expn1);

    // Free inner temporaries
    ALIGNED_FREE(phi_tmp);
}
