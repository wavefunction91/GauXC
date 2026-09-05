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

void gg_cca_cart_to_spherical_L0(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                 const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                 const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_00 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = cart[i];
    }
}
void gg_cca_cart_to_spherical_sum_L0(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart,
                                     const unsigned long ncart, double* PRAGMA_RESTRICT output,
                                     const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_00 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = cart[i];
        output[i] += tmp * vector[0];
    }
}
void gg_cca_cart_to_spherical_L1(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                 const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                 const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_10 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = cart[ncart + i];
    }

    // R_11c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = cart[2 * ncart + i];
    }
    // R_11s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = cart[i];
    }
}
void gg_cca_cart_to_spherical_sum_L1(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart,
                                     const unsigned long ncart, double* PRAGMA_RESTRICT output,
                                     const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_10 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = cart[ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_11c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = cart[2 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_11s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = cart[i];
        output[i] += tmp * vector[2];
    }
}
void gg_cca_cart_to_spherical_L2(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                 const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                 const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_20 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = 1.7320508075688772 * cart[ncart + i];
    }

    // R_21c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = 1.7320508075688772 * cart[4 * ncart + i];
    }
    // R_21s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = -0.5000000000000000 * cart[i];
        spherical[2 * nspherical + i] += -0.5000000000000000 * cart[3 * ncart + i];
        spherical[2 * nspherical + i] += cart[5 * ncart + i];
    }

    // R_22c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = 1.7320508075688772 * cart[2 * ncart + i];
    }
    // R_22s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = 0.8660254037844386 * cart[i];
        spherical[4 * nspherical + i] += -0.8660254037844386 * cart[3 * ncart + i];
    }
}
void gg_cca_cart_to_spherical_sum_L2(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart,
                                     const unsigned long ncart, double* PRAGMA_RESTRICT output,
                                     const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_20 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.7320508075688772 * cart[ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_21c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.7320508075688772 * cart[4 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_21s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.5000000000000000 * cart[i];
        tmp += -0.5000000000000000 * cart[3 * ncart + i];
        tmp += cart[5 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_22c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.7320508075688772 * cart[2 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_22s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.8660254037844386 * cart[i];
        tmp += -0.8660254037844386 * cart[3 * ncart + i];
        output[i] += tmp * vector[4];
    }
}
void gg_cca_cart_to_spherical_L3(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                 const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                 const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_30 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = 2.3717082451262845 * cart[ncart + i];
        spherical[i] += -0.7905694150420949 * cart[6 * ncart + i];
    }

    // R_31c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = 3.8729833462074170 * cart[4 * ncart + i];
    }
    // R_31s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = -0.6123724356957945 * cart[ncart + i];
        spherical[2 * nspherical + i] += -0.6123724356957945 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] += 2.4494897427831779 * cart[8 * ncart + i];
    }

    // R_32c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = -1.5000000000000000 * cart[2 * ncart + i];
        spherical[3 * nspherical + i] += -1.5000000000000000 * cart[7 * ncart + i];
        spherical[3 * nspherical + i] += cart[9 * ncart + i];
    }
    // R_32s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = -0.6123724356957945 * cart[i];
        spherical[4 * nspherical + i] += -0.6123724356957945 * cart[3 * ncart + i];
        spherical[4 * nspherical + i] += 2.4494897427831779 * cart[5 * ncart + i];
    }

    // R_33c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i] = 1.9364916731037085 * cart[2 * ncart + i];
        spherical[5 * nspherical + i] += -1.9364916731037085 * cart[7 * ncart + i];
    }
    // R_33s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i] = 0.7905694150420949 * cart[i];
        spherical[6 * nspherical + i] += -2.3717082451262845 * cart[3 * ncart + i];
    }
}
void gg_cca_cart_to_spherical_sum_L3(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart,
                                     const unsigned long ncart, double* PRAGMA_RESTRICT output,
                                     const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_30 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.3717082451262845 * cart[ncart + i];
        tmp += -0.7905694150420949 * cart[6 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_31c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 3.8729833462074170 * cart[4 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_31s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.6123724356957945 * cart[ncart + i];
        tmp += -0.6123724356957945 * cart[6 * ncart + i];
        tmp += 2.4494897427831779 * cart[8 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_32c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -1.5000000000000000 * cart[2 * ncart + i];
        tmp += -1.5000000000000000 * cart[7 * ncart + i];
        tmp += cart[9 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_32s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.6123724356957945 * cart[i];
        tmp += -0.6123724356957945 * cart[3 * ncart + i];
        tmp += 2.4494897427831779 * cart[5 * ncart + i];
        output[i] += tmp * vector[4];
    }

    // R_33c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.9364916731037085 * cart[2 * ncart + i];
        tmp += -1.9364916731037085 * cart[7 * ncart + i];
        output[i] += tmp * vector[5];
    }
    // R_33s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.7905694150420949 * cart[i];
        tmp += -2.3717082451262845 * cart[3 * ncart + i];
        output[i] += tmp * vector[6];
    }
}
void gg_cca_cart_to_spherical_L4(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                 const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                 const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_40 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = 2.9580398915498081 * cart[ncart + i];
        spherical[i] += -2.9580398915498081 * cart[6 * ncart + i];
    }

    // R_41c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = 6.2749501990055663 * cart[4 * ncart + i];
        spherical[nspherical + i] += -2.0916500663351889 * cart[11 * ncart + i];
    }
    // R_41s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = -1.1180339887498949 * cart[ncart + i];
        spherical[2 * nspherical + i] += -1.1180339887498949 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] += 6.7082039324993694 * cart[8 * ncart + i];
    }

    // R_42c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = -2.3717082451262845 * cart[4 * ncart + i];
        spherical[3 * nspherical + i] += -2.3717082451262845 * cart[11 * ncart + i];
        spherical[3 * nspherical + i] += 3.1622776601683795 * cart[13 * ncart + i];
    }
    // R_42s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = 0.3750000000000000 * cart[i];
        spherical[4 * nspherical + i] += 0.7500000000000000 * cart[3 * ncart + i];
        spherical[4 * nspherical + i] += 0.3750000000000000 * cart[10 * ncart + i];
        spherical[4 * nspherical + i] += -3.0000000000000000 * cart[5 * ncart + i];
        spherical[4 * nspherical + i] += -3.0000000000000000 * cart[12 * ncart + i];
        spherical[4 * nspherical + i] += cart[14 * ncart + i];
    }

    // R_43c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i] = -2.3717082451262845 * cart[2 * ncart + i];
        spherical[5 * nspherical + i] += -2.3717082451262845 * cart[7 * ncart + i];
        spherical[5 * nspherical + i] += 3.1622776601683795 * cart[9 * ncart + i];
    }
    // R_43s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i] = -0.5590169943749475 * cart[i];
        spherical[6 * nspherical + i] += 0.5590169943749475 * cart[10 * ncart + i];
        spherical[6 * nspherical + i] += 3.3541019662496847 * cart[5 * ncart + i];
        spherical[6 * nspherical + i] += -3.3541019662496847 * cart[12 * ncart + i];
    }

    // R_44c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i] = 2.0916500663351889 * cart[2 * ncart + i];
        spherical[7 * nspherical + i] += -6.2749501990055663 * cart[7 * ncart + i];
    }
    // R_44s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i] = 0.7395099728874520 * cart[i];
        spherical[8 * nspherical + i] += -4.4370598373247123 * cart[3 * ncart + i];
        spherical[8 * nspherical + i] += 0.7395099728874520 * cart[10 * ncart + i];
    }
}
void gg_cca_cart_to_spherical_sum_L4(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart,
                                     const unsigned long ncart, double* PRAGMA_RESTRICT output,
                                     const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_40 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.9580398915498081 * cart[ncart + i];
        tmp += -2.9580398915498081 * cart[6 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_41c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 6.2749501990055663 * cart[4 * ncart + i];
        tmp += -2.0916500663351889 * cart[11 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_41s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -1.1180339887498949 * cart[ncart + i];
        tmp += -1.1180339887498949 * cart[6 * ncart + i];
        tmp += 6.7082039324993694 * cart[8 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_42c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.3717082451262845 * cart[4 * ncart + i];
        tmp += -2.3717082451262845 * cart[11 * ncart + i];
        tmp += 3.1622776601683795 * cart[13 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_42s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.3750000000000000 * cart[i];
        tmp += 0.7500000000000000 * cart[3 * ncart + i];
        tmp += 0.3750000000000000 * cart[10 * ncart + i];
        tmp += -3.0000000000000000 * cart[5 * ncart + i];
        tmp += -3.0000000000000000 * cart[12 * ncart + i];
        tmp += cart[14 * ncart + i];
        output[i] += tmp * vector[4];
    }

    // R_43c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.3717082451262845 * cart[2 * ncart + i];
        tmp += -2.3717082451262845 * cart[7 * ncart + i];
        tmp += 3.1622776601683795 * cart[9 * ncart + i];
        output[i] += tmp * vector[5];
    }
    // R_43s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.5590169943749475 * cart[i];
        tmp += 0.5590169943749475 * cart[10 * ncart + i];
        tmp += 3.3541019662496847 * cart[5 * ncart + i];
        tmp += -3.3541019662496847 * cart[12 * ncart + i];
        output[i] += tmp * vector[6];
    }

    // R_44c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.0916500663351889 * cart[2 * ncart + i];
        tmp += -6.2749501990055663 * cart[7 * ncart + i];
        output[i] += tmp * vector[7];
    }
    // R_44s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.7395099728874520 * cart[i];
        tmp += -4.4370598373247123 * cart[3 * ncart + i];
        tmp += 0.7395099728874520 * cart[10 * ncart + i];
        output[i] += tmp * vector[8];
    }
}
void gg_cca_cart_to_spherical_L5(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                 const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                 const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_50 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = 3.5078038001005702 * cart[ncart + i];
        spherical[i] += -7.0156076002011405 * cart[6 * ncart + i];
        spherical[i] += 0.7015607600201140 * cart[15 * ncart + i];
    }

    // R_51c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = 8.8741196746494246 * cart[4 * ncart + i];
        spherical[nspherical + i] += -8.8741196746494246 * cart[11 * ncart + i];
    }
    // R_51s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = -1.5687375497513916 * cart[ncart + i];
        spherical[2 * nspherical + i] += -1.0458250331675945 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] += 0.5229125165837972 * cart[15 * ncart + i];
        spherical[2 * nspherical + i] += 12.5499003980111326 * cart[8 * ncart + i];
        spherical[2 * nspherical + i] += -4.1833001326703778 * cart[17 * ncart + i];
    }

    // R_52c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = -5.1234753829797990 * cart[4 * ncart + i];
        spherical[3 * nspherical + i] += -5.1234753829797990 * cart[11 * ncart + i];
        spherical[3 * nspherical + i] += 10.2469507659595980 * cart[13 * ncart + i];
    }
    // R_52s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = 0.4841229182759271 * cart[ncart + i];
        spherical[4 * nspherical + i] += 0.9682458365518543 * cart[6 * ncart + i];
        spherical[4 * nspherical + i] += 0.4841229182759271 * cart[15 * ncart + i];
        spherical[4 * nspherical + i] += -5.8094750193111251 * cart[8 * ncart + i];
        spherical[4 * nspherical + i] += -5.8094750193111251 * cart[17 * ncart + i];
        spherical[4 * nspherical + i] += 3.8729833462074170 * cart[19 * ncart + i];
    }

    // R_53c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i] = 1.8750000000000000 * cart[2 * ncart + i];
        spherical[5 * nspherical + i] += 3.7500000000000000 * cart[7 * ncart + i];
        spherical[5 * nspherical + i] += 1.8750000000000000 * cart[16 * ncart + i];
        spherical[5 * nspherical + i] += -5.0000000000000000 * cart[9 * ncart + i];
        spherical[5 * nspherical + i] += -5.0000000000000000 * cart[18 * ncart + i];
        spherical[5 * nspherical + i] += cart[20 * ncart + i];
    }
    // R_53s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i] = 0.4841229182759271 * cart[i];
        spherical[6 * nspherical + i] += 0.9682458365518543 * cart[3 * ncart + i];
        spherical[6 * nspherical + i] += 0.4841229182759271 * cart[10 * ncart + i];
        spherical[6 * nspherical + i] += -5.8094750193111251 * cart[5 * ncart + i];
        spherical[6 * nspherical + i] += -5.8094750193111251 * cart[12 * ncart + i];
        spherical[6 * nspherical + i] += 3.8729833462074170 * cart[14 * ncart + i];
    }

    // R_54c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i] = -2.5617376914898995 * cart[2 * ncart + i];
        spherical[7 * nspherical + i] += 2.5617376914898995 * cart[16 * ncart + i];
        spherical[7 * nspherical + i] += 5.1234753829797990 * cart[9 * ncart + i];
        spherical[7 * nspherical + i] += -5.1234753829797990 * cart[18 * ncart + i];
    }
    // R_54s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i] = -0.5229125165837972 * cart[i];
        spherical[8 * nspherical + i] += 1.0458250331675945 * cart[3 * ncart + i];
        spherical[8 * nspherical + i] += 1.5687375497513916 * cart[10 * ncart + i];
        spherical[8 * nspherical + i] += 4.1833001326703778 * cart[5 * ncart + i];
        spherical[8 * nspherical + i] += -12.5499003980111326 * cart[12 * ncart + i];
    }

    // R_55c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[9 * nspherical + i] = 2.2185299186623562 * cart[2 * ncart + i];
        spherical[9 * nspherical + i] += -13.3111795119741370 * cart[7 * ncart + i];
        spherical[9 * nspherical + i] += 2.2185299186623562 * cart[16 * ncart + i];
    }
    // R_55s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[10 * nspherical + i] = 0.7015607600201140 * cart[i];
        spherical[10 * nspherical + i] += -7.0156076002011405 * cart[3 * ncart + i];
        spherical[10 * nspherical + i] += 3.5078038001005702 * cart[10 * ncart + i];
    }
}
void gg_cca_cart_to_spherical_sum_L5(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart,
                                     const unsigned long ncart, double* PRAGMA_RESTRICT output,
                                     const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_50 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 3.5078038001005702 * cart[ncart + i];
        tmp += -7.0156076002011405 * cart[6 * ncart + i];
        tmp += 0.7015607600201140 * cart[15 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_51c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 8.8741196746494246 * cart[4 * ncart + i];
        tmp += -8.8741196746494246 * cart[11 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_51s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -1.5687375497513916 * cart[ncart + i];
        tmp += -1.0458250331675945 * cart[6 * ncart + i];
        tmp += 0.5229125165837972 * cart[15 * ncart + i];
        tmp += 12.5499003980111326 * cart[8 * ncart + i];
        tmp += -4.1833001326703778 * cart[17 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_52c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -5.1234753829797990 * cart[4 * ncart + i];
        tmp += -5.1234753829797990 * cart[11 * ncart + i];
        tmp += 10.2469507659595980 * cart[13 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_52s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.4841229182759271 * cart[ncart + i];
        tmp += 0.9682458365518543 * cart[6 * ncart + i];
        tmp += 0.4841229182759271 * cart[15 * ncart + i];
        tmp += -5.8094750193111251 * cart[8 * ncart + i];
        tmp += -5.8094750193111251 * cart[17 * ncart + i];
        tmp += 3.8729833462074170 * cart[19 * ncart + i];
        output[i] += tmp * vector[4];
    }

    // R_53c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.8750000000000000 * cart[2 * ncart + i];
        tmp += 3.7500000000000000 * cart[7 * ncart + i];
        tmp += 1.8750000000000000 * cart[16 * ncart + i];
        tmp += -5.0000000000000000 * cart[9 * ncart + i];
        tmp += -5.0000000000000000 * cart[18 * ncart + i];
        tmp += cart[20 * ncart + i];
        output[i] += tmp * vector[5];
    }
    // R_53s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.4841229182759271 * cart[i];
        tmp += 0.9682458365518543 * cart[3 * ncart + i];
        tmp += 0.4841229182759271 * cart[10 * ncart + i];
        tmp += -5.8094750193111251 * cart[5 * ncart + i];
        tmp += -5.8094750193111251 * cart[12 * ncart + i];
        tmp += 3.8729833462074170 * cart[14 * ncart + i];
        output[i] += tmp * vector[6];
    }

    // R_54c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.5617376914898995 * cart[2 * ncart + i];
        tmp += 2.5617376914898995 * cart[16 * ncart + i];
        tmp += 5.1234753829797990 * cart[9 * ncart + i];
        tmp += -5.1234753829797990 * cart[18 * ncart + i];
        output[i] += tmp * vector[7];
    }
    // R_54s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.5229125165837972 * cart[i];
        tmp += 1.0458250331675945 * cart[3 * ncart + i];
        tmp += 1.5687375497513916 * cart[10 * ncart + i];
        tmp += 4.1833001326703778 * cart[5 * ncart + i];
        tmp += -12.5499003980111326 * cart[12 * ncart + i];
        output[i] += tmp * vector[8];
    }

    // R_55c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.2185299186623562 * cart[2 * ncart + i];
        tmp += -13.3111795119741370 * cart[7 * ncart + i];
        tmp += 2.2185299186623562 * cart[16 * ncart + i];
        output[i] += tmp * vector[9];
    }
    // R_55s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.7015607600201140 * cart[i];
        tmp += -7.0156076002011405 * cart[3 * ncart + i];
        tmp += 3.5078038001005702 * cart[10 * ncart + i];
        output[i] += tmp * vector[10];
    }
}
void gg_cca_cart_to_spherical_L6(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                 const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                 const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_60 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = 4.0301597362883772 * cart[ncart + i];
        spherical[i] += -13.4338657876279228 * cart[6 * ncart + i];
        spherical[i] += 4.0301597362883772 * cart[15 * ncart + i];
    }

    // R_61c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = 11.6340690431164280 * cart[4 * ncart + i];
        spherical[nspherical + i] += -23.2681380862328560 * cart[11 * ncart + i];
        spherical[nspherical + i] += 2.3268138086232857 * cart[22 * ncart + i];
    }
    // R_61s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = -1.9843134832984430 * cart[ncart + i];
        spherical[2 * nspherical + i] += 1.9843134832984430 * cart[15 * ncart + i];
        spherical[2 * nspherical + i] += 19.8431348329844290 * cart[8 * ncart + i];
        spherical[2 * nspherical + i] += -19.8431348329844290 * cart[17 * ncart + i];
    }

    // R_62c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = -8.1513994197315593 * cart[4 * ncart + i];
        spherical[3 * nspherical + i] += -5.4342662798210393 * cart[11 * ncart + i];
        spherical[3 * nspherical + i] += 2.7171331399105196 * cart[22 * ncart + i];
        spherical[3 * nspherical + i] += 21.7370651192841571 * cart[13 * ncart + i];
        spherical[3 * nspherical + i] += -7.2456883730947190 * cart[24 * ncart + i];
    }
    // R_62s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = 0.9057110466368399 * cart[ncart + i];
        spherical[4 * nspherical + i] += 1.8114220932736798 * cart[6 * ncart + i];
        spherical[4 * nspherical + i] += 0.9057110466368399 * cart[15 * ncart + i];
        spherical[4 * nspherical + i] += -14.4913767461894381 * cart[8 * ncart + i];
        spherical[4 * nspherical + i] += -14.4913767461894381 * cart[17 * ncart + i];
        spherical[4 * nspherical + i] += 14.4913767461894381 * cart[19 * ncart + i];
    }

    // R_63c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i] = 2.8641098093473998 * cart[4 * ncart + i];
        spherical[5 * nspherical + i] += 5.7282196186947996 * cart[11 * ncart + i];
        spherical[5 * nspherical + i] += 2.8641098093473998 * cart[22 * ncart + i];
        spherical[5 * nspherical + i] += -11.4564392373895991 * cart[13 * ncart + i];
        spherical[5 * nspherical + i] += -11.4564392373895991 * cart[24 * ncart + i];
        spherical[5 * nspherical + i] += 4.5825756949558398 * cart[26 * ncart + i];
    }
    // R_63s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i] = -0.3125000000000000 * cart[i];
        spherical[6 * nspherical + i] += -0.9375000000000000 * cart[3 * ncart + i];
        spherical[6 * nspherical + i] += -0.9375000000000000 * cart[10 * ncart + i];
        spherical[6 * nspherical + i] += -0.3125000000000000 * cart[21 * ncart + i];
        spherical[6 * nspherical + i] += 5.6250000000000000 * cart[5 * ncart + i];
        spherical[6 * nspherical + i] += 11.2500000000000000 * cart[12 * ncart + i];
        spherical[6 * nspherical + i] += 5.6250000000000000 * cart[23 * ncart + i];
        spherical[6 * nspherical + i] += -7.5000000000000000 * cart[14 * ncart + i];
        spherical[6 * nspherical + i] += -7.5000000000000000 * cart[25 * ncart + i];
        spherical[6 * nspherical + i] += cart[27 * ncart + i];
    }

    // R_64c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i] = 2.8641098093473998 * cart[2 * ncart + i];
        spherical[7 * nspherical + i] += 5.7282196186947996 * cart[7 * ncart + i];
        spherical[7 * nspherical + i] += 2.8641098093473998 * cart[16 * ncart + i];
        spherical[7 * nspherical + i] += -11.4564392373895991 * cart[9 * ncart + i];
        spherical[7 * nspherical + i] += -11.4564392373895991 * cart[18 * ncart + i];
        spherical[7 * nspherical + i] += 4.5825756949558398 * cart[20 * ncart + i];
    }
    // R_64s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i] = 0.4528555233184199 * cart[i];
        spherical[8 * nspherical + i] += 0.4528555233184199 * cart[3 * ncart + i];
        spherical[8 * nspherical + i] += -0.4528555233184199 * cart[10 * ncart + i];
        spherical[8 * nspherical + i] += -0.4528555233184199 * cart[21 * ncart + i];
        spherical[8 * nspherical + i] += -7.2456883730947190 * cart[5 * ncart + i];
        spherical[8 * nspherical + i] += 7.2456883730947190 * cart[23 * ncart + i];
        spherical[8 * nspherical + i] += 7.2456883730947190 * cart[14 * ncart + i];
        spherical[8 * nspherical + i] += -7.2456883730947190 * cart[25 * ncart + i];
    }

    // R_65c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[9 * nspherical + i] = -2.7171331399105196 * cart[2 * ncart + i];
        spherical[9 * nspherical + i] += 5.4342662798210393 * cart[7 * ncart + i];
        spherical[9 * nspherical + i] += 8.1513994197315593 * cart[16 * ncart + i];
        spherical[9 * nspherical + i] += 7.2456883730947190 * cart[9 * ncart + i];
        spherical[9 * nspherical + i] += -21.7370651192841571 * cart[18 * ncart + i];
    }
    // R_65s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[10 * nspherical + i] = -0.4960783708246108 * cart[i];
        spherical[10 * nspherical + i] += 2.4803918541230536 * cart[3 * ncart + i];
        spherical[10 * nspherical + i] += 2.4803918541230536 * cart[10 * ncart + i];
        spherical[10 * nspherical + i] += -0.4960783708246108 * cart[21 * ncart + i];
        spherical[10 * nspherical + i] += 4.9607837082461073 * cart[5 * ncart + i];
        spherical[10 * nspherical + i] += -29.7647022494766453 * cart[12 * ncart + i];
        spherical[10 * nspherical + i] += 4.9607837082461073 * cart[23 * ncart + i];
    }

    // R_66c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[11 * nspherical + i] = 2.3268138086232857 * cart[2 * ncart + i];
        spherical[11 * nspherical + i] += -23.2681380862328560 * cart[7 * ncart + i];
        spherical[11 * nspherical + i] += 11.6340690431164280 * cart[16 * ncart + i];
    }
    // R_66s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[12 * nspherical + i] = 0.6716932893813962 * cart[i];
        spherical[12 * nspherical + i] += -10.0753993407209421 * cart[3 * ncart + i];
        spherical[12 * nspherical + i] += 10.0753993407209421 * cart[10 * ncart + i];
        spherical[12 * nspherical + i] += -0.6716932893813962 * cart[21 * ncart + i];
    }
}
void gg_cca_cart_to_spherical_sum_L6(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart,
                                     const unsigned long ncart, double* PRAGMA_RESTRICT output,
                                     const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_60 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 4.0301597362883772 * cart[ncart + i];
        tmp += -13.4338657876279228 * cart[6 * ncart + i];
        tmp += 4.0301597362883772 * cart[15 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_61c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 11.6340690431164280 * cart[4 * ncart + i];
        tmp += -23.2681380862328560 * cart[11 * ncart + i];
        tmp += 2.3268138086232857 * cart[22 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_61s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -1.9843134832984430 * cart[ncart + i];
        tmp += 1.9843134832984430 * cart[15 * ncart + i];
        tmp += 19.8431348329844290 * cart[8 * ncart + i];
        tmp += -19.8431348329844290 * cart[17 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_62c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -8.1513994197315593 * cart[4 * ncart + i];
        tmp += -5.4342662798210393 * cart[11 * ncart + i];
        tmp += 2.7171331399105196 * cart[22 * ncart + i];
        tmp += 21.7370651192841571 * cart[13 * ncart + i];
        tmp += -7.2456883730947190 * cart[24 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_62s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.9057110466368399 * cart[ncart + i];
        tmp += 1.8114220932736798 * cart[6 * ncart + i];
        tmp += 0.9057110466368399 * cart[15 * ncart + i];
        tmp += -14.4913767461894381 * cart[8 * ncart + i];
        tmp += -14.4913767461894381 * cart[17 * ncart + i];
        tmp += 14.4913767461894381 * cart[19 * ncart + i];
        output[i] += tmp * vector[4];
    }

    // R_63c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.8641098093473998 * cart[4 * ncart + i];
        tmp += 5.7282196186947996 * cart[11 * ncart + i];
        tmp += 2.8641098093473998 * cart[22 * ncart + i];
        tmp += -11.4564392373895991 * cart[13 * ncart + i];
        tmp += -11.4564392373895991 * cart[24 * ncart + i];
        tmp += 4.5825756949558398 * cart[26 * ncart + i];
        output[i] += tmp * vector[5];
    }
    // R_63s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.3125000000000000 * cart[i];
        tmp += -0.9375000000000000 * cart[3 * ncart + i];
        tmp += -0.9375000000000000 * cart[10 * ncart + i];
        tmp += -0.3125000000000000 * cart[21 * ncart + i];
        tmp += 5.6250000000000000 * cart[5 * ncart + i];
        tmp += 11.2500000000000000 * cart[12 * ncart + i];
        tmp += 5.6250000000000000 * cart[23 * ncart + i];
        tmp += -7.5000000000000000 * cart[14 * ncart + i];
        tmp += -7.5000000000000000 * cart[25 * ncart + i];
        tmp += cart[27 * ncart + i];
        output[i] += tmp * vector[6];
    }

    // R_64c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.8641098093473998 * cart[2 * ncart + i];
        tmp += 5.7282196186947996 * cart[7 * ncart + i];
        tmp += 2.8641098093473998 * cart[16 * ncart + i];
        tmp += -11.4564392373895991 * cart[9 * ncart + i];
        tmp += -11.4564392373895991 * cart[18 * ncart + i];
        tmp += 4.5825756949558398 * cart[20 * ncart + i];
        output[i] += tmp * vector[7];
    }
    // R_64s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.4528555233184199 * cart[i];
        tmp += 0.4528555233184199 * cart[3 * ncart + i];
        tmp += -0.4528555233184199 * cart[10 * ncart + i];
        tmp += -0.4528555233184199 * cart[21 * ncart + i];
        tmp += -7.2456883730947190 * cart[5 * ncart + i];
        tmp += 7.2456883730947190 * cart[23 * ncart + i];
        tmp += 7.2456883730947190 * cart[14 * ncart + i];
        tmp += -7.2456883730947190 * cart[25 * ncart + i];
        output[i] += tmp * vector[8];
    }

    // R_65c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.7171331399105196 * cart[2 * ncart + i];
        tmp += 5.4342662798210393 * cart[7 * ncart + i];
        tmp += 8.1513994197315593 * cart[16 * ncart + i];
        tmp += 7.2456883730947190 * cart[9 * ncart + i];
        tmp += -21.7370651192841571 * cart[18 * ncart + i];
        output[i] += tmp * vector[9];
    }
    // R_65s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.4960783708246108 * cart[i];
        tmp += 2.4803918541230536 * cart[3 * ncart + i];
        tmp += 2.4803918541230536 * cart[10 * ncart + i];
        tmp += -0.4960783708246108 * cart[21 * ncart + i];
        tmp += 4.9607837082461073 * cart[5 * ncart + i];
        tmp += -29.7647022494766453 * cart[12 * ncart + i];
        tmp += 4.9607837082461073 * cart[23 * ncart + i];
        output[i] += tmp * vector[10];
    }

    // R_66c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.3268138086232857 * cart[2 * ncart + i];
        tmp += -23.2681380862328560 * cart[7 * ncart + i];
        tmp += 11.6340690431164280 * cart[16 * ncart + i];
        output[i] += tmp * vector[11];
    }
    // R_66s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.6716932893813962 * cart[i];
        tmp += -10.0753993407209421 * cart[3 * ncart + i];
        tmp += 10.0753993407209421 * cart[10 * ncart + i];
        tmp += -0.6716932893813962 * cart[21 * ncart + i];
        output[i] += tmp * vector[12];
    }
}
void gg_cca_cart_to_spherical_L7(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                 const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                 const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_70 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = 4.5308189450142455 * cart[ncart + i];
        spherical[i] += -22.6540947250712286 * cart[6 * ncart + i];
        spherical[i] += 13.5924568350427357 * cart[15 * ncart + i];
        spherical[i] += -0.6472598492877494 * cart[28 * ncart + i];
    }

    // R_71c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = 14.5309475774981713 * cart[4 * ncart + i];
        spherical[nspherical + i] += -48.4364919249939092 * cart[11 * ncart + i];
        spherical[nspherical + i] += 14.5309475774981713 * cart[22 * ncart + i];
    }
    // R_71s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = -2.3747943989954163 * cart[ncart + i];
        spherical[2 * nspherical + i] += 2.3747943989954163 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] += 4.2746299181917493 * cart[15 * ncart + i];
        spherical[2 * nspherical + i] += -0.4749588797990832 * cart[28 * ncart + i];
        spherical[2 * nspherical + i] += 28.4975327879449942 * cart[8 * ncart + i];
        spherical[2 * nspherical + i] += -56.9950655758899885 * cart[17 * ncart + i];
        spherical[2 * nspherical + i] += 5.6995065575889985 * cart[30 * ncart + i];
    }

    // R_72c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = -11.3990131151779970 * cart[4 * ncart + i];
        spherical[3 * nspherical + i] += 11.3990131151779970 * cart[22 * ncart + i];
        spherical[3 * nspherical + i] += 37.9967103839266613 * cart[13 * ncart + i];
        spherical[3 * nspherical + i] += -37.9967103839266613 * cart[24 * ncart + i];
    }
    // R_72s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = 1.2888494142063300 * cart[ncart + i];
        spherical[4 * nspherical + i] += 2.1480823570105501 * cart[6 * ncart + i];
        spherical[4 * nspherical + i] += 0.4296164714021100 * cart[15 * ncart + i];
        spherical[4 * nspherical + i] += -0.4296164714021100 * cart[28 * ncart + i];
        spherical[4 * nspherical + i] += -25.7769882841265989 * cart[8 * ncart + i];
        spherical[4 * nspherical + i] += -17.1846588560844005 * cart[17 * ncart + i];
        spherical[4 * nspherical + i] += 8.5923294280422002 * cart[30 * ncart + i];
        spherical[4 * nspherical + i] += 34.3693177121688009 * cart[19 * ncart + i];
        spherical[4 * nspherical + i] += -11.4564392373895991 * cart[32 * ncart + i];
    }

    // R_73c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i] = 6.0756944047573693 * cart[4 * ncart + i];
        spherical[5 * nspherical + i] += 12.1513888095147387 * cart[11 * ncart + i];
        spherical[5 * nspherical + i] += 6.0756944047573693 * cart[22 * ncart + i];
        spherical[5 * nspherical + i] += -32.4037034920392983 * cart[13 * ncart + i];
        spherical[5 * nspherical + i] += -32.4037034920392983 * cart[24 * ncart + i];
        spherical[5 * nspherical + i] += 19.4422220952235811 * cart[26 * ncart + i];
    }
    // R_73s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i] = -0.4133986423538423 * cart[ncart + i];
        spherical[6 * nspherical + i] += -1.2401959270615268 * cart[6 * ncart + i];
        spherical[6 * nspherical + i] += -1.2401959270615268 * cart[15 * ncart + i];
        spherical[6 * nspherical + i] += -0.4133986423538423 * cart[28 * ncart + i];
        spherical[6 * nspherical + i] += 9.9215674164922145 * cart[8 * ncart + i];
        spherical[6 * nspherical + i] += 19.8431348329844290 * cart[17 * ncart + i];
        spherical[6 * nspherical + i] += 9.9215674164922145 * cart[30 * ncart + i];
        spherical[6 * nspherical + i] += -19.8431348329844290 * cart[19 * ncart + i];
        spherical[6 * nspherical + i] += -19.8431348329844290 * cart[32 * ncart + i];
        spherical[6 * nspherical + i] += 5.2915026221291814 * cart[34 * ncart + i];
    }

    // R_74c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i] = -2.1875000000000000 * cart[2 * ncart + i];
        spherical[7 * nspherical + i] += -6.5625000000000000 * cart[7 * ncart + i];
        spherical[7 * nspherical + i] += -6.5625000000000000 * cart[16 * ncart + i];
        spherical[7 * nspherical + i] += -2.1875000000000000 * cart[29 * ncart + i];
        spherical[7 * nspherical + i] += 13.1250000000000000 * cart[9 * ncart + i];
        spherical[7 * nspherical + i] += 26.2500000000000000 * cart[18 * ncart + i];
        spherical[7 * nspherical + i] += 13.1250000000000000 * cart[31 * ncart + i];
        spherical[7 * nspherical + i] += -10.5000000000000000 * cart[20 * ncart + i];
        spherical[7 * nspherical + i] += -10.5000000000000000 * cart[33 * ncart + i];
        spherical[7 * nspherical + i] += cart[35 * ncart + i];
    }
    // R_74s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i] = -0.4133986423538423 * cart[i];
        spherical[8 * nspherical + i] += -1.2401959270615268 * cart[3 * ncart + i];
        spherical[8 * nspherical + i] += -1.2401959270615268 * cart[10 * ncart + i];
        spherical[8 * nspherical + i] += -0.4133986423538423 * cart[21 * ncart + i];
        spherical[8 * nspherical + i] += 9.9215674164922145 * cart[5 * ncart + i];
        spherical[8 * nspherical + i] += 19.8431348329844290 * cart[12 * ncart + i];
        spherical[8 * nspherical + i] += 9.9215674164922145 * cart[23 * ncart + i];
        spherical[8 * nspherical + i] += -19.8431348329844290 * cart[14 * ncart + i];
        spherical[8 * nspherical + i] += -19.8431348329844290 * cart[25 * ncart + i];
        spherical[8 * nspherical + i] += 5.2915026221291814 * cart[27 * ncart + i];
    }

    // R_75c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[9 * nspherical + i] = 3.0378472023786847 * cart[2 * ncart + i];
        spherical[9 * nspherical + i] += 3.0378472023786847 * cart[7 * ncart + i];
        spherical[9 * nspherical + i] += -3.0378472023786847 * cart[16 * ncart + i];
        spherical[9 * nspherical + i] += -3.0378472023786847 * cart[29 * ncart + i];
        spherical[9 * nspherical + i] += -16.2018517460196492 * cart[9 * ncart + i];
        spherical[9 * nspherical + i] += 16.2018517460196492 * cart[31 * ncart + i];
        spherical[9 * nspherical + i] += 9.7211110476117906 * cart[20 * ncart + i];
        spherical[9 * nspherical + i] += -9.7211110476117906 * cart[33 * ncart + i];
    }
    // R_75s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[10 * nspherical + i] = 0.4296164714021100 * cart[i];
        spherical[10 * nspherical + i] += -0.4296164714021100 * cart[3 * ncart + i];
        spherical[10 * nspherical + i] += -2.1480823570105501 * cart[10 * ncart + i];
        spherical[10 * nspherical + i] += -1.2888494142063300 * cart[21 * ncart + i];
        spherical[10 * nspherical + i] += -8.5923294280422002 * cart[5 * ncart + i];
        spherical[10 * nspherical + i] += 17.1846588560844005 * cart[12 * ncart + i];
        spherical[10 * nspherical + i] += 25.7769882841265989 * cart[23 * ncart + i];
        spherical[10 * nspherical + i] += 11.4564392373895991 * cart[14 * ncart + i];
        spherical[10 * nspherical + i] += -34.3693177121688009 * cart[25 * ncart + i];
    }

    // R_76c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[11 * nspherical + i] = -2.8497532787944992 * cart[2 * ncart + i];
        spherical[11 * nspherical + i] += 14.2487663939724971 * cart[7 * ncart + i];
        spherical[11 * nspherical + i] += 14.2487663939724971 * cart[16 * ncart + i];
        spherical[11 * nspherical + i] += -2.8497532787944992 * cart[29 * ncart + i];
        spherical[11 * nspherical + i] += 9.4991775959816653 * cart[9 * ncart + i];
        spherical[11 * nspherical + i] += -56.9950655758899885 * cart[18 * ncart + i];
        spherical[11 * nspherical + i] += 9.4991775959816653 * cart[31 * ncart + i];
    }
    // R_76s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[12 * nspherical + i] = -0.4749588797990832 * cart[i];
        spherical[12 * nspherical + i] += 4.2746299181917493 * cart[3 * ncart + i];
        spherical[12 * nspherical + i] += 2.3747943989954163 * cart[10 * ncart + i];
        spherical[12 * nspherical + i] += -2.3747943989954163 * cart[21 * ncart + i];
        spherical[12 * nspherical + i] += 5.6995065575889985 * cart[5 * ncart + i];
        spherical[12 * nspherical + i] += -56.9950655758899885 * cart[12 * ncart + i];
        spherical[12 * nspherical + i] += 28.4975327879449942 * cart[23 * ncart + i];
    }

    // R_77c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[13 * nspherical + i] = 2.4218245962496954 * cart[2 * ncart + i];
        spherical[13 * nspherical + i] += -36.3273689437454337 * cart[7 * ncart + i];
        spherical[13 * nspherical + i] += 36.3273689437454337 * cart[16 * ncart + i];
        spherical[13 * nspherical + i] += -2.4218245962496954 * cart[29 * ncart + i];
    }
    // R_77s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[14 * nspherical + i] = 0.6472598492877494 * cart[i];
        spherical[14 * nspherical + i] += -13.5924568350427357 * cart[3 * ncart + i];
        spherical[14 * nspherical + i] += 22.6540947250712286 * cart[10 * ncart + i];
        spherical[14 * nspherical + i] += -4.5308189450142455 * cart[21 * ncart + i];
    }
}
void gg_cca_cart_to_spherical_sum_L7(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart,
                                     const unsigned long ncart, double* PRAGMA_RESTRICT output,
                                     const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_70 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 4.5308189450142455 * cart[ncart + i];
        tmp += -22.6540947250712286 * cart[6 * ncart + i];
        tmp += 13.5924568350427357 * cart[15 * ncart + i];
        tmp += -0.6472598492877494 * cart[28 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_71c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 14.5309475774981713 * cart[4 * ncart + i];
        tmp += -48.4364919249939092 * cart[11 * ncart + i];
        tmp += 14.5309475774981713 * cart[22 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_71s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.3747943989954163 * cart[ncart + i];
        tmp += 2.3747943989954163 * cart[6 * ncart + i];
        tmp += 4.2746299181917493 * cart[15 * ncart + i];
        tmp += -0.4749588797990832 * cart[28 * ncart + i];
        tmp += 28.4975327879449942 * cart[8 * ncart + i];
        tmp += -56.9950655758899885 * cart[17 * ncart + i];
        tmp += 5.6995065575889985 * cart[30 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_72c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -11.3990131151779970 * cart[4 * ncart + i];
        tmp += 11.3990131151779970 * cart[22 * ncart + i];
        tmp += 37.9967103839266613 * cart[13 * ncart + i];
        tmp += -37.9967103839266613 * cart[24 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_72s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.2888494142063300 * cart[ncart + i];
        tmp += 2.1480823570105501 * cart[6 * ncart + i];
        tmp += 0.4296164714021100 * cart[15 * ncart + i];
        tmp += -0.4296164714021100 * cart[28 * ncart + i];
        tmp += -25.7769882841265989 * cart[8 * ncart + i];
        tmp += -17.1846588560844005 * cart[17 * ncart + i];
        tmp += 8.5923294280422002 * cart[30 * ncart + i];
        tmp += 34.3693177121688009 * cart[19 * ncart + i];
        tmp += -11.4564392373895991 * cart[32 * ncart + i];
        output[i] += tmp * vector[4];
    }

    // R_73c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 6.0756944047573693 * cart[4 * ncart + i];
        tmp += 12.1513888095147387 * cart[11 * ncart + i];
        tmp += 6.0756944047573693 * cart[22 * ncart + i];
        tmp += -32.4037034920392983 * cart[13 * ncart + i];
        tmp += -32.4037034920392983 * cart[24 * ncart + i];
        tmp += 19.4422220952235811 * cart[26 * ncart + i];
        output[i] += tmp * vector[5];
    }
    // R_73s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.4133986423538423 * cart[ncart + i];
        tmp += -1.2401959270615268 * cart[6 * ncart + i];
        tmp += -1.2401959270615268 * cart[15 * ncart + i];
        tmp += -0.4133986423538423 * cart[28 * ncart + i];
        tmp += 9.9215674164922145 * cart[8 * ncart + i];
        tmp += 19.8431348329844290 * cart[17 * ncart + i];
        tmp += 9.9215674164922145 * cart[30 * ncart + i];
        tmp += -19.8431348329844290 * cart[19 * ncart + i];
        tmp += -19.8431348329844290 * cart[32 * ncart + i];
        tmp += 5.2915026221291814 * cart[34 * ncart + i];
        output[i] += tmp * vector[6];
    }

    // R_74c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.1875000000000000 * cart[2 * ncart + i];
        tmp += -6.5625000000000000 * cart[7 * ncart + i];
        tmp += -6.5625000000000000 * cart[16 * ncart + i];
        tmp += -2.1875000000000000 * cart[29 * ncart + i];
        tmp += 13.1250000000000000 * cart[9 * ncart + i];
        tmp += 26.2500000000000000 * cart[18 * ncart + i];
        tmp += 13.1250000000000000 * cart[31 * ncart + i];
        tmp += -10.5000000000000000 * cart[20 * ncart + i];
        tmp += -10.5000000000000000 * cart[33 * ncart + i];
        tmp += cart[35 * ncart + i];
        output[i] += tmp * vector[7];
    }
    // R_74s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.4133986423538423 * cart[i];
        tmp += -1.2401959270615268 * cart[3 * ncart + i];
        tmp += -1.2401959270615268 * cart[10 * ncart + i];
        tmp += -0.4133986423538423 * cart[21 * ncart + i];
        tmp += 9.9215674164922145 * cart[5 * ncart + i];
        tmp += 19.8431348329844290 * cart[12 * ncart + i];
        tmp += 9.9215674164922145 * cart[23 * ncart + i];
        tmp += -19.8431348329844290 * cart[14 * ncart + i];
        tmp += -19.8431348329844290 * cart[25 * ncart + i];
        tmp += 5.2915026221291814 * cart[27 * ncart + i];
        output[i] += tmp * vector[8];
    }

    // R_75c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 3.0378472023786847 * cart[2 * ncart + i];
        tmp += 3.0378472023786847 * cart[7 * ncart + i];
        tmp += -3.0378472023786847 * cart[16 * ncart + i];
        tmp += -3.0378472023786847 * cart[29 * ncart + i];
        tmp += -16.2018517460196492 * cart[9 * ncart + i];
        tmp += 16.2018517460196492 * cart[31 * ncart + i];
        tmp += 9.7211110476117906 * cart[20 * ncart + i];
        tmp += -9.7211110476117906 * cart[33 * ncart + i];
        output[i] += tmp * vector[9];
    }
    // R_75s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.4296164714021100 * cart[i];
        tmp += -0.4296164714021100 * cart[3 * ncart + i];
        tmp += -2.1480823570105501 * cart[10 * ncart + i];
        tmp += -1.2888494142063300 * cart[21 * ncart + i];
        tmp += -8.5923294280422002 * cart[5 * ncart + i];
        tmp += 17.1846588560844005 * cart[12 * ncart + i];
        tmp += 25.7769882841265989 * cart[23 * ncart + i];
        tmp += 11.4564392373895991 * cart[14 * ncart + i];
        tmp += -34.3693177121688009 * cart[25 * ncart + i];
        output[i] += tmp * vector[10];
    }

    // R_76c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.8497532787944992 * cart[2 * ncart + i];
        tmp += 14.2487663939724971 * cart[7 * ncart + i];
        tmp += 14.2487663939724971 * cart[16 * ncart + i];
        tmp += -2.8497532787944992 * cart[29 * ncart + i];
        tmp += 9.4991775959816653 * cart[9 * ncart + i];
        tmp += -56.9950655758899885 * cart[18 * ncart + i];
        tmp += 9.4991775959816653 * cart[31 * ncart + i];
        output[i] += tmp * vector[11];
    }
    // R_76s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.4749588797990832 * cart[i];
        tmp += 4.2746299181917493 * cart[3 * ncart + i];
        tmp += 2.3747943989954163 * cart[10 * ncart + i];
        tmp += -2.3747943989954163 * cart[21 * ncart + i];
        tmp += 5.6995065575889985 * cart[5 * ncart + i];
        tmp += -56.9950655758899885 * cart[12 * ncart + i];
        tmp += 28.4975327879449942 * cart[23 * ncart + i];
        output[i] += tmp * vector[12];
    }

    // R_77c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.4218245962496954 * cart[2 * ncart + i];
        tmp += -36.3273689437454337 * cart[7 * ncart + i];
        tmp += 36.3273689437454337 * cart[16 * ncart + i];
        tmp += -2.4218245962496954 * cart[29 * ncart + i];
        output[i] += tmp * vector[13];
    }
    // R_77s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.6472598492877494 * cart[i];
        tmp += -13.5924568350427357 * cart[3 * ncart + i];
        tmp += 22.6540947250712286 * cart[10 * ncart + i];
        tmp += -4.5308189450142455 * cart[21 * ncart + i];
        output[i] += tmp * vector[14];
    }
}
void gg_cca_cart_to_spherical_L8(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                 const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                 const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_80 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = 5.0136532339203512 * cart[ncart + i];
        spherical[i] += -35.0955726374424586 * cart[6 * ncart + i];
        spherical[i] += 35.0955726374424586 * cart[15 * ncart + i];
        spherical[i] += -5.0136532339203512 * cart[28 * ncart + i];
    }

    // R_81c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = 17.5477863187212293 * cart[4 * ncart + i];
        spherical[nspherical + i] += -87.7389315936061536 * cart[11 * ncart + i];
        spherical[nspherical + i] += 52.6433589561636950 * cart[22 * ncart + i];
        spherical[nspherical + i] += -2.5068266169601756 * cart[37 * ncart + i];
    }
    // R_81s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = -2.7460909717269018 * cart[ncart + i];
        spherical[2 * nspherical + i] += 6.4075456006961042 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] += 6.4075456006961042 * cart[15 * ncart + i];
        spherical[2 * nspherical + i] += -2.7460909717269018 * cart[28 * ncart + i];
        spherical[2 * nspherical + i] += 38.4452736041766272 * cart[8 * ncart + i];
        spherical[2 * nspherical + i] += -128.1509120139220954 * cart[17 * ncart + i];
        spherical[2 * nspherical + i] += 38.4452736041766272 * cart[30 * ncart + i];
    }

    // R_82c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = -14.8305862683341019 * cart[4 * ncart + i];
        spherical[3 * nspherical + i] += 14.8305862683341019 * cart[11 * ncart + i];
        spherical[3 * nspherical + i] += 26.6950552830013805 * cart[22 * ncart + i];
        spherical[3 * nspherical + i] += -2.9661172536668201 * cart[37 * ncart + i];
        spherical[3 * nspherical + i] += 59.3223450733364075 * cart[13 * ncart + i];
        spherical[3 * nspherical + i] += -118.6446901466728150 * cart[24 * ncart + i];
        spherical[3 * nspherical + i] += 11.8644690146672804 * cart[39 * ncart + i];
    }
    // R_82s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = 1.6453058226360229 * cart[ncart + i];
        spherical[4 * nspherical + i] += 1.6453058226360229 * cart[6 * ncart + i];
        spherical[4 * nspherical + i] += -1.6453058226360229 * cart[15 * ncart + i];
        spherical[4 * nspherical + i] += -1.6453058226360229 * cart[28 * ncart + i];
        spherical[4 * nspherical + i] += -39.4873397432645490 * cart[8 * ncart + i];
        spherical[4 * nspherical + i] += 39.4873397432645490 * cart[30 * ncart + i];
        spherical[4 * nspherical + i] += 65.8122329054409221 * cart[19 * ncart + i];
        spherical[4 * nspherical + i] += -65.8122329054409221 * cart[32 * ncart + i];
    }

    // R_83c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i] = 9.5583630757311155 * cart[4 * ncart + i];
        spherical[5 * nspherical + i] += 15.9306051262185271 * cart[11 * ncart + i];
        spherical[5 * nspherical + i] += 3.1861210252437053 * cart[22 * ncart + i];
        spherical[5 * nspherical + i] += -3.1861210252437053 * cart[37 * ncart + i];
        spherical[5 * nspherical + i] += -63.7224205048741084 * cart[13 * ncart + i];
        spherical[5 * nspherical + i] += -42.4816136699160722 * cart[24 * ncart + i];
        spherical[5 * nspherical + i] += 21.2408068349580361 * cart[39 * ncart + i];
        spherical[5 * nspherical + i] += 50.9779364038992853 * cart[26 * ncart + i];
        spherical[5 * nspherical + i] += -16.9926454679664296 * cart[41 * ncart + i];
    }
    // R_83s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i] = -0.7843687748756958 * cart[ncart + i];
        spherical[6 * nspherical + i] += -2.3531063246270874 * cart[6 * ncart + i];
        spherical[6 * nspherical + i] += -2.3531063246270874 * cart[15 * ncart + i];
        spherical[6 * nspherical + i] += -0.7843687748756958 * cart[28 * ncart + i];
        spherical[6 * nspherical + i] += 23.5310632462708753 * cart[8 * ncart + i];
        spherical[6 * nspherical + i] += 47.0621264925417506 * cart[17 * ncart + i];
        spherical[6 * nspherical + i] += 23.5310632462708753 * cart[30 * ncart + i];
        spherical[6 * nspherical + i] += -62.7495019900556628 * cart[19 * ncart + i];
        spherical[6 * nspherical + i] += -62.7495019900556628 * cart[32 * ncart + i];
        spherical[6 * nspherical + i] += 25.0998007960222651 * cart[34 * ncart + i];
    }

    // R_84c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i] = -3.2812500000000000 * cart[4 * ncart + i];
        spherical[7 * nspherical + i] += -9.8437500000000000 * cart[11 * ncart + i];
        spherical[7 * nspherical + i] += -9.8437500000000000 * cart[22 * ncart + i];
        spherical[7 * nspherical + i] += -3.2812500000000000 * cart[37 * ncart + i];
        spherical[7 * nspherical + i] += 26.2500000000000000 * cart[13 * ncart + i];
        spherical[7 * nspherical + i] += 52.5000000000000000 * cart[24 * ncart + i];
        spherical[7 * nspherical + i] += 26.2500000000000000 * cart[39 * ncart + i];
        spherical[7 * nspherical + i] += -31.5000000000000000 * cart[26 * ncart + i];
        spherical[7 * nspherical + i] += -31.5000000000000000 * cart[41 * ncart + i];
        spherical[7 * nspherical + i] += 6.0000000000000000 * cart[43 * ncart + i];
    }
    // R_84s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i] = 0.2734375000000000 * cart[i];
        spherical[8 * nspherical + i] += 1.0937500000000000 * cart[3 * ncart + i];
        spherical[8 * nspherical + i] += 1.6406250000000000 * cart[10 * ncart + i];
        spherical[8 * nspherical + i] += 1.0937500000000000 * cart[21 * ncart + i];
        spherical[8 * nspherical + i] += 0.2734375000000000 * cart[36 * ncart + i];
        spherical[8 * nspherical + i] += -8.7500000000000000 * cart[5 * ncart + i];
        spherical[8 * nspherical + i] += -26.2500000000000000 * cart[12 * ncart + i];
        spherical[8 * nspherical + i] += -26.2500000000000000 * cart[23 * ncart + i];
        spherical[8 * nspherical + i] += -8.7500000000000000 * cart[38 * ncart + i];
        spherical[8 * nspherical + i] += 26.2500000000000000 * cart[14 * ncart + i];
        spherical[8 * nspherical + i] += 52.5000000000000000 * cart[25 * ncart + i];
        spherical[8 * nspherical + i] += 26.2500000000000000 * cart[40 * ncart + i];
        spherical[8 * nspherical + i] += -14.0000000000000000 * cart[27 * ncart + i];
        spherical[8 * nspherical + i] += -14.0000000000000000 * cart[42 * ncart + i];
        spherical[8 * nspherical + i] += cart[44 * ncart + i];
    }

    // R_85c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[9 * nspherical + i] = -3.2812500000000000 * cart[2 * ncart + i];
        spherical[9 * nspherical + i] += -9.8437500000000000 * cart[7 * ncart + i];
        spherical[9 * nspherical + i] += -9.8437500000000000 * cart[16 * ncart + i];
        spherical[9 * nspherical + i] += -3.2812500000000000 * cart[29 * ncart + i];
        spherical[9 * nspherical + i] += 26.2500000000000000 * cart[9 * ncart + i];
        spherical[9 * nspherical + i] += 52.5000000000000000 * cart[18 * ncart + i];
        spherical[9 * nspherical + i] += 26.2500000000000000 * cart[31 * ncart + i];
        spherical[9 * nspherical + i] += -31.5000000000000000 * cart[20 * ncart + i];
        spherical[9 * nspherical + i] += -31.5000000000000000 * cart[33 * ncart + i];
        spherical[9 * nspherical + i] += 6.0000000000000000 * cart[35 * ncart + i];
    }
    // R_85s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[10 * nspherical + i] = -0.3921843874378479 * cart[i];
        spherical[10 * nspherical + i] += -0.7843687748756958 * cart[3 * ncart + i];
        spherical[10 * nspherical + i] += 0.7843687748756958 * cart[21 * ncart + i];
        spherical[10 * nspherical + i] += 0.3921843874378479 * cart[36 * ncart + i];
        spherical[10 * nspherical + i] += 11.7655316231354377 * cart[5 * ncart + i];
        spherical[10 * nspherical + i] += 11.7655316231354377 * cart[12 * ncart + i];
        spherical[10 * nspherical + i] += -11.7655316231354377 * cart[23 * ncart + i];
        spherical[10 * nspherical + i] += -11.7655316231354377 * cart[38 * ncart + i];
        spherical[10 * nspherical + i] += -31.3747509950278314 * cart[14 * ncart + i];
        spherical[10 * nspherical + i] += 31.3747509950278314 * cart[40 * ncart + i];
        spherical[10 * nspherical + i] += 12.5499003980111326 * cart[27 * ncart + i];
        spherical[10 * nspherical + i] += -12.5499003980111326 * cart[42 * ncart + i];
    }

    // R_86c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[11 * nspherical + i] = 3.1861210252437053 * cart[2 * ncart + i];
        spherical[11 * nspherical + i] += -3.1861210252437053 * cart[7 * ncart + i];
        spherical[11 * nspherical + i] += -15.9306051262185271 * cart[16 * ncart + i];
        spherical[11 * nspherical + i] += -9.5583630757311155 * cart[29 * ncart + i];
        spherical[11 * nspherical + i] += -21.2408068349580361 * cart[9 * ncart + i];
        spherical[11 * nspherical + i] += 42.4816136699160722 * cart[18 * ncart + i];
        spherical[11 * nspherical + i] += 63.7224205048741084 * cart[31 * ncart + i];
        spherical[11 * nspherical + i] += 16.9926454679664296 * cart[20 * ncart + i];
        spherical[11 * nspherical + i] += -50.9779364038992853 * cart[33 * ncart + i];
    }
    // R_86s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[12 * nspherical + i] = 0.4113264556590057 * cart[i];
        spherical[12 * nspherical + i] += -1.6453058226360229 * cart[3 * ncart + i];
        spherical[12 * nspherical + i] += -4.1132645565900576 * cart[10 * ncart + i];
        spherical[12 * nspherical + i] += -1.6453058226360229 * cart[21 * ncart + i];
        spherical[12 * nspherical + i] += 0.4113264556590057 * cart[36 * ncart + i];
        spherical[12 * nspherical + i] += -9.8718349358161372 * cart[5 * ncart + i];
        spherical[12 * nspherical + i] += 49.3591746790806880 * cart[12 * ncart + i];
        spherical[12 * nspherical + i] += 49.3591746790806880 * cart[23 * ncart + i];
        spherical[12 * nspherical + i] += -9.8718349358161372 * cart[38 * ncart + i];
        spherical[12 * nspherical + i] += 16.4530582263602305 * cart[14 * ncart + i];
        spherical[12 * nspherical + i] += -98.7183493581613760 * cart[25 * ncart + i];
        spherical[12 * nspherical + i] += 16.4530582263602305 * cart[40 * ncart + i];
    }

    // R_87c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[13 * nspherical + i] = -2.9661172536668201 * cart[2 * ncart + i];
        spherical[13 * nspherical + i] += 26.6950552830013805 * cart[7 * ncart + i];
        spherical[13 * nspherical + i] += 14.8305862683341019 * cart[16 * ncart + i];
        spherical[13 * nspherical + i] += -14.8305862683341019 * cart[29 * ncart + i];
        spherical[13 * nspherical + i] += 11.8644690146672804 * cart[9 * ncart + i];
        spherical[13 * nspherical + i] += -118.6446901466728150 * cart[18 * ncart + i];
        spherical[13 * nspherical + i] += 59.3223450733364075 * cart[31 * ncart + i];
    }
    // R_87s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[14 * nspherical + i] = -0.4576818286211503 * cart[i];
        spherical[14 * nspherical + i] += 6.4075456006961042 * cart[3 * ncart + i];
        spherical[14 * nspherical + i] += -6.4075456006961042 * cart[21 * ncart + i];
        spherical[14 * nspherical + i] += 0.4576818286211503 * cart[36 * ncart + i];
        spherical[14 * nspherical + i] += 6.4075456006961042 * cart[5 * ncart + i];
        spherical[14 * nspherical + i] += -96.1131840104415573 * cart[12 * ncart + i];
        spherical[14 * nspherical + i] += 96.1131840104415573 * cart[23 * ncart + i];
        spherical[14 * nspherical + i] += -6.4075456006961042 * cart[38 * ncart + i];
    }

    // R_88c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[15 * nspherical + i] = 2.5068266169601756 * cart[2 * ncart + i];
        spherical[15 * nspherical + i] += -52.6433589561636950 * cart[7 * ncart + i];
        spherical[15 * nspherical + i] += 87.7389315936061536 * cart[16 * ncart + i];
        spherical[15 * nspherical + i] += -17.5477863187212293 * cart[29 * ncart + i];
    }
    // R_88s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[16 * nspherical + i] = 0.6267066542400439 * cart[i];
        spherical[16 * nspherical + i] += -17.5477863187212293 * cart[3 * ncart + i];
        spherical[16 * nspherical + i] += 43.8694657968030768 * cart[10 * ncart + i];
        spherical[16 * nspherical + i] += -17.5477863187212293 * cart[21 * ncart + i];
        spherical[16 * nspherical + i] += 0.6267066542400439 * cart[36 * ncart + i];
    }
}
void gg_cca_cart_to_spherical_sum_L8(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart,
                                     const unsigned long ncart, double* PRAGMA_RESTRICT output,
                                     const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_80 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 5.0136532339203512 * cart[ncart + i];
        tmp += -35.0955726374424586 * cart[6 * ncart + i];
        tmp += 35.0955726374424586 * cart[15 * ncart + i];
        tmp += -5.0136532339203512 * cart[28 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_81c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 17.5477863187212293 * cart[4 * ncart + i];
        tmp += -87.7389315936061536 * cart[11 * ncart + i];
        tmp += 52.6433589561636950 * cart[22 * ncart + i];
        tmp += -2.5068266169601756 * cart[37 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_81s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.7460909717269018 * cart[ncart + i];
        tmp += 6.4075456006961042 * cart[6 * ncart + i];
        tmp += 6.4075456006961042 * cart[15 * ncart + i];
        tmp += -2.7460909717269018 * cart[28 * ncart + i];
        tmp += 38.4452736041766272 * cart[8 * ncart + i];
        tmp += -128.1509120139220954 * cart[17 * ncart + i];
        tmp += 38.4452736041766272 * cart[30 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_82c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -14.8305862683341019 * cart[4 * ncart + i];
        tmp += 14.8305862683341019 * cart[11 * ncart + i];
        tmp += 26.6950552830013805 * cart[22 * ncart + i];
        tmp += -2.9661172536668201 * cart[37 * ncart + i];
        tmp += 59.3223450733364075 * cart[13 * ncart + i];
        tmp += -118.6446901466728150 * cart[24 * ncart + i];
        tmp += 11.8644690146672804 * cart[39 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_82s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.6453058226360229 * cart[ncart + i];
        tmp += 1.6453058226360229 * cart[6 * ncart + i];
        tmp += -1.6453058226360229 * cart[15 * ncart + i];
        tmp += -1.6453058226360229 * cart[28 * ncart + i];
        tmp += -39.4873397432645490 * cart[8 * ncart + i];
        tmp += 39.4873397432645490 * cart[30 * ncart + i];
        tmp += 65.8122329054409221 * cart[19 * ncart + i];
        tmp += -65.8122329054409221 * cart[32 * ncart + i];
        output[i] += tmp * vector[4];
    }

    // R_83c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 9.5583630757311155 * cart[4 * ncart + i];
        tmp += 15.9306051262185271 * cart[11 * ncart + i];
        tmp += 3.1861210252437053 * cart[22 * ncart + i];
        tmp += -3.1861210252437053 * cart[37 * ncart + i];
        tmp += -63.7224205048741084 * cart[13 * ncart + i];
        tmp += -42.4816136699160722 * cart[24 * ncart + i];
        tmp += 21.2408068349580361 * cart[39 * ncart + i];
        tmp += 50.9779364038992853 * cart[26 * ncart + i];
        tmp += -16.9926454679664296 * cart[41 * ncart + i];
        output[i] += tmp * vector[5];
    }
    // R_83s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.7843687748756958 * cart[ncart + i];
        tmp += -2.3531063246270874 * cart[6 * ncart + i];
        tmp += -2.3531063246270874 * cart[15 * ncart + i];
        tmp += -0.7843687748756958 * cart[28 * ncart + i];
        tmp += 23.5310632462708753 * cart[8 * ncart + i];
        tmp += 47.0621264925417506 * cart[17 * ncart + i];
        tmp += 23.5310632462708753 * cart[30 * ncart + i];
        tmp += -62.7495019900556628 * cart[19 * ncart + i];
        tmp += -62.7495019900556628 * cart[32 * ncart + i];
        tmp += 25.0998007960222651 * cart[34 * ncart + i];
        output[i] += tmp * vector[6];
    }

    // R_84c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -3.2812500000000000 * cart[4 * ncart + i];
        tmp += -9.8437500000000000 * cart[11 * ncart + i];
        tmp += -9.8437500000000000 * cart[22 * ncart + i];
        tmp += -3.2812500000000000 * cart[37 * ncart + i];
        tmp += 26.2500000000000000 * cart[13 * ncart + i];
        tmp += 52.5000000000000000 * cart[24 * ncart + i];
        tmp += 26.2500000000000000 * cart[39 * ncart + i];
        tmp += -31.5000000000000000 * cart[26 * ncart + i];
        tmp += -31.5000000000000000 * cart[41 * ncart + i];
        tmp += 6.0000000000000000 * cart[43 * ncart + i];
        output[i] += tmp * vector[7];
    }
    // R_84s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.2734375000000000 * cart[i];
        tmp += 1.0937500000000000 * cart[3 * ncart + i];
        tmp += 1.6406250000000000 * cart[10 * ncart + i];
        tmp += 1.0937500000000000 * cart[21 * ncart + i];
        tmp += 0.2734375000000000 * cart[36 * ncart + i];
        tmp += -8.7500000000000000 * cart[5 * ncart + i];
        tmp += -26.2500000000000000 * cart[12 * ncart + i];
        tmp += -26.2500000000000000 * cart[23 * ncart + i];
        tmp += -8.7500000000000000 * cart[38 * ncart + i];
        tmp += 26.2500000000000000 * cart[14 * ncart + i];
        tmp += 52.5000000000000000 * cart[25 * ncart + i];
        tmp += 26.2500000000000000 * cart[40 * ncart + i];
        tmp += -14.0000000000000000 * cart[27 * ncart + i];
        tmp += -14.0000000000000000 * cart[42 * ncart + i];
        tmp += cart[44 * ncart + i];
        output[i] += tmp * vector[8];
    }

    // R_85c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -3.2812500000000000 * cart[2 * ncart + i];
        tmp += -9.8437500000000000 * cart[7 * ncart + i];
        tmp += -9.8437500000000000 * cart[16 * ncart + i];
        tmp += -3.2812500000000000 * cart[29 * ncart + i];
        tmp += 26.2500000000000000 * cart[9 * ncart + i];
        tmp += 52.5000000000000000 * cart[18 * ncart + i];
        tmp += 26.2500000000000000 * cart[31 * ncart + i];
        tmp += -31.5000000000000000 * cart[20 * ncart + i];
        tmp += -31.5000000000000000 * cart[33 * ncart + i];
        tmp += 6.0000000000000000 * cart[35 * ncart + i];
        output[i] += tmp * vector[9];
    }
    // R_85s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.3921843874378479 * cart[i];
        tmp += -0.7843687748756958 * cart[3 * ncart + i];
        tmp += 0.7843687748756958 * cart[21 * ncart + i];
        tmp += 0.3921843874378479 * cart[36 * ncart + i];
        tmp += 11.7655316231354377 * cart[5 * ncart + i];
        tmp += 11.7655316231354377 * cart[12 * ncart + i];
        tmp += -11.7655316231354377 * cart[23 * ncart + i];
        tmp += -11.7655316231354377 * cart[38 * ncart + i];
        tmp += -31.3747509950278314 * cart[14 * ncart + i];
        tmp += 31.3747509950278314 * cart[40 * ncart + i];
        tmp += 12.5499003980111326 * cart[27 * ncart + i];
        tmp += -12.5499003980111326 * cart[42 * ncart + i];
        output[i] += tmp * vector[10];
    }

    // R_86c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 3.1861210252437053 * cart[2 * ncart + i];
        tmp += -3.1861210252437053 * cart[7 * ncart + i];
        tmp += -15.9306051262185271 * cart[16 * ncart + i];
        tmp += -9.5583630757311155 * cart[29 * ncart + i];
        tmp += -21.2408068349580361 * cart[9 * ncart + i];
        tmp += 42.4816136699160722 * cart[18 * ncart + i];
        tmp += 63.7224205048741084 * cart[31 * ncart + i];
        tmp += 16.9926454679664296 * cart[20 * ncart + i];
        tmp += -50.9779364038992853 * cart[33 * ncart + i];
        output[i] += tmp * vector[11];
    }
    // R_86s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.4113264556590057 * cart[i];
        tmp += -1.6453058226360229 * cart[3 * ncart + i];
        tmp += -4.1132645565900576 * cart[10 * ncart + i];
        tmp += -1.6453058226360229 * cart[21 * ncart + i];
        tmp += 0.4113264556590057 * cart[36 * ncart + i];
        tmp += -9.8718349358161372 * cart[5 * ncart + i];
        tmp += 49.3591746790806880 * cart[12 * ncart + i];
        tmp += 49.3591746790806880 * cart[23 * ncart + i];
        tmp += -9.8718349358161372 * cart[38 * ncart + i];
        tmp += 16.4530582263602305 * cart[14 * ncart + i];
        tmp += -98.7183493581613760 * cart[25 * ncart + i];
        tmp += 16.4530582263602305 * cart[40 * ncart + i];
        output[i] += tmp * vector[12];
    }

    // R_87c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.9661172536668201 * cart[2 * ncart + i];
        tmp += 26.6950552830013805 * cart[7 * ncart + i];
        tmp += 14.8305862683341019 * cart[16 * ncart + i];
        tmp += -14.8305862683341019 * cart[29 * ncart + i];
        tmp += 11.8644690146672804 * cart[9 * ncart + i];
        tmp += -118.6446901466728150 * cart[18 * ncart + i];
        tmp += 59.3223450733364075 * cart[31 * ncart + i];
        output[i] += tmp * vector[13];
    }
    // R_87s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.4576818286211503 * cart[i];
        tmp += 6.4075456006961042 * cart[3 * ncart + i];
        tmp += -6.4075456006961042 * cart[21 * ncart + i];
        tmp += 0.4576818286211503 * cart[36 * ncart + i];
        tmp += 6.4075456006961042 * cart[5 * ncart + i];
        tmp += -96.1131840104415573 * cart[12 * ncart + i];
        tmp += 96.1131840104415573 * cart[23 * ncart + i];
        tmp += -6.4075456006961042 * cart[38 * ncart + i];
        output[i] += tmp * vector[14];
    }

    // R_88c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.5068266169601756 * cart[2 * ncart + i];
        tmp += -52.6433589561636950 * cart[7 * ncart + i];
        tmp += 87.7389315936061536 * cart[16 * ncart + i];
        tmp += -17.5477863187212293 * cart[29 * ncart + i];
        output[i] += tmp * vector[15];
    }
    // R_88s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.6267066542400439 * cart[i];
        tmp += -17.5477863187212293 * cart[3 * ncart + i];
        tmp += 43.8694657968030768 * cart[10 * ncart + i];
        tmp += -17.5477863187212293 * cart[21 * ncart + i];
        tmp += 0.6267066542400439 * cart[36 * ncart + i];
        output[i] += tmp * vector[16];
    }
}
void gg_gaussian_cart_to_spherical_L0(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                      const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                      const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_00 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = cart[i];
    }
}
void gg_gaussian_cart_to_spherical_sum_L0(const unsigned long size, const double* vector,
                                          const double* PRAGMA_RESTRICT cart, const unsigned long ncart,
                                          double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_00 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = cart[i];
        output[i] += tmp * vector[0];
    }
}
void gg_gaussian_cart_to_spherical_L1(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                      const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                      const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_10 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = cart[2 * ncart + i];
    }

    // R_11c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = cart[i];
    }
    // R_11s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = cart[ncart + i];
    }
}
void gg_gaussian_cart_to_spherical_sum_L1(const unsigned long size, const double* vector,
                                          const double* PRAGMA_RESTRICT cart, const unsigned long ncart,
                                          double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_10 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = cart[2 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_11c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = cart[i];
        output[i] += tmp * vector[1];
    }
    // R_11s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = cart[ncart + i];
        output[i] += tmp * vector[2];
    }
}
void gg_gaussian_cart_to_spherical_L2(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                      const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                      const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_20 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = -0.5000000000000000 * cart[i];
        spherical[i] += -0.5000000000000000 * cart[3 * ncart + i];
        spherical[i] += cart[5 * ncart + i];
    }

    // R_21c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = 1.7320508075688772 * cart[2 * ncart + i];
    }
    // R_21s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = 1.7320508075688772 * cart[4 * ncart + i];
    }

    // R_22c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = 0.8660254037844386 * cart[i];
        spherical[3 * nspherical + i] += -0.8660254037844386 * cart[3 * ncart + i];
    }
    // R_22s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = 1.7320508075688772 * cart[ncart + i];
    }
}
void gg_gaussian_cart_to_spherical_sum_L2(const unsigned long size, const double* vector,
                                          const double* PRAGMA_RESTRICT cart, const unsigned long ncart,
                                          double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_20 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.5000000000000000 * cart[i];
        tmp += -0.5000000000000000 * cart[3 * ncart + i];
        tmp += cart[5 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_21c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.7320508075688772 * cart[2 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_21s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.7320508075688772 * cart[4 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_22c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.8660254037844386 * cart[i];
        tmp += -0.8660254037844386 * cart[3 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_22s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.7320508075688772 * cart[ncart + i];
        output[i] += tmp * vector[4];
    }
}
void gg_gaussian_cart_to_spherical_L3(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                      const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                      const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_30 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = -1.5000000000000000 * cart[2 * ncart + i];
        spherical[i] += -1.5000000000000000 * cart[7 * ncart + i];
        spherical[i] += cart[9 * ncart + i];
    }

    // R_31c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = -0.6123724356957945 * cart[i];
        spherical[nspherical + i] += -0.6123724356957945 * cart[3 * ncart + i];
        spherical[nspherical + i] += 2.4494897427831779 * cart[5 * ncart + i];
    }
    // R_31s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = -0.6123724356957945 * cart[ncart + i];
        spherical[2 * nspherical + i] += -0.6123724356957945 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] += 2.4494897427831779 * cart[8 * ncart + i];
    }

    // R_32c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = 1.9364916731037085 * cart[2 * ncart + i];
        spherical[3 * nspherical + i] += -1.9364916731037085 * cart[7 * ncart + i];
    }
    // R_32s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = 3.8729833462074170 * cart[4 * ncart + i];
    }

    // R_33c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i] = 0.7905694150420949 * cart[i];
        spherical[5 * nspherical + i] += -2.3717082451262845 * cart[3 * ncart + i];
    }
    // R_33s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i] = 2.3717082451262845 * cart[ncart + i];
        spherical[6 * nspherical + i] += -0.7905694150420949 * cart[6 * ncart + i];
    }
}
void gg_gaussian_cart_to_spherical_sum_L3(const unsigned long size, const double* vector,
                                          const double* PRAGMA_RESTRICT cart, const unsigned long ncart,
                                          double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_30 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -1.5000000000000000 * cart[2 * ncart + i];
        tmp += -1.5000000000000000 * cart[7 * ncart + i];
        tmp += cart[9 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_31c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.6123724356957945 * cart[i];
        tmp += -0.6123724356957945 * cart[3 * ncart + i];
        tmp += 2.4494897427831779 * cart[5 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_31s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.6123724356957945 * cart[ncart + i];
        tmp += -0.6123724356957945 * cart[6 * ncart + i];
        tmp += 2.4494897427831779 * cart[8 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_32c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.9364916731037085 * cart[2 * ncart + i];
        tmp += -1.9364916731037085 * cart[7 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_32s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 3.8729833462074170 * cart[4 * ncart + i];
        output[i] += tmp * vector[4];
    }

    // R_33c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.7905694150420949 * cart[i];
        tmp += -2.3717082451262845 * cart[3 * ncart + i];
        output[i] += tmp * vector[5];
    }
    // R_33s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.3717082451262845 * cart[ncart + i];
        tmp += -0.7905694150420949 * cart[6 * ncart + i];
        output[i] += tmp * vector[6];
    }
}
void gg_gaussian_cart_to_spherical_L4(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                      const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                      const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_40 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = 0.3750000000000000 * cart[i];
        spherical[i] += 0.7500000000000000 * cart[3 * ncart + i];
        spherical[i] += 0.3750000000000000 * cart[10 * ncart + i];
        spherical[i] += -3.0000000000000000 * cart[5 * ncart + i];
        spherical[i] += -3.0000000000000000 * cart[12 * ncart + i];
        spherical[i] += cart[14 * ncart + i];
    }

    // R_41c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = -2.3717082451262845 * cart[2 * ncart + i];
        spherical[nspherical + i] += -2.3717082451262845 * cart[7 * ncart + i];
        spherical[nspherical + i] += 3.1622776601683795 * cart[9 * ncart + i];
    }
    // R_41s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = -2.3717082451262845 * cart[4 * ncart + i];
        spherical[2 * nspherical + i] += -2.3717082451262845 * cart[11 * ncart + i];
        spherical[2 * nspherical + i] += 3.1622776601683795 * cart[13 * ncart + i];
    }

    // R_42c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = -0.5590169943749475 * cart[i];
        spherical[3 * nspherical + i] += 0.5590169943749475 * cart[10 * ncart + i];
        spherical[3 * nspherical + i] += 3.3541019662496847 * cart[5 * ncart + i];
        spherical[3 * nspherical + i] += -3.3541019662496847 * cart[12 * ncart + i];
    }
    // R_42s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = -1.1180339887498949 * cart[ncart + i];
        spherical[4 * nspherical + i] += -1.1180339887498949 * cart[6 * ncart + i];
        spherical[4 * nspherical + i] += 6.7082039324993694 * cart[8 * ncart + i];
    }

    // R_43c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i] = 2.0916500663351889 * cart[2 * ncart + i];
        spherical[5 * nspherical + i] += -6.2749501990055663 * cart[7 * ncart + i];
    }
    // R_43s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i] = 6.2749501990055663 * cart[4 * ncart + i];
        spherical[6 * nspherical + i] += -2.0916500663351889 * cart[11 * ncart + i];
    }

    // R_44c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i] = 0.7395099728874520 * cart[i];
        spherical[7 * nspherical + i] += -4.4370598373247123 * cart[3 * ncart + i];
        spherical[7 * nspherical + i] += 0.7395099728874520 * cart[10 * ncart + i];
    }
    // R_44s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i] = 2.9580398915498081 * cart[ncart + i];
        spherical[8 * nspherical + i] += -2.9580398915498081 * cart[6 * ncart + i];
    }
}
void gg_gaussian_cart_to_spherical_sum_L4(const unsigned long size, const double* vector,
                                          const double* PRAGMA_RESTRICT cart, const unsigned long ncart,
                                          double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_40 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.3750000000000000 * cart[i];
        tmp += 0.7500000000000000 * cart[3 * ncart + i];
        tmp += 0.3750000000000000 * cart[10 * ncart + i];
        tmp += -3.0000000000000000 * cart[5 * ncart + i];
        tmp += -3.0000000000000000 * cart[12 * ncart + i];
        tmp += cart[14 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_41c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.3717082451262845 * cart[2 * ncart + i];
        tmp += -2.3717082451262845 * cart[7 * ncart + i];
        tmp += 3.1622776601683795 * cart[9 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_41s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.3717082451262845 * cart[4 * ncart + i];
        tmp += -2.3717082451262845 * cart[11 * ncart + i];
        tmp += 3.1622776601683795 * cart[13 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_42c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.5590169943749475 * cart[i];
        tmp += 0.5590169943749475 * cart[10 * ncart + i];
        tmp += 3.3541019662496847 * cart[5 * ncart + i];
        tmp += -3.3541019662496847 * cart[12 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_42s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -1.1180339887498949 * cart[ncart + i];
        tmp += -1.1180339887498949 * cart[6 * ncart + i];
        tmp += 6.7082039324993694 * cart[8 * ncart + i];
        output[i] += tmp * vector[4];
    }

    // R_43c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.0916500663351889 * cart[2 * ncart + i];
        tmp += -6.2749501990055663 * cart[7 * ncart + i];
        output[i] += tmp * vector[5];
    }
    // R_43s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 6.2749501990055663 * cart[4 * ncart + i];
        tmp += -2.0916500663351889 * cart[11 * ncart + i];
        output[i] += tmp * vector[6];
    }

    // R_44c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.7395099728874520 * cart[i];
        tmp += -4.4370598373247123 * cart[3 * ncart + i];
        tmp += 0.7395099728874520 * cart[10 * ncart + i];
        output[i] += tmp * vector[7];
    }
    // R_44s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.9580398915498081 * cart[ncart + i];
        tmp += -2.9580398915498081 * cart[6 * ncart + i];
        output[i] += tmp * vector[8];
    }
}
void gg_gaussian_cart_to_spherical_L5(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                      const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                      const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_50 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = 1.8750000000000000 * cart[2 * ncart + i];
        spherical[i] += 3.7500000000000000 * cart[7 * ncart + i];
        spherical[i] += 1.8750000000000000 * cart[16 * ncart + i];
        spherical[i] += -5.0000000000000000 * cart[9 * ncart + i];
        spherical[i] += -5.0000000000000000 * cart[18 * ncart + i];
        spherical[i] += cart[20 * ncart + i];
    }

    // R_51c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = 0.4841229182759271 * cart[i];
        spherical[nspherical + i] += 0.9682458365518543 * cart[3 * ncart + i];
        spherical[nspherical + i] += 0.4841229182759271 * cart[10 * ncart + i];
        spherical[nspherical + i] += -5.8094750193111251 * cart[5 * ncart + i];
        spherical[nspherical + i] += -5.8094750193111251 * cart[12 * ncart + i];
        spherical[nspherical + i] += 3.8729833462074170 * cart[14 * ncart + i];
    }
    // R_51s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = 0.4841229182759271 * cart[ncart + i];
        spherical[2 * nspherical + i] += 0.9682458365518543 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] += 0.4841229182759271 * cart[15 * ncart + i];
        spherical[2 * nspherical + i] += -5.8094750193111251 * cart[8 * ncart + i];
        spherical[2 * nspherical + i] += -5.8094750193111251 * cart[17 * ncart + i];
        spherical[2 * nspherical + i] += 3.8729833462074170 * cart[19 * ncart + i];
    }

    // R_52c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = -2.5617376914898995 * cart[2 * ncart + i];
        spherical[3 * nspherical + i] += 2.5617376914898995 * cart[16 * ncart + i];
        spherical[3 * nspherical + i] += 5.1234753829797990 * cart[9 * ncart + i];
        spherical[3 * nspherical + i] += -5.1234753829797990 * cart[18 * ncart + i];
    }
    // R_52s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = -5.1234753829797990 * cart[4 * ncart + i];
        spherical[4 * nspherical + i] += -5.1234753829797990 * cart[11 * ncart + i];
        spherical[4 * nspherical + i] += 10.2469507659595980 * cart[13 * ncart + i];
    }

    // R_53c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i] = -0.5229125165837972 * cart[i];
        spherical[5 * nspherical + i] += 1.0458250331675945 * cart[3 * ncart + i];
        spherical[5 * nspherical + i] += 1.5687375497513916 * cart[10 * ncart + i];
        spherical[5 * nspherical + i] += 4.1833001326703778 * cart[5 * ncart + i];
        spherical[5 * nspherical + i] += -12.5499003980111326 * cart[12 * ncart + i];
    }
    // R_53s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i] = -1.5687375497513916 * cart[ncart + i];
        spherical[6 * nspherical + i] += -1.0458250331675945 * cart[6 * ncart + i];
        spherical[6 * nspherical + i] += 0.5229125165837972 * cart[15 * ncart + i];
        spherical[6 * nspherical + i] += 12.5499003980111326 * cart[8 * ncart + i];
        spherical[6 * nspherical + i] += -4.1833001326703778 * cart[17 * ncart + i];
    }

    // R_54c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i] = 2.2185299186623562 * cart[2 * ncart + i];
        spherical[7 * nspherical + i] += -13.3111795119741370 * cart[7 * ncart + i];
        spherical[7 * nspherical + i] += 2.2185299186623562 * cart[16 * ncart + i];
    }
    // R_54s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i] = 8.8741196746494246 * cart[4 * ncart + i];
        spherical[8 * nspherical + i] += -8.8741196746494246 * cart[11 * ncart + i];
    }

    // R_55c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[9 * nspherical + i] = 0.7015607600201140 * cart[i];
        spherical[9 * nspherical + i] += -7.0156076002011405 * cart[3 * ncart + i];
        spherical[9 * nspherical + i] += 3.5078038001005702 * cart[10 * ncart + i];
    }
    // R_55s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[10 * nspherical + i] = 3.5078038001005702 * cart[ncart + i];
        spherical[10 * nspherical + i] += -7.0156076002011405 * cart[6 * ncart + i];
        spherical[10 * nspherical + i] += 0.7015607600201140 * cart[15 * ncart + i];
    }
}
void gg_gaussian_cart_to_spherical_sum_L5(const unsigned long size, const double* vector,
                                          const double* PRAGMA_RESTRICT cart, const unsigned long ncart,
                                          double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_50 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.8750000000000000 * cart[2 * ncart + i];
        tmp += 3.7500000000000000 * cart[7 * ncart + i];
        tmp += 1.8750000000000000 * cart[16 * ncart + i];
        tmp += -5.0000000000000000 * cart[9 * ncart + i];
        tmp += -5.0000000000000000 * cart[18 * ncart + i];
        tmp += cart[20 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_51c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.4841229182759271 * cart[i];
        tmp += 0.9682458365518543 * cart[3 * ncart + i];
        tmp += 0.4841229182759271 * cart[10 * ncart + i];
        tmp += -5.8094750193111251 * cart[5 * ncart + i];
        tmp += -5.8094750193111251 * cart[12 * ncart + i];
        tmp += 3.8729833462074170 * cart[14 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_51s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.4841229182759271 * cart[ncart + i];
        tmp += 0.9682458365518543 * cart[6 * ncart + i];
        tmp += 0.4841229182759271 * cart[15 * ncart + i];
        tmp += -5.8094750193111251 * cart[8 * ncart + i];
        tmp += -5.8094750193111251 * cart[17 * ncart + i];
        tmp += 3.8729833462074170 * cart[19 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_52c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.5617376914898995 * cart[2 * ncart + i];
        tmp += 2.5617376914898995 * cart[16 * ncart + i];
        tmp += 5.1234753829797990 * cart[9 * ncart + i];
        tmp += -5.1234753829797990 * cart[18 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_52s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -5.1234753829797990 * cart[4 * ncart + i];
        tmp += -5.1234753829797990 * cart[11 * ncart + i];
        tmp += 10.2469507659595980 * cart[13 * ncart + i];
        output[i] += tmp * vector[4];
    }

    // R_53c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.5229125165837972 * cart[i];
        tmp += 1.0458250331675945 * cart[3 * ncart + i];
        tmp += 1.5687375497513916 * cart[10 * ncart + i];
        tmp += 4.1833001326703778 * cart[5 * ncart + i];
        tmp += -12.5499003980111326 * cart[12 * ncart + i];
        output[i] += tmp * vector[5];
    }
    // R_53s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -1.5687375497513916 * cart[ncart + i];
        tmp += -1.0458250331675945 * cart[6 * ncart + i];
        tmp += 0.5229125165837972 * cart[15 * ncart + i];
        tmp += 12.5499003980111326 * cart[8 * ncart + i];
        tmp += -4.1833001326703778 * cart[17 * ncart + i];
        output[i] += tmp * vector[6];
    }

    // R_54c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.2185299186623562 * cart[2 * ncart + i];
        tmp += -13.3111795119741370 * cart[7 * ncart + i];
        tmp += 2.2185299186623562 * cart[16 * ncart + i];
        output[i] += tmp * vector[7];
    }
    // R_54s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 8.8741196746494246 * cart[4 * ncart + i];
        tmp += -8.8741196746494246 * cart[11 * ncart + i];
        output[i] += tmp * vector[8];
    }

    // R_55c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.7015607600201140 * cart[i];
        tmp += -7.0156076002011405 * cart[3 * ncart + i];
        tmp += 3.5078038001005702 * cart[10 * ncart + i];
        output[i] += tmp * vector[9];
    }
    // R_55s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 3.5078038001005702 * cart[ncart + i];
        tmp += -7.0156076002011405 * cart[6 * ncart + i];
        tmp += 0.7015607600201140 * cart[15 * ncart + i];
        output[i] += tmp * vector[10];
    }
}
void gg_gaussian_cart_to_spherical_L6(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                      const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                      const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_60 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = -0.3125000000000000 * cart[i];
        spherical[i] += -0.9375000000000000 * cart[3 * ncart + i];
        spherical[i] += -0.9375000000000000 * cart[10 * ncart + i];
        spherical[i] += -0.3125000000000000 * cart[21 * ncart + i];
        spherical[i] += 5.6250000000000000 * cart[5 * ncart + i];
        spherical[i] += 11.2500000000000000 * cart[12 * ncart + i];
        spherical[i] += 5.6250000000000000 * cart[23 * ncart + i];
        spherical[i] += -7.5000000000000000 * cart[14 * ncart + i];
        spherical[i] += -7.5000000000000000 * cart[25 * ncart + i];
        spherical[i] += cart[27 * ncart + i];
    }

    // R_61c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = 2.8641098093473998 * cart[2 * ncart + i];
        spherical[nspherical + i] += 5.7282196186947996 * cart[7 * ncart + i];
        spherical[nspherical + i] += 2.8641098093473998 * cart[16 * ncart + i];
        spherical[nspherical + i] += -11.4564392373895991 * cart[9 * ncart + i];
        spherical[nspherical + i] += -11.4564392373895991 * cart[18 * ncart + i];
        spherical[nspherical + i] += 4.5825756949558398 * cart[20 * ncart + i];
    }
    // R_61s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = 2.8641098093473998 * cart[4 * ncart + i];
        spherical[2 * nspherical + i] += 5.7282196186947996 * cart[11 * ncart + i];
        spherical[2 * nspherical + i] += 2.8641098093473998 * cart[22 * ncart + i];
        spherical[2 * nspherical + i] += -11.4564392373895991 * cart[13 * ncart + i];
        spherical[2 * nspherical + i] += -11.4564392373895991 * cart[24 * ncart + i];
        spherical[2 * nspherical + i] += 4.5825756949558398 * cart[26 * ncart + i];
    }

    // R_62c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = 0.4528555233184199 * cart[i];
        spherical[3 * nspherical + i] += 0.4528555233184199 * cart[3 * ncart + i];
        spherical[3 * nspherical + i] += -0.4528555233184199 * cart[10 * ncart + i];
        spherical[3 * nspherical + i] += -0.4528555233184199 * cart[21 * ncart + i];
        spherical[3 * nspherical + i] += -7.2456883730947190 * cart[5 * ncart + i];
        spherical[3 * nspherical + i] += 7.2456883730947190 * cart[23 * ncart + i];
        spherical[3 * nspherical + i] += 7.2456883730947190 * cart[14 * ncart + i];
        spherical[3 * nspherical + i] += -7.2456883730947190 * cart[25 * ncart + i];
    }
    // R_62s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = 0.9057110466368399 * cart[ncart + i];
        spherical[4 * nspherical + i] += 1.8114220932736798 * cart[6 * ncart + i];
        spherical[4 * nspherical + i] += 0.9057110466368399 * cart[15 * ncart + i];
        spherical[4 * nspherical + i] += -14.4913767461894381 * cart[8 * ncart + i];
        spherical[4 * nspherical + i] += -14.4913767461894381 * cart[17 * ncart + i];
        spherical[4 * nspherical + i] += 14.4913767461894381 * cart[19 * ncart + i];
    }

    // R_63c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i] = -2.7171331399105196 * cart[2 * ncart + i];
        spherical[5 * nspherical + i] += 5.4342662798210393 * cart[7 * ncart + i];
        spherical[5 * nspherical + i] += 8.1513994197315593 * cart[16 * ncart + i];
        spherical[5 * nspherical + i] += 7.2456883730947190 * cart[9 * ncart + i];
        spherical[5 * nspherical + i] += -21.7370651192841571 * cart[18 * ncart + i];
    }
    // R_63s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i] = -8.1513994197315593 * cart[4 * ncart + i];
        spherical[6 * nspherical + i] += -5.4342662798210393 * cart[11 * ncart + i];
        spherical[6 * nspherical + i] += 2.7171331399105196 * cart[22 * ncart + i];
        spherical[6 * nspherical + i] += 21.7370651192841571 * cart[13 * ncart + i];
        spherical[6 * nspherical + i] += -7.2456883730947190 * cart[24 * ncart + i];
    }

    // R_64c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i] = -0.4960783708246108 * cart[i];
        spherical[7 * nspherical + i] += 2.4803918541230536 * cart[3 * ncart + i];
        spherical[7 * nspherical + i] += 2.4803918541230536 * cart[10 * ncart + i];
        spherical[7 * nspherical + i] += -0.4960783708246108 * cart[21 * ncart + i];
        spherical[7 * nspherical + i] += 4.9607837082461073 * cart[5 * ncart + i];
        spherical[7 * nspherical + i] += -29.7647022494766453 * cart[12 * ncart + i];
        spherical[7 * nspherical + i] += 4.9607837082461073 * cart[23 * ncart + i];
    }
    // R_64s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i] = -1.9843134832984430 * cart[ncart + i];
        spherical[8 * nspherical + i] += 1.9843134832984430 * cart[15 * ncart + i];
        spherical[8 * nspherical + i] += 19.8431348329844290 * cart[8 * ncart + i];
        spherical[8 * nspherical + i] += -19.8431348329844290 * cart[17 * ncart + i];
    }

    // R_65c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[9 * nspherical + i] = 2.3268138086232857 * cart[2 * ncart + i];
        spherical[9 * nspherical + i] += -23.2681380862328560 * cart[7 * ncart + i];
        spherical[9 * nspherical + i] += 11.6340690431164280 * cart[16 * ncart + i];
    }
    // R_65s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[10 * nspherical + i] = 11.6340690431164280 * cart[4 * ncart + i];
        spherical[10 * nspherical + i] += -23.2681380862328560 * cart[11 * ncart + i];
        spherical[10 * nspherical + i] += 2.3268138086232857 * cart[22 * ncart + i];
    }

    // R_66c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[11 * nspherical + i] = 0.6716932893813962 * cart[i];
        spherical[11 * nspherical + i] += -10.0753993407209421 * cart[3 * ncart + i];
        spherical[11 * nspherical + i] += 10.0753993407209421 * cart[10 * ncart + i];
        spherical[11 * nspherical + i] += -0.6716932893813962 * cart[21 * ncart + i];
    }
    // R_66s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[12 * nspherical + i] = 4.0301597362883772 * cart[ncart + i];
        spherical[12 * nspherical + i] += -13.4338657876279228 * cart[6 * ncart + i];
        spherical[12 * nspherical + i] += 4.0301597362883772 * cart[15 * ncart + i];
    }
}
void gg_gaussian_cart_to_spherical_sum_L6(const unsigned long size, const double* vector,
                                          const double* PRAGMA_RESTRICT cart, const unsigned long ncart,
                                          double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_60 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.3125000000000000 * cart[i];
        tmp += -0.9375000000000000 * cart[3 * ncart + i];
        tmp += -0.9375000000000000 * cart[10 * ncart + i];
        tmp += -0.3125000000000000 * cart[21 * ncart + i];
        tmp += 5.6250000000000000 * cart[5 * ncart + i];
        tmp += 11.2500000000000000 * cart[12 * ncart + i];
        tmp += 5.6250000000000000 * cart[23 * ncart + i];
        tmp += -7.5000000000000000 * cart[14 * ncart + i];
        tmp += -7.5000000000000000 * cart[25 * ncart + i];
        tmp += cart[27 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_61c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.8641098093473998 * cart[2 * ncart + i];
        tmp += 5.7282196186947996 * cart[7 * ncart + i];
        tmp += 2.8641098093473998 * cart[16 * ncart + i];
        tmp += -11.4564392373895991 * cart[9 * ncart + i];
        tmp += -11.4564392373895991 * cart[18 * ncart + i];
        tmp += 4.5825756949558398 * cart[20 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_61s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.8641098093473998 * cart[4 * ncart + i];
        tmp += 5.7282196186947996 * cart[11 * ncart + i];
        tmp += 2.8641098093473998 * cart[22 * ncart + i];
        tmp += -11.4564392373895991 * cart[13 * ncart + i];
        tmp += -11.4564392373895991 * cart[24 * ncart + i];
        tmp += 4.5825756949558398 * cart[26 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_62c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.4528555233184199 * cart[i];
        tmp += 0.4528555233184199 * cart[3 * ncart + i];
        tmp += -0.4528555233184199 * cart[10 * ncart + i];
        tmp += -0.4528555233184199 * cart[21 * ncart + i];
        tmp += -7.2456883730947190 * cart[5 * ncart + i];
        tmp += 7.2456883730947190 * cart[23 * ncart + i];
        tmp += 7.2456883730947190 * cart[14 * ncart + i];
        tmp += -7.2456883730947190 * cart[25 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_62s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.9057110466368399 * cart[ncart + i];
        tmp += 1.8114220932736798 * cart[6 * ncart + i];
        tmp += 0.9057110466368399 * cart[15 * ncart + i];
        tmp += -14.4913767461894381 * cart[8 * ncart + i];
        tmp += -14.4913767461894381 * cart[17 * ncart + i];
        tmp += 14.4913767461894381 * cart[19 * ncart + i];
        output[i] += tmp * vector[4];
    }

    // R_63c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.7171331399105196 * cart[2 * ncart + i];
        tmp += 5.4342662798210393 * cart[7 * ncart + i];
        tmp += 8.1513994197315593 * cart[16 * ncart + i];
        tmp += 7.2456883730947190 * cart[9 * ncart + i];
        tmp += -21.7370651192841571 * cart[18 * ncart + i];
        output[i] += tmp * vector[5];
    }
    // R_63s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -8.1513994197315593 * cart[4 * ncart + i];
        tmp += -5.4342662798210393 * cart[11 * ncart + i];
        tmp += 2.7171331399105196 * cart[22 * ncart + i];
        tmp += 21.7370651192841571 * cart[13 * ncart + i];
        tmp += -7.2456883730947190 * cart[24 * ncart + i];
        output[i] += tmp * vector[6];
    }

    // R_64c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.4960783708246108 * cart[i];
        tmp += 2.4803918541230536 * cart[3 * ncart + i];
        tmp += 2.4803918541230536 * cart[10 * ncart + i];
        tmp += -0.4960783708246108 * cart[21 * ncart + i];
        tmp += 4.9607837082461073 * cart[5 * ncart + i];
        tmp += -29.7647022494766453 * cart[12 * ncart + i];
        tmp += 4.9607837082461073 * cart[23 * ncart + i];
        output[i] += tmp * vector[7];
    }
    // R_64s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -1.9843134832984430 * cart[ncart + i];
        tmp += 1.9843134832984430 * cart[15 * ncart + i];
        tmp += 19.8431348329844290 * cart[8 * ncart + i];
        tmp += -19.8431348329844290 * cart[17 * ncart + i];
        output[i] += tmp * vector[8];
    }

    // R_65c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.3268138086232857 * cart[2 * ncart + i];
        tmp += -23.2681380862328560 * cart[7 * ncart + i];
        tmp += 11.6340690431164280 * cart[16 * ncart + i];
        output[i] += tmp * vector[9];
    }
    // R_65s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 11.6340690431164280 * cart[4 * ncart + i];
        tmp += -23.2681380862328560 * cart[11 * ncart + i];
        tmp += 2.3268138086232857 * cart[22 * ncart + i];
        output[i] += tmp * vector[10];
    }

    // R_66c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.6716932893813962 * cart[i];
        tmp += -10.0753993407209421 * cart[3 * ncart + i];
        tmp += 10.0753993407209421 * cart[10 * ncart + i];
        tmp += -0.6716932893813962 * cart[21 * ncart + i];
        output[i] += tmp * vector[11];
    }
    // R_66s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 4.0301597362883772 * cart[ncart + i];
        tmp += -13.4338657876279228 * cart[6 * ncart + i];
        tmp += 4.0301597362883772 * cart[15 * ncart + i];
        output[i] += tmp * vector[12];
    }
}
void gg_gaussian_cart_to_spherical_L7(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                      const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                      const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_70 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = -2.1875000000000000 * cart[2 * ncart + i];
        spherical[i] += -6.5625000000000000 * cart[7 * ncart + i];
        spherical[i] += -6.5625000000000000 * cart[16 * ncart + i];
        spherical[i] += -2.1875000000000000 * cart[29 * ncart + i];
        spherical[i] += 13.1250000000000000 * cart[9 * ncart + i];
        spherical[i] += 26.2500000000000000 * cart[18 * ncart + i];
        spherical[i] += 13.1250000000000000 * cart[31 * ncart + i];
        spherical[i] += -10.5000000000000000 * cart[20 * ncart + i];
        spherical[i] += -10.5000000000000000 * cart[33 * ncart + i];
        spherical[i] += cart[35 * ncart + i];
    }

    // R_71c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = -0.4133986423538423 * cart[i];
        spherical[nspherical + i] += -1.2401959270615268 * cart[3 * ncart + i];
        spherical[nspherical + i] += -1.2401959270615268 * cart[10 * ncart + i];
        spherical[nspherical + i] += -0.4133986423538423 * cart[21 * ncart + i];
        spherical[nspherical + i] += 9.9215674164922145 * cart[5 * ncart + i];
        spherical[nspherical + i] += 19.8431348329844290 * cart[12 * ncart + i];
        spherical[nspherical + i] += 9.9215674164922145 * cart[23 * ncart + i];
        spherical[nspherical + i] += -19.8431348329844290 * cart[14 * ncart + i];
        spherical[nspherical + i] += -19.8431348329844290 * cart[25 * ncart + i];
        spherical[nspherical + i] += 5.2915026221291814 * cart[27 * ncart + i];
    }
    // R_71s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = -0.4133986423538423 * cart[ncart + i];
        spherical[2 * nspherical + i] += -1.2401959270615268 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] += -1.2401959270615268 * cart[15 * ncart + i];
        spherical[2 * nspherical + i] += -0.4133986423538423 * cart[28 * ncart + i];
        spherical[2 * nspherical + i] += 9.9215674164922145 * cart[8 * ncart + i];
        spherical[2 * nspherical + i] += 19.8431348329844290 * cart[17 * ncart + i];
        spherical[2 * nspherical + i] += 9.9215674164922145 * cart[30 * ncart + i];
        spherical[2 * nspherical + i] += -19.8431348329844290 * cart[19 * ncart + i];
        spherical[2 * nspherical + i] += -19.8431348329844290 * cart[32 * ncart + i];
        spherical[2 * nspherical + i] += 5.2915026221291814 * cart[34 * ncart + i];
    }

    // R_72c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = 3.0378472023786847 * cart[2 * ncart + i];
        spherical[3 * nspherical + i] += 3.0378472023786847 * cart[7 * ncart + i];
        spherical[3 * nspherical + i] += -3.0378472023786847 * cart[16 * ncart + i];
        spherical[3 * nspherical + i] += -3.0378472023786847 * cart[29 * ncart + i];
        spherical[3 * nspherical + i] += -16.2018517460196492 * cart[9 * ncart + i];
        spherical[3 * nspherical + i] += 16.2018517460196492 * cart[31 * ncart + i];
        spherical[3 * nspherical + i] += 9.7211110476117906 * cart[20 * ncart + i];
        spherical[3 * nspherical + i] += -9.7211110476117906 * cart[33 * ncart + i];
    }
    // R_72s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = 6.0756944047573693 * cart[4 * ncart + i];
        spherical[4 * nspherical + i] += 12.1513888095147387 * cart[11 * ncart + i];
        spherical[4 * nspherical + i] += 6.0756944047573693 * cart[22 * ncart + i];
        spherical[4 * nspherical + i] += -32.4037034920392983 * cart[13 * ncart + i];
        spherical[4 * nspherical + i] += -32.4037034920392983 * cart[24 * ncart + i];
        spherical[4 * nspherical + i] += 19.4422220952235811 * cart[26 * ncart + i];
    }

    // R_73c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i] = 0.4296164714021100 * cart[i];
        spherical[5 * nspherical + i] += -0.4296164714021100 * cart[3 * ncart + i];
        spherical[5 * nspherical + i] += -2.1480823570105501 * cart[10 * ncart + i];
        spherical[5 * nspherical + i] += -1.2888494142063300 * cart[21 * ncart + i];
        spherical[5 * nspherical + i] += -8.5923294280422002 * cart[5 * ncart + i];
        spherical[5 * nspherical + i] += 17.1846588560844005 * cart[12 * ncart + i];
        spherical[5 * nspherical + i] += 25.7769882841265989 * cart[23 * ncart + i];
        spherical[5 * nspherical + i] += 11.4564392373895991 * cart[14 * ncart + i];
        spherical[5 * nspherical + i] += -34.3693177121688009 * cart[25 * ncart + i];
    }
    // R_73s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i] = 1.2888494142063300 * cart[ncart + i];
        spherical[6 * nspherical + i] += 2.1480823570105501 * cart[6 * ncart + i];
        spherical[6 * nspherical + i] += 0.4296164714021100 * cart[15 * ncart + i];
        spherical[6 * nspherical + i] += -0.4296164714021100 * cart[28 * ncart + i];
        spherical[6 * nspherical + i] += -25.7769882841265989 * cart[8 * ncart + i];
        spherical[6 * nspherical + i] += -17.1846588560844005 * cart[17 * ncart + i];
        spherical[6 * nspherical + i] += 8.5923294280422002 * cart[30 * ncart + i];
        spherical[6 * nspherical + i] += 34.3693177121688009 * cart[19 * ncart + i];
        spherical[6 * nspherical + i] += -11.4564392373895991 * cart[32 * ncart + i];
    }

    // R_74c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i] = -2.8497532787944992 * cart[2 * ncart + i];
        spherical[7 * nspherical + i] += 14.2487663939724971 * cart[7 * ncart + i];
        spherical[7 * nspherical + i] += 14.2487663939724971 * cart[16 * ncart + i];
        spherical[7 * nspherical + i] += -2.8497532787944992 * cart[29 * ncart + i];
        spherical[7 * nspherical + i] += 9.4991775959816653 * cart[9 * ncart + i];
        spherical[7 * nspherical + i] += -56.9950655758899885 * cart[18 * ncart + i];
        spherical[7 * nspherical + i] += 9.4991775959816653 * cart[31 * ncart + i];
    }
    // R_74s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i] = -11.3990131151779970 * cart[4 * ncart + i];
        spherical[8 * nspherical + i] += 11.3990131151779970 * cart[22 * ncart + i];
        spherical[8 * nspherical + i] += 37.9967103839266613 * cart[13 * ncart + i];
        spherical[8 * nspherical + i] += -37.9967103839266613 * cart[24 * ncart + i];
    }

    // R_75c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[9 * nspherical + i] = -0.4749588797990832 * cart[i];
        spherical[9 * nspherical + i] += 4.2746299181917493 * cart[3 * ncart + i];
        spherical[9 * nspherical + i] += 2.3747943989954163 * cart[10 * ncart + i];
        spherical[9 * nspherical + i] += -2.3747943989954163 * cart[21 * ncart + i];
        spherical[9 * nspherical + i] += 5.6995065575889985 * cart[5 * ncart + i];
        spherical[9 * nspherical + i] += -56.9950655758899885 * cart[12 * ncart + i];
        spherical[9 * nspherical + i] += 28.4975327879449942 * cart[23 * ncart + i];
    }
    // R_75s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[10 * nspherical + i] = -2.3747943989954163 * cart[ncart + i];
        spherical[10 * nspherical + i] += 2.3747943989954163 * cart[6 * ncart + i];
        spherical[10 * nspherical + i] += 4.2746299181917493 * cart[15 * ncart + i];
        spherical[10 * nspherical + i] += -0.4749588797990832 * cart[28 * ncart + i];
        spherical[10 * nspherical + i] += 28.4975327879449942 * cart[8 * ncart + i];
        spherical[10 * nspherical + i] += -56.9950655758899885 * cart[17 * ncart + i];
        spherical[10 * nspherical + i] += 5.6995065575889985 * cart[30 * ncart + i];
    }

    // R_76c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[11 * nspherical + i] = 2.4218245962496954 * cart[2 * ncart + i];
        spherical[11 * nspherical + i] += -36.3273689437454337 * cart[7 * ncart + i];
        spherical[11 * nspherical + i] += 36.3273689437454337 * cart[16 * ncart + i];
        spherical[11 * nspherical + i] += -2.4218245962496954 * cart[29 * ncart + i];
    }
    // R_76s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[12 * nspherical + i] = 14.5309475774981713 * cart[4 * ncart + i];
        spherical[12 * nspherical + i] += -48.4364919249939092 * cart[11 * ncart + i];
        spherical[12 * nspherical + i] += 14.5309475774981713 * cart[22 * ncart + i];
    }

    // R_77c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[13 * nspherical + i] = 0.6472598492877494 * cart[i];
        spherical[13 * nspherical + i] += -13.5924568350427357 * cart[3 * ncart + i];
        spherical[13 * nspherical + i] += 22.6540947250712286 * cart[10 * ncart + i];
        spherical[13 * nspherical + i] += -4.5308189450142455 * cart[21 * ncart + i];
    }
    // R_77s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[14 * nspherical + i] = 4.5308189450142455 * cart[ncart + i];
        spherical[14 * nspherical + i] += -22.6540947250712286 * cart[6 * ncart + i];
        spherical[14 * nspherical + i] += 13.5924568350427357 * cart[15 * ncart + i];
        spherical[14 * nspherical + i] += -0.6472598492877494 * cart[28 * ncart + i];
    }
}
void gg_gaussian_cart_to_spherical_sum_L7(const unsigned long size, const double* vector,
                                          const double* PRAGMA_RESTRICT cart, const unsigned long ncart,
                                          double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_70 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.1875000000000000 * cart[2 * ncart + i];
        tmp += -6.5625000000000000 * cart[7 * ncart + i];
        tmp += -6.5625000000000000 * cart[16 * ncart + i];
        tmp += -2.1875000000000000 * cart[29 * ncart + i];
        tmp += 13.1250000000000000 * cart[9 * ncart + i];
        tmp += 26.2500000000000000 * cart[18 * ncart + i];
        tmp += 13.1250000000000000 * cart[31 * ncart + i];
        tmp += -10.5000000000000000 * cart[20 * ncart + i];
        tmp += -10.5000000000000000 * cart[33 * ncart + i];
        tmp += cart[35 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_71c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.4133986423538423 * cart[i];
        tmp += -1.2401959270615268 * cart[3 * ncart + i];
        tmp += -1.2401959270615268 * cart[10 * ncart + i];
        tmp += -0.4133986423538423 * cart[21 * ncart + i];
        tmp += 9.9215674164922145 * cart[5 * ncart + i];
        tmp += 19.8431348329844290 * cart[12 * ncart + i];
        tmp += 9.9215674164922145 * cart[23 * ncart + i];
        tmp += -19.8431348329844290 * cart[14 * ncart + i];
        tmp += -19.8431348329844290 * cart[25 * ncart + i];
        tmp += 5.2915026221291814 * cart[27 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_71s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.4133986423538423 * cart[ncart + i];
        tmp += -1.2401959270615268 * cart[6 * ncart + i];
        tmp += -1.2401959270615268 * cart[15 * ncart + i];
        tmp += -0.4133986423538423 * cart[28 * ncart + i];
        tmp += 9.9215674164922145 * cart[8 * ncart + i];
        tmp += 19.8431348329844290 * cart[17 * ncart + i];
        tmp += 9.9215674164922145 * cart[30 * ncart + i];
        tmp += -19.8431348329844290 * cart[19 * ncart + i];
        tmp += -19.8431348329844290 * cart[32 * ncart + i];
        tmp += 5.2915026221291814 * cart[34 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_72c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 3.0378472023786847 * cart[2 * ncart + i];
        tmp += 3.0378472023786847 * cart[7 * ncart + i];
        tmp += -3.0378472023786847 * cart[16 * ncart + i];
        tmp += -3.0378472023786847 * cart[29 * ncart + i];
        tmp += -16.2018517460196492 * cart[9 * ncart + i];
        tmp += 16.2018517460196492 * cart[31 * ncart + i];
        tmp += 9.7211110476117906 * cart[20 * ncart + i];
        tmp += -9.7211110476117906 * cart[33 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_72s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 6.0756944047573693 * cart[4 * ncart + i];
        tmp += 12.1513888095147387 * cart[11 * ncart + i];
        tmp += 6.0756944047573693 * cart[22 * ncart + i];
        tmp += -32.4037034920392983 * cart[13 * ncart + i];
        tmp += -32.4037034920392983 * cart[24 * ncart + i];
        tmp += 19.4422220952235811 * cart[26 * ncart + i];
        output[i] += tmp * vector[4];
    }

    // R_73c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.4296164714021100 * cart[i];
        tmp += -0.4296164714021100 * cart[3 * ncart + i];
        tmp += -2.1480823570105501 * cart[10 * ncart + i];
        tmp += -1.2888494142063300 * cart[21 * ncart + i];
        tmp += -8.5923294280422002 * cart[5 * ncart + i];
        tmp += 17.1846588560844005 * cart[12 * ncart + i];
        tmp += 25.7769882841265989 * cart[23 * ncart + i];
        tmp += 11.4564392373895991 * cart[14 * ncart + i];
        tmp += -34.3693177121688009 * cart[25 * ncart + i];
        output[i] += tmp * vector[5];
    }
    // R_73s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.2888494142063300 * cart[ncart + i];
        tmp += 2.1480823570105501 * cart[6 * ncart + i];
        tmp += 0.4296164714021100 * cart[15 * ncart + i];
        tmp += -0.4296164714021100 * cart[28 * ncart + i];
        tmp += -25.7769882841265989 * cart[8 * ncart + i];
        tmp += -17.1846588560844005 * cart[17 * ncart + i];
        tmp += 8.5923294280422002 * cart[30 * ncart + i];
        tmp += 34.3693177121688009 * cart[19 * ncart + i];
        tmp += -11.4564392373895991 * cart[32 * ncart + i];
        output[i] += tmp * vector[6];
    }

    // R_74c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.8497532787944992 * cart[2 * ncart + i];
        tmp += 14.2487663939724971 * cart[7 * ncart + i];
        tmp += 14.2487663939724971 * cart[16 * ncart + i];
        tmp += -2.8497532787944992 * cart[29 * ncart + i];
        tmp += 9.4991775959816653 * cart[9 * ncart + i];
        tmp += -56.9950655758899885 * cart[18 * ncart + i];
        tmp += 9.4991775959816653 * cart[31 * ncart + i];
        output[i] += tmp * vector[7];
    }
    // R_74s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -11.3990131151779970 * cart[4 * ncart + i];
        tmp += 11.3990131151779970 * cart[22 * ncart + i];
        tmp += 37.9967103839266613 * cart[13 * ncart + i];
        tmp += -37.9967103839266613 * cart[24 * ncart + i];
        output[i] += tmp * vector[8];
    }

    // R_75c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.4749588797990832 * cart[i];
        tmp += 4.2746299181917493 * cart[3 * ncart + i];
        tmp += 2.3747943989954163 * cart[10 * ncart + i];
        tmp += -2.3747943989954163 * cart[21 * ncart + i];
        tmp += 5.6995065575889985 * cart[5 * ncart + i];
        tmp += -56.9950655758899885 * cart[12 * ncart + i];
        tmp += 28.4975327879449942 * cart[23 * ncart + i];
        output[i] += tmp * vector[9];
    }
    // R_75s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.3747943989954163 * cart[ncart + i];
        tmp += 2.3747943989954163 * cart[6 * ncart + i];
        tmp += 4.2746299181917493 * cart[15 * ncart + i];
        tmp += -0.4749588797990832 * cart[28 * ncart + i];
        tmp += 28.4975327879449942 * cart[8 * ncart + i];
        tmp += -56.9950655758899885 * cart[17 * ncart + i];
        tmp += 5.6995065575889985 * cart[30 * ncart + i];
        output[i] += tmp * vector[10];
    }

    // R_76c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.4218245962496954 * cart[2 * ncart + i];
        tmp += -36.3273689437454337 * cart[7 * ncart + i];
        tmp += 36.3273689437454337 * cart[16 * ncart + i];
        tmp += -2.4218245962496954 * cart[29 * ncart + i];
        output[i] += tmp * vector[11];
    }
    // R_76s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 14.5309475774981713 * cart[4 * ncart + i];
        tmp += -48.4364919249939092 * cart[11 * ncart + i];
        tmp += 14.5309475774981713 * cart[22 * ncart + i];
        output[i] += tmp * vector[12];
    }

    // R_77c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.6472598492877494 * cart[i];
        tmp += -13.5924568350427357 * cart[3 * ncart + i];
        tmp += 22.6540947250712286 * cart[10 * ncart + i];
        tmp += -4.5308189450142455 * cart[21 * ncart + i];
        output[i] += tmp * vector[13];
    }
    // R_77s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 4.5308189450142455 * cart[ncart + i];
        tmp += -22.6540947250712286 * cart[6 * ncart + i];
        tmp += 13.5924568350427357 * cart[15 * ncart + i];
        tmp += -0.6472598492877494 * cart[28 * ncart + i];
        output[i] += tmp * vector[14];
    }
}
void gg_gaussian_cart_to_spherical_L8(const unsigned long size, const double* PRAGMA_RESTRICT cart,
                                      const unsigned long ncart, double* PRAGMA_RESTRICT spherical,
                                      const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_80 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i] = 0.2734375000000000 * cart[i];
        spherical[i] += 1.0937500000000000 * cart[3 * ncart + i];
        spherical[i] += 1.6406250000000000 * cart[10 * ncart + i];
        spherical[i] += 1.0937500000000000 * cart[21 * ncart + i];
        spherical[i] += 0.2734375000000000 * cart[36 * ncart + i];
        spherical[i] += -8.7500000000000000 * cart[5 * ncart + i];
        spherical[i] += -26.2500000000000000 * cart[12 * ncart + i];
        spherical[i] += -26.2500000000000000 * cart[23 * ncart + i];
        spherical[i] += -8.7500000000000000 * cart[38 * ncart + i];
        spherical[i] += 26.2500000000000000 * cart[14 * ncart + i];
        spherical[i] += 52.5000000000000000 * cart[25 * ncart + i];
        spherical[i] += 26.2500000000000000 * cart[40 * ncart + i];
        spherical[i] += -14.0000000000000000 * cart[27 * ncart + i];
        spherical[i] += -14.0000000000000000 * cart[42 * ncart + i];
        spherical[i] += cart[44 * ncart + i];
    }

    // R_81c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i] = -3.2812500000000000 * cart[2 * ncart + i];
        spherical[nspherical + i] += -9.8437500000000000 * cart[7 * ncart + i];
        spherical[nspherical + i] += -9.8437500000000000 * cart[16 * ncart + i];
        spherical[nspherical + i] += -3.2812500000000000 * cart[29 * ncart + i];
        spherical[nspherical + i] += 26.2500000000000000 * cart[9 * ncart + i];
        spherical[nspherical + i] += 52.5000000000000000 * cart[18 * ncart + i];
        spherical[nspherical + i] += 26.2500000000000000 * cart[31 * ncart + i];
        spherical[nspherical + i] += -31.5000000000000000 * cart[20 * ncart + i];
        spherical[nspherical + i] += -31.5000000000000000 * cart[33 * ncart + i];
        spherical[nspherical + i] += 6.0000000000000000 * cart[35 * ncart + i];
    }
    // R_81s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i] = -3.2812500000000000 * cart[4 * ncart + i];
        spherical[2 * nspherical + i] += -9.8437500000000000 * cart[11 * ncart + i];
        spherical[2 * nspherical + i] += -9.8437500000000000 * cart[22 * ncart + i];
        spherical[2 * nspherical + i] += -3.2812500000000000 * cart[37 * ncart + i];
        spherical[2 * nspherical + i] += 26.2500000000000000 * cart[13 * ncart + i];
        spherical[2 * nspherical + i] += 52.5000000000000000 * cart[24 * ncart + i];
        spherical[2 * nspherical + i] += 26.2500000000000000 * cart[39 * ncart + i];
        spherical[2 * nspherical + i] += -31.5000000000000000 * cart[26 * ncart + i];
        spherical[2 * nspherical + i] += -31.5000000000000000 * cart[41 * ncart + i];
        spherical[2 * nspherical + i] += 6.0000000000000000 * cart[43 * ncart + i];
    }

    // R_82c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i] = -0.3921843874378479 * cart[i];
        spherical[3 * nspherical + i] += -0.7843687748756958 * cart[3 * ncart + i];
        spherical[3 * nspherical + i] += 0.7843687748756958 * cart[21 * ncart + i];
        spherical[3 * nspherical + i] += 0.3921843874378479 * cart[36 * ncart + i];
        spherical[3 * nspherical + i] += 11.7655316231354377 * cart[5 * ncart + i];
        spherical[3 * nspherical + i] += 11.7655316231354377 * cart[12 * ncart + i];
        spherical[3 * nspherical + i] += -11.7655316231354377 * cart[23 * ncart + i];
        spherical[3 * nspherical + i] += -11.7655316231354377 * cart[38 * ncart + i];
        spherical[3 * nspherical + i] += -31.3747509950278314 * cart[14 * ncart + i];
        spherical[3 * nspherical + i] += 31.3747509950278314 * cart[40 * ncart + i];
        spherical[3 * nspherical + i] += 12.5499003980111326 * cart[27 * ncart + i];
        spherical[3 * nspherical + i] += -12.5499003980111326 * cart[42 * ncart + i];
    }
    // R_82s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i] = -0.7843687748756958 * cart[ncart + i];
        spherical[4 * nspherical + i] += -2.3531063246270874 * cart[6 * ncart + i];
        spherical[4 * nspherical + i] += -2.3531063246270874 * cart[15 * ncart + i];
        spherical[4 * nspherical + i] += -0.7843687748756958 * cart[28 * ncart + i];
        spherical[4 * nspherical + i] += 23.5310632462708753 * cart[8 * ncart + i];
        spherical[4 * nspherical + i] += 47.0621264925417506 * cart[17 * ncart + i];
        spherical[4 * nspherical + i] += 23.5310632462708753 * cart[30 * ncart + i];
        spherical[4 * nspherical + i] += -62.7495019900556628 * cart[19 * ncart + i];
        spherical[4 * nspherical + i] += -62.7495019900556628 * cart[32 * ncart + i];
        spherical[4 * nspherical + i] += 25.0998007960222651 * cart[34 * ncart + i];
    }

    // R_83c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i] = 3.1861210252437053 * cart[2 * ncart + i];
        spherical[5 * nspherical + i] += -3.1861210252437053 * cart[7 * ncart + i];
        spherical[5 * nspherical + i] += -15.9306051262185271 * cart[16 * ncart + i];
        spherical[5 * nspherical + i] += -9.5583630757311155 * cart[29 * ncart + i];
        spherical[5 * nspherical + i] += -21.2408068349580361 * cart[9 * ncart + i];
        spherical[5 * nspherical + i] += 42.4816136699160722 * cart[18 * ncart + i];
        spherical[5 * nspherical + i] += 63.7224205048741084 * cart[31 * ncart + i];
        spherical[5 * nspherical + i] += 16.9926454679664296 * cart[20 * ncart + i];
        spherical[5 * nspherical + i] += -50.9779364038992853 * cart[33 * ncart + i];
    }
    // R_83s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i] = 9.5583630757311155 * cart[4 * ncart + i];
        spherical[6 * nspherical + i] += 15.9306051262185271 * cart[11 * ncart + i];
        spherical[6 * nspherical + i] += 3.1861210252437053 * cart[22 * ncart + i];
        spherical[6 * nspherical + i] += -3.1861210252437053 * cart[37 * ncart + i];
        spherical[6 * nspherical + i] += -63.7224205048741084 * cart[13 * ncart + i];
        spherical[6 * nspherical + i] += -42.4816136699160722 * cart[24 * ncart + i];
        spherical[6 * nspherical + i] += 21.2408068349580361 * cart[39 * ncart + i];
        spherical[6 * nspherical + i] += 50.9779364038992853 * cart[26 * ncart + i];
        spherical[6 * nspherical + i] += -16.9926454679664296 * cart[41 * ncart + i];
    }

    // R_84c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i] = 0.4113264556590057 * cart[i];
        spherical[7 * nspherical + i] += -1.6453058226360229 * cart[3 * ncart + i];
        spherical[7 * nspherical + i] += -4.1132645565900576 * cart[10 * ncart + i];
        spherical[7 * nspherical + i] += -1.6453058226360229 * cart[21 * ncart + i];
        spherical[7 * nspherical + i] += 0.4113264556590057 * cart[36 * ncart + i];
        spherical[7 * nspherical + i] += -9.8718349358161372 * cart[5 * ncart + i];
        spherical[7 * nspherical + i] += 49.3591746790806880 * cart[12 * ncart + i];
        spherical[7 * nspherical + i] += 49.3591746790806880 * cart[23 * ncart + i];
        spherical[7 * nspherical + i] += -9.8718349358161372 * cart[38 * ncart + i];
        spherical[7 * nspherical + i] += 16.4530582263602305 * cart[14 * ncart + i];
        spherical[7 * nspherical + i] += -98.7183493581613760 * cart[25 * ncart + i];
        spherical[7 * nspherical + i] += 16.4530582263602305 * cart[40 * ncart + i];
    }
    // R_84s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i] = 1.6453058226360229 * cart[ncart + i];
        spherical[8 * nspherical + i] += 1.6453058226360229 * cart[6 * ncart + i];
        spherical[8 * nspherical + i] += -1.6453058226360229 * cart[15 * ncart + i];
        spherical[8 * nspherical + i] += -1.6453058226360229 * cart[28 * ncart + i];
        spherical[8 * nspherical + i] += -39.4873397432645490 * cart[8 * ncart + i];
        spherical[8 * nspherical + i] += 39.4873397432645490 * cart[30 * ncart + i];
        spherical[8 * nspherical + i] += 65.8122329054409221 * cart[19 * ncart + i];
        spherical[8 * nspherical + i] += -65.8122329054409221 * cart[32 * ncart + i];
    }

    // R_85c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[9 * nspherical + i] = -2.9661172536668201 * cart[2 * ncart + i];
        spherical[9 * nspherical + i] += 26.6950552830013805 * cart[7 * ncart + i];
        spherical[9 * nspherical + i] += 14.8305862683341019 * cart[16 * ncart + i];
        spherical[9 * nspherical + i] += -14.8305862683341019 * cart[29 * ncart + i];
        spherical[9 * nspherical + i] += 11.8644690146672804 * cart[9 * ncart + i];
        spherical[9 * nspherical + i] += -118.6446901466728150 * cart[18 * ncart + i];
        spherical[9 * nspherical + i] += 59.3223450733364075 * cart[31 * ncart + i];
    }
    // R_85s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[10 * nspherical + i] = -14.8305862683341019 * cart[4 * ncart + i];
        spherical[10 * nspherical + i] += 14.8305862683341019 * cart[11 * ncart + i];
        spherical[10 * nspherical + i] += 26.6950552830013805 * cart[22 * ncart + i];
        spherical[10 * nspherical + i] += -2.9661172536668201 * cart[37 * ncart + i];
        spherical[10 * nspherical + i] += 59.3223450733364075 * cart[13 * ncart + i];
        spherical[10 * nspherical + i] += -118.6446901466728150 * cart[24 * ncart + i];
        spherical[10 * nspherical + i] += 11.8644690146672804 * cart[39 * ncart + i];
    }

    // R_86c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[11 * nspherical + i] = -0.4576818286211503 * cart[i];
        spherical[11 * nspherical + i] += 6.4075456006961042 * cart[3 * ncart + i];
        spherical[11 * nspherical + i] += -6.4075456006961042 * cart[21 * ncart + i];
        spherical[11 * nspherical + i] += 0.4576818286211503 * cart[36 * ncart + i];
        spherical[11 * nspherical + i] += 6.4075456006961042 * cart[5 * ncart + i];
        spherical[11 * nspherical + i] += -96.1131840104415573 * cart[12 * ncart + i];
        spherical[11 * nspherical + i] += 96.1131840104415573 * cart[23 * ncart + i];
        spherical[11 * nspherical + i] += -6.4075456006961042 * cart[38 * ncart + i];
    }
    // R_86s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[12 * nspherical + i] = -2.7460909717269018 * cart[ncart + i];
        spherical[12 * nspherical + i] += 6.4075456006961042 * cart[6 * ncart + i];
        spherical[12 * nspherical + i] += 6.4075456006961042 * cart[15 * ncart + i];
        spherical[12 * nspherical + i] += -2.7460909717269018 * cart[28 * ncart + i];
        spherical[12 * nspherical + i] += 38.4452736041766272 * cart[8 * ncart + i];
        spherical[12 * nspherical + i] += -128.1509120139220954 * cart[17 * ncart + i];
        spherical[12 * nspherical + i] += 38.4452736041766272 * cart[30 * ncart + i];
    }

    // R_87c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[13 * nspherical + i] = 2.5068266169601756 * cart[2 * ncart + i];
        spherical[13 * nspherical + i] += -52.6433589561636950 * cart[7 * ncart + i];
        spherical[13 * nspherical + i] += 87.7389315936061536 * cart[16 * ncart + i];
        spherical[13 * nspherical + i] += -17.5477863187212293 * cart[29 * ncart + i];
    }
    // R_87s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[14 * nspherical + i] = 17.5477863187212293 * cart[4 * ncart + i];
        spherical[14 * nspherical + i] += -87.7389315936061536 * cart[11 * ncart + i];
        spherical[14 * nspherical + i] += 52.6433589561636950 * cart[22 * ncart + i];
        spherical[14 * nspherical + i] += -2.5068266169601756 * cart[37 * ncart + i];
    }

    // R_88c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[15 * nspherical + i] = 0.6267066542400439 * cart[i];
        spherical[15 * nspherical + i] += -17.5477863187212293 * cart[3 * ncart + i];
        spherical[15 * nspherical + i] += 43.8694657968030768 * cart[10 * ncart + i];
        spherical[15 * nspherical + i] += -17.5477863187212293 * cart[21 * ncart + i];
        spherical[15 * nspherical + i] += 0.6267066542400439 * cart[36 * ncart + i];
    }
    // R_88s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[16 * nspherical + i] = 5.0136532339203512 * cart[ncart + i];
        spherical[16 * nspherical + i] += -35.0955726374424586 * cart[6 * ncart + i];
        spherical[16 * nspherical + i] += 35.0955726374424586 * cart[15 * ncart + i];
        spherical[16 * nspherical + i] += -5.0136532339203512 * cart[28 * ncart + i];
    }
}
void gg_gaussian_cart_to_spherical_sum_L8(const unsigned long size, const double* vector,
                                          const double* PRAGMA_RESTRICT cart, const unsigned long ncart,
                                          double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_80 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.2734375000000000 * cart[i];
        tmp += 1.0937500000000000 * cart[3 * ncart + i];
        tmp += 1.6406250000000000 * cart[10 * ncart + i];
        tmp += 1.0937500000000000 * cart[21 * ncart + i];
        tmp += 0.2734375000000000 * cart[36 * ncart + i];
        tmp += -8.7500000000000000 * cart[5 * ncart + i];
        tmp += -26.2500000000000000 * cart[12 * ncart + i];
        tmp += -26.2500000000000000 * cart[23 * ncart + i];
        tmp += -8.7500000000000000 * cart[38 * ncart + i];
        tmp += 26.2500000000000000 * cart[14 * ncart + i];
        tmp += 52.5000000000000000 * cart[25 * ncart + i];
        tmp += 26.2500000000000000 * cart[40 * ncart + i];
        tmp += -14.0000000000000000 * cart[27 * ncart + i];
        tmp += -14.0000000000000000 * cart[42 * ncart + i];
        tmp += cart[44 * ncart + i];
        output[i] += tmp * vector[0];
    }

    // R_81c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -3.2812500000000000 * cart[2 * ncart + i];
        tmp += -9.8437500000000000 * cart[7 * ncart + i];
        tmp += -9.8437500000000000 * cart[16 * ncart + i];
        tmp += -3.2812500000000000 * cart[29 * ncart + i];
        tmp += 26.2500000000000000 * cart[9 * ncart + i];
        tmp += 52.5000000000000000 * cart[18 * ncart + i];
        tmp += 26.2500000000000000 * cart[31 * ncart + i];
        tmp += -31.5000000000000000 * cart[20 * ncart + i];
        tmp += -31.5000000000000000 * cart[33 * ncart + i];
        tmp += 6.0000000000000000 * cart[35 * ncart + i];
        output[i] += tmp * vector[1];
    }
    // R_81s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -3.2812500000000000 * cart[4 * ncart + i];
        tmp += -9.8437500000000000 * cart[11 * ncart + i];
        tmp += -9.8437500000000000 * cart[22 * ncart + i];
        tmp += -3.2812500000000000 * cart[37 * ncart + i];
        tmp += 26.2500000000000000 * cart[13 * ncart + i];
        tmp += 52.5000000000000000 * cart[24 * ncart + i];
        tmp += 26.2500000000000000 * cart[39 * ncart + i];
        tmp += -31.5000000000000000 * cart[26 * ncart + i];
        tmp += -31.5000000000000000 * cart[41 * ncart + i];
        tmp += 6.0000000000000000 * cart[43 * ncart + i];
        output[i] += tmp * vector[2];
    }

    // R_82c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.3921843874378479 * cart[i];
        tmp += -0.7843687748756958 * cart[3 * ncart + i];
        tmp += 0.7843687748756958 * cart[21 * ncart + i];
        tmp += 0.3921843874378479 * cart[36 * ncart + i];
        tmp += 11.7655316231354377 * cart[5 * ncart + i];
        tmp += 11.7655316231354377 * cart[12 * ncart + i];
        tmp += -11.7655316231354377 * cart[23 * ncart + i];
        tmp += -11.7655316231354377 * cart[38 * ncart + i];
        tmp += -31.3747509950278314 * cart[14 * ncart + i];
        tmp += 31.3747509950278314 * cart[40 * ncart + i];
        tmp += 12.5499003980111326 * cart[27 * ncart + i];
        tmp += -12.5499003980111326 * cart[42 * ncart + i];
        output[i] += tmp * vector[3];
    }
    // R_82s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.7843687748756958 * cart[ncart + i];
        tmp += -2.3531063246270874 * cart[6 * ncart + i];
        tmp += -2.3531063246270874 * cart[15 * ncart + i];
        tmp += -0.7843687748756958 * cart[28 * ncart + i];
        tmp += 23.5310632462708753 * cart[8 * ncart + i];
        tmp += 47.0621264925417506 * cart[17 * ncart + i];
        tmp += 23.5310632462708753 * cart[30 * ncart + i];
        tmp += -62.7495019900556628 * cart[19 * ncart + i];
        tmp += -62.7495019900556628 * cart[32 * ncart + i];
        tmp += 25.0998007960222651 * cart[34 * ncart + i];
        output[i] += tmp * vector[4];
    }

    // R_83c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 3.1861210252437053 * cart[2 * ncart + i];
        tmp += -3.1861210252437053 * cart[7 * ncart + i];
        tmp += -15.9306051262185271 * cart[16 * ncart + i];
        tmp += -9.5583630757311155 * cart[29 * ncart + i];
        tmp += -21.2408068349580361 * cart[9 * ncart + i];
        tmp += 42.4816136699160722 * cart[18 * ncart + i];
        tmp += 63.7224205048741084 * cart[31 * ncart + i];
        tmp += 16.9926454679664296 * cart[20 * ncart + i];
        tmp += -50.9779364038992853 * cart[33 * ncart + i];
        output[i] += tmp * vector[5];
    }
    // R_83s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 9.5583630757311155 * cart[4 * ncart + i];
        tmp += 15.9306051262185271 * cart[11 * ncart + i];
        tmp += 3.1861210252437053 * cart[22 * ncart + i];
        tmp += -3.1861210252437053 * cart[37 * ncart + i];
        tmp += -63.7224205048741084 * cart[13 * ncart + i];
        tmp += -42.4816136699160722 * cart[24 * ncart + i];
        tmp += 21.2408068349580361 * cart[39 * ncart + i];
        tmp += 50.9779364038992853 * cart[26 * ncart + i];
        tmp += -16.9926454679664296 * cart[41 * ncart + i];
        output[i] += tmp * vector[6];
    }

    // R_84c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.4113264556590057 * cart[i];
        tmp += -1.6453058226360229 * cart[3 * ncart + i];
        tmp += -4.1132645565900576 * cart[10 * ncart + i];
        tmp += -1.6453058226360229 * cart[21 * ncart + i];
        tmp += 0.4113264556590057 * cart[36 * ncart + i];
        tmp += -9.8718349358161372 * cart[5 * ncart + i];
        tmp += 49.3591746790806880 * cart[12 * ncart + i];
        tmp += 49.3591746790806880 * cart[23 * ncart + i];
        tmp += -9.8718349358161372 * cart[38 * ncart + i];
        tmp += 16.4530582263602305 * cart[14 * ncart + i];
        tmp += -98.7183493581613760 * cart[25 * ncart + i];
        tmp += 16.4530582263602305 * cart[40 * ncart + i];
        output[i] += tmp * vector[7];
    }
    // R_84s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 1.6453058226360229 * cart[ncart + i];
        tmp += 1.6453058226360229 * cart[6 * ncart + i];
        tmp += -1.6453058226360229 * cart[15 * ncart + i];
        tmp += -1.6453058226360229 * cart[28 * ncart + i];
        tmp += -39.4873397432645490 * cart[8 * ncart + i];
        tmp += 39.4873397432645490 * cart[30 * ncart + i];
        tmp += 65.8122329054409221 * cart[19 * ncart + i];
        tmp += -65.8122329054409221 * cart[32 * ncart + i];
        output[i] += tmp * vector[8];
    }

    // R_85c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.9661172536668201 * cart[2 * ncart + i];
        tmp += 26.6950552830013805 * cart[7 * ncart + i];
        tmp += 14.8305862683341019 * cart[16 * ncart + i];
        tmp += -14.8305862683341019 * cart[29 * ncart + i];
        tmp += 11.8644690146672804 * cart[9 * ncart + i];
        tmp += -118.6446901466728150 * cart[18 * ncart + i];
        tmp += 59.3223450733364075 * cart[31 * ncart + i];
        output[i] += tmp * vector[9];
    }
    // R_85s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -14.8305862683341019 * cart[4 * ncart + i];
        tmp += 14.8305862683341019 * cart[11 * ncart + i];
        tmp += 26.6950552830013805 * cart[22 * ncart + i];
        tmp += -2.9661172536668201 * cart[37 * ncart + i];
        tmp += 59.3223450733364075 * cart[13 * ncart + i];
        tmp += -118.6446901466728150 * cart[24 * ncart + i];
        tmp += 11.8644690146672804 * cart[39 * ncart + i];
        output[i] += tmp * vector[10];
    }

    // R_86c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -0.4576818286211503 * cart[i];
        tmp += 6.4075456006961042 * cart[3 * ncart + i];
        tmp += -6.4075456006961042 * cart[21 * ncart + i];
        tmp += 0.4576818286211503 * cart[36 * ncart + i];
        tmp += 6.4075456006961042 * cart[5 * ncart + i];
        tmp += -96.1131840104415573 * cart[12 * ncart + i];
        tmp += 96.1131840104415573 * cart[23 * ncart + i];
        tmp += -6.4075456006961042 * cart[38 * ncart + i];
        output[i] += tmp * vector[11];
    }
    // R_86s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = -2.7460909717269018 * cart[ncart + i];
        tmp += 6.4075456006961042 * cart[6 * ncart + i];
        tmp += 6.4075456006961042 * cart[15 * ncart + i];
        tmp += -2.7460909717269018 * cart[28 * ncart + i];
        tmp += 38.4452736041766272 * cart[8 * ncart + i];
        tmp += -128.1509120139220954 * cart[17 * ncart + i];
        tmp += 38.4452736041766272 * cart[30 * ncart + i];
        output[i] += tmp * vector[12];
    }

    // R_87c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 2.5068266169601756 * cart[2 * ncart + i];
        tmp += -52.6433589561636950 * cart[7 * ncart + i];
        tmp += 87.7389315936061536 * cart[16 * ncart + i];
        tmp += -17.5477863187212293 * cart[29 * ncart + i];
        output[i] += tmp * vector[13];
    }
    // R_87s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 17.5477863187212293 * cart[4 * ncart + i];
        tmp += -87.7389315936061536 * cart[11 * ncart + i];
        tmp += 52.6433589561636950 * cart[22 * ncart + i];
        tmp += -2.5068266169601756 * cart[37 * ncart + i];
        output[i] += tmp * vector[14];
    }

    // R_88c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 0.6267066542400439 * cart[i];
        tmp += -17.5477863187212293 * cart[3 * ncart + i];
        tmp += 43.8694657968030768 * cart[10 * ncart + i];
        tmp += -17.5477863187212293 * cart[21 * ncart + i];
        tmp += 0.6267066542400439 * cart[36 * ncart + i];
        output[i] += tmp * vector[15];
    }
    // R_88s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp = 5.0136532339203512 * cart[ncart + i];
        tmp += -35.0955726374424586 * cart[6 * ncart + i];
        tmp += 35.0955726374424586 * cart[15 * ncart + i];
        tmp += -5.0136532339203512 * cart[28 * ncart + i];
        output[i] += tmp * vector[16];
    }
}
void gg_cca_cart_copy_L0(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                         const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                         const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (0, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_cca_cart_sum_L0(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                        const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                        double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (0, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_cca_cart_copy_L1(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                         const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                         const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (1, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 1, 0)
    inp_shift = 1 * ncart_input;
    out_shift = 1 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 0, 1)
    inp_shift = 2 * ncart_input;
    out_shift = 2 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_cca_cart_sum_L1(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                        const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                        double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (1, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 1, 0)
    in_shift = 1 * ncart_input;
    coef = vector[1];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 0, 1)
    in_shift = 2 * ncart_input;
    coef = vector[2];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_cca_cart_copy_L2(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                         const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                         const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (2, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 1, 0)
    inp_shift = 1 * ncart_input;
    out_shift = 1 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 0, 1)
    inp_shift = 2 * ncart_input;
    out_shift = 2 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 2, 0)
    inp_shift = 3 * ncart_input;
    out_shift = 3 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 1, 1)
    inp_shift = 4 * ncart_input;
    out_shift = 4 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 0, 2)
    inp_shift = 5 * ncart_input;
    out_shift = 5 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_cca_cart_sum_L2(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                        const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                        double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (2, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 1, 0)
    in_shift = 1 * ncart_input;
    coef = vector[1];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 0, 1)
    in_shift = 2 * ncart_input;
    coef = vector[2];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 2, 0)
    in_shift = 3 * ncart_input;
    coef = vector[3];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 1, 1)
    in_shift = 4 * ncart_input;
    coef = vector[4];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 0, 2)
    in_shift = 5 * ncart_input;
    coef = vector[5];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_cca_cart_copy_L3(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                         const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                         const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (3, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 1, 0)
    inp_shift = 1 * ncart_input;
    out_shift = 1 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 0, 1)
    inp_shift = 2 * ncart_input;
    out_shift = 2 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 2, 0)
    inp_shift = 3 * ncart_input;
    out_shift = 3 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 1, 1)
    inp_shift = 4 * ncart_input;
    out_shift = 4 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 0, 2)
    inp_shift = 5 * ncart_input;
    out_shift = 5 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 3, 0)
    inp_shift = 6 * ncart_input;
    out_shift = 6 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 2, 1)
    inp_shift = 7 * ncart_input;
    out_shift = 7 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 1, 2)
    inp_shift = 8 * ncart_input;
    out_shift = 8 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 0, 3)
    inp_shift = 9 * ncart_input;
    out_shift = 9 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_cca_cart_sum_L3(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                        const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                        double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (3, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 1, 0)
    in_shift = 1 * ncart_input;
    coef = vector[1];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 0, 1)
    in_shift = 2 * ncart_input;
    coef = vector[2];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 2, 0)
    in_shift = 3 * ncart_input;
    coef = vector[3];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 1, 1)
    in_shift = 4 * ncart_input;
    coef = vector[4];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 0, 2)
    in_shift = 5 * ncart_input;
    coef = vector[5];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 3, 0)
    in_shift = 6 * ncart_input;
    coef = vector[6];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 2, 1)
    in_shift = 7 * ncart_input;
    coef = vector[7];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 1, 2)
    in_shift = 8 * ncart_input;
    coef = vector[8];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 0, 3)
    in_shift = 9 * ncart_input;
    coef = vector[9];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_cca_cart_copy_L4(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                         const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                         const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (4, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 1, 0)
    inp_shift = 1 * ncart_input;
    out_shift = 1 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 0, 1)
    inp_shift = 2 * ncart_input;
    out_shift = 2 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 2, 0)
    inp_shift = 3 * ncart_input;
    out_shift = 3 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 1, 1)
    inp_shift = 4 * ncart_input;
    out_shift = 4 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 0, 2)
    inp_shift = 5 * ncart_input;
    out_shift = 5 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 3, 0)
    inp_shift = 6 * ncart_input;
    out_shift = 6 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 2, 1)
    inp_shift = 7 * ncart_input;
    out_shift = 7 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 1, 2)
    inp_shift = 8 * ncart_input;
    out_shift = 8 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 0, 3)
    inp_shift = 9 * ncart_input;
    out_shift = 9 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 4, 0)
    inp_shift = 10 * ncart_input;
    out_shift = 10 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 3, 1)
    inp_shift = 11 * ncart_input;
    out_shift = 11 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 2, 2)
    inp_shift = 12 * ncart_input;
    out_shift = 12 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 1, 3)
    inp_shift = 13 * ncart_input;
    out_shift = 13 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 0, 4)
    inp_shift = 14 * ncart_input;
    out_shift = 14 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_cca_cart_sum_L4(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                        const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                        double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (4, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 1, 0)
    in_shift = 1 * ncart_input;
    coef = vector[1];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 0, 1)
    in_shift = 2 * ncart_input;
    coef = vector[2];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 2, 0)
    in_shift = 3 * ncart_input;
    coef = vector[3];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 1, 1)
    in_shift = 4 * ncart_input;
    coef = vector[4];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 0, 2)
    in_shift = 5 * ncart_input;
    coef = vector[5];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 3, 0)
    in_shift = 6 * ncart_input;
    coef = vector[6];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 2, 1)
    in_shift = 7 * ncart_input;
    coef = vector[7];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 1, 2)
    in_shift = 8 * ncart_input;
    coef = vector[8];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 0, 3)
    in_shift = 9 * ncart_input;
    coef = vector[9];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 4, 0)
    in_shift = 10 * ncart_input;
    coef = vector[10];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 3, 1)
    in_shift = 11 * ncart_input;
    coef = vector[11];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 2, 2)
    in_shift = 12 * ncart_input;
    coef = vector[12];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 1, 3)
    in_shift = 13 * ncart_input;
    coef = vector[13];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 0, 4)
    in_shift = 14 * ncart_input;
    coef = vector[14];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_cca_cart_copy_L5(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                         const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                         const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (5, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 1, 0)
    inp_shift = 1 * ncart_input;
    out_shift = 1 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 0, 1)
    inp_shift = 2 * ncart_input;
    out_shift = 2 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 2, 0)
    inp_shift = 3 * ncart_input;
    out_shift = 3 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 1, 1)
    inp_shift = 4 * ncart_input;
    out_shift = 4 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 0, 2)
    inp_shift = 5 * ncart_input;
    out_shift = 5 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 3, 0)
    inp_shift = 6 * ncart_input;
    out_shift = 6 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 2, 1)
    inp_shift = 7 * ncart_input;
    out_shift = 7 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 1, 2)
    inp_shift = 8 * ncart_input;
    out_shift = 8 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 0, 3)
    inp_shift = 9 * ncart_input;
    out_shift = 9 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 4, 0)
    inp_shift = 10 * ncart_input;
    out_shift = 10 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 3, 1)
    inp_shift = 11 * ncart_input;
    out_shift = 11 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 2, 2)
    inp_shift = 12 * ncart_input;
    out_shift = 12 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 1, 3)
    inp_shift = 13 * ncart_input;
    out_shift = 13 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 0, 4)
    inp_shift = 14 * ncart_input;
    out_shift = 14 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 5, 0)
    inp_shift = 15 * ncart_input;
    out_shift = 15 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 4, 1)
    inp_shift = 16 * ncart_input;
    out_shift = 16 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 3, 2)
    inp_shift = 17 * ncart_input;
    out_shift = 17 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 2, 3)
    inp_shift = 18 * ncart_input;
    out_shift = 18 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 1, 4)
    inp_shift = 19 * ncart_input;
    out_shift = 19 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 0, 5)
    inp_shift = 20 * ncart_input;
    out_shift = 20 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_cca_cart_sum_L5(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                        const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                        double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (5, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 1, 0)
    in_shift = 1 * ncart_input;
    coef = vector[1];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 0, 1)
    in_shift = 2 * ncart_input;
    coef = vector[2];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 2, 0)
    in_shift = 3 * ncart_input;
    coef = vector[3];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 1, 1)
    in_shift = 4 * ncart_input;
    coef = vector[4];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 0, 2)
    in_shift = 5 * ncart_input;
    coef = vector[5];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 3, 0)
    in_shift = 6 * ncart_input;
    coef = vector[6];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 2, 1)
    in_shift = 7 * ncart_input;
    coef = vector[7];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 1, 2)
    in_shift = 8 * ncart_input;
    coef = vector[8];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 0, 3)
    in_shift = 9 * ncart_input;
    coef = vector[9];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 4, 0)
    in_shift = 10 * ncart_input;
    coef = vector[10];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 3, 1)
    in_shift = 11 * ncart_input;
    coef = vector[11];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 2, 2)
    in_shift = 12 * ncart_input;
    coef = vector[12];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 1, 3)
    in_shift = 13 * ncart_input;
    coef = vector[13];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 0, 4)
    in_shift = 14 * ncart_input;
    coef = vector[14];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 5, 0)
    in_shift = 15 * ncart_input;
    coef = vector[15];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 4, 1)
    in_shift = 16 * ncart_input;
    coef = vector[16];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 3, 2)
    in_shift = 17 * ncart_input;
    coef = vector[17];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 2, 3)
    in_shift = 18 * ncart_input;
    coef = vector[18];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 1, 4)
    in_shift = 19 * ncart_input;
    coef = vector[19];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 0, 5)
    in_shift = 20 * ncart_input;
    coef = vector[20];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_cca_cart_copy_L6(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                         const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                         const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (6, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (5, 1, 0)
    inp_shift = 1 * ncart_input;
    out_shift = 1 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (5, 0, 1)
    inp_shift = 2 * ncart_input;
    out_shift = 2 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 2, 0)
    inp_shift = 3 * ncart_input;
    out_shift = 3 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 1, 1)
    inp_shift = 4 * ncart_input;
    out_shift = 4 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 0, 2)
    inp_shift = 5 * ncart_input;
    out_shift = 5 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 3, 0)
    inp_shift = 6 * ncart_input;
    out_shift = 6 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 2, 1)
    inp_shift = 7 * ncart_input;
    out_shift = 7 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 1, 2)
    inp_shift = 8 * ncart_input;
    out_shift = 8 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 0, 3)
    inp_shift = 9 * ncart_input;
    out_shift = 9 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 4, 0)
    inp_shift = 10 * ncart_input;
    out_shift = 10 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 3, 1)
    inp_shift = 11 * ncart_input;
    out_shift = 11 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 2, 2)
    inp_shift = 12 * ncart_input;
    out_shift = 12 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 1, 3)
    inp_shift = 13 * ncart_input;
    out_shift = 13 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 0, 4)
    inp_shift = 14 * ncart_input;
    out_shift = 14 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 5, 0)
    inp_shift = 15 * ncart_input;
    out_shift = 15 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 4, 1)
    inp_shift = 16 * ncart_input;
    out_shift = 16 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 3, 2)
    inp_shift = 17 * ncart_input;
    out_shift = 17 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 2, 3)
    inp_shift = 18 * ncart_input;
    out_shift = 18 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 1, 4)
    inp_shift = 19 * ncart_input;
    out_shift = 19 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 0, 5)
    inp_shift = 20 * ncart_input;
    out_shift = 20 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 6, 0)
    inp_shift = 21 * ncart_input;
    out_shift = 21 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 5, 1)
    inp_shift = 22 * ncart_input;
    out_shift = 22 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 4, 2)
    inp_shift = 23 * ncart_input;
    out_shift = 23 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 3, 3)
    inp_shift = 24 * ncart_input;
    out_shift = 24 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 2, 4)
    inp_shift = 25 * ncart_input;
    out_shift = 25 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 1, 5)
    inp_shift = 26 * ncart_input;
    out_shift = 26 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 0, 6)
    inp_shift = 27 * ncart_input;
    out_shift = 27 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_cca_cart_sum_L6(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                        const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                        double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (6, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (5, 1, 0)
    in_shift = 1 * ncart_input;
    coef = vector[1];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (5, 0, 1)
    in_shift = 2 * ncart_input;
    coef = vector[2];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 2, 0)
    in_shift = 3 * ncart_input;
    coef = vector[3];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 1, 1)
    in_shift = 4 * ncart_input;
    coef = vector[4];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 0, 2)
    in_shift = 5 * ncart_input;
    coef = vector[5];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 3, 0)
    in_shift = 6 * ncart_input;
    coef = vector[6];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 2, 1)
    in_shift = 7 * ncart_input;
    coef = vector[7];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 1, 2)
    in_shift = 8 * ncart_input;
    coef = vector[8];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 0, 3)
    in_shift = 9 * ncart_input;
    coef = vector[9];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 4, 0)
    in_shift = 10 * ncart_input;
    coef = vector[10];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 3, 1)
    in_shift = 11 * ncart_input;
    coef = vector[11];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 2, 2)
    in_shift = 12 * ncart_input;
    coef = vector[12];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 1, 3)
    in_shift = 13 * ncart_input;
    coef = vector[13];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 0, 4)
    in_shift = 14 * ncart_input;
    coef = vector[14];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 5, 0)
    in_shift = 15 * ncart_input;
    coef = vector[15];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 4, 1)
    in_shift = 16 * ncart_input;
    coef = vector[16];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 3, 2)
    in_shift = 17 * ncart_input;
    coef = vector[17];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 2, 3)
    in_shift = 18 * ncart_input;
    coef = vector[18];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 1, 4)
    in_shift = 19 * ncart_input;
    coef = vector[19];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 0, 5)
    in_shift = 20 * ncart_input;
    coef = vector[20];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 6, 0)
    in_shift = 21 * ncart_input;
    coef = vector[21];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 5, 1)
    in_shift = 22 * ncart_input;
    coef = vector[22];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 4, 2)
    in_shift = 23 * ncart_input;
    coef = vector[23];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 3, 3)
    in_shift = 24 * ncart_input;
    coef = vector[24];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 2, 4)
    in_shift = 25 * ncart_input;
    coef = vector[25];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 1, 5)
    in_shift = 26 * ncart_input;
    coef = vector[26];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 0, 6)
    in_shift = 27 * ncart_input;
    coef = vector[27];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_cca_cart_copy_L7(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                         const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                         const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (7, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (6, 1, 0)
    inp_shift = 1 * ncart_input;
    out_shift = 1 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (6, 0, 1)
    inp_shift = 2 * ncart_input;
    out_shift = 2 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (5, 2, 0)
    inp_shift = 3 * ncart_input;
    out_shift = 3 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (5, 1, 1)
    inp_shift = 4 * ncart_input;
    out_shift = 4 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (5, 0, 2)
    inp_shift = 5 * ncart_input;
    out_shift = 5 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 3, 0)
    inp_shift = 6 * ncart_input;
    out_shift = 6 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 2, 1)
    inp_shift = 7 * ncart_input;
    out_shift = 7 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 1, 2)
    inp_shift = 8 * ncart_input;
    out_shift = 8 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 0, 3)
    inp_shift = 9 * ncart_input;
    out_shift = 9 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 4, 0)
    inp_shift = 10 * ncart_input;
    out_shift = 10 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 3, 1)
    inp_shift = 11 * ncart_input;
    out_shift = 11 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 2, 2)
    inp_shift = 12 * ncart_input;
    out_shift = 12 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 1, 3)
    inp_shift = 13 * ncart_input;
    out_shift = 13 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 0, 4)
    inp_shift = 14 * ncart_input;
    out_shift = 14 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 5, 0)
    inp_shift = 15 * ncart_input;
    out_shift = 15 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 4, 1)
    inp_shift = 16 * ncart_input;
    out_shift = 16 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 3, 2)
    inp_shift = 17 * ncart_input;
    out_shift = 17 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 2, 3)
    inp_shift = 18 * ncart_input;
    out_shift = 18 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 1, 4)
    inp_shift = 19 * ncart_input;
    out_shift = 19 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 0, 5)
    inp_shift = 20 * ncart_input;
    out_shift = 20 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 6, 0)
    inp_shift = 21 * ncart_input;
    out_shift = 21 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 5, 1)
    inp_shift = 22 * ncart_input;
    out_shift = 22 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 4, 2)
    inp_shift = 23 * ncart_input;
    out_shift = 23 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 3, 3)
    inp_shift = 24 * ncart_input;
    out_shift = 24 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 2, 4)
    inp_shift = 25 * ncart_input;
    out_shift = 25 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 1, 5)
    inp_shift = 26 * ncart_input;
    out_shift = 26 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 0, 6)
    inp_shift = 27 * ncart_input;
    out_shift = 27 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 7, 0)
    inp_shift = 28 * ncart_input;
    out_shift = 28 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 6, 1)
    inp_shift = 29 * ncart_input;
    out_shift = 29 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 5, 2)
    inp_shift = 30 * ncart_input;
    out_shift = 30 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 4, 3)
    inp_shift = 31 * ncart_input;
    out_shift = 31 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 3, 4)
    inp_shift = 32 * ncart_input;
    out_shift = 32 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 2, 5)
    inp_shift = 33 * ncart_input;
    out_shift = 33 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 1, 6)
    inp_shift = 34 * ncart_input;
    out_shift = 34 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 0, 7)
    inp_shift = 35 * ncart_input;
    out_shift = 35 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_cca_cart_sum_L7(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                        const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                        double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (7, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (6, 1, 0)
    in_shift = 1 * ncart_input;
    coef = vector[1];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (6, 0, 1)
    in_shift = 2 * ncart_input;
    coef = vector[2];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (5, 2, 0)
    in_shift = 3 * ncart_input;
    coef = vector[3];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (5, 1, 1)
    in_shift = 4 * ncart_input;
    coef = vector[4];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (5, 0, 2)
    in_shift = 5 * ncart_input;
    coef = vector[5];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 3, 0)
    in_shift = 6 * ncart_input;
    coef = vector[6];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 2, 1)
    in_shift = 7 * ncart_input;
    coef = vector[7];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 1, 2)
    in_shift = 8 * ncart_input;
    coef = vector[8];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 0, 3)
    in_shift = 9 * ncart_input;
    coef = vector[9];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 4, 0)
    in_shift = 10 * ncart_input;
    coef = vector[10];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 3, 1)
    in_shift = 11 * ncart_input;
    coef = vector[11];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 2, 2)
    in_shift = 12 * ncart_input;
    coef = vector[12];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 1, 3)
    in_shift = 13 * ncart_input;
    coef = vector[13];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 0, 4)
    in_shift = 14 * ncart_input;
    coef = vector[14];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 5, 0)
    in_shift = 15 * ncart_input;
    coef = vector[15];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 4, 1)
    in_shift = 16 * ncart_input;
    coef = vector[16];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 3, 2)
    in_shift = 17 * ncart_input;
    coef = vector[17];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 2, 3)
    in_shift = 18 * ncart_input;
    coef = vector[18];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 1, 4)
    in_shift = 19 * ncart_input;
    coef = vector[19];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 0, 5)
    in_shift = 20 * ncart_input;
    coef = vector[20];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 6, 0)
    in_shift = 21 * ncart_input;
    coef = vector[21];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 5, 1)
    in_shift = 22 * ncart_input;
    coef = vector[22];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 4, 2)
    in_shift = 23 * ncart_input;
    coef = vector[23];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 3, 3)
    in_shift = 24 * ncart_input;
    coef = vector[24];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 2, 4)
    in_shift = 25 * ncart_input;
    coef = vector[25];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 1, 5)
    in_shift = 26 * ncart_input;
    coef = vector[26];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 0, 6)
    in_shift = 27 * ncart_input;
    coef = vector[27];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 7, 0)
    in_shift = 28 * ncart_input;
    coef = vector[28];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 6, 1)
    in_shift = 29 * ncart_input;
    coef = vector[29];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 5, 2)
    in_shift = 30 * ncart_input;
    coef = vector[30];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 4, 3)
    in_shift = 31 * ncart_input;
    coef = vector[31];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 3, 4)
    in_shift = 32 * ncart_input;
    coef = vector[32];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 2, 5)
    in_shift = 33 * ncart_input;
    coef = vector[33];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 1, 6)
    in_shift = 34 * ncart_input;
    coef = vector[34];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 0, 7)
    in_shift = 35 * ncart_input;
    coef = vector[35];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_cca_cart_copy_L8(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                         const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                         const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (8, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (7, 1, 0)
    inp_shift = 1 * ncart_input;
    out_shift = 1 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (7, 0, 1)
    inp_shift = 2 * ncart_input;
    out_shift = 2 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (6, 2, 0)
    inp_shift = 3 * ncart_input;
    out_shift = 3 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (6, 1, 1)
    inp_shift = 4 * ncart_input;
    out_shift = 4 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (6, 0, 2)
    inp_shift = 5 * ncart_input;
    out_shift = 5 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (5, 3, 0)
    inp_shift = 6 * ncart_input;
    out_shift = 6 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (5, 2, 1)
    inp_shift = 7 * ncart_input;
    out_shift = 7 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (5, 1, 2)
    inp_shift = 8 * ncart_input;
    out_shift = 8 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (5, 0, 3)
    inp_shift = 9 * ncart_input;
    out_shift = 9 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 4, 0)
    inp_shift = 10 * ncart_input;
    out_shift = 10 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 3, 1)
    inp_shift = 11 * ncart_input;
    out_shift = 11 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 2, 2)
    inp_shift = 12 * ncart_input;
    out_shift = 12 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 1, 3)
    inp_shift = 13 * ncart_input;
    out_shift = 13 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (4, 0, 4)
    inp_shift = 14 * ncart_input;
    out_shift = 14 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 5, 0)
    inp_shift = 15 * ncart_input;
    out_shift = 15 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 4, 1)
    inp_shift = 16 * ncart_input;
    out_shift = 16 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 3, 2)
    inp_shift = 17 * ncart_input;
    out_shift = 17 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 2, 3)
    inp_shift = 18 * ncart_input;
    out_shift = 18 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 1, 4)
    inp_shift = 19 * ncart_input;
    out_shift = 19 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 0, 5)
    inp_shift = 20 * ncart_input;
    out_shift = 20 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 6, 0)
    inp_shift = 21 * ncart_input;
    out_shift = 21 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 5, 1)
    inp_shift = 22 * ncart_input;
    out_shift = 22 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 4, 2)
    inp_shift = 23 * ncart_input;
    out_shift = 23 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 3, 3)
    inp_shift = 24 * ncart_input;
    out_shift = 24 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 2, 4)
    inp_shift = 25 * ncart_input;
    out_shift = 25 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 1, 5)
    inp_shift = 26 * ncart_input;
    out_shift = 26 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 0, 6)
    inp_shift = 27 * ncart_input;
    out_shift = 27 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 7, 0)
    inp_shift = 28 * ncart_input;
    out_shift = 28 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 6, 1)
    inp_shift = 29 * ncart_input;
    out_shift = 29 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 5, 2)
    inp_shift = 30 * ncart_input;
    out_shift = 30 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 4, 3)
    inp_shift = 31 * ncart_input;
    out_shift = 31 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 3, 4)
    inp_shift = 32 * ncart_input;
    out_shift = 32 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 2, 5)
    inp_shift = 33 * ncart_input;
    out_shift = 33 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 1, 6)
    inp_shift = 34 * ncart_input;
    out_shift = 34 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 0, 7)
    inp_shift = 35 * ncart_input;
    out_shift = 35 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 8, 0)
    inp_shift = 36 * ncart_input;
    out_shift = 36 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 7, 1)
    inp_shift = 37 * ncart_input;
    out_shift = 37 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 6, 2)
    inp_shift = 38 * ncart_input;
    out_shift = 38 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 5, 3)
    inp_shift = 39 * ncart_input;
    out_shift = 39 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 4, 4)
    inp_shift = 40 * ncart_input;
    out_shift = 40 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 3, 5)
    inp_shift = 41 * ncart_input;
    out_shift = 41 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 2, 6)
    inp_shift = 42 * ncart_input;
    out_shift = 42 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 1, 7)
    inp_shift = 43 * ncart_input;
    out_shift = 43 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 0, 8)
    inp_shift = 44 * ncart_input;
    out_shift = 44 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_cca_cart_sum_L8(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                        const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                        double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (8, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (7, 1, 0)
    in_shift = 1 * ncart_input;
    coef = vector[1];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (7, 0, 1)
    in_shift = 2 * ncart_input;
    coef = vector[2];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (6, 2, 0)
    in_shift = 3 * ncart_input;
    coef = vector[3];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (6, 1, 1)
    in_shift = 4 * ncart_input;
    coef = vector[4];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (6, 0, 2)
    in_shift = 5 * ncart_input;
    coef = vector[5];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (5, 3, 0)
    in_shift = 6 * ncart_input;
    coef = vector[6];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (5, 2, 1)
    in_shift = 7 * ncart_input;
    coef = vector[7];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (5, 1, 2)
    in_shift = 8 * ncart_input;
    coef = vector[8];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (5, 0, 3)
    in_shift = 9 * ncart_input;
    coef = vector[9];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 4, 0)
    in_shift = 10 * ncart_input;
    coef = vector[10];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 3, 1)
    in_shift = 11 * ncart_input;
    coef = vector[11];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 2, 2)
    in_shift = 12 * ncart_input;
    coef = vector[12];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 1, 3)
    in_shift = 13 * ncart_input;
    coef = vector[13];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (4, 0, 4)
    in_shift = 14 * ncart_input;
    coef = vector[14];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 5, 0)
    in_shift = 15 * ncart_input;
    coef = vector[15];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 4, 1)
    in_shift = 16 * ncart_input;
    coef = vector[16];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 3, 2)
    in_shift = 17 * ncart_input;
    coef = vector[17];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 2, 3)
    in_shift = 18 * ncart_input;
    coef = vector[18];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 1, 4)
    in_shift = 19 * ncart_input;
    coef = vector[19];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 0, 5)
    in_shift = 20 * ncart_input;
    coef = vector[20];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 6, 0)
    in_shift = 21 * ncart_input;
    coef = vector[21];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 5, 1)
    in_shift = 22 * ncart_input;
    coef = vector[22];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 4, 2)
    in_shift = 23 * ncart_input;
    coef = vector[23];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 3, 3)
    in_shift = 24 * ncart_input;
    coef = vector[24];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 2, 4)
    in_shift = 25 * ncart_input;
    coef = vector[25];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 1, 5)
    in_shift = 26 * ncart_input;
    coef = vector[26];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 0, 6)
    in_shift = 27 * ncart_input;
    coef = vector[27];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 7, 0)
    in_shift = 28 * ncart_input;
    coef = vector[28];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 6, 1)
    in_shift = 29 * ncart_input;
    coef = vector[29];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 5, 2)
    in_shift = 30 * ncart_input;
    coef = vector[30];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 4, 3)
    in_shift = 31 * ncart_input;
    coef = vector[31];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 3, 4)
    in_shift = 32 * ncart_input;
    coef = vector[32];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 2, 5)
    in_shift = 33 * ncart_input;
    coef = vector[33];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 1, 6)
    in_shift = 34 * ncart_input;
    coef = vector[34];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 0, 7)
    in_shift = 35 * ncart_input;
    coef = vector[35];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 8, 0)
    in_shift = 36 * ncart_input;
    coef = vector[36];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 7, 1)
    in_shift = 37 * ncart_input;
    coef = vector[37];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 6, 2)
    in_shift = 38 * ncart_input;
    coef = vector[38];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 5, 3)
    in_shift = 39 * ncart_input;
    coef = vector[39];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 4, 4)
    in_shift = 40 * ncart_input;
    coef = vector[40];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 3, 5)
    in_shift = 41 * ncart_input;
    coef = vector[41];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 2, 6)
    in_shift = 42 * ncart_input;
    coef = vector[42];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 1, 7)
    in_shift = 43 * ncart_input;
    coef = vector[43];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 0, 8)
    in_shift = 44 * ncart_input;
    coef = vector[44];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_molden_cart_copy_L0(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                            const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                            const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (0, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_molden_cart_sum_L0(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                           const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                           double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (0, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_molden_cart_copy_L1(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                            const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                            const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (1, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 1, 0)
    inp_shift = 1 * ncart_input;
    out_shift = 1 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 0, 1)
    inp_shift = 2 * ncart_input;
    out_shift = 2 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_molden_cart_sum_L1(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                           const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                           double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (1, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 1, 0)
    in_shift = 1 * ncart_input;
    coef = vector[1];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 0, 1)
    in_shift = 2 * ncart_input;
    coef = vector[2];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_molden_cart_copy_L2(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                            const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                            const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (2, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 1, 0)
    inp_shift = 1 * ncart_input;
    out_shift = 3 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 0, 1)
    inp_shift = 2 * ncart_input;
    out_shift = 4 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 2, 0)
    inp_shift = 3 * ncart_input;
    out_shift = 1 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 1, 1)
    inp_shift = 4 * ncart_input;
    out_shift = 5 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 0, 2)
    inp_shift = 5 * ncart_input;
    out_shift = 2 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_molden_cart_sum_L2(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                           const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                           double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (2, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 1, 0)
    in_shift = 1 * ncart_input;
    coef = vector[3];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 0, 1)
    in_shift = 2 * ncart_input;
    coef = vector[4];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 2, 0)
    in_shift = 3 * ncart_input;
    coef = vector[1];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 1, 1)
    in_shift = 4 * ncart_input;
    coef = vector[5];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 0, 2)
    in_shift = 5 * ncart_input;
    coef = vector[2];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_molden_cart_copy_L3(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                            const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                            const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (3, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 1, 0)
    inp_shift = 1 * ncart_input;
    out_shift = 4 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 0, 1)
    inp_shift = 2 * ncart_input;
    out_shift = 5 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 2, 0)
    inp_shift = 3 * ncart_input;
    out_shift = 3 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 1, 1)
    inp_shift = 4 * ncart_input;
    out_shift = 9 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 0, 2)
    inp_shift = 5 * ncart_input;
    out_shift = 6 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 3, 0)
    inp_shift = 6 * ncart_input;
    out_shift = 1 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 2, 1)
    inp_shift = 7 * ncart_input;
    out_shift = 8 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 1, 2)
    inp_shift = 8 * ncart_input;
    out_shift = 7 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 0, 3)
    inp_shift = 9 * ncart_input;
    out_shift = 2 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_molden_cart_sum_L3(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                           const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                           double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (3, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 1, 0)
    in_shift = 1 * ncart_input;
    coef = vector[4];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 0, 1)
    in_shift = 2 * ncart_input;
    coef = vector[5];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 2, 0)
    in_shift = 3 * ncart_input;
    coef = vector[3];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 1, 1)
    in_shift = 4 * ncart_input;
    coef = vector[9];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 0, 2)
    in_shift = 5 * ncart_input;
    coef = vector[6];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 3, 0)
    in_shift = 6 * ncart_input;
    coef = vector[1];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 2, 1)
    in_shift = 7 * ncart_input;
    coef = vector[8];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 1, 2)
    in_shift = 8 * ncart_input;
    coef = vector[7];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 0, 3)
    in_shift = 9 * ncart_input;
    coef = vector[2];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_molden_cart_copy_L4(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                            const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                            const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long inp_shift;
    unsigned long out_shift;

    // Copy (4, 0, 0)
    inp_shift = 0 * ncart_input;
    out_shift = 0 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 1, 0)
    inp_shift = 1 * ncart_input;
    out_shift = 3 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (3, 0, 1)
    inp_shift = 2 * ncart_input;
    out_shift = 4 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 2, 0)
    inp_shift = 3 * ncart_input;
    out_shift = 9 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 1, 1)
    inp_shift = 4 * ncart_input;
    out_shift = 12 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (2, 0, 2)
    inp_shift = 5 * ncart_input;
    out_shift = 10 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 3, 0)
    inp_shift = 6 * ncart_input;
    out_shift = 5 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 2, 1)
    inp_shift = 7 * ncart_input;
    out_shift = 13 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 1, 2)
    inp_shift = 8 * ncart_input;
    out_shift = 14 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (1, 0, 3)
    inp_shift = 9 * ncart_input;
    out_shift = 7 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 4, 0)
    inp_shift = 10 * ncart_input;
    out_shift = 1 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 3, 1)
    inp_shift = 11 * ncart_input;
    out_shift = 6 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 2, 2)
    inp_shift = 12 * ncart_input;
    out_shift = 11 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 1, 3)
    inp_shift = 13 * ncart_input;
    out_shift = 8 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }

    // Copy (0, 0, 4)
    inp_shift = 14 * ncart_input;
    out_shift = 2 * ncart_out;
    for (unsigned long i = 0; i < size; i++) {
        cart_out[out_shift + i] = cart_input[inp_shift + i];
    }
}
void gg_molden_cart_sum_L4(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                           const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                           double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
    ASSUME_ALIGNED(cart_input, 64);
    unsigned long in_shift;
    unsigned long out_shift;
    double coef;

    // Copy (4, 0, 0)
    in_shift = 0 * ncart_input;
    coef = vector[0];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 1, 0)
    in_shift = 1 * ncart_input;
    coef = vector[3];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (3, 0, 1)
    in_shift = 2 * ncart_input;
    coef = vector[4];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 2, 0)
    in_shift = 3 * ncart_input;
    coef = vector[9];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 1, 1)
    in_shift = 4 * ncart_input;
    coef = vector[12];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (2, 0, 2)
    in_shift = 5 * ncart_input;
    coef = vector[10];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 3, 0)
    in_shift = 6 * ncart_input;
    coef = vector[5];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 2, 1)
    in_shift = 7 * ncart_input;
    coef = vector[13];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 1, 2)
    in_shift = 8 * ncart_input;
    coef = vector[14];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (1, 0, 3)
    in_shift = 9 * ncart_input;
    coef = vector[7];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 4, 0)
    in_shift = 10 * ncart_input;
    coef = vector[1];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 3, 1)
    in_shift = 11 * ncart_input;
    coef = vector[6];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 2, 2)
    in_shift = 12 * ncart_input;
    coef = vector[11];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 1, 3)
    in_shift = 13 * ncart_input;
    coef = vector[8];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }

    // Copy (0, 0, 4)
    in_shift = 14 * ncart_input;
    coef = vector[2];
    for (unsigned long i = 0; i < size; i++) {
        cart_out[i] += coef * cart_input[in_shift + i];
    }
}
void gg_molden_cart_copy_L5(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                            const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                            const unsigned long ncart_out) {}
void gg_molden_cart_sum_L5(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                           const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                           double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {}
void gg_molden_cart_copy_L6(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                            const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                            const unsigned long ncart_out) {}
void gg_molden_cart_sum_L6(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                           const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                           double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {}
void gg_molden_cart_copy_L7(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                            const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                            const unsigned long ncart_out) {}
void gg_molden_cart_sum_L7(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                           const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                           double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {}
void gg_molden_cart_copy_L8(const unsigned long size, const double* PRAGMA_RESTRICT cart_input,
                            const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out,
                            const unsigned long ncart_out) {}
void gg_molden_cart_sum_L8(const unsigned long size, const double* PRAGMA_RESTRICT vector,
                           const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input,
                           double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {}
void gg_naive_transpose(unsigned long n, unsigned long m, const double* PRAGMA_RESTRICT input,
                        double* PRAGMA_RESTRICT output) {
    ASSUME_ALIGNED(input, 64);
    for (unsigned long i = 0; i < n; i++) {
        for (unsigned long j = 0; j < m; j++) {
            output[j * n + i] = input[i * m + j];
        }
    }
}
void gg_fast_transpose(unsigned long n, unsigned long m, const double* PRAGMA_RESTRICT input,
                       double* PRAGMA_RESTRICT output) {
// Temps
#ifdef _MSC_VER
    __declspec(align(64)) double tmp[64];
#else
    double tmp[64] __attribute__((aligned(64)));
#endif
    ASSUME_ALIGNED(input, 64);
    // Sizing
    unsigned long nblocks = n / 8;
    nblocks += (n % 8) ? 1 : 0;
    unsigned long mblocks = m / 8;
    mblocks += (m % 8) ? 1 : 0;
    // Outer blocks
    for (unsigned long nb = 0; nb < nblocks; nb++) {
        const unsigned long nstart = nb * 8;
        unsigned long nremain = ((nstart + 8) > n) ? (n - nstart) : 8;
        for (unsigned long mb = 0; mb < mblocks; mb++) {
            const unsigned long mstart = mb * 8;
            unsigned long mremain = ((mstart + 8) > m) ? (m - mstart) : 8;
            // Copy data to inner block
            for (unsigned long l = 0; l < nremain; l++) {
                const unsigned long start = (nstart + l) * m + mstart;
                for (unsigned long k = 0; k < mremain; k++) {
                    tmp[k * 8 + l] = input[start + k];
                }
            }
            // Copy data to inner block
            for (unsigned long k = 0; k < mremain; k++) {
                const unsigned long start = (mstart + k) * n + nstart;
                for (unsigned long l = 0; l < nremain; l++) {
                    output[start + l] = tmp[k * 8 + l];
                }
            }
        }
    }
}
void block_copy(unsigned long n, unsigned long m, const double* PRAGMA_RESTRICT input, unsigned long is,
                double* PRAGMA_RESTRICT output, unsigned long os, const int trans) {
    ASSUME_ALIGNED(input, 64);
    for (unsigned long i = 0; i < n; i++) {
        const unsigned long out_shift = i * os;
        const unsigned long inp_shift = i * is;

        for (unsigned long j = 0; j < m; j++) {
            output[out_shift + j] = input[inp_shift + j];
        }
    }
}
void block_matrix_vector(unsigned long n, unsigned long m, const double* vector, const double* PRAGMA_RESTRICT input,
                         unsigned long is, double* PRAGMA_RESTRICT output) {
    ASSUME_ALIGNED(input, 64);
    for (unsigned long i = 0; i < n; i++) {
        const unsigned long inp_shift = i * is;
        const double coef = vector[i];

        for (unsigned long j = 0; j < m; j++) {
            output[j] += coef * input[inp_shift + j];
        }
    }
}