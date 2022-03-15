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

void gg_cca_cart_to_spherical_L0(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_00 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  = cart[i];

    }

}
void gg_cca_cart_to_spherical_sum_L0(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_00 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = cart[i];
        output[i] += tmp * vector[0];

    }

}
void gg_cca_cart_to_spherical_L1(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_10 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  = cart[ncart + i];

    }

    // R_11c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i]  = cart[2 * ncart + i];

    }
    // R_11s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i]  = cart[i];

    }

}
void gg_cca_cart_to_spherical_sum_L1(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_10 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = cart[ncart + i];
        output[i] += tmp * vector[0];

    }

    // R_11c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = cart[2 * ncart + i];
        output[i] += tmp * vector[1];

    }
    // R_11s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = cart[i];
        output[i] += tmp * vector[2];

    }

}
void gg_cca_cart_to_spherical_L2(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_20 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  =  1.7320508075688772 * cart[ncart + i];

    }

    // R_21c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i]  =  1.7320508075688772 * cart[4 * ncart + i];

    }
    // R_21s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i]  = -0.5000000000000000 * cart[i];
        spherical[2 * nspherical + i] += -0.5000000000000000 * cart[3 * ncart + i];
        spherical[2 * nspherical + i] += cart[5 * ncart + i];

    }

    // R_22c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i]  =  1.7320508075688772 * cart[2 * ncart + i];

    }
    // R_22s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i]  =  0.8660254037844386 * cart[i];
        spherical[4 * nspherical + i] += -0.8660254037844386 * cart[3 * ncart + i];

    }

}
void gg_cca_cart_to_spherical_sum_L2(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_20 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  1.7320508075688772 * cart[ncart + i];
        output[i] += tmp * vector[0];

    }

    // R_21c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  1.7320508075688772 * cart[4 * ncart + i];
        output[i] += tmp * vector[1];

    }
    // R_21s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.5000000000000000 * cart[i];
        tmp += -0.5000000000000000 * cart[3 * ncart + i];
        tmp += cart[5 * ncart + i];
        output[i] += tmp * vector[2];

    }

    // R_22c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  1.7320508075688772 * cart[2 * ncart + i];
        output[i] += tmp * vector[3];

    }
    // R_22s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.8660254037844386 * cart[i];
        tmp += -0.8660254037844386 * cart[3 * ncart + i];
        output[i] += tmp * vector[4];

    }

}
void gg_cca_cart_to_spherical_L3(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_30 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  =  2.3717082451262845 * cart[ncart + i];
        spherical[i] += -0.7905694150420949 * cart[6 * ncart + i];

    }

    // R_31c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i]  =  3.8729833462074170 * cart[4 * ncart + i];

    }
    // R_31s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i]  = -0.6123724356957945 * cart[ncart + i];
        spherical[2 * nspherical + i] += -0.6123724356957945 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] +=  2.4494897427831779 * cart[8 * ncart + i];

    }

    // R_32c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i]  = -1.5000000000000000 * cart[2 * ncart + i];
        spherical[3 * nspherical + i] += -1.5000000000000000 * cart[7 * ncart + i];
        spherical[3 * nspherical + i] += cart[9 * ncart + i];

    }
    // R_32s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i]  = -0.6123724356957945 * cart[i];
        spherical[4 * nspherical + i] += -0.6123724356957945 * cart[3 * ncart + i];
        spherical[4 * nspherical + i] +=  2.4494897427831779 * cart[5 * ncart + i];

    }

    // R_33c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i]  =  1.9364916731037085 * cart[2 * ncart + i];
        spherical[5 * nspherical + i] += -1.9364916731037085 * cart[7 * ncart + i];

    }
    // R_33s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i]  =  0.7905694150420949 * cart[i];
        spherical[6 * nspherical + i] += -2.3717082451262845 * cart[3 * ncart + i];

    }

}
void gg_cca_cart_to_spherical_sum_L3(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_30 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.3717082451262845 * cart[ncart + i];
        tmp += -0.7905694150420949 * cart[6 * ncart + i];
        output[i] += tmp * vector[0];

    }

    // R_31c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  3.8729833462074170 * cart[4 * ncart + i];
        output[i] += tmp * vector[1];

    }
    // R_31s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.6123724356957945 * cart[ncart + i];
        tmp += -0.6123724356957945 * cart[6 * ncart + i];
        tmp +=  2.4494897427831779 * cart[8 * ncart + i];
        output[i] += tmp * vector[2];

    }

    // R_32c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -1.5000000000000000 * cart[2 * ncart + i];
        tmp += -1.5000000000000000 * cart[7 * ncart + i];
        tmp += cart[9 * ncart + i];
        output[i] += tmp * vector[3];

    }
    // R_32s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.6123724356957945 * cart[i];
        tmp += -0.6123724356957945 * cart[3 * ncart + i];
        tmp +=  2.4494897427831779 * cart[5 * ncart + i];
        output[i] += tmp * vector[4];

    }

    // R_33c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  1.9364916731037085 * cart[2 * ncart + i];
        tmp += -1.9364916731037085 * cart[7 * ncart + i];
        output[i] += tmp * vector[5];

    }
    // R_33s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.7905694150420949 * cart[i];
        tmp += -2.3717082451262845 * cart[3 * ncart + i];
        output[i] += tmp * vector[6];

    }

}
void gg_cca_cart_to_spherical_L4(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_40 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  =  2.9580398915498081 * cart[ncart + i];
        spherical[i] += -2.9580398915498081 * cart[6 * ncart + i];

    }

    // R_41c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i]  =  6.2749501990055663 * cart[4 * ncart + i];
        spherical[nspherical + i] += -2.0916500663351889 * cart[11 * ncart + i];

    }
    // R_41s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i]  = -1.1180339887498949 * cart[ncart + i];
        spherical[2 * nspherical + i] += -1.1180339887498949 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] +=  6.7082039324993694 * cart[8 * ncart + i];

    }

    // R_42c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i]  = -2.3717082451262845 * cart[4 * ncart + i];
        spherical[3 * nspherical + i] += -2.3717082451262845 * cart[11 * ncart + i];
        spherical[3 * nspherical + i] +=  3.1622776601683795 * cart[13 * ncart + i];

    }
    // R_42s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i]  =  0.3750000000000000 * cart[i];
        spherical[4 * nspherical + i] +=  0.7500000000000000 * cart[3 * ncart + i];
        spherical[4 * nspherical + i] +=  0.3750000000000000 * cart[10 * ncart + i];
        spherical[4 * nspherical + i] += -3.0000000000000000 * cart[5 * ncart + i];
        spherical[4 * nspherical + i] += -3.0000000000000000 * cart[12 * ncart + i];
        spherical[4 * nspherical + i] += cart[14 * ncart + i];

    }

    // R_43c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i]  = -2.3717082451262845 * cart[2 * ncart + i];
        spherical[5 * nspherical + i] += -2.3717082451262845 * cart[7 * ncart + i];
        spherical[5 * nspherical + i] +=  3.1622776601683795 * cart[9 * ncart + i];

    }
    // R_43s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i]  = -0.5590169943749475 * cart[i];
        spherical[6 * nspherical + i] +=  0.5590169943749475 * cart[10 * ncart + i];
        spherical[6 * nspherical + i] +=  3.3541019662496847 * cart[5 * ncart + i];
        spherical[6 * nspherical + i] += -3.3541019662496847 * cart[12 * ncart + i];

    }

    // R_44c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i]  =  2.0916500663351889 * cart[2 * ncart + i];
        spherical[7 * nspherical + i] += -6.2749501990055663 * cart[7 * ncart + i];

    }
    // R_44s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i]  =  0.7395099728874520 * cart[i];
        spherical[8 * nspherical + i] += -4.4370598373247123 * cart[3 * ncart + i];
        spherical[8 * nspherical + i] +=  0.7395099728874520 * cart[10 * ncart + i];

    }

}
void gg_cca_cart_to_spherical_sum_L4(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_40 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.9580398915498081 * cart[ncart + i];
        tmp += -2.9580398915498081 * cart[6 * ncart + i];
        output[i] += tmp * vector[0];

    }

    // R_41c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  6.2749501990055663 * cart[4 * ncart + i];
        tmp += -2.0916500663351889 * cart[11 * ncart + i];
        output[i] += tmp * vector[1];

    }
    // R_41s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -1.1180339887498949 * cart[ncart + i];
        tmp += -1.1180339887498949 * cart[6 * ncart + i];
        tmp +=  6.7082039324993694 * cart[8 * ncart + i];
        output[i] += tmp * vector[2];

    }

    // R_42c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -2.3717082451262845 * cart[4 * ncart + i];
        tmp += -2.3717082451262845 * cart[11 * ncart + i];
        tmp +=  3.1622776601683795 * cart[13 * ncart + i];
        output[i] += tmp * vector[3];

    }
    // R_42s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.3750000000000000 * cart[i];
        tmp +=  0.7500000000000000 * cart[3 * ncart + i];
        tmp +=  0.3750000000000000 * cart[10 * ncart + i];
        tmp += -3.0000000000000000 * cart[5 * ncart + i];
        tmp += -3.0000000000000000 * cart[12 * ncart + i];
        tmp += cart[14 * ncart + i];
        output[i] += tmp * vector[4];

    }

    // R_43c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -2.3717082451262845 * cart[2 * ncart + i];
        tmp += -2.3717082451262845 * cart[7 * ncart + i];
        tmp +=  3.1622776601683795 * cart[9 * ncart + i];
        output[i] += tmp * vector[5];

    }
    // R_43s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.5590169943749475 * cart[i];
        tmp +=  0.5590169943749475 * cart[10 * ncart + i];
        tmp +=  3.3541019662496847 * cart[5 * ncart + i];
        tmp += -3.3541019662496847 * cart[12 * ncart + i];
        output[i] += tmp * vector[6];

    }

    // R_44c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.0916500663351889 * cart[2 * ncart + i];
        tmp += -6.2749501990055663 * cart[7 * ncart + i];
        output[i] += tmp * vector[7];

    }
    // R_44s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.7395099728874520 * cart[i];
        tmp += -4.4370598373247123 * cart[3 * ncart + i];
        tmp +=  0.7395099728874520 * cart[10 * ncart + i];
        output[i] += tmp * vector[8];

    }

}
void gg_cca_cart_to_spherical_L5(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_50 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  =  3.5078038001005702 * cart[ncart + i];
        spherical[i] += -7.0156076002011405 * cart[6 * ncart + i];
        spherical[i] +=  0.7015607600201140 * cart[15 * ncart + i];

    }

    // R_51c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i]  =  8.8741196746494246 * cart[4 * ncart + i];
        spherical[nspherical + i] += -8.8741196746494246 * cart[11 * ncart + i];

    }
    // R_51s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i]  = -1.5687375497513916 * cart[ncart + i];
        spherical[2 * nspherical + i] += -1.0458250331675945 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] +=  0.5229125165837972 * cart[15 * ncart + i];
        spherical[2 * nspherical + i] +=  12.5499003980111326 * cart[8 * ncart + i];
        spherical[2 * nspherical + i] += -4.1833001326703778 * cart[17 * ncart + i];

    }

    // R_52c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i]  = -5.1234753829797990 * cart[4 * ncart + i];
        spherical[3 * nspherical + i] += -5.1234753829797990 * cart[11 * ncart + i];
        spherical[3 * nspherical + i] +=  10.2469507659595980 * cart[13 * ncart + i];

    }
    // R_52s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i]  =  0.4841229182759271 * cart[ncart + i];
        spherical[4 * nspherical + i] +=  0.9682458365518543 * cart[6 * ncart + i];
        spherical[4 * nspherical + i] +=  0.4841229182759271 * cart[15 * ncart + i];
        spherical[4 * nspherical + i] += -5.8094750193111251 * cart[8 * ncart + i];
        spherical[4 * nspherical + i] += -5.8094750193111251 * cart[17 * ncart + i];
        spherical[4 * nspherical + i] +=  3.8729833462074170 * cart[19 * ncart + i];

    }

    // R_53c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i]  =  1.8750000000000000 * cart[2 * ncart + i];
        spherical[5 * nspherical + i] +=  3.7500000000000000 * cart[7 * ncart + i];
        spherical[5 * nspherical + i] +=  1.8750000000000000 * cart[16 * ncart + i];
        spherical[5 * nspherical + i] += -5.0000000000000000 * cart[9 * ncart + i];
        spherical[5 * nspherical + i] += -5.0000000000000000 * cart[18 * ncart + i];
        spherical[5 * nspherical + i] += cart[20 * ncart + i];

    }
    // R_53s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i]  =  0.4841229182759271 * cart[i];
        spherical[6 * nspherical + i] +=  0.9682458365518543 * cart[3 * ncart + i];
        spherical[6 * nspherical + i] +=  0.4841229182759271 * cart[10 * ncart + i];
        spherical[6 * nspherical + i] += -5.8094750193111251 * cart[5 * ncart + i];
        spherical[6 * nspherical + i] += -5.8094750193111251 * cart[12 * ncart + i];
        spherical[6 * nspherical + i] +=  3.8729833462074170 * cart[14 * ncart + i];

    }

    // R_54c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i]  = -2.5617376914898995 * cart[2 * ncart + i];
        spherical[7 * nspherical + i] +=  2.5617376914898995 * cart[16 * ncart + i];
        spherical[7 * nspherical + i] +=  5.1234753829797990 * cart[9 * ncart + i];
        spherical[7 * nspherical + i] += -5.1234753829797990 * cart[18 * ncart + i];

    }
    // R_54s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i]  = -0.5229125165837972 * cart[i];
        spherical[8 * nspherical + i] +=  1.0458250331675945 * cart[3 * ncart + i];
        spherical[8 * nspherical + i] +=  1.5687375497513916 * cart[10 * ncart + i];
        spherical[8 * nspherical + i] +=  4.1833001326703778 * cart[5 * ncart + i];
        spherical[8 * nspherical + i] += -12.5499003980111326 * cart[12 * ncart + i];

    }

    // R_55c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[9 * nspherical + i]  =  2.2185299186623562 * cart[2 * ncart + i];
        spherical[9 * nspherical + i] += -13.3111795119741370 * cart[7 * ncart + i];
        spherical[9 * nspherical + i] +=  2.2185299186623562 * cart[16 * ncart + i];

    }
    // R_55s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[10 * nspherical + i]  =  0.7015607600201140 * cart[i];
        spherical[10 * nspherical + i] += -7.0156076002011405 * cart[3 * ncart + i];
        spherical[10 * nspherical + i] +=  3.5078038001005702 * cart[10 * ncart + i];

    }

}
void gg_cca_cart_to_spherical_sum_L5(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_50 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  3.5078038001005702 * cart[ncart + i];
        tmp += -7.0156076002011405 * cart[6 * ncart + i];
        tmp +=  0.7015607600201140 * cart[15 * ncart + i];
        output[i] += tmp * vector[0];

    }

    // R_51c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  8.8741196746494246 * cart[4 * ncart + i];
        tmp += -8.8741196746494246 * cart[11 * ncart + i];
        output[i] += tmp * vector[1];

    }
    // R_51s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -1.5687375497513916 * cart[ncart + i];
        tmp += -1.0458250331675945 * cart[6 * ncart + i];
        tmp +=  0.5229125165837972 * cart[15 * ncart + i];
        tmp +=  12.5499003980111326 * cart[8 * ncart + i];
        tmp += -4.1833001326703778 * cart[17 * ncart + i];
        output[i] += tmp * vector[2];

    }

    // R_52c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -5.1234753829797990 * cart[4 * ncart + i];
        tmp += -5.1234753829797990 * cart[11 * ncart + i];
        tmp +=  10.2469507659595980 * cart[13 * ncart + i];
        output[i] += tmp * vector[3];

    }
    // R_52s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.4841229182759271 * cart[ncart + i];
        tmp +=  0.9682458365518543 * cart[6 * ncart + i];
        tmp +=  0.4841229182759271 * cart[15 * ncart + i];
        tmp += -5.8094750193111251 * cart[8 * ncart + i];
        tmp += -5.8094750193111251 * cart[17 * ncart + i];
        tmp +=  3.8729833462074170 * cart[19 * ncart + i];
        output[i] += tmp * vector[4];

    }

    // R_53c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  1.8750000000000000 * cart[2 * ncart + i];
        tmp +=  3.7500000000000000 * cart[7 * ncart + i];
        tmp +=  1.8750000000000000 * cart[16 * ncart + i];
        tmp += -5.0000000000000000 * cart[9 * ncart + i];
        tmp += -5.0000000000000000 * cart[18 * ncart + i];
        tmp += cart[20 * ncart + i];
        output[i] += tmp * vector[5];

    }
    // R_53s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.4841229182759271 * cart[i];
        tmp +=  0.9682458365518543 * cart[3 * ncart + i];
        tmp +=  0.4841229182759271 * cart[10 * ncart + i];
        tmp += -5.8094750193111251 * cart[5 * ncart + i];
        tmp += -5.8094750193111251 * cart[12 * ncart + i];
        tmp +=  3.8729833462074170 * cart[14 * ncart + i];
        output[i] += tmp * vector[6];

    }

    // R_54c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -2.5617376914898995 * cart[2 * ncart + i];
        tmp +=  2.5617376914898995 * cart[16 * ncart + i];
        tmp +=  5.1234753829797990 * cart[9 * ncart + i];
        tmp += -5.1234753829797990 * cart[18 * ncart + i];
        output[i] += tmp * vector[7];

    }
    // R_54s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.5229125165837972 * cart[i];
        tmp +=  1.0458250331675945 * cart[3 * ncart + i];
        tmp +=  1.5687375497513916 * cart[10 * ncart + i];
        tmp +=  4.1833001326703778 * cart[5 * ncart + i];
        tmp += -12.5499003980111326 * cart[12 * ncart + i];
        output[i] += tmp * vector[8];

    }

    // R_55c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.2185299186623562 * cart[2 * ncart + i];
        tmp += -13.3111795119741370 * cart[7 * ncart + i];
        tmp +=  2.2185299186623562 * cart[16 * ncart + i];
        output[i] += tmp * vector[9];

    }
    // R_55s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.7015607600201140 * cart[i];
        tmp += -7.0156076002011405 * cart[3 * ncart + i];
        tmp +=  3.5078038001005702 * cart[10 * ncart + i];
        output[i] += tmp * vector[10];

    }

}
void gg_cca_cart_to_spherical_L6(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_60 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  =  4.0301597362883772 * cart[ncart + i];
        spherical[i] += -13.4338657876279228 * cart[6 * ncart + i];
        spherical[i] +=  4.0301597362883772 * cart[15 * ncart + i];

    }

    // R_61c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i]  =  11.6340690431164280 * cart[4 * ncart + i];
        spherical[nspherical + i] += -23.2681380862328560 * cart[11 * ncart + i];
        spherical[nspherical + i] +=  2.3268138086232857 * cart[22 * ncart + i];

    }
    // R_61s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i]  = -1.9843134832984430 * cart[ncart + i];
        spherical[2 * nspherical + i] +=  1.9843134832984430 * cart[15 * ncart + i];
        spherical[2 * nspherical + i] +=  19.8431348329844290 * cart[8 * ncart + i];
        spherical[2 * nspherical + i] += -19.8431348329844290 * cart[17 * ncart + i];

    }

    // R_62c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i]  = -8.1513994197315593 * cart[4 * ncart + i];
        spherical[3 * nspherical + i] += -5.4342662798210393 * cart[11 * ncart + i];
        spherical[3 * nspherical + i] +=  2.7171331399105196 * cart[22 * ncart + i];
        spherical[3 * nspherical + i] +=  21.7370651192841571 * cart[13 * ncart + i];
        spherical[3 * nspherical + i] += -7.2456883730947190 * cart[24 * ncart + i];

    }
    // R_62s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i]  =  0.9057110466368399 * cart[ncart + i];
        spherical[4 * nspherical + i] +=  1.8114220932736798 * cart[6 * ncart + i];
        spherical[4 * nspherical + i] +=  0.9057110466368399 * cart[15 * ncart + i];
        spherical[4 * nspherical + i] += -14.4913767461894381 * cart[8 * ncart + i];
        spherical[4 * nspherical + i] += -14.4913767461894381 * cart[17 * ncart + i];
        spherical[4 * nspherical + i] +=  14.4913767461894381 * cart[19 * ncart + i];

    }

    // R_63c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i]  =  2.8641098093473998 * cart[4 * ncart + i];
        spherical[5 * nspherical + i] +=  5.7282196186947996 * cart[11 * ncart + i];
        spherical[5 * nspherical + i] +=  2.8641098093473998 * cart[22 * ncart + i];
        spherical[5 * nspherical + i] += -11.4564392373895991 * cart[13 * ncart + i];
        spherical[5 * nspherical + i] += -11.4564392373895991 * cart[24 * ncart + i];
        spherical[5 * nspherical + i] +=  4.5825756949558398 * cart[26 * ncart + i];

    }
    // R_63s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i]  = -0.3125000000000000 * cart[i];
        spherical[6 * nspherical + i] += -0.9375000000000000 * cart[3 * ncart + i];
        spherical[6 * nspherical + i] += -0.9375000000000000 * cart[10 * ncart + i];
        spherical[6 * nspherical + i] += -0.3125000000000000 * cart[21 * ncart + i];
        spherical[6 * nspherical + i] +=  5.6250000000000000 * cart[5 * ncart + i];
        spherical[6 * nspherical + i] +=  11.2500000000000000 * cart[12 * ncart + i];
        spherical[6 * nspherical + i] +=  5.6250000000000000 * cart[23 * ncart + i];
        spherical[6 * nspherical + i] += -7.5000000000000000 * cart[14 * ncart + i];
        spherical[6 * nspherical + i] += -7.5000000000000000 * cart[25 * ncart + i];
        spherical[6 * nspherical + i] += cart[27 * ncart + i];

    }

    // R_64c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i]  =  2.8641098093473998 * cart[2 * ncart + i];
        spherical[7 * nspherical + i] +=  5.7282196186947996 * cart[7 * ncart + i];
        spherical[7 * nspherical + i] +=  2.8641098093473998 * cart[16 * ncart + i];
        spherical[7 * nspherical + i] += -11.4564392373895991 * cart[9 * ncart + i];
        spherical[7 * nspherical + i] += -11.4564392373895991 * cart[18 * ncart + i];
        spherical[7 * nspherical + i] +=  4.5825756949558398 * cart[20 * ncart + i];

    }
    // R_64s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i]  =  0.4528555233184199 * cart[i];
        spherical[8 * nspherical + i] +=  0.4528555233184199 * cart[3 * ncart + i];
        spherical[8 * nspherical + i] += -0.4528555233184199 * cart[10 * ncart + i];
        spherical[8 * nspherical + i] += -0.4528555233184199 * cart[21 * ncart + i];
        spherical[8 * nspherical + i] += -7.2456883730947190 * cart[5 * ncart + i];
        spherical[8 * nspherical + i] +=  7.2456883730947190 * cart[23 * ncart + i];
        spherical[8 * nspherical + i] +=  7.2456883730947190 * cart[14 * ncart + i];
        spherical[8 * nspherical + i] += -7.2456883730947190 * cart[25 * ncart + i];

    }

    // R_65c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[9 * nspherical + i]  = -2.7171331399105196 * cart[2 * ncart + i];
        spherical[9 * nspherical + i] +=  5.4342662798210393 * cart[7 * ncart + i];
        spherical[9 * nspherical + i] +=  8.1513994197315593 * cart[16 * ncart + i];
        spherical[9 * nspherical + i] +=  7.2456883730947190 * cart[9 * ncart + i];
        spherical[9 * nspherical + i] += -21.7370651192841571 * cart[18 * ncart + i];

    }
    // R_65s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[10 * nspherical + i]  = -0.4960783708246108 * cart[i];
        spherical[10 * nspherical + i] +=  2.4803918541230536 * cart[3 * ncart + i];
        spherical[10 * nspherical + i] +=  2.4803918541230536 * cart[10 * ncart + i];
        spherical[10 * nspherical + i] += -0.4960783708246108 * cart[21 * ncart + i];
        spherical[10 * nspherical + i] +=  4.9607837082461073 * cart[5 * ncart + i];
        spherical[10 * nspherical + i] += -29.7647022494766453 * cart[12 * ncart + i];
        spherical[10 * nspherical + i] +=  4.9607837082461073 * cart[23 * ncart + i];

    }

    // R_66c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[11 * nspherical + i]  =  2.3268138086232857 * cart[2 * ncart + i];
        spherical[11 * nspherical + i] += -23.2681380862328560 * cart[7 * ncart + i];
        spherical[11 * nspherical + i] +=  11.6340690431164280 * cart[16 * ncart + i];

    }
    // R_66s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[12 * nspherical + i]  =  0.6716932893813962 * cart[i];
        spherical[12 * nspherical + i] += -10.0753993407209421 * cart[3 * ncart + i];
        spherical[12 * nspherical + i] +=  10.0753993407209421 * cart[10 * ncart + i];
        spherical[12 * nspherical + i] += -0.6716932893813962 * cart[21 * ncart + i];

    }

}
void gg_cca_cart_to_spherical_sum_L6(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_60 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  4.0301597362883772 * cart[ncart + i];
        tmp += -13.4338657876279228 * cart[6 * ncart + i];
        tmp +=  4.0301597362883772 * cart[15 * ncart + i];
        output[i] += tmp * vector[0];

    }

    // R_61c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  11.6340690431164280 * cart[4 * ncart + i];
        tmp += -23.2681380862328560 * cart[11 * ncart + i];
        tmp +=  2.3268138086232857 * cart[22 * ncart + i];
        output[i] += tmp * vector[1];

    }
    // R_61s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -1.9843134832984430 * cart[ncart + i];
        tmp +=  1.9843134832984430 * cart[15 * ncart + i];
        tmp +=  19.8431348329844290 * cart[8 * ncart + i];
        tmp += -19.8431348329844290 * cart[17 * ncart + i];
        output[i] += tmp * vector[2];

    }

    // R_62c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -8.1513994197315593 * cart[4 * ncart + i];
        tmp += -5.4342662798210393 * cart[11 * ncart + i];
        tmp +=  2.7171331399105196 * cart[22 * ncart + i];
        tmp +=  21.7370651192841571 * cart[13 * ncart + i];
        tmp += -7.2456883730947190 * cart[24 * ncart + i];
        output[i] += tmp * vector[3];

    }
    // R_62s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.9057110466368399 * cart[ncart + i];
        tmp +=  1.8114220932736798 * cart[6 * ncart + i];
        tmp +=  0.9057110466368399 * cart[15 * ncart + i];
        tmp += -14.4913767461894381 * cart[8 * ncart + i];
        tmp += -14.4913767461894381 * cart[17 * ncart + i];
        tmp +=  14.4913767461894381 * cart[19 * ncart + i];
        output[i] += tmp * vector[4];

    }

    // R_63c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.8641098093473998 * cart[4 * ncart + i];
        tmp +=  5.7282196186947996 * cart[11 * ncart + i];
        tmp +=  2.8641098093473998 * cart[22 * ncart + i];
        tmp += -11.4564392373895991 * cart[13 * ncart + i];
        tmp += -11.4564392373895991 * cart[24 * ncart + i];
        tmp +=  4.5825756949558398 * cart[26 * ncart + i];
        output[i] += tmp * vector[5];

    }
    // R_63s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.3125000000000000 * cart[i];
        tmp += -0.9375000000000000 * cart[3 * ncart + i];
        tmp += -0.9375000000000000 * cart[10 * ncart + i];
        tmp += -0.3125000000000000 * cart[21 * ncart + i];
        tmp +=  5.6250000000000000 * cart[5 * ncart + i];
        tmp +=  11.2500000000000000 * cart[12 * ncart + i];
        tmp +=  5.6250000000000000 * cart[23 * ncart + i];
        tmp += -7.5000000000000000 * cart[14 * ncart + i];
        tmp += -7.5000000000000000 * cart[25 * ncart + i];
        tmp += cart[27 * ncart + i];
        output[i] += tmp * vector[6];

    }

    // R_64c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.8641098093473998 * cart[2 * ncart + i];
        tmp +=  5.7282196186947996 * cart[7 * ncart + i];
        tmp +=  2.8641098093473998 * cart[16 * ncart + i];
        tmp += -11.4564392373895991 * cart[9 * ncart + i];
        tmp += -11.4564392373895991 * cart[18 * ncart + i];
        tmp +=  4.5825756949558398 * cart[20 * ncart + i];
        output[i] += tmp * vector[7];

    }
    // R_64s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.4528555233184199 * cart[i];
        tmp +=  0.4528555233184199 * cart[3 * ncart + i];
        tmp += -0.4528555233184199 * cart[10 * ncart + i];
        tmp += -0.4528555233184199 * cart[21 * ncart + i];
        tmp += -7.2456883730947190 * cart[5 * ncart + i];
        tmp +=  7.2456883730947190 * cart[23 * ncart + i];
        tmp +=  7.2456883730947190 * cart[14 * ncart + i];
        tmp += -7.2456883730947190 * cart[25 * ncart + i];
        output[i] += tmp * vector[8];

    }

    // R_65c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -2.7171331399105196 * cart[2 * ncart + i];
        tmp +=  5.4342662798210393 * cart[7 * ncart + i];
        tmp +=  8.1513994197315593 * cart[16 * ncart + i];
        tmp +=  7.2456883730947190 * cart[9 * ncart + i];
        tmp += -21.7370651192841571 * cart[18 * ncart + i];
        output[i] += tmp * vector[9];

    }
    // R_65s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.4960783708246108 * cart[i];
        tmp +=  2.4803918541230536 * cart[3 * ncart + i];
        tmp +=  2.4803918541230536 * cart[10 * ncart + i];
        tmp += -0.4960783708246108 * cart[21 * ncart + i];
        tmp +=  4.9607837082461073 * cart[5 * ncart + i];
        tmp += -29.7647022494766453 * cart[12 * ncart + i];
        tmp +=  4.9607837082461073 * cart[23 * ncart + i];
        output[i] += tmp * vector[10];

    }

    // R_66c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.3268138086232857 * cart[2 * ncart + i];
        tmp += -23.2681380862328560 * cart[7 * ncart + i];
        tmp +=  11.6340690431164280 * cart[16 * ncart + i];
        output[i] += tmp * vector[11];

    }
    // R_66s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.6716932893813962 * cart[i];
        tmp += -10.0753993407209421 * cart[3 * ncart + i];
        tmp +=  10.0753993407209421 * cart[10 * ncart + i];
        tmp += -0.6716932893813962 * cart[21 * ncart + i];
        output[i] += tmp * vector[12];

    }

}
void gg_gaussian_cart_to_spherical_L0(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_00 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  = cart[i];

    }

}
void gg_gaussian_cart_to_spherical_sum_L0(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_00 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = cart[i];
        output[i] += tmp * vector[0];

    }

}
void gg_gaussian_cart_to_spherical_L1(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_10 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  = cart[2 * ncart + i];

    }

    // R_11c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i]  = cart[i];

    }
    // R_11s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i]  = cart[ncart + i];

    }

}
void gg_gaussian_cart_to_spherical_sum_L1(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_10 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = cart[2 * ncart + i];
        output[i] += tmp * vector[0];

    }

    // R_11c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = cart[i];
        output[i] += tmp * vector[1];

    }
    // R_11s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = cart[ncart + i];
        output[i] += tmp * vector[2];

    }

}
void gg_gaussian_cart_to_spherical_L2(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_20 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  = -0.5000000000000000 * cart[i];
        spherical[i] += -0.5000000000000000 * cart[3 * ncart + i];
        spherical[i] += cart[5 * ncart + i];

    }

    // R_21c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i]  =  1.7320508075688772 * cart[2 * ncart + i];

    }
    // R_21s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i]  =  1.7320508075688772 * cart[4 * ncart + i];

    }

    // R_22c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i]  =  0.8660254037844386 * cart[i];
        spherical[3 * nspherical + i] += -0.8660254037844386 * cart[3 * ncart + i];

    }
    // R_22s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i]  =  1.7320508075688772 * cart[ncart + i];

    }

}
void gg_gaussian_cart_to_spherical_sum_L2(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_20 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.5000000000000000 * cart[i];
        tmp += -0.5000000000000000 * cart[3 * ncart + i];
        tmp += cart[5 * ncart + i];
        output[i] += tmp * vector[0];

    }

    // R_21c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  1.7320508075688772 * cart[2 * ncart + i];
        output[i] += tmp * vector[1];

    }
    // R_21s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  1.7320508075688772 * cart[4 * ncart + i];
        output[i] += tmp * vector[2];

    }

    // R_22c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.8660254037844386 * cart[i];
        tmp += -0.8660254037844386 * cart[3 * ncart + i];
        output[i] += tmp * vector[3];

    }
    // R_22s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  1.7320508075688772 * cart[ncart + i];
        output[i] += tmp * vector[4];

    }

}
void gg_gaussian_cart_to_spherical_L3(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_30 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  = -1.5000000000000000 * cart[2 * ncart + i];
        spherical[i] += -1.5000000000000000 * cart[7 * ncart + i];
        spherical[i] += cart[9 * ncart + i];

    }

    // R_31c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i]  = -0.6123724356957945 * cart[i];
        spherical[nspherical + i] += -0.6123724356957945 * cart[3 * ncart + i];
        spherical[nspherical + i] +=  2.4494897427831779 * cart[5 * ncart + i];

    }
    // R_31s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i]  = -0.6123724356957945 * cart[ncart + i];
        spherical[2 * nspherical + i] += -0.6123724356957945 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] +=  2.4494897427831779 * cart[8 * ncart + i];

    }

    // R_32c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i]  =  1.9364916731037085 * cart[2 * ncart + i];
        spherical[3 * nspherical + i] += -1.9364916731037085 * cart[7 * ncart + i];

    }
    // R_32s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i]  =  3.8729833462074170 * cart[4 * ncart + i];

    }

    // R_33c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i]  =  0.7905694150420949 * cart[i];
        spherical[5 * nspherical + i] += -2.3717082451262845 * cart[3 * ncart + i];

    }
    // R_33s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i]  =  2.3717082451262845 * cart[ncart + i];
        spherical[6 * nspherical + i] += -0.7905694150420949 * cart[6 * ncart + i];

    }

}
void gg_gaussian_cart_to_spherical_sum_L3(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_30 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -1.5000000000000000 * cart[2 * ncart + i];
        tmp += -1.5000000000000000 * cart[7 * ncart + i];
        tmp += cart[9 * ncart + i];
        output[i] += tmp * vector[0];

    }

    // R_31c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.6123724356957945 * cart[i];
        tmp += -0.6123724356957945 * cart[3 * ncart + i];
        tmp +=  2.4494897427831779 * cart[5 * ncart + i];
        output[i] += tmp * vector[1];

    }
    // R_31s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.6123724356957945 * cart[ncart + i];
        tmp += -0.6123724356957945 * cart[6 * ncart + i];
        tmp +=  2.4494897427831779 * cart[8 * ncart + i];
        output[i] += tmp * vector[2];

    }

    // R_32c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  1.9364916731037085 * cart[2 * ncart + i];
        tmp += -1.9364916731037085 * cart[7 * ncart + i];
        output[i] += tmp * vector[3];

    }
    // R_32s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  3.8729833462074170 * cart[4 * ncart + i];
        output[i] += tmp * vector[4];

    }

    // R_33c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.7905694150420949 * cart[i];
        tmp += -2.3717082451262845 * cart[3 * ncart + i];
        output[i] += tmp * vector[5];

    }
    // R_33s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.3717082451262845 * cart[ncart + i];
        tmp += -0.7905694150420949 * cart[6 * ncart + i];
        output[i] += tmp * vector[6];

    }

}
void gg_gaussian_cart_to_spherical_L4(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_40 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  =  0.3750000000000000 * cart[i];
        spherical[i] +=  0.7500000000000000 * cart[3 * ncart + i];
        spherical[i] +=  0.3750000000000000 * cart[10 * ncart + i];
        spherical[i] += -3.0000000000000000 * cart[5 * ncart + i];
        spherical[i] += -3.0000000000000000 * cart[12 * ncart + i];
        spherical[i] += cart[14 * ncart + i];

    }

    // R_41c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i]  = -2.3717082451262845 * cart[2 * ncart + i];
        spherical[nspherical + i] += -2.3717082451262845 * cart[7 * ncart + i];
        spherical[nspherical + i] +=  3.1622776601683795 * cart[9 * ncart + i];

    }
    // R_41s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i]  = -2.3717082451262845 * cart[4 * ncart + i];
        spherical[2 * nspherical + i] += -2.3717082451262845 * cart[11 * ncart + i];
        spherical[2 * nspherical + i] +=  3.1622776601683795 * cart[13 * ncart + i];

    }

    // R_42c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i]  = -0.5590169943749475 * cart[i];
        spherical[3 * nspherical + i] +=  0.5590169943749475 * cart[10 * ncart + i];
        spherical[3 * nspherical + i] +=  3.3541019662496847 * cart[5 * ncart + i];
        spherical[3 * nspherical + i] += -3.3541019662496847 * cart[12 * ncart + i];

    }
    // R_42s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i]  = -1.1180339887498949 * cart[ncart + i];
        spherical[4 * nspherical + i] += -1.1180339887498949 * cart[6 * ncart + i];
        spherical[4 * nspherical + i] +=  6.7082039324993694 * cart[8 * ncart + i];

    }

    // R_43c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i]  =  2.0916500663351889 * cart[2 * ncart + i];
        spherical[5 * nspherical + i] += -6.2749501990055663 * cart[7 * ncart + i];

    }
    // R_43s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i]  =  6.2749501990055663 * cart[4 * ncart + i];
        spherical[6 * nspherical + i] += -2.0916500663351889 * cart[11 * ncart + i];

    }

    // R_44c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i]  =  0.7395099728874520 * cart[i];
        spherical[7 * nspherical + i] += -4.4370598373247123 * cart[3 * ncart + i];
        spherical[7 * nspherical + i] +=  0.7395099728874520 * cart[10 * ncart + i];

    }
    // R_44s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i]  =  2.9580398915498081 * cart[ncart + i];
        spherical[8 * nspherical + i] += -2.9580398915498081 * cart[6 * ncart + i];

    }

}
void gg_gaussian_cart_to_spherical_sum_L4(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_40 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.3750000000000000 * cart[i];
        tmp +=  0.7500000000000000 * cart[3 * ncart + i];
        tmp +=  0.3750000000000000 * cart[10 * ncart + i];
        tmp += -3.0000000000000000 * cart[5 * ncart + i];
        tmp += -3.0000000000000000 * cart[12 * ncart + i];
        tmp += cart[14 * ncart + i];
        output[i] += tmp * vector[0];

    }

    // R_41c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -2.3717082451262845 * cart[2 * ncart + i];
        tmp += -2.3717082451262845 * cart[7 * ncart + i];
        tmp +=  3.1622776601683795 * cart[9 * ncart + i];
        output[i] += tmp * vector[1];

    }
    // R_41s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -2.3717082451262845 * cart[4 * ncart + i];
        tmp += -2.3717082451262845 * cart[11 * ncart + i];
        tmp +=  3.1622776601683795 * cart[13 * ncart + i];
        output[i] += tmp * vector[2];

    }

    // R_42c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.5590169943749475 * cart[i];
        tmp +=  0.5590169943749475 * cart[10 * ncart + i];
        tmp +=  3.3541019662496847 * cart[5 * ncart + i];
        tmp += -3.3541019662496847 * cart[12 * ncart + i];
        output[i] += tmp * vector[3];

    }
    // R_42s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -1.1180339887498949 * cart[ncart + i];
        tmp += -1.1180339887498949 * cart[6 * ncart + i];
        tmp +=  6.7082039324993694 * cart[8 * ncart + i];
        output[i] += tmp * vector[4];

    }

    // R_43c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.0916500663351889 * cart[2 * ncart + i];
        tmp += -6.2749501990055663 * cart[7 * ncart + i];
        output[i] += tmp * vector[5];

    }
    // R_43s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  6.2749501990055663 * cart[4 * ncart + i];
        tmp += -2.0916500663351889 * cart[11 * ncart + i];
        output[i] += tmp * vector[6];

    }

    // R_44c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.7395099728874520 * cart[i];
        tmp += -4.4370598373247123 * cart[3 * ncart + i];
        tmp +=  0.7395099728874520 * cart[10 * ncart + i];
        output[i] += tmp * vector[7];

    }
    // R_44s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.9580398915498081 * cart[ncart + i];
        tmp += -2.9580398915498081 * cart[6 * ncart + i];
        output[i] += tmp * vector[8];

    }

}
void gg_gaussian_cart_to_spherical_L5(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_50 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  =  1.8750000000000000 * cart[2 * ncart + i];
        spherical[i] +=  3.7500000000000000 * cart[7 * ncart + i];
        spherical[i] +=  1.8750000000000000 * cart[16 * ncart + i];
        spherical[i] += -5.0000000000000000 * cart[9 * ncart + i];
        spherical[i] += -5.0000000000000000 * cart[18 * ncart + i];
        spherical[i] += cart[20 * ncart + i];

    }

    // R_51c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i]  =  0.4841229182759271 * cart[i];
        spherical[nspherical + i] +=  0.9682458365518543 * cart[3 * ncart + i];
        spherical[nspherical + i] +=  0.4841229182759271 * cart[10 * ncart + i];
        spherical[nspherical + i] += -5.8094750193111251 * cart[5 * ncart + i];
        spherical[nspherical + i] += -5.8094750193111251 * cart[12 * ncart + i];
        spherical[nspherical + i] +=  3.8729833462074170 * cart[14 * ncart + i];

    }
    // R_51s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i]  =  0.4841229182759271 * cart[ncart + i];
        spherical[2 * nspherical + i] +=  0.9682458365518543 * cart[6 * ncart + i];
        spherical[2 * nspherical + i] +=  0.4841229182759271 * cart[15 * ncart + i];
        spherical[2 * nspherical + i] += -5.8094750193111251 * cart[8 * ncart + i];
        spherical[2 * nspherical + i] += -5.8094750193111251 * cart[17 * ncart + i];
        spherical[2 * nspherical + i] +=  3.8729833462074170 * cart[19 * ncart + i];

    }

    // R_52c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i]  = -2.5617376914898995 * cart[2 * ncart + i];
        spherical[3 * nspherical + i] +=  2.5617376914898995 * cart[16 * ncart + i];
        spherical[3 * nspherical + i] +=  5.1234753829797990 * cart[9 * ncart + i];
        spherical[3 * nspherical + i] += -5.1234753829797990 * cart[18 * ncart + i];

    }
    // R_52s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i]  = -5.1234753829797990 * cart[4 * ncart + i];
        spherical[4 * nspherical + i] += -5.1234753829797990 * cart[11 * ncart + i];
        spherical[4 * nspherical + i] +=  10.2469507659595980 * cart[13 * ncart + i];

    }

    // R_53c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i]  = -0.5229125165837972 * cart[i];
        spherical[5 * nspherical + i] +=  1.0458250331675945 * cart[3 * ncart + i];
        spherical[5 * nspherical + i] +=  1.5687375497513916 * cart[10 * ncart + i];
        spherical[5 * nspherical + i] +=  4.1833001326703778 * cart[5 * ncart + i];
        spherical[5 * nspherical + i] += -12.5499003980111326 * cart[12 * ncart + i];

    }
    // R_53s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i]  = -1.5687375497513916 * cart[ncart + i];
        spherical[6 * nspherical + i] += -1.0458250331675945 * cart[6 * ncart + i];
        spherical[6 * nspherical + i] +=  0.5229125165837972 * cart[15 * ncart + i];
        spherical[6 * nspherical + i] +=  12.5499003980111326 * cart[8 * ncart + i];
        spherical[6 * nspherical + i] += -4.1833001326703778 * cart[17 * ncart + i];

    }

    // R_54c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i]  =  2.2185299186623562 * cart[2 * ncart + i];
        spherical[7 * nspherical + i] += -13.3111795119741370 * cart[7 * ncart + i];
        spherical[7 * nspherical + i] +=  2.2185299186623562 * cart[16 * ncart + i];

    }
    // R_54s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i]  =  8.8741196746494246 * cart[4 * ncart + i];
        spherical[8 * nspherical + i] += -8.8741196746494246 * cart[11 * ncart + i];

    }

    // R_55c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[9 * nspherical + i]  =  0.7015607600201140 * cart[i];
        spherical[9 * nspherical + i] += -7.0156076002011405 * cart[3 * ncart + i];
        spherical[9 * nspherical + i] +=  3.5078038001005702 * cart[10 * ncart + i];

    }
    // R_55s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[10 * nspherical + i]  =  3.5078038001005702 * cart[ncart + i];
        spherical[10 * nspherical + i] += -7.0156076002011405 * cart[6 * ncart + i];
        spherical[10 * nspherical + i] +=  0.7015607600201140 * cart[15 * ncart + i];

    }

}
void gg_gaussian_cart_to_spherical_sum_L5(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_50 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  1.8750000000000000 * cart[2 * ncart + i];
        tmp +=  3.7500000000000000 * cart[7 * ncart + i];
        tmp +=  1.8750000000000000 * cart[16 * ncart + i];
        tmp += -5.0000000000000000 * cart[9 * ncart + i];
        tmp += -5.0000000000000000 * cart[18 * ncart + i];
        tmp += cart[20 * ncart + i];
        output[i] += tmp * vector[0];

    }

    // R_51c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.4841229182759271 * cart[i];
        tmp +=  0.9682458365518543 * cart[3 * ncart + i];
        tmp +=  0.4841229182759271 * cart[10 * ncart + i];
        tmp += -5.8094750193111251 * cart[5 * ncart + i];
        tmp += -5.8094750193111251 * cart[12 * ncart + i];
        tmp +=  3.8729833462074170 * cart[14 * ncart + i];
        output[i] += tmp * vector[1];

    }
    // R_51s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.4841229182759271 * cart[ncart + i];
        tmp +=  0.9682458365518543 * cart[6 * ncart + i];
        tmp +=  0.4841229182759271 * cart[15 * ncart + i];
        tmp += -5.8094750193111251 * cart[8 * ncart + i];
        tmp += -5.8094750193111251 * cart[17 * ncart + i];
        tmp +=  3.8729833462074170 * cart[19 * ncart + i];
        output[i] += tmp * vector[2];

    }

    // R_52c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -2.5617376914898995 * cart[2 * ncart + i];
        tmp +=  2.5617376914898995 * cart[16 * ncart + i];
        tmp +=  5.1234753829797990 * cart[9 * ncart + i];
        tmp += -5.1234753829797990 * cart[18 * ncart + i];
        output[i] += tmp * vector[3];

    }
    // R_52s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -5.1234753829797990 * cart[4 * ncart + i];
        tmp += -5.1234753829797990 * cart[11 * ncart + i];
        tmp +=  10.2469507659595980 * cart[13 * ncart + i];
        output[i] += tmp * vector[4];

    }

    // R_53c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.5229125165837972 * cart[i];
        tmp +=  1.0458250331675945 * cart[3 * ncart + i];
        tmp +=  1.5687375497513916 * cart[10 * ncart + i];
        tmp +=  4.1833001326703778 * cart[5 * ncart + i];
        tmp += -12.5499003980111326 * cart[12 * ncart + i];
        output[i] += tmp * vector[5];

    }
    // R_53s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -1.5687375497513916 * cart[ncart + i];
        tmp += -1.0458250331675945 * cart[6 * ncart + i];
        tmp +=  0.5229125165837972 * cart[15 * ncart + i];
        tmp +=  12.5499003980111326 * cart[8 * ncart + i];
        tmp += -4.1833001326703778 * cart[17 * ncart + i];
        output[i] += tmp * vector[6];

    }

    // R_54c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.2185299186623562 * cart[2 * ncart + i];
        tmp += -13.3111795119741370 * cart[7 * ncart + i];
        tmp +=  2.2185299186623562 * cart[16 * ncart + i];
        output[i] += tmp * vector[7];

    }
    // R_54s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  8.8741196746494246 * cart[4 * ncart + i];
        tmp += -8.8741196746494246 * cart[11 * ncart + i];
        output[i] += tmp * vector[8];

    }

    // R_55c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.7015607600201140 * cart[i];
        tmp += -7.0156076002011405 * cart[3 * ncart + i];
        tmp +=  3.5078038001005702 * cart[10 * ncart + i];
        output[i] += tmp * vector[9];

    }
    // R_55s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  3.5078038001005702 * cart[ncart + i];
        tmp += -7.0156076002011405 * cart[6 * ncart + i];
        tmp +=  0.7015607600201140 * cart[15 * ncart + i];
        output[i] += tmp * vector[10];

    }

}
void gg_gaussian_cart_to_spherical_L6(const unsigned long size, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT spherical, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // R_60 Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[i]  = -0.3125000000000000 * cart[i];
        spherical[i] += -0.9375000000000000 * cart[3 * ncart + i];
        spherical[i] += -0.9375000000000000 * cart[10 * ncart + i];
        spherical[i] += -0.3125000000000000 * cart[21 * ncart + i];
        spherical[i] +=  5.6250000000000000 * cart[5 * ncart + i];
        spherical[i] +=  11.2500000000000000 * cart[12 * ncart + i];
        spherical[i] +=  5.6250000000000000 * cart[23 * ncart + i];
        spherical[i] += -7.5000000000000000 * cart[14 * ncart + i];
        spherical[i] += -7.5000000000000000 * cart[25 * ncart + i];
        spherical[i] += cart[27 * ncart + i];

    }

    // R_61c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[nspherical + i]  =  2.8641098093473998 * cart[2 * ncart + i];
        spherical[nspherical + i] +=  5.7282196186947996 * cart[7 * ncart + i];
        spherical[nspherical + i] +=  2.8641098093473998 * cart[16 * ncart + i];
        spherical[nspherical + i] += -11.4564392373895991 * cart[9 * ncart + i];
        spherical[nspherical + i] += -11.4564392373895991 * cart[18 * ncart + i];
        spherical[nspherical + i] +=  4.5825756949558398 * cart[20 * ncart + i];

    }
    // R_61s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[2 * nspherical + i]  =  2.8641098093473998 * cart[4 * ncart + i];
        spherical[2 * nspherical + i] +=  5.7282196186947996 * cart[11 * ncart + i];
        spherical[2 * nspherical + i] +=  2.8641098093473998 * cart[22 * ncart + i];
        spherical[2 * nspherical + i] += -11.4564392373895991 * cart[13 * ncart + i];
        spherical[2 * nspherical + i] += -11.4564392373895991 * cart[24 * ncart + i];
        spherical[2 * nspherical + i] +=  4.5825756949558398 * cart[26 * ncart + i];

    }

    // R_62c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[3 * nspherical + i]  =  0.4528555233184199 * cart[i];
        spherical[3 * nspherical + i] +=  0.4528555233184199 * cart[3 * ncart + i];
        spherical[3 * nspherical + i] += -0.4528555233184199 * cart[10 * ncart + i];
        spherical[3 * nspherical + i] += -0.4528555233184199 * cart[21 * ncart + i];
        spherical[3 * nspherical + i] += -7.2456883730947190 * cart[5 * ncart + i];
        spherical[3 * nspherical + i] +=  7.2456883730947190 * cart[23 * ncart + i];
        spherical[3 * nspherical + i] +=  7.2456883730947190 * cart[14 * ncart + i];
        spherical[3 * nspherical + i] += -7.2456883730947190 * cart[25 * ncart + i];

    }
    // R_62s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[4 * nspherical + i]  =  0.9057110466368399 * cart[ncart + i];
        spherical[4 * nspherical + i] +=  1.8114220932736798 * cart[6 * ncart + i];
        spherical[4 * nspherical + i] +=  0.9057110466368399 * cart[15 * ncart + i];
        spherical[4 * nspherical + i] += -14.4913767461894381 * cart[8 * ncart + i];
        spherical[4 * nspherical + i] += -14.4913767461894381 * cart[17 * ncart + i];
        spherical[4 * nspherical + i] +=  14.4913767461894381 * cart[19 * ncart + i];

    }

    // R_63c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[5 * nspherical + i]  = -2.7171331399105196 * cart[2 * ncart + i];
        spherical[5 * nspherical + i] +=  5.4342662798210393 * cart[7 * ncart + i];
        spherical[5 * nspherical + i] +=  8.1513994197315593 * cart[16 * ncart + i];
        spherical[5 * nspherical + i] +=  7.2456883730947190 * cart[9 * ncart + i];
        spherical[5 * nspherical + i] += -21.7370651192841571 * cart[18 * ncart + i];

    }
    // R_63s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[6 * nspherical + i]  = -8.1513994197315593 * cart[4 * ncart + i];
        spherical[6 * nspherical + i] += -5.4342662798210393 * cart[11 * ncart + i];
        spherical[6 * nspherical + i] +=  2.7171331399105196 * cart[22 * ncart + i];
        spherical[6 * nspherical + i] +=  21.7370651192841571 * cart[13 * ncart + i];
        spherical[6 * nspherical + i] += -7.2456883730947190 * cart[24 * ncart + i];

    }

    // R_64c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[7 * nspherical + i]  = -0.4960783708246108 * cart[i];
        spherical[7 * nspherical + i] +=  2.4803918541230536 * cart[3 * ncart + i];
        spherical[7 * nspherical + i] +=  2.4803918541230536 * cart[10 * ncart + i];
        spherical[7 * nspherical + i] += -0.4960783708246108 * cart[21 * ncart + i];
        spherical[7 * nspherical + i] +=  4.9607837082461073 * cart[5 * ncart + i];
        spherical[7 * nspherical + i] += -29.7647022494766453 * cart[12 * ncart + i];
        spherical[7 * nspherical + i] +=  4.9607837082461073 * cart[23 * ncart + i];

    }
    // R_64s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[8 * nspherical + i]  = -1.9843134832984430 * cart[ncart + i];
        spherical[8 * nspherical + i] +=  1.9843134832984430 * cart[15 * ncart + i];
        spherical[8 * nspherical + i] +=  19.8431348329844290 * cart[8 * ncart + i];
        spherical[8 * nspherical + i] += -19.8431348329844290 * cart[17 * ncart + i];

    }

    // R_65c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[9 * nspherical + i]  =  2.3268138086232857 * cart[2 * ncart + i];
        spherical[9 * nspherical + i] += -23.2681380862328560 * cart[7 * ncart + i];
        spherical[9 * nspherical + i] +=  11.6340690431164280 * cart[16 * ncart + i];

    }
    // R_65s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[10 * nspherical + i]  =  11.6340690431164280 * cart[4 * ncart + i];
        spherical[10 * nspherical + i] += -23.2681380862328560 * cart[11 * ncart + i];
        spherical[10 * nspherical + i] +=  2.3268138086232857 * cart[22 * ncart + i];

    }

    // R_66c Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[11 * nspherical + i]  =  0.6716932893813962 * cart[i];
        spherical[11 * nspherical + i] += -10.0753993407209421 * cart[3 * ncart + i];
        spherical[11 * nspherical + i] +=  10.0753993407209421 * cart[10 * ncart + i];
        spherical[11 * nspherical + i] += -0.6716932893813962 * cart[21 * ncart + i];

    }
    // R_66s Transform
    for (unsigned long i = 0; i < size; i++) {
        spherical[12 * nspherical + i]  =  4.0301597362883772 * cart[ncart + i];
        spherical[12 * nspherical + i] += -13.4338657876279228 * cart[6 * ncart + i];
        spherical[12 * nspherical + i] +=  4.0301597362883772 * cart[15 * ncart + i];

    }

}
void gg_gaussian_cart_to_spherical_sum_L6(const unsigned long size, const double* vector, const double* PRAGMA_RESTRICT cart, const unsigned long ncart, double* PRAGMA_RESTRICT output, const unsigned long nspherical) {
    ASSUME_ALIGNED(cart, 64);
    // temps
    double tmp;
    // R_60 Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.3125000000000000 * cart[i];
        tmp += -0.9375000000000000 * cart[3 * ncart + i];
        tmp += -0.9375000000000000 * cart[10 * ncart + i];
        tmp += -0.3125000000000000 * cart[21 * ncart + i];
        tmp +=  5.6250000000000000 * cart[5 * ncart + i];
        tmp +=  11.2500000000000000 * cart[12 * ncart + i];
        tmp +=  5.6250000000000000 * cart[23 * ncart + i];
        tmp += -7.5000000000000000 * cart[14 * ncart + i];
        tmp += -7.5000000000000000 * cart[25 * ncart + i];
        tmp += cart[27 * ncart + i];
        output[i] += tmp * vector[0];

    }

    // R_61c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.8641098093473998 * cart[2 * ncart + i];
        tmp +=  5.7282196186947996 * cart[7 * ncart + i];
        tmp +=  2.8641098093473998 * cart[16 * ncart + i];
        tmp += -11.4564392373895991 * cart[9 * ncart + i];
        tmp += -11.4564392373895991 * cart[18 * ncart + i];
        tmp +=  4.5825756949558398 * cart[20 * ncart + i];
        output[i] += tmp * vector[1];

    }
    // R_61s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.8641098093473998 * cart[4 * ncart + i];
        tmp +=  5.7282196186947996 * cart[11 * ncart + i];
        tmp +=  2.8641098093473998 * cart[22 * ncart + i];
        tmp += -11.4564392373895991 * cart[13 * ncart + i];
        tmp += -11.4564392373895991 * cart[24 * ncart + i];
        tmp +=  4.5825756949558398 * cart[26 * ncart + i];
        output[i] += tmp * vector[2];

    }

    // R_62c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.4528555233184199 * cart[i];
        tmp +=  0.4528555233184199 * cart[3 * ncart + i];
        tmp += -0.4528555233184199 * cart[10 * ncart + i];
        tmp += -0.4528555233184199 * cart[21 * ncart + i];
        tmp += -7.2456883730947190 * cart[5 * ncart + i];
        tmp +=  7.2456883730947190 * cart[23 * ncart + i];
        tmp +=  7.2456883730947190 * cart[14 * ncart + i];
        tmp += -7.2456883730947190 * cart[25 * ncart + i];
        output[i] += tmp * vector[3];

    }
    // R_62s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.9057110466368399 * cart[ncart + i];
        tmp +=  1.8114220932736798 * cart[6 * ncart + i];
        tmp +=  0.9057110466368399 * cart[15 * ncart + i];
        tmp += -14.4913767461894381 * cart[8 * ncart + i];
        tmp += -14.4913767461894381 * cart[17 * ncart + i];
        tmp +=  14.4913767461894381 * cart[19 * ncart + i];
        output[i] += tmp * vector[4];

    }

    // R_63c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -2.7171331399105196 * cart[2 * ncart + i];
        tmp +=  5.4342662798210393 * cart[7 * ncart + i];
        tmp +=  8.1513994197315593 * cart[16 * ncart + i];
        tmp +=  7.2456883730947190 * cart[9 * ncart + i];
        tmp += -21.7370651192841571 * cart[18 * ncart + i];
        output[i] += tmp * vector[5];

    }
    // R_63s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -8.1513994197315593 * cart[4 * ncart + i];
        tmp += -5.4342662798210393 * cart[11 * ncart + i];
        tmp +=  2.7171331399105196 * cart[22 * ncart + i];
        tmp +=  21.7370651192841571 * cart[13 * ncart + i];
        tmp += -7.2456883730947190 * cart[24 * ncart + i];
        output[i] += tmp * vector[6];

    }

    // R_64c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -0.4960783708246108 * cart[i];
        tmp +=  2.4803918541230536 * cart[3 * ncart + i];
        tmp +=  2.4803918541230536 * cart[10 * ncart + i];
        tmp += -0.4960783708246108 * cart[21 * ncart + i];
        tmp +=  4.9607837082461073 * cart[5 * ncart + i];
        tmp += -29.7647022494766453 * cart[12 * ncart + i];
        tmp +=  4.9607837082461073 * cart[23 * ncart + i];
        output[i] += tmp * vector[7];

    }
    // R_64s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  = -1.9843134832984430 * cart[ncart + i];
        tmp +=  1.9843134832984430 * cart[15 * ncart + i];
        tmp +=  19.8431348329844290 * cart[8 * ncart + i];
        tmp += -19.8431348329844290 * cart[17 * ncart + i];
        output[i] += tmp * vector[8];

    }

    // R_65c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  2.3268138086232857 * cart[2 * ncart + i];
        tmp += -23.2681380862328560 * cart[7 * ncart + i];
        tmp +=  11.6340690431164280 * cart[16 * ncart + i];
        output[i] += tmp * vector[9];

    }
    // R_65s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  11.6340690431164280 * cart[4 * ncart + i];
        tmp += -23.2681380862328560 * cart[11 * ncart + i];
        tmp +=  2.3268138086232857 * cart[22 * ncart + i];
        output[i] += tmp * vector[10];

    }

    // R_66c Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  0.6716932893813962 * cart[i];
        tmp += -10.0753993407209421 * cart[3 * ncart + i];
        tmp +=  10.0753993407209421 * cart[10 * ncart + i];
        tmp += -0.6716932893813962 * cart[21 * ncart + i];
        output[i] += tmp * vector[11];

    }
    // R_66s Transform
    for (unsigned long i = 0; i < size; i++) {
        tmp  =  4.0301597362883772 * cart[ncart + i];
        tmp += -13.4338657876279228 * cart[6 * ncart + i];
        tmp +=  4.0301597362883772 * cart[15 * ncart + i];
        output[i] += tmp * vector[12];

    }

}
void gg_cca_cart_copy_L0(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_sum_L0(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_copy_L1(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_sum_L1(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_copy_L2(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_sum_L2(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_copy_L3(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_sum_L3(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_copy_L4(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_sum_L4(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_copy_L5(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_sum_L5(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_copy_L6(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_cca_cart_sum_L6(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_molden_cart_copy_L0(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_molden_cart_sum_L0(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_molden_cart_copy_L1(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_molden_cart_sum_L1(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_molden_cart_copy_L2(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_molden_cart_sum_L2(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_molden_cart_copy_L3(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_molden_cart_sum_L3(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_molden_cart_copy_L4(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_molden_cart_sum_L4(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {

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
void gg_molden_cart_copy_L5(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
}
void gg_molden_cart_sum_L5(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
}
void gg_molden_cart_copy_L6(const unsigned long size, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
}
void gg_molden_cart_sum_L6(const unsigned long size, const double* PRAGMA_RESTRICT vector, const double* PRAGMA_RESTRICT cart_input, const unsigned long ncart_input, double* PRAGMA_RESTRICT cart_out, const unsigned long ncart_out) {
}
void gg_naive_transpose(unsigned long n, unsigned long m, const double* PRAGMA_RESTRICT input, double* PRAGMA_RESTRICT output) {
    ASSUME_ALIGNED(input, 64);
    for (unsigned long i = 0; i < n; i++) {
        for (unsigned long j = 0; j < m; j++) {
            output[j * n + i] = input[i * m + j];
        }
    }
}
void gg_fast_transpose(unsigned long n, unsigned long m, const double* PRAGMA_RESTRICT input, double* PRAGMA_RESTRICT output) {

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
void block_copy(unsigned long n, unsigned long m, const double* PRAGMA_RESTRICT input, unsigned long is, double* PRAGMA_RESTRICT output, unsigned long os, const int trans) {

    ASSUME_ALIGNED(input, 64);
    for (unsigned long i = 0; i < n; i++) {
        const unsigned long out_shift = i * os;
        const unsigned long inp_shift = i * is;

        for (unsigned long j = 0; j < m; j++) {
            output[out_shift + j] = input[inp_shift + j];
        }
    }
}
void block_matrix_vector(unsigned long n, unsigned long m, const double* vector, const double* PRAGMA_RESTRICT input, unsigned long is, double* PRAGMA_RESTRICT output) {

    ASSUME_ALIGNED(input, 64);
    for (unsigned long i = 0; i < n; i++) {
        const unsigned long inp_shift = i * is;
        const double coef = vector[i];

        for (unsigned long j = 0; j < m; j++) {
            output[j] += coef * input[inp_shift + j];
        }
    }
}