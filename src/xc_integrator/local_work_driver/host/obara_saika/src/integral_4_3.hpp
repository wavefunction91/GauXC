#ifndef __MY_INTEGRAL_4_3
#define __MY_INTEGRAL_4_3

#include "../include/cpu/integral_data_types.hpp"
namespace XCPU {
void integral_4_3(size_t npts,
                  double *points,
                  point rA,
                  point rB,
                  int nprim_pairs,
                  prim_pair *prim_pairs,
                  double *Xi,
                  double *Xj,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int ldG, 
                  double *weights, 
                  double *boys_table);
}

#endif
