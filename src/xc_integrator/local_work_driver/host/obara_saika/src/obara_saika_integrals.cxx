#include <stdio.h>
#include <stdlib.h>
#include "integral_data_types.h"
#include "obara_saika_integrals.h"
#include "integral_0_0.h"
#include "integral_1_0.h"
#include "integral_1_1.h"
#include "integral_2_0.h"
#include "integral_2_1.h"
#include "integral_2_2.h"
#include "integral_3_0.h"
#include "integral_3_1.h"
#include "integral_3_2.h"
#include "integral_3_3.h"
#include "integral_4_0.h"
#include "integral_4_1.h"
#include "integral_4_2.h"
#include "integral_4_3.h"
#include "integral_4_4.h"

void compute_integral_shell_pair(int npts,
                  shells shellA,
                  shells shellB,
                  point *points,
                  double *Xi,
                  double *Xj,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int ldG, 
                  double *weights) {
   int lA = shellA.L;
   int lB = shellB.L;

   if((lA == 0) && (lB == 0)) {
      integral_0_0(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   ldX,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 1) && (lB == 0)) {
      integral_1_0(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   Xj,
                   ldX,
                   Gi,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 1) && (lB == 1)) {
      integral_1_1(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   ldX,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 2) && (lB == 0)) {
      integral_2_0(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   Xj,
                   ldX,
                   Gi,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 2) && (lB == 1)) {
      integral_2_1(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   Xj,
                   ldX,
                   Gi,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 2) && (lB == 2)) {
      integral_2_2(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   ldX,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 3) && (lB == 0)) {
      integral_3_0(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   Xj,
                   ldX,
                   Gi,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 3) && (lB == 1)) {
      integral_3_1(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   Xj,
                   ldX,
                   Gi,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 3) && (lB == 2)) {
      integral_3_2(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   Xj,
                   ldX,
                   Gi,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 3) && (lB == 3)) {
      integral_3_3(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   ldX,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 4) && (lB == 0)) {
      integral_4_0(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   Xj,
                   ldX,
                   Gi,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 4) && (lB == 1)) {
      integral_4_1(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   Xj,
                   ldX,
                   Gi,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 4) && (lB == 2)) {
      integral_4_2(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   Xj,
                   ldX,
                   Gi,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 4) && (lB == 3)) {
      integral_4_3(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   Xj,
                   ldX,
                   Gi,
                   Gj,
                   ldG, 
                   weights);
   } else if((lA == 4) && (lB == 4)) {
      integral_4_4(npts,
                   shellA,
                   shellB,
                   points,
                   Xi,
                   ldX,
                   Gj,
                   ldG, 
                   weights);
   } else {
      printf("Type not defined!\n");
   }
}
