#include <stdio.h>
#include <stdlib.h>
#include "integral_data_types.h"
#include "obara_saika_integrals.h"
#include "integral_0.h"
#include "integral_1.h"
#include "integral_2.h"
#include "integral_3.h"
#include "integral_4.h"
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
                  int i,
                  int j,
                  shells *shell_list,
                  point *points,
                  double *Xi,
                  double *Xj,
                  int stX,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int stG, 
                  int ldG, 
                  double *weights) {
   if (i == j) {
      int lA = shell_list[i].L;

      if(lA == 0) {
         integral_0(npts,
                    shell_list[i],
                    points,
                    Xi,
                    stX,
                    ldX,
                    Gi,
                    stG, 
                    ldG, 
                    weights);
      } else if(lA == 1) {
        integral_1(npts,
                   shell_list[i],
                   points,
                   Xi,
                   stX,
                   ldX,
                   Gi,
                   stG, 
                   ldG, 
                   weights);
      } else if(lA == 2) {
        integral_2(npts,
                   shell_list[i],
                   points,
                   Xi,
                   stX,
                   ldX,
                   Gi,
                   stG, 
                   ldG, 
                   weights);
      } else if(lA == 3) {
        integral_3(npts,
                   shell_list[i],
                   points,
                   Xi,
                   stX,
                   ldX,
                   Gi,
                   stG, 
                   ldG, 
                   weights);
      } else if(lA == 4) {
        integral_4(npts,
                   shell_list[i],
                   points,
                   Xi,
                   stX,
                   ldX,
                   Gi,
                   stG, 
                   ldG, 
                   weights);
      } else {
         printf("Type not defined!\n");
      }
   } else {
      int lA = shell_list[i].L;
      int lB = shell_list[j].L;

      if((lA == 0) && (lB == 0)) {
         integral_0_0(npts,
                      shell_list[i],
                      shell_list[j],
                      points,
                      Xi,
                      Xj,
                      stX,
                      ldX,
                      Gi,
                      Gj,
                      stG, 
                      ldG, 
                      weights);
      } else if((lA == 1) && (lB == 0)) {
            integral_1_0(npts,
                         shell_list[i],
                         shell_list[j],
                         points,
                         Xi,
                         Xj,
                         stX,
                         ldX,
                         Gi,
                         Gj,
                         stG, 
                         ldG, 
                         weights);
      } else if((lA == 0) && (lB == 1)) {
         integral_1_0(npts,
                      shell_list[j],
                      shell_list[i],
                      points,
                      Xj,
                      Xi,
                      stX,
                      ldX,
                      Gj,
                      Gi,
                      stG, 
                      ldG, 
                      weights);
      } else if((lA == 1) && (lB == 1)) {
        integral_1_1(npts,
                     shell_list[i],
                     shell_list[j],
                     points,
                     Xi,
                     Xj,
                     stX,
                     ldX,
                     Gi,
                     Gj,
                     stG, 
                     ldG, 
                     weights);
      } else if((lA == 2) && (lB == 0)) {
            integral_2_0(npts,
                         shell_list[i],
                         shell_list[j],
                         points,
                         Xi,
                         Xj,
                         stX,
                         ldX,
                         Gi,
                         Gj,
                         stG, 
                         ldG, 
                         weights);
      } else if((lA == 0) && (lB == 2)) {
         integral_2_0(npts,
                      shell_list[j],
                      shell_list[i],
                      points,
                      Xj,
                      Xi,
                      stX,
                      ldX,
                      Gj,
                      Gi,
                      stG, 
                      ldG, 
                      weights);
      } else if((lA == 2) && (lB == 1)) {
            integral_2_1(npts,
                         shell_list[i],
                         shell_list[j],
                         points,
                         Xi,
                         Xj,
                         stX,
                         ldX,
                         Gi,
                         Gj,
                         stG, 
                         ldG, 
                         weights);
      } else if((lA == 1) && (lB == 2)) {
         integral_2_1(npts,
                      shell_list[j],
                      shell_list[i],
                      points,
                      Xj,
                      Xi,
                      stX,
                      ldX,
                      Gj,
                      Gi,
                      stG, 
                      ldG, 
                      weights);
      } else if((lA == 2) && (lB == 2)) {
        integral_2_2(npts,
                     shell_list[i],
                     shell_list[j],
                     points,
                     Xi,
                     Xj,
                     stX,
                     ldX,
                     Gi,
                     Gj,
                     stG, 
                     ldG, 
                     weights);
      } else if((lA == 3) && (lB == 0)) {
            integral_3_0(npts,
                         shell_list[i],
                         shell_list[j],
                         points,
                         Xi,
                         Xj,
                         stX,
                         ldX,
                         Gi,
                         Gj,
                         stG, 
                         ldG, 
                         weights);
      } else if((lA == 0) && (lB == 3)) {
         integral_3_0(npts,
                      shell_list[j],
                      shell_list[i],
                      points,
                      Xj,
                      Xi,
                      stX,
                      ldX,
                      Gj,
                      Gi,
                      stG, 
                      ldG, 
                      weights);
      } else if((lA == 3) && (lB == 1)) {
            integral_3_1(npts,
                         shell_list[i],
                         shell_list[j],
                         points,
                         Xi,
                         Xj,
                         stX,
                         ldX,
                         Gi,
                         Gj,
                         stG, 
                         ldG, 
                         weights);
      } else if((lA == 1) && (lB == 3)) {
         integral_3_1(npts,
                      shell_list[j],
                      shell_list[i],
                      points,
                      Xj,
                      Xi,
                      stX,
                      ldX,
                      Gj,
                      Gi,
                      stG, 
                      ldG, 
                      weights);
      } else if((lA == 3) && (lB == 2)) {
            integral_3_2(npts,
                         shell_list[i],
                         shell_list[j],
                         points,
                         Xi,
                         Xj,
                         stX,
                         ldX,
                         Gi,
                         Gj,
                         stG, 
                         ldG, 
                         weights);
      } else if((lA == 2) && (lB == 3)) {
         integral_3_2(npts,
                      shell_list[j],
                      shell_list[i],
                      points,
                      Xj,
                      Xi,
                      stX,
                      ldX,
                      Gj,
                      Gi,
                      stG, 
                      ldG, 
                      weights);
      } else if((lA == 3) && (lB == 3)) {
        integral_3_3(npts,
                     shell_list[i],
                     shell_list[j],
                     points,
                     Xi,
                     Xj,
                     stX,
                     ldX,
                     Gi,
                     Gj,
                     stG, 
                     ldG, 
                     weights);
      } else if((lA == 4) && (lB == 0)) {
            integral_4_0(npts,
                         shell_list[i],
                         shell_list[j],
                         points,
                         Xi,
                         Xj,
                         stX,
                         ldX,
                         Gi,
                         Gj,
                         stG, 
                         ldG, 
                         weights);
      } else if((lA == 0) && (lB == 4)) {
         integral_4_0(npts,
                      shell_list[j],
                      shell_list[i],
                      points,
                      Xj,
                      Xi,
                      stX,
                      ldX,
                      Gj,
                      Gi,
                      stG, 
                      ldG, 
                      weights);
      } else if((lA == 4) && (lB == 1)) {
            integral_4_1(npts,
                         shell_list[i],
                         shell_list[j],
                         points,
                         Xi,
                         Xj,
                         stX,
                         ldX,
                         Gi,
                         Gj,
                         stG, 
                         ldG, 
                         weights);
      } else if((lA == 1) && (lB == 4)) {
         integral_4_1(npts,
                      shell_list[j],
                      shell_list[i],
                      points,
                      Xj,
                      Xi,
                      stX,
                      ldX,
                      Gj,
                      Gi,
                      stG, 
                      ldG, 
                      weights);
      } else if((lA == 4) && (lB == 2)) {
            integral_4_2(npts,
                         shell_list[i],
                         shell_list[j],
                         points,
                         Xi,
                         Xj,
                         stX,
                         ldX,
                         Gi,
                         Gj,
                         stG, 
                         ldG, 
                         weights);
      } else if((lA == 2) && (lB == 4)) {
         integral_4_2(npts,
                      shell_list[j],
                      shell_list[i],
                      points,
                      Xj,
                      Xi,
                      stX,
                      ldX,
                      Gj,
                      Gi,
                      stG, 
                      ldG, 
                      weights);
      } else if((lA == 4) && (lB == 3)) {
            integral_4_3(npts,
                         shell_list[i],
                         shell_list[j],
                         points,
                         Xi,
                         Xj,
                         stX,
                         ldX,
                         Gi,
                         Gj,
                         stG, 
                         ldG, 
                         weights);
      } else if((lA == 3) && (lB == 4)) {
         integral_4_3(npts,
                      shell_list[j],
                      shell_list[i],
                      points,
                      Xj,
                      Xi,
                      stX,
                      ldX,
                      Gj,
                      Gi,
                      stG, 
                      ldG, 
                      weights);
      } else if((lA == 4) && (lB == 4)) {
        integral_4_4(npts,
                     shell_list[i],
                     shell_list[j],
                     points,
                     Xi,
                     Xj,
                     stX,
                     ldX,
                     Gi,
                     Gj,
                     stG, 
                     ldG, 
                     weights);
      } else {
         printf("Type not defined!\n");
      }
   }
}
