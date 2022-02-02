#include <stdio.h>
#include <stdlib.h>
#include "integral_data_types.hpp"
#include "obara_saika_integrals.hpp"
#include "integral_0.hpp"
#include "integral_1.hpp"
#include "integral_2.hpp"
#include "integral_3.hpp"
#include "integral_4.hpp"
#include "integral_0_0.hpp"
#include "integral_1_0.hpp"
#include "integral_1_1.hpp"
#include "integral_2_0.hpp"
#include "integral_2_1.hpp"
#include "integral_2_2.hpp"
#include "integral_3_0.hpp"
#include "integral_3_1.hpp"
#include "integral_3_2.hpp"
#include "integral_3_3.hpp"
#include "integral_4_0.hpp"
#include "integral_4_1.hpp"
#include "integral_4_2.hpp"
#include "integral_4_3.hpp"
#include "integral_4_4.hpp"
namespace XCPU {
void generate_shell_pair( const shells& A, const shells& B, shell_pair& AB) {
   // L Values
   AB.lA = A.L;
   AB.lB = B.L;

   AB.rA = A.origin;
   AB.rB = B.origin;

   const auto xA = A.origin.x;
   const auto yA = A.origin.y;
   const auto zA = A.origin.z;

   const auto xB = B.origin.x;
   const auto yB = B.origin.y;
   const auto zB = B.origin.z;

   AB.rAB.x = xA - xB;
   AB.rAB.y = yA - yB;
   AB.rAB.z = zA - zB;

   const double dAB = AB.rAB.x*AB.rAB.x + AB.rAB.y*AB.rAB.y + AB.rAB.z*AB.rAB.z;

   const int nprim_A = A.m;
   const int nprim_B = B.m;
   const int np = nprim_A * nprim_B;

   AB.nprim_pair = np;
   AB.prim_pairs = new prim_pair[np];
   for(int i = 0, ij = 0; i < nprim_A; ++i       )
   for(int j = 0        ; j < nprim_B; ++j, ++ij ) {
      auto& pair = AB.prim_pairs[ij];
      pair.coeff_prod = A.coeff[i].coeff * B.coeff[j].coeff;

      const auto alpha_A = A.coeff[i].alpha;
      const auto alpha_B = B.coeff[j].alpha;

      pair.gamma = alpha_A + alpha_B;
      const auto gamma_inv = 1. / pair.gamma;

      pair.P.x = (alpha_A * xA + alpha_B * xB) * gamma_inv;
      pair.P.y = (alpha_A * yA + alpha_B * yB) * gamma_inv;
      pair.P.z = (alpha_A * zA + alpha_B * zB) * gamma_inv;

      pair.PA.x = pair.P.x - xA;
      pair.PA.y = pair.P.y - yA;
      pair.PA.z = pair.P.z - zA;

      pair.PB.x = pair.P.x - xB;
      pair.PB.y = pair.P.y - yB;
      pair.PB.z = pair.P.z - zB;

      pair.K = 2 * M_PI * gamma_inv * std::exp( - alpha_A * alpha_B * dAB * gamma_inv );
   }
}

void compute_integral_shell_pair(size_t npts,
                  int i,
                  int j,
                  shells *shell_list,
                  double *points,
                  double *Xi,
                  double *Xj,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int ldG, 
                  double *weights, 
                  double *boys_table) {
   shell_pair shpair;
   // Account for permutational symmetry in kernels
   if( shell_list[i].L >= shell_list[j].L )
     generate_shell_pair(shell_list[i], shell_list[j], shpair);
   else
     generate_shell_pair(shell_list[j], shell_list[i], shpair);

   if (i == j) {
      int lA = shell_list[i].L;

      if(lA == 0) {
         integral_0(npts,
                    shpair,
                    points,
                    Xi,
                    ldX,
                    Gi,
                    ldG, 
                    weights, 
                    boys_table);
      } else if(lA == 1) {
        integral_1(npts,
                   shpair,
                   points,
                   Xi,
                   ldX,
                   Gi,
                   ldG, 
                   weights, 
                   boys_table);
      } else if(lA == 2) {
        integral_2(npts,
                   shpair,
                   points,
                   Xi,
                   ldX,
                   Gi,
                   ldG, 
                   weights, 
                   boys_table);
      } else if(lA == 3) {
        integral_3(npts,
                   shpair,
                   points,
                   Xi,
                   ldX,
                   Gi,
                   ldG, 
                   weights, 
                   boys_table);
      } else if(lA == 4) {
        integral_4(npts,
                   shpair,
                   points,
                   Xi,
                   ldX,
                   Gi,
                   ldG, 
                   weights, 
                   boys_table);
      } else {
         printf("Type not defined!\n");
      }
   } else {
      int lA = shell_list[i].L;
      int lB = shell_list[j].L;

      if((lA == 0) && (lB == 0)) {
         integral_0_0(npts,
                      shpair,
                      points,
                      Xi,
                      Xj,
                      ldX,
                      Gi,
                      Gj,
                      ldG, 
                      weights, 
                      boys_table);
      } else if((lA == 1) && (lB == 0)) {
            integral_1_0(npts,
                         shpair,
                         points,
                         Xi,
                         Xj,
                         ldX,
                         Gi,
                         Gj,
                         ldG, 
                         weights, 
                         boys_table);
      } else if((lA == 0) && (lB == 1)) {
         integral_1_0(npts,
                      shpair,
                      points,
                      Xj,
                      Xi,
                      ldX,
                      Gj,
                      Gi,
                      ldG, 
                      weights, 
                      boys_table);
      } else if((lA == 1) && (lB == 1)) {
        integral_1_1(npts,
                     shpair,
                     points,
                     Xi,
                     Xj,
                     ldX,
                     Gi,
                     Gj,
                     ldG, 
                     weights, 
                     boys_table);
      } else if((lA == 2) && (lB == 0)) {
            integral_2_0(npts,
                         shpair,
                         points,
                         Xi,
                         Xj,
                         ldX,
                         Gi,
                         Gj,
                         ldG, 
                         weights, 
                         boys_table);
      } else if((lA == 0) && (lB == 2)) {
         integral_2_0(npts,
                      shpair,
                      points,
                      Xj,
                      Xi,
                      ldX,
                      Gj,
                      Gi,
                      ldG, 
                      weights, 
                      boys_table);
      } else if((lA == 2) && (lB == 1)) {
            integral_2_1(npts,
                         shpair,
                         points,
                         Xi,
                         Xj,
                         ldX,
                         Gi,
                         Gj,
                         ldG, 
                         weights, 
                         boys_table);
      } else if((lA == 1) && (lB == 2)) {
         integral_2_1(npts,
                      shpair,
                      points,
                      Xj,
                      Xi,
                      ldX,
                      Gj,
                      Gi,
                      ldG, 
                      weights, 
                      boys_table);
      } else if((lA == 2) && (lB == 2)) {
        integral_2_2(npts,
                     shpair,
                     points,
                     Xi,
                     Xj,
                     ldX,
                     Gi,
                     Gj,
                     ldG, 
                     weights, 
                     boys_table);
      } else if((lA == 3) && (lB == 0)) {
            integral_3_0(npts,
                         shpair,
                         points,
                         Xi,
                         Xj,
                         ldX,
                         Gi,
                         Gj,
                         ldG, 
                         weights, 
                         boys_table);
      } else if((lA == 0) && (lB == 3)) {
         integral_3_0(npts,
                      shpair,
                      points,
                      Xj,
                      Xi,
                      ldX,
                      Gj,
                      Gi,
                      ldG, 
                      weights, 
                      boys_table);
      } else if((lA == 3) && (lB == 1)) {
            integral_3_1(npts,
                         shpair,
                         points,
                         Xi,
                         Xj,
                         ldX,
                         Gi,
                         Gj,
                         ldG, 
                         weights, 
                         boys_table);
      } else if((lA == 1) && (lB == 3)) {
         integral_3_1(npts,
                      shpair,
                      points,
                      Xj,
                      Xi,
                      ldX,
                      Gj,
                      Gi,
                      ldG, 
                      weights, 
                      boys_table);
      } else if((lA == 3) && (lB == 2)) {
            integral_3_2(npts,
                         shpair,
                         points,
                         Xi,
                         Xj,
                         ldX,
                         Gi,
                         Gj,
                         ldG, 
                         weights, 
                         boys_table);
      } else if((lA == 2) && (lB == 3)) {
         integral_3_2(npts,
                      shpair,
                      points,
                      Xj,
                      Xi,
                      ldX,
                      Gj,
                      Gi,
                      ldG, 
                      weights, 
                      boys_table);
      } else if((lA == 3) && (lB == 3)) {
        integral_3_3(npts,
                     shpair,
                     points,
                     Xi,
                     Xj,
                     ldX,
                     Gi,
                     Gj,
                     ldG, 
                     weights, 
                     boys_table);
      } else if((lA == 4) && (lB == 0)) {
            integral_4_0(npts,
                         shpair,
                         points,
                         Xi,
                         Xj,
                         ldX,
                         Gi,
                         Gj,
                         ldG, 
                         weights, 
                         boys_table);
      } else if((lA == 0) && (lB == 4)) {
         integral_4_0(npts,
                      shpair,
                      points,
                      Xj,
                      Xi,
                      ldX,
                      Gj,
                      Gi,
                      ldG, 
                      weights, 
                      boys_table);
      } else if((lA == 4) && (lB == 1)) {
            integral_4_1(npts,
                         shpair,
                         points,
                         Xi,
                         Xj,
                         ldX,
                         Gi,
                         Gj,
                         ldG, 
                         weights, 
                         boys_table);
      } else if((lA == 1) && (lB == 4)) {
         integral_4_1(npts,
                      shpair,
                      points,
                      Xj,
                      Xi,
                      ldX,
                      Gj,
                      Gi,
                      ldG, 
                      weights, 
                      boys_table);
      } else if((lA == 4) && (lB == 2)) {
            integral_4_2(npts,
                         shpair,
                         points,
                         Xi,
                         Xj,
                         ldX,
                         Gi,
                         Gj,
                         ldG, 
                         weights, 
                         boys_table);
      } else if((lA == 2) && (lB == 4)) {
         integral_4_2(npts,
                      shpair,
                      points,
                      Xj,
                      Xi,
                      ldX,
                      Gj,
                      Gi,
                      ldG, 
                      weights, 
                      boys_table);
      } else if((lA == 4) && (lB == 3)) {
            integral_4_3(npts,
                         shpair,
                         points,
                         Xi,
                         Xj,
                         ldX,
                         Gi,
                         Gj,
                         ldG, 
                         weights, 
                         boys_table);
      } else if((lA == 3) && (lB == 4)) {
         integral_4_3(npts,
                      shpair,
                      points,
                      Xj,
                      Xi,
                      ldX,
                      Gj,
                      Gi,
                      ldG, 
                      weights, 
                      boys_table);
      } else if((lA == 4) && (lB == 4)) {
        integral_4_4(npts,
                     shpair,
                     points,
                     Xi,
                     Xj,
                     ldX,
                     Gi,
                     Gj,
                     ldG, 
                     weights, 
                     boys_table);
      } else {
         printf("Type not defined!\n");
      }
   }
  delete shpair.prim_pairs;
}
}
