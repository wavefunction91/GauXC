#include <stdio.h>
#include <stdlib.h>
#include "integral_data_types.hpp"
#include "obara_saika_integrals.hpp"
#include "integral_0.hu"
#include "integral_1.hu"
#include "integral_2.hu"
#include "integral_0_0.hu"
#include "integral_1_0.hu"
#include "integral_1_1.hu"
#include "integral_2_0.hu"
#include "integral_2_1.hu"
#include "integral_2_2.hu"
namespace XGPU {

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

   prim_pair *prim_pairs = new prim_pair[np];
   for(int i = 0, ij = 0; i < nprim_A; ++i       )
   for(int j = 0        ; j < nprim_B; ++j, ++ij ) {
      auto& pair = prim_pairs[ij];
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

      pair.K_coeff_prod = 2 * M_PI * A.coeff[i].coeff * B.coeff[j].coeff * gamma_inv * std::exp( - alpha_A * alpha_B * dAB * gamma_inv );
   }

   AB.nprim_pair = np;
   AB.prim_pairs = prim_pairs;

}

void compute_integral_shell_pair(size_t npts,
                  int is_diag,
                  int lA,
                  int lB,
                  point rA,
                  point rB,
                  point rAB,
                  int nprim_pair,
                  prim_pair *ppair,
                  double *points,
                  double *Xi,
                  double *Xj,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int ldG, 
                  double *weights,
                  double *boys_table) {
   if (is_diag) {
      if(lA == 0) {
         integral_0<<<320, 128, 128 * 1 * sizeof(double)>>>(npts,
                                rA,
                                rB,
                                rAB,
                                nprim_pair,
                                ppair,
                                points,
                                Xi,
                                ldX,
                                Gi,
                                ldG, 
                                weights, 
                                boys_table);
      } else if(lA == 1) {
        integral_0<<<320, 128, 128 * 9 * sizeof(double)>>>(npts,
                               rA,
                               rB,
                               rAB,
                               nprim_pair,
                               ppair,
                               points,
                               Xi,
                               ldX,
                               Gi,
                               ldG, 
                               weights, 
                               boys_table);
      } else if(lA == 2) {
        integral_0<<<320, 128, 128 * 31 * sizeof(double)>>>(npts,
                               rA,
                               rB,
                               rAB,
                               nprim_pair,
                               ppair,
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
      if((lA == 0) && (lB == 0)) {
         integral_0_0<<<320, 128, 128 * 1 * sizeof(double)>>>(npts,
                                  rA,
                                  rB,
                                  rAB,
                                  nprim_pair,
                                  ppair,
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
         integral_1_0<<<320, 128, 128 * 3 * sizeof(double)>>>(npts,
                                  rA,
                                  rB,
                                  rAB,
                                  nprim_pair,
                                  ppair,
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
         integral_1_0<<<320, 128, 128 * 3 * sizeof(double)>>>(npts,
                                  rA,
                                  rB,
                                  rAB,
                                  nprim_pair,
                                  ppair,
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
        integral_1_1<<<320, 128, 128 * 9 * sizeof(double)>>>(npts,
                                 rA,
                                 rB,
                                 rAB,
                                 nprim_pair,
                                 ppair,
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
         integral_2_0<<<320, 128, 128 * 6 * sizeof(double)>>>(npts,
                                  rA,
                                  rB,
                                  rAB,
                                  nprim_pair,
                                  ppair,
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
         integral_2_0<<<320, 128, 128 * 6 * sizeof(double)>>>(npts,
                                  rA,
                                  rB,
                                  rAB,
                                  nprim_pair,
                                  ppair,
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
         integral_2_1<<<320, 128, 128 * 16 * sizeof(double)>>>(npts,
                                  rA,
                                  rB,
                                  rAB,
                                  nprim_pair,
                                  ppair,
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
         integral_2_1<<<320, 128, 128 * 16 * sizeof(double)>>>(npts,
                                  rA,
                                  rB,
                                  rAB,
                                  nprim_pair,
                                  ppair,
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
        integral_2_2<<<320, 128, 128 * 31 * sizeof(double)>>>(npts,
                                 rA,
                                 rB,
                                 rAB,
                                 nprim_pair,
                                 ppair,
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

}
}
