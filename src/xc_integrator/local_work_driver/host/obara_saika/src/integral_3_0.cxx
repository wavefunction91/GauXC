#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

void integral_3_0(size_t npts,
                  shell_pair shpair,
                  point *_points,
                  double *Xi,
                  double *Xj,
                  int stX,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int stG, 
                  int ldG, 
                  double *weights) {
   double temp[10 * NPTS_LOCAL];

   for(int i = 0; i < 10 * NPTS_LOCAL; ++i) {
      temp[i] = 0.0;
   }

   double X_AB = shpair.rAB.x;
   double Y_AB = shpair.rAB.y;
   double Z_AB = shpair.rAB.z;

   for(size_t p_outer = 0; p_outer < npts; p_outer += NPTS_LOCAL) {
      point *_point_outer = (_points + p_outer);

      for(int ij = 0; ij < shpair.nprim_pair; ++ij ) {
         double RHO = shpair.prim_pairs[ij].gamma;
         double RHO_INV = 1.0 / RHO;

         double xP = shpair.prim_pairs[ij].P.x;
         double yP = shpair.prim_pairs[ij].P.y;
         double zP = shpair.prim_pairs[ij].P.z;

         double X_PA = shpair.prim_pairs[ij].PA.x;
         double Y_PA = shpair.prim_pairs[ij].PA.y;
         double Z_PA = shpair.prim_pairs[ij].PA.z;

         double eval = shpair.prim_pairs[ij].coeff_prod * shpair.prim_pairs[ij].K;

         for(int p_inner = 0; p_inner < NPTS_LOCAL; ++p_inner) {
            point C = *(_point_outer + p_inner);

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            double X_PC = (xP - xC);
            double Y_PC = (yP - yC);
            double Z_PC = (zP - zC);

            double t00, t01, t02, t03, t10, t11, t12, t20, t21, t30;

            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

            t00 = eval * boys_function(0, tval);
            t01 = eval * boys_function(1, tval);
            t02 = eval * boys_function(2, tval);
            t03 = eval * boys_function(3, tval);
            t10 = X_PA * t00 - X_PC * t01;
            t11 = X_PA * t01 - X_PC * t02;
            t12 = X_PA * t02 - X_PC * t03;
            t20 = X_PA * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = X_PA * t11 - X_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t30 = X_PA * t20 - X_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            *(temp + 0 * NPTS_LOCAL + p_inner) += t30;
            t30 = Y_PA * t20 - Y_PC * t21;
            *(temp + 1 * NPTS_LOCAL + p_inner) += t30;
            t30 = Z_PA * t20 - Z_PC * t21;
            *(temp + 2 * NPTS_LOCAL + p_inner) += t30;
            t20 = Y_PA * t10 - Y_PC * t11;
            t21 = Y_PA * t11 - Y_PC * t12;
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            *(temp + 3 * NPTS_LOCAL + p_inner) += t30;
            t30 = Z_PA * t20 - Z_PC * t21;
            *(temp + 4 * NPTS_LOCAL + p_inner) += t30;
            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            *(temp + 5 * NPTS_LOCAL + p_inner) += t30;
            t10 = Y_PA * t00 - Y_PC * t01;
            t11 = Y_PA * t01 - Y_PC * t02;
            t12 = Y_PA * t02 - Y_PC * t03;
            t20 = Y_PA * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Y_PA * t11 - Y_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            *(temp + 6 * NPTS_LOCAL + p_inner) += t30;
            t30 = Z_PA * t20 - Z_PC * t21;
            *(temp + 7 * NPTS_LOCAL + p_inner) += t30;
            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            *(temp + 8 * NPTS_LOCAL + p_inner) += t30;
            t10 = Z_PA * t00 - Z_PC * t01;
            t11 = Z_PA * t01 - Z_PC * t02;
            t12 = Z_PA * t02 - Z_PC * t03;
            t20 = Z_PA * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Z_PA * t11 - Z_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            *(temp + 9 * NPTS_LOCAL + p_inner) += t30;
         }
      }

      for(int p_inner = 0; p_inner < NPTS_LOCAL; ++p_inner) {
         double *Xik = (Xi + (NPTS_LOCAL * p_outer + p_inner) * stX);
         double *Xjk = (Xj + (NPTS_LOCAL * p_outer + p_inner) * stX);
         double *Gik = (Gi + (NPTS_LOCAL * p_outer + p_inner) * stG);
         double *Gjk = (Gj + (NPTS_LOCAL * p_outer + p_inner) * stG);

         double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k, rcp_i, rcp_j, rcp_k;
         double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;

         X_ABp = 1.0; comb_m_i = 1.0;
         Y_ABp = 1.0; comb_n_j = 1.0;
         Z_ABp = 1.0; comb_p_k = 1.0;
         const_value = *(weights + p_outer * NPTS_LOCAL + p_inner) * comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
         t0 = *(temp + 0 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 0 * ldG) += *(Xjk + 0 * ldX) * t0;
         *(Gjk + 0 * ldG) += *(Xik + 0 * ldX) * t0;
         t1 = *(temp + 1 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 1 * ldG) += *(Xjk + 0 * ldX) * t1;
         *(Gjk + 0 * ldG) += *(Xik + 1 * ldX) * t1;
         t2 = *(temp + 2 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 2 * ldG) += *(Xjk + 0 * ldX) * t2;
         *(Gjk + 0 * ldG) += *(Xik + 2 * ldX) * t2;
         t3 = *(temp + 3 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 3 * ldG) += *(Xjk + 0 * ldX) * t3;
         *(Gjk + 0 * ldG) += *(Xik + 3 * ldX) * t3;
         t4 = *(temp + 4 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 4 * ldG) += *(Xjk + 0 * ldX) * t4;
         *(Gjk + 0 * ldG) += *(Xik + 4 * ldX) * t4;
         t5 = *(temp + 5 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 5 * ldG) += *(Xjk + 0 * ldX) * t5;
         *(Gjk + 0 * ldG) += *(Xik + 5 * ldX) * t5;
         t6 = *(temp + 6 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 6 * ldG) += *(Xjk + 0 * ldX) * t6;
         *(Gjk + 0 * ldG) += *(Xik + 6 * ldX) * t6;
         t7 = *(temp + 7 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 7 * ldG) += *(Xjk + 0 * ldX) * t7;
         *(Gjk + 0 * ldG) += *(Xik + 7 * ldX) * t7;
         t8 = *(temp + 8 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 8 * ldG) += *(Xjk + 0 * ldX) * t8;
         *(Gjk + 0 * ldG) += *(Xik + 8 * ldX) * t8;
         t9 = *(temp + 9 * NPTS_LOCAL + p_inner) * const_value;
         *(Gik + 9 * ldG) += *(Xjk + 0 * ldX) * t9;
         *(Gjk + 0 * ldG) += *(Xik + 9 * ldX) * t9;
      }
   }
}
