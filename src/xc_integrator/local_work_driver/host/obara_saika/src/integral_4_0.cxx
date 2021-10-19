#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

void integral_4_0(size_t npts,
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
   double temp[15];

   for(int i = 0; i < 15; ++i) {
      temp[i] = 0.0;
   }

   for(size_t point_idx = 0; point_idx < npts; ++point_idx) {
      point C = *(_points + point_idx);

      double xA = shpair.rA.x;
      double yA = shpair.rA.y;
      double zA = shpair.rA.z;

      double xB = shpair.rB.x;
      double yB = shpair.rB.y;
      double zB = shpair.rB.z;

      double X_AB = shpair.rAB.x;
      double Y_AB = shpair.rAB.y;
      double Z_AB = shpair.rAB.z;

      double beta_in = 0.0;
      for(int ij = 0; ij < shpair.nprim_pair; ++ij ) {
            double RHO = shpair.prim_pairs[ij].gamma;
            double RHO_INV = 1.0 / RHO;

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            double xP = shpair.prim_pairs[ij].P.x;
            double yP = shpair.prim_pairs[ij].P.y;
            double zP = shpair.prim_pairs[ij].P.z;

            double X_PA = shpair.prim_pairs[ij].PA.x;
            double Y_PA = shpair.prim_pairs[ij].PA.y;
            double Z_PA = shpair.prim_pairs[ij].PA.z;

            double X_PC = (xP - xC);
            double Y_PC = (yP - yC);
            double Z_PC = (zP - zC);

            double t00, t01, t02, t03, t04, t10, t11, t12, t13, t20, t21, t22, t30, t31, t40;

            double eval = shpair.prim_pairs[ij].coeff_prod * shpair.prim_pairs[ij].K;
            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

            t00 = eval * boys_function(0, tval);
            t01 = eval * boys_function(1, tval);
            t02 = eval * boys_function(2, tval);
            t03 = eval * boys_function(3, tval);
            t04 = eval * boys_function(4, tval);
            t10 = X_PA * t00 - X_PC * t01;
            t11 = X_PA * t01 - X_PC * t02;
            t12 = X_PA * t02 - X_PC * t03;
            t13 = X_PA * t03 - X_PC * t04;
            t20 = X_PA * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = X_PA * t11 - X_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = X_PA * t12 - X_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t30 = X_PA * t20 - X_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = X_PA * t21 - X_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            t40 = X_PA * t30 - X_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 0) = beta_in * (*(temp + 0)) + t40;

            t40 = Y_PA * t30 - Y_PC * t31;
            *(temp + 1) = beta_in * (*(temp + 1)) + t40;

            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 2) = beta_in * (*(temp + 2)) + t40;

            t30 = Y_PA * t20 - Y_PC * t21;
            t31 = Y_PA * t21 - Y_PC * t22;
            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 3) = beta_in * (*(temp + 3)) + t40;

            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 4) = beta_in * (*(temp + 4)) + t40;

            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 5) = beta_in * (*(temp + 5)) + t40;

            t20 = Y_PA * t10 - Y_PC * t11;
            t21 = Y_PA * t11 - Y_PC * t12;
            t22 = Y_PA * t12 - Y_PC * t13;
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Y_PA * t21 - Y_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 6) = beta_in * (*(temp + 6)) + t40;

            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 7) = beta_in * (*(temp + 7)) + t40;

            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 8) = beta_in * (*(temp + 8)) + t40;

            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t22 = Z_PA * t12 - Z_PC * t13;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 9) = beta_in * (*(temp + 9)) + t40;

            t10 = Y_PA * t00 - Y_PC * t01;
            t11 = Y_PA * t01 - Y_PC * t02;
            t12 = Y_PA * t02 - Y_PC * t03;
            t13 = Y_PA * t03 - Y_PC * t04;
            t20 = Y_PA * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Y_PA * t11 - Y_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Y_PA * t12 - Y_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Y_PA * t21 - Y_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 10) = beta_in * (*(temp + 10)) + t40;

            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 11) = beta_in * (*(temp + 11)) + t40;

            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 12) = beta_in * (*(temp + 12)) + t40;

            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t22 = Z_PA * t12 - Z_PC * t13;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 13) = beta_in * (*(temp + 13)) + t40;

            t10 = Z_PA * t00 - Z_PC * t01;
            t11 = Z_PA * t01 - Z_PC * t02;
            t12 = Z_PA * t02 - Z_PC * t03;
            t13 = Z_PA * t03 - Z_PC * t04;
            t20 = Z_PA * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Z_PA * t11 - Z_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Z_PA * t12 - Z_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 14) = beta_in * (*(temp + 14)) + t40;

            beta_in = 1.0;
      }

      double *Xik = (Xi + point_idx * stX);
      double *Xjk = (Xj + point_idx * stX);
      double *Gik = (Gi + point_idx * stG);
      double *Gjk = (Gj + point_idx * stG);

      double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k, rcp_i, rcp_j, rcp_k;
      double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value * (*(weights + point_idx));
      *(Gik + 0 * ldG) += *(Xjk + 0 * ldX) * t0;
      *(Gjk + 0 * ldG) += *(Xik + 0 * ldX) * t0;
      t1 = *(temp + 1) * const_value * (*(weights + point_idx));
      *(Gik + 1 * ldG) += *(Xjk + 0 * ldX) * t1;
      *(Gjk + 0 * ldG) += *(Xik + 1 * ldX) * t1;
      t2 = *(temp + 2) * const_value * (*(weights + point_idx));
      *(Gik + 2 * ldG) += *(Xjk + 0 * ldX) * t2;
      *(Gjk + 0 * ldG) += *(Xik + 2 * ldX) * t2;
      t3 = *(temp + 3) * const_value * (*(weights + point_idx));
      *(Gik + 3 * ldG) += *(Xjk + 0 * ldX) * t3;
      *(Gjk + 0 * ldG) += *(Xik + 3 * ldX) * t3;
      t4 = *(temp + 4) * const_value * (*(weights + point_idx));
      *(Gik + 4 * ldG) += *(Xjk + 0 * ldX) * t4;
      *(Gjk + 0 * ldG) += *(Xik + 4 * ldX) * t4;
      t5 = *(temp + 5) * const_value * (*(weights + point_idx));
      *(Gik + 5 * ldG) += *(Xjk + 0 * ldX) * t5;
      *(Gjk + 0 * ldG) += *(Xik + 5 * ldX) * t5;
      t6 = *(temp + 6) * const_value * (*(weights + point_idx));
      *(Gik + 6 * ldG) += *(Xjk + 0 * ldX) * t6;
      *(Gjk + 0 * ldG) += *(Xik + 6 * ldX) * t6;
      t7 = *(temp + 7) * const_value * (*(weights + point_idx));
      *(Gik + 7 * ldG) += *(Xjk + 0 * ldX) * t7;
      *(Gjk + 0 * ldG) += *(Xik + 7 * ldX) * t7;
      t8 = *(temp + 8) * const_value * (*(weights + point_idx));
      *(Gik + 8 * ldG) += *(Xjk + 0 * ldX) * t8;
      *(Gjk + 0 * ldG) += *(Xik + 8 * ldX) * t8;
      t9 = *(temp + 9) * const_value * (*(weights + point_idx));
      *(Gik + 9 * ldG) += *(Xjk + 0 * ldX) * t9;
      *(Gjk + 0 * ldG) += *(Xik + 9 * ldX) * t9;
      t10 = *(temp + 10) * const_value * (*(weights + point_idx));
      *(Gik + 10 * ldG) += *(Xjk + 0 * ldX) * t10;
      *(Gjk + 0 * ldG) += *(Xik + 10 * ldX) * t10;
      t11 = *(temp + 11) * const_value * (*(weights + point_idx));
      *(Gik + 11 * ldG) += *(Xjk + 0 * ldX) * t11;
      *(Gjk + 0 * ldG) += *(Xik + 11 * ldX) * t11;
      t12 = *(temp + 12) * const_value * (*(weights + point_idx));
      *(Gik + 12 * ldG) += *(Xjk + 0 * ldX) * t12;
      *(Gjk + 0 * ldG) += *(Xik + 12 * ldX) * t12;
      t13 = *(temp + 13) * const_value * (*(weights + point_idx));
      *(Gik + 13 * ldG) += *(Xjk + 0 * ldX) * t13;
      *(Gjk + 0 * ldG) += *(Xik + 13 * ldX) * t13;
      t14 = *(temp + 14) * const_value * (*(weights + point_idx));
      *(Gik + 14 * ldG) += *(Xjk + 0 * ldX) * t14;
      *(Gjk + 0 * ldG) += *(Xik + 14 * ldX) * t14;
   }
}
