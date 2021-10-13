#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

void integral_2_1(int npts,
                  shell shellA,
                  shell shellB,
                  point *_points,
                  double *Xi,
                  double *Xj,
                  int ldX,
                  double *Gi,
                  double *Gj,
                  int ldG, 
                  double *weights) {
   double temp[16];

   for(int point_idx = 0; point_idx < npts; ++point_idx) {
      point C = *(_points + point_idx);

      double xA = shellA.origin.x;
      double yA = shellA.origin.y;
      double zA = shellA.origin.z;

      double xB = shellB.origin.x;
      double yB = shellB.origin.y;
      double zB = shellB.origin.z;

      double X_AB = (xA - xB);
      double Y_AB = (yA - yB);
      double Z_AB = (zA - zB);

      double beta_in = 0.0;
      for(int i = 0; i < shellA.m; ++i) {
         for(int j = 0; j < shellB.m; ++j) {
            double aA = shellA.coeff[i].alpha;
            double cA = shellA.coeff[i].coeff;

            double aB = shellB.coeff[j].alpha;
            double cB = shellB.coeff[j].coeff;

            double RHO = aA + aB;
            double RHO_INV = 1.0 / RHO;

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            double xP = (aA * xA + aB * xB) * RHO_INV;
            double yP = (aA * yA + aB * yB) * RHO_INV;
            double zP = (aA * zA + aB * zB) * RHO_INV;

            double X_PB = (xP - xB);
            double Y_PB = (yP - yB);
            double Z_PB = (zP - zB);

            double X_PC = (xP - xC);
            double Y_PC = (yP - yC);
            double Z_PC = (zP - zC);

            double t00, t01, t02, t03, t10, t11, t12, t20, t21, t30;

            double eval = cA * cB * 2 * PI * RHO_INV * exp(-1.0 * (X_AB * X_AB + Y_AB * Y_AB + Z_AB * Z_AB) * aA * aB * RHO_INV);
            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

#ifdef BOYS_REFERENCE
            t00 = eval * boys_reference(0, tval);
            t01 = eval * boys_reference(1, tval);
            t02 = eval * boys_reference(2, tval);
            t03 = eval * boys_reference(3, tval);
#elif BOYS_ASYMP
            t00 = eval * boys_asymp(0, tval);
            t01 = eval * boys_asymp(1, tval);
            t02 = eval * boys_asymp(2, tval);
            t03 = eval * boys_asymp(3, tval);
#else
            #error "TYPE NOT DEFINED!"
#endif

            t10 = X_PB * t00 - X_PC * t01;
            t11 = X_PB * t01 - X_PC * t02;
            t12 = X_PB * t02 - X_PC * t03;
            t20 = X_PB * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = X_PB * t11 - X_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            *(temp + 0) = beta_in * (*(temp + 0)) + t20;

            t30 = X_PB * t20 - X_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            *(temp + 6) = beta_in * (*(temp + 6)) + t30;

            t30 = Y_PB * t20 - Y_PC * t21;
            *(temp + 7) = beta_in * (*(temp + 7)) + t30;

            t30 = Z_PB * t20 - Z_PC * t21;
            *(temp + 8) = beta_in * (*(temp + 8)) + t30;

            t20 = Y_PB * t10 - Y_PC * t11;
            t21 = Y_PB * t11 - Y_PC * t12;
            *(temp + 1) = beta_in * (*(temp + 1)) + t20;

            t30 = Y_PB * t20 - Y_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            *(temp + 9) = beta_in * (*(temp + 9)) + t30;

            t30 = Z_PB * t20 - Z_PC * t21;
            *(temp + 10) = beta_in * (*(temp + 10)) + t30;

            t20 = Z_PB * t10 - Z_PC * t11;
            t21 = Z_PB * t11 - Z_PC * t12;
            *(temp + 2) = beta_in * (*(temp + 2)) + t20;

            t30 = Z_PB * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            *(temp + 11) = beta_in * (*(temp + 11)) + t30;

            t10 = Y_PB * t00 - Y_PC * t01;
            t11 = Y_PB * t01 - Y_PC * t02;
            t12 = Y_PB * t02 - Y_PC * t03;
            t20 = Y_PB * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Y_PB * t11 - Y_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            *(temp + 3) = beta_in * (*(temp + 3)) + t20;

            t30 = Y_PB * t20 - Y_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            *(temp + 12) = beta_in * (*(temp + 12)) + t30;

            t30 = Z_PB * t20 - Z_PC * t21;
            *(temp + 13) = beta_in * (*(temp + 13)) + t30;

            t20 = Z_PB * t10 - Z_PC * t11;
            t21 = Z_PB * t11 - Z_PC * t12;
            *(temp + 4) = beta_in * (*(temp + 4)) + t20;

            t30 = Z_PB * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            *(temp + 14) = beta_in * (*(temp + 14)) + t30;

            t10 = Z_PB * t00 - Z_PC * t01;
            t11 = Z_PB * t01 - Z_PC * t02;
            t12 = Z_PB * t02 - Z_PC * t03;
            t20 = Z_PB * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Z_PB * t11 - Z_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            *(temp + 5) = beta_in * (*(temp + 5)) + t20;

            t30 = Z_PB * t20 - Z_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            *(temp + 15) = beta_in * (*(temp + 15)) + t30;

            beta_in = 1.0;
         }
      }

      double *Xik = (Xi + point_idx * ldX);
      double *Xjk = (Xj + point_idx * ldX);
      double *Gik = (Gi + point_idx * ldG);
      double *Gjk = (Gj + point_idx * ldG);

      double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k, rcp_i, rcp_j, rcp_k;
      double t0, t1, t2, t3, t4, t5;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 6) * const_value;
      *(Gjk + 0) += *(Xik + 0) * t0;
      *(Gik + 0) += *(Xjk + 0) * t0;
      t1 = *(temp + 7) * const_value;
      *(Gjk + 1) += *(Xik + 0) * t1;
      *(Gik + 0) += *(Xjk + 1) * t1;
      t2 = *(temp + 8) * const_value;
      *(Gjk + 2) += *(Xik + 0) * t2;
      *(Gik + 0) += *(Xjk + 2) * t2;
      t3 = *(temp + 9) * const_value;
      *(Gjk + 3) += *(Xik + 0) * t3;
      *(Gik + 0) += *(Xjk + 3) * t3;
      t4 = *(temp + 10) * const_value;
      *(Gjk + 4) += *(Xik + 0) * t4;
      *(Gik + 0) += *(Xjk + 4) * t4;
      t5 = *(temp + 11) * const_value;
      *(Gjk + 5) += *(Xik + 0) * t5;
      *(Gik + 0) += *(Xjk + 5) * t5;

      X_ABp *= (-1.0) * X_AB; rcp_i = 1.0 / (1.0 * 1); comb_m_i = (comb_m_i * 1) * rcp_i;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 0) * t0;
      *(Gik + 0) += *(Xjk + 0) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 0) * t1;
      *(Gik + 0) += *(Xjk + 1) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 0) * t2;
      *(Gik + 0) += *(Xjk + 2) * t2;
      t3 = *(temp + 3) * const_value;
      *(Gjk + 3) += *(Xik + 0) * t3;
      *(Gik + 0) += *(Xjk + 3) * t3;
      t4 = *(temp + 4) * const_value;
      *(Gjk + 4) += *(Xik + 0) * t4;
      *(Gik + 0) += *(Xjk + 4) * t4;
      t5 = *(temp + 5) * const_value;
      *(Gjk + 5) += *(Xik + 0) * t5;
      *(Gik + 0) += *(Xjk + 5) * t5;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 7) * const_value;
      *(Gjk + 0) += *(Xik + 1) * t0;
      *(Gik + 1) += *(Xjk + 0) * t0;
      t1 = *(temp + 9) * const_value;
      *(Gjk + 1) += *(Xik + 1) * t1;
      *(Gik + 1) += *(Xjk + 1) * t1;
      t2 = *(temp + 10) * const_value;
      *(Gjk + 2) += *(Xik + 1) * t2;
      *(Gik + 1) += *(Xjk + 2) * t2;
      t3 = *(temp + 12) * const_value;
      *(Gjk + 3) += *(Xik + 1) * t3;
      *(Gik + 1) += *(Xjk + 3) * t3;
      t4 = *(temp + 13) * const_value;
      *(Gjk + 4) += *(Xik + 1) * t4;
      *(Gik + 1) += *(Xjk + 4) * t4;
      t5 = *(temp + 14) * const_value;
      *(Gjk + 5) += *(Xik + 1) * t5;
      *(Gik + 1) += *(Xjk + 5) * t5;

      Y_ABp *= (-1.0) * Y_AB; rcp_j = 1.0 / (1.0 * 1); comb_n_j = (comb_n_j * 1) * rcp_j;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 1) * t0;
      *(Gik + 1) += *(Xjk + 0) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 1) * t1;
      *(Gik + 1) += *(Xjk + 1) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 1) * t2;
      *(Gik + 1) += *(Xjk + 2) * t2;
      t3 = *(temp + 3) * const_value;
      *(Gjk + 3) += *(Xik + 1) * t3;
      *(Gik + 1) += *(Xjk + 3) * t3;
      t4 = *(temp + 4) * const_value;
      *(Gjk + 4) += *(Xik + 1) * t4;
      *(Gik + 1) += *(Xjk + 4) * t4;
      t5 = *(temp + 5) * const_value;
      *(Gjk + 5) += *(Xik + 1) * t5;
      *(Gik + 1) += *(Xjk + 5) * t5;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 8) * const_value;
      *(Gjk + 0) += *(Xik + 2) * t0;
      *(Gik + 2) += *(Xjk + 0) * t0;
      t1 = *(temp + 10) * const_value;
      *(Gjk + 1) += *(Xik + 2) * t1;
      *(Gik + 2) += *(Xjk + 1) * t1;
      t2 = *(temp + 11) * const_value;
      *(Gjk + 2) += *(Xik + 2) * t2;
      *(Gik + 2) += *(Xjk + 2) * t2;
      t3 = *(temp + 13) * const_value;
      *(Gjk + 3) += *(Xik + 2) * t3;
      *(Gik + 2) += *(Xjk + 3) * t3;
      t4 = *(temp + 14) * const_value;
      *(Gjk + 4) += *(Xik + 2) * t4;
      *(Gik + 2) += *(Xjk + 4) * t4;
      t5 = *(temp + 15) * const_value;
      *(Gjk + 5) += *(Xik + 2) * t5;
      *(Gik + 2) += *(Xjk + 5) * t5;

      Z_ABp *= (-1.0) * Z_AB; rcp_k = 1.0 / (1.0 * 1); comb_p_k = (comb_p_k * 1) * rcp_k;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 2) * t0;
      *(Gik + 2) += *(Xjk + 0) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 2) * t1;
      *(Gik + 2) += *(Xjk + 1) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 2) * t2;
      *(Gik + 2) += *(Xjk + 2) * t2;
      t3 = *(temp + 3) * const_value;
      *(Gjk + 3) += *(Xik + 2) * t3;
      *(Gik + 2) += *(Xjk + 3) * t3;
      t4 = *(temp + 4) * const_value;
      *(Gjk + 4) += *(Xik + 2) * t4;
      *(Gik + 2) += *(Xjk + 4) * t4;
      t5 = *(temp + 5) * const_value;
      *(Gjk + 5) += *(Xik + 2) * t5;
      *(Gik + 2) += *(Xjk + 5) * t5;

      *(Gjk + 0) *= *(weights + point_idx);
      *(Gjk + 1) *= *(weights + point_idx);
      *(Gjk + 2) *= *(weights + point_idx);
      *(Gjk + 3) *= *(weights + point_idx);
      *(Gjk + 4) *= *(weights + point_idx);
      *(Gjk + 5) *= *(weights + point_idx);

      *(Gik + 0) *= *(weights + point_idx);
      *(Gik + 1) *= *(weights + point_idx);
      *(Gik + 2) *= *(weights + point_idx);
   }
}
