#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

void integral_2_2(int npts,
                  shell shellA,
                  shell shellB,
                  point *_points,
                  double *Xi,
                  int ldX,
                  double *Gj,
                  int ldG, 
                  double *weights) {
   double temp[31];

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

            double t00, t01, t02, t03, t04, t10, t11, t12, t13, t20, t21, t22, t30, t31, t40;

            double eval = cA * cB * 2 * PI * RHO_INV * exp(-1.0 * (X_AB * X_AB + Y_AB * Y_AB + Z_AB * Z_AB) * aA * aB * RHO_INV);
            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

#ifdef BOYS_REFERENCE
            t00 = eval * boys_reference(0, tval);
            t01 = eval * boys_reference(1, tval);
            t02 = eval * boys_reference(2, tval);
            t03 = eval * boys_reference(3, tval);
            t04 = eval * boys_reference(4, tval);
#elif BOYS_ASYMP
            t00 = eval * boys_asymp(0, tval);
            t01 = eval * boys_asymp(1, tval);
            t02 = eval * boys_asymp(2, tval);
            t03 = eval * boys_asymp(3, tval);
            t04 = eval * boys_asymp(4, tval);
#else
            #error "TYPE NOT DEFINED!"
#endif

            t10 = X_PB * t00 - X_PC * t01;
            t11 = X_PB * t01 - X_PC * t02;
            t12 = X_PB * t02 - X_PC * t03;
            t13 = X_PB * t03 - X_PC * t04;
            t20 = X_PB * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = X_PB * t11 - X_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = X_PB * t12 - X_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            *(temp + 0) = beta_in * (*(temp + 0)) + t20;

            t30 = X_PB * t20 - X_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = X_PB * t21 - X_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            *(temp + 6) = beta_in * (*(temp + 6)) + t30;

            t40 = X_PB * t30 - X_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 16) = beta_in * (*(temp + 16)) + t40;

            t40 = Y_PB * t30 - Y_PC * t31;
            *(temp + 17) = beta_in * (*(temp + 17)) + t40;

            t40 = Z_PB * t30 - Z_PC * t31;
            *(temp + 18) = beta_in * (*(temp + 18)) + t40;

            t30 = Y_PB * t20 - Y_PC * t21;
            t31 = Y_PB * t21 - Y_PC * t22;
            *(temp + 7) = beta_in * (*(temp + 7)) + t30;

            t40 = Y_PB * t30 - Y_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 19) = beta_in * (*(temp + 19)) + t40;

            t40 = Z_PB * t30 - Z_PC * t31;
            *(temp + 20) = beta_in * (*(temp + 20)) + t40;

            t30 = Z_PB * t20 - Z_PC * t21;
            t31 = Z_PB * t21 - Z_PC * t22;
            *(temp + 8) = beta_in * (*(temp + 8)) + t30;

            t40 = Z_PB * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 21) = beta_in * (*(temp + 21)) + t40;

            t20 = Y_PB * t10 - Y_PC * t11;
            t21 = Y_PB * t11 - Y_PC * t12;
            t22 = Y_PB * t12 - Y_PC * t13;
            *(temp + 1) = beta_in * (*(temp + 1)) + t20;

            t30 = Y_PB * t20 - Y_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Y_PB * t21 - Y_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            *(temp + 9) = beta_in * (*(temp + 9)) + t30;

            t40 = Y_PB * t30 - Y_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 22) = beta_in * (*(temp + 22)) + t40;

            t40 = Z_PB * t30 - Z_PC * t31;
            *(temp + 23) = beta_in * (*(temp + 23)) + t40;

            t30 = Z_PB * t20 - Z_PC * t21;
            t31 = Z_PB * t21 - Z_PC * t22;
            *(temp + 10) = beta_in * (*(temp + 10)) + t30;

            t40 = Z_PB * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 24) = beta_in * (*(temp + 24)) + t40;

            t20 = Z_PB * t10 - Z_PC * t11;
            t21 = Z_PB * t11 - Z_PC * t12;
            t22 = Z_PB * t12 - Z_PC * t13;
            *(temp + 2) = beta_in * (*(temp + 2)) + t20;

            t30 = Z_PB * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PB * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            *(temp + 11) = beta_in * (*(temp + 11)) + t30;

            t40 = Z_PB * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 25) = beta_in * (*(temp + 25)) + t40;

            t10 = Y_PB * t00 - Y_PC * t01;
            t11 = Y_PB * t01 - Y_PC * t02;
            t12 = Y_PB * t02 - Y_PC * t03;
            t13 = Y_PB * t03 - Y_PC * t04;
            t20 = Y_PB * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Y_PB * t11 - Y_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Y_PB * t12 - Y_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            *(temp + 3) = beta_in * (*(temp + 3)) + t20;

            t30 = Y_PB * t20 - Y_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Y_PB * t21 - Y_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            *(temp + 12) = beta_in * (*(temp + 12)) + t30;

            t40 = Y_PB * t30 - Y_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 26) = beta_in * (*(temp + 26)) + t40;

            t40 = Z_PB * t30 - Z_PC * t31;
            *(temp + 27) = beta_in * (*(temp + 27)) + t40;

            t30 = Z_PB * t20 - Z_PC * t21;
            t31 = Z_PB * t21 - Z_PC * t22;
            *(temp + 13) = beta_in * (*(temp + 13)) + t30;

            t40 = Z_PB * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 28) = beta_in * (*(temp + 28)) + t40;

            t20 = Z_PB * t10 - Z_PC * t11;
            t21 = Z_PB * t11 - Z_PC * t12;
            t22 = Z_PB * t12 - Z_PC * t13;
            *(temp + 4) = beta_in * (*(temp + 4)) + t20;

            t30 = Z_PB * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PB * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            *(temp + 14) = beta_in * (*(temp + 14)) + t30;

            t40 = Z_PB * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 29) = beta_in * (*(temp + 29)) + t40;

            t10 = Z_PB * t00 - Z_PC * t01;
            t11 = Z_PB * t01 - Z_PC * t02;
            t12 = Z_PB * t02 - Z_PC * t03;
            t13 = Z_PB * t03 - Z_PC * t04;
            t20 = Z_PB * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Z_PB * t11 - Z_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Z_PB * t12 - Z_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            *(temp + 5) = beta_in * (*(temp + 5)) + t20;

            t30 = Z_PB * t20 - Z_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Z_PB * t21 - Z_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            *(temp + 15) = beta_in * (*(temp + 15)) + t30;

            t40 = Z_PB * t30 - Z_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 30) = beta_in * (*(temp + 30)) + t40;

            beta_in = 1.0;
         }
      }

      double *Xik = (Xi + point_idx * ldX);
      double *Gjk = (Gj + point_idx * ldG);

      double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k, rcp_i, rcp_j, rcp_k;
      double t0, t1, t2, t3, t4, t5;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 16) * const_value;
      *(Gjk + 0) += *(Xik + 0) * t0;
      t1 = *(temp + 17) * const_value;
      *(Gjk + 1) += *(Xik + 0) * t1;
      t2 = *(temp + 18) * const_value;
      *(Gjk + 2) += *(Xik + 0) * t2;
      t3 = *(temp + 19) * const_value;
      *(Gjk + 3) += *(Xik + 0) * t3;
      t4 = *(temp + 20) * const_value;
      *(Gjk + 4) += *(Xik + 0) * t4;
      t5 = *(temp + 21) * const_value;
      *(Gjk + 5) += *(Xik + 0) * t5;

      X_ABp *= (-1.0) * X_AB; rcp_i = 1.0 / (1.0 * 1); comb_m_i = (comb_m_i * 2) * rcp_i;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 6) * const_value;
      *(Gjk + 0) += *(Xik + 0) * t0;
      t1 = *(temp + 7) * const_value;
      *(Gjk + 1) += *(Xik + 0) * t1;
      t2 = *(temp + 8) * const_value;
      *(Gjk + 2) += *(Xik + 0) * t2;
      t3 = *(temp + 9) * const_value;
      *(Gjk + 3) += *(Xik + 0) * t3;
      t4 = *(temp + 10) * const_value;
      *(Gjk + 4) += *(Xik + 0) * t4;
      t5 = *(temp + 11) * const_value;
      *(Gjk + 5) += *(Xik + 0) * t5;

      X_ABp *= (-1.0) * X_AB; rcp_i = 1.0 / (1.0 * 2); comb_m_i = (comb_m_i * 1) * rcp_i;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 0) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 0) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 0) * t2;
      t3 = *(temp + 3) * const_value;
      *(Gjk + 3) += *(Xik + 0) * t3;
      t4 = *(temp + 4) * const_value;
      *(Gjk + 4) += *(Xik + 0) * t4;
      t5 = *(temp + 5) * const_value;
      *(Gjk + 5) += *(Xik + 0) * t5;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 17) * const_value;
      *(Gjk + 0) += *(Xik + 1) * t0;
      t1 = *(temp + 19) * const_value;
      *(Gjk + 1) += *(Xik + 1) * t1;
      t2 = *(temp + 20) * const_value;
      *(Gjk + 2) += *(Xik + 1) * t2;
      t3 = *(temp + 22) * const_value;
      *(Gjk + 3) += *(Xik + 1) * t3;
      t4 = *(temp + 23) * const_value;
      *(Gjk + 4) += *(Xik + 1) * t4;
      t5 = *(temp + 24) * const_value;
      *(Gjk + 5) += *(Xik + 1) * t5;

      Y_ABp *= (-1.0) * Y_AB; rcp_j = 1.0 / (1.0 * 1); comb_n_j = (comb_n_j * 1) * rcp_j;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 6) * const_value;
      *(Gjk + 0) += *(Xik + 1) * t0;
      t1 = *(temp + 7) * const_value;
      *(Gjk + 1) += *(Xik + 1) * t1;
      t2 = *(temp + 8) * const_value;
      *(Gjk + 2) += *(Xik + 1) * t2;
      t3 = *(temp + 9) * const_value;
      *(Gjk + 3) += *(Xik + 1) * t3;
      t4 = *(temp + 10) * const_value;
      *(Gjk + 4) += *(Xik + 1) * t4;
      t5 = *(temp + 11) * const_value;
      *(Gjk + 5) += *(Xik + 1) * t5;

      X_ABp *= (-1.0) * X_AB; rcp_i = 1.0 / (1.0 * 1); comb_m_i = (comb_m_i * 1) * rcp_i;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 7) * const_value;
      *(Gjk + 0) += *(Xik + 1) * t0;
      t1 = *(temp + 9) * const_value;
      *(Gjk + 1) += *(Xik + 1) * t1;
      t2 = *(temp + 10) * const_value;
      *(Gjk + 2) += *(Xik + 1) * t2;
      t3 = *(temp + 12) * const_value;
      *(Gjk + 3) += *(Xik + 1) * t3;
      t4 = *(temp + 13) * const_value;
      *(Gjk + 4) += *(Xik + 1) * t4;
      t5 = *(temp + 14) * const_value;
      *(Gjk + 5) += *(Xik + 1) * t5;

      Y_ABp *= (-1.0) * Y_AB; rcp_j = 1.0 / (1.0 * 1); comb_n_j = (comb_n_j * 1) * rcp_j;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 1) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 1) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 1) * t2;
      t3 = *(temp + 3) * const_value;
      *(Gjk + 3) += *(Xik + 1) * t3;
      t4 = *(temp + 4) * const_value;
      *(Gjk + 4) += *(Xik + 1) * t4;
      t5 = *(temp + 5) * const_value;
      *(Gjk + 5) += *(Xik + 1) * t5;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 18) * const_value;
      *(Gjk + 0) += *(Xik + 2) * t0;
      t1 = *(temp + 20) * const_value;
      *(Gjk + 1) += *(Xik + 2) * t1;
      t2 = *(temp + 21) * const_value;
      *(Gjk + 2) += *(Xik + 2) * t2;
      t3 = *(temp + 23) * const_value;
      *(Gjk + 3) += *(Xik + 2) * t3;
      t4 = *(temp + 24) * const_value;
      *(Gjk + 4) += *(Xik + 2) * t4;
      t5 = *(temp + 25) * const_value;
      *(Gjk + 5) += *(Xik + 2) * t5;

      Z_ABp *= (-1.0) * Z_AB; rcp_k = 1.0 / (1.0 * 1); comb_p_k = (comb_p_k * 1) * rcp_k;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 6) * const_value;
      *(Gjk + 0) += *(Xik + 2) * t0;
      t1 = *(temp + 7) * const_value;
      *(Gjk + 1) += *(Xik + 2) * t1;
      t2 = *(temp + 8) * const_value;
      *(Gjk + 2) += *(Xik + 2) * t2;
      t3 = *(temp + 9) * const_value;
      *(Gjk + 3) += *(Xik + 2) * t3;
      t4 = *(temp + 10) * const_value;
      *(Gjk + 4) += *(Xik + 2) * t4;
      t5 = *(temp + 11) * const_value;
      *(Gjk + 5) += *(Xik + 2) * t5;

      X_ABp *= (-1.0) * X_AB; rcp_i = 1.0 / (1.0 * 1); comb_m_i = (comb_m_i * 1) * rcp_i;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 8) * const_value;
      *(Gjk + 0) += *(Xik + 2) * t0;
      t1 = *(temp + 10) * const_value;
      *(Gjk + 1) += *(Xik + 2) * t1;
      t2 = *(temp + 11) * const_value;
      *(Gjk + 2) += *(Xik + 2) * t2;
      t3 = *(temp + 13) * const_value;
      *(Gjk + 3) += *(Xik + 2) * t3;
      t4 = *(temp + 14) * const_value;
      *(Gjk + 4) += *(Xik + 2) * t4;
      t5 = *(temp + 15) * const_value;
      *(Gjk + 5) += *(Xik + 2) * t5;

      Z_ABp *= (-1.0) * Z_AB; rcp_k = 1.0 / (1.0 * 1); comb_p_k = (comb_p_k * 1) * rcp_k;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 2) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 2) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 2) * t2;
      t3 = *(temp + 3) * const_value;
      *(Gjk + 3) += *(Xik + 2) * t3;
      t4 = *(temp + 4) * const_value;
      *(Gjk + 4) += *(Xik + 2) * t4;
      t5 = *(temp + 5) * const_value;
      *(Gjk + 5) += *(Xik + 2) * t5;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 19) * const_value;
      *(Gjk + 0) += *(Xik + 3) * t0;
      t1 = *(temp + 22) * const_value;
      *(Gjk + 1) += *(Xik + 3) * t1;
      t2 = *(temp + 23) * const_value;
      *(Gjk + 2) += *(Xik + 3) * t2;
      t3 = *(temp + 26) * const_value;
      *(Gjk + 3) += *(Xik + 3) * t3;
      t4 = *(temp + 27) * const_value;
      *(Gjk + 4) += *(Xik + 3) * t4;
      t5 = *(temp + 28) * const_value;
      *(Gjk + 5) += *(Xik + 3) * t5;

      Y_ABp *= (-1.0) * Y_AB; rcp_j = 1.0 / (1.0 * 1); comb_n_j = (comb_n_j * 2) * rcp_j;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 7) * const_value;
      *(Gjk + 0) += *(Xik + 3) * t0;
      t1 = *(temp + 9) * const_value;
      *(Gjk + 1) += *(Xik + 3) * t1;
      t2 = *(temp + 10) * const_value;
      *(Gjk + 2) += *(Xik + 3) * t2;
      t3 = *(temp + 12) * const_value;
      *(Gjk + 3) += *(Xik + 3) * t3;
      t4 = *(temp + 13) * const_value;
      *(Gjk + 4) += *(Xik + 3) * t4;
      t5 = *(temp + 14) * const_value;
      *(Gjk + 5) += *(Xik + 3) * t5;

      Y_ABp *= (-1.0) * Y_AB; rcp_j = 1.0 / (1.0 * 2); comb_n_j = (comb_n_j * 1) * rcp_j;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 3) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 3) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 3) * t2;
      t3 = *(temp + 3) * const_value;
      *(Gjk + 3) += *(Xik + 3) * t3;
      t4 = *(temp + 4) * const_value;
      *(Gjk + 4) += *(Xik + 3) * t4;
      t5 = *(temp + 5) * const_value;
      *(Gjk + 5) += *(Xik + 3) * t5;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 20) * const_value;
      *(Gjk + 0) += *(Xik + 4) * t0;
      t1 = *(temp + 23) * const_value;
      *(Gjk + 1) += *(Xik + 4) * t1;
      t2 = *(temp + 24) * const_value;
      *(Gjk + 2) += *(Xik + 4) * t2;
      t3 = *(temp + 27) * const_value;
      *(Gjk + 3) += *(Xik + 4) * t3;
      t4 = *(temp + 28) * const_value;
      *(Gjk + 4) += *(Xik + 4) * t4;
      t5 = *(temp + 29) * const_value;
      *(Gjk + 5) += *(Xik + 4) * t5;

      Z_ABp *= (-1.0) * Z_AB; rcp_k = 1.0 / (1.0 * 1); comb_p_k = (comb_p_k * 1) * rcp_k;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 7) * const_value;
      *(Gjk + 0) += *(Xik + 4) * t0;
      t1 = *(temp + 9) * const_value;
      *(Gjk + 1) += *(Xik + 4) * t1;
      t2 = *(temp + 10) * const_value;
      *(Gjk + 2) += *(Xik + 4) * t2;
      t3 = *(temp + 12) * const_value;
      *(Gjk + 3) += *(Xik + 4) * t3;
      t4 = *(temp + 13) * const_value;
      *(Gjk + 4) += *(Xik + 4) * t4;
      t5 = *(temp + 14) * const_value;
      *(Gjk + 5) += *(Xik + 4) * t5;

      Y_ABp *= (-1.0) * Y_AB; rcp_j = 1.0 / (1.0 * 1); comb_n_j = (comb_n_j * 1) * rcp_j;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 8) * const_value;
      *(Gjk + 0) += *(Xik + 4) * t0;
      t1 = *(temp + 10) * const_value;
      *(Gjk + 1) += *(Xik + 4) * t1;
      t2 = *(temp + 11) * const_value;
      *(Gjk + 2) += *(Xik + 4) * t2;
      t3 = *(temp + 13) * const_value;
      *(Gjk + 3) += *(Xik + 4) * t3;
      t4 = *(temp + 14) * const_value;
      *(Gjk + 4) += *(Xik + 4) * t4;
      t5 = *(temp + 15) * const_value;
      *(Gjk + 5) += *(Xik + 4) * t5;

      Z_ABp *= (-1.0) * Z_AB; rcp_k = 1.0 / (1.0 * 1); comb_p_k = (comb_p_k * 1) * rcp_k;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 4) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 4) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 4) * t2;
      t3 = *(temp + 3) * const_value;
      *(Gjk + 3) += *(Xik + 4) * t3;
      t4 = *(temp + 4) * const_value;
      *(Gjk + 4) += *(Xik + 4) * t4;
      t5 = *(temp + 5) * const_value;
      *(Gjk + 5) += *(Xik + 4) * t5;

      X_ABp = 1.0; comb_m_i = 1.0;
      Y_ABp = 1.0; comb_n_j = 1.0;
      Z_ABp = 1.0; comb_p_k = 1.0;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 21) * const_value;
      *(Gjk + 0) += *(Xik + 5) * t0;
      t1 = *(temp + 24) * const_value;
      *(Gjk + 1) += *(Xik + 5) * t1;
      t2 = *(temp + 25) * const_value;
      *(Gjk + 2) += *(Xik + 5) * t2;
      t3 = *(temp + 28) * const_value;
      *(Gjk + 3) += *(Xik + 5) * t3;
      t4 = *(temp + 29) * const_value;
      *(Gjk + 4) += *(Xik + 5) * t4;
      t5 = *(temp + 30) * const_value;
      *(Gjk + 5) += *(Xik + 5) * t5;

      Z_ABp *= (-1.0) * Z_AB; rcp_k = 1.0 / (1.0 * 1); comb_p_k = (comb_p_k * 2) * rcp_k;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 8) * const_value;
      *(Gjk + 0) += *(Xik + 5) * t0;
      t1 = *(temp + 10) * const_value;
      *(Gjk + 1) += *(Xik + 5) * t1;
      t2 = *(temp + 11) * const_value;
      *(Gjk + 2) += *(Xik + 5) * t2;
      t3 = *(temp + 13) * const_value;
      *(Gjk + 3) += *(Xik + 5) * t3;
      t4 = *(temp + 14) * const_value;
      *(Gjk + 4) += *(Xik + 5) * t4;
      t5 = *(temp + 15) * const_value;
      *(Gjk + 5) += *(Xik + 5) * t5;

      Z_ABp *= (-1.0) * Z_AB; rcp_k = 1.0 / (1.0 * 2); comb_p_k = (comb_p_k * 1) * rcp_k;
      const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

      t0 = *(temp + 0) * const_value;
      *(Gjk + 0) += *(Xik + 5) * t0;
      t1 = *(temp + 1) * const_value;
      *(Gjk + 1) += *(Xik + 5) * t1;
      t2 = *(temp + 2) * const_value;
      *(Gjk + 2) += *(Xik + 5) * t2;
      t3 = *(temp + 3) * const_value;
      *(Gjk + 3) += *(Xik + 5) * t3;
      t4 = *(temp + 4) * const_value;
      *(Gjk + 4) += *(Xik + 5) * t4;
      t5 = *(temp + 5) * const_value;
      *(Gjk + 5) += *(Xik + 5) * t5;

      *(Gjk + 0) *= *(weights + point_idx);
      *(Gjk + 1) *= *(weights + point_idx);
      *(Gjk + 2) *= *(weights + point_idx);
      *(Gjk + 3) *= *(weights + point_idx);
      *(Gjk + 4) *= *(weights + point_idx);
      *(Gjk + 5) *= *(weights + point_idx);
   }
}
