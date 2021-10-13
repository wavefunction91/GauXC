#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

void integral_4_4(int npts,
                  shell shellA,
                  shell shellB,
                  point *_points,
                  double *Xi,
                  int ldX,
                  double *Gj,
                  int ldG, 
                  double *weights) {
   double temp[145];

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

            double t00, t01, t02, t03, t04, t05, t06, t07, t08, t10, t11, t12, t13, t14, t15, t16, t17, t20, t21, t22, t23, t24, t25, t26, t30, t31, t32, t33, t34, t35, t40, t41, t42, t43, t44, t50, t51, t52, t53, t60, t61, t62, t70, t71, t80;

            double eval = cA * cB * 2 * PI * RHO_INV * exp(-1.0 * (X_AB * X_AB + Y_AB * Y_AB + Z_AB * Z_AB) * aA * aB * RHO_INV);
            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

#ifdef BOYS_REFERENCE
            t00 = eval * boys_reference(0, tval);
            t01 = eval * boys_reference(1, tval);
            t02 = eval * boys_reference(2, tval);
            t03 = eval * boys_reference(3, tval);
            t04 = eval * boys_reference(4, tval);
            t05 = eval * boys_reference(5, tval);
            t06 = eval * boys_reference(6, tval);
            t07 = eval * boys_reference(7, tval);
            t08 = eval * boys_reference(8, tval);
#elif BOYS_ASYMP
            t00 = eval * boys_asymp(0, tval);
            t01 = eval * boys_asymp(1, tval);
            t02 = eval * boys_asymp(2, tval);
            t03 = eval * boys_asymp(3, tval);
            t04 = eval * boys_asymp(4, tval);
            t05 = eval * boys_asymp(5, tval);
            t06 = eval * boys_asymp(6, tval);
            t07 = eval * boys_asymp(7, tval);
            t08 = eval * boys_asymp(8, tval);
#else
            #error "TYPE NOT DEFINED!"
#endif

            t10 = X_PB * t00 - X_PC * t01;
            t11 = X_PB * t01 - X_PC * t02;
            t12 = X_PB * t02 - X_PC * t03;
            t13 = X_PB * t03 - X_PC * t04;
            t14 = X_PB * t04 - X_PC * t05;
            t15 = X_PB * t05 - X_PC * t06;
            t16 = X_PB * t06 - X_PC * t07;
            t17 = X_PB * t07 - X_PC * t08;
            t20 = X_PB * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = X_PB * t11 - X_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = X_PB * t12 - X_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t23 = X_PB * t13 - X_PC * t14 + 0.5 * RHO_INV * 1 * (t03 - t04);
            t24 = X_PB * t14 - X_PC * t15 + 0.5 * RHO_INV * 1 * (t04 - t05);
            t25 = X_PB * t15 - X_PC * t16 + 0.5 * RHO_INV * 1 * (t05 - t06);
            t26 = X_PB * t16 - X_PC * t17 + 0.5 * RHO_INV * 1 * (t06 - t07);
            t30 = X_PB * t20 - X_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = X_PB * t21 - X_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            t32 = X_PB * t22 - X_PC * t23 + 0.5 * RHO_INV * 2 * (t12 - t13);
            t33 = X_PB * t23 - X_PC * t24 + 0.5 * RHO_INV * 2 * (t13 - t14);
            t34 = X_PB * t24 - X_PC * t25 + 0.5 * RHO_INV * 2 * (t14 - t15);
            t35 = X_PB * t25 - X_PC * t26 + 0.5 * RHO_INV * 2 * (t15 - t16);
            t40 = X_PB * t30 - X_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            t41 = X_PB * t31 - X_PC * t32 + 0.5 * RHO_INV * 3 * (t21 - t22);
            t42 = X_PB * t32 - X_PC * t33 + 0.5 * RHO_INV * 3 * (t22 - t23);
            t43 = X_PB * t33 - X_PC * t34 + 0.5 * RHO_INV * 3 * (t23 - t24);
            t44 = X_PB * t34 - X_PC * t35 + 0.5 * RHO_INV * 3 * (t24 - t25);
            *(temp + 0) = beta_in * (*(temp + 0)) + t40;

            t50 = X_PB * t40 - X_PC * t41 + 0.5 * RHO_INV * 4 * (t30 - t31);
            t51 = X_PB * t41 - X_PC * t42 + 0.5 * RHO_INV * 4 * (t31 - t32);
            t52 = X_PB * t42 - X_PC * t43 + 0.5 * RHO_INV * 4 * (t32 - t33);
            t53 = X_PB * t43 - X_PC * t44 + 0.5 * RHO_INV * 4 * (t33 - t34);
            *(temp + 15) = beta_in * (*(temp + 15)) + t50;

            t60 = X_PB * t50 - X_PC * t51 + 0.5 * RHO_INV * 5 * (t40 - t41);
            t61 = X_PB * t51 - X_PC * t52 + 0.5 * RHO_INV * 5 * (t41 - t42);
            t62 = X_PB * t52 - X_PC * t53 + 0.5 * RHO_INV * 5 * (t42 - t43);
            *(temp + 36) = beta_in * (*(temp + 36)) + t60;

            t70 = X_PB * t60 - X_PC * t61 + 0.5 * RHO_INV * 6 * (t50 - t51);
            t71 = X_PB * t61 - X_PC * t62 + 0.5 * RHO_INV * 6 * (t51 - t52);
            *(temp + 64) = beta_in * (*(temp + 64)) + t70;

            t80 = X_PB * t70 - X_PC * t71 + 0.5 * RHO_INV * 7 * (t60 - t61);
            *(temp + 100) = beta_in * (*(temp + 100)) + t80;

            t80 = Y_PB * t70 - Y_PC * t71;
            *(temp + 101) = beta_in * (*(temp + 101)) + t80;

            t80 = Z_PB * t70 - Z_PC * t71;
            *(temp + 102) = beta_in * (*(temp + 102)) + t80;

            t70 = Y_PB * t60 - Y_PC * t61;
            t71 = Y_PB * t61 - Y_PC * t62;
            *(temp + 65) = beta_in * (*(temp + 65)) + t70;

            t80 = Y_PB * t70 - Y_PC * t71 + 0.5 * RHO_INV * 1 * (t60 - t61);
            *(temp + 103) = beta_in * (*(temp + 103)) + t80;

            t80 = Z_PB * t70 - Z_PC * t71;
            *(temp + 104) = beta_in * (*(temp + 104)) + t80;

            t70 = Z_PB * t60 - Z_PC * t61;
            t71 = Z_PB * t61 - Z_PC * t62;
            *(temp + 66) = beta_in * (*(temp + 66)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 1 * (t60 - t61);
            *(temp + 105) = beta_in * (*(temp + 105)) + t80;

            t60 = Y_PB * t50 - Y_PC * t51;
            t61 = Y_PB * t51 - Y_PC * t52;
            t62 = Y_PB * t52 - Y_PC * t53;
            *(temp + 37) = beta_in * (*(temp + 37)) + t60;

            t70 = Y_PB * t60 - Y_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            t71 = Y_PB * t61 - Y_PC * t62 + 0.5 * RHO_INV * 1 * (t51 - t52);
            *(temp + 67) = beta_in * (*(temp + 67)) + t70;

            t80 = Y_PB * t70 - Y_PC * t71 + 0.5 * RHO_INV * 2 * (t60 - t61);
            *(temp + 106) = beta_in * (*(temp + 106)) + t80;

            t80 = Z_PB * t70 - Z_PC * t71;
            *(temp + 107) = beta_in * (*(temp + 107)) + t80;

            t70 = Z_PB * t60 - Z_PC * t61;
            t71 = Z_PB * t61 - Z_PC * t62;
            *(temp + 68) = beta_in * (*(temp + 68)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 1 * (t60 - t61);
            *(temp + 108) = beta_in * (*(temp + 108)) + t80;

            t60 = Z_PB * t50 - Z_PC * t51;
            t61 = Z_PB * t51 - Z_PC * t52;
            t62 = Z_PB * t52 - Z_PC * t53;
            *(temp + 38) = beta_in * (*(temp + 38)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 1 * (t51 - t52);
            *(temp + 69) = beta_in * (*(temp + 69)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 2 * (t60 - t61);
            *(temp + 109) = beta_in * (*(temp + 109)) + t80;

            t50 = Y_PB * t40 - Y_PC * t41;
            t51 = Y_PB * t41 - Y_PC * t42;
            t52 = Y_PB * t42 - Y_PC * t43;
            t53 = Y_PB * t43 - Y_PC * t44;
            *(temp + 16) = beta_in * (*(temp + 16)) + t50;

            t60 = Y_PB * t50 - Y_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            t61 = Y_PB * t51 - Y_PC * t52 + 0.5 * RHO_INV * 1 * (t41 - t42);
            t62 = Y_PB * t52 - Y_PC * t53 + 0.5 * RHO_INV * 1 * (t42 - t43);
            *(temp + 39) = beta_in * (*(temp + 39)) + t60;

            t70 = Y_PB * t60 - Y_PC * t61 + 0.5 * RHO_INV * 2 * (t50 - t51);
            t71 = Y_PB * t61 - Y_PC * t62 + 0.5 * RHO_INV * 2 * (t51 - t52);
            *(temp + 70) = beta_in * (*(temp + 70)) + t70;

            t80 = Y_PB * t70 - Y_PC * t71 + 0.5 * RHO_INV * 3 * (t60 - t61);
            *(temp + 110) = beta_in * (*(temp + 110)) + t80;

            t80 = Z_PB * t70 - Z_PC * t71;
            *(temp + 111) = beta_in * (*(temp + 111)) + t80;

            t70 = Z_PB * t60 - Z_PC * t61;
            t71 = Z_PB * t61 - Z_PC * t62;
            *(temp + 71) = beta_in * (*(temp + 71)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 1 * (t60 - t61);
            *(temp + 112) = beta_in * (*(temp + 112)) + t80;

            t60 = Z_PB * t50 - Z_PC * t51;
            t61 = Z_PB * t51 - Z_PC * t52;
            t62 = Z_PB * t52 - Z_PC * t53;
            *(temp + 40) = beta_in * (*(temp + 40)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 1 * (t51 - t52);
            *(temp + 72) = beta_in * (*(temp + 72)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 2 * (t60 - t61);
            *(temp + 113) = beta_in * (*(temp + 113)) + t80;

            t50 = Z_PB * t40 - Z_PC * t41;
            t51 = Z_PB * t41 - Z_PC * t42;
            t52 = Z_PB * t42 - Z_PC * t43;
            t53 = Z_PB * t43 - Z_PC * t44;
            *(temp + 17) = beta_in * (*(temp + 17)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 1 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 1 * (t42 - t43);
            *(temp + 41) = beta_in * (*(temp + 41)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 2 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 2 * (t51 - t52);
            *(temp + 73) = beta_in * (*(temp + 73)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 3 * (t60 - t61);
            *(temp + 114) = beta_in * (*(temp + 114)) + t80;

            t40 = Y_PB * t30 - Y_PC * t31;
            t41 = Y_PB * t31 - Y_PC * t32;
            t42 = Y_PB * t32 - Y_PC * t33;
            t43 = Y_PB * t33 - Y_PC * t34;
            t44 = Y_PB * t34 - Y_PC * t35;
            *(temp + 1) = beta_in * (*(temp + 1)) + t40;

            t50 = Y_PB * t40 - Y_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Y_PB * t41 - Y_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            t52 = Y_PB * t42 - Y_PC * t43 + 0.5 * RHO_INV * 1 * (t32 - t33);
            t53 = Y_PB * t43 - Y_PC * t44 + 0.5 * RHO_INV * 1 * (t33 - t34);
            *(temp + 18) = beta_in * (*(temp + 18)) + t50;

            t60 = Y_PB * t50 - Y_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            t61 = Y_PB * t51 - Y_PC * t52 + 0.5 * RHO_INV * 2 * (t41 - t42);
            t62 = Y_PB * t52 - Y_PC * t53 + 0.5 * RHO_INV * 2 * (t42 - t43);
            *(temp + 42) = beta_in * (*(temp + 42)) + t60;

            t70 = Y_PB * t60 - Y_PC * t61 + 0.5 * RHO_INV * 3 * (t50 - t51);
            t71 = Y_PB * t61 - Y_PC * t62 + 0.5 * RHO_INV * 3 * (t51 - t52);
            *(temp + 74) = beta_in * (*(temp + 74)) + t70;

            t80 = Y_PB * t70 - Y_PC * t71 + 0.5 * RHO_INV * 4 * (t60 - t61);
            *(temp + 115) = beta_in * (*(temp + 115)) + t80;

            t80 = Z_PB * t70 - Z_PC * t71;
            *(temp + 116) = beta_in * (*(temp + 116)) + t80;

            t70 = Z_PB * t60 - Z_PC * t61;
            t71 = Z_PB * t61 - Z_PC * t62;
            *(temp + 75) = beta_in * (*(temp + 75)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 1 * (t60 - t61);
            *(temp + 117) = beta_in * (*(temp + 117)) + t80;

            t60 = Z_PB * t50 - Z_PC * t51;
            t61 = Z_PB * t51 - Z_PC * t52;
            t62 = Z_PB * t52 - Z_PC * t53;
            *(temp + 43) = beta_in * (*(temp + 43)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 1 * (t51 - t52);
            *(temp + 76) = beta_in * (*(temp + 76)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 2 * (t60 - t61);
            *(temp + 118) = beta_in * (*(temp + 118)) + t80;

            t50 = Z_PB * t40 - Z_PC * t41;
            t51 = Z_PB * t41 - Z_PC * t42;
            t52 = Z_PB * t42 - Z_PC * t43;
            t53 = Z_PB * t43 - Z_PC * t44;
            *(temp + 19) = beta_in * (*(temp + 19)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 1 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 1 * (t42 - t43);
            *(temp + 44) = beta_in * (*(temp + 44)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 2 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 2 * (t51 - t52);
            *(temp + 77) = beta_in * (*(temp + 77)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 3 * (t60 - t61);
            *(temp + 119) = beta_in * (*(temp + 119)) + t80;

            t40 = Z_PB * t30 - Z_PC * t31;
            t41 = Z_PB * t31 - Z_PC * t32;
            t42 = Z_PB * t32 - Z_PC * t33;
            t43 = Z_PB * t33 - Z_PC * t34;
            t44 = Z_PB * t34 - Z_PC * t35;
            *(temp + 2) = beta_in * (*(temp + 2)) + t40;

            t50 = Z_PB * t40 - Z_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Z_PB * t41 - Z_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            t52 = Z_PB * t42 - Z_PC * t43 + 0.5 * RHO_INV * 1 * (t32 - t33);
            t53 = Z_PB * t43 - Z_PC * t44 + 0.5 * RHO_INV * 1 * (t33 - t34);
            *(temp + 20) = beta_in * (*(temp + 20)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 2 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 2 * (t42 - t43);
            *(temp + 45) = beta_in * (*(temp + 45)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 3 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 3 * (t51 - t52);
            *(temp + 78) = beta_in * (*(temp + 78)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 4 * (t60 - t61);
            *(temp + 120) = beta_in * (*(temp + 120)) + t80;

            t30 = Y_PB * t20 - Y_PC * t21;
            t31 = Y_PB * t21 - Y_PC * t22;
            t32 = Y_PB * t22 - Y_PC * t23;
            t33 = Y_PB * t23 - Y_PC * t24;
            t34 = Y_PB * t24 - Y_PC * t25;
            t35 = Y_PB * t25 - Y_PC * t26;
            t40 = Y_PB * t30 - Y_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            t41 = Y_PB * t31 - Y_PC * t32 + 0.5 * RHO_INV * 1 * (t21 - t22);
            t42 = Y_PB * t32 - Y_PC * t33 + 0.5 * RHO_INV * 1 * (t22 - t23);
            t43 = Y_PB * t33 - Y_PC * t34 + 0.5 * RHO_INV * 1 * (t23 - t24);
            t44 = Y_PB * t34 - Y_PC * t35 + 0.5 * RHO_INV * 1 * (t24 - t25);
            *(temp + 3) = beta_in * (*(temp + 3)) + t40;

            t50 = Y_PB * t40 - Y_PC * t41 + 0.5 * RHO_INV * 2 * (t30 - t31);
            t51 = Y_PB * t41 - Y_PC * t42 + 0.5 * RHO_INV * 2 * (t31 - t32);
            t52 = Y_PB * t42 - Y_PC * t43 + 0.5 * RHO_INV * 2 * (t32 - t33);
            t53 = Y_PB * t43 - Y_PC * t44 + 0.5 * RHO_INV * 2 * (t33 - t34);
            *(temp + 21) = beta_in * (*(temp + 21)) + t50;

            t60 = Y_PB * t50 - Y_PC * t51 + 0.5 * RHO_INV * 3 * (t40 - t41);
            t61 = Y_PB * t51 - Y_PC * t52 + 0.5 * RHO_INV * 3 * (t41 - t42);
            t62 = Y_PB * t52 - Y_PC * t53 + 0.5 * RHO_INV * 3 * (t42 - t43);
            *(temp + 46) = beta_in * (*(temp + 46)) + t60;

            t70 = Y_PB * t60 - Y_PC * t61 + 0.5 * RHO_INV * 4 * (t50 - t51);
            t71 = Y_PB * t61 - Y_PC * t62 + 0.5 * RHO_INV * 4 * (t51 - t52);
            *(temp + 79) = beta_in * (*(temp + 79)) + t70;

            t80 = Y_PB * t70 - Y_PC * t71 + 0.5 * RHO_INV * 5 * (t60 - t61);
            *(temp + 121) = beta_in * (*(temp + 121)) + t80;

            t80 = Z_PB * t70 - Z_PC * t71;
            *(temp + 122) = beta_in * (*(temp + 122)) + t80;

            t70 = Z_PB * t60 - Z_PC * t61;
            t71 = Z_PB * t61 - Z_PC * t62;
            *(temp + 80) = beta_in * (*(temp + 80)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 1 * (t60 - t61);
            *(temp + 123) = beta_in * (*(temp + 123)) + t80;

            t60 = Z_PB * t50 - Z_PC * t51;
            t61 = Z_PB * t51 - Z_PC * t52;
            t62 = Z_PB * t52 - Z_PC * t53;
            *(temp + 47) = beta_in * (*(temp + 47)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 1 * (t51 - t52);
            *(temp + 81) = beta_in * (*(temp + 81)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 2 * (t60 - t61);
            *(temp + 124) = beta_in * (*(temp + 124)) + t80;

            t50 = Z_PB * t40 - Z_PC * t41;
            t51 = Z_PB * t41 - Z_PC * t42;
            t52 = Z_PB * t42 - Z_PC * t43;
            t53 = Z_PB * t43 - Z_PC * t44;
            *(temp + 22) = beta_in * (*(temp + 22)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 1 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 1 * (t42 - t43);
            *(temp + 48) = beta_in * (*(temp + 48)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 2 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 2 * (t51 - t52);
            *(temp + 82) = beta_in * (*(temp + 82)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 3 * (t60 - t61);
            *(temp + 125) = beta_in * (*(temp + 125)) + t80;

            t40 = Z_PB * t30 - Z_PC * t31;
            t41 = Z_PB * t31 - Z_PC * t32;
            t42 = Z_PB * t32 - Z_PC * t33;
            t43 = Z_PB * t33 - Z_PC * t34;
            t44 = Z_PB * t34 - Z_PC * t35;
            *(temp + 4) = beta_in * (*(temp + 4)) + t40;

            t50 = Z_PB * t40 - Z_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Z_PB * t41 - Z_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            t52 = Z_PB * t42 - Z_PC * t43 + 0.5 * RHO_INV * 1 * (t32 - t33);
            t53 = Z_PB * t43 - Z_PC * t44 + 0.5 * RHO_INV * 1 * (t33 - t34);
            *(temp + 23) = beta_in * (*(temp + 23)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 2 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 2 * (t42 - t43);
            *(temp + 49) = beta_in * (*(temp + 49)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 3 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 3 * (t51 - t52);
            *(temp + 83) = beta_in * (*(temp + 83)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 4 * (t60 - t61);
            *(temp + 126) = beta_in * (*(temp + 126)) + t80;

            t30 = Z_PB * t20 - Z_PC * t21;
            t31 = Z_PB * t21 - Z_PC * t22;
            t32 = Z_PB * t22 - Z_PC * t23;
            t33 = Z_PB * t23 - Z_PC * t24;
            t34 = Z_PB * t24 - Z_PC * t25;
            t35 = Z_PB * t25 - Z_PC * t26;
            t40 = Z_PB * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            t41 = Z_PB * t31 - Z_PC * t32 + 0.5 * RHO_INV * 1 * (t21 - t22);
            t42 = Z_PB * t32 - Z_PC * t33 + 0.5 * RHO_INV * 1 * (t22 - t23);
            t43 = Z_PB * t33 - Z_PC * t34 + 0.5 * RHO_INV * 1 * (t23 - t24);
            t44 = Z_PB * t34 - Z_PC * t35 + 0.5 * RHO_INV * 1 * (t24 - t25);
            *(temp + 5) = beta_in * (*(temp + 5)) + t40;

            t50 = Z_PB * t40 - Z_PC * t41 + 0.5 * RHO_INV * 2 * (t30 - t31);
            t51 = Z_PB * t41 - Z_PC * t42 + 0.5 * RHO_INV * 2 * (t31 - t32);
            t52 = Z_PB * t42 - Z_PC * t43 + 0.5 * RHO_INV * 2 * (t32 - t33);
            t53 = Z_PB * t43 - Z_PC * t44 + 0.5 * RHO_INV * 2 * (t33 - t34);
            *(temp + 24) = beta_in * (*(temp + 24)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 3 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 3 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 3 * (t42 - t43);
            *(temp + 50) = beta_in * (*(temp + 50)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 4 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 4 * (t51 - t52);
            *(temp + 84) = beta_in * (*(temp + 84)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 5 * (t60 - t61);
            *(temp + 127) = beta_in * (*(temp + 127)) + t80;

            t20 = Y_PB * t10 - Y_PC * t11;
            t21 = Y_PB * t11 - Y_PC * t12;
            t22 = Y_PB * t12 - Y_PC * t13;
            t23 = Y_PB * t13 - Y_PC * t14;
            t24 = Y_PB * t14 - Y_PC * t15;
            t25 = Y_PB * t15 - Y_PC * t16;
            t26 = Y_PB * t16 - Y_PC * t17;
            t30 = Y_PB * t20 - Y_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Y_PB * t21 - Y_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            t32 = Y_PB * t22 - Y_PC * t23 + 0.5 * RHO_INV * 1 * (t12 - t13);
            t33 = Y_PB * t23 - Y_PC * t24 + 0.5 * RHO_INV * 1 * (t13 - t14);
            t34 = Y_PB * t24 - Y_PC * t25 + 0.5 * RHO_INV * 1 * (t14 - t15);
            t35 = Y_PB * t25 - Y_PC * t26 + 0.5 * RHO_INV * 1 * (t15 - t16);
            t40 = Y_PB * t30 - Y_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            t41 = Y_PB * t31 - Y_PC * t32 + 0.5 * RHO_INV * 2 * (t21 - t22);
            t42 = Y_PB * t32 - Y_PC * t33 + 0.5 * RHO_INV * 2 * (t22 - t23);
            t43 = Y_PB * t33 - Y_PC * t34 + 0.5 * RHO_INV * 2 * (t23 - t24);
            t44 = Y_PB * t34 - Y_PC * t35 + 0.5 * RHO_INV * 2 * (t24 - t25);
            *(temp + 6) = beta_in * (*(temp + 6)) + t40;

            t50 = Y_PB * t40 - Y_PC * t41 + 0.5 * RHO_INV * 3 * (t30 - t31);
            t51 = Y_PB * t41 - Y_PC * t42 + 0.5 * RHO_INV * 3 * (t31 - t32);
            t52 = Y_PB * t42 - Y_PC * t43 + 0.5 * RHO_INV * 3 * (t32 - t33);
            t53 = Y_PB * t43 - Y_PC * t44 + 0.5 * RHO_INV * 3 * (t33 - t34);
            *(temp + 25) = beta_in * (*(temp + 25)) + t50;

            t60 = Y_PB * t50 - Y_PC * t51 + 0.5 * RHO_INV * 4 * (t40 - t41);
            t61 = Y_PB * t51 - Y_PC * t52 + 0.5 * RHO_INV * 4 * (t41 - t42);
            t62 = Y_PB * t52 - Y_PC * t53 + 0.5 * RHO_INV * 4 * (t42 - t43);
            *(temp + 51) = beta_in * (*(temp + 51)) + t60;

            t70 = Y_PB * t60 - Y_PC * t61 + 0.5 * RHO_INV * 5 * (t50 - t51);
            t71 = Y_PB * t61 - Y_PC * t62 + 0.5 * RHO_INV * 5 * (t51 - t52);
            *(temp + 85) = beta_in * (*(temp + 85)) + t70;

            t80 = Y_PB * t70 - Y_PC * t71 + 0.5 * RHO_INV * 6 * (t60 - t61);
            *(temp + 128) = beta_in * (*(temp + 128)) + t80;

            t80 = Z_PB * t70 - Z_PC * t71;
            *(temp + 129) = beta_in * (*(temp + 129)) + t80;

            t70 = Z_PB * t60 - Z_PC * t61;
            t71 = Z_PB * t61 - Z_PC * t62;
            *(temp + 86) = beta_in * (*(temp + 86)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 1 * (t60 - t61);
            *(temp + 130) = beta_in * (*(temp + 130)) + t80;

            t60 = Z_PB * t50 - Z_PC * t51;
            t61 = Z_PB * t51 - Z_PC * t52;
            t62 = Z_PB * t52 - Z_PC * t53;
            *(temp + 52) = beta_in * (*(temp + 52)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 1 * (t51 - t52);
            *(temp + 87) = beta_in * (*(temp + 87)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 2 * (t60 - t61);
            *(temp + 131) = beta_in * (*(temp + 131)) + t80;

            t50 = Z_PB * t40 - Z_PC * t41;
            t51 = Z_PB * t41 - Z_PC * t42;
            t52 = Z_PB * t42 - Z_PC * t43;
            t53 = Z_PB * t43 - Z_PC * t44;
            *(temp + 26) = beta_in * (*(temp + 26)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 1 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 1 * (t42 - t43);
            *(temp + 53) = beta_in * (*(temp + 53)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 2 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 2 * (t51 - t52);
            *(temp + 88) = beta_in * (*(temp + 88)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 3 * (t60 - t61);
            *(temp + 132) = beta_in * (*(temp + 132)) + t80;

            t40 = Z_PB * t30 - Z_PC * t31;
            t41 = Z_PB * t31 - Z_PC * t32;
            t42 = Z_PB * t32 - Z_PC * t33;
            t43 = Z_PB * t33 - Z_PC * t34;
            t44 = Z_PB * t34 - Z_PC * t35;
            *(temp + 7) = beta_in * (*(temp + 7)) + t40;

            t50 = Z_PB * t40 - Z_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Z_PB * t41 - Z_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            t52 = Z_PB * t42 - Z_PC * t43 + 0.5 * RHO_INV * 1 * (t32 - t33);
            t53 = Z_PB * t43 - Z_PC * t44 + 0.5 * RHO_INV * 1 * (t33 - t34);
            *(temp + 27) = beta_in * (*(temp + 27)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 2 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 2 * (t42 - t43);
            *(temp + 54) = beta_in * (*(temp + 54)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 3 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 3 * (t51 - t52);
            *(temp + 89) = beta_in * (*(temp + 89)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 4 * (t60 - t61);
            *(temp + 133) = beta_in * (*(temp + 133)) + t80;

            t30 = Z_PB * t20 - Z_PC * t21;
            t31 = Z_PB * t21 - Z_PC * t22;
            t32 = Z_PB * t22 - Z_PC * t23;
            t33 = Z_PB * t23 - Z_PC * t24;
            t34 = Z_PB * t24 - Z_PC * t25;
            t35 = Z_PB * t25 - Z_PC * t26;
            t40 = Z_PB * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            t41 = Z_PB * t31 - Z_PC * t32 + 0.5 * RHO_INV * 1 * (t21 - t22);
            t42 = Z_PB * t32 - Z_PC * t33 + 0.5 * RHO_INV * 1 * (t22 - t23);
            t43 = Z_PB * t33 - Z_PC * t34 + 0.5 * RHO_INV * 1 * (t23 - t24);
            t44 = Z_PB * t34 - Z_PC * t35 + 0.5 * RHO_INV * 1 * (t24 - t25);
            *(temp + 8) = beta_in * (*(temp + 8)) + t40;

            t50 = Z_PB * t40 - Z_PC * t41 + 0.5 * RHO_INV * 2 * (t30 - t31);
            t51 = Z_PB * t41 - Z_PC * t42 + 0.5 * RHO_INV * 2 * (t31 - t32);
            t52 = Z_PB * t42 - Z_PC * t43 + 0.5 * RHO_INV * 2 * (t32 - t33);
            t53 = Z_PB * t43 - Z_PC * t44 + 0.5 * RHO_INV * 2 * (t33 - t34);
            *(temp + 28) = beta_in * (*(temp + 28)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 3 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 3 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 3 * (t42 - t43);
            *(temp + 55) = beta_in * (*(temp + 55)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 4 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 4 * (t51 - t52);
            *(temp + 90) = beta_in * (*(temp + 90)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 5 * (t60 - t61);
            *(temp + 134) = beta_in * (*(temp + 134)) + t80;

            t20 = Z_PB * t10 - Z_PC * t11;
            t21 = Z_PB * t11 - Z_PC * t12;
            t22 = Z_PB * t12 - Z_PC * t13;
            t23 = Z_PB * t13 - Z_PC * t14;
            t24 = Z_PB * t14 - Z_PC * t15;
            t25 = Z_PB * t15 - Z_PC * t16;
            t26 = Z_PB * t16 - Z_PC * t17;
            t30 = Z_PB * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PB * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            t32 = Z_PB * t22 - Z_PC * t23 + 0.5 * RHO_INV * 1 * (t12 - t13);
            t33 = Z_PB * t23 - Z_PC * t24 + 0.5 * RHO_INV * 1 * (t13 - t14);
            t34 = Z_PB * t24 - Z_PC * t25 + 0.5 * RHO_INV * 1 * (t14 - t15);
            t35 = Z_PB * t25 - Z_PC * t26 + 0.5 * RHO_INV * 1 * (t15 - t16);
            t40 = Z_PB * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            t41 = Z_PB * t31 - Z_PC * t32 + 0.5 * RHO_INV * 2 * (t21 - t22);
            t42 = Z_PB * t32 - Z_PC * t33 + 0.5 * RHO_INV * 2 * (t22 - t23);
            t43 = Z_PB * t33 - Z_PC * t34 + 0.5 * RHO_INV * 2 * (t23 - t24);
            t44 = Z_PB * t34 - Z_PC * t35 + 0.5 * RHO_INV * 2 * (t24 - t25);
            *(temp + 9) = beta_in * (*(temp + 9)) + t40;

            t50 = Z_PB * t40 - Z_PC * t41 + 0.5 * RHO_INV * 3 * (t30 - t31);
            t51 = Z_PB * t41 - Z_PC * t42 + 0.5 * RHO_INV * 3 * (t31 - t32);
            t52 = Z_PB * t42 - Z_PC * t43 + 0.5 * RHO_INV * 3 * (t32 - t33);
            t53 = Z_PB * t43 - Z_PC * t44 + 0.5 * RHO_INV * 3 * (t33 - t34);
            *(temp + 29) = beta_in * (*(temp + 29)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 4 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 4 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 4 * (t42 - t43);
            *(temp + 56) = beta_in * (*(temp + 56)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 5 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 5 * (t51 - t52);
            *(temp + 91) = beta_in * (*(temp + 91)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 6 * (t60 - t61);
            *(temp + 135) = beta_in * (*(temp + 135)) + t80;

            t10 = Y_PB * t00 - Y_PC * t01;
            t11 = Y_PB * t01 - Y_PC * t02;
            t12 = Y_PB * t02 - Y_PC * t03;
            t13 = Y_PB * t03 - Y_PC * t04;
            t14 = Y_PB * t04 - Y_PC * t05;
            t15 = Y_PB * t05 - Y_PC * t06;
            t16 = Y_PB * t06 - Y_PC * t07;
            t17 = Y_PB * t07 - Y_PC * t08;
            t20 = Y_PB * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Y_PB * t11 - Y_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Y_PB * t12 - Y_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t23 = Y_PB * t13 - Y_PC * t14 + 0.5 * RHO_INV * 1 * (t03 - t04);
            t24 = Y_PB * t14 - Y_PC * t15 + 0.5 * RHO_INV * 1 * (t04 - t05);
            t25 = Y_PB * t15 - Y_PC * t16 + 0.5 * RHO_INV * 1 * (t05 - t06);
            t26 = Y_PB * t16 - Y_PC * t17 + 0.5 * RHO_INV * 1 * (t06 - t07);
            t30 = Y_PB * t20 - Y_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Y_PB * t21 - Y_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            t32 = Y_PB * t22 - Y_PC * t23 + 0.5 * RHO_INV * 2 * (t12 - t13);
            t33 = Y_PB * t23 - Y_PC * t24 + 0.5 * RHO_INV * 2 * (t13 - t14);
            t34 = Y_PB * t24 - Y_PC * t25 + 0.5 * RHO_INV * 2 * (t14 - t15);
            t35 = Y_PB * t25 - Y_PC * t26 + 0.5 * RHO_INV * 2 * (t15 - t16);
            t40 = Y_PB * t30 - Y_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            t41 = Y_PB * t31 - Y_PC * t32 + 0.5 * RHO_INV * 3 * (t21 - t22);
            t42 = Y_PB * t32 - Y_PC * t33 + 0.5 * RHO_INV * 3 * (t22 - t23);
            t43 = Y_PB * t33 - Y_PC * t34 + 0.5 * RHO_INV * 3 * (t23 - t24);
            t44 = Y_PB * t34 - Y_PC * t35 + 0.5 * RHO_INV * 3 * (t24 - t25);
            *(temp + 10) = beta_in * (*(temp + 10)) + t40;

            t50 = Y_PB * t40 - Y_PC * t41 + 0.5 * RHO_INV * 4 * (t30 - t31);
            t51 = Y_PB * t41 - Y_PC * t42 + 0.5 * RHO_INV * 4 * (t31 - t32);
            t52 = Y_PB * t42 - Y_PC * t43 + 0.5 * RHO_INV * 4 * (t32 - t33);
            t53 = Y_PB * t43 - Y_PC * t44 + 0.5 * RHO_INV * 4 * (t33 - t34);
            *(temp + 30) = beta_in * (*(temp + 30)) + t50;

            t60 = Y_PB * t50 - Y_PC * t51 + 0.5 * RHO_INV * 5 * (t40 - t41);
            t61 = Y_PB * t51 - Y_PC * t52 + 0.5 * RHO_INV * 5 * (t41 - t42);
            t62 = Y_PB * t52 - Y_PC * t53 + 0.5 * RHO_INV * 5 * (t42 - t43);
            *(temp + 57) = beta_in * (*(temp + 57)) + t60;

            t70 = Y_PB * t60 - Y_PC * t61 + 0.5 * RHO_INV * 6 * (t50 - t51);
            t71 = Y_PB * t61 - Y_PC * t62 + 0.5 * RHO_INV * 6 * (t51 - t52);
            *(temp + 92) = beta_in * (*(temp + 92)) + t70;

            t80 = Y_PB * t70 - Y_PC * t71 + 0.5 * RHO_INV * 7 * (t60 - t61);
            *(temp + 136) = beta_in * (*(temp + 136)) + t80;

            t80 = Z_PB * t70 - Z_PC * t71;
            *(temp + 137) = beta_in * (*(temp + 137)) + t80;

            t70 = Z_PB * t60 - Z_PC * t61;
            t71 = Z_PB * t61 - Z_PC * t62;
            *(temp + 93) = beta_in * (*(temp + 93)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 1 * (t60 - t61);
            *(temp + 138) = beta_in * (*(temp + 138)) + t80;

            t60 = Z_PB * t50 - Z_PC * t51;
            t61 = Z_PB * t51 - Z_PC * t52;
            t62 = Z_PB * t52 - Z_PC * t53;
            *(temp + 58) = beta_in * (*(temp + 58)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 1 * (t51 - t52);
            *(temp + 94) = beta_in * (*(temp + 94)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 2 * (t60 - t61);
            *(temp + 139) = beta_in * (*(temp + 139)) + t80;

            t50 = Z_PB * t40 - Z_PC * t41;
            t51 = Z_PB * t41 - Z_PC * t42;
            t52 = Z_PB * t42 - Z_PC * t43;
            t53 = Z_PB * t43 - Z_PC * t44;
            *(temp + 31) = beta_in * (*(temp + 31)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 1 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 1 * (t42 - t43);
            *(temp + 59) = beta_in * (*(temp + 59)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 2 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 2 * (t51 - t52);
            *(temp + 95) = beta_in * (*(temp + 95)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 3 * (t60 - t61);
            *(temp + 140) = beta_in * (*(temp + 140)) + t80;

            t40 = Z_PB * t30 - Z_PC * t31;
            t41 = Z_PB * t31 - Z_PC * t32;
            t42 = Z_PB * t32 - Z_PC * t33;
            t43 = Z_PB * t33 - Z_PC * t34;
            t44 = Z_PB * t34 - Z_PC * t35;
            *(temp + 11) = beta_in * (*(temp + 11)) + t40;

            t50 = Z_PB * t40 - Z_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Z_PB * t41 - Z_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            t52 = Z_PB * t42 - Z_PC * t43 + 0.5 * RHO_INV * 1 * (t32 - t33);
            t53 = Z_PB * t43 - Z_PC * t44 + 0.5 * RHO_INV * 1 * (t33 - t34);
            *(temp + 32) = beta_in * (*(temp + 32)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 2 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 2 * (t42 - t43);
            *(temp + 60) = beta_in * (*(temp + 60)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 3 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 3 * (t51 - t52);
            *(temp + 96) = beta_in * (*(temp + 96)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 4 * (t60 - t61);
            *(temp + 141) = beta_in * (*(temp + 141)) + t80;

            t30 = Z_PB * t20 - Z_PC * t21;
            t31 = Z_PB * t21 - Z_PC * t22;
            t32 = Z_PB * t22 - Z_PC * t23;
            t33 = Z_PB * t23 - Z_PC * t24;
            t34 = Z_PB * t24 - Z_PC * t25;
            t35 = Z_PB * t25 - Z_PC * t26;
            t40 = Z_PB * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            t41 = Z_PB * t31 - Z_PC * t32 + 0.5 * RHO_INV * 1 * (t21 - t22);
            t42 = Z_PB * t32 - Z_PC * t33 + 0.5 * RHO_INV * 1 * (t22 - t23);
            t43 = Z_PB * t33 - Z_PC * t34 + 0.5 * RHO_INV * 1 * (t23 - t24);
            t44 = Z_PB * t34 - Z_PC * t35 + 0.5 * RHO_INV * 1 * (t24 - t25);
            *(temp + 12) = beta_in * (*(temp + 12)) + t40;

            t50 = Z_PB * t40 - Z_PC * t41 + 0.5 * RHO_INV * 2 * (t30 - t31);
            t51 = Z_PB * t41 - Z_PC * t42 + 0.5 * RHO_INV * 2 * (t31 - t32);
            t52 = Z_PB * t42 - Z_PC * t43 + 0.5 * RHO_INV * 2 * (t32 - t33);
            t53 = Z_PB * t43 - Z_PC * t44 + 0.5 * RHO_INV * 2 * (t33 - t34);
            *(temp + 33) = beta_in * (*(temp + 33)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 3 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 3 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 3 * (t42 - t43);
            *(temp + 61) = beta_in * (*(temp + 61)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 4 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 4 * (t51 - t52);
            *(temp + 97) = beta_in * (*(temp + 97)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 5 * (t60 - t61);
            *(temp + 142) = beta_in * (*(temp + 142)) + t80;

            t20 = Z_PB * t10 - Z_PC * t11;
            t21 = Z_PB * t11 - Z_PC * t12;
            t22 = Z_PB * t12 - Z_PC * t13;
            t23 = Z_PB * t13 - Z_PC * t14;
            t24 = Z_PB * t14 - Z_PC * t15;
            t25 = Z_PB * t15 - Z_PC * t16;
            t26 = Z_PB * t16 - Z_PC * t17;
            t30 = Z_PB * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PB * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            t32 = Z_PB * t22 - Z_PC * t23 + 0.5 * RHO_INV * 1 * (t12 - t13);
            t33 = Z_PB * t23 - Z_PC * t24 + 0.5 * RHO_INV * 1 * (t13 - t14);
            t34 = Z_PB * t24 - Z_PC * t25 + 0.5 * RHO_INV * 1 * (t14 - t15);
            t35 = Z_PB * t25 - Z_PC * t26 + 0.5 * RHO_INV * 1 * (t15 - t16);
            t40 = Z_PB * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            t41 = Z_PB * t31 - Z_PC * t32 + 0.5 * RHO_INV * 2 * (t21 - t22);
            t42 = Z_PB * t32 - Z_PC * t33 + 0.5 * RHO_INV * 2 * (t22 - t23);
            t43 = Z_PB * t33 - Z_PC * t34 + 0.5 * RHO_INV * 2 * (t23 - t24);
            t44 = Z_PB * t34 - Z_PC * t35 + 0.5 * RHO_INV * 2 * (t24 - t25);
            *(temp + 13) = beta_in * (*(temp + 13)) + t40;

            t50 = Z_PB * t40 - Z_PC * t41 + 0.5 * RHO_INV * 3 * (t30 - t31);
            t51 = Z_PB * t41 - Z_PC * t42 + 0.5 * RHO_INV * 3 * (t31 - t32);
            t52 = Z_PB * t42 - Z_PC * t43 + 0.5 * RHO_INV * 3 * (t32 - t33);
            t53 = Z_PB * t43 - Z_PC * t44 + 0.5 * RHO_INV * 3 * (t33 - t34);
            *(temp + 34) = beta_in * (*(temp + 34)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 4 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 4 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 4 * (t42 - t43);
            *(temp + 62) = beta_in * (*(temp + 62)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 5 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 5 * (t51 - t52);
            *(temp + 98) = beta_in * (*(temp + 98)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 6 * (t60 - t61);
            *(temp + 143) = beta_in * (*(temp + 143)) + t80;

            t10 = Z_PB * t00 - Z_PC * t01;
            t11 = Z_PB * t01 - Z_PC * t02;
            t12 = Z_PB * t02 - Z_PC * t03;
            t13 = Z_PB * t03 - Z_PC * t04;
            t14 = Z_PB * t04 - Z_PC * t05;
            t15 = Z_PB * t05 - Z_PC * t06;
            t16 = Z_PB * t06 - Z_PC * t07;
            t17 = Z_PB * t07 - Z_PC * t08;
            t20 = Z_PB * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Z_PB * t11 - Z_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Z_PB * t12 - Z_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t23 = Z_PB * t13 - Z_PC * t14 + 0.5 * RHO_INV * 1 * (t03 - t04);
            t24 = Z_PB * t14 - Z_PC * t15 + 0.5 * RHO_INV * 1 * (t04 - t05);
            t25 = Z_PB * t15 - Z_PC * t16 + 0.5 * RHO_INV * 1 * (t05 - t06);
            t26 = Z_PB * t16 - Z_PC * t17 + 0.5 * RHO_INV * 1 * (t06 - t07);
            t30 = Z_PB * t20 - Z_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Z_PB * t21 - Z_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            t32 = Z_PB * t22 - Z_PC * t23 + 0.5 * RHO_INV * 2 * (t12 - t13);
            t33 = Z_PB * t23 - Z_PC * t24 + 0.5 * RHO_INV * 2 * (t13 - t14);
            t34 = Z_PB * t24 - Z_PC * t25 + 0.5 * RHO_INV * 2 * (t14 - t15);
            t35 = Z_PB * t25 - Z_PC * t26 + 0.5 * RHO_INV * 2 * (t15 - t16);
            t40 = Z_PB * t30 - Z_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            t41 = Z_PB * t31 - Z_PC * t32 + 0.5 * RHO_INV * 3 * (t21 - t22);
            t42 = Z_PB * t32 - Z_PC * t33 + 0.5 * RHO_INV * 3 * (t22 - t23);
            t43 = Z_PB * t33 - Z_PC * t34 + 0.5 * RHO_INV * 3 * (t23 - t24);
            t44 = Z_PB * t34 - Z_PC * t35 + 0.5 * RHO_INV * 3 * (t24 - t25);
            *(temp + 14) = beta_in * (*(temp + 14)) + t40;

            t50 = Z_PB * t40 - Z_PC * t41 + 0.5 * RHO_INV * 4 * (t30 - t31);
            t51 = Z_PB * t41 - Z_PC * t42 + 0.5 * RHO_INV * 4 * (t31 - t32);
            t52 = Z_PB * t42 - Z_PC * t43 + 0.5 * RHO_INV * 4 * (t32 - t33);
            t53 = Z_PB * t43 - Z_PC * t44 + 0.5 * RHO_INV * 4 * (t33 - t34);
            *(temp + 35) = beta_in * (*(temp + 35)) + t50;

            t60 = Z_PB * t50 - Z_PC * t51 + 0.5 * RHO_INV * 5 * (t40 - t41);
            t61 = Z_PB * t51 - Z_PC * t52 + 0.5 * RHO_INV * 5 * (t41 - t42);
            t62 = Z_PB * t52 - Z_PC * t53 + 0.5 * RHO_INV * 5 * (t42 - t43);
            *(temp + 63) = beta_in * (*(temp + 63)) + t60;

            t70 = Z_PB * t60 - Z_PC * t61 + 0.5 * RHO_INV * 6 * (t50 - t51);
            t71 = Z_PB * t61 - Z_PC * t62 + 0.5 * RHO_INV * 6 * (t51 - t52);
            *(temp + 99) = beta_in * (*(temp + 99)) + t70;

            t80 = Z_PB * t70 - Z_PC * t71 + 0.5 * RHO_INV * 7 * (t60 - t61);
            *(temp + 144) = beta_in * (*(temp + 144)) + t80;

            beta_in = 1.0;
         }
      }

      double *Xik = (Xi + point_idx * ldX);
      double *Gjk = (Gj + point_idx * ldG);

      for(int c0 = 0; c0 <= 4; ++c0) {
         for(int c1 = 0; c1 <= c0; ++c1) {
            int m = 4 - c0;
            int n = c0 - c1;
            int p = c1;

            int idxB = (((4 - m) * (4 - m + 1)) >> 1) + p;

            double X_ABp = 1.0, comb_m_i = 1.0;
            for(int i = 0; i <= m; ++i) {
               double rcp_i;

               double Y_ABp = 1.0, comb_n_j = 1.0;
               for(int j = 0; j <= n; ++j) {
                  double rcp_j;

                  double Z_ABp = 1.0, comb_p_k = 1.0;
                  for(int k = 0; k <= p; ++k) {
                     double rcp_k;
                     int mv, pv, Lv = 8 - i - j - k;

                     int offset = (Lv * (Lv + 1) * (Lv + 2) - 120) / 6;
                     double const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;

                     mv = 4 + m - i; pv = 0 + p - k;
                     double t0 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 0) += *(Xik + idxB) * t0;

                     mv = 3 + m - i; pv = 0 + p - k;
                     double t1 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 1) += *(Xik + idxB) * t1;

                     mv = 3 + m - i; pv = 1 + p - k;
                     double t2 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 2) += *(Xik + idxB) * t2;

                     mv = 2 + m - i; pv = 0 + p - k;
                     double t3 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 3) += *(Xik + idxB) * t3;

                     mv = 2 + m - i; pv = 1 + p - k;
                     double t4 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 4) += *(Xik + idxB) * t4;

                     mv = 2 + m - i; pv = 2 + p - k;
                     double t5 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 5) += *(Xik + idxB) * t5;

                     mv = 1 + m - i; pv = 0 + p - k;
                     double t6 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 6) += *(Xik + idxB) * t6;

                     mv = 1 + m - i; pv = 1 + p - k;
                     double t7 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 7) += *(Xik + idxB) * t7;

                     mv = 1 + m - i; pv = 2 + p - k;
                     double t8 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 8) += *(Xik + idxB) * t8;

                     mv = 1 + m - i; pv = 3 + p - k;
                     double t9 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 9) += *(Xik + idxB) * t9;

                     mv = 0 + m - i; pv = 0 + p - k;
                     double t10 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 10) += *(Xik + idxB) * t10;

                     mv = 0 + m - i; pv = 1 + p - k;
                     double t11 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 11) += *(Xik + idxB) * t11;

                     mv = 0 + m - i; pv = 2 + p - k;
                     double t12 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 12) += *(Xik + idxB) * t12;

                     mv = 0 + m - i; pv = 3 + p - k;
                     double t13 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 13) += *(Xik + idxB) * t13;

                     mv = 0 + m - i; pv = 4 + p - k;
                     double t14 = *(temp + offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * const_value;
                     *(Gjk + 14) += *(Xik + idxB) * t14;

                     Z_ABp *= (-1.0) * Z_AB; rcp_k = 1.0 / (1.0 * (k + 1)); comb_p_k = (comb_p_k * (p - k)) * rcp_k;
                  }

                  Y_ABp *= (-1.0) * Y_AB; rcp_j = 1.0 / (1.0 * (j + 1)); comb_n_j = (comb_n_j * (n - j)) * rcp_j;
               }

               X_ABp *= (-1.0) * X_AB; rcp_i = 1.0 / (1.0 * (i + 1)); comb_m_i = (comb_m_i * (m - i)) * rcp_i;
            }
         }
      }

      *(Gjk + 0) *= *(weights + point_idx);
      *(Gjk + 1) *= *(weights + point_idx);
      *(Gjk + 2) *= *(weights + point_idx);
      *(Gjk + 3) *= *(weights + point_idx);
      *(Gjk + 4) *= *(weights + point_idx);
      *(Gjk + 5) *= *(weights + point_idx);
      *(Gjk + 6) *= *(weights + point_idx);
      *(Gjk + 7) *= *(weights + point_idx);
      *(Gjk + 8) *= *(weights + point_idx);
      *(Gjk + 9) *= *(weights + point_idx);
      *(Gjk + 10) *= *(weights + point_idx);
      *(Gjk + 11) *= *(weights + point_idx);
      *(Gjk + 12) *= *(weights + point_idx);
      *(Gjk + 13) *= *(weights + point_idx);
      *(Gjk + 14) *= *(weights + point_idx);
   }
}
