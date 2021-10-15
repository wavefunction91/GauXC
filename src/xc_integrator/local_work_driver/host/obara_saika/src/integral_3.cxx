#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

void integral_3(int npts,
               shells shellA,
               point *_points,
               double *Xi,
               int stX,
               int ldX,
               double *Gi,
               int stG, 
               int ldG, 
               double *weights) {
   double temp[74];

   for(int point_idx = 0; point_idx < npts; ++point_idx) {
      point C = *(_points + point_idx);

      double xA = shellA.origin.x;
      double yA = shellA.origin.y;
      double zA = shellA.origin.z;

      double beta_in = 0.0;
      for(int i = 0; i < shellA.m; ++i) {
         for(int j = 0; j < shellA.m; ++j) {
            double aA = shellA.coeff[i].alpha;
            double cA = shellA.coeff[i].coeff;

            double aB = shellA.coeff[j].alpha;
            double cB = shellA.coeff[j].coeff;

            double RHO = aA + aB;
            double RHO_INV = 1.0 / RHO;

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            double X_PA = 0.0;
            double Y_PA = 0.0;
            double Z_PA = 0.0;

            double X_PC = (xA - xC);
            double Y_PC = (yA - yC);
            double Z_PC = (zA - zC);

            double t00, t01, t02, t03, t04, t05, t06, t10, t11, t12, t13, t14, t15, t20, t21, t22, t23, t24, t30, t31, t32, t33, t40, t41, t42, t50, t51, t60;

            double eval = cA * cB * 2 * PI * RHO_INV;
            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

#ifdef BOYS_REFERENCE
            t00 = eval * boys_reference(0, tval);
            t01 = eval * boys_reference(1, tval);
            t02 = eval * boys_reference(2, tval);
            t03 = eval * boys_reference(3, tval);
            t04 = eval * boys_reference(4, tval);
            t05 = eval * boys_reference(5, tval);
            t06 = eval * boys_reference(6, tval);
#elif BOYS_ASYMP
            t00 = eval * boys_asymp(0, tval);
            t01 = eval * boys_asymp(1, tval);
            t02 = eval * boys_asymp(2, tval);
            t03 = eval * boys_asymp(3, tval);
            t04 = eval * boys_asymp(4, tval);
            t05 = eval * boys_asymp(5, tval);
            t06 = eval * boys_asymp(6, tval);
#else
            #error "TYPE NOT DEFINED!"
#endif

            t10 = X_PA * t00 - X_PC * t01;
            t11 = X_PA * t01 - X_PC * t02;
            t12 = X_PA * t02 - X_PC * t03;
            t13 = X_PA * t03 - X_PC * t04;
            t14 = X_PA * t04 - X_PC * t05;
            t15 = X_PA * t05 - X_PC * t06;
            t20 = X_PA * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = X_PA * t11 - X_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = X_PA * t12 - X_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t23 = X_PA * t13 - X_PC * t14 + 0.5 * RHO_INV * 1 * (t03 - t04);
            t24 = X_PA * t14 - X_PC * t15 + 0.5 * RHO_INV * 1 * (t04 - t05);
            t30 = X_PA * t20 - X_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = X_PA * t21 - X_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            t32 = X_PA * t22 - X_PC * t23 + 0.5 * RHO_INV * 2 * (t12 - t13);
            t33 = X_PA * t23 - X_PC * t24 + 0.5 * RHO_INV * 2 * (t13 - t14);
            *(temp + 0) = beta_in * (*(temp + 0)) + t30;

            t40 = X_PA * t30 - X_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            t41 = X_PA * t31 - X_PC * t32 + 0.5 * RHO_INV * 3 * (t21 - t22);
            t42 = X_PA * t32 - X_PC * t33 + 0.5 * RHO_INV * 3 * (t22 - t23);
            *(temp + 10) = beta_in * (*(temp + 10)) + t40;

            t50 = X_PA * t40 - X_PC * t41 + 0.5 * RHO_INV * 4 * (t30 - t31);
            t51 = X_PA * t41 - X_PC * t42 + 0.5 * RHO_INV * 4 * (t31 - t32);
            *(temp + 25) = beta_in * (*(temp + 25)) + t50;

            t60 = X_PA * t50 - X_PC * t51 + 0.5 * RHO_INV * 5 * (t40 - t41);
            *(temp + 46) = beta_in * (*(temp + 46)) + t60;

            t60 = Y_PA * t50 - Y_PC * t51;
            *(temp + 47) = beta_in * (*(temp + 47)) + t60;

            t60 = Z_PA * t50 - Z_PC * t51;
            *(temp + 48) = beta_in * (*(temp + 48)) + t60;

            t50 = Y_PA * t40 - Y_PC * t41;
            t51 = Y_PA * t41 - Y_PC * t42;
            *(temp + 26) = beta_in * (*(temp + 26)) + t50;

            t60 = Y_PA * t50 - Y_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            *(temp + 49) = beta_in * (*(temp + 49)) + t60;

            t60 = Z_PA * t50 - Z_PC * t51;
            *(temp + 50) = beta_in * (*(temp + 50)) + t60;

            t50 = Z_PA * t40 - Z_PC * t41;
            t51 = Z_PA * t41 - Z_PC * t42;
            *(temp + 27) = beta_in * (*(temp + 27)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            *(temp + 51) = beta_in * (*(temp + 51)) + t60;

            t40 = Y_PA * t30 - Y_PC * t31;
            t41 = Y_PA * t31 - Y_PC * t32;
            t42 = Y_PA * t32 - Y_PC * t33;
            *(temp + 11) = beta_in * (*(temp + 11)) + t40;

            t50 = Y_PA * t40 - Y_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Y_PA * t41 - Y_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            *(temp + 28) = beta_in * (*(temp + 28)) + t50;

            t60 = Y_PA * t50 - Y_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            *(temp + 52) = beta_in * (*(temp + 52)) + t60;

            t60 = Z_PA * t50 - Z_PC * t51;
            *(temp + 53) = beta_in * (*(temp + 53)) + t60;

            t50 = Z_PA * t40 - Z_PC * t41;
            t51 = Z_PA * t41 - Z_PC * t42;
            *(temp + 29) = beta_in * (*(temp + 29)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            *(temp + 54) = beta_in * (*(temp + 54)) + t60;

            t40 = Z_PA * t30 - Z_PC * t31;
            t41 = Z_PA * t31 - Z_PC * t32;
            t42 = Z_PA * t32 - Z_PC * t33;
            *(temp + 12) = beta_in * (*(temp + 12)) + t40;

            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            *(temp + 30) = beta_in * (*(temp + 30)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            *(temp + 55) = beta_in * (*(temp + 55)) + t60;

            t30 = Y_PA * t20 - Y_PC * t21;
            t31 = Y_PA * t21 - Y_PC * t22;
            t32 = Y_PA * t22 - Y_PC * t23;
            t33 = Y_PA * t23 - Y_PC * t24;
            *(temp + 1) = beta_in * (*(temp + 1)) + t30;

            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            t41 = Y_PA * t31 - Y_PC * t32 + 0.5 * RHO_INV * 1 * (t21 - t22);
            t42 = Y_PA * t32 - Y_PC * t33 + 0.5 * RHO_INV * 1 * (t22 - t23);
            *(temp + 13) = beta_in * (*(temp + 13)) + t40;

            t50 = Y_PA * t40 - Y_PC * t41 + 0.5 * RHO_INV * 2 * (t30 - t31);
            t51 = Y_PA * t41 - Y_PC * t42 + 0.5 * RHO_INV * 2 * (t31 - t32);
            *(temp + 31) = beta_in * (*(temp + 31)) + t50;

            t60 = Y_PA * t50 - Y_PC * t51 + 0.5 * RHO_INV * 3 * (t40 - t41);
            *(temp + 56) = beta_in * (*(temp + 56)) + t60;

            t60 = Z_PA * t50 - Z_PC * t51;
            *(temp + 57) = beta_in * (*(temp + 57)) + t60;

            t50 = Z_PA * t40 - Z_PC * t41;
            t51 = Z_PA * t41 - Z_PC * t42;
            *(temp + 32) = beta_in * (*(temp + 32)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            *(temp + 58) = beta_in * (*(temp + 58)) + t60;

            t40 = Z_PA * t30 - Z_PC * t31;
            t41 = Z_PA * t31 - Z_PC * t32;
            t42 = Z_PA * t32 - Z_PC * t33;
            *(temp + 14) = beta_in * (*(temp + 14)) + t40;

            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            *(temp + 33) = beta_in * (*(temp + 33)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            *(temp + 59) = beta_in * (*(temp + 59)) + t60;

            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            t32 = Z_PA * t22 - Z_PC * t23;
            t33 = Z_PA * t23 - Z_PC * t24;
            *(temp + 2) = beta_in * (*(temp + 2)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            t41 = Z_PA * t31 - Z_PC * t32 + 0.5 * RHO_INV * 1 * (t21 - t22);
            t42 = Z_PA * t32 - Z_PC * t33 + 0.5 * RHO_INV * 1 * (t22 - t23);
            *(temp + 15) = beta_in * (*(temp + 15)) + t40;

            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 2 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 2 * (t31 - t32);
            *(temp + 34) = beta_in * (*(temp + 34)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 3 * (t40 - t41);
            *(temp + 60) = beta_in * (*(temp + 60)) + t60;

            t20 = Y_PA * t10 - Y_PC * t11;
            t21 = Y_PA * t11 - Y_PC * t12;
            t22 = Y_PA * t12 - Y_PC * t13;
            t23 = Y_PA * t13 - Y_PC * t14;
            t24 = Y_PA * t14 - Y_PC * t15;
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Y_PA * t21 - Y_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            t32 = Y_PA * t22 - Y_PC * t23 + 0.5 * RHO_INV * 1 * (t12 - t13);
            t33 = Y_PA * t23 - Y_PC * t24 + 0.5 * RHO_INV * 1 * (t13 - t14);
            *(temp + 3) = beta_in * (*(temp + 3)) + t30;

            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            t41 = Y_PA * t31 - Y_PC * t32 + 0.5 * RHO_INV * 2 * (t21 - t22);
            t42 = Y_PA * t32 - Y_PC * t33 + 0.5 * RHO_INV * 2 * (t22 - t23);
            *(temp + 16) = beta_in * (*(temp + 16)) + t40;

            t50 = Y_PA * t40 - Y_PC * t41 + 0.5 * RHO_INV * 3 * (t30 - t31);
            t51 = Y_PA * t41 - Y_PC * t42 + 0.5 * RHO_INV * 3 * (t31 - t32);
            *(temp + 35) = beta_in * (*(temp + 35)) + t50;

            t60 = Y_PA * t50 - Y_PC * t51 + 0.5 * RHO_INV * 4 * (t40 - t41);
            *(temp + 61) = beta_in * (*(temp + 61)) + t60;

            t60 = Z_PA * t50 - Z_PC * t51;
            *(temp + 62) = beta_in * (*(temp + 62)) + t60;

            t50 = Z_PA * t40 - Z_PC * t41;
            t51 = Z_PA * t41 - Z_PC * t42;
            *(temp + 36) = beta_in * (*(temp + 36)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            *(temp + 63) = beta_in * (*(temp + 63)) + t60;

            t40 = Z_PA * t30 - Z_PC * t31;
            t41 = Z_PA * t31 - Z_PC * t32;
            t42 = Z_PA * t32 - Z_PC * t33;
            *(temp + 17) = beta_in * (*(temp + 17)) + t40;

            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            *(temp + 37) = beta_in * (*(temp + 37)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            *(temp + 64) = beta_in * (*(temp + 64)) + t60;

            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            t32 = Z_PA * t22 - Z_PC * t23;
            t33 = Z_PA * t23 - Z_PC * t24;
            *(temp + 4) = beta_in * (*(temp + 4)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            t41 = Z_PA * t31 - Z_PC * t32 + 0.5 * RHO_INV * 1 * (t21 - t22);
            t42 = Z_PA * t32 - Z_PC * t33 + 0.5 * RHO_INV * 1 * (t22 - t23);
            *(temp + 18) = beta_in * (*(temp + 18)) + t40;

            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 2 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 2 * (t31 - t32);
            *(temp + 38) = beta_in * (*(temp + 38)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 3 * (t40 - t41);
            *(temp + 65) = beta_in * (*(temp + 65)) + t60;

            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t22 = Z_PA * t12 - Z_PC * t13;
            t23 = Z_PA * t13 - Z_PC * t14;
            t24 = Z_PA * t14 - Z_PC * t15;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            t32 = Z_PA * t22 - Z_PC * t23 + 0.5 * RHO_INV * 1 * (t12 - t13);
            t33 = Z_PA * t23 - Z_PC * t24 + 0.5 * RHO_INV * 1 * (t13 - t14);
            *(temp + 5) = beta_in * (*(temp + 5)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            t41 = Z_PA * t31 - Z_PC * t32 + 0.5 * RHO_INV * 2 * (t21 - t22);
            t42 = Z_PA * t32 - Z_PC * t33 + 0.5 * RHO_INV * 2 * (t22 - t23);
            *(temp + 19) = beta_in * (*(temp + 19)) + t40;

            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 3 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 3 * (t31 - t32);
            *(temp + 39) = beta_in * (*(temp + 39)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 4 * (t40 - t41);
            *(temp + 66) = beta_in * (*(temp + 66)) + t60;

            t10 = Y_PA * t00 - Y_PC * t01;
            t11 = Y_PA * t01 - Y_PC * t02;
            t12 = Y_PA * t02 - Y_PC * t03;
            t13 = Y_PA * t03 - Y_PC * t04;
            t14 = Y_PA * t04 - Y_PC * t05;
            t15 = Y_PA * t05 - Y_PC * t06;
            t20 = Y_PA * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Y_PA * t11 - Y_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Y_PA * t12 - Y_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t23 = Y_PA * t13 - Y_PC * t14 + 0.5 * RHO_INV * 1 * (t03 - t04);
            t24 = Y_PA * t14 - Y_PC * t15 + 0.5 * RHO_INV * 1 * (t04 - t05);
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Y_PA * t21 - Y_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            t32 = Y_PA * t22 - Y_PC * t23 + 0.5 * RHO_INV * 2 * (t12 - t13);
            t33 = Y_PA * t23 - Y_PC * t24 + 0.5 * RHO_INV * 2 * (t13 - t14);
            *(temp + 6) = beta_in * (*(temp + 6)) + t30;

            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            t41 = Y_PA * t31 - Y_PC * t32 + 0.5 * RHO_INV * 3 * (t21 - t22);
            t42 = Y_PA * t32 - Y_PC * t33 + 0.5 * RHO_INV * 3 * (t22 - t23);
            *(temp + 20) = beta_in * (*(temp + 20)) + t40;

            t50 = Y_PA * t40 - Y_PC * t41 + 0.5 * RHO_INV * 4 * (t30 - t31);
            t51 = Y_PA * t41 - Y_PC * t42 + 0.5 * RHO_INV * 4 * (t31 - t32);
            *(temp + 40) = beta_in * (*(temp + 40)) + t50;

            t60 = Y_PA * t50 - Y_PC * t51 + 0.5 * RHO_INV * 5 * (t40 - t41);
            *(temp + 67) = beta_in * (*(temp + 67)) + t60;

            t60 = Z_PA * t50 - Z_PC * t51;
            *(temp + 68) = beta_in * (*(temp + 68)) + t60;

            t50 = Z_PA * t40 - Z_PC * t41;
            t51 = Z_PA * t41 - Z_PC * t42;
            *(temp + 41) = beta_in * (*(temp + 41)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            *(temp + 69) = beta_in * (*(temp + 69)) + t60;

            t40 = Z_PA * t30 - Z_PC * t31;
            t41 = Z_PA * t31 - Z_PC * t32;
            t42 = Z_PA * t32 - Z_PC * t33;
            *(temp + 21) = beta_in * (*(temp + 21)) + t40;

            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            *(temp + 42) = beta_in * (*(temp + 42)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            *(temp + 70) = beta_in * (*(temp + 70)) + t60;

            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            t32 = Z_PA * t22 - Z_PC * t23;
            t33 = Z_PA * t23 - Z_PC * t24;
            *(temp + 7) = beta_in * (*(temp + 7)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            t41 = Z_PA * t31 - Z_PC * t32 + 0.5 * RHO_INV * 1 * (t21 - t22);
            t42 = Z_PA * t32 - Z_PC * t33 + 0.5 * RHO_INV * 1 * (t22 - t23);
            *(temp + 22) = beta_in * (*(temp + 22)) + t40;

            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 2 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 2 * (t31 - t32);
            *(temp + 43) = beta_in * (*(temp + 43)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 3 * (t40 - t41);
            *(temp + 71) = beta_in * (*(temp + 71)) + t60;

            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t22 = Z_PA * t12 - Z_PC * t13;
            t23 = Z_PA * t13 - Z_PC * t14;
            t24 = Z_PA * t14 - Z_PC * t15;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            t32 = Z_PA * t22 - Z_PC * t23 + 0.5 * RHO_INV * 1 * (t12 - t13);
            t33 = Z_PA * t23 - Z_PC * t24 + 0.5 * RHO_INV * 1 * (t13 - t14);
            *(temp + 8) = beta_in * (*(temp + 8)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            t41 = Z_PA * t31 - Z_PC * t32 + 0.5 * RHO_INV * 2 * (t21 - t22);
            t42 = Z_PA * t32 - Z_PC * t33 + 0.5 * RHO_INV * 2 * (t22 - t23);
            *(temp + 23) = beta_in * (*(temp + 23)) + t40;

            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 3 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 3 * (t31 - t32);
            *(temp + 44) = beta_in * (*(temp + 44)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 4 * (t40 - t41);
            *(temp + 72) = beta_in * (*(temp + 72)) + t60;

            t10 = Z_PA * t00 - Z_PC * t01;
            t11 = Z_PA * t01 - Z_PC * t02;
            t12 = Z_PA * t02 - Z_PC * t03;
            t13 = Z_PA * t03 - Z_PC * t04;
            t14 = Z_PA * t04 - Z_PC * t05;
            t15 = Z_PA * t05 - Z_PC * t06;
            t20 = Z_PA * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Z_PA * t11 - Z_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Z_PA * t12 - Z_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t23 = Z_PA * t13 - Z_PC * t14 + 0.5 * RHO_INV * 1 * (t03 - t04);
            t24 = Z_PA * t14 - Z_PC * t15 + 0.5 * RHO_INV * 1 * (t04 - t05);
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            t32 = Z_PA * t22 - Z_PC * t23 + 0.5 * RHO_INV * 2 * (t12 - t13);
            t33 = Z_PA * t23 - Z_PC * t24 + 0.5 * RHO_INV * 2 * (t13 - t14);
            *(temp + 9) = beta_in * (*(temp + 9)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            t41 = Z_PA * t31 - Z_PC * t32 + 0.5 * RHO_INV * 3 * (t21 - t22);
            t42 = Z_PA * t32 - Z_PC * t33 + 0.5 * RHO_INV * 3 * (t22 - t23);
            *(temp + 24) = beta_in * (*(temp + 24)) + t40;

            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 4 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 4 * (t31 - t32);
            *(temp + 45) = beta_in * (*(temp + 45)) + t50;

            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 5 * (t40 - t41);
            *(temp + 73) = beta_in * (*(temp + 73)) + t60;

            beta_in = 1.0;
         }
      }

      double *Xik = (Xi + point_idx * stX);
      double *Gik = (Gi + point_idx * stG);

      for(int c0 = 0; c0 <= 3; ++c0) {
         for(int c1 = 0; c1 <= c0; ++c1) {
            int m = 3 - c0;
            int n = c0 - c1;
            int p = c1;

            int idxB = (((3 - m) * (3 - m + 1)) >> 1) + p;

            int mv, pv;

            mv = 3 + m; pv = 0 + p;
            double t0 = *(temp + 46 + (((6 - mv) * (6 - mv + 1)) >> 1) + pv);
            *(Gik + 0 * ldG) += *(Xik + idxB * ldX) * t0;

            mv = 2 + m; pv = 0 + p;
            double t1 = *(temp + 46 + (((6 - mv) * (6 - mv + 1)) >> 1) + pv);
            *(Gik + 1 * ldG) += *(Xik + idxB * ldX) * t1;

            mv = 2 + m; pv = 1 + p;
            double t2 = *(temp + 46 + (((6 - mv) * (6 - mv + 1)) >> 1) + pv);
            *(Gik + 2 * ldG) += *(Xik + idxB * ldX) * t2;

            mv = 1 + m; pv = 0 + p;
            double t3 = *(temp + 46 + (((6 - mv) * (6 - mv + 1)) >> 1) + pv);
            *(Gik + 3 * ldG) += *(Xik + idxB * ldX) * t3;

            mv = 1 + m; pv = 1 + p;
            double t4 = *(temp + 46 + (((6 - mv) * (6 - mv + 1)) >> 1) + pv);
            *(Gik + 4 * ldG) += *(Xik + idxB * ldX) * t4;

            mv = 1 + m; pv = 2 + p;
            double t5 = *(temp + 46 + (((6 - mv) * (6 - mv + 1)) >> 1) + pv);
            *(Gik + 5 * ldG) += *(Xik + idxB * ldX) * t5;

            mv = 0 + m; pv = 0 + p;
            double t6 = *(temp + 46 + (((6 - mv) * (6 - mv + 1)) >> 1) + pv);
            *(Gik + 6 * ldG) += *(Xik + idxB * ldX) * t6;

            mv = 0 + m; pv = 1 + p;
            double t7 = *(temp + 46 + (((6 - mv) * (6 - mv + 1)) >> 1) + pv);
            *(Gik + 7 * ldG) += *(Xik + idxB * ldX) * t7;

            mv = 0 + m; pv = 2 + p;
            double t8 = *(temp + 46 + (((6 - mv) * (6 - mv + 1)) >> 1) + pv);
            *(Gik + 8 * ldG) += *(Xik + idxB * ldX) * t8;

            mv = 0 + m; pv = 3 + p;
            double t9 = *(temp + 46 + (((6 - mv) * (6 - mv + 1)) >> 1) + pv);
            *(Gik + 9 * ldG) += *(Xik + idxB * ldX) * t9;
         }
      }

      *(Gik + 0 * ldG) *= *(weights + point_idx);
      *(Gik + 1 * ldG) *= *(weights + point_idx);
      *(Gik + 2 * ldG) *= *(weights + point_idx);
      *(Gik + 3 * ldG) *= *(weights + point_idx);
      *(Gik + 4 * ldG) *= *(weights + point_idx);
      *(Gik + 5 * ldG) *= *(weights + point_idx);
      *(Gik + 6 * ldG) *= *(weights + point_idx);
      *(Gik + 7 * ldG) *= *(weights + point_idx);
      *(Gik + 8 * ldG) *= *(weights + point_idx);
      *(Gik + 9 * ldG) *= *(weights + point_idx);
   }
}
