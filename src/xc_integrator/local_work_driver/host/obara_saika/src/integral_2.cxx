#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

void integral_2(int npts,
               shells shellA,
               point *_points,
               double *Xi,
               int stX,
               int ldX,
               double *Gi,
               int stG, 
               int ldG, 
               double *weights) {
   double temp[31];

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

            double t00, t01, t02, t03, t04, t10, t11, t12, t13, t20, t21, t22, t30, t31, t40;

            double eval = cA * cB * 2 * PI * RHO_INV;
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

            t10 = X_PA * t00 - X_PC * t01;
            t11 = X_PA * t01 - X_PC * t02;
            t12 = X_PA * t02 - X_PC * t03;
            t13 = X_PA * t03 - X_PC * t04;
            t20 = X_PA * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = X_PA * t11 - X_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = X_PA * t12 - X_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            *(temp + 0) = beta_in * (*(temp + 0)) + t20;

            t30 = X_PA * t20 - X_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = X_PA * t21 - X_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            *(temp + 6) = beta_in * (*(temp + 6)) + t30;

            t40 = X_PA * t30 - X_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 16) = beta_in * (*(temp + 16)) + t40;

            t40 = Y_PA * t30 - Y_PC * t31;
            *(temp + 17) = beta_in * (*(temp + 17)) + t40;

            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 18) = beta_in * (*(temp + 18)) + t40;

            t30 = Y_PA * t20 - Y_PC * t21;
            t31 = Y_PA * t21 - Y_PC * t22;
            *(temp + 7) = beta_in * (*(temp + 7)) + t30;

            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 19) = beta_in * (*(temp + 19)) + t40;

            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 20) = beta_in * (*(temp + 20)) + t40;

            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            *(temp + 8) = beta_in * (*(temp + 8)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 21) = beta_in * (*(temp + 21)) + t40;

            t20 = Y_PA * t10 - Y_PC * t11;
            t21 = Y_PA * t11 - Y_PC * t12;
            t22 = Y_PA * t12 - Y_PC * t13;
            *(temp + 1) = beta_in * (*(temp + 1)) + t20;

            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Y_PA * t21 - Y_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            *(temp + 9) = beta_in * (*(temp + 9)) + t30;

            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 22) = beta_in * (*(temp + 22)) + t40;

            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 23) = beta_in * (*(temp + 23)) + t40;

            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            *(temp + 10) = beta_in * (*(temp + 10)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 24) = beta_in * (*(temp + 24)) + t40;

            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t22 = Z_PA * t12 - Z_PC * t13;
            *(temp + 2) = beta_in * (*(temp + 2)) + t20;

            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            *(temp + 11) = beta_in * (*(temp + 11)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 25) = beta_in * (*(temp + 25)) + t40;

            t10 = Y_PA * t00 - Y_PC * t01;
            t11 = Y_PA * t01 - Y_PC * t02;
            t12 = Y_PA * t02 - Y_PC * t03;
            t13 = Y_PA * t03 - Y_PC * t04;
            t20 = Y_PA * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Y_PA * t11 - Y_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Y_PA * t12 - Y_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            *(temp + 3) = beta_in * (*(temp + 3)) + t20;

            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Y_PA * t21 - Y_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            *(temp + 12) = beta_in * (*(temp + 12)) + t30;

            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 26) = beta_in * (*(temp + 26)) + t40;

            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 27) = beta_in * (*(temp + 27)) + t40;

            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            *(temp + 13) = beta_in * (*(temp + 13)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 28) = beta_in * (*(temp + 28)) + t40;

            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t22 = Z_PA * t12 - Z_PC * t13;
            *(temp + 4) = beta_in * (*(temp + 4)) + t20;

            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            *(temp + 14) = beta_in * (*(temp + 14)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 29) = beta_in * (*(temp + 29)) + t40;

            t10 = Z_PA * t00 - Z_PC * t01;
            t11 = Z_PA * t01 - Z_PC * t02;
            t12 = Z_PA * t02 - Z_PC * t03;
            t13 = Z_PA * t03 - Z_PC * t04;
            t20 = Z_PA * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Z_PA * t11 - Z_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Z_PA * t12 - Z_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            *(temp + 5) = beta_in * (*(temp + 5)) + t20;

            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            *(temp + 15) = beta_in * (*(temp + 15)) + t30;

            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 30) = beta_in * (*(temp + 30)) + t40;

            beta_in = 1.0;
         }
      }

      double *Xik = (Xi + point_idx * stX);
      double *Gik = (Gi + point_idx * stG);

      double t0, t1, t2, t3, t4, t5;

      t0 = *(temp + 16);
      *(Gik + 0 * ldG) += *(Xik + 0 * ldX) * t0;
      t1 = *(temp + 17);
      *(Gik + 1 * ldG) += *(Xik + 0 * ldX) * t1;
      t2 = *(temp + 18);
      *(Gik + 2 * ldG) += *(Xik + 0 * ldX) * t2;
      t3 = *(temp + 19);
      *(Gik + 3 * ldG) += *(Xik + 0 * ldX) * t3;
      t4 = *(temp + 20);
      *(Gik + 4 * ldG) += *(Xik + 0 * ldX) * t4;
      t5 = *(temp + 21);
      *(Gik + 5 * ldG) += *(Xik + 0 * ldX) * t5;

      t0 = *(temp + 17);
      *(Gik + 0 * ldG) += *(Xik + 1 * ldX) * t0;
      t1 = *(temp + 19);
      *(Gik + 1 * ldG) += *(Xik + 1 * ldX) * t1;
      t2 = *(temp + 20);
      *(Gik + 2 * ldG) += *(Xik + 1 * ldX) * t2;
      t3 = *(temp + 22);
      *(Gik + 3 * ldG) += *(Xik + 1 * ldX) * t3;
      t4 = *(temp + 23);
      *(Gik + 4 * ldG) += *(Xik + 1 * ldX) * t4;
      t5 = *(temp + 24);
      *(Gik + 5 * ldG) += *(Xik + 1 * ldX) * t5;

      t0 = *(temp + 18);
      *(Gik + 0 * ldG) += *(Xik + 2 * ldX) * t0;
      t1 = *(temp + 20);
      *(Gik + 1 * ldG) += *(Xik + 2 * ldX) * t1;
      t2 = *(temp + 21);
      *(Gik + 2 * ldG) += *(Xik + 2 * ldX) * t2;
      t3 = *(temp + 23);
      *(Gik + 3 * ldG) += *(Xik + 2 * ldX) * t3;
      t4 = *(temp + 24);
      *(Gik + 4 * ldG) += *(Xik + 2 * ldX) * t4;
      t5 = *(temp + 25);
      *(Gik + 5 * ldG) += *(Xik + 2 * ldX) * t5;

      t0 = *(temp + 19);
      *(Gik + 0 * ldG) += *(Xik + 3 * ldX) * t0;
      t1 = *(temp + 22);
      *(Gik + 1 * ldG) += *(Xik + 3 * ldX) * t1;
      t2 = *(temp + 23);
      *(Gik + 2 * ldG) += *(Xik + 3 * ldX) * t2;
      t3 = *(temp + 26);
      *(Gik + 3 * ldG) += *(Xik + 3 * ldX) * t3;
      t4 = *(temp + 27);
      *(Gik + 4 * ldG) += *(Xik + 3 * ldX) * t4;
      t5 = *(temp + 28);
      *(Gik + 5 * ldG) += *(Xik + 3 * ldX) * t5;

      t0 = *(temp + 20);
      *(Gik + 0 * ldG) += *(Xik + 4 * ldX) * t0;
      t1 = *(temp + 23);
      *(Gik + 1 * ldG) += *(Xik + 4 * ldX) * t1;
      t2 = *(temp + 24);
      *(Gik + 2 * ldG) += *(Xik + 4 * ldX) * t2;
      t3 = *(temp + 27);
      *(Gik + 3 * ldG) += *(Xik + 4 * ldX) * t3;
      t4 = *(temp + 28);
      *(Gik + 4 * ldG) += *(Xik + 4 * ldX) * t4;
      t5 = *(temp + 29);
      *(Gik + 5 * ldG) += *(Xik + 4 * ldX) * t5;

      t0 = *(temp + 21);
      *(Gik + 0 * ldG) += *(Xik + 5 * ldX) * t0;
      t1 = *(temp + 24);
      *(Gik + 1 * ldG) += *(Xik + 5 * ldX) * t1;
      t2 = *(temp + 25);
      *(Gik + 2 * ldG) += *(Xik + 5 * ldX) * t2;
      t3 = *(temp + 28);
      *(Gik + 3 * ldG) += *(Xik + 5 * ldX) * t3;
      t4 = *(temp + 29);
      *(Gik + 4 * ldG) += *(Xik + 5 * ldX) * t4;
      t5 = *(temp + 30);
      *(Gik + 5 * ldG) += *(Xik + 5 * ldX) * t5;

      *(Gik + 0 * ldG) *= *(weights + point_idx);
      *(Gik + 1 * ldG) *= *(weights + point_idx);
      *(Gik + 2 * ldG) *= *(weights + point_idx);
      *(Gik + 3 * ldG) *= *(weights + point_idx);
      *(Gik + 4 * ldG) *= *(weights + point_idx);
      *(Gik + 5 * ldG) *= *(weights + point_idx);
   }
}
