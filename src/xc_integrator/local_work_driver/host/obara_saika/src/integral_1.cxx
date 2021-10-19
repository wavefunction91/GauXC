#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

void integral_1(size_t npts,
               shells shellA,
               point *_points,
               double *Xi,
               int stX,
               int ldX,
               double *Gi,
               int stG, 
               int ldG, 
               double *weights) {
   double temp[9];

   for(int i = 0; i < 9; ++i) {
      temp[i] = 0.0;
   }

   for(size_t point_idx = 0; point_idx < npts; ++point_idx) {
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

            double t00, t01, t02, t10, t11, t20;

            double eval = cA * cB * 2 * PI * RHO_INV;
            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

            t00 = eval * boys_function(0, tval);
            t01 = eval * boys_function(1, tval);
            t02 = eval * boys_function(2, tval);
            t10 = X_PA * t00 - X_PC * t01;
            t11 = X_PA * t01 - X_PC * t02;
            *(temp + 0) = beta_in * (*(temp + 0)) + t10;

            t20 = X_PA * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            *(temp + 3) = beta_in * (*(temp + 3)) + t20;

            t20 = Y_PA * t10 - Y_PC * t11;
            *(temp + 4) = beta_in * (*(temp + 4)) + t20;

            t20 = Z_PA * t10 - Z_PC * t11;
            *(temp + 5) = beta_in * (*(temp + 5)) + t20;

            t10 = Y_PA * t00 - Y_PC * t01;
            t11 = Y_PA * t01 - Y_PC * t02;
            *(temp + 1) = beta_in * (*(temp + 1)) + t10;

            t20 = Y_PA * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            *(temp + 6) = beta_in * (*(temp + 6)) + t20;

            t20 = Z_PA * t10 - Z_PC * t11;
            *(temp + 7) = beta_in * (*(temp + 7)) + t20;

            t10 = Z_PA * t00 - Z_PC * t01;
            t11 = Z_PA * t01 - Z_PC * t02;
            *(temp + 2) = beta_in * (*(temp + 2)) + t10;

            t20 = Z_PA * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            *(temp + 8) = beta_in * (*(temp + 8)) + t20;

            beta_in = 1.0;
         }
      }

      double *Xik = (Xi + point_idx * stX);
      double *Gik = (Gi + point_idx * stG);

      double t0, t1, t2;

      t0 = *(temp + 3) * (*(weights + point_idx));
      *(Gik + 0 * ldG) += *(Xik + 0 * ldX) * t0;
      t1 = *(temp + 4) * (*(weights + point_idx));
      *(Gik + 1 * ldG) += *(Xik + 0 * ldX) * t1;
      t2 = *(temp + 5) * (*(weights + point_idx));
      *(Gik + 2 * ldG) += *(Xik + 0 * ldX) * t2;

      t0 = *(temp + 4) * (*(weights + point_idx));
      *(Gik + 0 * ldG) += *(Xik + 1 * ldX) * t0;
      t1 = *(temp + 6) * (*(weights + point_idx));
      *(Gik + 1 * ldG) += *(Xik + 1 * ldX) * t1;
      t2 = *(temp + 7) * (*(weights + point_idx));
      *(Gik + 2 * ldG) += *(Xik + 1 * ldX) * t2;

      t0 = *(temp + 5) * (*(weights + point_idx));
      *(Gik + 0 * ldG) += *(Xik + 2 * ldX) * t0;
      t1 = *(temp + 7) * (*(weights + point_idx));
      *(Gik + 1 * ldG) += *(Xik + 2 * ldX) * t1;
      t2 = *(temp + 8) * (*(weights + point_idx));
      *(Gik + 2 * ldG) += *(Xik + 2 * ldX) * t2;
   }
}
