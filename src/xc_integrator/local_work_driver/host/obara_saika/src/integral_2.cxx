#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

void integral_2(size_t npts,
               shell_pair shpair,
               point *_points,
               double *Xi,
               int stX,
               int ldX,
               double *Gi,
               int stG, 
               int ldG, 
               double *weights) {
   double temp[31 * NPTS_LOCAL];

   for(int i = 0; i < 31 * NPTS_LOCAL; ++i) {
      temp[i] = 0.0;
   }

   for(size_t p_outer = 0; p_outer < npts; p_outer += NPTS_LOCAL) {
      size_t npts_inner = MIN(NPTS_LOCAL, npts - p_outer);
      point *_point_outer = (_points + p_outer);

      double xA = shpair.rA.x;
      double yA = shpair.rA.y;
      double zA = shpair.rA.z;

      for( int ij = 0; ij < shpair.nprim_pair; ++ij ) {
         double RHO = shpair.prim_pairs[ij].gamma;
         double RHO_INV = 1.0 / RHO;

         constexpr double X_PA = 0.0;
         constexpr double Y_PA = 0.0;
         constexpr double Z_PA = 0.0;

         double eval = shpair.prim_pairs[ij].coeff_prod * 2 * PI * RHO_INV;

         for(size_t p_inner = 0; p_inner < npts_inner; ++p_inner) {
            point C = *(_point_outer + p_inner);

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            double X_PC = (xA - xC);
            double Y_PC = (yA - yC);
            double Z_PC = (zA - zC);

            double t00, t01, t02, t03, t04, t10, t11, t12, t13, t20, t21, t22, t30, t31, t40;

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
            *(temp + 0 * NPTS_LOCAL + p_inner) += t20;
            t30 = X_PA * t20 - X_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = X_PA * t21 - X_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            *(temp + 6 * NPTS_LOCAL + p_inner) += t30;
            t40 = X_PA * t30 - X_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 16 * NPTS_LOCAL + p_inner) += t40;
            t40 = Y_PA * t30 - Y_PC * t31;
            *(temp + 17 * NPTS_LOCAL + p_inner) += t40;
            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 18 * NPTS_LOCAL + p_inner) += t40;
            t30 = Y_PA * t20 - Y_PC * t21;
            t31 = Y_PA * t21 - Y_PC * t22;
            *(temp + 7 * NPTS_LOCAL + p_inner) += t30;
            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 19 * NPTS_LOCAL + p_inner) += t40;
            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 20 * NPTS_LOCAL + p_inner) += t40;
            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            *(temp + 8 * NPTS_LOCAL + p_inner) += t30;
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 21 * NPTS_LOCAL + p_inner) += t40;
            t20 = Y_PA * t10 - Y_PC * t11;
            t21 = Y_PA * t11 - Y_PC * t12;
            t22 = Y_PA * t12 - Y_PC * t13;
            *(temp + 1 * NPTS_LOCAL + p_inner) += t20;
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Y_PA * t21 - Y_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            *(temp + 9 * NPTS_LOCAL + p_inner) += t30;
            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 22 * NPTS_LOCAL + p_inner) += t40;
            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 23 * NPTS_LOCAL + p_inner) += t40;
            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            *(temp + 10 * NPTS_LOCAL + p_inner) += t30;
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 24 * NPTS_LOCAL + p_inner) += t40;
            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t22 = Z_PA * t12 - Z_PC * t13;
            *(temp + 2 * NPTS_LOCAL + p_inner) += t20;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            *(temp + 11 * NPTS_LOCAL + p_inner) += t30;
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 25 * NPTS_LOCAL + p_inner) += t40;
            t10 = Y_PA * t00 - Y_PC * t01;
            t11 = Y_PA * t01 - Y_PC * t02;
            t12 = Y_PA * t02 - Y_PC * t03;
            t13 = Y_PA * t03 - Y_PC * t04;
            t20 = Y_PA * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Y_PA * t11 - Y_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Y_PA * t12 - Y_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            *(temp + 3 * NPTS_LOCAL + p_inner) += t20;
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Y_PA * t21 - Y_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            *(temp + 12 * NPTS_LOCAL + p_inner) += t30;
            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 26 * NPTS_LOCAL + p_inner) += t40;
            t40 = Z_PA * t30 - Z_PC * t31;
            *(temp + 27 * NPTS_LOCAL + p_inner) += t40;
            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            *(temp + 13 * NPTS_LOCAL + p_inner) += t30;
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            *(temp + 28 * NPTS_LOCAL + p_inner) += t40;
            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t22 = Z_PA * t12 - Z_PC * t13;
            *(temp + 4 * NPTS_LOCAL + p_inner) += t20;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            *(temp + 14 * NPTS_LOCAL + p_inner) += t30;
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            *(temp + 29 * NPTS_LOCAL + p_inner) += t40;
            t10 = Z_PA * t00 - Z_PC * t01;
            t11 = Z_PA * t01 - Z_PC * t02;
            t12 = Z_PA * t02 - Z_PC * t03;
            t13 = Z_PA * t03 - Z_PC * t04;
            t20 = Z_PA * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Z_PA * t11 - Z_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Z_PA * t12 - Z_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            *(temp + 5 * NPTS_LOCAL + p_inner) += t20;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            *(temp + 15 * NPTS_LOCAL + p_inner) += t30;
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            *(temp + 30 * NPTS_LOCAL + p_inner) += t40;
         }
      }

      for(size_t p_inner = 0; p_inner < npts_inner; ++p_inner) {;
         double *Xik = (Xi + (p_outer + p_inner) * stX);
         double *Gik = (Gi + (p_outer + p_inner) * stG);

         *(Gik + 0 * ldG) += *(Xik + 0 * ldX) * (*(temp + 16 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 1 * ldG) += *(Xik + 0 * ldX) * (*(temp + 17 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 2 * ldG) += *(Xik + 0 * ldX) * (*(temp + 18 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 3 * ldG) += *(Xik + 0 * ldX) * (*(temp + 19 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 4 * ldG) += *(Xik + 0 * ldX) * (*(temp + 20 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 5 * ldG) += *(Xik + 0 * ldX) * (*(temp + 21 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 0 * ldG) += *(Xik + 1 * ldX) * (*(temp + 17 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 1 * ldG) += *(Xik + 1 * ldX) * (*(temp + 19 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 2 * ldG) += *(Xik + 1 * ldX) * (*(temp + 20 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 3 * ldG) += *(Xik + 1 * ldX) * (*(temp + 22 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 4 * ldG) += *(Xik + 1 * ldX) * (*(temp + 23 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 5 * ldG) += *(Xik + 1 * ldX) * (*(temp + 24 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 0 * ldG) += *(Xik + 2 * ldX) * (*(temp + 18 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 1 * ldG) += *(Xik + 2 * ldX) * (*(temp + 20 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 2 * ldG) += *(Xik + 2 * ldX) * (*(temp + 21 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 3 * ldG) += *(Xik + 2 * ldX) * (*(temp + 23 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 4 * ldG) += *(Xik + 2 * ldX) * (*(temp + 24 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 5 * ldG) += *(Xik + 2 * ldX) * (*(temp + 25 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 0 * ldG) += *(Xik + 3 * ldX) * (*(temp + 19 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 1 * ldG) += *(Xik + 3 * ldX) * (*(temp + 22 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 2 * ldG) += *(Xik + 3 * ldX) * (*(temp + 23 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 3 * ldG) += *(Xik + 3 * ldX) * (*(temp + 26 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 4 * ldG) += *(Xik + 3 * ldX) * (*(temp + 27 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 5 * ldG) += *(Xik + 3 * ldX) * (*(temp + 28 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 0 * ldG) += *(Xik + 4 * ldX) * (*(temp + 20 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 1 * ldG) += *(Xik + 4 * ldX) * (*(temp + 23 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 2 * ldG) += *(Xik + 4 * ldX) * (*(temp + 24 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 3 * ldG) += *(Xik + 4 * ldX) * (*(temp + 27 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 4 * ldG) += *(Xik + 4 * ldX) * (*(temp + 28 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 5 * ldG) += *(Xik + 4 * ldX) * (*(temp + 29 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 0 * ldG) += *(Xik + 5 * ldX) * (*(temp + 21 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 1 * ldG) += *(Xik + 5 * ldX) * (*(temp + 24 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 2 * ldG) += *(Xik + 5 * ldX) * (*(temp + 25 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 3 * ldG) += *(Xik + 5 * ldX) * (*(temp + 28 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 4 * ldG) += *(Xik + 5 * ldX) * (*(temp + 29 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 5 * ldG) += *(Xik + 5 * ldX) * (*(temp + 30 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
      }
   }
}
