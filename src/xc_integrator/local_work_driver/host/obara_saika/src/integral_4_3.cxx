#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

void integral_4_3(size_t npts,
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
   double temp[100 * NPTS_LOCAL];

   for(int i = 0; i < 100 * NPTS_LOCAL; ++i) {
      temp[i] = 0.0;
   }

   double X_AB = shpair.rAB.x;
   double Y_AB = shpair.rAB.y;
   double Z_AB = shpair.rAB.z;

   for(size_t p_outer = 0; p_outer < npts; p_outer += NPTS_LOCAL) {
      size_t npts_inner = MIN(NPTS_LOCAL, npts - p_outer);
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

         for(int p_inner = 0; p_inner < npts_inner; ++p_inner) {
            point C = *(_point_outer + p_inner);

            double xC = C.x;
            double yC = C.y;
            double zC = C.z;

            double X_PC = (xP - xC);
            double Y_PC = (yP - yC);
            double Z_PC = (zP - zC);

            double t00, t01, t02, t03, t04, t05, t06, t07, t10, t11, t12, t13, t14, t15, t16, t20, t21, t22, t23, t24, t25, t30, t31, t32, t33, t34, t40, t41, t42, t43, t50, t51, t52, t60, t61, t70;

            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

            t00 = eval * boys_function(0, tval);
            t01 = eval * boys_function(1, tval);
            t02 = eval * boys_function(2, tval);
            t03 = eval * boys_function(3, tval);
            t04 = eval * boys_function(4, tval);
            t05 = eval * boys_function(5, tval);
            t06 = eval * boys_function(6, tval);
            t07 = eval * boys_function(7, tval);
            t10 = X_PA * t00 - X_PC * t01;
            t11 = X_PA * t01 - X_PC * t02;
            t12 = X_PA * t02 - X_PC * t03;
            t13 = X_PA * t03 - X_PC * t04;
            t14 = X_PA * t04 - X_PC * t05;
            t15 = X_PA * t05 - X_PC * t06;
            t16 = X_PA * t06 - X_PC * t07;
            t20 = X_PA * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = X_PA * t11 - X_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = X_PA * t12 - X_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t23 = X_PA * t13 - X_PC * t14 + 0.5 * RHO_INV * 1 * (t03 - t04);
            t24 = X_PA * t14 - X_PC * t15 + 0.5 * RHO_INV * 1 * (t04 - t05);
            t25 = X_PA * t15 - X_PC * t16 + 0.5 * RHO_INV * 1 * (t05 - t06);
            t30 = X_PA * t20 - X_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = X_PA * t21 - X_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            t32 = X_PA * t22 - X_PC * t23 + 0.5 * RHO_INV * 2 * (t12 - t13);
            t33 = X_PA * t23 - X_PC * t24 + 0.5 * RHO_INV * 2 * (t13 - t14);
            t34 = X_PA * t24 - X_PC * t25 + 0.5 * RHO_INV * 2 * (t14 - t15);
            t40 = X_PA * t30 - X_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            t41 = X_PA * t31 - X_PC * t32 + 0.5 * RHO_INV * 3 * (t21 - t22);
            t42 = X_PA * t32 - X_PC * t33 + 0.5 * RHO_INV * 3 * (t22 - t23);
            t43 = X_PA * t33 - X_PC * t34 + 0.5 * RHO_INV * 3 * (t23 - t24);
            *(temp + 0 * NPTS_LOCAL + p_inner) += t40;
            t50 = X_PA * t40 - X_PC * t41 + 0.5 * RHO_INV * 4 * (t30 - t31);
            t51 = X_PA * t41 - X_PC * t42 + 0.5 * RHO_INV * 4 * (t31 - t32);
            t52 = X_PA * t42 - X_PC * t43 + 0.5 * RHO_INV * 4 * (t32 - t33);
            *(temp + 15 * NPTS_LOCAL + p_inner) += t50;
            t60 = X_PA * t50 - X_PC * t51 + 0.5 * RHO_INV * 5 * (t40 - t41);
            t61 = X_PA * t51 - X_PC * t52 + 0.5 * RHO_INV * 5 * (t41 - t42);
            *(temp + 36 * NPTS_LOCAL + p_inner) += t60;
            t70 = X_PA * t60 - X_PC * t61 + 0.5 * RHO_INV * 6 * (t50 - t51);
            *(temp + 64 * NPTS_LOCAL + p_inner) += t70;
            t70 = Y_PA * t60 - Y_PC * t61;
            *(temp + 65 * NPTS_LOCAL + p_inner) += t70;
            t70 = Z_PA * t60 - Z_PC * t61;
            *(temp + 66 * NPTS_LOCAL + p_inner) += t70;
            t60 = Y_PA * t50 - Y_PC * t51;
            t61 = Y_PA * t51 - Y_PC * t52;
            *(temp + 37 * NPTS_LOCAL + p_inner) += t60;
            t70 = Y_PA * t60 - Y_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            *(temp + 67 * NPTS_LOCAL + p_inner) += t70;
            t70 = Z_PA * t60 - Z_PC * t61;
            *(temp + 68 * NPTS_LOCAL + p_inner) += t70;
            t60 = Z_PA * t50 - Z_PC * t51;
            t61 = Z_PA * t51 - Z_PC * t52;
            *(temp + 38 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            *(temp + 69 * NPTS_LOCAL + p_inner) += t70;
            t50 = Y_PA * t40 - Y_PC * t41;
            t51 = Y_PA * t41 - Y_PC * t42;
            t52 = Y_PA * t42 - Y_PC * t43;
            *(temp + 16 * NPTS_LOCAL + p_inner) += t50;
            t60 = Y_PA * t50 - Y_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            t61 = Y_PA * t51 - Y_PC * t52 + 0.5 * RHO_INV * 1 * (t41 - t42);
            *(temp + 39 * NPTS_LOCAL + p_inner) += t60;
            t70 = Y_PA * t60 - Y_PC * t61 + 0.5 * RHO_INV * 2 * (t50 - t51);
            *(temp + 70 * NPTS_LOCAL + p_inner) += t70;
            t70 = Z_PA * t60 - Z_PC * t61;
            *(temp + 71 * NPTS_LOCAL + p_inner) += t70;
            t60 = Z_PA * t50 - Z_PC * t51;
            t61 = Z_PA * t51 - Z_PC * t52;
            *(temp + 40 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            *(temp + 72 * NPTS_LOCAL + p_inner) += t70;
            t50 = Z_PA * t40 - Z_PC * t41;
            t51 = Z_PA * t41 - Z_PC * t42;
            t52 = Z_PA * t42 - Z_PC * t43;
            *(temp + 17 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 1 * (t41 - t42);
            *(temp + 41 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 2 * (t50 - t51);
            *(temp + 73 * NPTS_LOCAL + p_inner) += t70;
            t40 = Y_PA * t30 - Y_PC * t31;
            t41 = Y_PA * t31 - Y_PC * t32;
            t42 = Y_PA * t32 - Y_PC * t33;
            t43 = Y_PA * t33 - Y_PC * t34;
            *(temp + 1 * NPTS_LOCAL + p_inner) += t40;
            t50 = Y_PA * t40 - Y_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Y_PA * t41 - Y_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            t52 = Y_PA * t42 - Y_PC * t43 + 0.5 * RHO_INV * 1 * (t32 - t33);
            *(temp + 18 * NPTS_LOCAL + p_inner) += t50;
            t60 = Y_PA * t50 - Y_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            t61 = Y_PA * t51 - Y_PC * t52 + 0.5 * RHO_INV * 2 * (t41 - t42);
            *(temp + 42 * NPTS_LOCAL + p_inner) += t60;
            t70 = Y_PA * t60 - Y_PC * t61 + 0.5 * RHO_INV * 3 * (t50 - t51);
            *(temp + 74 * NPTS_LOCAL + p_inner) += t70;
            t70 = Z_PA * t60 - Z_PC * t61;
            *(temp + 75 * NPTS_LOCAL + p_inner) += t70;
            t60 = Z_PA * t50 - Z_PC * t51;
            t61 = Z_PA * t51 - Z_PC * t52;
            *(temp + 43 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            *(temp + 76 * NPTS_LOCAL + p_inner) += t70;
            t50 = Z_PA * t40 - Z_PC * t41;
            t51 = Z_PA * t41 - Z_PC * t42;
            t52 = Z_PA * t42 - Z_PC * t43;
            *(temp + 19 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 1 * (t41 - t42);
            *(temp + 44 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 2 * (t50 - t51);
            *(temp + 77 * NPTS_LOCAL + p_inner) += t70;
            t40 = Z_PA * t30 - Z_PC * t31;
            t41 = Z_PA * t31 - Z_PC * t32;
            t42 = Z_PA * t32 - Z_PC * t33;
            t43 = Z_PA * t33 - Z_PC * t34;
            *(temp + 2 * NPTS_LOCAL + p_inner) += t40;
            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            t52 = Z_PA * t42 - Z_PC * t43 + 0.5 * RHO_INV * 1 * (t32 - t33);
            *(temp + 20 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 2 * (t41 - t42);
            *(temp + 45 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 3 * (t50 - t51);
            *(temp + 78 * NPTS_LOCAL + p_inner) += t70;
            t30 = Y_PA * t20 - Y_PC * t21;
            t31 = Y_PA * t21 - Y_PC * t22;
            t32 = Y_PA * t22 - Y_PC * t23;
            t33 = Y_PA * t23 - Y_PC * t24;
            t34 = Y_PA * t24 - Y_PC * t25;
            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            t41 = Y_PA * t31 - Y_PC * t32 + 0.5 * RHO_INV * 1 * (t21 - t22);
            t42 = Y_PA * t32 - Y_PC * t33 + 0.5 * RHO_INV * 1 * (t22 - t23);
            t43 = Y_PA * t33 - Y_PC * t34 + 0.5 * RHO_INV * 1 * (t23 - t24);
            *(temp + 3 * NPTS_LOCAL + p_inner) += t40;
            t50 = Y_PA * t40 - Y_PC * t41 + 0.5 * RHO_INV * 2 * (t30 - t31);
            t51 = Y_PA * t41 - Y_PC * t42 + 0.5 * RHO_INV * 2 * (t31 - t32);
            t52 = Y_PA * t42 - Y_PC * t43 + 0.5 * RHO_INV * 2 * (t32 - t33);
            *(temp + 21 * NPTS_LOCAL + p_inner) += t50;
            t60 = Y_PA * t50 - Y_PC * t51 + 0.5 * RHO_INV * 3 * (t40 - t41);
            t61 = Y_PA * t51 - Y_PC * t52 + 0.5 * RHO_INV * 3 * (t41 - t42);
            *(temp + 46 * NPTS_LOCAL + p_inner) += t60;
            t70 = Y_PA * t60 - Y_PC * t61 + 0.5 * RHO_INV * 4 * (t50 - t51);
            *(temp + 79 * NPTS_LOCAL + p_inner) += t70;
            t70 = Z_PA * t60 - Z_PC * t61;
            *(temp + 80 * NPTS_LOCAL + p_inner) += t70;
            t60 = Z_PA * t50 - Z_PC * t51;
            t61 = Z_PA * t51 - Z_PC * t52;
            *(temp + 47 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            *(temp + 81 * NPTS_LOCAL + p_inner) += t70;
            t50 = Z_PA * t40 - Z_PC * t41;
            t51 = Z_PA * t41 - Z_PC * t42;
            t52 = Z_PA * t42 - Z_PC * t43;
            *(temp + 22 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 1 * (t41 - t42);
            *(temp + 48 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 2 * (t50 - t51);
            *(temp + 82 * NPTS_LOCAL + p_inner) += t70;
            t40 = Z_PA * t30 - Z_PC * t31;
            t41 = Z_PA * t31 - Z_PC * t32;
            t42 = Z_PA * t32 - Z_PC * t33;
            t43 = Z_PA * t33 - Z_PC * t34;
            *(temp + 4 * NPTS_LOCAL + p_inner) += t40;
            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            t52 = Z_PA * t42 - Z_PC * t43 + 0.5 * RHO_INV * 1 * (t32 - t33);
            *(temp + 23 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 2 * (t41 - t42);
            *(temp + 49 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 3 * (t50 - t51);
            *(temp + 83 * NPTS_LOCAL + p_inner) += t70;
            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            t32 = Z_PA * t22 - Z_PC * t23;
            t33 = Z_PA * t23 - Z_PC * t24;
            t34 = Z_PA * t24 - Z_PC * t25;
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            t41 = Z_PA * t31 - Z_PC * t32 + 0.5 * RHO_INV * 1 * (t21 - t22);
            t42 = Z_PA * t32 - Z_PC * t33 + 0.5 * RHO_INV * 1 * (t22 - t23);
            t43 = Z_PA * t33 - Z_PC * t34 + 0.5 * RHO_INV * 1 * (t23 - t24);
            *(temp + 5 * NPTS_LOCAL + p_inner) += t40;
            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 2 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 2 * (t31 - t32);
            t52 = Z_PA * t42 - Z_PC * t43 + 0.5 * RHO_INV * 2 * (t32 - t33);
            *(temp + 24 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 3 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 3 * (t41 - t42);
            *(temp + 50 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 4 * (t50 - t51);
            *(temp + 84 * NPTS_LOCAL + p_inner) += t70;
            t20 = Y_PA * t10 - Y_PC * t11;
            t21 = Y_PA * t11 - Y_PC * t12;
            t22 = Y_PA * t12 - Y_PC * t13;
            t23 = Y_PA * t13 - Y_PC * t14;
            t24 = Y_PA * t14 - Y_PC * t15;
            t25 = Y_PA * t15 - Y_PC * t16;
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Y_PA * t21 - Y_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            t32 = Y_PA * t22 - Y_PC * t23 + 0.5 * RHO_INV * 1 * (t12 - t13);
            t33 = Y_PA * t23 - Y_PC * t24 + 0.5 * RHO_INV * 1 * (t13 - t14);
            t34 = Y_PA * t24 - Y_PC * t25 + 0.5 * RHO_INV * 1 * (t14 - t15);
            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            t41 = Y_PA * t31 - Y_PC * t32 + 0.5 * RHO_INV * 2 * (t21 - t22);
            t42 = Y_PA * t32 - Y_PC * t33 + 0.5 * RHO_INV * 2 * (t22 - t23);
            t43 = Y_PA * t33 - Y_PC * t34 + 0.5 * RHO_INV * 2 * (t23 - t24);
            *(temp + 6 * NPTS_LOCAL + p_inner) += t40;
            t50 = Y_PA * t40 - Y_PC * t41 + 0.5 * RHO_INV * 3 * (t30 - t31);
            t51 = Y_PA * t41 - Y_PC * t42 + 0.5 * RHO_INV * 3 * (t31 - t32);
            t52 = Y_PA * t42 - Y_PC * t43 + 0.5 * RHO_INV * 3 * (t32 - t33);
            *(temp + 25 * NPTS_LOCAL + p_inner) += t50;
            t60 = Y_PA * t50 - Y_PC * t51 + 0.5 * RHO_INV * 4 * (t40 - t41);
            t61 = Y_PA * t51 - Y_PC * t52 + 0.5 * RHO_INV * 4 * (t41 - t42);
            *(temp + 51 * NPTS_LOCAL + p_inner) += t60;
            t70 = Y_PA * t60 - Y_PC * t61 + 0.5 * RHO_INV * 5 * (t50 - t51);
            *(temp + 85 * NPTS_LOCAL + p_inner) += t70;
            t70 = Z_PA * t60 - Z_PC * t61;
            *(temp + 86 * NPTS_LOCAL + p_inner) += t70;
            t60 = Z_PA * t50 - Z_PC * t51;
            t61 = Z_PA * t51 - Z_PC * t52;
            *(temp + 52 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            *(temp + 87 * NPTS_LOCAL + p_inner) += t70;
            t50 = Z_PA * t40 - Z_PC * t41;
            t51 = Z_PA * t41 - Z_PC * t42;
            t52 = Z_PA * t42 - Z_PC * t43;
            *(temp + 26 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 1 * (t41 - t42);
            *(temp + 53 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 2 * (t50 - t51);
            *(temp + 88 * NPTS_LOCAL + p_inner) += t70;
            t40 = Z_PA * t30 - Z_PC * t31;
            t41 = Z_PA * t31 - Z_PC * t32;
            t42 = Z_PA * t32 - Z_PC * t33;
            t43 = Z_PA * t33 - Z_PC * t34;
            *(temp + 7 * NPTS_LOCAL + p_inner) += t40;
            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            t52 = Z_PA * t42 - Z_PC * t43 + 0.5 * RHO_INV * 1 * (t32 - t33);
            *(temp + 27 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 2 * (t41 - t42);
            *(temp + 54 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 3 * (t50 - t51);
            *(temp + 89 * NPTS_LOCAL + p_inner) += t70;
            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            t32 = Z_PA * t22 - Z_PC * t23;
            t33 = Z_PA * t23 - Z_PC * t24;
            t34 = Z_PA * t24 - Z_PC * t25;
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            t41 = Z_PA * t31 - Z_PC * t32 + 0.5 * RHO_INV * 1 * (t21 - t22);
            t42 = Z_PA * t32 - Z_PC * t33 + 0.5 * RHO_INV * 1 * (t22 - t23);
            t43 = Z_PA * t33 - Z_PC * t34 + 0.5 * RHO_INV * 1 * (t23 - t24);
            *(temp + 8 * NPTS_LOCAL + p_inner) += t40;
            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 2 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 2 * (t31 - t32);
            t52 = Z_PA * t42 - Z_PC * t43 + 0.5 * RHO_INV * 2 * (t32 - t33);
            *(temp + 28 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 3 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 3 * (t41 - t42);
            *(temp + 55 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 4 * (t50 - t51);
            *(temp + 90 * NPTS_LOCAL + p_inner) += t70;
            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t22 = Z_PA * t12 - Z_PC * t13;
            t23 = Z_PA * t13 - Z_PC * t14;
            t24 = Z_PA * t14 - Z_PC * t15;
            t25 = Z_PA * t15 - Z_PC * t16;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            t32 = Z_PA * t22 - Z_PC * t23 + 0.5 * RHO_INV * 1 * (t12 - t13);
            t33 = Z_PA * t23 - Z_PC * t24 + 0.5 * RHO_INV * 1 * (t13 - t14);
            t34 = Z_PA * t24 - Z_PC * t25 + 0.5 * RHO_INV * 1 * (t14 - t15);
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            t41 = Z_PA * t31 - Z_PC * t32 + 0.5 * RHO_INV * 2 * (t21 - t22);
            t42 = Z_PA * t32 - Z_PC * t33 + 0.5 * RHO_INV * 2 * (t22 - t23);
            t43 = Z_PA * t33 - Z_PC * t34 + 0.5 * RHO_INV * 2 * (t23 - t24);
            *(temp + 9 * NPTS_LOCAL + p_inner) += t40;
            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 3 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 3 * (t31 - t32);
            t52 = Z_PA * t42 - Z_PC * t43 + 0.5 * RHO_INV * 3 * (t32 - t33);
            *(temp + 29 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 4 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 4 * (t41 - t42);
            *(temp + 56 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 5 * (t50 - t51);
            *(temp + 91 * NPTS_LOCAL + p_inner) += t70;
            t10 = Y_PA * t00 - Y_PC * t01;
            t11 = Y_PA * t01 - Y_PC * t02;
            t12 = Y_PA * t02 - Y_PC * t03;
            t13 = Y_PA * t03 - Y_PC * t04;
            t14 = Y_PA * t04 - Y_PC * t05;
            t15 = Y_PA * t05 - Y_PC * t06;
            t16 = Y_PA * t06 - Y_PC * t07;
            t20 = Y_PA * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Y_PA * t11 - Y_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Y_PA * t12 - Y_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t23 = Y_PA * t13 - Y_PC * t14 + 0.5 * RHO_INV * 1 * (t03 - t04);
            t24 = Y_PA * t14 - Y_PC * t15 + 0.5 * RHO_INV * 1 * (t04 - t05);
            t25 = Y_PA * t15 - Y_PC * t16 + 0.5 * RHO_INV * 1 * (t05 - t06);
            t30 = Y_PA * t20 - Y_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Y_PA * t21 - Y_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            t32 = Y_PA * t22 - Y_PC * t23 + 0.5 * RHO_INV * 2 * (t12 - t13);
            t33 = Y_PA * t23 - Y_PC * t24 + 0.5 * RHO_INV * 2 * (t13 - t14);
            t34 = Y_PA * t24 - Y_PC * t25 + 0.5 * RHO_INV * 2 * (t14 - t15);
            t40 = Y_PA * t30 - Y_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            t41 = Y_PA * t31 - Y_PC * t32 + 0.5 * RHO_INV * 3 * (t21 - t22);
            t42 = Y_PA * t32 - Y_PC * t33 + 0.5 * RHO_INV * 3 * (t22 - t23);
            t43 = Y_PA * t33 - Y_PC * t34 + 0.5 * RHO_INV * 3 * (t23 - t24);
            *(temp + 10 * NPTS_LOCAL + p_inner) += t40;
            t50 = Y_PA * t40 - Y_PC * t41 + 0.5 * RHO_INV * 4 * (t30 - t31);
            t51 = Y_PA * t41 - Y_PC * t42 + 0.5 * RHO_INV * 4 * (t31 - t32);
            t52 = Y_PA * t42 - Y_PC * t43 + 0.5 * RHO_INV * 4 * (t32 - t33);
            *(temp + 30 * NPTS_LOCAL + p_inner) += t50;
            t60 = Y_PA * t50 - Y_PC * t51 + 0.5 * RHO_INV * 5 * (t40 - t41);
            t61 = Y_PA * t51 - Y_PC * t52 + 0.5 * RHO_INV * 5 * (t41 - t42);
            *(temp + 57 * NPTS_LOCAL + p_inner) += t60;
            t70 = Y_PA * t60 - Y_PC * t61 + 0.5 * RHO_INV * 6 * (t50 - t51);
            *(temp + 92 * NPTS_LOCAL + p_inner) += t70;
            t70 = Z_PA * t60 - Z_PC * t61;
            *(temp + 93 * NPTS_LOCAL + p_inner) += t70;
            t60 = Z_PA * t50 - Z_PC * t51;
            t61 = Z_PA * t51 - Z_PC * t52;
            *(temp + 58 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 1 * (t50 - t51);
            *(temp + 94 * NPTS_LOCAL + p_inner) += t70;
            t50 = Z_PA * t40 - Z_PC * t41;
            t51 = Z_PA * t41 - Z_PC * t42;
            t52 = Z_PA * t42 - Z_PC * t43;
            *(temp + 31 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 1 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 1 * (t41 - t42);
            *(temp + 59 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 2 * (t50 - t51);
            *(temp + 95 * NPTS_LOCAL + p_inner) += t70;
            t40 = Z_PA * t30 - Z_PC * t31;
            t41 = Z_PA * t31 - Z_PC * t32;
            t42 = Z_PA * t32 - Z_PC * t33;
            t43 = Z_PA * t33 - Z_PC * t34;
            *(temp + 11 * NPTS_LOCAL + p_inner) += t40;
            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 1 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 1 * (t31 - t32);
            t52 = Z_PA * t42 - Z_PC * t43 + 0.5 * RHO_INV * 1 * (t32 - t33);
            *(temp + 32 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 2 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 2 * (t41 - t42);
            *(temp + 60 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 3 * (t50 - t51);
            *(temp + 96 * NPTS_LOCAL + p_inner) += t70;
            t30 = Z_PA * t20 - Z_PC * t21;
            t31 = Z_PA * t21 - Z_PC * t22;
            t32 = Z_PA * t22 - Z_PC * t23;
            t33 = Z_PA * t23 - Z_PC * t24;
            t34 = Z_PA * t24 - Z_PC * t25;
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 1 * (t20 - t21);
            t41 = Z_PA * t31 - Z_PC * t32 + 0.5 * RHO_INV * 1 * (t21 - t22);
            t42 = Z_PA * t32 - Z_PC * t33 + 0.5 * RHO_INV * 1 * (t22 - t23);
            t43 = Z_PA * t33 - Z_PC * t34 + 0.5 * RHO_INV * 1 * (t23 - t24);
            *(temp + 12 * NPTS_LOCAL + p_inner) += t40;
            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 2 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 2 * (t31 - t32);
            t52 = Z_PA * t42 - Z_PC * t43 + 0.5 * RHO_INV * 2 * (t32 - t33);
            *(temp + 33 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 3 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 3 * (t41 - t42);
            *(temp + 61 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 4 * (t50 - t51);
            *(temp + 97 * NPTS_LOCAL + p_inner) += t70;
            t20 = Z_PA * t10 - Z_PC * t11;
            t21 = Z_PA * t11 - Z_PC * t12;
            t22 = Z_PA * t12 - Z_PC * t13;
            t23 = Z_PA * t13 - Z_PC * t14;
            t24 = Z_PA * t14 - Z_PC * t15;
            t25 = Z_PA * t15 - Z_PC * t16;
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 1 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 1 * (t11 - t12);
            t32 = Z_PA * t22 - Z_PC * t23 + 0.5 * RHO_INV * 1 * (t12 - t13);
            t33 = Z_PA * t23 - Z_PC * t24 + 0.5 * RHO_INV * 1 * (t13 - t14);
            t34 = Z_PA * t24 - Z_PC * t25 + 0.5 * RHO_INV * 1 * (t14 - t15);
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 2 * (t20 - t21);
            t41 = Z_PA * t31 - Z_PC * t32 + 0.5 * RHO_INV * 2 * (t21 - t22);
            t42 = Z_PA * t32 - Z_PC * t33 + 0.5 * RHO_INV * 2 * (t22 - t23);
            t43 = Z_PA * t33 - Z_PC * t34 + 0.5 * RHO_INV * 2 * (t23 - t24);
            *(temp + 13 * NPTS_LOCAL + p_inner) += t40;
            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 3 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 3 * (t31 - t32);
            t52 = Z_PA * t42 - Z_PC * t43 + 0.5 * RHO_INV * 3 * (t32 - t33);
            *(temp + 34 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 4 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 4 * (t41 - t42);
            *(temp + 62 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 5 * (t50 - t51);
            *(temp + 98 * NPTS_LOCAL + p_inner) += t70;
            t10 = Z_PA * t00 - Z_PC * t01;
            t11 = Z_PA * t01 - Z_PC * t02;
            t12 = Z_PA * t02 - Z_PC * t03;
            t13 = Z_PA * t03 - Z_PC * t04;
            t14 = Z_PA * t04 - Z_PC * t05;
            t15 = Z_PA * t05 - Z_PC * t06;
            t16 = Z_PA * t06 - Z_PC * t07;
            t20 = Z_PA * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            t21 = Z_PA * t11 - Z_PC * t12 + 0.5 * RHO_INV * 1 * (t01 - t02);
            t22 = Z_PA * t12 - Z_PC * t13 + 0.5 * RHO_INV * 1 * (t02 - t03);
            t23 = Z_PA * t13 - Z_PC * t14 + 0.5 * RHO_INV * 1 * (t03 - t04);
            t24 = Z_PA * t14 - Z_PC * t15 + 0.5 * RHO_INV * 1 * (t04 - t05);
            t25 = Z_PA * t15 - Z_PC * t16 + 0.5 * RHO_INV * 1 * (t05 - t06);
            t30 = Z_PA * t20 - Z_PC * t21 + 0.5 * RHO_INV * 2 * (t10 - t11);
            t31 = Z_PA * t21 - Z_PC * t22 + 0.5 * RHO_INV * 2 * (t11 - t12);
            t32 = Z_PA * t22 - Z_PC * t23 + 0.5 * RHO_INV * 2 * (t12 - t13);
            t33 = Z_PA * t23 - Z_PC * t24 + 0.5 * RHO_INV * 2 * (t13 - t14);
            t34 = Z_PA * t24 - Z_PC * t25 + 0.5 * RHO_INV * 2 * (t14 - t15);
            t40 = Z_PA * t30 - Z_PC * t31 + 0.5 * RHO_INV * 3 * (t20 - t21);
            t41 = Z_PA * t31 - Z_PC * t32 + 0.5 * RHO_INV * 3 * (t21 - t22);
            t42 = Z_PA * t32 - Z_PC * t33 + 0.5 * RHO_INV * 3 * (t22 - t23);
            t43 = Z_PA * t33 - Z_PC * t34 + 0.5 * RHO_INV * 3 * (t23 - t24);
            *(temp + 14 * NPTS_LOCAL + p_inner) += t40;
            t50 = Z_PA * t40 - Z_PC * t41 + 0.5 * RHO_INV * 4 * (t30 - t31);
            t51 = Z_PA * t41 - Z_PC * t42 + 0.5 * RHO_INV * 4 * (t31 - t32);
            t52 = Z_PA * t42 - Z_PC * t43 + 0.5 * RHO_INV * 4 * (t32 - t33);
            *(temp + 35 * NPTS_LOCAL + p_inner) += t50;
            t60 = Z_PA * t50 - Z_PC * t51 + 0.5 * RHO_INV * 5 * (t40 - t41);
            t61 = Z_PA * t51 - Z_PC * t52 + 0.5 * RHO_INV * 5 * (t41 - t42);
            *(temp + 63 * NPTS_LOCAL + p_inner) += t60;
            t70 = Z_PA * t60 - Z_PC * t61 + 0.5 * RHO_INV * 6 * (t50 - t51);
            *(temp + 99 * NPTS_LOCAL + p_inner) += t70;
         }
      }

      for(int p_inner = 0; p_inner < npts_inner; ++p_inner) {
         double *Xik = (Xi + (NPTS_LOCAL * p_outer + p_inner) * stX);
         double *Xjk = (Xj + (NPTS_LOCAL * p_outer + p_inner) * stX);
         double *Gik = (Gi + (NPTS_LOCAL * p_outer + p_inner) * stG);
         double *Gjk = (Gj + (NPTS_LOCAL * p_outer + p_inner) * stG);

         for(int c0 = 0; c0 <= 3; ++c0) {
            for(int c1 = 0; c1 <= c0; ++c1) {
               int m = 3 - c0;
               int n = c0 - c1;
               int p = c1;

               int idxB = (((3 - m) * (3 - m + 1)) >> 1) + p;

               double X_ABp = 1.0, comb_m_i = 1.0;
               for(int i = 0; i <= m; ++i) {
                  double rcp_i;

                  double Y_ABp = 1.0, comb_n_j = 1.0;
                  for(int j = 0; j <= n; ++j) {
                     double rcp_j;

                     double Z_ABp = 1.0, comb_p_k = 1.0;
                     for(int k = 0; k <= p; ++k) {
                        double rcp_k;
                        int mv, pv, Lv = 7 - i - j - k;

                        int offset = (Lv * (Lv + 1) * (Lv + 2) - 120) / 6;
                        double const_value = *(weights + NPTS_LOCAL * p_outer + p_inner) * comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
                        mv = 4 + m - i; pv = 0 + p - k;
                        double t0 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 0 * ldG) += *(Xjk + idxB * ldX) * t0;
                        *(Gjk + idxB * ldG) += *(Xik + 0 * ldX) * t0;
                        mv = 3 + m - i; pv = 0 + p - k;
                        double t1 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 1 * ldG) += *(Xjk + idxB * ldX) * t1;
                        *(Gjk + idxB * ldG) += *(Xik + 1 * ldX) * t1;
                        mv = 3 + m - i; pv = 1 + p - k;
                        double t2 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 2 * ldG) += *(Xjk + idxB * ldX) * t2;
                        *(Gjk + idxB * ldG) += *(Xik + 2 * ldX) * t2;
                        mv = 2 + m - i; pv = 0 + p - k;
                        double t3 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 3 * ldG) += *(Xjk + idxB * ldX) * t3;
                        *(Gjk + idxB * ldG) += *(Xik + 3 * ldX) * t3;
                        mv = 2 + m - i; pv = 1 + p - k;
                        double t4 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 4 * ldG) += *(Xjk + idxB * ldX) * t4;
                        *(Gjk + idxB * ldG) += *(Xik + 4 * ldX) * t4;
                        mv = 2 + m - i; pv = 2 + p - k;
                        double t5 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 5 * ldG) += *(Xjk + idxB * ldX) * t5;
                        *(Gjk + idxB * ldG) += *(Xik + 5 * ldX) * t5;
                        mv = 1 + m - i; pv = 0 + p - k;
                        double t6 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 6 * ldG) += *(Xjk + idxB * ldX) * t6;
                        *(Gjk + idxB * ldG) += *(Xik + 6 * ldX) * t6;
                        mv = 1 + m - i; pv = 1 + p - k;
                        double t7 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 7 * ldG) += *(Xjk + idxB * ldX) * t7;
                        *(Gjk + idxB * ldG) += *(Xik + 7 * ldX) * t7;
                        mv = 1 + m - i; pv = 2 + p - k;
                        double t8 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 8 * ldG) += *(Xjk + idxB * ldX) * t8;
                        *(Gjk + idxB * ldG) += *(Xik + 8 * ldX) * t8;
                        mv = 1 + m - i; pv = 3 + p - k;
                        double t9 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 9 * ldG) += *(Xjk + idxB * ldX) * t9;
                        *(Gjk + idxB * ldG) += *(Xik + 9 * ldX) * t9;
                        mv = 0 + m - i; pv = 0 + p - k;
                        double t10 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 10 * ldG) += *(Xjk + idxB * ldX) * t10;
                        *(Gjk + idxB * ldG) += *(Xik + 10 * ldX) * t10;
                        mv = 0 + m - i; pv = 1 + p - k;
                        double t11 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 11 * ldG) += *(Xjk + idxB * ldX) * t11;
                        *(Gjk + idxB * ldG) += *(Xik + 11 * ldX) * t11;
                        mv = 0 + m - i; pv = 2 + p - k;
                        double t12 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 12 * ldG) += *(Xjk + idxB * ldX) * t12;
                        *(Gjk + idxB * ldG) += *(Xik + 12 * ldX) * t12;
                        mv = 0 + m - i; pv = 3 + p - k;
                        double t13 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 13 * ldG) += *(Xjk + idxB * ldX) * t13;
                        *(Gjk + idxB * ldG) += *(Xik + 13 * ldX) * t13;
                        mv = 0 + m - i; pv = 4 + p - k;
                        double t14 = *(temp + (offset + (((Lv - mv) * (Lv - mv + 1)) >> 1) + pv) * NPTS_LOCAL + p_inner) * const_value;
                        *(Gik + 14 * ldG) += *(Xjk + idxB * ldX) * t14;
                        *(Gjk + idxB * ldG) += *(Xik + 14 * ldX) * t14;

                        Z_ABp *= Z_AB; rcp_k = 1.0 / (1.0 * (k + 1)); comb_p_k = (comb_p_k * (p - k)) * rcp_k;
                     }

                     Y_ABp *= Y_AB; rcp_j = 1.0 / (1.0 * (j + 1)); comb_n_j = (comb_n_j * (n - j)) * rcp_j;
                  }

                  X_ABp *= X_AB; rcp_i = 1.0 / (1.0 * (i + 1)); comb_m_i = (comb_m_i * (m - i)) * rcp_i;
               }
            }
         }
      }
   }
}
