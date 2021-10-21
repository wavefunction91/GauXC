#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

void integral_1(size_t npts,
               shell_pair shpair,
               point *_points,
               double *Xi,
               int stX,
               int ldX,
               double *Gi,
               int stG, 
               int ldG, 
               double *weights) {
   double temp[9 * NPTS_LOCAL];

   for(size_t p_outer = 0; p_outer < npts; p_outer += NPTS_LOCAL) {
      size_t npts_inner = MIN(NPTS_LOCAL, npts - p_outer);
      point *_point_outer = (_points + p_outer);

      double xA = shpair.rA.x;
      double yA = shpair.rA.y;
      double zA = shpair.rA.z;

      for(int i = 0; i < 9 * NPTS_LOCAL; ++i) temp[i] = 0.0;

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

            double t00, t01, t02, t10, t11, t20;

            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

            t00 = eval * boys_function(0, tval);
            t01 = eval * boys_function(1, tval);
            t02 = eval * boys_function(2, tval);
            t10 = X_PA * t00 - X_PC * t01;
            t11 = X_PA * t01 - X_PC * t02;
            *(temp + 0 * NPTS_LOCAL + p_inner) += t10;
            t20 = X_PA * t10 - X_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            *(temp + 3 * NPTS_LOCAL + p_inner) += t20;
            t20 = Y_PA * t10 - Y_PC * t11;
            *(temp + 4 * NPTS_LOCAL + p_inner) += t20;
            t20 = Z_PA * t10 - Z_PC * t11;
            *(temp + 5 * NPTS_LOCAL + p_inner) += t20;
            t10 = Y_PA * t00 - Y_PC * t01;
            t11 = Y_PA * t01 - Y_PC * t02;
            *(temp + 1 * NPTS_LOCAL + p_inner) += t10;
            t20 = Y_PA * t10 - Y_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            *(temp + 6 * NPTS_LOCAL + p_inner) += t20;
            t20 = Z_PA * t10 - Z_PC * t11;
            *(temp + 7 * NPTS_LOCAL + p_inner) += t20;
            t10 = Z_PA * t00 - Z_PC * t01;
            t11 = Z_PA * t01 - Z_PC * t02;
            *(temp + 2 * NPTS_LOCAL + p_inner) += t10;
            t20 = Z_PA * t10 - Z_PC * t11 + 0.5 * RHO_INV * 1 * (t00 - t01);
            *(temp + 8 * NPTS_LOCAL + p_inner) += t20;
         }
      }

      for(size_t p_inner = 0; p_inner < npts_inner; ++p_inner) {;
         double *Xik = (Xi + (p_outer + p_inner) * stX);
         double *Gik = (Gi + (p_outer + p_inner) * stG);

         *(Gik + 0 * ldG) += *(Xik + 0 * ldX) * (*(temp + 3 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 1 * ldG) += *(Xik + 0 * ldX) * (*(temp + 4 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 2 * ldG) += *(Xik + 0 * ldX) * (*(temp + 5 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 0 * ldG) += *(Xik + 1 * ldX) * (*(temp + 4 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 1 * ldG) += *(Xik + 1 * ldX) * (*(temp + 6 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 2 * ldG) += *(Xik + 1 * ldX) * (*(temp + 7 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 0 * ldG) += *(Xik + 2 * ldX) * (*(temp + 5 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 1 * ldG) += *(Xik + 2 * ldX) * (*(temp + 7 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
         *(Gik + 2 * ldG) += *(Xik + 2 * ldX) * (*(temp + 8 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
      }
   }
}
