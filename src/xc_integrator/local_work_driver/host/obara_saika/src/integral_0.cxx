#include <math.h>
#include "boys_computation.h"
#include "integral_data_types.h"

#define PI 3.14159265358979323846

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

void integral_0(size_t npts,
               shell_pair shpair,
               point *_points,
               double *Xi,
               int stX,
               int ldX,
               double *Gi,
               int stG, 
               int ldG, 
               double *weights) {
   double temp[1 * NPTS_LOCAL];

   for(size_t p_outer = 0; p_outer < npts; p_outer += NPTS_LOCAL) {
      size_t npts_inner = MIN(NPTS_LOCAL, npts - p_outer);
      point *_point_outer = (_points + p_outer);

      double xA = shpair.rA.x;
      double yA = shpair.rA.y;
      double zA = shpair.rA.z;

      for(int i = 0; i < 1 * NPTS_LOCAL; ++i) temp[i] = 0.0;

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

            double t00;

            double tval = RHO * (X_PC * X_PC + Y_PC * Y_PC + Z_PC * Z_PC);

            t00 = eval * boys_function(0, tval);
            *(temp + 0 * NPTS_LOCAL + p_inner) += t00;
         }
      }

      for(size_t p_inner = 0; p_inner < npts_inner; ++p_inner) {;
         double *Xik = (Xi + (p_outer + p_inner) * stX);
         double *Gik = (Gi + (p_outer + p_inner) * stG);

         *(Gik + 0 * ldG) += *(Xik + 0 * ldX) * (*(temp + 0 * NPTS_LOCAL + p_inner)) * (*(weights + p_outer + p_inner));
      }
   }
}
