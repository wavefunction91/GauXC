#include <math.h>
#include "chebyshev_boys_computation.hpp"
#include "integral_data_types.hpp"
#include "config_obara_saika.hpp"
#include "integral_0.hu"

#define PI 3.14159265358979323846

#define MIN(a,b)			\
  ({ __typeof__ (a) _a = (a);	        \
  __typeof__ (b) _b = (b);		\
  _a < _b ? _a : _b; })

namespace XGPU {
__global__ void integral_0(size_t npts,
                          point rA,
                          point rB,
                          int nprim_pairs,
                          prim_pair *prim_pairs,
                          double *_points,
                          double *Xi,
                          int ldX,
                          double *Gi,
                          int ldG, 
                          double *weights,
                          double *boys_table) {
   __shared__ double *temp;
   for(size_t p_outer = blockIdx.x * blockDim.x; p_outer < npts; p_outer += gridDim.x * blockDim.x) {
      double *_point_outer = (_points + p_outer);

      size_t p_inner = (threadIdx.x < (npts - p_outer)) ? threadIdx.x : (npts - p_outer);

      double xA = rA.x;
      double yA = rA.y;
      double zA = rA.z;

      for(int i = 0; i < 1; ++i) SCALAR_STORE((temp + i * blockDim.x + threadIdx.x), SCALAR_ZERO());

      for(int ij = 0; ij < nprim_pairs; ++ij) {
         double RHO = prim_pairs[ij].gamma;

         double eval = prim_pairs[ij].K_coeff_prod;

         // Evaluate T Values
         SCALAR_TYPE xC = SCALAR_LOAD((_point_outer + p_inner + 0 * npts));
         SCALAR_TYPE yC = SCALAR_LOAD((_point_outer + p_inner + 1 * npts));
         SCALAR_TYPE zC = SCALAR_LOAD((_point_outer + p_inner + 2 * npts));

         SCALAR_TYPE X_PC = SCALAR_SUB(xA, xC);
         SCALAR_TYPE Y_PC = SCALAR_SUB(yA, yC);
         SCALAR_TYPE Z_PC = SCALAR_SUB(zA, zC);

         X_PC = SCALAR_MUL(X_PC, X_PC);
         X_PC = SCALAR_FMA(Y_PC, Y_PC, X_PC);
         X_PC = SCALAR_FMA(Z_PC, Z_PC, X_PC);
         SCALAR_TYPE TVAL = SCALAR_MUL(RHO, X_PC);

         SCALAR_TYPE t00, TVAL_inv_e;

         // Evaluate Boys function
         boys_element<0>(&TVAL, &TVAL_inv_e, &t00, boys_table);

         // Evaluate VRR Buffer
         SCALAR_TYPE tx;


         t00 = SCALAR_MUL(eval, t00);
         tx = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
         tx = SCALAR_ADD(tx, t00);
         SCALAR_STORE((temp + 0 * blockDim.x + threadIdx.x), tx);
      }

      double *Xik = (Xi + p_outer + p_inner);
      double *Gik = (Gi + p_outer + p_inner);

      SCALAR_TYPE tx, wg, xik, gik;
      tx  = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
      wg  = SCALAR_LOAD((weights + p_outer + p_inner));

      xik = SCALAR_LOAD((Xik + 0 * ldX));
      gik = SCALAR_LOAD((Gik + 0 * ldG));

      tx = SCALAR_MUL(tx, wg);
      gik = SCALAR_FMA(tx, xik, gik);
      SCALAR_STORE((Gik + 0 * ldG), gik);
   }
}
}
