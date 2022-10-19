#include <math.h>
#include "../include/gpu/chebyshev_boys_computation.hpp"
#include "config_obara_saika.hpp"
#include "integral_1_1.hu"

namespace XGPU {
  __inline__ __device__ void dev_integral_1_1_driver(double X_AB,
				   double Y_AB,
				   double Z_AB,
				   size_t npts,
				   double *_points_x,
				   double *_points_y,
				   double *_points_z,
           shell_pair* sp,
				   double *Xi,
				   double *Xj,
				   int ldX,
				   double *Gi,
				   double *Gj,
				   int ldG, 
				   double *weights, 
				   double *boys_table) {
    //__shared__ double temp[128 * 9];
    double temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8;
    const auto nprim_pairs = sp->nprim_pairs();
    const auto prim_pairs  = sp->prim_pairs();
    
    for(size_t p_outer = blockIdx.x * blockDim.x; p_outer < npts; p_outer += gridDim.x * blockDim.x) {
      double *_point_outer_x = (_points_x + p_outer);
      double *_point_outer_y = (_points_y + p_outer);
      double *_point_outer_z = (_points_z + p_outer);

      size_t p_inner = (threadIdx.x < (npts - p_outer)) ? threadIdx.x : (npts - p_outer);

      //for(int i = 0; i < 9; ++i) SCALAR_STORE((temp + i * blockDim.x + threadIdx.x), SCALAR_ZERO());
      temp_0 = SCALAR_ZERO();
      temp_1 = SCALAR_ZERO();
      temp_2 = SCALAR_ZERO();
      temp_3 = SCALAR_ZERO();
      temp_4 = SCALAR_ZERO();
      temp_5 = SCALAR_ZERO();
      temp_6 = SCALAR_ZERO();
      temp_7 = SCALAR_ZERO();
      temp_8 = SCALAR_ZERO();


      for(int ij = 0; ij < nprim_pairs; ++ij) {
	double RHO = prim_pairs[ij].gamma;
	double RHO_INV = prim_pairs[ij].gamma_inv;
	double X_PA = prim_pairs[ij].PA.x;
	double Y_PA = prim_pairs[ij].PA.y;
	double Z_PA = prim_pairs[ij].PA.z;

	double xP = prim_pairs[ij].P.x;
	double yP = prim_pairs[ij].P.y;
	double zP = prim_pairs[ij].P.z;

	double eval = prim_pairs[ij].K_coeff_prod;
    if( std::abs(eval) < 1e-12 ) continue;

	// Evaluate T Values
	SCALAR_TYPE xC = SCALAR_LOAD((_point_outer_x + p_inner));
	SCALAR_TYPE yC = SCALAR_LOAD((_point_outer_y + p_inner));
	SCALAR_TYPE zC = SCALAR_LOAD((_point_outer_z + p_inner));

	SCALAR_TYPE X_PC = SCALAR_SUB(xP, xC);
	SCALAR_TYPE Y_PC = SCALAR_SUB(yP, yC);
	SCALAR_TYPE Z_PC = SCALAR_SUB(zP, zC);

	SCALAR_TYPE TVAL = SCALAR_MUL(X_PC, X_PC);
	TVAL = SCALAR_FMA(Y_PC, Y_PC, TVAL);
	TVAL = SCALAR_FMA(Z_PC, Z_PC, TVAL);
	TVAL = SCALAR_MUL(RHO, TVAL);

	SCALAR_TYPE t00, t01, t02, TVAL_inv_e;

	// Evaluate Boys function
	boys_element<2>(&TVAL, &TVAL_inv_e, &t02, boys_table);

	// Evaluate VRR Buffer
	SCALAR_TYPE t10, t11, t20, tx, ty;

	t01 = SCALAR_MUL(SCALAR_ADD(SCALAR_MUL(TVAL, t02), TVAL_inv_e), SCALAR_SET1(0.66666666666666662966));
	t00 = SCALAR_MUL(SCALAR_ADD(SCALAR_MUL(TVAL, t01), TVAL_inv_e), SCALAR_SET1(2.00000000000000000000));

	t00 = SCALAR_MUL(eval, t00);
	t01 = SCALAR_MUL(eval, t01);
	t02 = SCALAR_MUL(eval, t02);
	t10 = SCALAR_MUL(X_PA, t00);
	t10 = SCALAR_FNMA(X_PC, t01, t10);
	t11 = SCALAR_MUL(X_PA, t01);
	t11 = SCALAR_FNMA(X_PC, t02, t11);
	//tx = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
  tx = temp_0;
	tx = SCALAR_ADD(tx, t10);
	//SCALAR_STORE((temp + 0 * blockDim.x + threadIdx.x), tx);
  temp_0 = tx;
	t20 = SCALAR_MUL(X_PA, t10);
	t20 = SCALAR_FNMA(X_PC, t11, t20);
	tx = SCALAR_SUB(t00, t01);
	ty = SCALAR_SET1(0.5 * 1);
	ty = SCALAR_MUL(ty, RHO_INV);
	t20 = SCALAR_FMA(tx, ty, t20);
	//tx = SCALAR_LOAD((temp + 3 * blockDim.x + threadIdx.x));
  tx = temp_3;
	tx = SCALAR_ADD(tx, t20);
	//SCALAR_STORE((temp + 3 * blockDim.x + threadIdx.x), tx);
  temp_3 = tx;
	t20 = SCALAR_MUL(Y_PA, t10);
	t20 = SCALAR_FNMA(Y_PC, t11, t20);
	//tx = SCALAR_LOAD((temp + 4 * blockDim.x + threadIdx.x));
  tx = temp_4;
	tx = SCALAR_ADD(tx, t20);
	//SCALAR_STORE((temp + 4 * blockDim.x + threadIdx.x), tx);
  temp_4 = tx;
	t20 = SCALAR_MUL(Z_PA, t10);
	t20 = SCALAR_FNMA(Z_PC, t11, t20);
	//tx = SCALAR_LOAD((temp + 5 * blockDim.x + threadIdx.x));
  tx = temp_5;
	tx = SCALAR_ADD(tx, t20);
	//SCALAR_STORE((temp + 5 * blockDim.x + threadIdx.x), tx);
  temp_5 = tx;
	t10 = SCALAR_MUL(Y_PA, t00);
	t10 = SCALAR_FNMA(Y_PC, t01, t10);
	t11 = SCALAR_MUL(Y_PA, t01);
	t11 = SCALAR_FNMA(Y_PC, t02, t11);
	//tx = SCALAR_LOAD((temp + 1 * blockDim.x + threadIdx.x));
  tx = temp_1;
	tx = SCALAR_ADD(tx, t10);
	//SCALAR_STORE((temp + 1 * blockDim.x + threadIdx.x), tx);
  temp_1 = tx;
	t20 = SCALAR_MUL(Y_PA, t10);
	t20 = SCALAR_FNMA(Y_PC, t11, t20);
	tx = SCALAR_SUB(t00, t01);
	ty = SCALAR_SET1(0.5 * 1);
	ty = SCALAR_MUL(ty, RHO_INV);
	t20 = SCALAR_FMA(tx, ty, t20);
	//tx = SCALAR_LOAD((temp + 6 * blockDim.x + threadIdx.x));
  tx = temp_6;
	tx = SCALAR_ADD(tx, t20);
	//SCALAR_STORE((temp + 6 * blockDim.x + threadIdx.x), tx);
  temp_6 = tx;
	t20 = SCALAR_MUL(Z_PA, t10);
	t20 = SCALAR_FNMA(Z_PC, t11, t20);
	//tx = SCALAR_LOAD((temp + 7 * blockDim.x + threadIdx.x));
  tx = temp_7;
	tx = SCALAR_ADD(tx, t20);
	//SCALAR_STORE((temp + 7 * blockDim.x + threadIdx.x), tx);
  temp_7 = tx;
	t10 = SCALAR_MUL(Z_PA, t00);
	t10 = SCALAR_FNMA(Z_PC, t01, t10);
	t11 = SCALAR_MUL(Z_PA, t01);
	t11 = SCALAR_FNMA(Z_PC, t02, t11);
	//tx = SCALAR_LOAD((temp + 2 * blockDim.x + threadIdx.x));
  tx = temp_2;
	tx = SCALAR_ADD(tx, t10);
	//SCALAR_STORE((temp + 2 * blockDim.x + threadIdx.x), tx);
  temp_2 = tx;
	t20 = SCALAR_MUL(Z_PA, t10);
	t20 = SCALAR_FNMA(Z_PC, t11, t20);
	tx = SCALAR_SUB(t00, t01);
	ty = SCALAR_SET1(0.5 * 1);
	ty = SCALAR_MUL(ty, RHO_INV);
	t20 = SCALAR_FMA(tx, ty, t20);
	//tx = SCALAR_LOAD((temp + 8 * blockDim.x + threadIdx.x));
  tx = temp_8;
	tx = SCALAR_ADD(tx, t20);
	//SCALAR_STORE((temp + 8 * blockDim.x + threadIdx.x), tx);
  temp_8 = tx;
      }

      if(threadIdx.x < npts - p_outer) {
	double *Xik = (Xi + p_outer + p_inner);
	double *Xjk = (Xj + p_outer + p_inner);
	double *Gik = (Gi + p_outer + p_inner);
	double *Gjk = (Gj + p_outer + p_inner);

	SCALAR_TYPE const_value_v = SCALAR_LOAD((weights + p_outer + p_inner));

	double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k;
	SCALAR_TYPE const_value_w;
	SCALAR_TYPE tx, ty, tz, tw, t0, t1, t2;

  #if 0
	X_ABp = 1.0; comb_m_i = 1.0;
	Y_ABp = 1.0; comb_n_j = 1.0;
	Z_ABp = 1.0; comb_p_k = 1.0;
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);
	tx = SCALAR_LOAD((Xik + 0 * ldX));
	ty = SCALAR_LOAD((Xjk + 0 * ldX));
	tz = SCALAR_LOAD((Gik + 0 * ldG));
	tw = SCALAR_LOAD((Gjk + 0 * ldG));
	//t0 = SCALAR_LOAD((temp + 3 * blockDim.x + threadIdx.x));
  t0 = temp_3;
	t0 = SCALAR_MUL(t0, const_value_w);
	tz = SCALAR_FMA(ty, t0, tz);
	tw = SCALAR_FMA(tx, t0, tw);
	SCALAR_STORE((Gik + 0 * ldG), tz);
	SCALAR_STORE((Gjk + 0 * ldG), tw);
	tx = SCALAR_LOAD((Xik + 1 * ldX));
	ty = SCALAR_LOAD((Xjk + 0 * ldX));
	tz = SCALAR_LOAD((Gik + 1 * ldG));
	tw = SCALAR_LOAD((Gjk + 0 * ldG));
	//t1 = SCALAR_LOAD((temp + 4 * blockDim.x + threadIdx.x));
  t1 = temp_4;
	t1 = SCALAR_MUL(t1, const_value_w);
	tz = SCALAR_FMA(ty, t1, tz);
	tw = SCALAR_FMA(tx, t1, tw);
	SCALAR_STORE((Gik + 1 * ldG), tz);
	SCALAR_STORE((Gjk + 0 * ldG), tw);
	tx = SCALAR_LOAD((Xik + 2 * ldX));
	ty = SCALAR_LOAD((Xjk + 0 * ldX));
	tz = SCALAR_LOAD((Gik + 2 * ldG));
	tw = SCALAR_LOAD((Gjk + 0 * ldG));
	//t2 = SCALAR_LOAD((temp + 5 * blockDim.x + threadIdx.x));
  t2 = temp_5;
	t2 = SCALAR_MUL(t2, const_value_w);
	tz = SCALAR_FMA(ty, t2, tz);
	tw = SCALAR_FMA(tx, t2, tw);
	SCALAR_STORE((Gik + 2 * ldG), tz);
	SCALAR_STORE((Gjk + 0 * ldG), tw);
	X_ABp = SCALAR_MUL(X_ABp, X_AB); comb_m_i = SCALAR_MUL(comb_m_i * 1, SCALAR_RECIPROCAL(1));
	Y_ABp = 1.0; comb_n_j = 1.0;
	Z_ABp = 1.0; comb_p_k = 1.0;
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);
	tx = SCALAR_LOAD((Xik + 0 * ldX));
	ty = SCALAR_LOAD((Xjk + 0 * ldX));
	tz = SCALAR_LOAD((Gik + 0 * ldG));
	tw = SCALAR_LOAD((Gjk + 0 * ldG));
	//t0 = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
  t0 = temp_0;
	t0 = SCALAR_MUL(t0, const_value_w);
	tz = SCALAR_FMA(ty, t0, tz);
	tw = SCALAR_FMA(tx, t0, tw);
	SCALAR_STORE((Gik + 0 * ldG), tz);
	SCALAR_STORE((Gjk + 0 * ldG), tw);
	tx = SCALAR_LOAD((Xik + 1 * ldX));
	ty = SCALAR_LOAD((Xjk + 0 * ldX));
	tz = SCALAR_LOAD((Gik + 1 * ldG));
	tw = SCALAR_LOAD((Gjk + 0 * ldG));
	//t1 = SCALAR_LOAD((temp + 1 * blockDim.x + threadIdx.x));
  t1 = temp_1;
	t1 = SCALAR_MUL(t1, const_value_w);
	tz = SCALAR_FMA(ty, t1, tz);
	tw = SCALAR_FMA(tx, t1, tw);
	SCALAR_STORE((Gik + 1 * ldG), tz);
	SCALAR_STORE((Gjk + 0 * ldG), tw);
	tx = SCALAR_LOAD((Xik + 2 * ldX));
	ty = SCALAR_LOAD((Xjk + 0 * ldX));
	tz = SCALAR_LOAD((Gik + 2 * ldG));
	tw = SCALAR_LOAD((Gjk + 0 * ldG));
	//t2 = SCALAR_LOAD((temp + 2 * blockDim.x + threadIdx.x));
  t2 = temp_2;
	t2 = SCALAR_MUL(t2, const_value_w);
	tz = SCALAR_FMA(ty, t2, tz);
	tw = SCALAR_FMA(tx, t2, tw);
	SCALAR_STORE((Gik + 2 * ldG), tz);
	SCALAR_STORE((Gjk + 0 * ldG), tw);
	X_ABp = 1.0; comb_m_i = 1.0;
	Y_ABp = 1.0; comb_n_j = 1.0;
	Z_ABp = 1.0; comb_p_k = 1.0;
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);
	tx = SCALAR_LOAD((Xik + 0 * ldX));
	ty = SCALAR_LOAD((Xjk + 1 * ldX));
	tz = SCALAR_LOAD((Gik + 0 * ldG));
	tw = SCALAR_LOAD((Gjk + 1 * ldG));
	//t0 = SCALAR_LOAD((temp + 4 * blockDim.x + threadIdx.x));
  t0 = temp_4;
	t0 = SCALAR_MUL(t0, const_value_w);
	tz = SCALAR_FMA(ty, t0, tz);
	tw = SCALAR_FMA(tx, t0, tw);
	SCALAR_STORE((Gik + 0 * ldG), tz);
	SCALAR_STORE((Gjk + 1 * ldG), tw);
	tx = SCALAR_LOAD((Xik + 1 * ldX));
	ty = SCALAR_LOAD((Xjk + 1 * ldX));
	tz = SCALAR_LOAD((Gik + 1 * ldG));
	tw = SCALAR_LOAD((Gjk + 1 * ldG));
	//t1 = SCALAR_LOAD((temp + 6 * blockDim.x + threadIdx.x));
  t1 = temp_6;
	t1 = SCALAR_MUL(t1, const_value_w);
	tz = SCALAR_FMA(ty, t1, tz);
	tw = SCALAR_FMA(tx, t1, tw);
	SCALAR_STORE((Gik + 1 * ldG), tz);
	SCALAR_STORE((Gjk + 1 * ldG), tw);
	tx = SCALAR_LOAD((Xik + 2 * ldX));
	ty = SCALAR_LOAD((Xjk + 1 * ldX));
	tz = SCALAR_LOAD((Gik + 2 * ldG));
	tw = SCALAR_LOAD((Gjk + 1 * ldG));
	//t2 = SCALAR_LOAD((temp + 7 * blockDim.x + threadIdx.x));
  t2 = temp_7;
	t2 = SCALAR_MUL(t2, const_value_w);
	tz = SCALAR_FMA(ty, t2, tz);
	tw = SCALAR_FMA(tx, t2, tw);
	SCALAR_STORE((Gik + 2 * ldG), tz);
	SCALAR_STORE((Gjk + 1 * ldG), tw);
	Y_ABp = SCALAR_MUL(Y_ABp, Y_AB); comb_n_j = SCALAR_MUL(comb_n_j * 1, SCALAR_RECIPROCAL(1));
	Z_ABp = 1.0; comb_p_k = 1.0;
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);
	tx = SCALAR_LOAD((Xik + 0 * ldX));
	ty = SCALAR_LOAD((Xjk + 1 * ldX));
	tz = SCALAR_LOAD((Gik + 0 * ldG));
	tw = SCALAR_LOAD((Gjk + 1 * ldG));
	//t0 = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
  t0 = temp_0;
	t0 = SCALAR_MUL(t0, const_value_w);
	tz = SCALAR_FMA(ty, t0, tz);
	tw = SCALAR_FMA(tx, t0, tw);
	SCALAR_STORE((Gik + 0 * ldG), tz);
	SCALAR_STORE((Gjk + 1 * ldG), tw);
	tx = SCALAR_LOAD((Xik + 1 * ldX));
	ty = SCALAR_LOAD((Xjk + 1 * ldX));
	tz = SCALAR_LOAD((Gik + 1 * ldG));
	tw = SCALAR_LOAD((Gjk + 1 * ldG));
	//t1 = SCALAR_LOAD((temp + 1 * blockDim.x + threadIdx.x));
  t1 = temp_1;
	t1 = SCALAR_MUL(t1, const_value_w);
	tz = SCALAR_FMA(ty, t1, tz);
	tw = SCALAR_FMA(tx, t1, tw);
	SCALAR_STORE((Gik + 1 * ldG), tz);
	SCALAR_STORE((Gjk + 1 * ldG), tw);
	tx = SCALAR_LOAD((Xik + 2 * ldX));
	ty = SCALAR_LOAD((Xjk + 1 * ldX));
	tz = SCALAR_LOAD((Gik + 2 * ldG));
	tw = SCALAR_LOAD((Gjk + 1 * ldG));
	//t2 = SCALAR_LOAD((temp + 2 * blockDim.x + threadIdx.x));
  t2 = temp_2;
	t2 = SCALAR_MUL(t2, const_value_w);
	tz = SCALAR_FMA(ty, t2, tz);
	tw = SCALAR_FMA(tx, t2, tw);
	SCALAR_STORE((Gik + 2 * ldG), tz);
	SCALAR_STORE((Gjk + 1 * ldG), tw);
	X_ABp = 1.0; comb_m_i = 1.0;
	Y_ABp = 1.0; comb_n_j = 1.0;
	Z_ABp = 1.0; comb_p_k = 1.0;
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);
	tx = SCALAR_LOAD((Xik + 0 * ldX));
	ty = SCALAR_LOAD((Xjk + 2 * ldX));
	tz = SCALAR_LOAD((Gik + 0 * ldG));
	tw = SCALAR_LOAD((Gjk + 2 * ldG));
	//t0 = SCALAR_LOAD((temp + 5 * blockDim.x + threadIdx.x));
  t0 = temp_5;
	t0 = SCALAR_MUL(t0, const_value_w);
	tz = SCALAR_FMA(ty, t0, tz);
	tw = SCALAR_FMA(tx, t0, tw);
	SCALAR_STORE((Gik + 0 * ldG), tz);
	SCALAR_STORE((Gjk + 2 * ldG), tw);
	tx = SCALAR_LOAD((Xik + 1 * ldX));
	ty = SCALAR_LOAD((Xjk + 2 * ldX));
	tz = SCALAR_LOAD((Gik + 1 * ldG));
	tw = SCALAR_LOAD((Gjk + 2 * ldG));
	//t1 = SCALAR_LOAD((temp + 7 * blockDim.x + threadIdx.x));
  t1 = temp_7;
	t1 = SCALAR_MUL(t1, const_value_w);
	tz = SCALAR_FMA(ty, t1, tz);
	tw = SCALAR_FMA(tx, t1, tw);
	SCALAR_STORE((Gik + 1 * ldG), tz);
	SCALAR_STORE((Gjk + 2 * ldG), tw);
	tx = SCALAR_LOAD((Xik + 2 * ldX));
	ty = SCALAR_LOAD((Xjk + 2 * ldX));
	tz = SCALAR_LOAD((Gik + 2 * ldG));
	tw = SCALAR_LOAD((Gjk + 2 * ldG));
	//t2 = SCALAR_LOAD((temp + 8 * blockDim.x + threadIdx.x));
  t2 = temp_8;
	t2 = SCALAR_MUL(t2, const_value_w);
	tz = SCALAR_FMA(ty, t2, tz);
	tw = SCALAR_FMA(tx, t2, tw);
	SCALAR_STORE((Gik + 2 * ldG), tz);
	SCALAR_STORE((Gjk + 2 * ldG), tw);
	Z_ABp = SCALAR_MUL(Z_ABp, Z_AB); comb_p_k = SCALAR_MUL(comb_p_k * 1, SCALAR_RECIPROCAL(1));
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);
	tx = SCALAR_LOAD((Xik + 0 * ldX));
	ty = SCALAR_LOAD((Xjk + 2 * ldX));
	tz = SCALAR_LOAD((Gik + 0 * ldG));
	tw = SCALAR_LOAD((Gjk + 2 * ldG));
	//t0 = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
  t0 = temp_0;
	t0 = SCALAR_MUL(t0, const_value_w);
	tz = SCALAR_FMA(ty, t0, tz);
	tw = SCALAR_FMA(tx, t0, tw);
	SCALAR_STORE((Gik + 0 * ldG), tz);
	SCALAR_STORE((Gjk + 2 * ldG), tw);
	tx = SCALAR_LOAD((Xik + 1 * ldX));
	ty = SCALAR_LOAD((Xjk + 2 * ldX));
	tz = SCALAR_LOAD((Gik + 1 * ldG));
	tw = SCALAR_LOAD((Gjk + 2 * ldG));
	//t1 = SCALAR_LOAD((temp + 1 * blockDim.x + threadIdx.x));
  t1 = temp_1;
	t1 = SCALAR_MUL(t1, const_value_w);
	tz = SCALAR_FMA(ty, t1, tz);
	tw = SCALAR_FMA(tx, t1, tw);
	SCALAR_STORE((Gik + 1 * ldG), tz);
	SCALAR_STORE((Gjk + 2 * ldG), tw);
	tx = SCALAR_LOAD((Xik + 2 * ldX));
	ty = SCALAR_LOAD((Xjk + 2 * ldX));
	tz = SCALAR_LOAD((Gik + 2 * ldG));
	tw = SCALAR_LOAD((Gjk + 2 * ldG));
	//t2 = SCALAR_LOAD((temp + 2 * blockDim.x + threadIdx.x));
  t2 = temp_2;
	t2 = SCALAR_MUL(t2, const_value_w);
	tz = SCALAR_FMA(ty, t2, tz);
	tw = SCALAR_FMA(tx, t2, tw);
	SCALAR_STORE((Gik + 2 * ldG), tz);
	SCALAR_STORE((Gjk + 2 * ldG), tw);
  #else

  /**** j = 0 ****/
	X_ABp = 1.0; comb_m_i = 1.0;
	Y_ABp = 1.0; comb_n_j = 1.0;
	Z_ABp = 1.0; comb_p_k = 1.0;
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);

	tx = SCALAR_LOAD((Xik + 0 * ldX));
	ty = SCALAR_LOAD((Xjk + 0 * ldX));
	t0 = SCALAR_MUL(temp_3, const_value_w);
	tz = SCALAR_MUL(ty, t0);
	tw = SCALAR_MUL(tx, t0);
	atomicAdd((Gik + 0 * ldG), tz);

	tx = SCALAR_LOAD((Xik + 1 * ldX));
	t1 = SCALAR_MUL(temp_4, const_value_w);
	tz = SCALAR_MUL(ty, t1);
	tw = SCALAR_FMA(tx, t1, tw);
	atomicAdd((Gik + 1 * ldG), tz);

	tx = SCALAR_LOAD((Xik + 2 * ldX));
	t2 = SCALAR_MUL(temp_5, const_value_w);
	tz = SCALAR_MUL(ty, t2);
	tw = SCALAR_FMA(tx, t2, tw);
	atomicAdd((Gik + 2 * ldG), tz);

	X_ABp = SCALAR_MUL(X_ABp, X_AB); comb_m_i = SCALAR_MUL(comb_m_i * 1, SCALAR_RECIPROCAL(1));
	Y_ABp = 1.0; comb_n_j = 1.0;
	Z_ABp = 1.0; comb_p_k = 1.0;
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);

	tx = SCALAR_LOAD((Xik + 0 * ldX));
	t0 = SCALAR_MUL(temp_0, const_value_w);
	tz = SCALAR_MUL(ty, t0);
	tw = SCALAR_FMA(tx, t0, tw);
	atomicAdd((Gik + 0 * ldG), tz);

	tx = SCALAR_LOAD((Xik + 1 * ldX));
	t1 = SCALAR_MUL(temp_1, const_value_w);
	tz = SCALAR_MUL(ty, t1);
	tw = SCALAR_FMA(tx, t1, tw);
	atomicAdd((Gik + 1 * ldG), tz);

	tx = SCALAR_LOAD((Xik + 2 * ldX));
	t2 = SCALAR_MUL(temp_2, const_value_w);
	tz = SCALAR_MUL(ty, t2);
	tw = SCALAR_FMA(tx, t2, tw);
	atomicAdd((Gik + 2 * ldG), tz);
	atomicAdd((Gjk + 0 * ldG), tw);



  /**** j = 1 ****/
	X_ABp = 1.0; comb_m_i = 1.0;
	Y_ABp = 1.0; comb_n_j = 1.0;
	Z_ABp = 1.0; comb_p_k = 1.0;
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);

	tx = SCALAR_LOAD((Xik + 0 * ldX));
	ty = SCALAR_LOAD((Xjk + 1 * ldX));
	t0 = SCALAR_MUL(temp_4, const_value_w);
	tz = SCALAR_MUL(ty, t0);
	tw = SCALAR_MUL(tx, t0);
	atomicAdd((Gik + 0 * ldG), tz);

	tx = SCALAR_LOAD((Xik + 1 * ldX));
	t1 = SCALAR_MUL(temp_6, const_value_w);
	tz = SCALAR_MUL(ty, t1);
	tw = SCALAR_FMA(tx, t1, tw);
	atomicAdd((Gik + 1 * ldG), tz);

	tx = SCALAR_LOAD((Xik + 2 * ldX));
	t2 = SCALAR_MUL(temp_7, const_value_w);
	tz = SCALAR_MUL(ty, t2);
	tw = SCALAR_FMA(tx, t2, tw);
	atomicAdd((Gik + 2 * ldG), tz);

	Y_ABp = SCALAR_MUL(Y_ABp, Y_AB); comb_n_j = SCALAR_MUL(comb_n_j * 1, SCALAR_RECIPROCAL(1));
	Z_ABp = 1.0; comb_p_k = 1.0;
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);

	tx = SCALAR_LOAD((Xik + 0 * ldX));
	t0 = SCALAR_MUL(temp_0, const_value_w);
	tz = SCALAR_MUL(ty, t0);
	tw = SCALAR_FMA(tx, t0, tw);
	atomicAdd((Gik + 0 * ldG), tz);

	tx = SCALAR_LOAD((Xik + 1 * ldX));
	t1 = SCALAR_MUL(temp_1, const_value_w);
	tz = SCALAR_MUL(ty, t1);
	tw = SCALAR_FMA(tx, t1, tw);
	atomicAdd((Gik + 1 * ldG), tz);

	tx = SCALAR_LOAD((Xik + 2 * ldX));
	t2 = SCALAR_MUL(temp_2, const_value_w);
	tz = SCALAR_MUL(ty, t2);
	tw = SCALAR_FMA(tx, t2, tw);
	atomicAdd((Gik + 2 * ldG), tz);
	atomicAdd((Gjk + 1 * ldG), tw);


  /**** j = 2 ****/
	X_ABp = 1.0; comb_m_i = 1.0;
	Y_ABp = 1.0; comb_n_j = 1.0;
	Z_ABp = 1.0; comb_p_k = 1.0;
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);

	tx = SCALAR_LOAD((Xik + 0 * ldX));
	ty = SCALAR_LOAD((Xjk + 2 * ldX));
	t0 = SCALAR_MUL(temp_5, const_value_w);
	tz = SCALAR_MUL(ty, t0);
	tw = SCALAR_MUL(tx, t0);
	atomicAdd((Gik + 0 * ldG), tz);

	tx = SCALAR_LOAD((Xik + 1 * ldX));
	t1 = SCALAR_MUL(temp_7, const_value_w);
	tz = SCALAR_MUL(ty, t1);
	tw = SCALAR_FMA(tx, t1, tw);
	atomicAdd((Gik + 1 * ldG), tz);

	tx = SCALAR_LOAD((Xik + 2 * ldX));
	t2 = SCALAR_MUL(temp_8, const_value_w);
	tz = SCALAR_MUL(ty, t2);
	tw = SCALAR_FMA(tx, t2, tw);
	atomicAdd((Gik + 2 * ldG), tz);

	Z_ABp = SCALAR_MUL(Z_ABp, Z_AB); comb_p_k = SCALAR_MUL(comb_p_k * 1, SCALAR_RECIPROCAL(1));
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);

	tx = SCALAR_LOAD((Xik + 0 * ldX));
	t0 = SCALAR_MUL(temp_0, const_value_w);
	tz = SCALAR_MUL(ty, t0);
	tw = SCALAR_FMA(tx, t0, tw);
	atomicAdd((Gik + 0 * ldG), tz);

	tx = SCALAR_LOAD((Xik + 1 * ldX));
	t1 = SCALAR_MUL(temp_1, const_value_w);
	tz = SCALAR_MUL(ty, t1);
	tw = SCALAR_FMA(tx, t1, tw);
	atomicAdd((Gik + 1 * ldG), tz);

	tx = SCALAR_LOAD((Xik + 2 * ldX));
	t2 = SCALAR_MUL(temp_2, const_value_w);
	tz = SCALAR_MUL(ty, t2);
	tw = SCALAR_FMA(tx, t2, tw);
	atomicAdd((Gik + 2 * ldG), tz);
	atomicAdd((Gjk + 2 * ldG), tw);
  #endif
      }
    }
  }

  __global__ void dev_integral_1_1(
           double X_AB,
				   double Y_AB,
				   double Z_AB,
           size_t npts,
				   double *points_x,
				   double *points_y,
				   double *points_z,
           shell_pair* sp,
				   double *Xi,
				   double *Xj,
				   int ldX,
				   double *Gi,
				   double *Gj,
				   int ldG, 
				   double *weights, 
				   double *boys_table) {
    dev_integral_1_1_driver( X_AB, Y_AB, Z_AB, npts, points_x, points_y, 
      points_z, sp, Xi, Xj, ldX, Gi, Gj, ldG, weights, boys_table );
  }

  void integral_1_1(double X_AB,
		    double Y_AB,
		    double Z_AB,
		    size_t npts,
		    double *points_x,
		    double *points_y,
		    double *points_z,
        shell_pair* sp,
		    double *Xi,
		    double *Xj,
		    int ldX,
		    double *Gi,
		    double *Gj,
		    int ldG, 
		    double *weights, 
		  double *boys_table,
      cudaStream_t stream) {
    dev_integral_1_1<<<320, 128, 0, stream>>>(X_AB,
				   Y_AB,
				   Z_AB,
				   npts,
				   points_x,
				   points_y,
				   points_z,
           sp,
				   Xi,
				   Xj,
				   ldX,
				   Gi,
				   Gj,
				   ldG, 
				   weights,
				   boys_table);
  }

  __inline__ __device__ void dev_integral_1_1_batched_driver(
           double X_AB,
				   double Y_AB,
				   double Z_AB,
           const GauXC::ShellPairToTaskDevice* sp2task,
           GauXC::XCDeviceTask*                device_tasks,
				   double *boys_table) {

    const int ntask = sp2task->ntask;
    #pragma unroll 1
    for( int i_task = blockIdx.y; i_task < ntask; i_task += gridDim.y ) {
    
      const auto iT = sp2task->task_idx_device[i_task];
      const auto* task  = device_tasks + iT;
      const auto  npts  = task->npts;

      const auto  i_off = sp2task->task_shell_off_row_device[i_task]*npts;
      const auto  j_off = sp2task->task_shell_off_col_device[i_task]*npts;


      dev_integral_1_1_driver( 
        X_AB, Y_AB, Z_AB,
        npts,
        task->points_x,
        task->points_y,
        task->points_z,
        sp2task->shell_pair_device,
        task->fmat + i_off,
        task->fmat + j_off,
        npts,
        task->gmat + i_off,
        task->gmat + j_off,
        npts,
        task->weights, boys_table );
    }

  }
  __global__ void dev_integral_1_1_batched(
           double X_AB,
				   double Y_AB,
				   double Z_AB,
           const GauXC::ShellPairToTaskDevice* sp2task,
           GauXC::XCDeviceTask*                device_tasks,
				   double *boys_table) {
   dev_integral_1_1_batched_driver(X_AB,Y_AB,Z_AB,sp2task,device_tasks,boys_table);
 }



  void integral_1_1_batched(size_t ntask_sp,
        double X_AB,
				double Y_AB,
				double Z_AB,
        const GauXC::ShellPairToTaskDevice* sp2task,
        GauXC::XCDeviceTask*                device_tasks,
		    double *boys_table,
        cudaStream_t stream) {

    int nthreads = 128;
    int nblocks_x = 160;
    int nblocks_y = ntask_sp;
    dim3 nblocks(nblocks_x, nblocks_y);

    dev_integral_1_1_batched<<<nblocks,nthreads,0,stream>>>(
      X_AB, Y_AB, Z_AB, sp2task, device_tasks, boys_table );

  }



  __global__ void dev_integral_1_1_shell_batched(
           int nsp,
           const GauXC::ShellPairToTaskDevice* sp2task,
           GauXC::XCDeviceTask*                device_tasks,
				   double *boys_table) {

   for( int i = blockIdx.z; i < nsp; i += gridDim.z ) {
     auto sp = sp2task + i;
     const auto X_AB = sp->X_AB;
     const auto Y_AB = sp->Y_AB;
     const auto Z_AB = sp->Z_AB;
     dev_integral_1_1_batched_driver(X_AB,Y_AB,Z_AB,sp,device_tasks,boys_table);
   }
 }

  void integral_1_1_shell_batched(
        size_t nsp,
        size_t max_ntask,
        const GauXC::ShellPairToTaskDevice* sp2task,
        GauXC::XCDeviceTask*                device_tasks,
		    double *boys_table,
        cudaStream_t stream) {

    size_t xy_max = (1ul << 16) - 1;
    int nthreads = 128;
    int nblocks_x = 1;
    int nblocks_y = std::min(max_ntask, xy_max);
    int nblocks_z = std::min(nsp,  xy_max);
    dim3 nblocks(nblocks_x, nblocks_y, nblocks_z);
    dev_integral_1_1_shell_batched<<<nblocks,nthreads,0,stream>>>(
      nsp, sp2task, device_tasks, boys_table );

  }
}
