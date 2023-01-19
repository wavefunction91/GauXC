#include <math.h>
#include "../include/gpu/chebyshev_boys_computation.hpp"
#include "config_obara_saika.hpp"
#include "integral_1.hu"

namespace XGPU {
  __inline__ __device__ void dev_integral_1_driver(size_t npts,
				 double *points_x,
				 double *points_y,
				 double *points_z,
         shell_pair* sp,
				 double *Xi,
				 int ldX,
				 double *Gi,
				 int ldG,
				 double *weights,
				 double *boys_table) {
    __shared__ double temp[128 * 9];
    const auto nprim_pairs = sp->nprim_pairs();
    const auto prim_pairs  = sp->prim_pairs();
    
    for(size_t p_outer = blockIdx.x * blockDim.x; p_outer < npts; p_outer += gridDim.x * blockDim.x) {
      double *_point_outer_x = (points_x + p_outer);
      double *_point_outer_y = (points_y + p_outer);
      double *_point_outer_z = (points_z + p_outer);

      size_t p_inner = (threadIdx.x < (npts - p_outer)) ? threadIdx.x : (npts - p_outer);

      for(int i = 0; i < 9; ++i) SCALAR_STORE((temp + i * blockDim.x + threadIdx.x), SCALAR_ZERO());

      for(int ij = 0; ij < nprim_pairs; ++ij) {
	double RHO = prim_pairs[ij].gamma;
	double RHO_INV = prim_pairs[ij].gamma_inv;

	double xA = prim_pairs[ij].P.x;
	double yA = prim_pairs[ij].P.y;
	double zA = prim_pairs[ij].P.z;
	
	constexpr double X_PA = 0.0;
	constexpr double Y_PA = 0.0;
	constexpr double Z_PA = 0.0;

	double eval = prim_pairs[ij].K_coeff_prod;

	// Evaluate T Values
	SCALAR_TYPE xC = SCALAR_LOAD((_point_outer_x + p_inner));
	SCALAR_TYPE yC = SCALAR_LOAD((_point_outer_y + p_inner));
	SCALAR_TYPE zC = SCALAR_LOAD((_point_outer_z + p_inner));

	SCALAR_TYPE X_PC = SCALAR_SUB(xA, xC);
	SCALAR_TYPE Y_PC = SCALAR_SUB(yA, yC);
	SCALAR_TYPE Z_PC = SCALAR_SUB(zA, zC);

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
	tx = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
	tx = SCALAR_ADD(tx, t10);
	SCALAR_STORE((temp + 0 * blockDim.x + threadIdx.x), tx);
	t20 = SCALAR_MUL(X_PA, t10);
	t20 = SCALAR_FNMA(X_PC, t11, t20);
	tx = SCALAR_SUB(t00, t01);
	ty = SCALAR_SET1(0.5 * 1);
	ty = SCALAR_MUL(ty, RHO_INV);
	t20 = SCALAR_FMA(tx, ty, t20);
	tx = SCALAR_LOAD((temp + 3 * blockDim.x + threadIdx.x));
	tx = SCALAR_ADD(tx, t20);
	SCALAR_STORE((temp + 3 * blockDim.x + threadIdx.x), tx);
	t20 = SCALAR_MUL(Y_PA, t10);
	t20 = SCALAR_FNMA(Y_PC, t11, t20);
	tx = SCALAR_LOAD((temp + 4 * blockDim.x + threadIdx.x));
	tx = SCALAR_ADD(tx, t20);
	SCALAR_STORE((temp + 4 * blockDim.x + threadIdx.x), tx);
	t20 = SCALAR_MUL(Z_PA, t10);
	t20 = SCALAR_FNMA(Z_PC, t11, t20);
	tx = SCALAR_LOAD((temp + 5 * blockDim.x + threadIdx.x));
	tx = SCALAR_ADD(tx, t20);
	SCALAR_STORE((temp + 5 * blockDim.x + threadIdx.x), tx);
	t10 = SCALAR_MUL(Y_PA, t00);
	t10 = SCALAR_FNMA(Y_PC, t01, t10);
	t11 = SCALAR_MUL(Y_PA, t01);
	t11 = SCALAR_FNMA(Y_PC, t02, t11);
	tx = SCALAR_LOAD((temp + 1 * blockDim.x + threadIdx.x));
	tx = SCALAR_ADD(tx, t10);
	SCALAR_STORE((temp + 1 * blockDim.x + threadIdx.x), tx);
	t20 = SCALAR_MUL(Y_PA, t10);
	t20 = SCALAR_FNMA(Y_PC, t11, t20);
	tx = SCALAR_SUB(t00, t01);
	ty = SCALAR_SET1(0.5 * 1);
	ty = SCALAR_MUL(ty, RHO_INV);
	t20 = SCALAR_FMA(tx, ty, t20);
	tx = SCALAR_LOAD((temp + 6 * blockDim.x + threadIdx.x));
	tx = SCALAR_ADD(tx, t20);
	SCALAR_STORE((temp + 6 * blockDim.x + threadIdx.x), tx);
	t20 = SCALAR_MUL(Z_PA, t10);
	t20 = SCALAR_FNMA(Z_PC, t11, t20);
	tx = SCALAR_LOAD((temp + 7 * blockDim.x + threadIdx.x));
	tx = SCALAR_ADD(tx, t20);
	SCALAR_STORE((temp + 7 * blockDim.x + threadIdx.x), tx);
	t10 = SCALAR_MUL(Z_PA, t00);
	t10 = SCALAR_FNMA(Z_PC, t01, t10);
	t11 = SCALAR_MUL(Z_PA, t01);
	t11 = SCALAR_FNMA(Z_PC, t02, t11);
	tx = SCALAR_LOAD((temp + 2 * blockDim.x + threadIdx.x));
	tx = SCALAR_ADD(tx, t10);
	SCALAR_STORE((temp + 2 * blockDim.x + threadIdx.x), tx);
	t20 = SCALAR_MUL(Z_PA, t10);
	t20 = SCALAR_FNMA(Z_PC, t11, t20);
	tx = SCALAR_SUB(t00, t01);
	ty = SCALAR_SET1(0.5 * 1);
	ty = SCALAR_MUL(ty, RHO_INV);
	t20 = SCALAR_FMA(tx, ty, t20);
	tx = SCALAR_LOAD((temp + 8 * blockDim.x + threadIdx.x));
	tx = SCALAR_ADD(tx, t20);
	SCALAR_STORE((temp + 8 * blockDim.x + threadIdx.x), tx);
      }

      if(threadIdx.x < npts - p_outer) {
	double *Xik = (Xi + p_outer + p_inner);
	double *Gik = (Gi + p_outer + p_inner);

	SCALAR_TYPE tx, wg, xik, gik;
	tx  = SCALAR_LOAD((temp + 3 * blockDim.x + threadIdx.x));
	wg  = SCALAR_LOAD((weights + p_outer + p_inner));

	xik = SCALAR_LOAD((Xik + 0 * ldX));
	gik = SCALAR_LOAD((Gik + 0 * ldG));

	tx = SCALAR_MUL(tx, wg);
	gik = SCALAR_FMA(tx, xik, gik);
	SCALAR_STORE((Gik + 0 * ldG), gik);
	tx  = SCALAR_LOAD((temp + 4 * blockDim.x + threadIdx.x));
	wg  = SCALAR_LOAD((weights + p_outer + p_inner));

	xik = SCALAR_LOAD((Xik + 0 * ldX));
	gik = SCALAR_LOAD((Gik + 1 * ldG));

	tx = SCALAR_MUL(tx, wg);
	gik = SCALAR_FMA(tx, xik, gik);
	SCALAR_STORE((Gik + 1 * ldG), gik);
	tx  = SCALAR_LOAD((temp + 5 * blockDim.x + threadIdx.x));
	wg  = SCALAR_LOAD((weights + p_outer + p_inner));

	xik = SCALAR_LOAD((Xik + 0 * ldX));
	gik = SCALAR_LOAD((Gik + 2 * ldG));

	tx = SCALAR_MUL(tx, wg);
	gik = SCALAR_FMA(tx, xik, gik);
	SCALAR_STORE((Gik + 2 * ldG), gik);
	tx  = SCALAR_LOAD((temp + 4 * blockDim.x + threadIdx.x));
	wg  = SCALAR_LOAD((weights + p_outer + p_inner));

	xik = SCALAR_LOAD((Xik + 1 * ldX));
	gik = SCALAR_LOAD((Gik + 0 * ldG));

	tx = SCALAR_MUL(tx, wg);
	gik = SCALAR_FMA(tx, xik, gik);
	SCALAR_STORE((Gik + 0 * ldG), gik);
	tx  = SCALAR_LOAD((temp + 6 * blockDim.x + threadIdx.x));
	wg  = SCALAR_LOAD((weights + p_outer + p_inner));

	xik = SCALAR_LOAD((Xik + 1 * ldX));
	gik = SCALAR_LOAD((Gik + 1 * ldG));

	tx = SCALAR_MUL(tx, wg);
	gik = SCALAR_FMA(tx, xik, gik);
	SCALAR_STORE((Gik + 1 * ldG), gik);
	tx  = SCALAR_LOAD((temp + 7 * blockDim.x + threadIdx.x));
	wg  = SCALAR_LOAD((weights + p_outer + p_inner));

	xik = SCALAR_LOAD((Xik + 1 * ldX));
	gik = SCALAR_LOAD((Gik + 2 * ldG));

	tx = SCALAR_MUL(tx, wg);
	gik = SCALAR_FMA(tx, xik, gik);
	SCALAR_STORE((Gik + 2 * ldG), gik);
	tx  = SCALAR_LOAD((temp + 5 * blockDim.x + threadIdx.x));
	wg  = SCALAR_LOAD((weights + p_outer + p_inner));

	xik = SCALAR_LOAD((Xik + 2 * ldX));
	gik = SCALAR_LOAD((Gik + 0 * ldG));

	tx = SCALAR_MUL(tx, wg);
	gik = SCALAR_FMA(tx, xik, gik);
	SCALAR_STORE((Gik + 0 * ldG), gik);
	tx  = SCALAR_LOAD((temp + 7 * blockDim.x + threadIdx.x));
	wg  = SCALAR_LOAD((weights + p_outer + p_inner));

	xik = SCALAR_LOAD((Xik + 2 * ldX));
	gik = SCALAR_LOAD((Gik + 1 * ldG));

	tx = SCALAR_MUL(tx, wg);
	gik = SCALAR_FMA(tx, xik, gik);
	SCALAR_STORE((Gik + 1 * ldG), gik);
	tx  = SCALAR_LOAD((temp + 8 * blockDim.x + threadIdx.x));
	wg  = SCALAR_LOAD((weights + p_outer + p_inner));

	xik = SCALAR_LOAD((Xik + 2 * ldX));
	gik = SCALAR_LOAD((Gik + 2 * ldG));

	tx = SCALAR_MUL(tx, wg);
	gik = SCALAR_FMA(tx, xik, gik);
	SCALAR_STORE((Gik + 2 * ldG), gik);
      }
    }
  }

  __global__ void dev_integral_1(size_t npts,
				   double *points_x,
				   double *points_y,
				   double *points_z,
           shell_pair* sp,
				   double *Xi,
				   int ldX,
				   double *Gi,
				   int ldG, 
				   double *weights, 
				   double *boys_table) {
    dev_integral_1_driver( npts, points_x, points_y, points_z, sp, Xi, ldX,
      Gi, ldG, weights, boys_table );
  }

  void integral_1(size_t npts,
		  double *_points_x,	
		  double *_points_y,	
		  double *_points_z,	
      shell_pair* sp,
		  double *Xi,
		  int ldX,
		  double *Gi,
		  int ldG, 
		  double *weights,
		  double *boys_table,
      cudaStream_t stream) {
    dev_integral_1<<<320, 128, 0, stream>>>(npts,
				 _points_x,
				 _points_y,
				 _points_z,
         sp,
				 Xi,
				 ldX,
				 Gi,
				 ldG, 
				 weights, 
				 boys_table);
  }

  __global__ void dev_integral_1_batched(
           const GauXC::ShellPairToTaskDevice* sp2task,
           GauXC::XCDeviceTask*                device_tasks,
				   double *boys_table) {

    const int ntask = sp2task->ntask;
    for( int i_task = blockIdx.y; i_task < ntask; i_task += gridDim.y ) {
    
      const auto iT = sp2task->task_idx_device[i_task];
      const auto* task  = device_tasks + iT;
      const auto  npts  = task->npts;

      const auto  i_off = sp2task->task_shell_off_row_device[i_task]*npts;

      dev_integral_1_driver( 
        npts,
        task->points_x,
        task->points_y,
        task->points_z,
        sp2task->shell_pair_device,
        task->fmat + i_off,
        npts,
        task->gmat + i_off,
        npts,
        task->weights, boys_table );
    }

  }


  void integral_1_batched(size_t ntask_sp,
        const GauXC::ShellPairToTaskDevice* sp2task,
        GauXC::XCDeviceTask*                device_tasks,
		    double *boys_table,
        cudaStream_t stream) {

    int nthreads = 128;
    int nblocks_x = 160;
    int nblocks_y = ntask_sp;
    dim3 nblocks(nblocks_x, nblocks_y);
    dev_integral_1_batched<<<nblocks,nthreads,0,stream>>>(
      sp2task, device_tasks, boys_table );

  }

#define K1_WARPSIZE 32
#define K1_NUMWARPS 8
#define K1_NUMTHREADS (K1_WARPSIZE * K1_NUMWARPS)
#define K1_USESHARED 1
#define K1_MAX_PRIMPAIRS 32

__inline__ __device__ void dev_integral_1_task(
  const int i,
  const int npts,
  const int nprim_pairs,
  // Data
  double X_AB,
  double Y_AB,
  double Z_AB,
  // Point data
  double4 (&s_task_data)[K1_NUMTHREADS],
  // Shell Pair Data
  const shell_pair* sp,
  // Output Data
  const double *Xi,
  const double *Xj,
  int ldX,
  double *Gi,
  double *Gj,
  int ldG, 
  // Other
  double *boys_table) {

  const int laneId = threadIdx.x % K1_WARPSIZE;

  const auto& prim_pairs = sp->prim_pairs();

#if K1_USESHARED
  // Load Primpairs to shared
  __shared__ GauXC::PrimitivePair<double> s_prim_pairs[K1_NUMWARPS][K1_MAX_PRIMPAIRS];

  const int warpId = (threadIdx.x / K1_WARPSIZE);
  const int32_t* src = (int32_t*) &(prim_pairs[0]);
  int32_t* dst = (int32_t*) &(s_prim_pairs[warpId][0]);

  for (int i = laneId; i < nprim_pairs * sizeof(GauXC::PrimitivePair<double>) / sizeof(int32_t); i+=K1_WARPSIZE) {
    dst[i] = src[i]; 
  }
  __syncwarp();
#endif


  __shared__ double outBuffer[K1_NUMTHREADS][3];

  // Loop over points in shared in batches of 32
  for (int i = 0; i < K1_NUMTHREADS / K1_WARPSIZE; i++) {
    for (int j = 0; j < 3; j++) {
      outBuffer[threadIdx.x][j] = 0.0;
    }

    double temp_0, temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8;
    temp_0 = SCALAR_ZERO();
    temp_1 = SCALAR_ZERO();
    temp_2 = SCALAR_ZERO();
    temp_3 = SCALAR_ZERO();
    temp_4 = SCALAR_ZERO();
    temp_5 = SCALAR_ZERO();
    temp_6 = SCALAR_ZERO();
    temp_7 = SCALAR_ZERO();
    temp_8 = SCALAR_ZERO();

    const int pointIndex = i * K1_WARPSIZE + laneId;

    if (pointIndex < npts) {
      const double point_x = s_task_data[pointIndex].x;
      const double point_y = s_task_data[pointIndex].y;
      const double point_z = s_task_data[pointIndex].z;
      const double weight = s_task_data[pointIndex].w;

      for(int ij = 0; ij < nprim_pairs; ++ij) {
#if K1_USESHARED
        double RHO = s_prim_pairs[warpId][ij].gamma;
        double RHO_INV = s_prim_pairs[warpId][ij].gamma_inv;
        double X_PA = s_prim_pairs[warpId][ij].PA.x;
        double Y_PA = s_prim_pairs[warpId][ij].PA.y;
        double Z_PA = s_prim_pairs[warpId][ij].PA.z;

        double xP = s_prim_pairs[warpId][ij].P.x;
        double yP = s_prim_pairs[warpId][ij].P.y;
        double zP = s_prim_pairs[warpId][ij].P.z;

        double eval = s_prim_pairs[warpId][ij].K_coeff_prod;
#else
        double RHO = prim_pairs[ij].gamma;
        double RHO_INV = prim_pairs[ij].gamma_inv;
        double X_PA = prim_pairs[ij].PA.x;
        double Y_PA = prim_pairs[ij].PA.y;
        double Z_PA = prim_pairs[ij].PA.z;

        double xP = prim_pairs[ij].P.x;
        double yP = prim_pairs[ij].P.y;
        double zP = prim_pairs[ij].P.z;

        double eval = prim_pairs[ij].K_coeff_prod;
#endif

        // Evaluate T Values
        SCALAR_TYPE X_PC = SCALAR_SUB(xP, point_x);
        SCALAR_TYPE Y_PC = SCALAR_SUB(yP, point_y);
        SCALAR_TYPE Z_PC = SCALAR_SUB(zP, point_z);


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
        tx = temp_0;
        tx = SCALAR_ADD(tx, t10);
        temp_0 = tx;
        t20 = SCALAR_MUL(X_PA, t10);
        t20 = SCALAR_FNMA(X_PC, t11, t20);
        tx = SCALAR_SUB(t00, t01);
        ty = SCALAR_SET1(0.5 * 1);
        ty = SCALAR_MUL(ty, RHO_INV);
        t20 = SCALAR_FMA(tx, ty, t20);
        tx = temp_3;
        tx = SCALAR_ADD(tx, t20);
        temp_3 = tx;
        t20 = SCALAR_MUL(Y_PA, t10);
        t20 = SCALAR_FNMA(Y_PC, t11, t20);
        tx = temp_4;
        tx = SCALAR_ADD(tx, t20);
        temp_4 = tx;
        t20 = SCALAR_MUL(Z_PA, t10);
        t20 = SCALAR_FNMA(Z_PC, t11, t20);
        tx = temp_5;
        tx = SCALAR_ADD(tx, t20);
        temp_5 = tx;
        t10 = SCALAR_MUL(Y_PA, t00);
        t10 = SCALAR_FNMA(Y_PC, t01, t10);
        t11 = SCALAR_MUL(Y_PA, t01);
        t11 = SCALAR_FNMA(Y_PC, t02, t11);
        tx = temp_1;
        tx = SCALAR_ADD(tx, t10);
        temp_1 = tx;
        t20 = SCALAR_MUL(Y_PA, t10);
        t20 = SCALAR_FNMA(Y_PC, t11, t20);
        tx = SCALAR_SUB(t00, t01);
        ty = SCALAR_SET1(0.5 * 1);
        ty = SCALAR_MUL(ty, RHO_INV);
        t20 = SCALAR_FMA(tx, ty, t20);
        tx = temp_6;
        tx = SCALAR_ADD(tx, t20);
        temp_6 = tx;
        t20 = SCALAR_MUL(Z_PA, t10);
        t20 = SCALAR_FNMA(Z_PC, t11, t20);
        tx = temp_7;
        tx = SCALAR_ADD(tx, t20);
        temp_7 = tx;
        t10 = SCALAR_MUL(Z_PA, t00);
        t10 = SCALAR_FNMA(Z_PC, t01, t10);
        t11 = SCALAR_MUL(Z_PA, t01);
        t11 = SCALAR_FNMA(Z_PC, t02, t11);
        tx = temp_2;
        tx = SCALAR_ADD(tx, t10);
        temp_2 = tx;
        t20 = SCALAR_MUL(Z_PA, t10);
        t20 = SCALAR_FNMA(Z_PC, t11, t20);
        tx = SCALAR_SUB(t00, t01);
        ty = SCALAR_SET1(0.5 * 1);
        ty = SCALAR_MUL(ty, RHO_INV);
        t20 = SCALAR_FMA(tx, ty, t20);
        tx = temp_8;
        tx = SCALAR_ADD(tx, t20);
        temp_8 = tx;
      }

#if 0
      if (
        abs(temp_0) > 1e-12 || abs(temp_1) > 1e-12 || abs(temp_2) > 1e-12 ||
        abs(temp_3) > 1e-12 || abs(temp_4) > 1e-12 || abs(temp_5) > 1e-12 ||
        abs(temp_6) > 1e-12 || abs(temp_7) > 1e-12 || abs(temp_8) > 1e-12
      ) 
#endif
      {
        const double * __restrict__ Xik = (Xi + pointIndex);
        const double * __restrict__ Xjk = (Xj + pointIndex);
        double * __restrict__ Gik = (Gi + pointIndex);
        double * __restrict__ Gjk = (Gj + pointIndex);

        SCALAR_TYPE const_value_v = weight;

        double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k;
        SCALAR_TYPE const_value_w;
        SCALAR_TYPE tx, ty, tz, tw, t0, t1, t2;

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
        outBuffer[threadIdx.x][0] += tz;

        tx = SCALAR_LOAD((Xik + 1 * ldX));
        t1 = SCALAR_MUL(temp_4, const_value_w);
        tz = SCALAR_MUL(ty, t1);
        tw = SCALAR_FMA(tx, t1, tw);
        outBuffer[threadIdx.x][1] += tz;

        tx = SCALAR_LOAD((Xik + 2 * ldX));
        t2 = SCALAR_MUL(temp_5, const_value_w);
        tz = SCALAR_MUL(ty, t2);
        tw = SCALAR_FMA(tx, t2, tw);
        outBuffer[threadIdx.x][2] += tz;

        X_ABp = SCALAR_MUL(X_ABp, X_AB); comb_m_i = SCALAR_MUL(comb_m_i * 1, SCALAR_RECIPROCAL(1));
        Y_ABp = 1.0; comb_n_j = 1.0;
        Z_ABp = 1.0; comb_p_k = 1.0;
        const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
        const_value_w = SCALAR_MUL(const_value_v, const_value);

        tx = SCALAR_LOAD((Xik + 0 * ldX));
        t0 = SCALAR_MUL(temp_0, const_value_w);
        tz = SCALAR_MUL(ty, t0);
        tw = SCALAR_FMA(tx, t0, tw);
        outBuffer[threadIdx.x][0] += tz;

        tx = SCALAR_LOAD((Xik + 1 * ldX));
        t1 = SCALAR_MUL(temp_1, const_value_w);
        tz = SCALAR_MUL(ty, t1);
        tw = SCALAR_FMA(tx, t1, tw);
        outBuffer[threadIdx.x][1] += tz;

        tx = SCALAR_LOAD((Xik + 2 * ldX));
        t2 = SCALAR_MUL(temp_2, const_value_w);
        tz = SCALAR_MUL(ty, t2);
        tw = SCALAR_FMA(tx, t2, tw);
        outBuffer[threadIdx.x][2] += tz;
        //atomicAdd((Gjk + 0 * ldG), tw);

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
        outBuffer[threadIdx.x][0] += tz;

        tx = SCALAR_LOAD((Xik + 1 * ldX));
        t1 = SCALAR_MUL(temp_6, const_value_w);
        tz = SCALAR_MUL(ty, t1);
        tw = SCALAR_FMA(tx, t1, tw);
        outBuffer[threadIdx.x][1] += tz;

        tx = SCALAR_LOAD((Xik + 2 * ldX));
        t2 = SCALAR_MUL(temp_7, const_value_w);
        tz = SCALAR_MUL(ty, t2);
        tw = SCALAR_FMA(tx, t2, tw);
        outBuffer[threadIdx.x][2] += tz;

        Y_ABp = SCALAR_MUL(Y_ABp, Y_AB); comb_n_j = SCALAR_MUL(comb_n_j * 1, SCALAR_RECIPROCAL(1));
        Z_ABp = 1.0; comb_p_k = 1.0;
        const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
        const_value_w = SCALAR_MUL(const_value_v, const_value);

        tx = SCALAR_LOAD((Xik + 0 * ldX));
        t0 = SCALAR_MUL(temp_0, const_value_w);
        tz = SCALAR_MUL(ty, t0);
        tw = SCALAR_FMA(tx, t0, tw);
        outBuffer[threadIdx.x][0] += tz;

        tx = SCALAR_LOAD((Xik + 1 * ldX));
        t1 = SCALAR_MUL(temp_1, const_value_w);
        tz = SCALAR_MUL(ty, t1);
        tw = SCALAR_FMA(tx, t1, tw);
        outBuffer[threadIdx.x][1] += tz;

        tx = SCALAR_LOAD((Xik + 2 * ldX));
        t2 = SCALAR_MUL(temp_2, const_value_w);
        tz = SCALAR_MUL(ty, t2);
        tw = SCALAR_FMA(tx, t2, tw);
        outBuffer[threadIdx.x][2] += tz;
        //atomicAdd((Gjk + 1 * ldG), tw);

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
        outBuffer[threadIdx.x][0] += tz;

        tx = SCALAR_LOAD((Xik + 1 * ldX));
        t1 = SCALAR_MUL(temp_7, const_value_w);
        tz = SCALAR_MUL(ty, t1);
        tw = SCALAR_FMA(tx, t1, tw);
        outBuffer[threadIdx.x][1] += tz;

        tx = SCALAR_LOAD((Xik + 2 * ldX));
        t2 = SCALAR_MUL(temp_8, const_value_w);
        tz = SCALAR_MUL(ty, t2);
        tw = SCALAR_FMA(tx, t2, tw);
        outBuffer[threadIdx.x][2] += tz;

        Z_ABp = SCALAR_MUL(Z_ABp, Z_AB); comb_p_k = SCALAR_MUL(comb_p_k * 1, SCALAR_RECIPROCAL(1));
        const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
        const_value_w = SCALAR_MUL(const_value_v, const_value);

        tx = SCALAR_LOAD((Xik + 0 * ldX));
        t0 = SCALAR_MUL(temp_0, const_value_w);
        tz = SCALAR_MUL(ty, t0);
        tw = SCALAR_FMA(tx, t0, tw);
        outBuffer[threadIdx.x][0] += tz;

        tx = SCALAR_LOAD((Xik + 1 * ldX));
        t1 = SCALAR_MUL(temp_1, const_value_w);
        tz = SCALAR_MUL(ty, t1);
        tw = SCALAR_FMA(tx, t1, tw);
        outBuffer[threadIdx.x][1] += tz;

        tx = SCALAR_LOAD((Xik + 2 * ldX));
        t2 = SCALAR_MUL(temp_2, const_value_w);
        tz = SCALAR_MUL(ty, t2);
        tw = SCALAR_FMA(tx, t2, tw);
        outBuffer[threadIdx.x][2] += tz;
        //atomicAdd((Gjk + 2 * ldG), tw);

        atomicAdd((Gik + 0 * ldG), outBuffer[threadIdx.x][0]);
        atomicAdd((Gik + 1 * ldG), outBuffer[threadIdx.x][1]);
        atomicAdd((Gik + 2 * ldG), outBuffer[threadIdx.x][2]);

      }
    }
  }
  __syncwarp();
}

__global__ void 
__launch_bounds__(K1_NUMTHREADS, 1)
dev_integral_1_task_batched(
  int ntask, int nsubtask,
  GauXC::XCDeviceTask*                device_tasks,
  const GauXC::TaskToShellPairDevice* task2sp,
  const int4* subtasks,
  const int32_t* nprim_pairs_device,
  shell_pair** sp_ptr_device,
  double* sp_X_AB_device,
  double* sp_Y_AB_device,
  double* sp_Z_AB_device,
  double *boys_table) {

  __shared__ double4 s_task_data[K1_NUMTHREADS];

  const int warpId = threadIdx.x / K1_WARPSIZE;
  
  const int i_subtask = blockIdx.x;
  const int i_task = subtasks[i_subtask].x;
  const int point_start = subtasks[i_subtask].y;
  const int point_end = subtasks[i_subtask].z;
  const int point_count = point_end - point_start;

  const auto* task = device_tasks + i_task;

  const int npts = task->npts;

  const auto* points_x = task->points_x;
  const auto* points_y = task->points_y;
  const auto* points_z = task->points_z;
  const auto* weights = task->weights;

  const auto nsp = task2sp[i_task].nsp;

  const int npts_block = (point_count + blockDim.x - 1) / blockDim.x;

  for (int i_block = 0; i_block < npts_block; i_block++) {
    const int i = point_start + i_block * blockDim.x;

    // load point into registers
    const double point_x = points_x[i + threadIdx.x];
    const double point_y = points_y[i + threadIdx.x];
    const double point_z = points_z[i + threadIdx.x];
    const double weight = weights[i + threadIdx.x];

    s_task_data[threadIdx.x].x = point_x;
    s_task_data[threadIdx.x].y = point_y;
    s_task_data[threadIdx.x].z = point_z;
    s_task_data[threadIdx.x].w = weight;
    __syncthreads();

    for (int j = K1_NUMWARPS*blockIdx.y+warpId; j < nsp; j+=K1_NUMWARPS*gridDim.y) {
      const auto i_off = task2sp[i_task].task_shell_off_row_device[j];

      const auto index =  task2sp[i_task].shell_pair_linear_idx_device[j];
      const auto* sp = sp_ptr_device[index];
      const auto nprim_pairs = nprim_pairs_device[index];
      const double X_AB = 0.0;
      const double Y_AB = 0.0;
      const double Z_AB = 0.0;

      dev_integral_1_task(
        i, point_count, nprim_pairs,
        X_AB, Y_AB, Z_AB,
        s_task_data,
        sp,
        task->fmat + i_off + i,
        task->fmat + i_off + i,
        npts,
        task->gmat + i_off + i,
        task->gmat + i_off + i,
        npts,
        boys_table);
    }
    __syncthreads();
  }
}

  void integral_1_task_batched(
    size_t ntasks, size_t nsubtask,
    int max_primpairs, size_t max_nsp,
    GauXC::XCDeviceTask*                device_tasks,
    const GauXC::TaskToShellPairDevice* task2sp,
    const std::array<int32_t, 4>*  subtasks,
    const int32_t* nprim_pairs_device,
    shell_pair** sp_ptr_device,
    double* sp_X_AB_device,
    double* sp_Y_AB_device,
    double* sp_Z_AB_device,
    double *boys_table,
    cudaStream_t stream) {

    size_t xy_max = (1ul << 16) - 1;
    int nthreads = K1_NUMTHREADS;
    int nblocks_x = nsubtask;
    int nblocks_y = 8; //std::min(max_nsp,  xy_max);
    int nblocks_z = 1;
    dim3 nblocks(nblocks_x, nblocks_y, nblocks_z);

    dev_integral_1_task_batched<<<nblocks, nthreads, 0, stream>>>(
      ntasks, nsubtask, 
      device_tasks, task2sp, 
      (int4*) subtasks, nprim_pairs_device, sp_ptr_device,
      sp_X_AB_device, sp_Y_AB_device, sp_Z_AB_device,
      boys_table );
  }
}
