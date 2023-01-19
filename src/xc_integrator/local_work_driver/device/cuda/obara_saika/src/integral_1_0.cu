#include <math.h>
#include "../include/gpu/chebyshev_boys_computation.hpp"
#include "config_obara_saika.hpp"
#include "integral_1_0.hu"

namespace XGPU {
  __inline__ __device__ void dev_integral_1_0_driver(size_t npts,
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
    //__shared__ double temp[128 * 3];
    double temp_0, temp_1, temp_2;
    const auto nprim_pairs = sp->nprim_pairs();
    const auto prim_pairs  = sp->prim_pairs();


    const int npts_int = (int) npts;
    
    for(int p_outer = blockIdx.x * blockDim.x; p_outer < npts_int; p_outer += gridDim.x * blockDim.x) {
      double *_point_outer_x = (_points_x + p_outer);
      double *_point_outer_y = (_points_y + p_outer);
      double *_point_outer_z = (_points_z + p_outer);

      int p_inner = threadIdx.x;
      if (threadIdx.x < npts_int - p_outer) {

      //for(int i = 0; i < 3; ++i) SCALAR_STORE((temp + i * blockDim.x + threadIdx.x), SCALAR_ZERO());
      temp_0 = SCALAR_ZERO();
      temp_1 = SCALAR_ZERO();
      temp_2 = SCALAR_ZERO();

      for(int ij = 0; ij < nprim_pairs; ++ij) {
	double RHO = prim_pairs[ij].gamma;
	double X_PA = prim_pairs[ij].PA.x;
	double Y_PA = prim_pairs[ij].PA.y;
	double Z_PA = prim_pairs[ij].PA.z;

	double xP = prim_pairs[ij].P.x;
	double yP = prim_pairs[ij].P.y;
	double zP = prim_pairs[ij].P.z;

	double eval = prim_pairs[ij].K_coeff_prod;

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
	
	SCALAR_TYPE t00, t01, TVAL_inv_e;

	// Evaluate Boys function
	boys_element<1>(&TVAL, &TVAL_inv_e, &t01, boys_table);

	// Evaluate VRR Buffer
	SCALAR_TYPE t10, tx;

	t00 = SCALAR_MUL(SCALAR_ADD(SCALAR_MUL(TVAL, t01), TVAL_inv_e), SCALAR_SET1(2.00000000000000000000));

	t00 = SCALAR_MUL(eval, t00);
	t01 = SCALAR_MUL(eval, t01);
	t10 = SCALAR_MUL(X_PA, t00);
	t10 = SCALAR_FNMA(X_PC, t01, t10);
	//tx = SCALAR_LOAD((temp + 0 * blockDim.x + threadIdx.x));
	tx = temp_0;
	tx = SCALAR_ADD(tx, t10);
	//SCALAR_STORE((temp + 0 * blockDim.x + threadIdx.x), tx);
  temp_0 = tx;
	t10 = SCALAR_MUL(Y_PA, t00);
	t10 = SCALAR_FNMA(Y_PC, t01, t10);
	//tx = SCALAR_LOAD((temp + 1 * blockDim.x + threadIdx.x));
	tx = temp_1;
	tx = SCALAR_ADD(tx, t10);
	//SCALAR_STORE((temp + 1 * blockDim.x + threadIdx.x), tx);
  temp_1 = tx;
	t10 = SCALAR_MUL(Z_PA, t00);
	t10 = SCALAR_FNMA(Z_PC, t01, t10);
	//tx = SCALAR_LOAD((temp + 2 * blockDim.x + threadIdx.x));
	tx = temp_2;
	tx = SCALAR_ADD(tx, t10);
	//SCALAR_STORE((temp + 2 * blockDim.x + threadIdx.x), tx);
  temp_2 = tx;
      }

  
      if (abs(temp_0) > 1e-12 || abs(temp_1) > 1e-12 || abs(temp_2) > 1e-12) {
	double *Xik = (Xi + p_outer + p_inner);
	double *Xjk = (Xj + p_outer + p_inner);
	double *Gik = (Gi + p_outer + p_inner);
	double *Gjk = (Gj + p_outer + p_inner);

	SCALAR_TYPE const_value_v = SCALAR_LOAD((weights + p_outer + p_inner));

	double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k;
	SCALAR_TYPE const_value_w;
	SCALAR_TYPE tx, ty, tz, tw, t0, t1, t2;

	X_ABp = 1.0; comb_m_i = 1.0;
	Y_ABp = 1.0; comb_n_j = 1.0;
	Z_ABp = 1.0; comb_p_k = 1.0;
	const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
	const_value_w = SCALAR_MUL(const_value_v, const_value);

  #if 0
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
  #else
	tx = SCALAR_LOAD((Xik + 0 * ldX));
	ty = SCALAR_LOAD((Xjk + 0 * ldX));

	t0 = SCALAR_MUL(temp_0, const_value_w);
	tz = SCALAR_MUL(ty, t0);
	tw = SCALAR_MUL(tx, t0);
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
  #endif
      }
      }
    }
  }

  __global__ void dev_integral_1_0(size_t npts,
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
    dev_integral_1_0_driver( npts, points_x, points_y, points_z, sp, Xi, Xj, ldX,
      Gi, Gj, ldG, weights, boys_table );
  }

    void integral_1_0(size_t npts,
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
      dev_integral_1_0<<<320, 128, 0, stream>>>(npts,
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




  template <bool swap>
  __inline__ __device__ void dev_integral_1_0_batched_driver(
           const GauXC::ShellPairToTaskDevice* sp2task,
           GauXC::XCDeviceTask*                device_tasks,
				   double *boys_table) {

    //if (sp2task->shell_pair_device->nprim_pairs() == 0) return;
    const int ntask = sp2task->ntask;
    #pragma unroll 1
    for( int i_task = blockIdx.y; i_task < ntask; i_task += gridDim.y ) {
    
      const auto iT = sp2task->task_idx_device[i_task];
      const auto* task  = device_tasks + iT;
      const auto  npts  = task->npts;

      int i_off, j_off;
      if constexpr ( swap ) {
        j_off = sp2task->task_shell_off_row_device[i_task]*npts;
        i_off = sp2task->task_shell_off_col_device[i_task]*npts;
      } else {
        i_off = sp2task->task_shell_off_row_device[i_task]*npts;
        j_off = sp2task->task_shell_off_col_device[i_task]*npts;
      }


      dev_integral_1_0_driver( 
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

  template <bool swap>
  __global__ void dev_integral_1_0_batched(
           const GauXC::ShellPairToTaskDevice* sp2task,
           GauXC::XCDeviceTask*                device_tasks,
				   double *boys_table) {
    dev_integral_1_0_batched_driver<swap>(sp2task,device_tasks,boys_table);
  }


  void integral_1_0_batched(bool swap, size_t ntask_sp,
        const GauXC::ShellPairToTaskDevice* sp2task,
        GauXC::XCDeviceTask*                device_tasks,
		    double *boys_table,
        cudaStream_t stream) {

    int nthreads = 128;
    int nblocks_x = 160;
    int nblocks_y = ntask_sp;
    dim3 nblocks(nblocks_x, nblocks_y);

    if(swap)
      dev_integral_1_0_batched<true><<<nblocks,nthreads,0,stream>>>(
        sp2task, device_tasks, boys_table );
    else
      dev_integral_1_0_batched<false><<<nblocks,nthreads,0,stream>>>(
        sp2task, device_tasks, boys_table );

  }

  template <bool swap>
  __global__ void dev_integral_1_0_shell_batched(
           int nsp,
           const GauXC::ShellPairToTaskDevice* sp2task,
           GauXC::XCDeviceTask*                device_tasks,
				   double *boys_table) {
    for(int i = blockIdx.z; i < nsp; i+= gridDim.z ) {
      dev_integral_1_0_batched_driver<swap>(sp2task+i,device_tasks,boys_table);
    }
  }

  void integral_1_0_shell_batched(
        bool swap,
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
    if(swap)
      dev_integral_1_0_shell_batched<true><<<nblocks,nthreads,0,stream>>>(
        nsp, sp2task, device_tasks, boys_table );
    else
      dev_integral_1_0_shell_batched<false><<<nblocks,nthreads,0,stream>>>(
        nsp, sp2task, device_tasks, boys_table );

  }

#define K10_WARPSIZE 32
#define K10_NUMWARPS 8
#define K10_NUMTHREADS (K10_WARPSIZE * K10_NUMWARPS)
#define K10_USESHARED 1
#define K10_MAX_PRIMPAIRS 32

__inline__ __device__ void dev_integral_1_0_task(
  const int i,
  const int npts,
  const int nprim_pairs,
  // Point data
  double4 (&s_task_data)[K10_NUMTHREADS],
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

  const int laneId = threadIdx.x % K10_WARPSIZE;

  const auto& prim_pairs = sp->prim_pairs();

#if K10_USESHARED
  // Load Primpairs to shared
  __shared__ GauXC::PrimitivePair<double> s_prim_pairs[K10_NUMWARPS][K10_MAX_PRIMPAIRS];

  const int warpId = (threadIdx.x / K10_WARPSIZE);
  const int32_t* src = (int32_t*) &(prim_pairs[0]);
  int32_t* dst = (int32_t*) &(s_prim_pairs[warpId][0]);

  for (int i = laneId; i < nprim_pairs * sizeof(GauXC::PrimitivePair<double>) / sizeof(int32_t); i+=K10_WARPSIZE) {
    dst[i] = src[i]; 
  }
  __syncwarp();
#endif


  // Loop over points in shared in batches of 32
  for (int i = 0; i < K10_NUMTHREADS / K10_WARPSIZE; i++) {
    double temp_0 = SCALAR_ZERO();
    double temp_1 = SCALAR_ZERO();
    double temp_2 = SCALAR_ZERO();

    const int pointIndex = i * K10_WARPSIZE + laneId;

    if (pointIndex < npts) {
      const double point_x = s_task_data[pointIndex].x;
      const double point_y = s_task_data[pointIndex].y;
      const double point_z = s_task_data[pointIndex].z;
      const double weight = s_task_data[pointIndex].w;

      for(int ij = 0; ij < nprim_pairs; ++ij) {
#if K10_USESHARED
        double RHO = s_prim_pairs[warpId][ij].gamma;
        double X_PA = s_prim_pairs[warpId][ij].PA.x;
        double Y_PA = s_prim_pairs[warpId][ij].PA.y;
        double Z_PA = s_prim_pairs[warpId][ij].PA.z;

        double xP = s_prim_pairs[warpId][ij].P.x;
        double yP = s_prim_pairs[warpId][ij].P.y;
        double zP = s_prim_pairs[warpId][ij].P.z;

        double eval = s_prim_pairs[warpId][ij].K_coeff_prod;
#else
        double RHO = prim_pairs[ij].gamma;
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
        
        SCALAR_TYPE t00, t01, TVAL_inv_e;

        // Evaluate Boys function
        boys_element<1>(&TVAL, &TVAL_inv_e, &t01, boys_table);

        // Evaluate VRR Buffer
        SCALAR_TYPE t10, tx;

        t00 = SCALAR_MUL(SCALAR_ADD(SCALAR_MUL(TVAL, t01), TVAL_inv_e), SCALAR_SET1(2.00000000000000000000));

        t00 = SCALAR_MUL(eval, t00);
        t01 = SCALAR_MUL(eval, t01);
        t10 = SCALAR_MUL(X_PA, t00);
        t10 = SCALAR_FNMA(X_PC, t01, t10);
        tx = temp_0;
        tx = SCALAR_ADD(tx, t10);
        temp_0 = tx;
        t10 = SCALAR_MUL(Y_PA, t00);
        t10 = SCALAR_FNMA(Y_PC, t01, t10);
        tx = temp_1;
        tx = SCALAR_ADD(tx, t10);
        temp_1 = tx;
        t10 = SCALAR_MUL(Z_PA, t00);
        t10 = SCALAR_FNMA(Z_PC, t01, t10);
        tx = temp_2;
        tx = SCALAR_ADD(tx, t10);
        temp_2 = tx;
      }

  
      if (abs(temp_0) > 1e-12 || abs(temp_1) > 1e-12 || abs(temp_2) > 1e-12) {
        const double * __restrict__ Xik = (Xi + pointIndex);
        const double * __restrict__ Xjk = (Xj + pointIndex);
        double * __restrict__ Gik = (Gi + pointIndex);
        double * __restrict__ Gjk = (Gj + pointIndex);

        SCALAR_TYPE const_value_v = weight;

        double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k;
        SCALAR_TYPE const_value_w;
        SCALAR_TYPE tx, ty, tz, tw, t0, t1, t2;

        X_ABp = 1.0; comb_m_i = 1.0;
        Y_ABp = 1.0; comb_n_j = 1.0;
        Z_ABp = 1.0; comb_p_k = 1.0;
        const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
        const_value_w = SCALAR_MUL(const_value_v, const_value);

        tx = SCALAR_LOAD((Xik + 0 * ldX));
        ty = SCALAR_LOAD((Xjk + 0 * ldX));

        t0 = SCALAR_MUL(temp_0, const_value_w);
        tz = SCALAR_MUL(ty, t0);
        tw = SCALAR_MUL(tx, t0);
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
      }
    }
  }
  __syncwarp();
}


template <bool swap>
__global__ void 
__launch_bounds__(K10_NUMTHREADS, 1)
dev_integral_1_0_task_batched(
  int ntask, int nsubtask,
  GauXC::XCDeviceTask*                device_tasks,
  const GauXC::TaskToShellPairDevice* task2sp,
  const int4* subtasks,
  const int32_t* nprim_pairs_device,
  shell_pair** sp_ptr_device,
  double *boys_table) {

  __shared__ double4 s_task_data[K10_NUMTHREADS];

  const int warpId = threadIdx.x / K10_WARPSIZE;
  
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

    for (int j = K10_NUMWARPS*blockIdx.y+warpId; j < nsp; j+=K10_NUMWARPS*gridDim.y) {
      const auto i_off = swap ? task2sp[i_task].task_shell_off_col_device[j] :
                                task2sp[i_task].task_shell_off_row_device[j];
      const auto j_off = swap ? task2sp[i_task].task_shell_off_row_device[j] :
                                task2sp[i_task].task_shell_off_col_device[j];


      const auto index =  task2sp[i_task].shell_pair_linear_idx_device[j];
      const auto* sp = sp_ptr_device[index];
      const auto nprim_pairs = nprim_pairs_device[index];

      dev_integral_1_0_task(
        i, point_count, nprim_pairs,
        s_task_data,
        sp,
        task->fmat + i_off + i,
        task->fmat + j_off + i,
        npts,
        task->gmat + i_off + i,
        task->gmat + j_off + i,
        npts,
        boys_table);
    }
    __syncthreads();
  }
}

  void integral_1_0_task_batched(
    bool swap,
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
    int nthreads = K10_NUMTHREADS;
    int nblocks_x = nsubtask;
    int nblocks_y = 8; //std::min(max_nsp,  xy_max);
    int nblocks_z = 1;
    dim3 nblocks(nblocks_x, nblocks_y, nblocks_z);

    if (swap) {
      dev_integral_1_0_task_batched<true><<<nblocks, nthreads, 0, stream>>>(
        ntasks, nsubtask, 
        device_tasks, task2sp, 
        (int4*) subtasks, nprim_pairs_device, sp_ptr_device,
        boys_table );
    } else {
      dev_integral_1_0_task_batched<false><<<nblocks, nthreads, 0, stream>>>(
        ntasks, nsubtask, 
        device_tasks, task2sp, 
        (int4*) subtasks, nprim_pairs_device, sp_ptr_device,
        boys_table );
    }
  }
}
