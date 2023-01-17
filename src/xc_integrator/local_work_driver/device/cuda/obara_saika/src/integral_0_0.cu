#include <math.h>
#include "../include/gpu/chebyshev_boys_computation.hpp"
#include "config_obara_saika.hpp"
#include "integral_0_0.hu"

namespace XGPU {
  __inline__ __device__ void dev_integral_0_0_driver(size_t npts, 
				   const double *points_x,
				   const double *points_y,
				   const double *points_z,
           const shell_pair* sp,
				   const double *Xi,
				   const double *Xj,
				   int ldX,
				   double *Gi,
				   double *Gj,
				   int ldG, 
				   const double *weights, 
				   const double *boys_table) {

    double temp;

    // Load PrimPairs into shared mem
    const int nprim_pairs = sp->nprim_pairs();
    #if 1
    __shared__ GauXC::PrimitivePair<double> prim_pairs[GauXC::detail::nprim_pair_max];
    __syncthreads();
    if( threadIdx.x < 32 ) {
      const auto pp = sp->prim_pairs();
      for(int ij = threadIdx.x; ij < nprim_pairs; ij += 32) {
        prim_pairs[ij] = pp[ij];
      }
    }
    __syncthreads();
    #else
    const auto& prim_pairs = sp->prim_pairs();
    #endif

    const int npts_int = (int) npts;

    #pragma unroll(1)
    for(int p_outer = blockIdx.x * 128; p_outer < npts_int; p_outer += gridDim.x * 128) {

      const double * __restrict__ _point_outer_x = (points_x + p_outer);
      const double * __restrict__ _point_outer_y = (points_y + p_outer);
      const double * __restrict__ _point_outer_z = (points_z + p_outer);

      int p_inner = threadIdx.x;
      if (threadIdx.x < npts_int - p_outer) {

      temp = SCALAR_ZERO();
	    const SCALAR_TYPE xC = SCALAR_LOAD((_point_outer_x + p_inner));
	    const SCALAR_TYPE yC = SCALAR_LOAD((_point_outer_y + p_inner));
	    const SCALAR_TYPE zC = SCALAR_LOAD((_point_outer_z + p_inner));

      for(int ij = 0; ij < nprim_pairs; ++ij) {
        double RHO = prim_pairs[ij].gamma;
      
        double xP = prim_pairs[ij].P.x;
        double yP = prim_pairs[ij].P.y;
        double zP = prim_pairs[ij].P.z;
      
        double eval = prim_pairs[ij].K_coeff_prod;
      
        // Evaluate T Values
        const SCALAR_TYPE X_PC = SCALAR_SUB(xP, xC);
        const SCALAR_TYPE Y_PC = SCALAR_SUB(yP, yC);
        const SCALAR_TYPE Z_PC = SCALAR_SUB(zP, zC);
      
        SCALAR_TYPE TVAL = SCALAR_MUL(X_PC, X_PC);
        TVAL = SCALAR_FMA(Y_PC, Y_PC, TVAL);
        TVAL = SCALAR_FMA(Z_PC, Z_PC, TVAL);
        TVAL = SCALAR_MUL(RHO, TVAL);
      
        // Evaluate VRR Buffer
        const SCALAR_TYPE t00 = boys_element_0(TVAL);
        temp = SCALAR_FMA( eval, t00, temp );
      }
      if (abs(temp) > 1e-12) {
        const double * __restrict__ Xik = (Xi + p_outer + p_inner);
        const double * __restrict__ Xjk = (Xj + p_outer + p_inner);
        double * __restrict__ Gik = (Gi + p_outer + p_inner);
        double * __restrict__ Gjk = (Gj + p_outer + p_inner);
      
        SCALAR_TYPE const_value_v = SCALAR_LOAD((weights + p_outer + p_inner));
      
        double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k;
        SCALAR_TYPE const_value_w;
        SCALAR_TYPE tx, ty, tz, tw, t0;
      
        X_ABp = 1.0; comb_m_i = 1.0;
        Y_ABp = 1.0; comb_n_j = 1.0;
        Z_ABp = 1.0; comb_p_k = 1.0;
        const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
        const_value_w = SCALAR_MUL(const_value_v, const_value);
        tx = SCALAR_LOAD(Xik);
        ty = SCALAR_LOAD(Xjk);
        t0 = SCALAR_MUL(temp, const_value_w);
        tz = SCALAR_MUL(ty, t0);
        tw = SCALAR_MUL(tx, t0);
        atomicAdd(Gik, tz);
        atomicAdd(Gjk, tw);
      }
      }
    }
  }





  __global__ void dev_integral_0_0(size_t npts,
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
    dev_integral_0_0_driver( npts, points_x, points_y, points_z, sp, Xi, Xj, ldX,
      Gi, Gj, ldG, weights, boys_table );
  }



  void integral_0_0(size_t npts,
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
    int nthreads = 128;
    int nblocks = std::min(intmax_t(320), GauXC::util::div_ceil(npts,nthreads));
    dev_integral_0_0<<<nblocks, nthreads,0,stream>>>(npts,
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





  __inline__ __device__ void dev_integral_0_0_batched_driver(
           const GauXC::ShellPairToTaskDevice* sp2task,
           GauXC::XCDeviceTask*                device_tasks,
				   double *boys_table) {

    //if (sp2task->shell_pair_device->nprim_pairs() == 0) return;
    const int ntask = sp2task->ntask;

    for( int i_task = blockIdx.y; i_task < ntask; i_task += gridDim.y ) {
    
      const auto iT = sp2task->task_idx_device[i_task];
      const auto* task  = device_tasks + iT;
      const auto  npts  = task->npts;

      const auto  i_off = sp2task->task_shell_off_row_device[i_task]*npts;
      const auto  j_off = sp2task->task_shell_off_col_device[i_task]*npts;

      dev_integral_0_0_driver( 
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

  __global__ void dev_integral_0_0_batched(
           const GauXC::ShellPairToTaskDevice* sp2task,
           GauXC::XCDeviceTask*                device_tasks,
				   double *boys_table) {
    dev_integral_0_0_batched_driver( sp2task, device_tasks, boys_table );
  }

  void integral_0_0_batched(size_t ntask_sp,
        const GauXC::ShellPairToTaskDevice* sp2task,
        GauXC::XCDeviceTask*                device_tasks,
		    double *boys_table,
        cudaStream_t stream) {

    int nthreads = 128;
    int nblocks_x = 160;
    int nblocks_y = ntask_sp;
    dim3 nblocks(nblocks_x, nblocks_y);
    dev_integral_0_0_batched<<<nblocks,nthreads,0,stream>>>(
      sp2task, device_tasks, boys_table );

  }





  __inline__ __device__ void dev_integral_0_0_soa_batched_driver(
           int32_t                         ntask,
           const int32_t*                  sp2task_idx_device,
           const int32_t*                  sp2task_shell_off_row_device,
           const int32_t*                  sp2task_shell_off_col_device,
           const GauXC::ShellPair<double>* shell_pair_device,
           const int32_t*                  task_npts,
           const double**                  task_points_x,
           const double**                  task_points_y,
           const double**                  task_points_z,
           const double**                  task_weights,
           const double**                  task_fmat,
           double**                        task_gmat,
				   double *                        boys_table) {

    for( int i_task = blockIdx.y; i_task < ntask; i_task += gridDim.y ) {
    
      const auto iT   = sp2task_idx_device[i_task];
      const auto npts = task_npts[iT];

      const auto  i_off = sp2task_shell_off_row_device[i_task] * npts;
      const auto  j_off = sp2task_shell_off_col_device[i_task] * npts;

      dev_integral_0_0_driver( 
        npts,
        task_points_x[iT],
        task_points_y[iT],
        task_points_z[iT],
        shell_pair_device,
        task_fmat[iT] + i_off,
        task_fmat[iT] + j_off,
        npts,
        task_gmat[iT] + i_off,
        task_gmat[iT] + j_off,
        npts,
        task_weights[iT], boys_table );
    }

  }

  __global__ void dev_integral_0_0_soa_batched(
           int32_t                         ntask,
           const int32_t*                  sp2task_idx_device,
           const int32_t*                  sp2task_shell_off_row_device,
           const int32_t*                  sp2task_shell_off_col_device,
           const GauXC::ShellPair<double>* shell_pair_device,
           const int32_t*                  task_npts,
           const double**                   task_points_x,
           const double**                   task_points_y,
           const double**                   task_points_z,
           const double**                   task_weights,
           const double**                   task_fmat,
           double**                         task_gmat,
				   double *boys_table) {
    dev_integral_0_0_soa_batched_driver( ntask, sp2task_idx_device, 
      sp2task_shell_off_row_device, sp2task_shell_off_col_device, shell_pair_device,
      task_npts, task_points_x, task_points_y, task_points_z, task_weights,
      task_fmat, task_gmat, boys_table );
  }


  __global__ void 
  __launch_bounds__(128, 16)
  dev_integral_0_0_shell_batched(
           int nsp,
           const GauXC::ShellPairToTaskDevice* sp2task,
           GauXC::XCDeviceTask*                device_tasks,
				   double *boys_table) {

    for( int i = blockIdx.z; i < nsp; i += gridDim.z ) {
      dev_integral_0_0_batched_driver( sp2task + i, device_tasks, boys_table );
    }

  }

  void integral_0_0_shell_batched(
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

    dev_integral_0_0_shell_batched<<<nblocks,nthreads,0,stream>>>(
      nsp, sp2task, device_tasks, boys_table );

  }
   
#define K00_WARPSIZE 32
#define K00_NUMWARPS 8
#define K00_NUMTHREADS (K00_WARPSIZE * K00_NUMWARPS)
#define K00_USESHARED 1
#define K00_MAX_PRIMPAIRS 32

__inline__ __device__ void dev_integral_0_0_task(
  const int i,
  const int npts,
  const int nprim_pairs,
  // Point data
  double4 (&s_task_data)[K00_NUMTHREADS],
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
  const double *boys_table) {


  const auto& prim_pairs = sp->prim_pairs();

  const int laneId = threadIdx.x % K00_WARPSIZE;

#if K00_USESHARED
  // Load Primpairs to shared
  __shared__ GauXC::PrimitivePair<double> s_prim_pairs[K00_NUMWARPS][K00_MAX_PRIMPAIRS];

  const int warpId = (threadIdx.x / K00_WARPSIZE);
  const int32_t* src = (int32_t*) &(prim_pairs[0]);
  int32_t* dst = (int32_t*) &(s_prim_pairs[warpId][0]);

  for (int i = laneId; i < nprim_pairs * sizeof(GauXC::PrimitivePair<double>) / sizeof(int32_t); i+=K00_WARPSIZE) {
    dst[i] = src[i]; 
  }
  __syncwarp();
#endif

  // Loop over points in shared in batches of 32
  for (int i = 0; i < K00_NUMTHREADS / K00_WARPSIZE; i++) {
    double temp = SCALAR_ZERO();

    const int pointIndex = i * K00_WARPSIZE + laneId;

    if (pointIndex < npts) {

      const double point_x = s_task_data[pointIndex].x;
      const double point_y = s_task_data[pointIndex].y;
      const double point_z = s_task_data[pointIndex].z;
      const double weight = s_task_data[pointIndex].w;

      for (int ij = 0; ij < nprim_pairs; ij++) {
#if K00_USESHARED
        double RHO = s_prim_pairs[warpId][ij].gamma;
        double xP = s_prim_pairs[warpId][ij].P.x;
        double yP = s_prim_pairs[warpId][ij].P.y;
        double zP = s_prim_pairs[warpId][ij].P.z;
        double eval = s_prim_pairs[warpId][ij].K_coeff_prod;
#else
        double RHO = prim_pairs[ij].gamma;
        double xP = prim_pairs[ij].P.x;
        double yP = prim_pairs[ij].P.y;
        double zP = prim_pairs[ij].P.z;
        double eval = prim_pairs[ij].K_coeff_prod;
#endif
      
        // Evaluate T Values
        const SCALAR_TYPE X_PC = SCALAR_SUB(xP, point_x);
        const SCALAR_TYPE Y_PC = SCALAR_SUB(yP, point_y);
        const SCALAR_TYPE Z_PC = SCALAR_SUB(zP, point_z);
      
        SCALAR_TYPE TVAL = SCALAR_MUL(X_PC, X_PC);
        TVAL = SCALAR_FMA(Y_PC, Y_PC, TVAL);
        TVAL = SCALAR_FMA(Z_PC, Z_PC, TVAL);
        TVAL = SCALAR_MUL(RHO, TVAL);
      
        // Evaluate VRR Buffer
        const SCALAR_TYPE t00 = boys_element_0(TVAL);
        temp = SCALAR_FMA( eval, t00, temp );
      }

      // Output
      if (abs(temp) > 1e-12) {
        const double * __restrict__ Xik = (Xi + pointIndex);
        const double * __restrict__ Xjk = (Xj + pointIndex);
        double * __restrict__ Gik = (Gi + pointIndex);
        double * __restrict__ Gjk = (Gj + pointIndex);

        SCALAR_TYPE const_value_v = weight;
      
        double const_value, X_ABp, Y_ABp, Z_ABp, comb_m_i, comb_n_j, comb_p_k;
        SCALAR_TYPE const_value_w;
        SCALAR_TYPE tx, ty, tz, tw, t0;
      
        X_ABp = 1.0; comb_m_i = 1.0;
        Y_ABp = 1.0; comb_n_j = 1.0;
        Z_ABp = 1.0; comb_p_k = 1.0;
        const_value = comb_m_i * comb_n_j * comb_p_k * X_ABp * Y_ABp * Z_ABp;
        const_value_w = SCALAR_MUL(const_value_v, const_value);
        tx = SCALAR_LOAD(Xik);
        ty = SCALAR_LOAD(Xjk);
        t0 = SCALAR_MUL(temp, const_value_w);
        tz = SCALAR_MUL(ty, t0);
        tw = SCALAR_MUL(tx, t0);
        atomicAdd(Gik, tz);
        atomicAdd(Gjk, tw);
      }
    }
  }
  __syncwarp();
}




__global__ void 
__launch_bounds__(K00_NUMTHREADS, 1)
dev_integral_0_0_task_batched(
  int ntask, int nsubtask,
  GauXC::XCDeviceTask*                device_tasks,
  const GauXC::TaskToShellPairDevice* task2sp,
  const int4* subtasks,
  const int32_t* nprim_pairs_device,
  shell_pair** sp_ptr_device,
  double *boys_table) {

  __shared__ double4 s_task_data[K00_NUMTHREADS];

  const int warpId = threadIdx.x / K00_WARPSIZE;
  
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

    for (int j = K00_NUMWARPS*blockIdx.y+warpId; j < nsp; j+=K00_NUMWARPS*gridDim.y) {
      const auto i_off = task2sp[i_task].task_shell_off_row_device[j];
      const auto j_off = task2sp[i_task].task_shell_off_col_device[j];

      const auto index =  task2sp[i_task].shell_pair_linear_idx_device[j];
      const auto* sp = sp_ptr_device[index];
      const auto nprim_pairs = nprim_pairs_device[index];

      dev_integral_0_0_task(
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

  void integral_0_0_task_batched(
    size_t ntasks, size_t nsubtask,
    size_t max_nsp,
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
    int nthreads = K00_NUMTHREADS;
    int nblocks_x = nsubtask;
    int nblocks_y = 8; //std::min(max_nsp,  xy_max);
    int nblocks_z = 1;
    dim3 nblocks(nblocks_x, nblocks_y, nblocks_z);

    dev_integral_0_0_task_batched<<<nblocks, nthreads, 0, stream>>>(
      ntasks, nsubtask,
      device_tasks, task2sp, 
      (int4*) subtasks, nprim_pairs_device, sp_ptr_device,
      boys_table );

  }

}
