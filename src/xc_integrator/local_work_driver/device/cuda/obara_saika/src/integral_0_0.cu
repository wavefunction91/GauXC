#include <math.h>
#include "../include/gpu/chebyshev_boys_computation.hpp"
#include "config_obara_saika.hpp" /* KXI Most macros are here*/
#include "integral_0_0.hu"


namespace XGPU {
  __inline__ __device__ void dev_integral_0_0_driver(size_t npts, int init,
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
#if OPTIMIZE_ACCESS_0_0
    int block_tid = threadIdx.x + threadIdx.y * blockDim.x;
    int block_nthreads = blockDim.x * blockDim.y;
    __shared__ double prim_pairs_P_x[GauXC::detail::nprim_pair_max+1];
    __shared__ double prim_pairs_P_y[GauXC::detail::nprim_pair_max];
    __shared__ double prim_pairs_P_z[GauXC::detail::nprim_pair_max];
    __shared__ double prim_pairs_K_coeff_prod[GauXC::detail::nprim_pair_max];
    __shared__ double prim_pairs_gamma[GauXC::detail::nprim_pair_max];

    if(init) {
    GauXC::PrimitivePair<double> tmp;
    const auto pp = sp->prim_pairs();
    for(int ij = block_tid; ij < nprim_pairs; ij += block_nthreads) {
        tmp =  pp[ij];
        prim_pairs_P_x[ij] = tmp.P.x;
        prim_pairs_P_y[ij] = tmp.P.y;
        prim_pairs_P_z[ij] = tmp.P.z;
        prim_pairs_gamma[ij] = tmp.gamma;
        prim_pairs_K_coeff_prod[ij] = tmp.K_coeff_prod;
    }
    __syncthreads();
   }
#else
    #if 1
    __shared__ GauXC::PrimitivePair<double> prim_pairs[GauXC::detail::nprim_pair_max];
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
#endif    
    #pragma unroll(1)
  for(size_t p_outer = blockIdx.x * 128; p_outer < npts; p_outer += gridDim.x * 128) {
      const double * __restrict__ _point_outer_x = (points_x + p_outer);
      const double * __restrict__ _point_outer_y = (points_y + p_outer);
      const double * __restrict__ _point_outer_z = (points_z + p_outer);
      size_t p_inner = (threadIdx.x < (npts - p_outer)) ? threadIdx.x : (npts - p_outer);
      temp = SCALAR_ZERO();
      const SCALAR_TYPE xC = SCALAR_LOAD((_point_outer_x + p_inner));
      const SCALAR_TYPE yC = SCALAR_LOAD((_point_outer_y + p_inner));
      const SCALAR_TYPE zC = SCALAR_LOAD((_point_outer_z + p_inner));

      #if OPTIMIZE_ACCESS_0_0
      int ij = block_tid % nprim_pairs;
      for(int s = 0; s < nprim_pairs; ++s, ++ij) {
        ij = (ij != nprim_pairs) * ij;
        double RHO = prim_pairs_gamma[ij];

        double xP = prim_pairs_P_x[ij];
        double yP = prim_pairs_P_y[ij];
        double zP = prim_pairs_P_z[ij];

        double eval = prim_pairs_K_coeff_prod[ij];
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
     #else
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
     #endif

      if(threadIdx.x < npts - p_outer) {
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
    dev_integral_0_0_driver( npts, 1, points_x, points_y, points_z, sp, Xi, Xj, ldX,
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
    int init = 1;
    const int ntask = sp2task->ntask;
    for( int i_task = blockIdx.y; i_task < ntask; i_task += gridDim.y ) {
      const auto iT = sp2task->task_idx_device[i_task];
      const auto* task  = device_tasks + iT;
      const auto  npts  = task->npts;

      const auto  i_off = sp2task->task_shell_off_row_device[i_task]*npts;
      const auto  j_off = sp2task->task_shell_off_col_device[i_task]*npts;

      dev_integral_0_0_driver( 
        npts, init,
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
      init = 0;
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

    int init = 1;
    for( int i_task = blockIdx.y; i_task < ntask; i_task += gridDim.y ) {
      const auto iT   = sp2task_idx_device[i_task];
      const auto npts = task_npts[iT];

      const auto  i_off = sp2task_shell_off_row_device[i_task] * npts;
      const auto  j_off = sp2task_shell_off_col_device[i_task] * npts;

      dev_integral_0_0_driver( 
        npts, init,
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
      init = 0;
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

  __global__ void dev_integral_0_0_shell_batched(
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

    int nthreads = 128;
    int nblocks_x = 80;
    int nblocks_y = max_ntask/8;
    int nblocks_z = nsp/8;
    dim3 nblocks(nblocks_x, nblocks_y, nblocks_z);
    dev_integral_0_0_shell_batched<<<nblocks,nthreads,0,stream>>>(
      nsp, sp2task, device_tasks, boys_table );

  }
      
}
