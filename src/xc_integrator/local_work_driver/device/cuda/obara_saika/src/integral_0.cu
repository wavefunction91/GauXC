#include <math.h>
#include "../include/gpu/chebyshev_boys_computation.hpp"
#include "config_obara_saika.hpp"
#include "integral_0.hu"

#include "device_specific/cuda_device_constants.hpp"

namespace XGPU {

using namespace GauXC;

  __inline__ __device__ void dev_integral_0_driver(size_t npts,
				 double *_points_x,
				 double *_points_y,
				 double *_points_z,
         shell_pair* sp,
				 double *Xi,
				 int ldX,
				 double *Gi,
				 int ldG, 
				 double *weights,
				 double *boys_table) {
    __shared__ double temp[128 * 1];
    const auto nprim_pairs = sp->nprim_pairs();
    const auto prim_pairs  = sp->prim_pairs();
    
    for(size_t p_outer = blockIdx.x * blockDim.x; p_outer < npts; p_outer += gridDim.x * blockDim.x) {
      double *_point_outer_x = (_points_x + p_outer);
      double *_point_outer_y = (_points_y + p_outer);
      double *_point_outer_z = (_points_z + p_outer);
      
      size_t p_inner = (threadIdx.x < (npts - p_outer)) ? threadIdx.x : (npts - p_outer);

      for(int i = 0; i < 1; ++i) SCALAR_STORE((temp + i * blockDim.x + threadIdx.x), SCALAR_ZERO());

      for(int ij = 0; ij < nprim_pairs; ++ij) {
	double RHO = prim_pairs[ij].gamma;

	double xA = prim_pairs[ij].P.x;
	double yA = prim_pairs[ij].P.y;
	double zA = prim_pairs[ij].P.z;
	
	double eval = prim_pairs[ij].K_coeff_prod;

	// Evaluate T Values
	SCALAR_TYPE xC = SCALAR_LOAD((_point_outer_x + p_inner));
	SCALAR_TYPE yC = SCALAR_LOAD((_point_outer_y + p_inner));
	SCALAR_TYPE zC = SCALAR_LOAD((_point_outer_z + p_inner));

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

      if(threadIdx.x < npts - p_outer) {
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

  __global__ void dev_integral_0(size_t npts,
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
    dev_integral_0_driver( npts, points_x, points_y, points_z, sp, Xi, ldX,
      Gi, ldG, weights, boys_table );
  }

  void integral_0(size_t npts,
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
    dev_integral_0<<<320, 128, 0, stream>>>(npts,
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

  __global__ void dev_integral_0_batched(
           const GauXC::ShellPairToTaskDevice* sp2task,
           GauXC::XCDeviceTask*                device_tasks,
				   double *boys_table) {

    const int ntask = sp2task->ntask;
    for( int i_task = blockIdx.y; i_task < ntask; i_task += gridDim.y ) {
    
      const auto iT = sp2task->task_idx_device[i_task];
      const auto* task  = device_tasks + iT;
      const auto  npts  = task->npts;

      const auto  i_off = sp2task->task_shell_off_row_device[i_task]*npts;

      dev_integral_0_driver( 
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


  void integral_0_batched(size_t ntask_sp,
        const GauXC::ShellPairToTaskDevice* sp2task,
        GauXC::XCDeviceTask*                device_tasks,
		    double *boys_table,
        cudaStream_t stream) {

    int nthreads = 128;
    int nblocks_x = 160;
    int nblocks_y = ntask_sp;
    dim3 nblocks(nblocks_x, nblocks_y);
    dev_integral_0_batched<<<nblocks,nthreads,0,stream>>>(
      sp2task, device_tasks, boys_table );

  }

template<bool use_shared, int max_primpairs, int points_per_subtask>
__inline__ __device__ void dev_integral_0_task(
  const int i,
  const int npts,
  const int nprim_pairs,
  // Point data
  double4 (&s_task_data)[points_per_subtask],
  // Shell Pair Data
  const shell_pair* sp,
  // Output Data
  const double *Xi,
  int ldX,
  double *Gi,
  int ldG, 
  // Other
  const double *boys_table) {

  static constexpr int num_warps = points_per_subtask / cuda::warp_size;

  const int laneId = threadIdx.x % cuda::warp_size;
  const int warpId = threadIdx.x / cuda::warp_size;

  const auto& prim_pairs = sp->prim_pairs();
  __shared__ GauXC::PrimitivePair<double> s_prim_pairs[num_warps * max_primpairs];


  if constexpr (use_shared) {
      // Load Primpairs to shared
      const int32_t* src = (int32_t*) &(prim_pairs[0]);
      int32_t* dst = (int32_t*) &(s_prim_pairs[warpId * max_primpairs]);
      const int num_transfers = nprim_pairs * sizeof(GauXC::PrimitivePair<double>) / sizeof(int32_t);

      for (int i = laneId; i < num_transfers; i += cuda::warp_size) {
        dst[i] = src[i]; 
      }
      __syncwarp();
  }

  // Loop over points in shared in batches of 32
  for (int i = 0; i <  num_warps; i++) {
    double temp = SCALAR_ZERO();

    const int pointIndex = i * cuda::warp_size + laneId;

    if (pointIndex < npts) {

      const double point_x = s_task_data[pointIndex].x;
      const double point_y = s_task_data[pointIndex].y;
      const double point_z = s_task_data[pointIndex].z;
      const double weight = s_task_data[pointIndex].w;

      for (int ij = 0; ij < nprim_pairs; ij++) {
        const GauXC::PrimitivePair<double>* prim_pairs_use = nullptr; 
        if constexpr (use_shared) prim_pairs_use = &(s_prim_pairs[warpId * max_primpairs]);
        else                      prim_pairs_use = &(prim_pairs[0]);

        double RHO = prim_pairs_use[ij].gamma;
        double xP = prim_pairs_use[ij].P.x;
        double yP = prim_pairs_use[ij].P.y;
        double zP = prim_pairs_use[ij].P.z;
        double eval = prim_pairs_use[ij].K_coeff_prod;

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
     // if (abs(temp) > 1e-12)
      {
        const double * __restrict__ Xik = (Xi + pointIndex);
        double * __restrict__ Gik = (Gi + pointIndex);

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
        t0 = SCALAR_MUL(temp, const_value_w);
        tw = SCALAR_MUL(tx, t0);
        atomicAdd(Gik, tw);
      }
    }
  }
  __syncwarp();
}

template<bool use_shared, int max_primpairs, int points_per_subtask>
__global__ void 
__launch_bounds__(points_per_subtask, 1)
dev_integral_0_task_batched(
  int ntask, int nsubtask,
  GauXC::XCDeviceTask*                device_tasks,
  const GauXC::TaskToShellPairDevice* task2sp,
  const int4* subtasks,
  const int32_t* nprim_pairs_device,
  shell_pair** sp_ptr_device,
  double *boys_table) {

  static constexpr int num_warps = points_per_subtask / cuda::warp_size;

  __shared__ double4 s_task_data[points_per_subtask];

  const int warpId = threadIdx.x / cuda::warp_size;

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

  const int npts_block = util::div_ceil(point_count, blockDim.x);

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

    for (int j = num_warps*blockIdx.y+warpId; j < nsp; j+=num_warps*gridDim.y) {
      const auto i_off = task2sp[i_task].task_shell_off_row_device[j];

      const auto index =  task2sp[i_task].shell_pair_linear_idx_device[j];
      const auto* sp = sp_ptr_device[index];
      const auto nprim_pairs = nprim_pairs_device[index];

      dev_integral_0_task<use_shared, max_primpairs, points_per_subtask>(
        i, point_count, nprim_pairs,
        s_task_data,
        sp,
        task->fmat + i_off + i,
        npts,
        task->gmat + i_off + i,
        npts,
        boys_table);
    }
    __syncthreads();
  }
}

  void integral_0_task_batched(
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
    int nthreads = cuda::obara_saika::points_per_subtask;

    int nblocks_x = nsubtask;
    int nblocks_y = 8; //std::min(max_nsp,  xy_max);
    int nblocks_z = 1;
    dim3 nblocks(nblocks_x, nblocks_y, nblocks_z);

    if (max_primpairs <= cuda::obara_saika::max_primpairs) {
      dev_integral_0_task_batched<
        cuda::obara_saika::k0_use_shared, 
        cuda::obara_saika::max_primpairs, 
        cuda::obara_saika::points_per_subtask
      ><<<nblocks, nthreads, 0, stream>>>(
        ntasks, nsubtask,
        device_tasks, task2sp, 
        (int4*) subtasks, nprim_pairs_device, sp_ptr_device,
        boys_table );
    } else {
      dev_integral_0_task_batched<
        false,
        cuda::obara_saika::max_primpairs, 
        cuda::obara_saika::points_per_subtask
      ><<<nblocks, nthreads, 0, stream>>>(
        ntasks, nsubtask,
        device_tasks, task2sp, 
        (int4*) subtasks, nprim_pairs_device, sp_ptr_device,
        boys_table );
    }
  }
}
