#pragma once

namespace XGPU {
  
  void generate_shell_pair( const shells& A, const shells& B, prim_pair *prim_pairs);
  void compute_integral_shell_pair(int is_diag,
				   size_t npts,
				   double *points_x,
				   double *points_y,
				   double *points_z,
				   int lA,
				   int lB,
				   point rA,
				   point rB,
           shell_pair* sp,
				   double *Xi,
				   double *Xj,
				   int ldX,
				   double *Gi,
				   double *Gj,
				   int ldG, 
				   double *weights,
				   double *boys_table,
           cudaStream_t stream);

  void compute_integral_shell_pair_batched( int is_diag,
    size_t ntask_sp,
    int lA, int lB, 
    double X_AB,
		double Y_AB,
		double Z_AB,
    const GauXC::ShellPairToTaskDevice* sp2task,
    GauXC::XCDeviceTask*                device_tasks,
		double *boys_table,
    cudaStream_t stream); 
}
