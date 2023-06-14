#include <iostream>
#include <cassert>

#include <gauxc/shell.hpp>
#include <gauxc/xc_task.hpp>

#include "cuda/collocation/collocation_angular_cartesian.hpp"
#include "cuda/collocation/collocation_angular_spherical_unnorm.hpp"
#include "cuda/cuda_alg_variant_control.hpp"

namespace GauXC      {
namespace integrator {
namespace cuda       {

using namespace GauXC::cuda;

template <typename T>
__global__
void collocation_device_masked_combined_kernel(
  size_t                        ntasks,
  Shell<T>*        __restrict__ shells_device,
  XCTaskDevice<T>* __restrict__ device_tasks
) {

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( blockIdx.z < ntasks ) {

    auto& task = device_tasks[ blockIdx.z ];
  
    const auto               nshells     = task.nshells;
    const auto               nbf         = task.nbe;
    const auto               npts        = task.npts;
    const auto* __restrict__ pts_device  = task.points;
    const auto* __restrict__ mask_device = task.shell_list;
    const auto* __restrict__ offs_device = task.shell_offs;

    auto* __restrict__ eval_device    = task.bf;

  if( tid_x < npts and tid_y < nshells ) {

    const uint32_t ipt = tid_x;
    const uint32_t ish = tid_y;
    const uint32_t ibf = offs_device[ish];

    const auto& shell = shells_device[mask_device[ish]];
    const auto* pt    = pts_device + 3*ipt;
  

    const auto* __restrict__ O     = shell.O_data();
    const auto* __restrict__ alpha = shell.alpha_data();
    const auto* __restrict__ coeff = shell.coeff_data();

    const auto xc = pt[0] - O[0];
    const auto yc = pt[1] - O[1];
    const auto zc = pt[2] - O[2];
  
    const auto rsq = xc*xc + yc*yc + zc*zc;
  
    const uint32_t nprim = shell.nprim(); 
    auto tmp = 0.;
    for( uint32_t i = 0; i < nprim; ++i )
      tmp += coeff[i] * std::exp( - alpha[i] * rsq );

    auto * bf_eval = eval_device + ibf*npts + ipt;

    const bool do_sph = shell.pure();
    if( do_sph )
      collocation_spherical_unnorm_angular( npts, shell.l(), tmp, xc, yc, zc, bf_eval );
    else
      collocation_cartesian_angular( npts, shell.l(), tmp, xc, yc, zc, bf_eval );

  } // shell / point idx check

  } // Batch idx check

}















template <typename T>
__global__
void collocation_device_masked_combined_kernel_deriv1(
  size_t                        ntasks,
  Shell<T>*        __restrict__ shells_device,
  XCTaskDevice<T>* __restrict__ device_tasks
) {

  // DBWY: These are factored into the loop for this optimization
  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( blockIdx.z < ntasks ) {

    auto& task = device_tasks[ blockIdx.z ];
  
    const auto               nshells     = task.nshells;
    const auto               nbf         = task.nbe;
    const auto               npts        = task.npts;
    const auto* __restrict__ pts_device  = task.points;
    const auto* __restrict__ mask_device = task.shell_list;
    const auto* __restrict__ offs_device = task.shell_offs;

    auto* __restrict__ eval_device    = task.bf;
    auto* __restrict__ deval_device_x = task.dbfx;
    auto* __restrict__ deval_device_y = task.dbfy;
    auto* __restrict__ deval_device_z = task.dbfz;

  if( tid_y < nshells and tid_x < npts ) {

    const uint32_t ish = tid_y;
    const uint32_t ipt = tid_x;
    const uint32_t ibf = offs_device[ish];

    const auto& shell = shells_device[mask_device[ish]];

    const auto* pt    = pts_device + 3*ipt;
  

    const auto* __restrict__ O     = shell.O_data();
    const auto* __restrict__ alpha = shell.alpha_data();
    const auto* __restrict__ coeff = shell.coeff_data();

    const auto xc = pt[0] - O[0];
    const auto yc = pt[1] - O[1];
    const auto zc = pt[2] - O[2];
  
    const auto rsq = xc*xc + yc*yc + zc*zc;
  
    const uint32_t nprim = shell.nprim(); 
    auto tmp = 0.;
    auto tmp_x = 0., tmp_y = 0., tmp_z = 0.;
    for( uint32_t i = 0; i < nprim; ++i ) {

      const auto a = alpha[i];
      const auto e = coeff[i] * std::exp( - a * rsq );

      const auto ae = 2. * a * e;

      tmp   += e;
      tmp_x -= ae * xc;
      tmp_y -= ae * yc;
      tmp_z -= ae * zc;

    }

    auto * bf_eval = eval_device    + ibf*npts + ipt;
    auto * dx_eval = deval_device_x + ibf*npts + ipt;
    auto * dy_eval = deval_device_y + ibf*npts + ipt;
    auto * dz_eval = deval_device_z + ibf*npts + ipt;

    const bool do_sph = shell.pure();
    if( do_sph ) 
      collocation_spherical_unnorm_angular_deriv1( npts, shell.l(), tmp, tmp_x, tmp_y, 
                                               tmp_z, xc, yc, zc, bf_eval, dx_eval, 
                                               dy_eval, dz_eval );
    else
      collocation_cartesian_angular_deriv1( npts, shell.l(), tmp, tmp_x, tmp_y, tmp_z, 
                                        xc, yc, zc, bf_eval, dx_eval, 
                                        dy_eval, dz_eval );

  } // shell / point idx check
  } // Batch idx check


}

} // namespace cuda
} // namespace integrator
} // namespace GauXC
