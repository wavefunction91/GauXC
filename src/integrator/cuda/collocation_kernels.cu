#include <iostream>
#include <cassert>

#include "collocation_kernels.hpp"
#include "collocation_angular_cartesian.hpp"
#include "collocation_angular_spherical_unnorm.hpp"

namespace GauXC      {
namespace integrator {
namespace cuda       {



template <typename T>
__global__
void collocation_device_petite_kernel(
  size_t          nshells,
  size_t          nbf,
  size_t          npts,
  const Shell<T>* shells_device,
  const size_t*   offs_device,
  const T*        pts_device,
  T*              eval_device
) {

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < npts and tid_y < nshells ) {

    const size_t ipt = tid_x;
    const size_t ish = tid_y;

    const size_t ibf = offs_device[ish];

    const auto& shell = shells_device[ish];
    const auto* pt    = pts_device + 3*ipt;
  

    const auto* O     = shell.O_data();
    const auto* alpha = shell.alpha_data();
    const auto* coeff = shell.coeff_data();

    const auto xc = pt[0] - O[0];
    const auto yc = pt[1] - O[1];
    const auto zc = pt[2] - O[2];
  
    const auto rsq = xc*xc + yc*yc + zc*zc;
  
    const size_t nprim = shell.nprim(); 
    auto tmp = 0.;
    for( size_t i = 0; i < nprim; ++i )
      tmp += coeff[i] * std::exp( - alpha[i] * rsq );

    auto * bf_eval = eval_device + ibf + ipt*nbf;

    const bool do_sph = shell.pure();
    if( do_sph )
      collocation_spherical_unnorm_angular( shell.l(), tmp, xc, yc, zc, bf_eval );
    else
      collocation_cartesian_angular( shell.l(), tmp, xc, yc, zc, bf_eval );

  }

}

template
__global__
void collocation_device_petite_kernel(
  size_t               nshells,
  size_t               nbf,
  size_t               npts,
  const Shell<double>* shells_device,
  const size_t*        offs_device,
  const double*        pts_device,
  double*              eval_device
); 














template <typename T>
__global__
void collocation_device_masked_kernel(
  size_t          nshells,
  size_t          nbf,
  size_t          npts,
  const Shell<T>* shells_device,
  const size_t*   mask_device,
  const size_t*   offs_device,
  const T*        pts_device,
  T*              eval_device
) {

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < npts and tid_y < nshells ) {

    const size_t ipt = tid_x;
    const size_t ish = tid_y;

    const size_t ibf = offs_device[ish];

    const auto& shell = shells_device[mask_device[ish]];
    const auto* pt    = pts_device + 3*ipt;
  

    const auto* O     = shell.O_data();
    const auto* alpha = shell.alpha_data();
    const auto* coeff = shell.coeff_data();

    const auto xc = pt[0] - O[0];
    const auto yc = pt[1] - O[1];
    const auto zc = pt[2] - O[2];
  
    const auto rsq = xc*xc + yc*yc + zc*zc;
  
    const size_t nprim = shell.nprim(); 
    auto tmp = 0.;
    for( size_t i = 0; i < nprim; ++i )
      tmp += coeff[i] * std::exp( - alpha[i] * rsq );

    auto * bf_eval = eval_device + ibf + ipt*nbf;

    const bool do_sph = shell.pure();
    if( do_sph )
      collocation_spherical_unnorm_angular( shell.l(), tmp, xc, yc, zc, bf_eval );
    else
      collocation_cartesian_angular( shell.l(), tmp, xc, yc, zc, bf_eval );

  }

}

template
__global__
void collocation_device_masked_kernel(
  size_t               nshells,
  size_t               nbf,
  size_t               npts,
  const Shell<double>* shells_device,
  const size_t*        mask_device,
  const size_t*        offs_device,
  const double*        pts_device,
  double*              eval_device
); 















template <typename T>
__global__
void collocation_device_petite_kernel_deriv1(
  size_t          nshells,
  size_t          nbf,
  size_t          npts,
  const Shell<T>* shells_device,
  const size_t*   offs_device,
  const T*        pts_device,
  T*              eval_device,
  T*              deval_device_x,
  T*              deval_device_y,
  T*              deval_device_z
) {

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < npts and tid_y < nshells ) {

    const size_t ipt = tid_x;
    const size_t ish = tid_y;

    const size_t ibf = offs_device[ish];

    const auto& shell = shells_device[ish];
    const auto* pt    = pts_device + 3*ipt;
  

    const auto* O     = shell.O_data();
    const auto* alpha = shell.alpha_data();
    const auto* coeff = shell.coeff_data();

    const auto xc = pt[0] - O[0];
    const auto yc = pt[1] - O[1];
    const auto zc = pt[2] - O[2];
  
    const auto rsq = xc*xc + yc*yc + zc*zc;
  
    const size_t nprim = shell.nprim(); 
    auto tmp = 0.;
    auto tmp_x = 0., tmp_y = 0., tmp_z = 0.;
    for( size_t i = 0; i < nprim; ++i ) {

      const auto a = alpha[i];
      const auto e = coeff[i] * std::exp( - a * rsq );

      const auto ae = 2. * a * e;

      tmp   += e;
      tmp_x -= ae * xc;
      tmp_y -= ae * yc;
      tmp_z -= ae * zc;

    }

    auto * bf_eval = eval_device    + ibf + ipt*nbf;
    auto * dx_eval = deval_device_x + ibf + ipt*nbf;
    auto * dy_eval = deval_device_y + ibf + ipt*nbf;
    auto * dz_eval = deval_device_z + ibf + ipt*nbf;

    const bool do_sph = shell.pure();
    if( do_sph ) 
      collocation_spherical_unnorm_angular_deriv1( shell.l(), tmp, tmp_x, tmp_y, tmp_z, 
                                               xc, yc, zc, bf_eval, dx_eval, 
                                               dy_eval, dz_eval );
    else
      collocation_cartesian_angular_deriv1( shell.l(), tmp, tmp_x, tmp_y, tmp_z, 
                                        xc, yc, zc, bf_eval, dx_eval, 
                                        dy_eval, dz_eval );

  }


}

template
__global__
void collocation_device_petite_kernel_deriv1(
  size_t               nshells,
  size_t               nbf,
  size_t               npts,
  const Shell<double>* shells_device,
  const size_t*        offs_device,
  const double*        pts_device,
  double*              eval_device,
  double*              deval_device_x,
  double*              deval_device_y,
  double*              deval_device_z
);












template <typename T>
__global__
void collocation_device_masked_kernel_deriv1(
  size_t          nshells,
  size_t          nbf,
  size_t          npts,
  const Shell<T>* shells_device,
  const size_t*   mask_device,
  const size_t*   offs_device,
  const T*        pts_device,
  T*              eval_device,
  T*              deval_device_x,
  T*              deval_device_y,
  T*              deval_device_z
) {

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < npts and tid_y < nshells ) {

    const size_t ipt = tid_x;
    const size_t ish = tid_y;

    const size_t ibf = offs_device[ish];

    const auto& shell = shells_device[mask_device[ish]];
    const auto* pt    = pts_device + 3*ipt;
  

    const auto* O     = shell.O_data();
    const auto* alpha = shell.alpha_data();
    const auto* coeff = shell.coeff_data();

    const auto xc = pt[0] - O[0];
    const auto yc = pt[1] - O[1];
    const auto zc = pt[2] - O[2];
  
    const auto rsq = xc*xc + yc*yc + zc*zc;
  
    const size_t nprim = shell.nprim(); 
    auto tmp = 0.;
    auto tmp_x = 0., tmp_y = 0., tmp_z = 0.;
    for( size_t i = 0; i < nprim; ++i ) {

      const auto a = alpha[i];
      const auto e = coeff[i] * std::exp( - a * rsq );

      const auto ae = 2. * a * e;

      tmp   += e;
      tmp_x -= ae * xc;
      tmp_y -= ae * yc;
      tmp_z -= ae * zc;

    }

    auto * bf_eval = eval_device    + ibf + ipt*nbf;
    auto * dx_eval = deval_device_x + ibf + ipt*nbf;
    auto * dy_eval = deval_device_y + ibf + ipt*nbf;
    auto * dz_eval = deval_device_z + ibf + ipt*nbf;

    const bool do_sph = shell.pure();
    if( do_sph ) 
      collocation_spherical_unnorm_angular_deriv1( shell.l(), tmp, tmp_x, tmp_y, tmp_z, 
                                               xc, yc, zc, bf_eval, dx_eval, 
                                               dy_eval, dz_eval );
    else
      collocation_cartesian_angular_deriv1( shell.l(), tmp, tmp_x, tmp_y, tmp_z, 
                                        xc, yc, zc, bf_eval, dx_eval, 
                                        dy_eval, dz_eval );

  }


}

template
__global__
void collocation_device_masked_kernel_deriv1(
  size_t               nshells,
  size_t               nbf,
  size_t               npts,
  const Shell<double>* shells_device,
  const size_t*        mask_device,
  const size_t*        offs_device,
  const double*        pts_device,
  double*              eval_device,
  double*              deval_device_x,
  double*              deval_device_y,
  double*              deval_device_z
);

} // namespace cuda
} // namespace integrator
} // namespace GauXC
