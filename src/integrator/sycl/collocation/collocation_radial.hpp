#include <CL/sycl.hpp>
#include <iostream>
#include <cassert>

#include <gauxc/shell.hpp>


namespace GauXC      {
namespace integrator {
namespace sycl       {

__inline__  void collocation_device_radial_eval(
  const Shell&   shell,
  const double*  pt,
  double*        eval_device) {

  const auto* O     = shell.O_data();
  const auto* alpha = shell.alpha_data();
  const auto* coeff = shell.coeff_data();

  const double xc = pt[0] - O[0];
  const double yc = pt[1] - O[1];
  const double zc = pt[2] - O[2];

  const double rsq = xc*xc + yc*yc + zc*zc;

  const size_t nprim = shell.nprim;
  double tmp = 0.;
  for( size_t i = 0; i < nprim; ++i )
    tmp += coeff[i] * cl::sycl::exp( - alpha[i] * rsq );

  *eval_device = tmp;
}



__inline__ void collocation_device_radial_eval_deriv1(
  const Shell&   shell,
  const double*  pt,
  double*        eval_device,
  double*        deval_device_x,
  double*        deval_device_y,
  double*        deval_device_z) {

  const auto* O     = shell.O_data();
  const auto* alpha = shell.alpha_data();
  const auto* coeff = shell.coeff_data();

  const double xc = pt[0] - O[0];
  const double yc = pt[1] - O[1];
  const double zc = pt[2] - O[2];

  const double rsq = xc*xc + yc*yc + zc*zc;

  const size_t nprim = shell.nprim;
  double tmp = 0.;
  double tmp_x = 0., tmp_y = 0., tmp_z = 0.;
  for( size_t i = 0; i < nprim; ++i ) {

    const double a = alpha[i];
    const double e = coeff[i] * cl::sycl::exp( - a * rsq );

    const double ae = 2. * a * e;

    tmp   += e;
    tmp_x -= ae * xc;
    tmp_y -= ae * yc;
    tmp_z -= ae * zc;

  }

  *eval_device    = tmp;
  *deval_device_x = tmp_x;
  *deval_device_y = tmp_y;
  *deval_device_z = tmp_z;
}

} // namespace sycl
} // namespace integrator
} // namespace GauXC
