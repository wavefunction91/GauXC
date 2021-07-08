#pragma once


namespace GauXC      {
namespace integrator {
namespace cuda       {

template <typename T>
__inline__ __device__ void collocation_device_radial_eval(
  uint32_t       nprim,
  const T*       alpha,
  const T*       coeff,
  const T        xc,
  const T        yc,
  const T        zc,
  T*             eval_device
) {

  const T rsq = xc*xc + yc*yc + zc*zc;
  
  T tmp = 0.;
  for( uint32_t i = 0; i < nprim; ++i )
    tmp += coeff[i] * std::exp( - alpha[i] * rsq );

  *eval_device = tmp;

}



template <typename T>
__inline__ __device__ void collocation_device_radial_eval_deriv1(
  uint32_t       nprim,
  const T*       alpha,
  const T*       coeff,
  const T        xc,
  const T        yc,
  const T        zc,
  T*             eval_device,
  T*             deval_device_x,
  T*             deval_device_y,
  T*             deval_device_z
) {

  
  const T rsq = xc*xc + yc*yc + zc*zc;
  
  T tmp = 0.;
  T tmp_x = 0., tmp_y = 0., tmp_z = 0.;
  for( uint32_t i = 0; i < nprim; ++i ) {

    const T a = alpha[i];
    const T e = coeff[i] * std::exp( - a * rsq );

    const T ae = 2. * a * e;

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

} // namespace cuda
} // namespace integrator
} // namespace GauXC


