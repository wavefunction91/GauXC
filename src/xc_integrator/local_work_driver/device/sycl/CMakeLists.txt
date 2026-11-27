#
# GauXC Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of
# any required approvals from the U.S. Dept. of Energy).
#
# (c) 2024-2025, Microsoft Corporation
#
# All rights reserved.
#
# See LICENSE.txt for details
#
set( GAUXC_SYCL_LWD_SOURCES
  sycl_aos_scheme1_data.cxx
  sycl_aos_scheme1.cxx
  sycl_aos_scheme1_weights.cxx

  xc_functional_eval_wrapper.cxx

  kernels/collocation_device.cxx
  kernels/grid_to_center.cxx
  kernels/sycl_ssf_1d.cxx
  kernels/sycl_ssf_2d.cxx
  kernels/pack_submat.cxx
  kernels/onemkl_extensions.cxx
  kernels/uvvars.cxx
  kernels/zmat_vxc.cxx
  kernels/zmat_fxc.cxx
  kernels/sycl_inc_potential.cxx
  kernels/symmetrize_mat.cxx
  kernels/increment_exc_grad.cxx
)

target_sources(gauxc PRIVATE ${GAUXC_SYCL_LWD_SOURCES})
gauxc_sycl_device_sources( ${GAUXC_SYCL_LWD_SOURCES} )
