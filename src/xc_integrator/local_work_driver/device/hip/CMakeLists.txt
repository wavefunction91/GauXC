target_sources(gauxc PRIVATE
  hip_aos_scheme1_data.cxx
  hip_aos_scheme1.cxx

  xc_functional_eval_wrapper.cxx

  kernels/collocation_device.hip
  kernels/grid_to_center.hip
  kernels/hip_ssf_1d.hip
  #cuda_aos_scheme1_weights.cu # cuda_ssf_2d w/ CudaAoSScheme1 constants
  kernels/pack_submat.hip
  kernels/hipblas_extensions.hip
  kernels/uvvars.hip
  kernels/zmat_vxc.hip
  kernels/hip_inc_potential.hip
  kernels/symmetrize_mat.hip

)

