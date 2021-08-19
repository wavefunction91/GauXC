target_sources(gauxc PRIVATE
  cuda_aos_scheme1_data.cxx
  cuda_aos_scheme1.cxx

  kernels/collocation_device.cu
  #kernels/cuda_weights.cu
  kernels/grid_to_center.cu
  kernels/cuda_ssf_1d.cu
  cuda_aos_scheme1_weights.cu # cuda_ssf_2d w/ CudaAoSScheme1 constants
  kernels/pack_submat.cu
  kernels/cublas_extensions.cu
  kernels/uvvars.cu
  kernels/zmat_vxc.cu
  kernels/cuda_inc_potential.cu

)
