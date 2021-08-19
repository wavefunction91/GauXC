target_sources(gauxc PRIVATE
  cuda_aos_scheme1_data.cxx
  cuda_aos_scheme1.cxx
  cuda_aos_scheme1_weights.cu

  kernels/collocation_device.cu
  #kernels/cuda_weights.cu
  kernels/grid_to_center.cu
  kernels/cuda_ssf_1d.cu
)
