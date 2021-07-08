#!/bin/sh
hipify-perl ../cuda/collocation/collocation_angular_cartesian.hpp        > collocation/collocation_angular_cartesian.hpp
hipify-perl ../cuda/collocation/collocation_angular_spherical_unnorm.hpp > collocation/collocation_angular_spherical_unnorm.hpp
hipify-perl ../cuda/collocation/collocation_device_constants.hpp         > collocation/collocation_device_constants.hpp
hipify-perl ../cuda/collocation/collocation_radial.hpp                   > collocation/collocation_radial.hpp

#hipify-perl ../cuda/collocation_device.cu                    > collocation_device.cxx
hipify-perl ../cuda/collocation_device.hpp                   > collocation_device.hpp
hipify-perl ../cuda/collocation_masked_combined_kernels.hpp  > collocation_masked_combined_kernels.hpp
hipify-perl ../cuda/collocation_masked_kernels.hpp           > collocation_masked_kernels.hpp
hipify-perl ../cuda/collocation_petite_combined_kernels.hpp  > collocation_petite_combined_kernels.hpp
hipify-perl ../cuda/collocation_petite_kernels.hpp           > collocation_petite_kernels.hpp
hipify-perl ../cuda/cublas_extensions.cu                     > hipblas_extensions.cxx
hipify-perl ../cuda/cublas_extensions.hpp                    > hipblas_extensions.hpp
#hipify-perl ../cuda/cuda_eval_denvars.cu                     > hip_eval_denvars.cxx
hipify-perl ../cuda/cuda_eval_denvars.hpp                    > hip_eval_denvars.hpp
hipify-perl ../cuda/cuda_extensions.hpp                      > hip_extensions.hpp
hipify-perl ../cuda/cuda_alg_variant_control.hpp             > hip_alg_variant_control.hpp
hipify-perl ../cuda/cuda_inc_potential.cu                    > hip_inc_potential.cxx
hipify-perl ../cuda/cuda_inc_potential.hpp                   > hip_inc_potential.hpp
hipify-perl ../cuda/cuda_pack_density.cu                     > hip_pack_density.cxx
hipify-perl ../cuda/cuda_pack_density.hpp                    > hip_pack_density.hpp
hipify-perl ../cuda/cuda_weights.cu                          > hip_weights.cxx
hipify-perl ../cuda/cuda_weights.hpp                         > hip_weights.hpp
hipify-perl ../cuda/cuda_zmat.cu                             > hip_zmat.cxx
hipify-perl ../cuda/cuda_zmat.hpp                            > hip_zmat.hpp
hipify-perl ../cuda/xc_cuda_data.cxx                         > xc_hip_data.cxx
hipify-perl ../cuda/xc_cuda_util.cxx                         > xc_hip_util.cxx

sed -i -e "s/cuda/hip/g" *.cxx *.hpp collocation/*.hpp
sed -i -e "s/CUDA/HIP/g" *.cxx *.hpp collocation/*.hpp
sed -i -e "s/Cuda/Hip/g" *.cxx *.hpp collocation/*.hpp
sed -i -e "s/cublas/hipblas/g" *.cxx *.hpp collocation/*.hpp
sed -i -e "s/CUBLAS/HIPBLAS/g" *.cxx *.hpp collocation/*.hpp

sed -i -e "s/__global__/__global__\n__launch_bounds__(1024,1)/g" \
	collocation_masked_combined_kernels.hpp collocation_masked_kernels.hpp collocation_petite_combined_kernels.hpp collocation_petite_kernels.hpp

sed -i -e "s/register //g" *.cxx *.hpp
sed -i -e "s/#define GAUXC_HIP_ENABLE_WARP_REDUCTIONS/\/\/#define GAUXC_HIP_ENABLE_WARP_REDUCTIONS/g" hip_alg_variant_control.hpp
