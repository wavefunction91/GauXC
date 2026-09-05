#/bin/bash

if [ ! -d kernels ]
then
  mkdir kernels
fi

if [ ! -d kernels/collocation ]
then
  mkdir -p kernels/collocation
fi

export CUDA_PREFIX=$PWD/../cuda/kernels
export HIP_PREFIX=$PWD/kernels

# Generate collocation kernels
hipify-perl $CUDA_PREFIX/collocation/collocation_angular_cartesian.hpp > \
            $HIP_PREFIX/collocation/collocation_angular_cartesian.hpp
hipify-perl $CUDA_PREFIX/collocation/collocation_angular_spherical_unnorm.hpp > \
            $HIP_PREFIX/collocation/collocation_angular_spherical_unnorm.hpp
hipify-perl $CUDA_PREFIX/collocation/collocation_device_constants.hpp > \
            $HIP_PREFIX/collocation/collocation_device_constants.hpp
hipify-perl $CUDA_PREFIX/collocation_masked_combined_kernels.hpp > \
            $HIP_PREFIX/collocation_masked_combined_kernels.hpp
hipify-perl $CUDA_PREFIX/collocation_masked_kernels.hpp > \
            $HIP_PREFIX/collocation_masked_kernels.hpp
#hipify-perl $CUDA_PREFIX/collocation_device.hpp > \
#            $HIP_PREFIX/collocation_device.hpp
hipify-perl $CUDA_PREFIX/collocation_device.cu > \
            $HIP_PREFIX/collocation_device.hip


# Generate Weights Kernels
#hipify-perl $CUDA_PREFIX/grid_to_center.hpp > $HIP_PREFIX/grid_to_center.hpp
hipify-perl $CUDA_PREFIX/grid_to_center.cu  > $HIP_PREFIX/grid_to_center.hip
#hipify-perl $CUDA_PREFIX/cuda_ssf_1d.hpp > $HIP_PREFIX/hip_ssf_1d.hpp
hipify-perl $CUDA_PREFIX/cuda_ssf_1d.cu  > $HIP_PREFIX/hip_ssf_1d.hip


# cuBLAS -> hipBLAS
#hipify-perl $CUDA_PREFIX/cublas_extensions.hpp > $HIP_PREFIX/hipblas_extensions.hpp
hipify-perl $CUDA_PREFIX/cublas_extensions.cu > $HIP_PREFIX/hipblas_extensions.hip

# Z Matrix
#hipify-perl $CUDA_PREFIX/zmat_vxc.hpp > $HIP_PREFIX/zmat_vxc.hpp
hipify-perl $CUDA_PREFIX/zmat_vxc.cu > $HIP_PREFIX/zmat_vxc.hip


#hipify-perl $CUDA_PREFIX/../cuda_aos_scheme1.cxx > $HIP_PREFIX/../hip_aos_scheme1.cxx

sed -i -e "s/cuda/hip/g" kernels/{,*/}*.hpp *.{cxx,hpp}
sed -i -e "s/cuda/hip/g" kernels/*.hip
sed -i -e "s/CUDA/HIP/g" kernels/*.hip
sed -i -e "s/cublas/hipblas/g" kernels/*.hip
sed -i -e "s/CUBLAS/HIPBLAS/g" kernels/*.hip
sed -i -e "s/register //g" kernels/*.hip

#sed -i -e "s/Cuda/Hip/g" *.{cxx,hpp}
