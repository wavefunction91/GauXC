#pragma once
#include <gauxc/gauxc_config.hpp>

namespace GauXC {

// FWD decl all exception types for optional handling

#ifdef GAUXC_ENABLE_CUDA
class cuda_exception;
class cublas_exception;
#endif

#ifdef GAUXC_ENABLE_HIP
class hip_exception;
class hipblas_exception;
#endif

#ifdef GAUXC_ENABLE_MAGMA
class magma_exception;
#endif

}
