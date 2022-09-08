#pragma once

#include "cutlass/cutlass.h"
#include "cutlass/gemm/gemm.h"

#include "device/device_queue.hpp"

namespace GauXC {

void cutlass_gemm(
  cutlass::gemm::GemmCoord* problem_sizes_device,
  cutlass::gemm::GemmCoord* problem_sizes_host,
  const int problem_count,
  double ** ptr_A,
  double ** ptr_B,
  double ** ptr_C,
  double ** ptr_D,
  int64_t* lda,
  int64_t* ldb,
  int64_t* ldc,
  int64_t* ldd,
  const double alpha,
  const double beta,
  device_queue queue
);

void cutlass_syr2k(
  cutlass::gemm::GemmCoord* problem_sizes_device,
  cutlass::gemm::GemmCoord* problem_sizes_host,
  const int problem_count,
  double ** ptr_A,
  double ** ptr_B,
  double ** ptr_C,
  double ** ptr_D,
  int64_t* lda,
  int64_t* ldb,
  int64_t* ldc,
  int64_t* ldd,
  const double alpha,
  const double beta,
  device_queue queue
);

}
