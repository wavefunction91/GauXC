#include <iostream>
#include <algorithm>

#include "cutlass/cutlass.h"
#include "cutlass/gemm/gemm.h"
#include "cutlass/gemm/kernel/gemm_grouped.h"
#include "cutlass/gemm/kernel/default_gemm_grouped.h"
#include "cutlass/gemm/device/gemm_grouped.h"
#include "cutlass/gemm/device/gemm_universal.h"

#include "cutlass/util/device_memory.h"

#include "device_specific/cuda_device_constants.hpp"
#include "device_specific/cuda_util.hpp"
#include "device/device_queue.hpp"

namespace GauXC {

__global__ void symmetrize_batch_matrix_device( 
  size_t num_matrices, 
  double** As, 
  int64_t *LDAs, 
  cutlass::gemm::GemmCoord* problem_sizes_device
) {
  constexpr uint32_t block_size = cuda::warp_size;

  __shared__ double buffer[block_size][block_size+1];  // Pad shared memory to resolve shared memory

  for (int n = blockIdx.x; n < num_matrices; n += gridDim.x) {
    const int64_t LDA  = LDAs[n];
    const int32_t N     = problem_sizes_device[n].n();
    double*       A    = As[n];

    const size_t num_blocks = ((N + block_size - 1) / block_size);

    for (int i = 0; i < num_blocks; i++) {
      const int i_coord = i * block_size;
      for (int j = i; j < num_blocks; j++) {
        const int j_coord = j * block_size;

        if (i_coord + threadIdx.y < N && j_coord + threadIdx.x < N) {
          buffer[threadIdx.y][threadIdx.x] = A[(i_coord + threadIdx.y) * LDA + j_coord + threadIdx.x];
        }
        __syncthreads();

        if (j_coord + threadIdx.y < N && i_coord + threadIdx.x < N) {
          buffer[threadIdx.x][threadIdx.y] += A[(j_coord + threadIdx.y) * LDA + i_coord + threadIdx.x];
        }
        __syncthreads();

        if (j_coord + threadIdx.y < N && i_coord + threadIdx.x < N) {
          A[(j_coord + threadIdx.y) * LDA + i_coord + threadIdx.x] = buffer[threadIdx.x][threadIdx.y];
        }
        if (i_coord + threadIdx.y < N && j_coord + threadIdx.x < N) {
          A[(i_coord + threadIdx.y) * LDA + j_coord + threadIdx.x] = buffer[threadIdx.y][threadIdx.x];
        }
      }
      __syncthreads();
    }
  }
}


template<typename CutlassGemm>
void cutlass_wrapper_launch(
  cutlass::gemm::GemmCoord* problem_sizes_device,
  const int problem_count,
  const double alpha,
  const double beta,
  const int threadblock_count,
  double ** ptr_A,
  double ** ptr_B,
  double ** ptr_C,
  double ** ptr_D,
  int64_t* lda,
  int64_t* ldb,
  int64_t* ldc,
  int64_t* ldd,
  cudaStream_t stream
) {
  cutlass::Status status;

  typename CutlassGemm::GemmKernel::Epilogue::OutputOp::Params epilogue_op(alpha, beta);

  // Configure GEMM arguments
  typename CutlassGemm::Arguments args(
    problem_sizes_device,
    problem_count,
    threadblock_count,
    epilogue_op,
    ptr_A,
    ptr_B,
    ptr_C,
    ptr_D,
    lda,
    ldb,
    ldc,
    ldd
  );

  // Initialize the GEMM object
  CutlassGemm gemm;

  status = gemm.initialize(args);

  if (status != cutlass::Status::kSuccess) {
    std::cerr << "Failed to initialize CUTLASS Grouped GEMM kernel." << std::endl;
    return;
  }

  // Run the grouped GEMM object
  status = gemm.run(stream);

  if (status != cutlass::Status::kSuccess) {
    std::cerr << "Failed to run CUTLASS Grouped GEMM kernel." << std::endl;
    return;
  }
}


template<typename CutlassGemm>
int getThreadBlockCount() {
  int smem_size = int(sizeof(typename CutlassGemm::GemmKernel::SharedStorage));
  cudaDeviceProp properties;

  cudaGetDeviceProperties(&properties, 0);
  int occupancy = std::min(2, int(properties.sharedMemPerMultiprocessor / smem_size));

  return properties.multiProcessorCount * occupancy;
}


void cutlass_gemm(
  cutlass::gemm::GemmCoord* problem_sizes_device,
  const int problem_count,
  double ** ptr_A,
  double ** ptr_B,
  double ** ptr_C,
  double ** ptr_D,
  int64_t* lda,
  int64_t* ldb,
  int64_t* ldc,
  int64_t* ldd,
  device_queue queue
) {
  using ElementOutput = double;
  using ElementAccumulator = double;
  using ElementA = double; 
  using ElementB = double; 
  using ElementC = double; 

  static int const kAlignmentA = 1;
  static int const kAlignmentB = 1;
  
  using ThreadblockShape = cutlass::gemm::GemmShape<64, 64, 16>;
  using WarpShape = cutlass::gemm::GemmShape<32, 32, 16>;
  using InstructionShape = cutlass::gemm::GemmShape<8, 8, 4>;
  static int const kStages = 4;

  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    ElementC,  1,
    ElementAccumulator,
    ElementAccumulator>; 

  using GemmKernel = typename cutlass::gemm::kernel::DefaultGemmGrouped<
    ElementA, 
    cutlass::layout::ColumnMajor, 
    cutlass::ComplexTransform::kNone,
    kAlignmentA,
    ElementB,
    cutlass::layout::ColumnMajor, 
    cutlass::ComplexTransform::kNone,
    kAlignmentB,
    ElementOutput, cutlass::layout::ColumnMajor,
    ElementAccumulator, 
    cutlass::arch::OpClassTensorOp, 
    cutlass::arch::Sm80,
    ThreadblockShape,
    WarpShape,
    InstructionShape,
    EpilogueOutputOp,
    cutlass::gemm::threadblock::GemmBatchedIdentityThreadblockSwizzle, 
    kStages>::GemmKernel;

  using GemmGrouped = cutlass::gemm::device::GemmGrouped<GemmKernel>; 

  const int threadblock_count = getThreadBlockCount<GemmGrouped>();

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();

  const double alpha = 1.0;
  const double beta  = 0.0;

  cutlass_wrapper_launch<GemmGrouped>(
    problem_sizes_device,
    problem_count,
    alpha, beta,
    threadblock_count,
    ptr_A,
    ptr_B,
    ptr_C,
    ptr_D,
    lda,
    ldb,
    ldc,
    ldd,
    stream
  );
}


void cutlass_syr2k(
  cutlass::gemm::GemmCoord* problem_sizes_device,
  const int problem_count,
  double ** ptr_A,
  double ** ptr_B,
  double ** ptr_C,
  double ** ptr_D,
  int64_t* lda,
  int64_t* ldb,
  int64_t* ldc,
  int64_t* ldd,
  device_queue queue
) {
  using ElementOutput = double;
  using ElementAccumulator = double;
  using ElementA = double; 
  using ElementB = double; 
  using ElementC = double; 

  static int const kAlignmentA = 1;
  static int const kAlignmentB = 1;
  
  using ThreadblockShape = cutlass::gemm::GemmShape<64, 64, 16>;
  using WarpShape = cutlass::gemm::GemmShape<32, 32, 16>;
  using InstructionShape = cutlass::gemm::GemmShape<8, 8, 4>;
  static int const kStages = 4;

  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    ElementC,  1,
    ElementAccumulator,
    ElementAccumulator>; 

  using GemmKernel = typename cutlass::gemm::kernel::DefaultGemmGrouped<
    ElementA, 
    cutlass::layout::RowMajor, 
    cutlass::ComplexTransform::kNone,
    kAlignmentA,
    ElementB,
    cutlass::layout::ColumnMajor, 
    cutlass::ComplexTransform::kNone,
    kAlignmentB,
    ElementOutput, cutlass::layout::ColumnMajor,
    ElementAccumulator, 
    cutlass::arch::OpClassTensorOp, 
    cutlass::arch::Sm80,
    ThreadblockShape,
    WarpShape,
    InstructionShape,
    EpilogueOutputOp,
    cutlass::gemm::threadblock::GemmBatchedIdentityThreadblockSwizzle, 
    kStages>::GemmKernel;

  using GemmGrouped = cutlass::gemm::device::GemmGrouped<GemmKernel>; 

  const int threadblock_count = getThreadBlockCount<GemmGrouped>();

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();

  cutlass_wrapper_launch<GemmGrouped>(
    problem_sizes_device,
    problem_count,
    1.0, 0.0,
    threadblock_count,
    ptr_A,
    ptr_B,
    ptr_C,
    ptr_D,
    lda,
    ldb,
    ldc,
    ldd,
    stream
  );

  const size_t num_blocks = ((problem_count + cuda::warp_size - 1) / cuda::warp_size);
  // Warp size must equal max_warps_per_thread_block must equal 32
  dim3 threads(cuda::warp_size, cuda::max_warps_per_thread_block), blocks(num_blocks);
  symmetrize_batch_matrix_device<<<blocks, threads, 0, stream>>>(problem_count, ptr_D, ldd, problem_sizes_device);
}

}
