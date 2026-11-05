#!/bin/bash
#
# Regenerate the HIP kernels from their CUDA counterparts.
#
# The HIP backend is a machine translation of the CUDA backend.  This
# script MUST be re-run whenever a CUDA kernel is added or changed,
# otherwise the shared (backend-agnostic) driver in ../scheme1_base.cxx
# calls into kernels that have no HIP definition and the HIP build fails
# to link.  That is exactly what happened between the 2022 port and the
# 1.0 release: mGGA, FXC contraction, EXC gradients, sn-LinK screening
# and the shell-to-task collocation were added on the CUDA side only.
#
# The translation is performed with explicit sed rules rather than
# hipify-perl so that it is reproducible on machines without a ROCm
# installation, and so that the handful of NON-mechanical decisions are
# documented in one place:
#
#   * Where possible, warp/wavefront size is NOT hardcoded; kernels take
#     it from GauXC::cuda::warp_size -> GauXC::hip::warp_size (32 vs 64),
#     so the launch geometry adapts. A handful of kernels below instead
#     hardcode a fixed 32 for HIP -- these are cases where 32 is really a
#     tile/launch WIDTH that must stay in lockstep with an array extent
#     sized against CUDA's warp_size, not a stand-in for the true
#     wavefront width; see symmetrize_mat.cu, uvvars_gga/mgga.hpp and
#     exx_ek_screening_bfn_stats.cu below.  Kernels that assume a 32-wide
#     reduction must otherwise be reviewed by hand -- see CHECK_WAVEFRONT
#     below.
#   * __syncwarp() has no HIP equivalent; wavefronts execute in lockstep
#     on AMD, so it is commented out (the convention already used by the
#     2022 port in this directory).
#   * __shfl_*_sync(mask, ...) -> __shfl_*(...): HIP shuffles take no
#     mask argument.
#   * cub -> hipcub.
#   * nvcc implicitly converts a __global__ function pointer to
#     const void* (used to pass a kernel to cudaFuncSetAttribute);
#     clang/HIP does not, so the call needs an explicit
#     reinterpret_cast<const void*>(...).
#   * symmetrize_mat.cu's block-transpose kernels use a SQUARE
#     block_size x block_size shared tile and launch a matching square
#     block_size x block_size thread block, hardcoding
#     block_size == warp_size == max_warps_per_thread_block (true on
#     CUDA: 32 == 32 == 32). On HIP that triple equality is impossible
#     (warp_size 64 forces max_warps_per_thread_block down to 16, and
#     64x64 would exceed the 1024 threads/block hardware limit anyway),
#     so block_size is decoupled to a fixed 32 for HIP -- independent of
#     hip::warp_size, still square, still <=1024 threads/block and well
#     under the LDS budget.
#   * uvvars_gga.hpp/uvvars_mgga.hpp's eval_vvar_{gga,mgga}_kern build a
#     __shared__ den_shared[rows][warp_size][SM_BLOCK+1] tile that is
#     WRITTEN with threadIdx.x indexing the warp_size dimension and READ
#     BACK with threadIdx.x indexing the SM_BLOCK dimension (a transpose).
#     That is only a valid transpose -- and only avoids reading
#     uninitialised/aliased entries -- when blockDim.x == SM_BLOCK ==
#     warp_size, which is how CUDA sizes it (32 == 32 == 32) and how it is
#     launched (uvvars.cu: dim3 threads(cuda::warp_size, ...)). This local
#     "warp_size" is therefore the tile/launch width, not the true
#     warp/wavefront width (contrast the direct-atomicAdd reduction in
#     uvvars_lda.hpp, which has no tile and is correct at the true 64-wide
#     hip::warp_size). Naively carrying hip::warp_size (64) into it, as the
#     generic cuda::->hip:: rule does, keeps SM_BLOCK/the launch at 32 and
#     breaks the invariant: reads alias into neighbouring point rows, and
#     the reduction sums 64 lanes of which only 32 were ever written --
#     silently wrong densities, not a crash. So this local warp_size, and
#     the matching GGA/MGGA vvars launch geometry in uvvars.hip, are
#     pinned to a fixed 32 for HIP, independent of hip::warp_size --
#     matching CUDA's numbers exactly and keeping the tile
#     ([rows][32][33]) well under the 64KB LDS limit without needing to
#     shrink SM_BLOCK at all.
#   * <cuda.h> (the CUDA driver API header) is dropped rather than mapped
#     to a HIP equivalent; nothing in the translated files uses driver-API
#     symbols, and hip/hip_runtime.h (inserted below) already covers the
#     runtime-API surface these files actually rely on.
#   * CUTLASS has no HIP counterpart; cutlass_wrapper is CUDA-only and is
#     deliberately not translated (it is guarded by GAUXC_ENABLE_CUTLASS,
#     which is a CUDA-only dependent option in the top-level CMakeLists).
#   * `register` is accepted (as a no-op) by nvcc but is ISO-illegal in
#     C++17, which the HIP/clang compiler enforces; the keyword is
#     stripped rather than translated.
#   * The shell-to-task collocation kernels load one primitive
#     coefficient/exponent per lane, in lockstep, on the assumption that
#     shell_nprim_max (32, a Shell data-layout constant, see shell.hpp)
#     equals the warp width -- true for cuda::warp_size (32) but not for
#     hip::warp_size (64, an AMD wavefront). The static_assert is
#     loosened to <= and the two loads are guarded so lanes >=
#     shell_nprim_max don't read/write out of bounds (those lanes are NOT
#     idle on AMD -- they go on to process points in the ipt = threadIdx.x
#     % hip::warp_size loop below; only the primitive load is guarded).
#     That guard is necessary but not sufficient: my_alpha/my_coeff pick
#     their per-wavefront __shared__ row via threadIdx.x/32, an
#     untranslated literal left over from CUDA (where 32 coincides with
#     cuda::warp_size). On HIP the true wavefront is 64 lanes, so each
#     64-lane wavefront spans TWO threadIdx.x/32 groups/rows; the guard
#     above (warp_rank = threadIdx.x % hip::warp_size < shell_nprim_max)
#     is true only for the low 32 lanes of the wavefront, so the row
#     belonging to the high 32 lanes is never written and is read back as
#     uninitialised garbage. threadIdx.x/32 is changed to
#     threadIdx.x/hip::warp_size so the row groups line up with actual
#     wavefronts (and, correspondingly, with warp_rank's own use of
#     hip::warp_size); the alpha/coeff shared arrays have 16 rows, so at
#     64 lanes/row this uses at most 8 of them -- no resize needed.
#   * pack_submat.cu/cuda_inc_potential.cu's submatrix scatter/gather
#     kernels (pack density into per-task submatrices; gather task-local
#     VXC/K contributions back into the global matrix) iterate over
#     CUT_Y/CUT_X-sized blocks and WARP_X-wide lanes -- WARP_X (16),
#     CUT_X, CUT_Y (8, 8) are plain tuning constants, untouched by any
#     cuda::/hip:: translation, on the assumption that the launch
#     geometry equals them exactly: cuda::warp_size/2 == 16 == WARP_X and
#     cuda::max_warps_per_thread_block*2 == 64 == CUT_Y*CUT_X. On HIP the
#     generic translation instead gives (32, 32): tid_yy (=threadIdx.y/
#     CUT_X) only reaches half of CUT_Y's range, so every other block of
#     rows is never visited (undercounting -- explains e.g. EXC/N_EL
#     coming out systematically low), while tid_xy (=threadIdx.x/WARP_X,
#     always 0 on CUDA since blockDim.x==WARP_X there) now also takes 1,
#     redundantly re-covering columns already handled by tid_xy==0 and
#     double-atomicAdd-ing them (undefined accumulation order -- explains
#     the run-to-run non-determinism in repeated identical evaluations).
#     Neither kernel body calls a warp-synchronous primitive (only
#     atomicAdd), so there is no real dependence on the physical
#     wavefront width here at all; the launch geometry is pinned to a
#     fixed (16, 64), matching CUDA's numbers exactly.
#
# Usage:  ./hipify.sh          (from this directory)

set -euo pipefail

CUDA_PREFIX=$PWD/../cuda/kernels
HIP_PREFIX=$PWD/kernels

mkdir -p "$HIP_PREFIX/collocation"

hipify_file() {
  local src=$1 dst=$2
  sed \
    -e 's|device_specific/cuda_util\.hpp|device_specific/hip_util.hpp|g' \
    -e 's|device_specific/cuda_device_constants\.hpp|device_specific/hip_device_constants.hpp|g' \
    -e 's|device_specific/cublas_util\.hpp|device_specific/hipblas_util.hpp|g' \
    -e 's|cuda_extensions\.hpp|hip_extensions.hpp|g' \
    -e 's|cuda_aos_scheme1\.hpp|hip_aos_scheme1.hpp|g' \
    -e 's|cuda_ssf_1d\.hpp|hip_ssf_1d.hpp|g' \
    -e 's|\bregister[[:space:]]\+||g' \
    -e 's|#include <cub/\(.*\)\.cuh>|#include <hipcub/\1.hpp>|g' \
    -e '/^#include <cuda\.h>$/d' \
    -e 's|\bcub::|hipcub::|g' \
    -e 's|\bcuda::|hip::|g' \
    -e 's|\bnamespace cuda\b|namespace hip|g' \
    -e 's|\bcudaStream_t\b|hipStream_t|g' \
    -e 's|\bcuda_stream\b|hip_stream|g' \
    -e 's|\bcudaError_t\b|hipError_t|g' \
    -e 's|\bcudaSuccess\b|hipSuccess|g' \
    -e 's|\bcudaGetErrorString\b|hipGetErrorString|g' \
    -e 's|\bcudaGetLastError\b|hipGetLastError|g' \
    -e 's|\bcudaDeviceSynchronize\b|hipDeviceSynchronize|g' \
    -e 's|\bcudaMalloc\b|hipMalloc|g' \
    -e 's|\bcudaFree\b|hipFree|g' \
    -e 's|\bcudaMemcpy|hipMemcpy|g' \
    -e 's|\bcudaMemset|hipMemset|g' \
    -e 's|device/cuda/kernels|device/hip/kernels|g' \
    -e 's|util::cuda_|util::hip_|g' \
    -e 's|\bcuda_kernel_max_threads_per_block\b|hip_kernel_max_threads_per_block|g' \
    -e 's|\bcudaDeviceGetAttribute\b|hipDeviceGetAttribute|g' \
    -e 's|\bcudaDevAttrMaxSharedMemoryPerBlockOptin\b|hipDeviceAttributeMaxSharedMemoryPerBlock|g' \
    -e 's|\bcudaDevAttrMaxSharedMemoryPerBlock\b|hipDeviceAttributeMaxSharedMemoryPerBlock|g' \
    -e 's|\bcudaFuncSetAttribute\b|hipFuncSetAttribute|g' \
    -e 's|hipFuncSetAttribute(&\([A-Za-z_][A-Za-z0-9_]*\),|hipFuncSetAttribute(reinterpret_cast<const void*>(\&\1),|g' \
    -e 's|\bcudaFuncAttribute|hipFuncAttribute|g' \
    -e 's|\bcuda_exception\b|hip_exception|g' \
    -e 's|GAUXC_CUDA|GAUXC_HIP|g' \
    -e 's|GAUXC_CUBLAS|GAUXC_HIPBLAS|g' \
    -e 's|cuda_exception\.hpp|hip_exception.hpp|g' \
    -e 's|\bCUDA_|HIP_|g' \
    -e 's|\bCUBLAS_|HIPBLAS_|g' \
    -e 's|\bcublas|hipblas|g' \
    -e 's|__shfl_\([a-z]*\)_sync *( *[^,]*, *|__shfl_\1(|g' \
    -e 's|\(^[[:space:]]*\)__syncwarp();|\1// __syncwarp();  // lockstep wavefronts on AMD|g' \
    -e 's|static_assert( detail::shell_nprim_max == hip::warp_size );|static_assert( detail::shell_nprim_max <= hip::warp_size );|g' \
    -e 's|my_alpha\[warp_rank\] = alpha_gm\[warp_rank\];|if( warp_rank < detail::shell_nprim_max ) my_alpha[warp_rank] = alpha_gm[warp_rank];|g' \
    -e 's|my_coeff\[warp_rank\] = coeff_gm\[warp_rank\];|if( warp_rank < detail::shell_nprim_max ) my_coeff[warp_rank] = coeff_gm[warp_rank];|g' \
    -e 's|alpha\[threadIdx\.x/32\]|alpha[threadIdx.x/hip::warp_size]|g' \
    -e 's|coeff\[threadIdx\.x/32\]|coeff[threadIdx.x/hip::warp_size]|g' \
    -e 's|constexpr uint32_t block_size = hip::warp_size;|constexpr uint32_t block_size = 32; // fixed, independent of hip::warp_size -- see hipify.sh|g' \
    -e 's|// Warp size must equal max_warps_per_thread_block must equal 32|// block_size is fixed at 32 for HIP, independent of hip::warp_size -- see hipify.sh|g' \
    -e 's|dim3 threads(hip::warp_size, hip::max_warps_per_thread_block), blocks(num_blocks);|dim3 threads(32, 32), blocks(num_blocks);|g' \
    -e 's|const size_t num_blocks = ((N + hip::warp_size - 1) / hip::warp_size);|const size_t num_blocks = ((N + 32 - 1) / 32);|g' \
    -e 's|dim3 threads( hip::warp_size/2, hip::max_warps_per_thread_block \* 2, 1 );|dim3 threads( 16, 64, 1 ); // matches WARP_X=16, CUT_Y*CUT_X=64 below, independent of hip::warp_size -- see hipify.sh|g' \
    "$src" > "$dst"

  # CUDA cache-hint stores have no HIP counterpart: __stcs() is an NVIDIA
  # intrinsic and the pre-CUDA-11 fallback is inline PTX.  CUDART_VERSION
  # is undefined under HIP, so the preprocessor would otherwise select the
  # PTX branch and fail to compile.  Keep the intrinsic branch and lower
  # the store to a plain one (the hint is an optimization, not semantics).
  if grep -q 'CUDART_VERSION' "$dst"; then
    awk '
      /^#if \(CUDART_VERSION/ { skipelse = 1; next }
      /^#else/ && skipelse    { drop = 1; next }
      /^#endif/ && skipelse   { skipelse = 0; drop = 0; next }
      drop                    { next }
      { print }
    ' "$dst" > "$dst.tmp"
    sed -e 's|__stcs( *\([^,]*\), *\([^)]*\));|*(\1) = \2;|g' \
        "$dst.tmp" > "$dst"
    rm -f "$dst.tmp"
  fi

  # bitvector_to_position_list_{shellpair,shells} in
  # exx_ek_screening_bfn_stats.cu declare a local "warp_size" that is
  # really just the (hardcoded, architecture-independent) 32x32 launch
  # geometry from dim3 threads(32,32) at the call site -- the name is a
  # coincidence of CUDA's warp being 32 wide, not a dependence on the
  # true warp/wavefront width (unlike e.g. the warp_reduce_* callers in
  # uvvars_*.hpp, which this pattern does not match). Sizing its
  # collisions_buffer[warp_size][warp_size][...] off hip::warp_size (64)
  # quadruples the tile and overflows AMD's 64KB LDS; pin it back to the
  # actual launch geometry instead.
  if grep -q 'collisions_buffer\[warp_size\]\[warp_size\]\[buffer_size\]' "$dst"; then
    sed -i -e 's|constexpr auto warp_size = hip::warp_size;|constexpr auto warp_size = 32; // matches dim3 threads(32,32) below, independent of hip::warp_size -- see hipify.sh|g' \
        "$dst"
  fi

  # eval_vvar_{gga,mgga}_kern in uvvars_gga.hpp/uvvars_mgga.hpp (marked by
  # the den_shared[] tile they declare) use "warp_size" as a tile/launch
  # WIDTH that must stay 32 in lockstep with SM_BLOCK and the kernel's
  # launch geometry, not the true wavefront width -- see the header
  # comment above. Pin it, matching the exx_ek_screening treatment above.
  if grep -q 'den_shared\[' "$dst"; then
    sed -i -e 's|constexpr auto warp_size = hip::warp_size;|constexpr auto warp_size = 32; // tile/launch width, independent of hip::warp_size -- see hipify.sh|g' \
        "$dst"
  fi

  # eval_vvars_gga_impl/eval_vvars_mgga_impl in uvvars.hip launch
  # eval_vvar_{gga,mgga}_kern with dim3 threads( hip::warp_size, ... ) --
  # the generic cuda::->hip:: translation of CUDA's dim3
  # threads(cuda::warp_size, ...). blockDim.x must match the tile/launch
  # width pinned to 32 above (not the true 64-wide hip::warp_size), so
  # this launch geometry is fixed to 32 for those two impls only;
  # eval_vvars_lda_impl's identical-looking launch has no tile and is
  # correctly left at hip::warp_size (uvvars_lda.hpp's reduction is a
  # direct atomicAdd with no shared-memory transpose to keep in step).
  if grep -q 'eval_vvars_gga_impl\|eval_vvars_mgga_impl' "$dst"; then
    awk '
      /^void eval_vvars_gga_impl\(/  { fn = "gga" }
      /^void eval_vvars_mgga_impl\(/ { fn = "mgga" }
      /^void eval_vvars_lda_impl\(/  { fn = "" }
      /^void eval_tmat_/             { fn = "" }
      /dim3 threads\( hip::warp_size, hip::max_warps_per_thread_block, 1 \);/ && (fn == "gga" || fn == "mgga") {
        print "  dim3 threads( 32, hip::max_warps_per_thread_block, 1 ); // tile/launch width pinned to 32, independent of hip::warp_size -- see hipify.sh"
        next
      }
      { print }
    ' "$dst" > "$dst.tmp"
    mv "$dst.tmp" "$dst"
  fi

  # exx_ek_screening_bfn_stats_kernel (same source file as the
  # bitvector_to_position_list_* kernels above, but a different kernel)
  # sizes its bf_shared[32][32+1]/bfn_sum_shared[32] tile off the literal
  # 32 that CUDA's warp_lane/warp_id/nwarp/warp_reduce_* all coincide
  # with, but otherwise consistently uses hip::warp_size (64) for those
  # same quantities. That mismatch is an out-of-bounds LDS write/read
  # (warp_lane reaches 63), not merely wrong numbers. Rather than resize
  # the tile to the true wavefront width (blockDim.x is hardcoded to 1024
  # at the call site, so the row count is a fixed constant either way),
  # pin warp_lane/warp_id/nwarp and the warp_reduce_* widths to 32 for
  # this one kernel -- each real 64-lane wavefront then behaves as two
  # independent 32-lane groups, matching CUDA's numbers exactly and the
  # tile's existing 32x33 sizing, consistent with the [32]-pinning already
  # used for this file's other kernel above. Scoped to this kernel's body
  # only (brace-depth tracked) so the true-64-wide uses of hip::warp_size
  # elsewhere in this file (exx_ek_collapse_fmax_to_shells_kernel, which
  # has no shared-memory tile, and the host-side launch geometry) are
  # untouched.
  if grep -q 'exx_ek_screening_bfn_stats_kernel' "$dst"; then
    awk '
      {
        line = $0
        if (!infn && line ~ /^__global__ void exx_ek_screening_bfn_stats_kernel\(/) {
          infn = 1; depth = 0; entered = 0
        }
        if (infn) {
          if (line ~ /hip::warp_size/)
            gsub(/hip::warp_size/, "32 /* pinned, independent of hip::warp_size -- see hipify.sh */", line)
          opens  = gsub(/\{/, "{", line)
          closes = gsub(/\}/, "}", line)
          if (opens > 0) entered = 1
          depth += opens - closes
        }
        print line
        if (infn && entered && depth <= 0) infn = 0
      }
    ' "$dst" > "$dst.tmp"
    mv "$dst.tmp" "$dst"
  fi

  # compute_grid_to_center_dist's inner loop stride is hardcoded to
  # hip::warp_size/2 (32), on the assumption that blockDim.y ==
  # warp_size/2 -- true on CUDA (distance_thread_y == cuda::warp_size/2
  # == 16 == blockDim.y) but not on HIP, where
  # distance_thread_y == hip::max_warps_per_thread_block/2 == 8 while the
  # stride is still 32: threadIdx.y in [0,8) only reaches k in
  # {0..7} u {32..39}, leaving 48 of every 64 points' dist[] entries
  # whatever was there before. Stride by the actual blockDim.y instead of
  # a hardcoded value derived from warp_size -- warp-width-agnostic and
  # correct on both backends (it reduces to the same 16 on CUDA).
  sed -i -e 's|for (int k = threadIdx\.y; k < hip::warp_size; k+=hip::warp_size/2) {|for (int k = threadIdx.y; k < hip::warp_size; k+=blockDim.y) { // stride by actual blockDim.y, not a hardcoded warp_size/2 -- see hipify.sh|' \
      "$dst"

  # HIP needs its runtime header; insert before the first #include (i.e.
  # after the license comment block).
  if ! grep -q 'hip/hip_runtime.h' "$dst"; then
    awk '
      BEGIN { done = 0 }
      /^#include/ && !done { print "#include \"hip/hip_runtime.h\""; done = 1 }
      { print }
      END { if (!done) print "#include \"hip/hip_runtime.h\"" }
    ' "$dst" > "$dst.tmp"
    mv "$dst.tmp" "$dst"
  fi
}

# ---- collocation ------------------------------------------------------
# Every header directly under collocation/ (per-l cartesian/spherical
# kernels for values/gradient/hessian/laplacian/lapgrad, plus the shared
# angular/radial/constants headers); deprecated/, scripts/ and templates/
# are CUDA-generator-only and are not part of the build.
for f in "$CUDA_PREFIX"/collocation/*.hpp; do
  hipify_file "$f" "$HIP_PREFIX/collocation/$(basename "$f")"
done

hipify_file "$CUDA_PREFIX/collocation_masked_combined_kernels.hpp" \
            "$HIP_PREFIX/collocation_masked_combined_kernels.hpp"
hipify_file "$CUDA_PREFIX/collocation_masked_kernels.hpp" \
            "$HIP_PREFIX/collocation_masked_kernels.hpp"
hipify_file "$CUDA_PREFIX/collocation_shell_to_task_kernels.hpp" \
            "$HIP_PREFIX/collocation_shell_to_task_kernels.hpp"
hipify_file "$CUDA_PREFIX/collocation_device.cu" \
            "$HIP_PREFIX/collocation_device.hip"

# ---- weights ----------------------------------------------------------
hipify_file "$CUDA_PREFIX/grid_to_center.cu"  "$HIP_PREFIX/grid_to_center.hip"
hipify_file "$CUDA_PREFIX/grid_to_center.hpp" "$HIP_PREFIX/grid_to_center.hpp"
hipify_file "$CUDA_PREFIX/cuda_ssf_1d.cu"     "$HIP_PREFIX/hip_ssf_1d.hip"
hipify_file "$CUDA_PREFIX/cuda_ssf_1d.hpp"    "$HIP_PREFIX/hip_ssf_1d.hpp"

# ---- BLAS extensions --------------------------------------------------
hipify_file "$CUDA_PREFIX/cublas_extensions.cu" "$HIP_PREFIX/hipblas_extensions.hip"
hipify_file "$CUDA_PREFIX/cuda_extensions.hpp"  "$HIP_PREFIX/hip_extensions.hpp"

# ---- density / potential / XC assembly --------------------------------
hipify_file "$CUDA_PREFIX/uvvars.cu"          "$HIP_PREFIX/uvvars.hip"
hipify_file "$CUDA_PREFIX/uvvars_lda.hpp"     "$HIP_PREFIX/uvvars_lda.hpp"
hipify_file "$CUDA_PREFIX/uvvars_gga.hpp"     "$HIP_PREFIX/uvvars_gga.hpp"
hipify_file "$CUDA_PREFIX/uvvars_mgga.hpp"    "$HIP_PREFIX/uvvars_mgga.hpp"
hipify_file "$CUDA_PREFIX/zmat_vxc.cu"        "$HIP_PREFIX/zmat_vxc.hip"
hipify_file "$CUDA_PREFIX/zmat_fxc.cu"        "$HIP_PREFIX/zmat_fxc.hip"
hipify_file "$CUDA_PREFIX/pack_submat.cu"     "$HIP_PREFIX/pack_submat.hip"
hipify_file "$CUDA_PREFIX/symmetrize_mat.cu"  "$HIP_PREFIX/symmetrize_mat.hip"
hipify_file "$CUDA_PREFIX/cuda_inc_potential.cu" "$HIP_PREFIX/hip_inc_potential.hip"

# ---- gradients and sn-LinK screening ----------------------------------
hipify_file "$CUDA_PREFIX/increment_exc_grad.cu" "$HIP_PREFIX/increment_exc_grad.hip"
hipify_file "$CUDA_PREFIX/exx_ek_screening_bfn_stats.cu" \
            "$HIP_PREFIX/exx_ek_screening_bfn_stats.hip"

echo "hipify: regenerated $(ls "$HIP_PREFIX"/*.hip "$HIP_PREFIX"/*.hpp | wc -l) files"
echo
echo "CHECK_WAVEFRONT: kernels performing intra-warp reductions were written"
echo "against a 32-lane warp.  On AMD the wavefront is 64 lanes and"
echo "GauXC::hip::warp_size reflects that, but any reduction whose trip"
echo "count is written as a literal must be reviewed.  Grep for '16;' '8;'"
echo "'4;' '2;' '1;' shuffle ladders in the generated files before trusting"
echo "numerical results on AMD hardware."
