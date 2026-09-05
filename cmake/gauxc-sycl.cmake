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
find_package( SYCL REQUIRED )

# Link against the SYCL runtime, but do NOT inherit its -fsycl interface
# compile option: that would run the SYCL frontend over every GauXC
# translation unit, including the large generated host-only sources where it
# costs a substantial fraction of their compile time and emits no device code.
# The flag is applied per source by gauxc_sycl_device_sources() instead, while
# the link still needs it so the device images are gathered. The imported
# target itself is left untouched -- ExchCXX consumes it and does want the
# target-wide behavior.
# $<LINK_ONLY:> only suppresses a dependency's usage requirements for the
# purposes of linking; INTERFACE_COMPILE_OPTIONS still propagate, and they
# reach GauXC both directly and transitively through ExchCXX. Clear -fsycl
# from the imported target itself, once ExchCXX has finished configuring and
# compiling its own sources with it. Without this the SYCL frontend runs over
# every GauXC translation unit, including the large generated host-only
# sources, adding a substantial fraction to their compile time while emitting
# no device code. gauxc_sycl_device_sources() applies the flag per source, and
# the link keeps it so the device images are still gathered.
get_target_property( _sycl_copts SYCL::SYCL INTERFACE_COMPILE_OPTIONS )
if( _sycl_copts )
  set( _sycl_copts_keep )
  foreach( _opt ${_sycl_copts} )
    if( NOT _opt MATCHES "(^|[:;])-fsycl($|[;>])" )
      list( APPEND _sycl_copts_keep ${_opt} )
    endif()
  endforeach()
  set_target_properties( SYCL::SYCL PROPERTIES
    INTERFACE_COMPILE_OPTIONS "${_sycl_copts_keep}" )
endif()

target_link_libraries( gauxc PUBLIC SYCL::SYCL )
target_link_options( gauxc PUBLIC -fsycl )

# oneMKL supplies the device BLAS for the SYCL backend. Only the BLAS domain
# is needed -- MKL::MKL_SYCL is an umbrella over every SYCL domain (LAPACK,
# DFT, RNG, sparse, stats, VM, data fitting), all of which would be dragged
# onto the link line. This mirrors how the CUDA backend takes CUDA::cublas
# rather than the whole toolkit.
#
# The host-side BLAS is a separate concern: it is located by
# gauxc-linalg-modules (BLAS::BLAS) when src/.../host is added, and must not be
# resolved here.
#
# Both halves must agree on the integer interface. oneMKL defaults to ILP64,
# and mixing the two pulls in libmkl_intel_lp64 and libmkl_intel_ilp64 at once;
# the duplicate dgemm_ then binds to whichever copy the linker sees first, so a
# 32-bit call site is read as 64-bit arguments. That surfaces at runtime as
# "Intel oneMKL ERROR: Parameter 10 was incorrect on entry to DGEMM".
# MKLConfig keeps a separate cache entry for the SYCL domains, so both the
# CPU-side (MKL_INTERFACE_FULL) and device-side (MKL_SYCL_INTERFACE_FULL)
# selections have to be pinned. Keyed off GAUXC_BLAS_PREFER_ILP64 rather than
# the discovered GAUXC_BLAS_IS_LP64, which is not set until later.
if( GAUXC_BLAS_PREFER_ILP64 )
  set( _gauxc_mkl_interface intel_ilp64 )
else()
  set( _gauxc_mkl_interface intel_lp64 )
endif()
set( MKL_INTERFACE_FULL      ${_gauxc_mkl_interface} CACHE STRING "" FORCE )
set( MKL_SYCL_INTERFACE_FULL ${_gauxc_mkl_interface} CACHE STRING "" FORCE )

if( NOT TARGET MKL::MKL_SYCL::BLAS )
  find_package( MKL CONFIG REQUIRED )
endif()

if( NOT TARGET MKL::MKL_SYCL::BLAS )
  message( FATAL_ERROR
    "oneMKL was found but does not provide MKL::MKL_SYCL::BLAS, which the "
    "GauXC SYCL backend requires for its device BLAS." )
endif()

# LINK_ONLY: MKLConfig attaches -fsycl to its SYCL targets as an interface
# compile option, which would push the SYCL frontend onto every consuming
# translation unit, including the large generated host-only sources where it
# adds substantial compile time and emits no device code. The flag is applied
# per source by gauxc_sycl_device_sources() below.
# $<LINK_ONLY:> suppresses a dependency's usage requirements for linking only;
# INTERFACE_COMPILE_OPTIONS still propagate. MKLConfig puts -fsycl there, so
# strip it from a private copy of the target rather than inherit it: the flag
# would otherwise reach every GauXC translation unit, including the large
# generated host-only sources, where it adds a substantial fraction to their
# compile time and produces no device code. The entry is wrapped in a
# COMPILE_LANGUAGE generator expression, hence the match rather than a string
# compare. gauxc_sycl_device_sources() applies -fsycl per source instead.
get_target_property( _mkl_blas_copts MKL::MKL_SYCL::BLAS INTERFACE_COMPILE_OPTIONS )
if( _mkl_blas_copts )
  set( _mkl_blas_copts_keep )
  foreach( _opt ${_mkl_blas_copts} )
    if( NOT _opt MATCHES "(^|[:;])-fsycl($|[;>])" )
      list( APPEND _mkl_blas_copts_keep ${_opt} )
    endif()
  endforeach()
  set_target_properties( MKL::MKL_SYCL::BLAS PROPERTIES
    INTERFACE_COMPILE_OPTIONS "${_mkl_blas_copts_keep}" )
endif()

target_link_libraries( gauxc PUBLIC $<LINK_ONLY:MKL::MKL_SYCL::BLAS> )

if( GAUXC_ENABLE_ONECCL )
  find_package( oneCCL REQUIRED )
  target_link_libraries( gauxc PUBLIC oneCCL )
  set( GAUXC_HAS_ONECCL TRUE CACHE BOOL "GauXC has oneCCL" FORCE )
endif()


# --- AoT SYCL target alias pass-through, mirrors ExchCXX's mapping ---
set( _GAUXC_SYCL_ALLOWED
  intel_gpu_pvc
  spir64_x86_64
  nvidia_gpu_sm_80
  nvidia_gpu_sm_90
  amd_gpu_gfx90a
  amd_gpu_gfx942
)

if( DEFINED GAUXC_SYCL_TARGET AND NOT GAUXC_SYCL_TARGET STREQUAL "" )
  list( FIND _GAUXC_SYCL_ALLOWED "${GAUXC_SYCL_TARGET}" _gauxc_sycl_idx )
  if( _gauxc_sycl_idx EQUAL -1 )
    message( FATAL_ERROR
      "Invalid GAUXC_SYCL_TARGET='${GAUXC_SYCL_TARGET}'. "
      "Allowed values: ${_GAUXC_SYCL_ALLOWED}" )
  endif()

  unset( _gauxc_sycl_compile_opts )
  unset( _gauxc_sycl_link_opts )

  if( GAUXC_SYCL_TARGET STREQUAL "intel_gpu_pvc" )
    # These land on a SOURCE COMPILE_OPTIONS property, which does not honor
    # the SHELL: prefix -- the whole string would be passed as one argument.
    # Keep each token separate so -Xsycl-target-backend gets its operand.
    list( APPEND _gauxc_sycl_compile_opts
      -fsycl-targets=spir64_gen
      -Xsycl-target-backend=spir64_gen
      "-device pvc"
    )
    list( APPEND _gauxc_sycl_link_opts
      "SHELL:-ftarget-register-alloc-mode=pvc:large"
      "SHELL:-fsycl-targets=spir64_gen"
      "SHELL:-Xsycl-target-backend \"-device pvc\""
    )
  elseif( GAUXC_SYCL_TARGET STREQUAL "spir64_x86_64" )
    list( APPEND _gauxc_sycl_compile_opts -fsycl-targets=spir64_x86_64 )
    list( APPEND _gauxc_sycl_link_opts "SHELL:-fsycl-targets=spir64_x86_64" )
  elseif( GAUXC_SYCL_TARGET MATCHES "^nvidia_gpu_" OR
          GAUXC_SYCL_TARGET MATCHES "^amd_gpu_" )
    list( APPEND _gauxc_sycl_compile_opts -fsycl-targets=${GAUXC_SYCL_TARGET} )
    list( APPEND _gauxc_sycl_link_opts "SHELL:-fsycl-targets=${GAUXC_SYCL_TARGET}" )
  endif()

  if( NOT _gauxc_sycl_compile_opts )
    message( FATAL_ERROR
      "GAUXC_SYCL_TARGET='${GAUXC_SYCL_TARGET}' is allowed but has no AoT flag "
      "mapping in ${CMAKE_CURRENT_LIST_FILE}; refusing to fall back to a JIT "
      "build that would be reported as AoT." )
  endif()

  # AoT compilation is applied per-source rather than to the whole gauxc
  # target: only the SYCL sources carry device code, and running the AoT
  # backend over the host-only translation units (notably the generated
  # Obara-Saika integrals) would cost a full rebuild for no device image.
  # gauxc_sycl_aot_sources() below stamps the flags onto the SYCL sources.
  # The link options are PUBLIC: gauxc is a static archive, so device linking
  # happens when the consuming executable is linked, and that link needs the
  # same target selection.
  set( GAUXC_SYCL_AOT_COMPILE_OPTIONS ${_gauxc_sycl_compile_opts} )
  target_link_options( gauxc PUBLIC ${_gauxc_sycl_link_opts} )

  message( STATUS "GauXC SYCL AoT enabled for target: ${GAUXC_SYCL_TARGET}" )
endif()

# Mark the given sources (relative to the calling CMakeLists) as carrying SYCL
# device code. Only these translation units are compiled with -fsycl and the
# device-code flags -- running the SYCL frontend over the host-only sources
# costs a substantial fraction of their compile time and produces no device
# image. When GAUXC_SYCL_TARGET requests it, the ahead-of-time flags are
# stamped here too.
function( gauxc_sycl_device_sources )
  # TARGET_DIRECTORY is required: SOURCE properties are scoped to the directory
  # that sets them, and these sources are added from subdirectories while the
  # gauxc target itself is defined in src/. Without it the property is set in a
  # scope the target never reads, and the flags silently never appear.
  foreach( _src ${ARGN} )
    get_filename_component( _abs ${_src} ABSOLUTE )
    set_property( SOURCE ${_abs} TARGET_DIRECTORY gauxc
      APPEND PROPERTY COMPILE_OPTIONS
        ${SYCL_DEVICE_COMPILE_OPTIONS}
        ${GAUXC_SYCL_DEVICE_COMPILE_OPTIONS}
        ${GAUXC_SYCL_AOT_COMPILE_OPTIONS} )
  endforeach()
endfunction()


include( CheckCXXCompilerFlag )
check_cxx_compiler_flag( "-fno-sycl-id-queries-fit-in-int"     GAUXC_SYCL_ID_QUERIES_FIT_IN_INT )
check_cxx_compiler_flag( "-fsycl-device-code-split=per_kernel" GAUXC_SYCL_DEVICE_CODE_SPLIT_PER_KERNEL )
check_cxx_compiler_flag( "-Xsycl-target-frontend \"-fp-model=precise\"" GAUXC_HAVE_SYCL_TARGET_FRONTEND_FP_MODEL_PRECISE )

include( CheckLinkerFlag )
check_linker_flag( CXX "-flink-huge-device-code"         GAUXC_SYCL_LINK_HUGE_DEVICE_CODE )
check_linker_flag( CXX "-fsycl-max-parallel-link-jobs=4" GAUXC_SYCL_MAX_PARALLEL_LINK_JOBS )

# GauXC indexes into task batches which routinely exceed INT_MAX
if( GAUXC_SYCL_ID_QUERIES_FIT_IN_INT )
  list( APPEND GAUXC_SYCL_DEVICE_COMPILE_OPTIONS -fno-sycl-id-queries-fit-in-int )
endif()

# The collocation kernels alone are ~50 distinct kernels; splitting keeps JIT
# time and device image size manageable
if( GAUXC_SYCL_DEVICE_CODE_SPLIT_PER_KERNEL )
  list( APPEND GAUXC_SYCL_DEVICE_COMPILE_OPTIONS -fsycl-device-code-split=per_kernel )
  target_link_options( gauxc PUBLIC
    $<$<LINK_LANGUAGE:CXX>:-fsycl-device-code-split=per_kernel>
  )
endif()

# Device-side FP model. The matching host-side -fp-model=precise is applied
# project-wide in the top-level CMakeLists.txt, so host and device agree.
if( GAUXC_HAVE_SYCL_TARGET_FRONTEND_FP_MODEL_PRECISE )
  # SHELL: is only honored for target-level COMPILE_OPTIONS; on a SOURCE
  # property the whole string is passed as one argument. Pass the two tokens
  # separately so -Xsycl-target-frontend picks up its operand.
  list( APPEND GAUXC_SYCL_DEVICE_COMPILE_OPTIONS
    -Xsycl-target-frontend -fp-model=precise )
endif()

if( GAUXC_SYCL_LINK_HUGE_DEVICE_CODE )
  target_link_options( gauxc PUBLIC
    $<$<LINK_LANGUAGE:CXX>:-flink-huge-device-code>
  )
endif()

if( GAUXC_SYCL_MAX_PARALLEL_LINK_JOBS )
  target_link_options( gauxc PUBLIC
    $<$<LINK_LANGUAGE:CXX>:-fsycl-max-parallel-link-jobs=4>
  )
endif()
