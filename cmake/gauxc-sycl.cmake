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

target_link_libraries( gauxc PUBLIC SYCL::SYCL )
target_link_options( gauxc PUBLIC -fsycl )

# oneMKL supplies the device BLAS; the host BLAS is BLAS::BLAS, found
# separately. Host and device must agree on the integer interface or the
# duplicate dgemm_ binds to the wrong copy. Re-pinned from the host driver
# once BLAS_IS_LP64 is discovered.
function( gauxc_sycl_pin_mkl_interface _is_lp64 )
  if( _is_lp64 )
    set( _iface intel_lp64 )
  else()
    set( _iface intel_ilp64 )
  endif()
  set( MKL_INTERFACE_FULL      ${_iface} CACHE STRING "" FORCE )
  set( MKL_SYCL_INTERFACE_FULL ${_iface} CACHE STRING "" FORCE )
endfunction()

# Provisional; re-pinned from the discovered value in the host driver.
if( GAUXC_BLAS_PREFER_ILP64 )
  gauxc_sycl_pin_mkl_interface( FALSE )
else()
  gauxc_sycl_pin_mkl_interface( TRUE )
endif()

if( NOT TARGET MKL::MKL_SYCL::BLAS )
  find_package( MKL CONFIG REQUIRED )
endif()

if( NOT TARGET MKL::MKL_SYCL::BLAS )
  message( FATAL_ERROR
    "oneMKL was found but does not provide MKL::MKL_SYCL::BLAS, which the "
    "GauXC SYCL backend requires for its device BLAS." )
endif()

# MKLConfig also puts -fsycl in INTERFACE_COMPILE_OPTIONS; strip it for the
# same reason as SYCL::SYCL above. Wrapped in a generator expression, hence
# the regex match.
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


# AoT target aliases, mirroring ExchCXX's mapping
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
    # SOURCE COMPILE_OPTIONS does not honor SHELL:, so keep tokens separate.
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

  # Compile options are applied per source; link options are PUBLIC because
  # gauxc is a static archive and device linking happens in the consumer.
  set( GAUXC_SYCL_AOT_COMPILE_OPTIONS ${_gauxc_sycl_compile_opts} )
  target_link_options( gauxc PUBLIC ${_gauxc_sycl_link_opts} )

  message( STATUS "GauXC SYCL AoT enabled for target: ${GAUXC_SYCL_TARGET}" )
endif()

# Mark sources as carrying SYCL device code: only these get -fsycl and the
# device-code flags.
function( gauxc_sycl_device_sources )
  # TARGET_DIRECTORY is required: SOURCE properties are directory-scoped and
  # these sources are added from subdirectories, gauxc is defined in src/.
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

# ~50 collocation kernels alone; splitting keeps device image size manageable
if( GAUXC_SYCL_DEVICE_CODE_SPLIT_PER_KERNEL )
  list( APPEND GAUXC_SYCL_DEVICE_COMPILE_OPTIONS -fsycl-device-code-split=per_kernel )
  target_link_options( gauxc PUBLIC
    $<$<LINK_LANGUAGE:CXX>:-fsycl-device-code-split=per_kernel>
  )
endif()

# Device-side FP model; host side is set project-wide in the top-level file.
if( GAUXC_HAVE_SYCL_TARGET_FRONTEND_FP_MODEL_PRECISE )
  # Tokens kept separate: SOURCE properties do not honor SHELL:.
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
