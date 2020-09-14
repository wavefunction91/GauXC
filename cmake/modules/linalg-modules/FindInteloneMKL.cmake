#    FindInteloneMKL.cmake
#
#    Finds Intel(R) MKL and exports its linkange as
#    CMake TARGETS
#
#    This module is meant to serve as part of FindLinAlg.
#    It can also be used by itself.
#
#    The module will define the following variables:
#
#      InteloneMKL_FOUND       - Found oneMKL installation
#      InteloneMKL_INCLUDE_DIR - Location of oneMKL headers (mkl_sycl.hpp)
#      InteloneMKL_LIBRARIES   - oneMKL libraries
#
#    This module will export the following CMake TARGETS if possible
#
#      InteloneMKL::mkl
#
#      intelonemkl_PREFERS_STATIC          - default OFF
#      intelonemkl_PREFERED_THREAD_LEVEL   - ( sequential, openmp, tbb ) default: openmp
#      intelonemkl_PREFERED_THREAD_LIBRARY - ( intel, gnu, pgi )         default: depends on compiler


# SANITY CHECK
if( "ilp64" IN_LIST Intelonemkl_FIND_COMPONENTS AND "lp64" IN_LIST Intelonemkl_FIND_COMPONENTS )
  message( FATAL_ERROR "InteloneMKL cannot link to both ILP64 and LP64 iterfaces" )
endif()

if( "scalapack" IN_LIST Intelonemkl_FIND_COMPONENTS AND NOT ("blacs" IN_LIST Intelonemkl_FIND_COMPONENTS) )
  list(APPEND Intelonemkl_FIND_COMPONENTS "blacs" )
endif()

# MKL lib names
if( intelonemkl_PREFERS_STATIC )
  set( intelonemkl_SYCL_LIBRARY_NAME       "libmkl_sycl.a"         )
  set( intelonemkl_LP64_LIBRARY_NAME       "libmkl_intel_lp64.a"   )
  set( intelonemkl_ILP64_LIBRARY_NAME      "libmkl_intel_ilp64.a"  )
  set( intelonemkl_SEQUENTIAL_LIBRARY_NAME "libmkl_sequential.a"   )
  set( intelonemkl_OMP_INTEL_LIBRARY_NAME  "libmkl_intel_thread.a" )
  set( intelonemkl_OMP_GNU_LIBRARY_NAME    "libmkl_gnu_thread.a"   )
  set( intelonemkl_OMP_PGI_LIBRARY_NAME    "libmkl_pgi_thread.a"   )
  set( intelonemkl_TBB_LIBRARY_NAME        "libmkl_tbb_thread.a"   )
  set( intelonemkl_CORE_LIBRARY_NAME       "libmkl_core.a"         )

  set( intelonemkl_LP64_SCALAPACK_LIBRARY_NAME  "libmkl_scalapack_lp64.a"  )
  set( intelonemkl_ILP64_SCALAPACK_LIBRARY_NAME "libmkl_scalapack_ilp64.a" )

  set( intelonemkl_LP64_INTELMPI_BLACS_LIBRARY_NAME  "libmkl_blacs_intelmpi_lp64.a"  )
  set( intelonemkl_LP64_OPENMPI_BLACS_LIBRARY_NAME   "libmkl_blacs_openmpi_lp64.a"   )
  set( intelonemkl_LP64_SGIMPT_BLACS_LIBRARY_NAME    "libmkl_blacs_sgimpt_lp64.a"    )
  set( intelonemkl_ILP64_INTELMPI_BLACS_LIBRARY_NAME "libmkl_blacs_intelmpi_ilp64.a" )
  set( intelonemkl_ILP64_OPENMPI_BLACS_LIBRARY_NAME  "libmkl_blacs_openmpi_ilp64.a"  )
  set( intelonemkl_ILP64_SGIMPT_BLACS_LIBRARY_NAME   "libmkl_blacs_sgimpt_ilp64.a"   )
else()
  set( intelonemkl_SYCL_LIBRARY_NAME       "mkl_sycl"         )
  set( intelonemkl_LP64_LIBRARY_NAME       "mkl_intel_lp64"   )
  set( intelonemkl_ILP64_LIBRARY_NAME      "mkl_intel_ilp64"  )
  set( intelonemkl_SEQUENTIAL_LIBRARY_NAME "mkl_sequential"   )
  set( intelonemkl_OMP_INTEL_LIBRARY_NAME  "mkl_intel_thread" )
  set( intelonemkl_OMP_GNU_LIBRARY_NAME    "mkl_gnu_thread"   )
  set( intelonemkl_OMP_PGI_LIBRARY_NAME    "mkl_pgi_thread"   )
  set( intelonemkl_TBB_LIBRARY_NAME        "mkl_tbb_thread"   )
  set( intelonemkl_CORE_LIBRARY_NAME       "mkl_core"         )

  set( intelonemkl_LP64_SCALAPACK_LIBRARY_NAME  "mkl_scalapack_lp64"  )
  set( intelonemkl_ILP64_SCALAPACK_LIBRARY_NAME "mkl_scalapack_ilp64" )

  set( intelonemkl_LP64_INTELMPI_BLACS_LIBRARY_NAME  "mkl_blacs_intelmpi_lp64"  )
  set( intelonemkl_LP64_OPENMPI_BLACS_LIBRARY_NAME   "mkl_blacs_openmpi_lp64"   )
  set( intelonemkl_LP64_SGIMPT_BLACS_LIBRARY_NAME    "mkl_blacs_sgimpt_lp64"    )
  set( intelonemkl_ILP64_INTELMPI_BLACS_LIBRARY_NAME "mkl_blacs_intelmpi_ilp64" )
  set( intelonemkl_ILP64_OPENMPI_BLACS_LIBRARY_NAME  "mkl_blacs_openmpi_ilp64"  )
  set( intelonemkl_ILP64_SGIMPT_BLACS_LIBRARY_NAME   "mkl_blacs_sgimpt_ilp64"   )
endif()


# Defaults
if( NOT intelonemkl_PREFERED_THREAD_LEVEL )
  set( intelonemkl_PREFERED_THREAD_LEVEL "openmp" )
endif()

if( NOT intelonemkl_PREFERED_MPI_LIBRARY )
  set( intelonemkl_PREFERED_MPI_LIBRARY "intelmpi" )
endif()

if( NOT intelonemkl_PREFIX )
  set( intelonemkl_PREFIX $ENV{MKLROOT} )
endif()



# MKL Threads
if( intelonemkl_PREFERED_THREAD_LEVEL MATCHES "sequential" )
  set( intelonemkl_THREAD_LIBRARY_NAME ${intelonemkl_SEQUENTIAL_LIBRARY_NAME} )
elseif( intelonemkl_PREFERED_THREAD_LEVEL MATCHES "tbb" )
  set( intelonemkl_THREAD_LIBRARY_NAME ${intelonemkl_TBB_LIBRARY_NAME} )
else() # OpenMP
  if( CMAKE_C_COMPILER_ID MATCHES "Intel" )
    set( intelonemkl_THREAD_LIBRARY_NAME ${intelonemkl_OMP_INTEL_LIBRARY_NAME} )
  elseif( CMAKE_C_COMPILER_ID MATCHES "PGI" )
    set( intelonemkl_THREAD_LIBRARY_NAME ${intelonemkl_OMP_PGI_LIBRARY_NAME} )
  else()
    set( intelonemkl_THREAD_LIBRARY_NAME ${intelonemkl_OMP_GNU_LIBRARY_NAME} )
  endif()
endif()


# MKL MPI for BLACS
if( intelonemkl_PREFERED_MPI_LIBRARY MATCHES "openmpi" )
  set( intelonemkl_LP64_BLACS_LIBRARY_NAME  ${intelonemkl_LP64_OPENMPI_BLACS_LIBRARY_NAME}  )
  set( intelonemkl_ILP64_BLACS_LIBRARY_NAME ${intelonemkl_ILP64_OPENMPI_BLACS_LIBRARY_NAME} )
elseif( intelonemkl_PREFERED_MPI_LIBRARY MATCHES "sgimpt" )
  set( intelonemkl_LP64_BLACS_LIBRARY_NAME  ${intelonemkl_LP64_SGIMPT_BLACS_LIBRARY_NAME}  )
  set( intelonemkl_ILP64_BLACS_LIBRARY_NAME ${intelonemkl_ILP64_SGIMPT_BLACS_LIBRARY_NAME} )
else() # Intel MPI
  set( intelonemkl_LP64_BLACS_LIBRARY_NAME  ${intelonemkl_LP64_INTELMPI_BLACS_LIBRARY_NAME}  )
  set( intelonemkl_ILP64_BLACS_LIBRARY_NAME ${intelonemkl_ILP64_INTELMPI_BLACS_LIBRARY_NAME} )
endif()


# Header
find_path( InteloneMKL_INCLUDE_DIR
  NAMES mkl_sycl.hpp
  HINTS ${intelonemkl_PREFIX}
  PATHS ${InteloneMKL_INCLUDE_DIR}
  PATH_SUFFIXES include
  DOC "Intel(R) MKL header"
)

find_library( intelonemkl_THREAD_LIBRARY
  NAMES ${intelonemkl_THREAD_LIBRARY_NAME}
  HINTS ${intelonemkl_PREFIX}
  PATHS ${intelonemkl_LIBRARY_DIR} ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
  PATH_SUFFIXES lib/intel64 lib/ia32
  DOC "Intel(R) MKL THREAD Library"
)

find_library( intelonemkl_CORE_LIBRARY
  NAMES ${intelonemkl_CORE_LIBRARY_NAME}
  HINTS ${intelonemkl_PREFIX}
  PATHS ${intelonemkl_LIBRARY_DIR} ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
  PATH_SUFFIXES lib/intel64 lib/ia32
  DOC "Intel(R) MKL CORE Library"
)



# Check version
if( EXISTS ${InteloneMKL_INCLUDE_DIR}/mkl_version.h )
  set( version_pattern
  "^#define[\t ]+__INTEL_MKL(|_MINOR|_UPDATE)__[\t ]+([0-9\\.]+)$"
  )
  file( STRINGS ${InteloneMKL_INCLUDE_DIR}/mkl_version.h mkl_version
        REGEX ${version_pattern} )

  foreach( match ${mkl_version} )

    if(Intelonemkl_VERSION_STRING)
      set(Intelonemkl_VERSION_STRING "${Intelonemkl_VERSION_STRING}.")
    endif()

    string(REGEX REPLACE ${version_pattern}
      "${Intelonemkl_VERSION_STRING}\\2"
      Intelonemkl_VERSION_STRING ${match}
    )

    set(Intelonemkl_VERSION_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})

  endforeach()

  unset( mkl_version )
  unset( version_pattern )
endif()



if( InteloneMKL_INCLUDE_DIR )
  set( Intelonemkl_INCLUDE_DIR ${InteloneMKL_INCLUDE_DIR} )
endif()


# Handle LP64 / ILP64
find_library( intelonemkl_ILP64_LIBRARY
  NAMES ${intelonemkl_ILP64_LIBRARY_NAME}
  HINTS ${intelonemkl_PREFIX}
  PATHS ${intelonemkl_LIBRARY_DIR} ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
  PATH_SUFFIXES lib/intel64 lib/ia32
  DOC "Intel(R) ILP64 MKL Library"
)

find_library( intelonemkl_LP64_LIBRARY
  NAMES ${intelonemkl_LP64_LIBRARY_NAME}
  HINTS ${intelonemkl_PREFIX}
  PATHS ${intelonemkl_LIBRARY_DIR} ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
  PATH_SUFFIXES lib/intel64 lib/ia32
  DOC "Intel(R) LP64 MKL Library"
)

if( intelonemkl_ILP64_LIBRARY )
  set( Intelonemkl_ilp64_FOUND TRUE )
endif()

if( intelonemkl_LP64_LIBRARY )
  set( Intelonemkl_lp64_FOUND TRUE )
endif()



# BLACS / ScaLAPACK

find_library( intelonemkl_ILP64_BLACS_LIBRARY
  NAMES ${intelonemkl_ILP64_BLACS_LIBRARY_NAME}
  HINTS ${intelonemkl_PREFIX}
  PATHS ${intelonemkl_LIBRARY_DIR} ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
  PATH_SUFFIXES lib/intel64 lib/ia32
  DOC "Intel(R) ILP64 MKL BLACS Library"
)

find_library( intelonemkl_LP64_BLACS_LIBRARY
  NAMES ${intelonemkl_LP64_BLACS_LIBRARY_NAME}
  HINTS ${intelonemkl_PREFIX}
  PATHS ${intelonemkl_LIBRARY_DIR} ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
  PATH_SUFFIXES lib/intel64 lib/ia32
  DOC "Intel(R) LP64 MKL BLACS Library"
)

find_library( intelonemkl_ILP64_SCALAPACK_LIBRARY
  NAMES ${intelonemkl_ILP64_SCALAPACK_LIBRARY_NAME}
  HINTS ${intelonemkl_PREFIX}
  PATHS ${intelonemkl_LIBRARY_DIR} ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
  PATH_SUFFIXES lib/intel64 lib/ia32
  DOC "Intel(R) ILP64 MKL SCALAPACK Library"
)

find_library( intelonemkl_LP64_SCALAPACK_LIBRARY
  NAMES ${intelonemkl_LP64_SCALAPACK_LIBRARY_NAME}
  HINTS ${intelonemkl_PREFIX}
  PATHS ${intelonemkl_LIBRARY_DIR} ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
  PATH_SUFFIXES lib/intel64 lib/ia32
  DOC "Intel(R) LP64 MKL SCALAPACK Library"
)



# Default to LP64
if( "ilp64" IN_LIST Intelonemkl_FIND_COMPONENTS )
  set( Intelonemkl_COMPILE_DEFINITIONS "MKL_ILP64" )
  if( CMAKE_C_COMPILER_ID MATCHES "GNU" )
    set( Intelonemkl_C_COMPILE_FLAGS        "-m64" )
    set( Intelonemkl_Fortran_COMPILE_FLAGS  "-m64" "-fdefault-integer-8" )
  elseif( CMAKE_C_COMPILER_ID MATCHES "PGI" )
    set( Intelonemkl_Fortran_COMPILE_FLAGS "-i8" )
  endif()
  set( intelonemkl_LIBRARY ${intelonemkl_ILP64_LIBRARY} )

  if( intelonemkl_ILP64_BLACS_LIBRARY )
    set( intelonemkl_BLACS_LIBRARY ${intelonemkl_ILP64_BLACS_LIBRARY} )
    set( Intelonemkl_blacs_FOUND TRUE )
  endif()

  if( intelonemkl_ILP64_SCALAPACK_LIBRARY )
    set( intelonemkl_SCALAPACK_LIBRARY ${intelonemkl_ILP64_SCALAPACK_LIBRARY} )
    set( Intelonemkl_scalapack_FOUND TRUE )
  endif()

else()
  set( intelonemkl_LIBRARY ${intelonemkl_LP64_LIBRARY} )

  if( intelonemkl_LP64_BLACS_LIBRARY )
    set( intelonemkl_BLACS_LIBRARY ${intelonemkl_LP64_BLACS_LIBRARY} )
    set( Intelonemkl_blacs_FOUND TRUE )
  endif()

  if( intelonemkl_LP64_SCALAPACK_LIBRARY )
    set( intelonemkl_SCALAPACK_LIBRARY ${intelonemkl_LP64_SCALAPACK_LIBRARY} )
    set( Intelonemkl_scalapack_FOUND TRUE )
  endif()
endif()





# Check if found library is actually static
if( intelonemkl_CORE_LIBRARY MATCHES ".+libmkl_core.a" )
  set( intelonemkl_PREFERS_STATIC TRUE )
endif()




if( intelonemkl_LIBRARY AND intelonemkl_THREAD_LIBRARY AND intelonemkl_CORE_LIBRARY )

  if( intelonemkl_PREFERS_STATIC )

    if( "scalapack" IN_LIST Intelonemkl_FIND_COMPONENTS )
      set( InteloneMKL_LIBRARIES ${intelonemkl_SCALAPACK_LIBRARY} )
    endif()

    list( APPEND InteloneMKL_LIBRARIES  "-Wl,--start-group" ${intelonemkl_LIBRARY} ${intelonemkl_THREAD_LIBRARY} ${intelonemkl_CORE_LIBRARY} )

    if( "blacs" IN_LIST Intelonemkl_FIND_COMPONENTS )
      list( APPEND InteloneMKL_LIBRARIES ${intelonemkl_BLACS_LIBRARY} )
    endif()

    list( APPEND InteloneMKL_LIBRARIES "-Wl,--end-group" )

  else()

    set( InteloneMKL_LIBRARIES "-Wl,--no-as-needed" )
    if( "scalapack" IN_LIST Intelonemkl_FIND_COMPONENTS )
      list( APPEND InteloneMKL_LIBRARIES ${intelonemkl_SCALAPACK_LIBRARY} )
    endif()

    list( APPEND InteloneMKL_LIBRARIES  ${intelonemkl_LIBRARY} ${intelonemkl_THREAD_LIBRARY} ${intelonemkl_CORE_LIBRARY} )

    if( "blacs" IN_LIST Intelonemkl_FIND_COMPONENTS )
      list( APPEND InteloneMKL_LIBRARIES ${intelonemkl_BLACS_LIBRARY} )
    endif()

  endif()


  if( intelonemkl_PREFERED_THREAD_LEVEL MATCHES "openmp" )
    find_package( OpenMP QUIET )
    set( InteloneMKL_LIBRARIES ${InteloneMKL_LIBRARIES} OpenMP::OpenMP_C )
  elseif( intelonemkl_PREFERED_THREAD_LEVEL MATCHES "tbb" )
    find_package( TBB QUIET )
    set( InteloneMKL_LIBRARIES ${InteloneMKL_LIBRARIES} tbb )
  endif()
  set( InteloneMKL_LIBRARIES ${InteloneMKL_LIBRARIES} "m" "dl" )

  if( "scalapack" IN_LIST Intelonemkl_FIND_COMPONENTS OR
      "blacs"     IN_LIST Intelonemkl_FIND_COMPONENTS )

    if( NOT TARGET MPI::MPI_C )
      find_dependency( MPI )
    endif()
    list( APPEND InteloneMKL_LIBRARIES MPI::MPI_C )

  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( InteloneMKL
  REQUIRED_VARS InteloneMKL_LIBRARIES Intelonemkl_INCLUDE_DIR
  VERSION_VAR Intelonemkl_VERSION_STRING
  HANDLE_COMPONENTS
)

if( InteloneMKL_FOUND AND NOT TARGET InteloneMKL::mkl )

  add_library( InteloneMKL::mkl INTERFACE IMPORTED )
  set_target_properties( InteloneMKL::mkl PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${InteloneMKL_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES      "${InteloneMKL_LIBRARIES}"
    INTERFACE_COMPILE_OPTIONS     "${InteloneMKL_C_COMPILE_FLAGS}"
    INTERFACE_COMPILE_DEFINITIONS "${InteloneMKL_COMPILE_DEFINITIONS}"
  )

  if( "scalapack" IN_LIST Intelonemkl_FIND_COMPONENTS AND NOT scalapack_LIBRARIES )
    set( scalapack_LIBRARIES InteloneMKL::mkl )
  endif()

  if( "blacs" IN_LIST Intelonemkl_FIND_COMPONENTS AND NOT blacs_LIBRARIES )
    set( blacs_LIBRARIES InteloneMKL::mkl )
  endif()

endif()
