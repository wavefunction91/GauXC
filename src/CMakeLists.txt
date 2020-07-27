# Parallelism
find_package( OpenMP  REQUIRED )
find_package( MPI     REQUIRED )

# Required Dependencies

find_package( IntegratorXX QUIET )
if( NOT ${IntegratorXX_FOUND} )

  set( INTEGRATORXX_ENABLE_TESTS OFF CACHE BOOL "" )
  FetchContent_Declare(
    integratorxx
    GIT_REPOSITORY 
      https://github.com/wavefunction91/IntegratorXX.git
    GIT_TAG 
      redo
  )

  FetchContent_MakeAvailable( integratorxx )

endif()

find_package( ExchCXX QUIET )
if( NOT ${ExchCXX_FOUND} )

  set( EXCHCXX_ENABLE_CUDA  ${GAUXC_ENABLE_CUDA} CACHE BOOL "" )
  set( EXCHCXX_ENABLE_TESTS OFF                  CACHE BOOL "" )

  FetchContent_Declare(
    exchcxx
    GIT_REPOSITORY 
      https://github.com/wavefunction91/ExchCXX.git
    GIT_TAG 
      master
  )

  FetchContent_MakeAvailable( exchcxx )


else()

  if( ${GAUXC_ENABLE_CUDA} AND NOT ${EXCHCXX_ENABLE_CUDA} )
    message( FATAL_ERROR "GauXC CUDA BINDINGS REQUIRE ExchCXX CUDA Bindings" )
  endif()

endif()

if( GAUXC_ENABLE_CEREAL )
  include( BuildFindCereal )
endif()


# Optional Dependencies
find_package( gau2grid QUIET ) 
set( GAUXC_HAS_GAU2GRID ${gau2grid_FOUND} )


add_library( gauxc STATIC grid.cxx grid_impl.cxx molmeta.cxx molgrid.cxx molgrid_impl.cxx molgrid_defaults.cxx atomic_radii.cxx load_balancer_impl.cxx load_balancer.cxx default_load_balancer.cxx )

target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

target_link_libraries( gauxc PUBLIC ExchCXX::exchcxx IntegratorXX::integratorxx MPI::MPI_C)
if( GAUXC_ENABLE_CEREAL )
  target_link_libraries( gauxc PUBLIC cereal )
endif()
