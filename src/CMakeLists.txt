# Parallelism
find_package( OpenMP  REQUIRED )
find_package( MPI     REQUIRED )
find_package( LAPACK  REQUIRED )

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



find_package( Eigen3 CONFIG QUIET )
if( NOT Eigen3_FOUND )
  
  FetchContent_Declare(
    eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
  )

  FetchContent_GetPropertieS( eigen3 )
  if( NOT eigen3_POPULATED )
    FetchContent_Populate( eigen3 )
  endif()

  message( FATAL_ERROR "Eigen3 Pull Not Yet Configured" )

endif()



if( GAUXC_ENABLE_CEREAL )
  include( BuildFindCereal )
endif()


set( FETCHCONTENT_UPDATES_DISCONNECTED ON )
add_library( gauxc_gau2grid INTERFACE )

find_package( gau2grid CONFIG QUIET ) 
if( NOT gau2grid_FOUND )

  FetchContent_Declare(
    gau2grid
    GIT_REPOSITORY https://github.com/dgasmith/gau2grid.git
  )

  set( MAX_AM 6 CACHE STRING "" )
  FetchContent_GetProperties( gau2grid )
  if( NOT gau2grid_POPULATED )
    FetchContent_Populate( gau2grid )
    file( READ ${gau2grid_SOURCE_DIR}/CMakeLists.txt GAU2GRID_CMAKE_LISTS )
    string( REPLACE "CMAKE_SOURCE_DIR" "PROJECT_SOURCE_DIR" GAU2GRID_CORRECT "${GAU2GRID_CMAKE_LISTS}" )
    file( WRITE ${gau2grid_SOURCE_DIR}/CMakeLists.txt "${GAU2GRID_CORRECT}" )
    add_subdirectory( ${gau2grid_SOURCE_DIR} ${gau2grid_BINARY_DIR} )
  endif()



  target_link_libraries( gauxc_gau2grid INTERFACE gg )
  target_include_directories( gauxc_gau2grid INTERFACE ${gau2grid_BINARY_DIR} )

else()

  target_link_libraries( gauxc_gau2grid INTERFACE gau2grid::gg )

endif()



add_library( gauxc STATIC 
  grid.cxx 
  grid_impl.cxx 
  molmeta.cxx 
  molgrid.cxx 
  molgrid_impl.cxx 
  molgrid_defaults.cxx 
  atomic_radii.cxx 
  load_balancer_impl.cxx 
  load_balancer.cxx 
  replicated_load_balancer.cxx 
  integrator_common.cxx
  host_integrator/xc_host_util.cxx
  host_integrator/host_weights.cxx
  host_integrator/host_collocation.cxx
  host_integrator/host_zmat.cxx
  host_integrator/blas.cxx
)

target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
)

configure_file( 
  ${PROJECT_SOURCE_DIR}/include/gauxc/gauxc_config.hpp.in
  ${PROJECT_BINARY_DIR}/include/gauxc/gauxc_config.hpp
)

target_link_libraries( gauxc PUBLIC 
  ExchCXX::exchcxx 
  IntegratorXX::integratorxx 
  MPI::MPI_C 
  Eigen3::Eigen 
  LAPACK::lapack
  gauxc_gau2grid
)
if( GAUXC_ENABLE_CEREAL )
  target_link_libraries( gauxc PUBLIC cereal )
endif()
