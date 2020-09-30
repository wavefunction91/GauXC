# Parallelism
find_package( OpenMP  REQUIRED )
find_package( MPI     REQUIRED )

# Required Dependencies
include( gauxc-integratorxx )
include( gauxc-exchcxx      )
include( gauxc-cereal       )
include( gauxc-eigen3       )


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
)

target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
)

target_compile_options( gauxc
  PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>: -Wall -Wextra -Wpedantic -Wnon-virtual-dtor>
)



configure_file( 
  ${PROJECT_SOURCE_DIR}/include/gauxc/gauxc_config.hpp.in
  ${PROJECT_BINARY_DIR}/include/gauxc/gauxc_config.hpp
)

target_link_libraries( gauxc PUBLIC 
  ExchCXX::exchcxx 
  IntegratorXX::integratorxx 
  MPI::MPI_C 
  OpenMP::OpenMP_CXX
  Eigen3::Eigen 
)

if( GAUXC_ENABLE_CEREAL )
  target_link_libraries( gauxc PUBLIC cereal )
endif()

add_subdirectory( integrator )
