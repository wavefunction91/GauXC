# Parallelism
find_package( OpenMP  REQUIRED )
find_package( MPI     REQUIRED )
find_package( LAPACK  REQUIRED )

# Required Dependencies
include( gauxc-integratorxx )
include( gauxc-exchcxx      )
include( gauxc-cereal       )
include( gauxc-eigen3       )
include( gauxc-gau2grid     )


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
