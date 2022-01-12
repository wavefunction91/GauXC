find_package( LAPACK  REQUIRED )
include( gauxc-gau2grid     )

target_sources( gauxc PRIVATE 
  local_host_work_driver.cxx
  local_host_work_driver_pimpl.cxx
  reference_local_host_work_driver.cxx

  reference/weights.cxx
  reference/gau2grid_collocation.cxx

  blas.cxx
)

target_link_libraries( gauxc PUBLIC LAPACK::LAPACK )

if( GAUXC_ENABLE_GAU2GRID AND TARGET gau2grid::gg )
  target_link_libraries( gauxc PUBLIC gau2grid::gg )
endif()

add_subdirectory(rys)
add_subdirectory(obara_saika)
