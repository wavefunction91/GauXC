find_package( LAPACK  REQUIRED )
include( gauxc-gau2grid     )
target_sources( gauxc PRIVATE 
  # Common Host Utilities
  host/host_weights.cxx
  host/host_collocation.cxx
  host/blas.cxx
  
  # XC Specific
  host/host_exc_vxc_zmat.cxx
  host/local_work_replicated_exc_vxc.cxx
  
  # Interfaces
  host/reference_xc_host_integrator.cxx
)

target_link_libraries( gauxc PUBLIC LAPACK::LAPACK )

if( GAUXC_ENABLE_GAU2GRID )
  target_link_libraries( gauxc PUBLIC gau2grid::gg )
endif()


