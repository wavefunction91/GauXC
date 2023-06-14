# Implementations of generic interfaces
target_sources( gauxc PRIVATE replicated/replicated_xc_integrator_impl.cxx )

if( GAUXC_ENABLE_HOST )
  target_sources( gauxc PRIVATE replicated/reference_xc_host_integrator.cxx )
endif()

