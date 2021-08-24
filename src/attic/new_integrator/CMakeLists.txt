# Implementations of generic interfaces
target_sources( gauxc PRIVATE replicated_xc_integrator_impl.cxx )

target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/new_integrator>
)

# Common Utilities
include( common/gauxc-common.cmake )

# Host Integrator Utilities
if( GAUXC_ENABLE_HOST )
  include( host/gauxc-host.cmake )
endif()

# Device Integrator Utilities
if( GAUXC_ENABLE_DEVICE )
  include( device/gauxc-device.cmake )
endif()
