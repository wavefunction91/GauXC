##target_sources( gauxc
##  PRIVATE
##    integrator_common.cxx
##    host/xc_host_util.cxx
##    host/host_weights.cxx
##    host/host_collocation.cxx
##    host/host_zmat.cxx
##    host/blas.cxx
##)
##
##target_include_directories( gauxc
##  PUBLIC
##    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator>
##)
##
##if( GAUXC_ENABLE_DEVICE )
##
##  if( GAUXC_ENABLE_CUDA )
##    include( cuda/gauxc-cuda_integrator.cmake )
##  endif()
##
##  get_target_property( DEVICE_TARGET_DEPENDS gauxc_device_integrator
##    INTERFACE_LINK_LIBRARIES ) 
##
##  target_link_libraries( gauxc PUBLIC 
##    $<TARGET_OBJECTS:gauxc_device_integrator> 
##    ${DEVICE_TARGET_DEPENDS}
##  )
##
##  add_dependencies( gauxc gauxc_device_integrator )
##
##endif()
#
#
#
#set( GAUXC_COMMON_INTEGRATOR_INCLUDES
#     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
#     $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
#     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
#     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator> )
#
#get_target_property( GAUXC_COMMON_INTEGRATOR_LIBRARIES 
#  gauxc INTERFACE_LINK_LIBRARIES )
#
#
#
#add_library( gauxc_integrator_common STATIC integrator_common.cxx )
#
#target_link_libraries( gauxc_integrator_common PUBLIC 
#                       "${GAUXC_COMMON_INTEGRATOR_LIBRARIES}" )
#target_include_directories( gauxc_integrator_common PUBLIC 
#                       "${GAUXC_COMMON_INTEGRATOR_INCLUDES}" )
#
#target_compile_features( gauxc_integrator_common PUBLIC cxx_std_17 )
#
#add_library( gauxc_integrator_host  STATIC  host/xc_host_util.cxx
#                                            host/host_weights.cxx
#                                            host/host_collocation.cxx
#                                            host/host_zmat.cxx
#                                            host/blas.cxx              )
#
#target_link_libraries( gauxc_integrator_host PUBLIC gauxc_integrator_common )
#
#set( GAUXC_INTEGRATOR_LIBRARIES gauxc_integrator_host )
#             
#
#target_link_libraries( gauxc PUBLIC gauxc_host_integrator )
#
#
#
#if( GAUXC_ENABLE_DEVICE )
#
#  if( GAUXC_ENABLE_CUDA )
#    include( cuda/gauxc-cuda_integrator.cmake )
#  endif()
#
#  #get_target_property( DEVICE_TARGET_DEPENDS gauxc_device_integrator
#  #  INTERFACE_LINK_LIBRARIES ) 
#
#  ##target_link_libraries( gauxc PUBLIC 
#  ##  $<TARGET_OBJECTS:gauxc_device_integrator> 
#  ##  ${DEVICE_TARGET_DEPENDS}
#  ##)
#
#  #list( PREPEND GAUXC_INTEGRATOR_LIBRARIES 
#  #      $<TARGET_OBJECTS:gauxc_device_integrator> )
#
#  #list( APPEND GAUXC_INTEGRATOR_LIBRARIES ${DEVICE_TARGET_DEPENDS} )
#  target_link_libraries( gauxc PUBLIC gauxc_device_integrator )
#
#
#endif()


# Common Integrator Utilities
target_sources( gauxc PRIVATE integrator_common.cxx )
target_include_directories( gauxc
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/integrator>
)

# Host Integrator Utilities
if( GAUXC_ENABLE_HOST )
  include( host/gauxc-host_integrator.cmake )
endif()

if( GAUXC_ENABLE_CUDA )
  include( cuda/gauxc-cuda_integrator.cmake )
endif()

if( GAUXC_ENABLE_HIP )
  include( hip/gauxc-hip_integrator.cmake )
endif()

if( GAUXC_ENABLE_SYCL )
  include( sycl/gauxc-sycl_integrator.cmake )
endif()
