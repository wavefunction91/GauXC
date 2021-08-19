target_sources( gauxc PRIVATE host/replicated_load_balancer.cxx )

#if( GAUXC_ENABLE_CUDA )
#  include( cuda/gauxc-cuda_load_balancer.cmake )
#endif()
