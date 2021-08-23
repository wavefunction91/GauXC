
find_package( hip REQUIRED )
find_package( hipblas REQUIRED )

target_sources( gauxc PRIVATE hip_backend.cxx )
target_link_libraries( gauxc PUBLIC hip::host roc::hipblas )
