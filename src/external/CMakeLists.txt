include(FetchContent)
FetchContent_Declare( HighFive
  GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
)

set(HIGHFIVE_USE_BOOST OFF CACHE BOOL "" )
set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "" )
set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "" )
#set(HIGHFIVE_PARALLEL_HDF5 ON CACHE BOOL "" )
set(HIGHFIVE_BUILD_DOCS OFF CACHE BOOL "" )
FetchContent_MakeAvailable( HighFive )

target_sources( gauxc PRIVATE hdf5.cxx )
target_link_libraries( gauxc PUBLIC HighFive )
