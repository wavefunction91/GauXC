set( GAUXC_RYS_HOST_SRC 
     src/boys_table.c
     src/jacobi_table.c
     src/rys_1rw.c
     src/rys_2rw.c
     src/rys_3rw.c
     src/rys_4rw.c
     src/rys_5rw.c
     src/rys_integral.c
     src/rys_rw.c
     src/rys_xrw.c )

target_sources( gauxc PRIVATE ${GAUXC_RYS_HOST_SRC} )
target_include_directories( gauxc PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
)
