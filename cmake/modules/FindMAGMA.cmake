if( NOT DEFINED MAGMA_ROOT_DIR )
  find_package(PkgConfig)
  pkg_check_modules( PC_MAGMA magma )
endif()

if( NOT MAGMA_INCLUDE_DIR )
find_path( MAGMA_INCLUDE_DIR magma.h
           HINTS ${PC_MAGMA_INCLUDEDIR}
                 ${PC_MAGMA_INCLUDE_DIRS}
           PATHS ${MAGMA_ROOT_DIR}
           PATH_SUFFIXES include
)
endif()

if(NOT MAGMA_LIBRARY) 
find_library( MAGMA_LIBRARY NAMES magma
              HINTS ${PC_MAGMA_LIBDIR}
                    ${PC_MAGMA_LIBRARY_DIRS}
              PATHS ${MAGMA_ROOT_DIR}
              PATH_SUFFIXES lib lib64 lib32
)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( 
  MAGMA DEFAULT_MSG
  MAGMA_LIBRARY
  MAGMA_INCLUDE_DIR
)

if( MAGMA_FOUND AND NOT TARGET MAGMA::magma )

  set( MAGMA_INCLUDE_DIRS ${MAGMA_INCLUDE_DIR} )
  set( MAGMA_LIBRARIES    ${MAGMA_LIBRARY}     )

  add_library( MAGMA::magma INTERFACE IMPORTED )
  set_target_properties( MAGMA::magma PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${MAGMA_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES      "${MAGMA_LIBRARIES}"
  )

endif()

