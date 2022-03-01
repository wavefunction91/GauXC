cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(gau2grid
        VERSION 2.0.5
        LANGUAGES C)
set(gau2grid_AUTHORS      "Daniel G. A. Smith")
set(gau2grid_DESCRIPTION  "Fast computation of a gaussian and its derivative on a grid")
set(gau2grid_URL          "https://github.com/dgasmith/gau2grid")
set(gau2grid_LICENSE      "BSD 3-clause")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


#############################  Options: Build How?  #############################
include(psi4OptionsTools)
option_with_default(MAX_AM "The maximum gaussian angular momentum to compile" 8)
option_with_default(CMAKE_BUILD_TYPE "Build type (Release or Debug)" Release)
if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
option_with_flags(ENABLE_XHOST "Enables processor-specific optimization (with MSVC, it enables AVX2 instructions)" ON
                  "-xHost" "-march=native" "/arch:AVX2")
else()
option_with_flags(ENABLE_XHOST "Enables processor-specific optimization (with MSVC, it enables AVX2 instructions)" ON
                  "-march=native" "-xHost" "/arch:AVX2")
endif()
option_with_default(BUILD_FPIC "Libraries will be compiled with position independent code" ON)
option_with_print(BUILD_SHARED_LIBS "Build final library as shared, not static" ON)
option_with_default(ENABLE_GENERIC "Enables mostly static linking of system libraries for shared library" OFF)
option_with_default(DISABLE_PRAGMA "Disable certain pragma optimizations, appends _GG_NO_PRAGMA to compile flags" OFF )

# Warnings
if((${BUILD_SHARED_LIBS}) AND NOT ${BUILD_FPIC})
    message(FATAL_ERROR "BUILD_SHARED_LIBS ON and BUILD_FPIC OFF are incompatible, as shared library requires position independent code")
endif()

# Install
option_with_default(CMAKE_INSTALL_LIBDIR "Directory to which libraries installed" lib)
option_with_default(PYMOD_INSTALL_LIBDIR "Location within CMAKE_INSTALL_LIBDIR to which python modules are installed
                                          Must start with: / . Used to imitate python install: /python3.6/site-packages ." /)
option_with_print(INSTALL_PYMOD "Additionally installs as independent python module in PYMOD_INSTALL_LIBDIR" OFF)
option_with_default(NATIVE_PYTHON_INSTALL "For INSTALL_PYMOD=ON, install in Python manner to PYTHON_EXECUTABLE's site-packages rather than Linux manner to prefix. Overrides CMAKE_INSTALL_PREFIX, CMAKE_INSTALL_LIBDIR, PYMOD_INSTALL_LIBDIR. Only Py module installed." OFF)
option_with_print(NATIVE_PYTHON_INSTALL_WITH_LIB "Same as NATIVE_PYTHON_INSTALL except installs library, too, _without_ overriding CMAKE_INSTALL_* options." OFF)

########################  Process & Validate Options  ##########################
include(autocmake_safeguards)
include(custom_color_messages)
include(custom_static_library)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/usr/local/gau2grid" CACHE PATH "Install path" FORCE)
endif()
message(STATUS "gau2grid install: ${CMAKE_INSTALL_PREFIX}")

#  <<  Python  >>
set(Python_ADDITIONAL_VERSIONS 3.9 3.8 3.7 3.6 3.5)  # adjust with CMake minimum FindPythonInterp
find_package(PythonLibsNew 3.6 REQUIRED)
message(STATUS "${Cyan}Found Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}${ColourReset}: ${PYTHON_EXECUTABLE} (found version ${PYTHON_VERSION_STRING})")


################################  Main Project  ################################
add_custom_command(
    OUTPUT  gau2grid/gau2grid.h gau2grid_orbital.c gau2grid_phi.c gau2grid_deriv1.c gau2grid_deriv2.c gau2grid_deriv3.c gau2grid_transform.c gau2grid_helper.c
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys; \
                                     sys.path.append('${PROJECT_SOURCE_DIR}'); \
                                     import gau2grid as gg; \
                                     gg.c_gen.generate_c_gau2grid(${MAX_AM}, path='${CMAKE_CURRENT_BINARY_DIR}')"
    DEPENDS gau2grid/c_generator.py
            gau2grid/c_generator.py
            gau2grid/codegen.py
            gau2grid/c_pragma.py
            gau2grid/c_util_generator.py
            gau2grid/c_wrapper.py
            gau2grid/docs_generator.py
            gau2grid/order.py
            gau2grid/python_reference.py
            gau2grid/RSH.py
            gau2grid/utility.py
    VERBATIM)

set(sources_list ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_phi.c
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_orbital.c
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_deriv1.c
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_deriv2.c
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_deriv3.c
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_transform.c
                 ${CMAKE_CURRENT_BINARY_DIR}/gau2grid_helper.c)

add_library(gg ${sources_list})
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    set_target_properties(gg PROPERTIES COMPILE_FLAGS "-c11")
else()
    set_target_properties(gg PROPERTIES COMPILE_FLAGS "-std=c11")
endif()
set_target_properties(gg PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC}
                                    SOVERSION 2)  # bump whenever interface has changes or removals

if( DISABLE_PRAGMA )
  target_compile_definitions( gg PRIVATE $<BUILD_INTERFACE:__GG_NO_PRAGMA> )
endif()

find_package(StandardMathLibraryC)
target_link_libraries(gg PRIVATE ${STANDARD_MATH_LIBRARY})

if(${BUILD_SHARED_LIBS})
    target_link_libraries(gg PRIVATE ${LIBC_INTERJECT})
endif()


###################################  Install  ##################################
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(PN ${PROJECT_NAME})

# Alias to allow for consistent manipulation as a subproject
add_library( ${PN}::gg ALIAS gg )

target_include_directories(gg PUBLIC
                           $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# GNUInstallDirs "DATADIR" wrong here; CMake search path wants "share".
set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}")
configure_package_config_file(cmake/${PN}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
                                 VERSION ${${PN}_VERSION}
                                 COMPATIBILITY SameMajorVersion)

# Install our files
if(${NATIVE_PYTHON_INSTALL_WITH_LIB} OR (NOT(${INSTALL_PYMOD} AND ${NATIVE_PYTHON_INSTALL})))
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gau2grid/gau2grid.h
                  ${CMAKE_CURRENT_BINARY_DIR}/gau2grid/gau2grid_pragma.h
                  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PN})

    install(TARGETS gg
            EXPORT "${PN}Targets"
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
                  ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
            DESTINATION ${CMAKECONFIG_INSTALL_DIR})
    install(EXPORT "${PN}Targets"
            NAMESPACE "${PN}::"
            DESTINATION ${CMAKECONFIG_INSTALL_DIR})
    export(EXPORT "${PN}Targets"
           NAMESPACE "${PN}::"
           FILE "${PROJECT_BINARY_DIR}/${PN}Targets.cmake")
endif()

if(${INSTALL_PYMOD})
    if(${NATIVE_PYTHON_INSTALL})
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                        "import sys; print(sys.prefix);"
                        OUTPUT_VARIABLE CMAKE_INSTALL_PREFIX
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                        "from distutils import sysconfig as s; import os; import sys; cmake_install_prefix = sys.prefix; prefix_lib = s.get_config_var('LIBDIR'); print(prefix_lib.replace(os.path.commonpath([prefix_lib, cmake_install_prefix]), '').strip('/'));"
                        OUTPUT_VARIABLE CMAKE_INSTALL_LIBDIR
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                        "from distutils import sysconfig as s; import os; prefix_lib = s.get_config_var('LIBDIR'); spdir = s.get_python_lib(plat_specific=True); print(spdir.replace(os.path.commonpath([prefix_lib, spdir]), ''));"
                        OUTPUT_VARIABLE PYMOD_INSTALL_LIBDIR
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()

    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
                    "from numpy import distutils; print(distutils.misc_util.get_shared_lib_extension(is_python_ext=False))"
                    OUTPUT_VARIABLE PYLIB_EXTENSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    install(DIRECTORY gau2grid
            DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}
            USE_SOURCE_PERMISSIONS
            FILES_MATCHING PATTERN "*.py")

    install(FILES $<TARGET_FILE:gg>
            DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/gau2grid
            RENAME "gg${PYLIB_EXTENSION}")

    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE
            DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/gau2grid)
endif()
