#
# GauXC Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of
# any required approvals from the U.S. Dept. of Energy).
#
# (c) 2024-2025, Microsoft Corporation
#
# All rights reserved.
#
# See LICENSE.txt for details
#

target_sources( gauxc PRIVATE
  status.f90
  types.f90
  enums.f90
  atom.f90
  molecule.f90
  shell.f90
  basisset.f90
  molgrid.f90
  runtime_environment.F90
  load_balancer.f90
  molecular_weights.f90
  functional.f90
  xc_integrator.f90
)
set_target_properties( gauxc PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules
)
if(NOT EXISTS ${PROJECT_BINARY_DIR}/modules)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
endif()
target_include_directories(
  gauxc
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/gauxc/modules>
)
if(GAUXC_ENABLE_MPI)
  target_link_libraries( gauxc PUBLIC MPI::MPI_Fortran )
endif()

# Generate config file
configure_file( 
  ${PROJECT_SOURCE_DIR}/include/gauxc/gauxc_config.f.in
  ${PROJECT_BINARY_DIR}/include/gauxc/gauxc_config.f
)

# Install Fortran modules
install(
  DIRECTORY ${PROJECT_BINARY_DIR}/modules/
  DESTINATION include/gauxc/modules
  FILES_MATCHING PATTERN "*.mod"
)

install(
  FILES       ${PROJECT_BINARY_DIR}/include/gauxc/gauxc_config.f
  DESTINATION include/gauxc
)
