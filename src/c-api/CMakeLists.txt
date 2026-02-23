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
  c_status.cxx
  c_types.cxx
  c_molecule.cxx
  c_basisset.cxx
  c_molgrid.cxx
  c_runtime_environment.cxx
  c_load_balancer.cxx
  c_molecular_weights.cxx
  c_functional.cxx
  c_xc_integrator.cxx
)
target_include_directories( gauxc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})