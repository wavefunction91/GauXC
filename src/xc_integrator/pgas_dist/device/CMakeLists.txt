#
# GauXC Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of
# any required approvals from the U.S. Dept. of Energy). All rights reserved.
#
# See LICENSE.txt for details
#
target_sources( gauxc PRIVATE 
  pgas_dist_xc_device_integrator.cxx
  shellbatched_pgas_dist_xc_device_integrator.cxx
)

