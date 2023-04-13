/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#ifdef GAUXC_ENABLE_MPI
  #define GAUXC_MPI_CODE(...) __VA_ARGS__
#else
  #define GAUXC_MPI_CODE(...) 
#endif

