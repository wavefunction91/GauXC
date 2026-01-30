/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/gauxc_config.h>

#ifdef GAUXC_HAS_MPI
  #define GAUXC_MPI_CODE(...) __VA_ARGS__
#else
  #define GAUXC_MPI_CODE(...) 
#endif

#ifdef GAUXC_HAS_MPI
#include <mpi.h>
#endif