#pragma once

#include <iostream>

#define DEFAULT_NCHEB  7
#define DEFAULT_MAX_M  8
#define DEFAULT_MAX_T 30

#define DEFAULT_NSEGMENT ((DEFAULT_MAX_T * DEFAULT_NCHEB) / 2)
#define DEFAULT_LD_TABLE (DEFAULT_NCHEB + 1)

namespace XGPU {
  // create tables
  double *boys_init();
  void boys_finalize(double *boys_table);
}

