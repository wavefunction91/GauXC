#pragma once

#ifdef GAUXC_ENABLE_MPI
  #define GAUXC_MPI_CODE(...) __VA_ARGS__
#else
  #define GAUXC_MPI_CODE(...) 
#endif

