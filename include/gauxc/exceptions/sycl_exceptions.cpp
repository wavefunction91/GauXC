#pragma once

#include <gauxc/gauxc_config.hpp>
#include <CL/sycl.hpp>
#include <stdexcept>
#include <sstream>

#ifdef GAUXC_ENABLE_SYCL

#define GAUXC_SYCL_ERROR( EXPR )                                 \
  try {                                                          \
    EXPR;                                                        \
  }                                                              \
  catch (cl::sycl::exception const &ex) {                        \
    std::stringstream msg;                                       \
    msg << "SYCL Exception at " << __FILE__ << " : " << __LINE__ \
        << std::endl;                                            \
    throw(std::runtime_error( ex.what() ));                      \
  }

#endif // GAUXC_ENABLE_SYCL
