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

#include <gauxc/exceptions.hpp>
#include <stdexcept>
#include <string>
#include <sstream>

#ifdef GAUXC_HAS_SYCL
#include <sycl/sycl.hpp>
// Only the exception type is needed here. The umbrella <oneapi/mkl.hpp> pulls
// in the C LAPACK prototypes, whose parameter named "csl" collides with the
// csl macro from buffer_adaptor.hpp
#include <oneapi/mkl/exceptions.hpp>

namespace GauXC {

/**
 *  @brief A class to handle excecptions arising from oneMKL operations
 */
class onemkl_exception : public std::exception {

  std::string file_;       ///< File which contains the code that threw the exception
  int         line_;       ///< Line number of file_ that threw exception
  std::string msg_prefix_; ///< General descriptor of task which threw exception
  std::string err_str_;    ///< oneMKL error string pertaining to the thrown exception
  std::string what_msg_;

  /**
   *  @brief Get a descriptive message pertaining to the thrown oneMKL error
   *
   *  @returns a descritive message pertaining to the oneMKL error represented
   *  by the internal state of the exception object.
   */
  const char* what() const noexcept override {
     return what_msg_.c_str();
  }

public:

  /**
   *  @brief Construct a onemkl_exception object
   *
   *  @param[in] file File which contains the code that threw the exception
   *  @param[in] line Line number of file that threw exception
   *  @param[in] msg  General descriptor of task which threw exception
   *  @param[in] err  oneMKL exception pertaining to the thrown exception
   */
  onemkl_exception( std::string file, int line, std::string msg,
    const std::exception& err ) :
    file_(file), line_(line), msg_prefix_(msg), err_str_(err.what()) {
    std::stringstream ss;
    ss << "oneMKL Exception (" << msg_prefix_ << ")" << std::endl
       << "  Error      \"" << err_str_ << "\"" << std::endl
       << "  File       " << file_ << std::endl
       << "  Line       " << line_ << std::endl;
    what_msg_ = ss.str();
  }

}; // class onemkl_exception

} // namespace GauXC

// Macro to wrap oneMKL error handling
#define GAUXC_ONEMKL_ERROR( MSG, EXPR ) \
  try { EXPR; } \
  catch( const oneapi::mkl::exception& e ) { \
    throw onemkl_exception( __FILE__, __LINE__, MSG, e ); \
  }

#endif
