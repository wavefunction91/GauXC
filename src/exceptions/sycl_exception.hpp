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

namespace GauXC {

/**
 *  @brief A class to handle excecptions arising from SYCL operations
 */
class sycl_exception : public std::exception {

  std::string file_;       ///< File which contains the code that threw the exception
  int         line_;       ///< Line number of file_ that threw exception
  std::string msg_prefix_; ///< General descriptor of task which threw exception
  std::string err_str_;    ///< SYCL error string pertaining to the thrown exception
  std::string what_msg_;

  /**
   *  @brief Get a descriptive message pertaining to the thrown SYCL error
   *
   *  @returns a descritive message pertaining to the SYCL error represented by
   *  the internal state of the exception object.
   */
  const char* what() const noexcept override {
     return what_msg_.c_str();
  }

public:

  /**
   *  @brief Construct a sycl_exception object
   *
   *  @param[in] file File which contains the code that threw the exception
   *  @param[in] line Line number of file that threw exception
   *  @param[in] msg  General descriptor of task which threw exception
   *  @param[in] err  SYCL exception pertaining to the thrown exception
   */
  sycl_exception( std::string file, int line, std::string msg,
    const ::sycl::exception& err ) :
    file_(file), line_(line), msg_prefix_(msg), err_str_(err.what()) {
    std::stringstream ss;
    ss << "SYCL Exception (" << msg_prefix_ << ")" << std::endl
       << "  Error      \"" << err_str_ << "\"" << std::endl
       << "  File       " << file_ << std::endl
       << "  Line       " << line_ << std::endl;
    what_msg_ = ss.str();
  }

}; // class sycl_exception

} // namespace GauXC

// Macro to wrap SYCL error handling
#define GAUXC_SYCL_ERROR( MSG, EXPR ) \
  try { EXPR; } \
  catch( const ::sycl::exception& e ) { \
    throw sycl_exception( __FILE__, __LINE__, MSG, e ); \
  }

#endif
