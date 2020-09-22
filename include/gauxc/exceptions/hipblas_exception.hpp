#pragma once

#include <gauxc/gauxc_config.hpp>
#include <stdexcept>
#include <string>
#include <sstream>

#ifdef GAUXC_ENABLE_HIP

namespace GauXC {

namespace detail {

static std::string hipblasGetErrorString(hipblasStatus_t error) {
    switch (error)
    {
        case HIPBLAS_STATUS_SUCCESS:
            return "HIPBLAS_STATUS_SUCCESS";

        case HIPBLAS_STATUS_NOT_INITIALIZED:
            return "HIPBLAS_STATUS_NOT_INITIALIZED";

        case HIPBLAS_STATUS_ALLOC_FAILED:
            return "HIPBLAS_STATUS_ALLOC_FAILED";

        case HIPBLAS_STATUS_INVALID_VALUE:
            return "HIPBLAS_STATUS_INVALID_VALUE";

        case HIPBLAS_STATUS_ARCH_MISMATCH:
            return "HIPBLAS_STATUS_ARCH_MISMATCH";

        case HIPBLAS_STATUS_MAPPING_ERROR:
            return "HIPBLAS_STATUS_MAPPING_ERROR";

        case HIPBLAS_STATUS_EXECUTION_FAILED:
            return "HIPBLAS_STATUS_EXECUTION_FAILED";

        case HIPBLAS_STATUS_INTERNAL_ERROR:
            return "HIPBLAS_STATUS_INTERNAL_ERROR";

        case HIPBLAS_STATUS_NOT_SUPPORTED:
            return "HIPBLAS_STATUS_NOT_SUPPORTED";

    }

    return "<unknown>";
}

}

class hipblas_exception : public std::exception {

  std::string file_;
  int         line_;
  std::string msg_prefix_;
  hipblasStatus_t err_code_;

  const char* what() const noexcept override {
     std::stringstream ss;
     ss << "HIPBLAS Exception (" << msg_prefix_ << ")" << std::endl
        << "  Error Code " << int(err_code_) << ": \"" 
                           << detail::hipblasGetErrorString( err_code_ ) 
                           << "\"" << std::endl
        << "  File       " << file_ << std::endl
        << "  Line       " << line_ << std::endl;

     auto msg = ss.str();

     return strdup( msg.c_str() );
  }

public:

  hipblas_exception( std::string file, int line, std::string msg, 
                    hipblasStatus_t err ) :
    file_(file), line_(line), msg_prefix_(msg), err_code_(err) { }

};

}


#define GAUXC_HIPBLAS_ERROR( MSG, ERR ) \
  if( ERR != HIPBLAS_STATUS_SUCCESS ) \
    throw hipblas_exception( __FILE__, __LINE__, MSG, ERR );

#endif
