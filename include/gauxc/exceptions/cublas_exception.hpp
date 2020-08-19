#pragma once

#include <gauxc/gauxc_config.hpp>
#include <stdexcept>
#include <string>
#include <sstream>

#ifdef GAUXC_ENABLE_CUDA

namespace GauXC {

namespace detail {

static std::string cublasGetErrorString(cublasStatus_t error) {
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";

        case CUBLAS_STATUS_NOT_SUPPORTED:
            return "CUBLAS_STATUS_NOT_SUPPORTED";

        case CUBLAS_STATUS_LICENSE_ERROR:
            return "CUBLAS_STATUS_LICENSE_ERROR";
    }

    return "<unknown>";
}

}

class cublas_exception : public std::exception {

  std::string file_;
  int         line_;
  std::string msg_prefix_;
  cublasStatus_t err_code_;

  const char* what() const noexcept override {
     std::stringstream ss;
     ss << "CUBLAS Exception (" << msg_prefix_ << ")" << std::endl
        << "  Error Code " << int(err_code_) << ": \"" 
                           << detail::cublasGetErrorString( err_code_ ) 
                           << "\"" << std::endl
        << "  File       " << file_ << std::endl
        << "  Line       " << line_ << std::endl;

     auto msg = ss.str();

     return strdup( msg.c_str() );
  }

public:

  cublas_exception( std::string file, int line, std::string msg, 
                    cublasStatus_t err ) :
    file_(file), line_(line), msg_prefix_(msg), err_code_(err) { }

};

}


#define GAUXC_CUBLAS_ERROR( MSG, ERR ) \
  if( ERR != CUBLAS_STATUS_SUCCESS ) \
    throw cublas_exception( __FILE__, __LINE__, MSG, ERR );

#endif
