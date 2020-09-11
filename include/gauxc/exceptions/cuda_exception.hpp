#pragma once

#include <gauxc/gauxc_config.hpp>
#include <stdexcept>
#include <string>
#include <sstream>

#ifdef GAUXC_ENABLE_CUDA

namespace GauXC {

class cuda_exception : public std::exception {

  std::string file_;
  int         line_;
  std::string msg_prefix_;
  cudaError_t err_code_;

  const char* what() const noexcept override {
     std::stringstream ss;
     ss << "CUDA Exception (" << msg_prefix_ << ")" << std::endl
        << "  Error Code " << int(err_code_) << ": \"" 
                           << cudaGetErrorString( err_code_ ) << "\"" << std::endl
        << "  File       " << file_ << std::endl
        << "  Line       " << line_ << std::endl;

     auto msg = ss.str();

     return strdup( msg.c_str() );
  }

public:

  cuda_exception( std::string file, int line, std::string msg, cudaError_t err ) :
    file_(file), line_(line), msg_prefix_(msg), err_code_(err) { }

};

}

#define GAUXC_CUDA_ERROR( MSG, ERR ) \
  if( ERR != cudaSuccess ) \
    throw cuda_exception( __FILE__, __LINE__, MSG, ERR );

#endif
