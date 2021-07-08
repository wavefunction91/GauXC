#pragma once

#include <gauxc/gauxc_config.hpp>
#include <stdexcept>
#include <string>
#include <sstream>

#ifdef GAUXC_ENABLE_HIP

namespace GauXC {

class hip_exception : public std::exception {

  std::string file_;
  int         line_;
  std::string msg_prefix_;
  hipError_t err_code_;

  const char* what() const noexcept override {
     std::stringstream ss;
     ss << "HIP Exception (" << msg_prefix_ << ")" << std::endl
        << "  Error Code " << int(err_code_) << ": \"" 
                           << hipGetErrorString( err_code_ ) << "\"" << std::endl
        << "  File       " << file_ << std::endl
        << "  Line       " << line_ << std::endl;

     auto msg = ss.str();

     return strdup( msg.c_str() );
  }

public:

  hip_exception( std::string file, int line, std::string msg, hipError_t err ) :
    file_(file), line_(line), msg_prefix_(msg), err_code_(err) { }

};

}

#define GAUXC_HIP_ERROR( MSG, ERR ) \
  if( ERR != hipSuccess ) \
    throw hip_exception( __FILE__, __LINE__, MSG, ERR );

#endif
