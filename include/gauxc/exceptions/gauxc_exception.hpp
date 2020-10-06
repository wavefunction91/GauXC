#pragma once

#include <stdexcept>
#include <string>
#include <string.h>
#include <sstream>

namespace GauXC {

class exchcxx_exception : public std::exception {

  std::string file_;
  int         line_;
  std::string msg_prefix_;

  const char* what() const noexcept override {
     std::stringstream ss;
     ss << "GAUXC Exception (" << msg_prefix_ << ")" << std::endl
        << "  File       " << file_ << std::endl
        << "  Line       " << line_ << std::endl;

     auto msg = ss.str();

     return strdup( msg.c_str() );
  }

public:

  exchcxx_exception( std::string file, int line, std::string msg) :
    file_(file), line_(line), msg_prefix_(msg) { }

};

#define GAUXC_BOOL_CHECK( MSG, V ) \
  if( not (V) ) \
    throw exchcxx_exception( __FILE__, __LINE__, MSG );

}

