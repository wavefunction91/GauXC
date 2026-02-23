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

#include <gauxc/functional.h>
#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/enums/spin.hpp>
#include <exchcxx/enums/backend.hpp>

#include "c_functional.hpp"
#include "c_status.hpp"

namespace GauXC::detail {
/**
 *  Splits a string into tokens  based on a delimiter
 *
 *  \param [out] tokens     std::vector of std::string objects which hold
 *                          the split tokens
 *  \param [in]  str        std::string to split
 *  \param [in]  delimiters Delimiters on which to split str
 */
static inline void split(std::vector<std::string>& tokens, 
  const std::string& str, const std::string& delimiters = " ") {

    tokens.clear();
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}; // split
} // namespace GauXC::detail

namespace GauXC::C {
extern "C" {

GauXCFunctional gauxc_functional_from_string(
  GauXCStatus* status,
  const char* functional_spec,
  bool polarized
) {
  detail::gauxc_status_init(status);
  GauXCFunctional functional{};
  functional.hdr = GauXCHeader{GauXC_Type_Functional};
  functional.ptr = nullptr;

  try {
    auto polar = polarized ? ExchCXX::Spin::Polarized : ExchCXX::Spin::Unpolarized;
    functional_type func;
    if(ExchCXX::functional_map.key_exists(functional_spec)) {
      func = functional_type( ExchCXX::Backend::builtin, ExchCXX::functional_map.value(functional_spec), 
        polar );
    } 
#ifdef EXCHCXX_ENABLE_LIBXC
    else { 
      std::vector<std::pair<double, ExchCXX::XCKernel>> funcs;
      std::vector<std::string> libxc_names;
      detail::split(libxc_names, functional_spec, ",");
      for( auto n : libxc_names ) {
        funcs.push_back( {1.0, ExchCXX::XCKernel(ExchCXX::libxc_name_string(n), polar)} );
      }
      func = functional_type(funcs);
    }
#endif

    functional.ptr = new functional_type( std::move(func) );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return functional;
}

GauXCFunctional gauxc_functional_from_enum(
  GauXCStatus* status,
  enum GauXC_Functional functional_enum,
  bool polarized
) {
  detail::gauxc_status_init(status);
  GauXCFunctional functional{};
  functional.hdr = GauXCHeader{GauXC_Type_Functional};
  functional.ptr = nullptr;

  try {
    auto polar = polarized ? ExchCXX::Spin::Polarized : ExchCXX::Spin::Unpolarized;
    functional_type func( ExchCXX::Backend::builtin, 
      static_cast<ExchCXX::Functional>(functional_enum), polar );

    functional.ptr = new functional_type( std::move(func) );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return functional;
}

void gauxc_functional_delete(
  GauXCStatus* status,
  GauXCFunctional* functional
) {
  detail::gauxc_status_init(status);
  if(functional == nullptr) return;
  if(functional->ptr != nullptr)
    delete detail::get_functional_ptr(*functional);
  functional->ptr = nullptr;
}

} // extern "C"
} // namespace GauXC::C