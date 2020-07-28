#pragma once

#include <vector>
#include <numeric>

#include <gauxc/shell.hpp>

namespace GauXC {

template <typename F>
class BasisSet : public std::vector<Shell<F>> {

  using ao_range = std::pair<int32_t,int32_t>;

  std::vector<int32_t>   shell_to_first_ao_;
  std::vector<ao_range> shell_to_ao_range_;


public:

  template <typename... Args>
  BasisSet( Args&&... args ) :
    std::vector<Shell<F>>( std::forward<Args>(args)... )  { }

  BasisSet( const BasisSet& ) = default;
  BasisSet( BasisSet&& ) noexcept = default;

  inline int32_t nshells() const { return this->size(); }; 
  inline int32_t nbf()     const {
    return std::accumulate( this->cbegin(), this->cend(), 0ul,
      [](const auto& a, const auto& b) { 
        return a + b.size();
      } );
  };

  void generate_shell_to_ao() {

    size_t st_idx = 0;
    for( const auto& shell : (*this) ) {
      size_t range_end = st_idx + shell.size();
      shell_to_first_ao_.emplace_back( st_idx );
      shell_to_ao_range_.push_back({ st_idx, range_end });
      st_idx = range_end;
    }

  }

  const auto& shell_to_first_ao() const { return shell_to_first_ao_; }
        auto& shell_to_first_ao()       { return shell_to_first_ao_; }
  const auto& shell_to_ao_range() const { return shell_to_ao_range_; }
        auto& shell_to_ao_range()       { return shell_to_ao_range_; }

  auto shell_to_first_ao(int32_t i) const { return shell_to_first_ao_.at(i); }
  auto shell_to_ao_range(int32_t i) const { return shell_to_ao_range_.at(i); }

};

}
