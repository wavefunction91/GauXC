#pragma once

#include <vector>
#include <numeric>

#include <gauxc/shell.hpp>

namespace GauXC {

/**
 *  @brief A class to manage a Gaussian type orbital (GTO) basis set
 *
 *  Extends std::vector<Shell<F>>
 *
 *  @tparam F Datatype representing the internal basis set storage
 */
template <typename F>
class BasisSet : public std::vector<Shell<F>> {

  // TODO: Move these to BasisSetMap
  using ao_range = std::pair<int32_t,int32_t>;

  std::vector<int32_t>   shell_to_first_ao_; ///< Map from shell index to first basis function of that shell
  std::vector<ao_range> shell_to_ao_range_; ///< Map from shell index to range of basis functions for that shell

public:

  /**
   *  @brief Construct a BasisSet object
   *
   *  Delegates to std::vector<Shell<F>>::vector
   *
   *  @tparam Args Parameter pack for arguements that are passed to
   *  base constructor
   */
  template <typename... Args>
  BasisSet( Args&&... args ) :
    std::vector<Shell<F>>( std::forward<Args>(args)... )  { }

  /// Copy a BasisSet object
  BasisSet( const BasisSet& )     = default;

  /// Move a BasisSet object
  BasisSet( BasisSet&& ) noexcept = default;

  /// Copy-assign BasisSet object
  BasisSet& operator=( const BasisSet& ) = default;

  /// Move-assign BasisSet object
  BasisSet& operator=( BasisSet&& ) noexcept = default;

  /**
   *  @brief Return the number of GTO shells which comprise the BasisSet object
   *
   *  Delegates to std::vector<Shell<F>>::size
   *
   *  @returns the number of GTO shells which comprise the BasisSet object
   */
  inline int32_t nshells() const { return this->size(); }; 

  /**
   *  @brief Return the number of GTO basis functions which comprise the 
   *  BasisSet object.
   *
   *  This routine accumulates the shell sizes (accounting for Cart/Sph angular
   *  factors) for each shell in the basis set.
   *
   *  @returns the number of GTO basis functions which comprise the BasisSet
   *  object.
   */
  inline int32_t nbf()     const {
    return std::accumulate( this->cbegin(), this->cend(), 0ul,
      [](const auto& a, const auto& b) { 
        return a + b.size();
      } );
  };

  /**
   *  @brief Determine the number of basis functions contained in a
   *  specified subset of the BasisSet object.
   *
   *  Performs the following operation:
   *    for( i in shell_list ) nbf += size of shell i
   *
   *  @tparam IntegralIterator Iterator type representing the list of
   *  shell indices.
   *
   *  @param[in] shell_list_begin Start iterator for shell list
   *  @param[in] shell_list_end   End iterator for shell_list
   *  @returns   Number of basis functions in the specified shell subset.
   */
  template <typename IntegralIterator>
  inline int32_t nbf_subset( IntegralIterator shell_list_begin,
                             IntegralIterator shell_list_end ) const {
    int32_t _nbf = 0;
    for( auto it = shell_list_begin; it != shell_list_end; ++it )
      _nbf += std::vector<Shell<F>>::at(*it).size();
    return _nbf;
  }

  /**
   *  @brief Generate the maps from shell indices to basis function indices
   *  TODO: Move this to BasisSetMap
   */
  void generate_shell_to_ao() {

    size_t st_idx = 0;
    for( const auto& shell : (*this) ) {
      size_t range_end = st_idx + shell.size();
      shell_to_first_ao_.emplace_back( st_idx );
      shell_to_ao_range_.push_back({ st_idx, range_end });
      st_idx = range_end;
    }

  }


  // TODO: Move these to BasisSetMap

  /**
   *  @brief Return the map from shell indicies to starting basis functions
   *
   *  Const variant
   *
   *  @returns The map from shell indices to starting basis functions
   */
  const auto& shell_to_first_ao() const { return shell_to_first_ao_; }

  /**
   *  @brief Return the map from shell indicies to starting basis functions
   *
   *  Non-const variant
   *
   *  @returns The map from shell indices to starting basis functions
   */
        auto& shell_to_first_ao()       { return shell_to_first_ao_; }

  /**
   *  @brief Return the map from shell indicies to basis function ranges
   *
   *  Const variant
   *
   *  @returns The map from shell indices to basis function ranges
   */
  const auto& shell_to_ao_range() const { return shell_to_ao_range_; }

  /**
   *  @brief Return the map from shell indicies to basis function ranges
   *
   *  Non-const variant
   *
   *  @returns The map from shell indices to basis function ranges
   */
        auto& shell_to_ao_range()       { return shell_to_ao_range_; }



  /**
   *  @brief Get first basis function index for a specified shell
   *
   *  @param[in] i Shell index
   *  @returns   Basis function index for shell "i"
   */
  auto shell_to_first_ao(int32_t i) const { return shell_to_first_ao_.at(i); }

  /**
   *  @brief Get first basis function range for a specified shell
   *
   *  @param[in] i Shell index
   *  @returns   Basis function range for shell "i"
   */
  auto shell_to_ao_range(int32_t i) const { return shell_to_ao_range_.at(i); }

}; // class BasisSet

} // namespace GauXC
