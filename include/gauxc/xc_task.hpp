#pragma once

#include <array>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <gauxc/gauxc_config.hpp>
#include <gauxc/shell.hpp>
#include <gauxc/exceptions.hpp>

namespace GauXC {

struct XCTask {

  int32_t                              iParent = -1;
  std::vector< std::array<double,3> >  points;
  std::vector< double  >               weights;
  int32_t                              npts = 0;

  double                               dist_nearest;
  double                               max_weight = std::numeric_limits<double>::infinity();

  struct screening_data {
    std::vector<int32_t>               shell_list;
    std::vector<int32_t>               submat_block;
    std::vector<std::array<int32_t,3>> submat_map;
    int32_t                            nbe = 0;

    bool equiv_with( const screening_data& other ) const {
      return shell_list == other.shell_list;
    }
  };

  screening_data bfn_screening;
  screening_data cou_screening;

  void merge_with( const XCTask& other ) {
    if( !equiv_with(other) )
      GAUXC_GENERIC_EXCEPTION("Cannot Perform Requested Merge: Incompatible Tasks");
    points.insert( points.end(), other.points.begin(), other.points.end() );
    weights.insert( weights.end(), other.weights.begin(), other.weights.end() );
    npts = points.size();
  }

  template <typename TaskIt>
  void merge_with( TaskIt begin, TaskIt end ) {

    size_t old_sz = points.size();
    size_t pts_add = std::accumulate( begin, end, 0ul,
      []( const auto &a, const auto &t ) {
        return a + t.points.size();
      });

    size_t new_sz = old_sz + pts_add;
    points.resize( new_sz );
    weights.resize( new_sz );

    auto points_it  = points.begin()  + old_sz;
    auto weights_it = weights.begin() + old_sz;
    for( auto it = begin; it != end; ++it ) {
      if( !equiv_with(*it) )
        GAUXC_GENERIC_EXCEPTION("Cannot Perform Requested Task Merge");
      points_it  = std::copy( it->points.begin(), it->points.end(), points_it );
      weights_it = std::copy( it->weights.begin(), it->weights.end(), weights_it );
    }

    npts = points.size();
  }


  inline bool equiv_with( const XCTask& other ) const {
    return iParent == other.iParent and 
      bfn_screening.equiv_with(other.bfn_screening);
  }

  template <typename Archive>
  void serialize( Archive& ar ) {
    ar( iParent, bfn_screening.nbe, npts, dist_nearest, max_weight, 
      bfn_screening.shell_list, points, weights );  
  }


  inline size_t cost(size_t n_deriv, size_t natoms) const {
    return (bfn_screening.nbe * ( 1 + bfn_screening.nbe + n_deriv ) + natoms * natoms) * npts;
  }
};


}
