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

  int32_t                              iParent;
  std::vector< std::array<double,3> >  points;
  std::vector< double  >               weights;
  std::vector< int32_t >               shell_list;
  int32_t                              nbe;
  int32_t                              npts;

  double                               dist_nearest;
  double                               max_weight;

  std::vector< std::array<int32_t,3> > submat_map;
  std::vector< int32_t >               submat_block;

  void merge_with( const XCTask& other ) {
    if( shell_list != other.shell_list or iParent != other.iParent )
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
      if( shell_list != it->shell_list or iParent != it->iParent )
        GAUXC_GENERIC_EXCEPTION("Cannot Perform Requested Task Merge");
      points_it  = std::copy( it->points.begin(), it->points.end(), points_it );
      weights_it = std::copy( it->weights.begin(), it->weights.end(), weights_it );
    }

    npts = points.size();
  }


  inline bool equiv_with( const XCTask& other ) const {
    return iParent == other.iParent and 
           shell_list == other.shell_list;
  }

  template <typename Archive>
  void serialize( Archive& ar ) {
    ar( iParent, nbe, npts, dist_nearest, max_weight, 
      shell_list, points, weights );  
  }


  inline size_t cost(size_t n_deriv, size_t natoms) const {
    return (nbe * ( 1 + nbe + n_deriv ) + natoms * natoms) * npts;
  }
};


}
