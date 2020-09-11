#pragma once
#include "blas.hpp"
#include <vector>
#include <tuple>
#include <cstdint>

namespace GauXC  {
namespace detail {

template <typename _F1, typename _F2>
void submat_set(int32_t M, int32_t N, int32_t MSub, 
  int32_t NSub, _F1 *ABig, int32_t LDAB, _F2 *ASmall, 
  int32_t LDAS, 
  std::vector<std::pair<int32_t,int32_t>> &submat_map) {

  (void)(M);
  (void)(N);
  (void)(MSub);
  (void)(NSub);

  int32_t i(0);
  for( auto& iCut : submat_map ) {
    int32_t deltaI = iCut.second - iCut.first;
    int32_t j(0);
  for( auto& jCut : submat_map ) {
    int32_t deltaJ = jCut.second - jCut.first;
  
    auto* ABig_use   = ABig   + iCut.first + jCut.first * LDAB;
    auto* ASmall_use = ASmall + i          + j          * LDAS;


    GauXC::blas::lacpy( 'A', deltaI, deltaJ, ABig_use, LDAB, 
                         ASmall_use, LDAS );

  
    j += deltaJ;
  }
    i += deltaI;
  }
  

}

template <typename _F1, typename _F2>
void inc_by_submat(int32_t M, int32_t N, int32_t MSub, 
  int32_t NSub, _F1 *ABig, int32_t LDAB, _F2 *ASmall, 
  int32_t LDAS, 
  std::vector<std::pair<int32_t,int32_t>> &submat_map) {

  (void)(M);
  (void)(N);
  (void)(MSub);
  (void)(NSub);

  int32_t i(0);
  for( auto& iCut : submat_map ) {
    int32_t deltaI = iCut.second - iCut.first;
    int32_t j(0);
  for( auto& jCut : submat_map ) {
    int32_t deltaJ = jCut.second - jCut.first;
  
    auto* ABig_use   = ABig   + iCut.first + jCut.first * LDAB;
    auto* ASmall_use = ASmall + i          + j          * LDAS;


    for( int32_t jj = 0; jj < deltaJ; ++jj )
    for( int32_t ii = 0; ii < deltaI; ++ii )
      ABig_use[ ii + jj * LDAB ] += ASmall_use[ ii + jj * LDAS ];

  
    j += deltaJ;
  }
    i += deltaI;
  }
  

}

}
}
