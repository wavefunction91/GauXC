#pragma once
#include <memory>
#include <vector>

namespace GauXC {

template <uint32_t M, uint32_t NCheb, uint32_t MaxT, uint32_t NSegment, uint32_t LDTable>
void boys_chebyshev( int npts, const double* T, const double* boys_table, double* eval ); 

template <uint32_t NCheb, uint32_t MaxM, uint32_t MaxT, 
  uint32_t NSegment = (MaxT*NCheb)/2, 
  uint32_t LDTable = NCheb+1>
class ChebyshevBoysEvaluator {

  std::vector<double> table_;

public:

  ChebyshevBoysEvaluator();
  template <uint32_t M>
  void eval( size_t npts, const double* T, double* FmT ) {
    boys_chebyshev<M,NCheb,MaxT,NSegment,LDTable>( npts, T, table_.data(), FmT );
  }

};

void gauxc_boys_init();
void gauxc_boys_finalize();

namespace detail {
  constexpr uint32_t default_ncheb = 7;
  constexpr uint32_t default_max_m = 16;
  constexpr uint32_t default_max_t = 30;
  using default_chebyshev_type = ChebyshevBoysEvaluator< default_ncheb, default_max_m, default_max_t >;
}
extern std::unique_ptr<detail::default_chebyshev_type> chebyshev_boys_instance;
extern template class ChebyshevBoysEvaluator< detail::default_ncheb, detail::default_max_m, detail::default_max_t >;

}
