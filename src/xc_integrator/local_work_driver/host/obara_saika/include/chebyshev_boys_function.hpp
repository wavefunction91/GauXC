#pragma once
#include <memory>
#include <vector>

namespace GauXC {

class ChebyshevBoysEvaluator {

  int ncheb_;
  int maxM_;
  int nseg_;
  int ldtable_;

  double min_t_thresh_;
  double max_t_thresh_;

  std::vector<double> table_;

public:

  ChebyshevBoysEvaluator( int ncheb, int maxM, int nseg, double minT, double maxT );
  void eval( size_t npts, int m, const double* T, double* FmT );

};

void gauxc_boys_init(int ncheb, int maxM, int nseg, double minT, double maxT);
void gauxc_boys_finalize();
extern std::unique_ptr<ChebyshevBoysEvaluator> chebyshev_boys_instance;

}
