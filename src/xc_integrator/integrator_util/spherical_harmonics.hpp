#include <cmath>
#include <vector>
#include <array>
#include <cassert>


void scaled_ylm_matrix(const int lmax, const double* points, const int32_t  npts, const std::array<double, 3> center, const double radius, double* ylm_matrix);