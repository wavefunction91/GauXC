namespace GauXC {

void cuda_aos_scheme1_weights_wrapper( int32_t npts, int32_t natoms,
  const double* points, const double* RAB, int32_t ldRAB, const double* coords, 
  double* dist, int32_t lddist, const int32_t* iparent,
  const double* dist_nearest, double* weights, cudaStream_t stream );

}
