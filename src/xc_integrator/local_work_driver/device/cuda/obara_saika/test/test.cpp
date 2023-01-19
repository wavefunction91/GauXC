#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <libint2.hpp>
#include <gpu/integral_data_types.hpp>
#include <gpu/obara_saika_integrals.hpp>
#include <gpu/chebyshev_boys_computation.hpp>
#include <vector>
#include <iostream>
#include <algorithm>
#include <random>
#include <sys/time.h>

int main(int argc, char** argv) {
  libint2::initialize();

  // Benzene
  std::vector<libint2::Atom> atoms = {
    libint2::Atom{ 6,  6.92768e-01,  -1.77656e+00,   1.40218e-03},
    libint2::Atom{ 6,  3.35108e+00,  -1.77668e+00,   2.21098e-03},
    libint2::Atom{ 6,  4.68035e+00,   5.25219e-01,   1.22454e-03},
    libint2::Atom{ 6,  3.35121e+00,   2.82744e+00,  -7.02978e-04},
    libint2::Atom{ 6,  6.93087e-01,   2.82756e+00,  -1.55902e-03},
    libint2::Atom{ 6, -6.36278e-01,   5.25491e-01,  -4.68652e-04},
    libint2::Atom{ 1, -3.41271e-01,  -3.56759e+00,   2.21287e-03},
    libint2::Atom{ 1,  4.38492e+00,  -3.56783e+00,   3.73599e-03},
    libint2::Atom{ 1,  6.74844e+00,   5.25274e-01,   1.88028e-03},
    libint2::Atom{ 1,  4.38551e+00,   4.61832e+00,  -1.48721e-03},
    libint2::Atom{ 1, -3.41001e-01,   4.61857e+00,  -3.05569e-03},
    libint2::Atom{ 1, -2.70437e+00,   5.25727e-01,  -1.09793e-03} 
  };

  // Create cc-pVDZ BasisSet
  const std::string basis_name = "cc-pVDZ";
  libint2::BasisSet basis( basis_name, atoms );
  basis.set_pure(false); // Reset to Cartesian
  auto shell2bf = basis.shell2bf();

  auto [min_x, max_x] = std::minmax_element( atoms.begin(), atoms.end(), 
    []( const auto& a, const auto& b) { return a.x < b.x; } );
  auto [min_y, max_y] = std::minmax_element( atoms.begin(), atoms.end(), 
    []( const auto& a, const auto& b) { return a.y < b.y; } );
  auto [min_z, max_z] = std::minmax_element( atoms.begin(), atoms.end(), 
    []( const auto& a, const auto& b) { return a.z < b.z; } );

  std::array<double,3> box_lo = { min_x->x, min_y->y, min_z->z };
  std::array<double,3> box_hi = { max_x->x, max_y->y, max_z->z };

  std::default_random_engine gen;
  std::uniform_real_distribution<double> 
    dist_x( box_lo[0], box_hi[0] ),
    dist_y( box_lo[1], box_hi[1] ),
    dist_z( box_lo[2], box_hi[2] );

  auto gen_grid_point = [&]() {
    return std::array<double,3>{ dist_x(gen), dist_y(gen), dist_z(gen) };
  };

  if( argc != 2 ) throw std::runtime_error("Must Specify NGrid");
  
  const int ngrid = std::stoll( std::string(argv[1]) );
  
  std::vector< std::array<double,3> > grid_points( ngrid );
  std::generate( grid_points.begin(), grid_points.end(), gen_grid_point );

  const size_t nbf = basis.nbf();
  std::cout << "Running sn-LinK Proxy App with Settings:" << std::endl
	    << "  * NBF   = " << nbf << std::endl
	    << "  * NGRID = " << ngrid << std::endl
	    << std::endl;

  std::vector<libint2::Engine> engines;
  engines.reserve(ngrid);
  for( const auto& g : grid_points ) {
    engines.emplace_back( libint2::Operator::nuclear, basis.max_nprim(),
		          basis.max_l(), 0 );
    std::vector< std::pair<double, std::array<double,3>> > q = { {-1., g} }; 
    engines.back().set_params(q);
  }

  // Generate a random F matrix
  std::vector<double> F( ngrid * nbf );
  std::generate( F.begin(), F.end(), [&](){ return dist_x(gen); } );
  
  // Generate random grid weights
  std::vector<double> w( ngrid );
  std::generate( w.begin(), w.end(), [&](){ return dist_x(gen); } );

  // Compute A
  std::vector<double> A( nbf * nbf * ngrid );
  
  using row_major_mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using const_row_major_map = Eigen::Map< const row_major_mat >;
  
  using col_major_mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  using col_major_map = Eigen::Map< col_major_mat >;

  // correctness - libint implementation
  
  for( int k = 0; k < ngrid; ++k ) {
    auto& engine = engines.at(k);
    const auto& engine_buf = engine.results();

    col_major_map A_k( A.data() + nbf * nbf * k, nbf, nbf );

    for( int j = 0; j < basis.size(); ++j) {
      auto bf_j = shell2bf[j];
      auto nj   = basis[j].size();
      
      for( int i = 0; i < basis.size(); ++i) {
        auto bf_i = shell2bf[i];
        auto ni   = basis[i].size();

	engine.compute( basis[i], basis[j] );
	  
	const_row_major_map buf_map( engine_buf[0], ni, nj );
	A_k.block( bf_i, bf_j, ni, nj ) = buf_map;
      }
    }
  }

  std::vector<double> G_libint( ngrid * nbf );
  for(int i = 0; i < ngrid * nbf; ++i) {
    G_libint[i] = 0.0;
  }
  for( int k = 0; k < ngrid; ++k ) {
    for( int i = 0; i < nbf; ++i ) {
      double tmp = 0.0;

      for( int j = 0; j < nbf; ++j )
        tmp += A[i + j * nbf + k * nbf * nbf] * F[j * ngrid + k];

      G_libint[ i * ngrid + k ] = w[k] * tmp;
    }
  }

  // correctness - own implementation

  struct timeval start, end;

  int nshells = basis.size();
  
  double *dev_boys_table = XGPU::boys_init(); 

  std::vector<double> GPU_G_own( ngrid * nbf );
  for(int i = 0; i < ngrid * nbf; ++i) {
    GPU_G_own[i] = 0.0;
  }
  
  std::vector<XGPU::shells> _shells;
  std::vector<double> _points_transposed(3 * ngrid);

  _shells.resize(nshells);
  
  for( int i = 0; i < ngrid; ++i ){
    _points_transposed[i + 0 * ngrid] = grid_points[i][0];
    _points_transposed[i + 1 * ngrid] = grid_points[i][1];
    _points_transposed[i + 2 * ngrid] = grid_points[i][2];
      }

  for( int i = 0; i < nshells; ++i ) {
    _shells[i].origin.x = basis[i].O[0];
    _shells[i].origin.y = basis[i].O[1];
    _shells[i].origin.z = basis[i].O[2];

    _shells[i].m = basis[i].alpha.size();
    _shells[i].L = basis[i].contr[0].l;
    
    _shells[i].coeff = new XGPU::coefficients[_shells[i].m];
    for( int j = 0; j < _shells[i].m; ++j ) {
      _shells[i].coeff[j].alpha = basis[i].alpha[j];
      _shells[i].coeff[j].coeff = basis[i].contr[0].coeff[j];
    }
  }

  int total_prim_pairs = 0;
  for( int i = 0; i < nshells; ++i) {
    for( int j = 0; j <= i; ++j) {
      total_prim_pairs += (_shells[i].m * _shells[j].m);
    }
  }

  XGPU::prim_pair *prim_pairs = new XGPU::prim_pair[total_prim_pairs];

  int offset = 0;
  for( int i = 0; i < nshells; ++i) {
    for( int j = 0; j <= i; ++j) {
      if( _shells[i].L >= _shells[j].L )
	XGPU::generate_shell_pair(_shells[i], _shells[j], (prim_pairs + offset));
      else
	XGPU::generate_shell_pair(_shells[j], _shells[i], (prim_pairs + offset));

      offset += (_shells[i].m * _shells[j].m);
    }
  }

  double *dev_points_transposed, *dev_X, *dev_G, *dev_weights;
  XGPU::prim_pair *dev_prim_pairs;
  
  cudaMalloc((void**) &dev_points_transposed, 3 * ngrid * sizeof(double));
  cudaMalloc((void**) &dev_X, ngrid * nbf * sizeof(double));
  cudaMalloc((void**) &dev_G, ngrid * nbf * sizeof(double));
  cudaMalloc((void**) &dev_weights, ngrid * sizeof(double));

  cudaMalloc((void**) &dev_prim_pairs, total_prim_pairs * sizeof(XGPU::prim_pair));
  
  cudaMemcpy(dev_points_transposed, _points_transposed.data(), 3 * ngrid * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_X, F.data(), ngrid * nbf * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_G, GPU_G_own.data(), ngrid * nbf * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_weights, w.data(), ngrid * sizeof(double), cudaMemcpyHostToDevice);

  cudaMemcpy(dev_prim_pairs, prim_pairs, total_prim_pairs * sizeof(XGPU::prim_pair), cudaMemcpyHostToDevice);
    
  double *Xi = dev_X;
  double *Xj = dev_X;

  double *Gi = dev_G;
  double *Gj = dev_G;

  gettimeofday(&start, NULL);
  offset = 0;
  int ioff_cart = 0;
  for( int i = 0; i < nshells; ++i) {
    XGPU::shells bra_shell = _shells[i];
    int bra_cart_size = (bra_shell.L + 1) * (bra_shell.L + 2) / 2;

    int joff_cart = 0;
    for( int j = 0; j <= i; ++j) {
      XGPU::shells ket_shell = _shells[j];
      int ket_cart_size = (ket_shell.L + 1) * (ket_shell.L + 2) / 2;

      XGPU::compute_integral_shell_pair(i == j,
					ngrid,
					dev_points_transposed,
					_shells[i].L,
					_shells[j].L,
					_shells[i].origin,
					_shells[j].origin,
					(_shells[i].m * _shells[j].m),
					(dev_prim_pairs + offset),
					(Xi + ioff_cart * ngrid),
					(Xj + joff_cart * ngrid),
					ngrid,
					(Gi + ioff_cart * ngrid),
					(Gj + joff_cart * ngrid),
					ngrid,
					dev_weights,
					dev_boys_table);
      
      offset += (_shells[i].m * _shells[j].m);
      
      joff_cart += ket_cart_size;
    }

    ioff_cart += bra_cart_size;
  }

  cudaDeviceSynchronize();
  
  gettimeofday(&end, NULL);

  cudaMemcpy(GPU_G_own.data(), dev_G, ngrid * nbf * sizeof(double), cudaMemcpyDeviceToHost);

  cudaFree(dev_X);
  cudaFree(dev_G);
  cudaFree(dev_points_transposed);
  cudaFree(dev_weights);
  cudaFree(dev_prim_pairs);
  
  XGPU::boys_finalize(dev_boys_table);

  int correct = 1;
  
  for( int i = 0; i < nbf * ngrid; ++i) {	
    if((fabs(G_libint[i] - GPU_G_own[i]) > 1e-2) || std::isnan(GPU_G_own[i])) {
      printf("%lf %lf\n", G_libint[i], GPU_G_own[i]);
      correct = 0;
    }
  }

  std::cout << "Correctness: " << correct << "\tExecution: "<< 1000000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) << std::endl;

  cudaFree(dev_X);
  cudaFree(dev_G);
  cudaFree(dev_points_transposed);
  cudaFree(dev_weights);
  
  libint2::finalize();  // done with libint
}
