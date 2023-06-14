#if 0
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "gpu/integral_data_types.hpp"
#include "gpu/obara_saika_integrals.hpp"
#include "gpu/chebyshev_boys_computation.hpp"

#include "cpu/integral_data_types.hpp"
#include "cpu/obara_saika_integrals.hpp"
#include "cpu/chebyshev_boys_computation.hpp"

#include <vector>
#include <iostream>
#include <algorithm>
#include <random>
#include <sys/time.h>

#include <gauxc/molecule.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/external/hdf5.hpp>
#include <highfive/H5File.hpp>

int main(int argc, char* argv[]) {

  if( argc < 3 ) throw std::runtime_error("NOT VALID INPUT");

  struct timeval cpu_start, cpu_end, gpu_start, gpu_end;
  
  std::string data_file = argv[1];

  const int ngrid = std::stoll( std::string(argv[2]) );
  
  GauXC::Molecule mol;
  GauXC::read_hdf5_record(mol, data_file, "/MOLECULE");

  GauXC::BasisSet<double> basis;
  GauXC::read_hdf5_record(basis, data_file, "/BASIS");
  for( auto& sh : basis ) sh.set_pure(false); // Reset to cartesian

  std::cout << mol.size() << std::endl; 
  std::cout << basis.size() << std::endl;
  std::cout << basis.nbf() << std::endl;
  
  auto [min_x, max_x] = std::minmax_element( mol.begin(), mol.end(), []( const auto& a, const auto& b) { return a.x < b.x; } );
  auto [min_y, max_y] = std::minmax_element( mol.begin(), mol.end(), []( const auto& a, const auto& b) { return a.y < b.y; } );
  auto [min_z, max_z] = std::minmax_element( mol.begin(), mol.end(), []( const auto& a, const auto& b) { return a.z < b.z; } );

  std::array<double,3> box_lo = { min_x->x, min_y->y, min_z->z };
  std::array<double,3> box_hi = { max_x->x, max_y->y, max_z->z };

  std::default_random_engine gen;
  std::uniform_real_distribution<double>
    dist_x( box_lo[0], box_hi[0] ),
    dist_y( box_lo[1], box_hi[1] ),
    dist_z( box_lo[2], box_hi[2] );

  auto gen_grid_point = [&]() { return std::array<double,3>{ dist_x(gen), dist_y(gen), dist_z(gen) }; };

  std::vector< std::array<double,3> > grid_points( ngrid );
  std::generate( grid_points.begin(), grid_points.end(), gen_grid_point );

  int nbf = basis.nbf();
  std::vector<double> F( ngrid * nbf );
  std::generate( F.begin(), F.end(), [&](){ return dist_x(gen); } );

  std::vector<double> w( ngrid );
  std::generate( w.begin(), w.end(), [&](){ return dist_x(gen); } );

  int nshells = basis.size();
  
  std::vector<double> CPU_G_own( ngrid * nbf );
  for(int i = 0; i < ngrid * nbf; ++i) {
    CPU_G_own[i] = 0.0;
  }

  {
    double *cpu_boys_table = XCPU::boys_init();
    
    std::vector<XCPU::shells> _shells;
    std::vector<double> _points_transposed(3 * ngrid);

    _shells.resize(nshells);
  
    for( int i = 0; i < ngrid; ++i ){
      _points_transposed[i + 0 * ngrid] = grid_points[i][0];
      _points_transposed[i + 1 * ngrid] = grid_points[i][1];
      _points_transposed[i + 2 * ngrid] = grid_points[i][2];
    }

    for( int i = 0; i < nshells; ++i ) {
      _shells[i].origin.x = basis[i].O()[0];
      _shells[i].origin.y = basis[i].O()[1];
      _shells[i].origin.z = basis[i].O()[2];

      _shells[i].m = basis[i].nprim();
      _shells[i].L = basis[i].l();

      _shells[i].coeff = new XCPU::coefficients[_shells[i].m];
      for( int j = 0; j < _shells[i].m; ++j ) {
	_shells[i].coeff[j].alpha = basis[i].alpha()[j];
	_shells[i].coeff[j].coeff = basis[i].coeff()[j];
      }
    }

    int total_prim_pairs = 0;
    for( int i = 0; i < nshells; ++i) {
      for( int j = 0; j <= i; ++j) {
	total_prim_pairs += (_shells[i].m * _shells[j].m);
      }
    }

    XCPU::prim_pair *prim_pairs = new XCPU::prim_pair[total_prim_pairs];

    int offset = 0;
    for( int i = 0; i < nshells; ++i) {
      for( int j = 0; j <= i; ++j) {
	if( _shells[i].L >= _shells[j].L )
	  XCPU::generate_shell_pair(_shells[i], _shells[j], (prim_pairs + offset));
	else
	  XCPU::generate_shell_pair(_shells[j], _shells[i], (prim_pairs + offset));

	offset += (_shells[i].m * _shells[j].m);
      }
    }
  
    // CPU implementation
    double *Xi = F.data();
    double *Xj = F.data();

    double *Gi = CPU_G_own.data();
    double *Gj = CPU_G_own.data();

    gettimeofday(&cpu_start, NULL);
    offset = 0;
    int ioff_cart = 0;
    for( int i = 0; i < nshells; ++i) {
      XCPU::shells bra_shell = _shells[i];
      int bra_cart_size = (bra_shell.L + 1) * (bra_shell.L + 2) / 2;

      int joff_cart = 0;
      for( int j = 0; j <= i; ++j) {
	XCPU::shells ket_shell = _shells[j];
	int ket_cart_size = (ket_shell.L + 1) * (ket_shell.L + 2) / 2;

	XCPU::compute_integral_shell_pair(i == j,
					  ngrid,
					  _points_transposed.data(),
					  _shells[i].L,
					  _shells[j].L,
					  _shells[i].origin,
					  _shells[j].origin,
					  (_shells[i].m * _shells[j].m),
					  (prim_pairs + offset),
					  (Xi + ioff_cart * ngrid),
					  (Xj + joff_cart * ngrid),
					  ngrid,
					  (Gi + ioff_cart * ngrid),
					  (Gj + joff_cart * ngrid),
					  ngrid,
					  w.data(),
					  cpu_boys_table);

	offset += (_shells[i].m * _shells[j].m);
      
	joff_cart += ket_cart_size;
      }

      ioff_cart += bra_cart_size;
    }

    gettimeofday(&cpu_end, NULL);
    
    XCPU::boys_finalize(cpu_boys_table);
  }
 
  // GPU implementation

  std::vector<double> GPU_G_own( ngrid * nbf );
  for(int i = 0; i < ngrid * nbf; ++i) {
    GPU_G_own[i] = 0.0;
  }

  {
    double *dev_boys_table = XGPU::boys_init(); 

    std::vector<double> _points_transposed(3 * ngrid);
    for( int i = 0; i < ngrid; ++i ){
      _points_transposed[i + 0 * ngrid] = grid_points[i][0];
      _points_transposed[i + 1 * ngrid] = grid_points[i][1];
      _points_transposed[i + 2 * ngrid] = grid_points[i][2];
    }

    
    #if 0
    std::vector<XGPU::shells> _shells;
    _shells.resize(nshells);
    for( int i = 0; i < nshells; ++i ) {
      _shells[i].origin.x = basis[i].O()[0];
      _shells[i].origin.y = basis[i].O()[1];
      _shells[i].origin.z = basis[i].O()[2];

      _shells[i].m = basis[i].nprim();
      _shells[i].L = basis[i].l();

      _shells[i].coeff = new XGPU::coefficients[_shells[i].m];
      for( int j = 0; j < _shells[i].m; ++j ) {
	_shells[i].coeff[j].alpha = basis[i].alpha()[j];
	_shells[i].coeff[j].coeff = basis[i].coeff()[j];
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
    #else
    GauXC::ShellPairCollection shell_pairs(basis);
    #endif

    double *dev_points_transposed, *dev_X, *dev_G, *dev_weights;
    //XGPU::prim_pair *dev_prim_pairs;
    XGPU::shell_pair* dev_shell_pairs;
  
    cudaMalloc((void**) &dev_points_transposed, 3 * ngrid * sizeof(double));
    cudaMalloc((void**) &dev_X, ngrid * nbf * sizeof(double));
    cudaMalloc((void**) &dev_G, ngrid * nbf * sizeof(double));
    cudaMalloc((void**) &dev_weights, ngrid * sizeof(double));

    #if 0
    cudaMalloc((void**) &dev_prim_pairs, total_prim_pairs * sizeof(XGPU::prim_pair));
    #else
    cudaMalloc((void**) &dev_shell_pairs, shell_pairs.npairs() * sizeof(XGPU::shell_pair));
    #endif
  
    cudaMemcpy(dev_points_transposed, _points_transposed.data(), 3 * ngrid * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_X, F.data(), ngrid * nbf * sizeof(double), cudaMemcpyHostToDevice);
    //cudaMemcpy(dev_G, GPU_G_own.data(), ngrid * nbf * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset( dev_G, 0, ngrid*nbf*sizeof(double));
    cudaMemcpy(dev_weights, w.data(), ngrid * sizeof(double), cudaMemcpyHostToDevice);

    #if 0
    cudaMemcpy(dev_prim_pairs, prim_pairs, total_prim_pairs * sizeof(XGPU::prim_pair), cudaMemcpyHostToDevice);
    #else
    cudaMemcpy(dev_shell_pairs, shell_pairs.shell_pairs(), shell_pairs.npairs() * sizeof(XGPU::shell_pair), cudaMemcpyHostToDevice);
    #endif
    
    double *Xi = dev_X;
    double *Xj = dev_X;

    double *Gi = dev_G;
    double *Gj = dev_G;

    gettimeofday(&gpu_start, NULL);
    //offset = 0;
    int ioff_cart = 0;
    for( int i = 0; i < nshells; ++i) {
      #if 0
      XGPU::shells bra_shell = _shells[i];
      int bra_cart_size = (bra_shell.L + 1) * (bra_shell.L + 2) / 2;
      #else
      auto& bra_shell = basis[i];
      auto bra_cart_size = bra_shell.size();
      #endif

      int joff_cart = 0;
      for( int j = 0; j <= i; ++j) {
        #if 0
        XGPU::shells ket_shell = _shells[i];
        int ket_cart_size = (ket_shell.L + 1) * (ket_shell.L + 2) / 2;
        #else
        auto& ket_shell = basis[j];
        auto ket_cart_size = ket_shell.size();
        #endif
        
        XGPU::point A{bra_shell.O()[0], bra_shell.O()[1], bra_shell.O()[2]};
        XGPU::point B{ket_shell.O()[0], ket_shell.O()[1], ket_shell.O()[2]};
        auto sp = dev_shell_pairs + GauXC::detail::packed_lt_index(i,j,nshells);
        XGPU::compute_integral_shell_pair(i == j,
        				  ngrid,
        				  dev_points_transposed + 0*ngrid,
        				  dev_points_transposed + 1*ngrid,
        				  dev_points_transposed + 2*ngrid,
        				  bra_shell.l(), ket_shell.l(),
                  A, B,
        				  sp,
        				  (Xi + ioff_cart * ngrid),
        				  (Xj + joff_cart * ngrid),
        				  ngrid,
        				  (Gi + ioff_cart * ngrid),
        				  (Gj + joff_cart * ngrid),
        				  ngrid,
        				  dev_weights,
        				  dev_boys_table,0 );
        
        //offset += (_shells[i].m * _shells[j].m);
            
        joff_cart += ket_cart_size;
      }

      ioff_cart += bra_cart_size;
    }

    cudaDeviceSynchronize();
  
    gettimeofday(&gpu_end, NULL);

    cudaMemcpy(GPU_G_own.data(), dev_G, ngrid * nbf * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(dev_X);
    cudaFree(dev_G);
    cudaFree(dev_points_transposed);
    cudaFree(dev_weights);
    cudaFree(dev_shell_pairs);
  
    XGPU::boys_finalize(dev_boys_table);
  }

  int correct = 1;
  
  for( int i = 0; i < nbf * ngrid; ++i) {
    if((fabs(CPU_G_own[i] - GPU_G_own[i]) > 1e-10) || std::isnan(GPU_G_own[i]) || std::isnan(GPU_G_own[i])) {
      correct = 0;
    }
  }

  std::cout << "Correctness: " << correct << "\tCPU Execution: "<< 1000000 * (cpu_end.tv_sec - cpu_start.tv_sec) + (cpu_end.tv_usec - cpu_start.tv_usec) << "\tGPU Execution: "<< 1000000 * (gpu_end.tv_sec - gpu_start.tv_sec) + (gpu_end.tv_usec - gpu_start.tv_usec) << std::endl;
  
  return 0;
}
#else
int main() {}
#endif
