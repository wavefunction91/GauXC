#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <integral_data_types.hpp>
#include <obara_saika_integrals.hpp>
#include <chebyshev_boys_computation.hpp>
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

  if( argc < 2 ) throw std::runtime_error("NOT VALID INPUT");

  std::string data_file = argv[1];

  GauXC::Molecule mol;
  GauXC::read_hdf5_record(mol, data_file, "/MOLECULE");

  GauXC::BasisSet<double> basis;
  GauXC::read_hdf5_record(basis, data_file, "/BASIS");
  for( auto& sh : basis ) sh.set_pure(false); // Reset to cartesian

  std::cout << mol.size() << std::endl; 
  std::cout << basis.size() << std::endl;
  std::cout << basis.nbf() << std::endl;
}
