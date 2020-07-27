#include "catch2/catch.hpp"
#include "standards.hpp"

#include <gauxc/load_balancer.hpp>

#include <random>
#include <algorithm>

using namespace GauXC;

TEST_CASE( "DefaultLoadBalancer", "[load_balancer]" ) {

  MPI_Comm comm          = MPI_COMM_WORLD;
  Molecule mol           = make_benzene();
  BasisSet<double> basis = make_ccpvdz( mol, SphericalType(true) );

  for( auto& sh : basis ) sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );

  for( auto i = 0; i < basis.size(); ++i ) {
    const auto& sh = basis[i];
    std::cout << "CEN = " << sh.O()[0] << ", " << sh.O()[1] << ", " << sh.O()[2] << std::endl;
    std::cout << "L = " << sh.l() << std::endl;
    std::cout << "CR = " << sh.cutoff_radius() << std::endl;
    std::cout << "PRIMS" << std::endl;
    for( auto p = 0; p < sh.nprim(); ++p )
      std::cout << "  " << sh.alpha()[p] << ", " << sh.coeff()[p] << std::endl;
    std::cout << std::endl;
  }

  std::cout << basis.nbf() << std::endl;

  MolGrid mg(AtomicGridSizeDefault::UltraFineGrid, mol);

  LoadBalancer lb(comm, mol, mg, basis);

  auto& tasks = lb.get_tasks();

  std::cout << "TASKS" << std::endl;
  for( auto& task : tasks ) {
    std::cout << task.iParent << ", ";
    std::cout << task.nbe << ", ";
    std::cout << task.dist_nearest << ", ";
    for( auto& sh : task.shell_list )
      std::cout << sh << ", ";
    std::cout << std::endl;
  }

}
