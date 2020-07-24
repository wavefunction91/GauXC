#include "catch2/catch.hpp"
#include "standards.hpp"

#include <gauxc/load_balancer.hpp>

#include <random>
#include <algorithm>

using namespace GauXC;

TEST_CASE( "DefaultLoadBalancer", "[load_balancer]" ) {

  MPI_Comm comm          = MPI_COMM_WORLD;
  Molecule mol           = make_water();
  BasisSet<double> basis = make_631Gd( mol, SphericalType(true) );

  MolGrid mg(AtomicGridSizeDefault::UltraFineGrid, mol);

  LoadBalancer lb(comm, mol, mg, basis);

  auto& tasks = lb.get_tasks();

}
