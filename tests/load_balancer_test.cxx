#include "ut_common.hpp"
#include <gauxc/load_balancer.hpp>

using namespace GauXC;

TEST_CASE( "DefaultLoadBalancer", "[load_balancer]" ) {

  MPI_Comm comm          = MPI_COMM_WORLD;
  Molecule mol           = make_benzene();
  BasisSet<double> basis = make_ccpvdz( mol, SphericalType(true) );

  for( auto& sh : basis ) 
    sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );

  MolGrid mg(AtomicGridSizeDefault::UltraFineGrid, mol);

  LoadBalancer lb(comm, mol, mg, basis);

  auto& tasks = lb.get_tasks();

  std::string ref_file = GAUXC_REF_DATA_PATH "/benzene_cc-pvdz_ufg_tasks.bin";
#ifdef GAUXC_GEN_TESTS

  {
    // Points / Weights not stored in reference data to 
    // save space
    for( auto& t : tasks ) {
      t.points.clear();
      t.weights.clear();
    }

    std::ofstream of( ref_file, std::ios::binary );
    cereal::BinaryOutputArchive ar(of);
    ar( tasks );
  }

#else

  std::vector<XCTask> ref_tasks;
  {
    std::ifstream ifile( ref_file, std::ios::binary );
    cereal::BinaryInputArchive ar(ifile);
    ar( ref_tasks );
  }

  REQUIRE( tasks.size() == ref_tasks.size() );

  size_t ntasks = tasks.size();
  for( size_t i = 0; i < ntasks; ++i ) {

    const auto& t  = tasks[i]; 
    const auto& rt = ref_tasks[i];
    CHECK( t.iParent == rt.iParent );
    CHECK( t.shell_list == rt.shell_list );
    CHECK( t.nbe == rt.nbe );
    CHECK( t.dist_nearest == Approx(rt.dist_nearest) );

    /* 
    // Points / Weights not stored in reference data to 
    // save space
    REQUIRE( t.points.size() == rt.points.size() );
    size_t npts = t.points.size();
    for( size_t j = 0; j < npts; ++j ) {
      CHECK( t.points[j][0] == Approx(rt.points[j][0]) );
      CHECK( t.points[j][1] == Approx(rt.points[j][1]) );
      CHECK( t.points[j][2] == Approx(rt.points[j][2]) );
      CHECK( t.weights[j] == Approx(rt.weights[j]) );
    }
    */

  }

#endif


}
