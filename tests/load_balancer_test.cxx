#include "ut_common.hpp"
#include <gauxc/load_balancer.hpp>

using namespace GauXC;


void gen_ref_lb_data( std::vector<XCTask>& tasks ) {

  int world_size;
  int world_rank;
#ifdef GAUXC_ENABLE_MPI
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
  MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
#else
  world_size = 1;
  world_rank = 0;
#endif
  std::string ref_file = GAUXC_REF_DATA_PATH "/benzene_cc-pvdz_ufg_tasks_" + std::to_string(world_size) + "mpi_rank" + std::to_string(world_rank) + ".bin";

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

void check_lb_data( const std::vector<XCTask>& tasks ) {

  int world_size;
  int world_rank;
#ifdef GAUXC_ENABLE_MPI
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
  MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
#else
  world_size = 1;
  world_rank = 0;
#endif
  std::string ref_file = GAUXC_REF_DATA_PATH "/benzene_cc-pvdz_ufg_tasks_" + std::to_string(world_size) + "mpi_rank" + std::to_string(world_rank) + ".bin";

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

}


//#define GAUXC_GEN_TESTS
TEST_CASE( "DefaultLoadBalancer", "[load_balancer]" ) {

#ifdef GAUXC_ENABLE_MPI
  MPI_Comm comm          = MPI_COMM_WORLD;
  #define GAUXC_MPI_ARG comm,
#else
  #define GAUXC_MPI_ARG 
#endif

  Molecule mol           = make_benzene();
  BasisSet<double> basis = make_ccpvdz( mol, SphericalType(true) );

  for( auto& sh : basis ) 
    sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );

  MolGrid mg(AtomicGridSizeDefault::UltraFineGrid, mol);

  auto meta = std::make_shared<MolMeta>( mol );

#ifdef GAUXC_GEN_TESTS

  LoadBalancerFactory lb_factory( ExecutionSpace::Host, "Default" );
  auto lb = lb_factory.get_instance( GAUXC_MPI_ARG mol, mg, basis);
  auto& tasks = lb.get_tasks();
  gen_ref_lb_data(tasks);

#else

  SECTION("Default Host") {

    LoadBalancerFactory lb_factory( ExecutionSpace::Host, "Default" );
    auto lb = lb_factory.get_instance( GAUXC_MPI_ARG mol, mg, basis);
    auto& tasks = lb.get_tasks();
    check_lb_data( tasks );

  }

#ifdef GAUXC_ENABLE_DEVICE
  SECTION("Default Device") {

    LoadBalancerFactory lb_factory( ExecutionSpace::Device, "Default" );
    auto lb = lb_factory.get_instance( GAUXC_MPI_ARG mol, mg, basis);
    auto& tasks = lb.get_tasks();
    check_lb_data( tasks );

  }
#endif

#endif


}
