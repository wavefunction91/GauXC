#include "ut_common.hpp"
#include <gauxc/molgrid.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/load_balancer.hpp>
#include <fstream>
#include <mpi.h>

#include "host/host_collocation.hpp"

using namespace GauXC;


struct ref_collocation_data {
  std::vector<int32_t>              mask;
  std::vector<std::array<double,3>> pts;
  std::vector<double>               eval;
  std::vector<double>               deval_x;
  std::vector<double>               deval_y;
  std::vector<double>               deval_z;

  template <typename Archive>
  void serialize( Archive& ar ) {
    ar( mask, pts, eval, deval_x, deval_y, deval_z );
  }

};

void generate_collocation_data( const Molecule& mol, const BasisSet<double>& basis,
                                std::ofstream& out_file, size_t ntask_save = 10 ) {


  MolGrid mg(AtomicGridSizeDefault::FineGrid, mol);
  LoadBalancer lb(MPI_COMM_WORLD, mol, mg, basis); 
  auto& tasks = lb.get_tasks();


  std::vector< ref_collocation_data > ref_data;
  
  for( size_t i = 0; i < ntask_save; ++i ) {
    auto& task = tasks[i];

    auto& pts  = task.points;
    auto& mask = task.shell_list;

    // Only keep first 67 points to save on space
    if( task.points.size() > 67 )
      task.points.erase( task.points.begin() + 67, task.points.end() );

    const auto npts = task.points.size();
    const auto nbf  = task.nbe;

    std::vector<double> eval   ( nbf * npts ),
                        deval_x( nbf * npts ),
                        deval_y( nbf * npts ),
                        deval_z( nbf * npts );

    integrator::host::eval_collocation_deriv1( npts, mask.size(), nbf,
                                               pts.data()->data(), basis, 
                                               mask.data(), 
                                               eval.data(), deval_x.data(),
                                               deval_y.data(), deval_z.data() );

    auto max_abs = *std::max_element( eval.begin(), eval.end(), 
                   [](auto a, auto b){ return std::abs(a) < std::abs(b); } );
    if( std::abs(max_abs) < 1e-9 ) continue;

    ref_collocation_data d{ std::move(mask), std::move(pts), std::move(eval), 
                            std::move(deval_x), std::move(deval_y), 
                            std::move(deval_z) };

    ref_data.emplace_back( std::move(d) );

  }

  {
    cereal::BinaryOutputArchive ar( out_file );
    ar( ref_data );
  }

}



void test_host_collocation( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  for( auto& d : ref_data ) {
  
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    std::vector<double> eval   ( nbf * npts ),
                        deval_x( nbf * npts ),
                        deval_y( nbf * npts ),
                        deval_z( nbf * npts );


    integrator::host::eval_collocation_deriv1( npts, mask.size(), nbf,
                                               pts.data()->data(), basis, 
                                               mask.data(), 
                                               eval.data(), deval_x.data(),
                                               deval_y.data(), deval_z.data() );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_x[i] == Approx( d.deval_x[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_y[i] == Approx( d.deval_y[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_z[i] == Approx( d.deval_z[i] ) );
  }

}

//#define GENERATE_TESTS

TEST_CASE( "Water / cc-pVDZ", "[collocation]" ) {

#ifdef GENERATE_TESTS
  int world_size;
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
  if( world_size > 1 ) return;
#endif

  Molecule mol           = make_water();
  BasisSet<double> basis = make_ccpvdz( mol, SphericalType(true) );

  for( auto& sh : basis ) sh.set_shell_tolerance( 1e-6 );

#ifdef GENERATE_TESTS
  std::ofstream ref_data( "water_cc-pVDZ_collocation.bin", std::ios::binary );
  generate_collocation_data( mol, basis, ref_data );  
#else
  std::ifstream ref_data( GAUXC_REF_DATA_PATH "/water_cc-pVDZ_collocation.bin", 
                          std::ios::binary );

  SECTION( "Host Eval" ) {
    test_host_collocation( basis, ref_data );
  }
#endif

}
