#include "ut_common.hpp"
#include <gauxc/molgrid.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/load_balancer.hpp>
#include <fstream>
#include <string>

#ifdef GAUXC_ENABLE_HOST
#include "host/host_weights.hpp"
#endif

#ifdef GAUXC_ENABLE_CUDA
#include <gauxc/exceptions/cuda_exception.hpp>
#include <gauxc/util/cuda_util.hpp>
#include "cuda/cuda_weights.hpp"
#endif

using namespace GauXC;

struct ref_weights_data {
  Molecule                  mol;
  std::shared_ptr<MolMeta>  meta;
  std::vector< XCTask > tasks_unm;
  std::vector< XCTask > tasks_mod; // This is only the weights

  template <typename Archive>
  void load( Archive& ar ) {
    ar( mol, tasks_unm, tasks_mod );
    meta = std::make_shared<MolMeta>(mol);
  }
  template <typename Archive>
  void save( Archive& ar ) const {
    ar( mol, tasks_unm, tasks_mod );
  }
};


#ifdef GAUXC_ENABLE_HOST
void generate_weights_data( const Molecule& mol, const BasisSet<double>& basis,
                                std::ofstream& out_file, size_t ntask_save = 15 ) {


  MolGrid mg(AtomicGridSizeDefault::FineGrid, mol);
#ifdef GAUXC_ENABLE_MPI
  LoadBalancer lb(MPI_COMM_WORLD, mol, mg, basis);
#else
  LoadBalancer lb(mol, mg, basis);
#endif
  auto& tasks = lb.get_tasks();

  ref_weights_data   ref_data;
  ref_data.mol       = mol;

  auto abs_comparator = []( const auto& a, const auto& b ) {
    return std::abs(a) < std::abs(b);
  };

  std::sort( tasks.begin(), tasks.end(),
    [&]( const auto& a, const auto& b ) {
      auto a_max =
        *std::max_element( a.weights.begin(), a.weights.end(),
                           abs_comparator );
      auto b_max =
        *std::max_element( b.weights.begin(), b.weights.end(),
                           abs_comparator );

      return a_max < b_max;
    });

  if( tasks.size() > ntask_save )
    tasks.erase( tasks.begin() + ntask_save, tasks.end() );

  ref_data.tasks_unm = tasks; // Make a copy of un modified tasks

  integrator::host::partition_weights_host( XCWeightAlg::SSF,
    mol, lb.molmeta(), tasks );

  // Clear out unneeded data
  for( auto& task : tasks ) {
    task.points.clear();
    task.shell_list.clear();
  }
  ref_data.tasks_mod = tasks;

  {
    cereal::BinaryOutputArchive ar( out_file );
    ar( ref_data );
  }

}


void test_host_weights( std::ifstream& in_file ) {

  ref_weights_data ref_data;
  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  integrator::host::partition_weights_host( XCWeightAlg::SSF,
    ref_data.mol, *ref_data.meta, ref_data.tasks_unm );


  size_t ntasks = ref_data.tasks_unm.size();
  for( size_t itask = 0; itask < ntasks; ++itask ) {
    auto& task     = ref_data.tasks_unm.at(itask);
    auto& ref_task = ref_data.tasks_mod.at(itask);

    size_t npts = task.weights.size();
    for( size_t i = 0; i < npts; ++i ) {
      CHECK( task.weights.at(i) ==
             Approx(ref_task.weights.at(i)) );
    }
  }

}
#endif

#ifdef GAUXC_ENABLE_CUDA
void test_cuda_weights( std::ifstream& in_file ) {

  ref_weights_data ref_data;
  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  std::vector< std::array<double,3> > points;
  std::vector< double >               weights, weights_ref;
  std::vector< double >               dist_nearest;
  std::vector< int32_t >              iparent;

  for( auto& task : ref_data.tasks_unm ) {
    points.insert( points.end(),
                   task.points.begin(),
                   task.points.end() );
    weights.insert( weights.end(),
                    task.weights.begin(),
                    task.weights.end() );

    size_t npts = task.points.size();
    dist_nearest.insert( dist_nearest.end(), npts,
                         task.dist_nearest );
    iparent.insert( iparent.end(), npts, task.iParent );
  }

  for( auto& task : ref_data.tasks_mod ) {
    weights_ref.insert( weights_ref.end(),
                        task.weights.begin(),
                        task.weights.end() );
  }

  size_t npts   = points.size();
  size_t natoms = ref_data.mol.natoms();

  std::vector< double >  coords( 3 * natoms );
  for( auto iat = 0 ; iat < natoms; ++iat ) {
    coords[ 3*iat + 0 ] = ref_data.mol.at(iat).x;
    coords[ 3*iat + 1 ] = ref_data.mol.at(iat).y;
    coords[ 3*iat + 2 ] = ref_data.mol.at(iat).z;
  }


  auto* points_d  = util::cuda_malloc<double>( 3*npts );
  auto* weights_d = util::cuda_malloc<double>( npts   );
  auto* iparent_d = util::cuda_malloc<int32_t>( npts  );
  auto* distnea_d = util::cuda_malloc<double>( npts   );
  auto* rab_d     = util::cuda_malloc<double>( natoms*natoms );
  auto* coords_d  = util::cuda_malloc<double>( 3*natoms );
  auto* dist_scr_d= util::cuda_malloc<double>( npts*natoms );

  util::cuda_copy( 3*npts, points_d,  points.data()->data() );
  util::cuda_copy( npts,   weights_d, weights.data() );
  util::cuda_copy( npts,   iparent_d, iparent.data() );
  util::cuda_copy( npts,   distnea_d, dist_nearest.data() );
  util::cuda_copy( natoms*natoms, rab_d,
                   ref_data.meta->rab().data() );
  util::cuda_copy( 3*natoms, coords_d, coords.data() );

  cudaStream_t stream = 0;
  integrator::cuda::partition_weights_cuda_SoA(
    XCWeightAlg::SSF, npts, natoms, points_d,
    iparent_d, distnea_d, rab_d, coords_d,
    weights_d, dist_scr_d, stream );

  util::cuda_device_sync();
  util::cuda_copy( npts, weights.data(), weights_d );
  util::cuda_free( points_d, weights_d, iparent_d, distnea_d,
                   rab_d, coords_d, dist_scr_d );

  for( auto i = 0ul; i < npts; ++i )
    CHECK( weights.at(i) == Approx( weights_ref.at(i) ) );

}
#endif

//#define GENERATE_TESTS
TEST_CASE( "Benzene", "[weights]" ) {

#ifdef GENERATE_TESTS
#ifdef GAUXC_ENABLE_MPI
  int world_size;
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
  if( world_size > 1 ) return;
#endif
#endif

  Molecule mol = make_benzene();

#ifdef GENERATE_TESTS
  BasisSet<double> basis = make_631Gd( mol, SphericalType(true) );
  for( auto& sh : basis ) sh.set_shell_tolerance( 1e-6 );

  std::ofstream ref_data( "benzene_weights_ssf.bin", std::ios::binary );
  generate_weights_data( mol, basis, ref_data );  
#else

  std::ifstream ref_data( GAUXC_REF_DATA_PATH "/benzene_weights_ssf.bin", 
                          std::ios::binary );

#ifdef GAUXC_ENABLE_HOST
  SECTION( "Host Weights" ) {
    test_host_weights( ref_data );
  }
#endif

#ifdef GAUXC_ENABLE_CUDA
  SECTION( "Device Weights" ) {
    test_cuda_weights( ref_data );
  }
#endif

#endif
}


