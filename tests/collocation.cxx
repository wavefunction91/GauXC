#include "ut_common.hpp"
#include <gauxc/molgrid.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/load_balancer.hpp>
#include <fstream>
#include <string>


#define MAX_NPTS_CHECK 67


#include "host/host_collocation.hpp"



#ifdef GAUXC_ENABLE_CUDA

#include <gauxc/exceptions/cuda_exception.hpp>
#include <gauxc/util/cuda_util.hpp>
#include "cuda/collocation_device.hpp"

#endif

#ifdef GAUXC_ENABLE_HIP

#include <gauxc/exceptions/hip_exception.hpp>
#include <gauxc/util/hip_util.hpp>
#include "hip/collocation_device.hpp"

#endif

#ifdef GAUXC_ENABLE_SYCL

#include <gauxc/exceptions/sycl_exception.hpp>
#include <gauxc/util/sycl_util.hpp>
#include "sycl/collocation_device.hpp"

#endif

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
#ifdef GAUXC_ENABLE_MPI
  LoadBalancer lb(MPI_COMM_WORLD, mol, mg, basis);
#else
  LoadBalancer lb(mol, mg, basis);
#endif
  auto& tasks = lb.get_tasks();


  std::vector< ref_collocation_data > ref_data;

  for( size_t i = 0; i < ntask_save; ++i ) {
    auto& task = tasks[i];

    auto& pts  = task.points;
    auto& mask = task.shell_list;

    // Only keep first MAX_NPTS_CHECK points to save on space
    if( task.points.size() > MAX_NPTS_CHECK )
      task.points.erase( task.points.begin() + MAX_NPTS_CHECK, task.points.end() );

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

    std::vector<double> eval( nbf * npts );


    integrator::host::eval_collocation( npts, mask.size(), nbf,
                                        pts.data()->data(), basis,
                                        mask.data(),
                                        eval.data() );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );

  }

}

void test_host_collocation_deriv1( const BasisSet<double>& basis, std::ifstream& in_file) {



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


#ifdef GAUXC_ENABLE_CUDA



void test_cuda_collocation_petite( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  auto shells_device  = util::cuda_malloc<Shell<double>>( basis.size() );
  auto offs_device    = util::cuda_malloc<size_t>( basis.size() );
  auto pts_device     = util::cuda_malloc<double>( 3 * MAX_NPTS_CHECK );
  auto eval_device    = util::cuda_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );



  cudaStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    util::cuda_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::cuda_copy( offs.size(), offs_device, offs.data()  );

    std::vector<Shell<double>> shells;
    for( auto idx : mask ) shells.emplace_back(basis[idx]);
    util::cuda_copy( shells.size(), shells_device, shells.data() );

    integrator::cuda::eval_collocation_petite( shells.size(), nbf, npts,
                                               shells_device, offs_device,
                                               pts_device,
                                               eval_device, stream );

    std::vector<double> eval( nbf * npts );

    util::cuda_copy( nbf * npts, eval.data(),    eval_device    );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );

  }
  util::cuda_device_sync();
  util::cuda_free(shells_device, offs_device, pts_device, eval_device );

}




void test_cuda_collocation_masked( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  auto shells_device  = util::cuda_malloc<Shell<double>>( basis.size() );
  auto offs_device    = util::cuda_malloc<size_t>( basis.size() );
  auto mask_device    = util::cuda_malloc<size_t>( basis.size() );
  auto pts_device     = util::cuda_malloc<double>( 3 * MAX_NPTS_CHECK );
  auto eval_device    = util::cuda_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );


  std::vector<Shell<double>> shells( basis );
  util::cuda_copy( basis.size(), shells_device, shells.data() );


  cudaStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    util::cuda_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::cuda_copy( mask.size(), mask_device, mask_ul.data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::cuda_copy( offs.size(), offs_device, offs.data()  );


    integrator::cuda::eval_collocation_masked( mask.size(), nbf, npts,
                                               shells_device, mask_device,
                                               offs_device, pts_device,
                                               eval_device, stream );

    std::vector<double> eval( nbf * npts );

    util::cuda_copy( nbf * npts, eval.data(),    eval_device    );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );

  }
  util::cuda_device_sync();
  util::cuda_free(shells_device, offs_device, pts_device, eval_device );

}















void test_cuda_collocation_petite_combined( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  std::vector< cuda::XCTaskDevice<double> > tasks;


  cudaStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.nbe     = nbf;
    task.npts    = npts;
    task.nshells = mask.size();

    task.points     = util::cuda_malloc<double>( 3 * npts );
    task.shell_offs = util::cuda_malloc<size_t>( mask.size() );
    task.shells     = util::cuda_malloc<Shell<double>>(mask.size());
    task.bf         = util::cuda_malloc<double>( nbf * npts );

    auto* pts_device = task.points;
    auto* offs_device = task.shell_offs;
    auto* shells_device = task.shells;


    util::cuda_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::cuda_copy( offs.size(), offs_device, offs.data()  );

    std::vector<Shell<double>> shells;
    for( auto idx : mask ) shells.emplace_back(basis[idx]);
    util::cuda_copy( shells.size(), shells_device, shells.data() );


  }


  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.nshells < b.nshells;
    })->nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::cuda_malloc<cuda::XCTaskDevice<double>>( tasks.size() );
  util::cuda_copy( tasks.size(), tasks_device, tasks.data() );

  integrator::cuda::eval_collocation_petite_combined( tasks.size(), npts_max,
    nshells_max, tasks_device, stream );

  util::cuda_device_sync();


  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    util::cuda_copy( eval.size(), eval.data(), tasks[i].bf );

    for( auto i = 0; i < eval.size(); ++i )
      CHECK( eval[i] == Approx( ref_eval[i] ) );
  }


  for( auto& t : tasks ) {
    util::cuda_free( t.points, t.shell_offs, t.shells, t.bf );
  }
  util::cuda_free( tasks_device );
}


void test_cuda_collocation_masked_combined( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  std::vector< cuda::XCTaskDevice<double> > tasks;

  auto shells_device  = util::cuda_malloc<Shell<double>>( basis.size() );
  std::vector<Shell<double>> shells( basis );
  util::cuda_copy( basis.size(), shells_device, shells.data() );

  cudaStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.nbe     = nbf;
    task.npts    = npts;
    task.nshells = mask.size();

    task.points     = util::cuda_malloc<double>( 3 * npts );
    task.shell_offs = util::cuda_malloc<size_t>( mask.size() );
    task.shell_list = util::cuda_malloc<size_t>( mask.size() );
    task.bf         = util::cuda_malloc<double>( nbf * npts );

    auto* pts_device = task.points;
    auto* offs_device = task.shell_offs;
    auto* mask_device = task.shell_list;


    util::cuda_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::cuda_copy( mask.size(), mask_device, mask_ul.data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::cuda_copy( offs.size(), offs_device, offs.data()  );


  }


  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.nshells < b.nshells;
    })->nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::cuda_malloc<cuda::XCTaskDevice<double>>( tasks.size() );
  util::cuda_copy( tasks.size(), tasks_device, tasks.data() );

  integrator::cuda::eval_collocation_masked_combined( tasks.size(), npts_max,
    nshells_max, shells_device, tasks_device, stream );

  util::cuda_device_sync();


  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    util::cuda_copy( eval.size(), eval.data(), tasks[i].bf );

    for( auto i = 0; i < eval.size(); ++i )
      CHECK( eval[i] == Approx( ref_eval[i] ) );
  }


  for( auto& t : tasks ) {
    util::cuda_free( t.points, t.shell_offs, t.shell_list, t.bf );
  }
  util::cuda_free( tasks_device, shells_device );
}



















void test_cuda_collocation_deriv1_petite( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  auto shells_device  = util::cuda_malloc<Shell<double>>( basis.size() );
  auto offs_device    = util::cuda_malloc<size_t>( basis.size() );
  auto pts_device     = util::cuda_malloc<double>( 3 * MAX_NPTS_CHECK );
  auto eval_device    = util::cuda_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );
  auto deval_device_x = util::cuda_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );
  auto deval_device_y = util::cuda_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );
  auto deval_device_z = util::cuda_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );



  cudaStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    util::cuda_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::cuda_copy( offs.size(), offs_device, offs.data()  );

    std::vector<Shell<double>> shells;
    for( auto idx : mask ) shells.emplace_back(basis[idx]);
    util::cuda_copy( shells.size(), shells_device, shells.data() );

    integrator::cuda::eval_collocation_petite_deriv1( shells.size(), nbf, npts,
                                                      shells_device, offs_device,
                                                      pts_device,
                                                      eval_device, deval_device_x,
                                                      deval_device_y, deval_device_z,
                                                      stream );

    std::vector<double> eval   ( nbf * npts ),
                        deval_x( nbf * npts ),
                        deval_y( nbf * npts ),
                        deval_z( nbf * npts );

    util::cuda_copy( nbf * npts, eval.data(),    eval_device    );
    util::cuda_copy( nbf * npts, deval_x.data(), deval_device_x );
    util::cuda_copy( nbf * npts, deval_y.data(), deval_device_y );
    util::cuda_copy( nbf * npts, deval_z.data(), deval_device_z );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_x[i] == Approx( d.deval_x[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_y[i] == Approx( d.deval_y[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_z[i] == Approx( d.deval_z[i] ) );

  }
  util::cuda_device_sync();
  util::cuda_free(shells_device, offs_device, pts_device, eval_device,
                 deval_device_x, deval_device_y, deval_device_z);
}




void test_cuda_collocation_deriv1_masked( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  auto shells_device  = util::cuda_malloc<Shell<double>>( basis.size() );
  auto offs_device    = util::cuda_malloc<size_t>( basis.size() );
  auto mask_device    = util::cuda_malloc<size_t>( basis.size() );
  auto pts_device     = util::cuda_malloc<double>( 3 * MAX_NPTS_CHECK );
  auto eval_device    = util::cuda_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );
  auto deval_device_x = util::cuda_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );
  auto deval_device_y = util::cuda_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );
  auto deval_device_z = util::cuda_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );


  std::vector<Shell<double>> shells( basis );
  util::cuda_copy( basis.size(), shells_device, shells.data() );

  cudaStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    util::cuda_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::cuda_copy( mask.size(), mask_device, mask_ul.data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::cuda_copy( offs.size(), offs_device, offs.data()  );


    integrator::cuda::eval_collocation_masked_deriv1( mask.size(), nbf, npts,
                                                      shells_device, mask_device,
                                                      offs_device, pts_device,
                                                      eval_device, deval_device_x,
                                                      deval_device_y, deval_device_z,
                                                      stream );

    std::vector<double> eval   ( nbf * npts ),
                        deval_x( nbf * npts ),
                        deval_y( nbf * npts ),
                        deval_z( nbf * npts );

    util::cuda_copy( nbf * npts, eval.data(),    eval_device    );
    util::cuda_copy( nbf * npts, deval_x.data(), deval_device_x );
    util::cuda_copy( nbf * npts, deval_y.data(), deval_device_y );
    util::cuda_copy( nbf * npts, deval_z.data(), deval_device_z );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_x[i] == Approx( d.deval_x[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_y[i] == Approx( d.deval_y[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_z[i] == Approx( d.deval_z[i] ) );

  }
  util::cuda_device_sync();
  util::cuda_free(shells_device, offs_device, pts_device, eval_device,
                 deval_device_x, deval_device_y, deval_device_z);
}







void test_cuda_collocation_petite_combined_deriv1( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  std::vector< cuda::XCTaskDevice<double> > tasks;


  cudaStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.nbe     = nbf;
    task.npts    = npts;
    task.nshells = mask.size();

    task.points     = util::cuda_malloc<double>( 3 * npts );
    task.shell_offs = util::cuda_malloc<size_t>( mask.size() );
    task.shells     = util::cuda_malloc<Shell<double>>(mask.size());
    task.bf         = util::cuda_malloc<double>( nbf * npts );
    task.dbfx       = util::cuda_malloc<double>( nbf * npts );
    task.dbfy       = util::cuda_malloc<double>( nbf * npts );
    task.dbfz       = util::cuda_malloc<double>( nbf * npts );

    auto* pts_device = task.points;
    auto* offs_device = task.shell_offs;
    auto* shells_device = task.shells;


    util::cuda_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::cuda_copy( offs.size(), offs_device, offs.data()  );

    std::vector<Shell<double>> shells;
    for( auto idx : mask ) shells.emplace_back(basis[idx]);
    util::cuda_copy( shells.size(), shells_device, shells.data() );


  }


  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.nshells < b.nshells;
    })->nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::cuda_malloc<cuda::XCTaskDevice<double>>( tasks.size() );
  util::cuda_copy( tasks.size(), tasks_device, tasks.data() );

  integrator::cuda::eval_collocation_petite_combined_deriv1( tasks.size(), npts_max,
    nshells_max, tasks_device, stream );

  util::cuda_device_sync();


  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    auto* ref_deval_x = ref_data[i].deval_x.data();
    auto* ref_deval_y = ref_data[i].deval_y.data();
    auto* ref_deval_z = ref_data[i].deval_z.data();

    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_x (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_y (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_z (tasks[i].nbe * tasks[i].npts);

    util::cuda_copy( eval.size(), eval.data(), tasks[i].bf );
    util::cuda_copy( eval.size(), deval_x.data(), tasks[i].dbfx );
    util::cuda_copy( eval.size(), deval_y.data(), tasks[i].dbfy );
    util::cuda_copy( eval.size(), deval_z.data(), tasks[i].dbfz );

    for( auto i = 0; i < eval.size(); ++i )
      CHECK( eval[i] == Approx( ref_eval[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_x[i] == Approx( ref_deval_x[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_y[i] == Approx( ref_deval_y[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_z[i] == Approx( ref_deval_z[i] ) );
  }


  for( auto& t : tasks ) {
    util::cuda_free( t.points, t.shell_offs, t.shells, t.bf, t.dbfx, t.dbfy,
      t.dbfz );
  }
  util::cuda_free( tasks_device );
}


void test_cuda_collocation_masked_combined_deriv1( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  std::vector< cuda::XCTaskDevice<double> > tasks;

  auto shells_device  = util::cuda_malloc<Shell<double>>( basis.size() );
  std::vector<Shell<double>> shells( basis );
  util::cuda_copy( basis.size(), shells_device, shells.data() );

  cudaStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.nbe     = nbf;
    task.npts    = npts;
    task.nshells = mask.size();

    task.points     = util::cuda_malloc<double>( 3 * npts );
    task.shell_offs = util::cuda_malloc<size_t>( mask.size() );
    task.shell_list = util::cuda_malloc<size_t>( mask.size() );
    task.bf         = util::cuda_malloc<double>( nbf * npts );
    task.dbfx       = util::cuda_malloc<double>( nbf * npts );
    task.dbfy       = util::cuda_malloc<double>( nbf * npts );
    task.dbfz       = util::cuda_malloc<double>( nbf * npts );


    auto* pts_device = task.points;
    auto* offs_device = task.shell_offs;
    auto* mask_device = task.shell_list;


    util::cuda_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::cuda_copy( mask.size(), mask_device, mask_ul.data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::cuda_copy( offs.size(), offs_device, offs.data()  );


  }


  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.nshells < b.nshells;
    })->nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::cuda_malloc<cuda::XCTaskDevice<double>>( tasks.size() );
  util::cuda_copy( tasks.size(), tasks_device, tasks.data() );

  integrator::cuda::eval_collocation_masked_combined_deriv1( tasks.size(), npts_max,
    nshells_max, shells_device, tasks_device, stream );

  util::cuda_device_sync();


  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    auto* ref_deval_x = ref_data[i].deval_x.data();
    auto* ref_deval_y = ref_data[i].deval_y.data();
    auto* ref_deval_z = ref_data[i].deval_z.data();

    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_x (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_y (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_z (tasks[i].nbe * tasks[i].npts);

    util::cuda_copy( eval.size(), eval.data(), tasks[i].bf );
    util::cuda_copy( eval.size(), deval_x.data(), tasks[i].dbfx );
    util::cuda_copy( eval.size(), deval_y.data(), tasks[i].dbfy );
    util::cuda_copy( eval.size(), deval_z.data(), tasks[i].dbfz );

    for( auto i = 0; i < eval.size(); ++i )
      CHECK( eval[i] == Approx( ref_eval[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_x[i] == Approx( ref_deval_x[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_y[i] == Approx( ref_deval_y[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_z[i] == Approx( ref_deval_z[i] ) );
  }


  for( auto& t : tasks ) {
    util::cuda_free( t.points, t.shell_offs, t.shell_list, t.bf, t.dbfx, t.dbfy,
      t.dbfz );
  }
  util::cuda_free( tasks_device, shells_device );
}
#endif

#ifdef GAUXC_ENABLE_HIP



void test_hip_collocation_petite( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  auto shells_device  = util::hip_malloc<Shell<double>>( basis.size() );
  auto offs_device    = util::hip_malloc<size_t>( basis.size() );
  auto pts_device     = util::hip_malloc<double>( 3 * MAX_NPTS_CHECK );
  auto eval_device    = util::hip_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );



  hipStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    util::hip_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::hip_copy( offs.size(), offs_device, offs.data()  );

    std::vector<Shell<double>> shells;
    for( auto idx : mask ) shells.emplace_back(basis[idx]);
    util::hip_copy( shells.size(), shells_device, shells.data() );

    integrator::hip::eval_collocation_petite( shells.size(), nbf, npts,
                                               shells_device, offs_device,
                                               pts_device,
                                               eval_device, stream );

    std::vector<double> eval( nbf * npts );

    util::hip_copy( nbf * npts, eval.data(),    eval_device    );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );

  }
  util::hip_device_sync();
  util::hip_free(shells_device, offs_device, pts_device, eval_device );

}




void test_hip_collocation_masked( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  auto shells_device  = util::hip_malloc<Shell<double>>( basis.size() );
  auto offs_device    = util::hip_malloc<size_t>( basis.size() );
  auto mask_device    = util::hip_malloc<size_t>( basis.size() );
  auto pts_device     = util::hip_malloc<double>( 3 * MAX_NPTS_CHECK );
  auto eval_device    = util::hip_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );


  std::vector<Shell<double>> shells( basis );
  util::hip_copy( basis.size(), shells_device, shells.data() );


  hipStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    util::hip_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::hip_copy( mask.size(), mask_device, mask_ul.data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::hip_copy( offs.size(), offs_device, offs.data()  );


    integrator::hip::eval_collocation_masked( mask.size(), nbf, npts,
                                               shells_device, mask_device,
                                               offs_device, pts_device,
                                               eval_device, stream );

    std::vector<double> eval( nbf * npts );

    util::hip_copy( nbf * npts, eval.data(),    eval_device    );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );

  }
  util::hip_device_sync();
  util::hip_free(shells_device, offs_device, pts_device, eval_device );

}















void test_hip_collocation_petite_combined( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  std::vector< hip::XCTaskDevice<double> > tasks;


  hipStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.nbe     = nbf;
    task.npts    = npts;
    task.nshells = mask.size();

    task.points     = util::hip_malloc<double>( 3 * npts );
    task.shell_offs = util::hip_malloc<size_t>( mask.size() );
    task.shells     = util::hip_malloc<Shell<double>>(mask.size());
    task.bf         = util::hip_malloc<double>( nbf * npts );

    auto* pts_device = task.points;
    auto* offs_device = task.shell_offs;
    auto* shells_device = task.shells;


    util::hip_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::hip_copy( offs.size(), offs_device, offs.data()  );

    std::vector<Shell<double>> shells;
    for( auto idx : mask ) shells.emplace_back(basis[idx]);
    util::hip_copy( shells.size(), shells_device, shells.data() );


  }


  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.nshells < b.nshells;
    })->nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::hip_malloc<hip::XCTaskDevice<double>>( tasks.size() );
  util::hip_copy( tasks.size(), tasks_device, tasks.data() );

  integrator::hip::eval_collocation_petite_combined( tasks.size(), npts_max,
    nshells_max, tasks_device, stream );

  util::hip_device_sync();


  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    util::hip_copy( eval.size(), eval.data(), tasks[i].bf );

    for( auto i = 0; i < eval.size(); ++i )
      CHECK( eval[i] == Approx( ref_eval[i] ) );
  }


  for( auto& t : tasks ) {
    util::hip_free( t.points, t.shell_offs, t.shells, t.bf );
  }
  util::hip_free( tasks_device );
}


void test_hip_collocation_masked_combined( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  std::vector< hip::XCTaskDevice<double> > tasks;

  auto shells_device  = util::hip_malloc<Shell<double>>( basis.size() );
  std::vector<Shell<double>> shells( basis );
  util::hip_copy( basis.size(), shells_device, shells.data() );

  hipStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.nbe     = nbf;
    task.npts    = npts;
    task.nshells = mask.size();

    task.points     = util::hip_malloc<double>( 3 * npts );
    task.shell_offs = util::hip_malloc<size_t>( mask.size() );
    task.shell_list = util::hip_malloc<size_t>( mask.size() );
    task.bf         = util::hip_malloc<double>( nbf * npts );

    auto* pts_device = task.points;
    auto* offs_device = task.shell_offs;
    auto* mask_device = task.shell_list;


    util::hip_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::hip_copy( mask.size(), mask_device, mask_ul.data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::hip_copy( offs.size(), offs_device, offs.data()  );


  }


  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.nshells < b.nshells;
    })->nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::hip_malloc<hip::XCTaskDevice<double>>( tasks.size() );
  util::hip_copy( tasks.size(), tasks_device, tasks.data() );

  integrator::hip::eval_collocation_masked_combined( tasks.size(), npts_max,
    nshells_max, shells_device, tasks_device, stream );

  util::hip_device_sync();


  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    util::hip_copy( eval.size(), eval.data(), tasks[i].bf );

    for( auto i = 0; i < eval.size(); ++i )
      CHECK( eval[i] == Approx( ref_eval[i] ) );
  }


  for( auto& t : tasks ) {
    util::hip_free( t.points, t.shell_offs, t.shell_list, t.bf );
  }
  util::hip_free( tasks_device, shells_device );
}



















void test_hip_collocation_deriv1_petite( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  auto shells_device  = util::hip_malloc<Shell<double>>( basis.size() );
  auto offs_device    = util::hip_malloc<size_t>( basis.size() );
  auto pts_device     = util::hip_malloc<double>( 3 * MAX_NPTS_CHECK );
  auto eval_device    = util::hip_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );
  auto deval_device_x = util::hip_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );
  auto deval_device_y = util::hip_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );
  auto deval_device_z = util::hip_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );



  hipStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    util::hip_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::hip_copy( offs.size(), offs_device, offs.data()  );

    std::vector<Shell<double>> shells;
    for( auto idx : mask ) shells.emplace_back(basis[idx]);
    util::hip_copy( shells.size(), shells_device, shells.data() );

    integrator::hip::eval_collocation_petite_deriv1( shells.size(), nbf, npts,
                                                      shells_device, offs_device,
                                                      pts_device,
                                                      eval_device, deval_device_x,
                                                      deval_device_y, deval_device_z,
                                                      stream );

    std::vector<double> eval   ( nbf * npts ),
                        deval_x( nbf * npts ),
                        deval_y( nbf * npts ),
                        deval_z( nbf * npts );

    util::hip_copy( nbf * npts, eval.data(),    eval_device    );
    util::hip_copy( nbf * npts, deval_x.data(), deval_device_x );
    util::hip_copy( nbf * npts, deval_y.data(), deval_device_y );
    util::hip_copy( nbf * npts, deval_z.data(), deval_device_z );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_x[i] == Approx( d.deval_x[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_y[i] == Approx( d.deval_y[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_z[i] == Approx( d.deval_z[i] ) );

  }
  util::hip_device_sync();
  util::hip_free(shells_device, offs_device, pts_device, eval_device,
                 deval_device_x, deval_device_y, deval_device_z);
}




void test_hip_collocation_deriv1_masked( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  auto shells_device  = util::hip_malloc<Shell<double>>( basis.size() );
  auto offs_device    = util::hip_malloc<size_t>( basis.size() );
  auto mask_device    = util::hip_malloc<size_t>( basis.size() );
  auto pts_device     = util::hip_malloc<double>( 3 * MAX_NPTS_CHECK );
  auto eval_device    = util::hip_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );
  auto deval_device_x = util::hip_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );
  auto deval_device_y = util::hip_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );
  auto deval_device_z = util::hip_malloc<double>( basis.nbf() * MAX_NPTS_CHECK );


  std::vector<Shell<double>> shells( basis );
  util::hip_copy( basis.size(), shells_device, shells.data() );

  hipStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    util::hip_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::hip_copy( mask.size(), mask_device, mask_ul.data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::hip_copy( offs.size(), offs_device, offs.data()  );


    integrator::hip::eval_collocation_masked_deriv1( mask.size(), nbf, npts,
                                                      shells_device, mask_device,
                                                      offs_device, pts_device,
                                                      eval_device, deval_device_x,
                                                      deval_device_y, deval_device_z,
                                                      stream );

    std::vector<double> eval   ( nbf * npts ),
                        deval_x( nbf * npts ),
                        deval_y( nbf * npts ),
                        deval_z( nbf * npts );

    util::hip_copy( nbf * npts, eval.data(),    eval_device    );
    util::hip_copy( nbf * npts, deval_x.data(), deval_device_x );
    util::hip_copy( nbf * npts, deval_y.data(), deval_device_y );
    util::hip_copy( nbf * npts, deval_z.data(), deval_device_z );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_x[i] == Approx( d.deval_x[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_y[i] == Approx( d.deval_y[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_z[i] == Approx( d.deval_z[i] ) );

  }
  util::hip_device_sync();
  util::hip_free(shells_device, offs_device, pts_device, eval_device,
                 deval_device_x, deval_device_y, deval_device_z);
}







void test_hip_collocation_petite_combined_deriv1( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  std::vector< hip::XCTaskDevice<double> > tasks;


  hipStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.nbe     = nbf;
    task.npts    = npts;
    task.nshells = mask.size();

    task.points     = util::hip_malloc<double>( 3 * npts );
    task.shell_offs = util::hip_malloc<size_t>( mask.size() );
    task.shells     = util::hip_malloc<Shell<double>>(mask.size());
    task.bf         = util::hip_malloc<double>( nbf * npts );
    task.dbfx       = util::hip_malloc<double>( nbf * npts );
    task.dbfy       = util::hip_malloc<double>( nbf * npts );
    task.dbfz       = util::hip_malloc<double>( nbf * npts );

    auto* pts_device = task.points;
    auto* offs_device = task.shell_offs;
    auto* shells_device = task.shells;


    util::hip_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::hip_copy( offs.size(), offs_device, offs.data()  );

    std::vector<Shell<double>> shells;
    for( auto idx : mask ) shells.emplace_back(basis[idx]);
    util::hip_copy( shells.size(), shells_device, shells.data() );


  }


  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.nshells < b.nshells;
    })->nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::hip_malloc<hip::XCTaskDevice<double>>( tasks.size() );
  util::hip_copy( tasks.size(), tasks_device, tasks.data() );

  integrator::hip::eval_collocation_petite_combined_deriv1( tasks.size(), npts_max,
    nshells_max, tasks_device, stream );

  util::hip_device_sync();


  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    auto* ref_deval_x = ref_data[i].deval_x.data();
    auto* ref_deval_y = ref_data[i].deval_y.data();
    auto* ref_deval_z = ref_data[i].deval_z.data();

    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_x (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_y (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_z (tasks[i].nbe * tasks[i].npts);

    util::hip_copy( eval.size(), eval.data(), tasks[i].bf );
    util::hip_copy( eval.size(), deval_x.data(), tasks[i].dbfx );
    util::hip_copy( eval.size(), deval_y.data(), tasks[i].dbfy );
    util::hip_copy( eval.size(), deval_z.data(), tasks[i].dbfz );

    for( auto i = 0; i < eval.size(); ++i )
      CHECK( eval[i] == Approx( ref_eval[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_x[i] == Approx( ref_deval_x[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_y[i] == Approx( ref_deval_y[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_z[i] == Approx( ref_deval_z[i] ) );
  }


  for( auto& t : tasks ) {
    util::hip_free( t.points, t.shell_offs, t.shells, t.bf, t.dbfx, t.dbfy,
      t.dbfz );
  }
  util::hip_free( tasks_device );
}


void test_hip_collocation_masked_combined_deriv1( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  std::vector< hip::XCTaskDevice<double> > tasks;

  auto shells_device  = util::hip_malloc<Shell<double>>( basis.size() );
  std::vector<Shell<double>> shells( basis );
  util::hip_copy( basis.size(), shells_device, shells.data() );

  hipStream_t stream = 0;
  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.nbe     = nbf;
    task.npts    = npts;
    task.nshells = mask.size();

    task.points     = util::hip_malloc<double>( 3 * npts );
    task.shell_offs = util::hip_malloc<size_t>( mask.size() );
    task.shell_list = util::hip_malloc<size_t>( mask.size() );
    task.bf         = util::hip_malloc<double>( nbf * npts );
    task.dbfx       = util::hip_malloc<double>( nbf * npts );
    task.dbfy       = util::hip_malloc<double>( nbf * npts );
    task.dbfz       = util::hip_malloc<double>( nbf * npts );


    auto* pts_device = task.points;
    auto* offs_device = task.shell_offs;
    auto* mask_device = task.shell_list;


    util::hip_copy( 3*npts, pts_device, pts.data()->data() );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::hip_copy( mask.size(), mask_device, mask_ul.data() );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::hip_copy( offs.size(), offs_device, offs.data()  );


  }


  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.nshells < b.nshells;
    })->nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::hip_malloc<hip::XCTaskDevice<double>>( tasks.size() );
  util::hip_copy( tasks.size(), tasks_device, tasks.data() );

  integrator::hip::eval_collocation_masked_combined_deriv1( tasks.size(), npts_max,
    nshells_max, shells_device, tasks_device, stream );

  util::hip_device_sync();


  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    auto* ref_deval_x = ref_data[i].deval_x.data();
    auto* ref_deval_y = ref_data[i].deval_y.data();
    auto* ref_deval_z = ref_data[i].deval_z.data();

    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_x (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_y (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_z (tasks[i].nbe * tasks[i].npts);

    util::hip_copy( eval.size(), eval.data(), tasks[i].bf );
    util::hip_copy( eval.size(), deval_x.data(), tasks[i].dbfx );
    util::hip_copy( eval.size(), deval_y.data(), tasks[i].dbfy );
    util::hip_copy( eval.size(), deval_z.data(), tasks[i].dbfz );

    for( auto i = 0; i < eval.size(); ++i )
      CHECK( eval[i] == Approx( ref_eval[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_x[i] == Approx( ref_deval_x[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_y[i] == Approx( ref_deval_y[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_z[i] == Approx( ref_deval_z[i] ) );
  }


  for( auto& t : tasks ) {
    util::hip_free( t.points, t.shell_offs, t.shell_list, t.bf, t.dbfx, t.dbfy,
      t.dbfz );
  }
  util::hip_free( tasks_device, shells_device );
}
#endif



#ifdef GAUXC_ENABLE_SYCL

void test_sycl_collocation_petite( const BasisSet<double>& basis, std::ifstream& in_file,
                                   cl::sycl::queue& syclQueue) {

  std::vector<ref_collocation_data> ref_data;
  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  auto shells_device  = util::sycl_malloc<Shell<double>>( basis.size(), syclQueue );
  auto offs_device    = util::sycl_malloc<size_t>( basis.size(), syclQueue );
  auto pts_device     = util::sycl_malloc<double>( 3 * MAX_NPTS_CHECK, syclQueue );
  auto eval_device    = util::sycl_malloc<double>( basis.nbf() * MAX_NPTS_CHECK, syclQueue );

  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    util::sycl_copy( 3*npts, pts_device, pts.data()->data(), syclQueue );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::sycl_copy( offs.size(), offs_device, offs.data(), syclQueue );

    std::vector<Shell<double>> shells;
    for( auto idx : mask ) shells.emplace_back(basis[idx]);
    util::sycl_copy( shells.size(), shells_device, shells.data(), syclQueue );

    integrator::sycl::eval_collocation_petite( shells.size(), nbf, npts,
                                               shells_device, offs_device,
                                               pts_device,
                                               eval_device, &syclQueue );

    std::vector<double> eval( nbf * npts );

    util::sycl_copy( nbf * npts, eval.data(), eval_device, syclQueue );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );

  }
  util::sycl_device_sync(syclQueue);
  util::sycl_free(shells_device, syclQueue);
  util::sycl_free(offs_device, syclQueue);
  util::sycl_free(pts_device, syclQueue);
  util::sycl_free(eval_device, syclQueue);
}


void test_sycl_collocation_masked( const BasisSet<double>& basis, std::ifstream& in_file,
                                   cl::sycl::queue& syclQueue) {

  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  auto shells_device  = util::sycl_malloc<Shell<double>>( basis.size(), syclQueue );
  auto offs_device    = util::sycl_malloc<size_t>( basis.size(), syclQueue );
  auto mask_device    = util::sycl_malloc<size_t>( basis.size(), syclQueue );
  auto pts_device     = util::sycl_malloc<double>( 3 * MAX_NPTS_CHECK, syclQueue );
  auto eval_device    = util::sycl_malloc<double>( basis.nbf() * MAX_NPTS_CHECK, syclQueue );

  std::vector<Shell<double>> shells( basis );
  util::sycl_copy( basis.size(), shells_device, shells.data(), syclQueue );

  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    util::sycl_copy( 3*npts, pts_device, pts.data()->data(), syclQueue );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::sycl_copy( mask.size(), mask_device, mask_ul.data(), syclQueue );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::sycl_copy( offs.size(), offs_device, offs.data(), syclQueue );


    integrator::sycl::eval_collocation_masked( mask.size(), nbf, npts,
                                               shells_device, mask_device,
                                               offs_device, pts_device,
                                               eval_device, &syclQueue );

    std::vector<double> eval( nbf * npts );

    util::sycl_copy( nbf * npts, eval.data(), eval_device, syclQueue );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );
  }

  util::sycl_device_sync(syclQueue);
  util::sycl_free(shells_device, syclQueue);
  util::sycl_free(offs_device, syclQueue);
  util::sycl_free(pts_device, syclQueue);
  util::sycl_free(eval_device, syclQueue);

}

void test_sycl_collocation_petite_combined( const BasisSet<double>& basis, std::ifstream& in_file,
                                            cl::sycl::queue& syclQueue) {

  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  std::vector< GauXC::sycl::XCTaskDevice<double> > tasks;

  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.nbe     = nbf;
    task.npts    = npts;
    task.nshells = mask.size();

    task.points     = util::sycl_malloc<double>( 3 * npts, syclQueue );
    task.shell_offs = util::sycl_malloc<size_t>( mask.size(), syclQueue );
    task.shells     = util::sycl_malloc<Shell<double>>(mask.size(), syclQueue);
    task.bf         = util::sycl_malloc<double>( nbf * npts, syclQueue );

    auto* pts_device = task.points;
    auto* offs_device = task.shell_offs;
    auto* shells_device = task.shells;

    util::sycl_copy( 3*npts, pts_device, pts.data()->data(), syclQueue );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::sycl_copy( offs.size(), offs_device, offs.data(), syclQueue );

    std::vector<Shell<double>> shells;
    for( auto idx : mask ) shells.emplace_back(basis[idx]);
    util::sycl_copy( shells.size(), shells_device, shells.data(), syclQueue );
  }

  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.nshells < b.nshells;
    })->nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::sycl_malloc<GauXC::sycl::XCTaskDevice<double>>( tasks.size(), syclQueue );
  util::sycl_copy( tasks.size(), tasks_device, tasks.data(), syclQueue );

  integrator::sycl::eval_collocation_petite_combined( tasks.size(), npts_max,
    nshells_max, tasks_device, &syclQueue );

  util::sycl_device_sync(syclQueue);


  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    util::sycl_copy( eval.size(), eval.data(), tasks[i].bf, syclQueue );

    for( auto i = 0; i < eval.size(); ++i )
      CHECK( eval[i] == Approx( ref_eval[i] ) );
  }


  for( auto& t : tasks ) {
      util::sycl_free( t.points, syclQueue);
      util::sycl_free( t.shell_offs, syclQueue);
      util::sycl_free( t.shells, syclQueue);
      util::sycl_free( t.bf, syclQueue );
  }
  util::sycl_free( tasks_device, syclQueue );
}

void test_sycl_collocation_masked_combined( const BasisSet<double>& basis, std::ifstream& in_file,
                                            cl::sycl::queue& syclQueue) {

  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  std::vector< GauXC::sycl::XCTaskDevice<double> > tasks;

  auto shells_device  = util::sycl_malloc<Shell<double>>( basis.size(), syclQueue );
  std::vector<Shell<double>> shells( basis );
  util::sycl_copy( basis.size(), shells_device, shells.data(), syclQueue );

  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.nbe     = nbf;
    task.npts    = npts;
    task.nshells = mask.size();

    task.points     = util::sycl_malloc<double>( 3 * npts, syclQueue );
    task.shell_offs = util::sycl_malloc<size_t>( mask.size(), syclQueue );
    task.shell_list = util::sycl_malloc<size_t>( mask.size(), syclQueue );
    task.bf         = util::sycl_malloc<double>( nbf * npts, syclQueue );

    auto* pts_device = task.points;
    auto* offs_device = task.shell_offs;
    auto* mask_device = task.shell_list;

    util::sycl_copy( 3*npts, pts_device, pts.data()->data(), syclQueue );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::sycl_copy( mask.size(), mask_device, mask_ul.data(), syclQueue );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::sycl_copy( offs.size(), offs_device, offs.data(), syclQueue );
  }

  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.nshells < b.nshells;
    })->nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::sycl_malloc<GauXC::sycl::XCTaskDevice<double>>( tasks.size(), syclQueue );
  util::sycl_copy( tasks.size(), tasks_device, tasks.data(), syclQueue );

  integrator::sycl::eval_collocation_masked_combined( tasks.size(), npts_max,
    nshells_max, shells_device, tasks_device, &syclQueue );

  util::sycl_device_sync(syclQueue);

  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    util::sycl_copy( eval.size(), eval.data(), tasks[i].bf, syclQueue );

    for( auto i = 0; i < eval.size(); ++i )
      CHECK( eval[i] == Approx( ref_eval[i] ) );
  }

  for( auto& t : tasks ) {
      util::sycl_free( t.points, syclQueue );
      util::sycl_free( t.shell_offs, syclQueue );
      util::sycl_free( t.shell_list, syclQueue );
      util::sycl_free( t.bf, syclQueue );
  }
  util::sycl_free( tasks_device, syclQueue );
  util::sycl_free( shells_device, syclQueue );
}

void test_sycl_collocation_deriv1_petite( const BasisSet<double>& basis, std::ifstream& in_file,
                                          cl::sycl::queue& syclQueue) {

  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  auto shells_device  = util::sycl_malloc<Shell<double>>( basis.size(), syclQueue );
  auto offs_device    = util::sycl_malloc<size_t>( basis.size(), syclQueue );
  auto pts_device     = util::sycl_malloc<double>( 3 * MAX_NPTS_CHECK, syclQueue );
  auto eval_device    = util::sycl_malloc<double>( basis.nbf() * MAX_NPTS_CHECK, syclQueue );
  auto deval_device_x = util::sycl_malloc<double>( basis.nbf() * MAX_NPTS_CHECK, syclQueue );
  auto deval_device_y = util::sycl_malloc<double>( basis.nbf() * MAX_NPTS_CHECK, syclQueue );
  auto deval_device_z = util::sycl_malloc<double>( basis.nbf() * MAX_NPTS_CHECK, syclQueue );

  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    util::sycl_copy( 3*npts, pts_device, pts.data()->data(), syclQueue );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::sycl_copy( offs.size(), offs_device, offs.data(), syclQueue );

    std::vector<Shell<double>> shells;
    for( auto idx : mask ) shells.emplace_back(basis[idx]);
    util::sycl_copy( shells.size(), shells_device, shells.data(), syclQueue );

    integrator::sycl::eval_collocation_petite_deriv1( shells.size(), nbf, npts,
                                                      shells_device, offs_device,
                                                      pts_device,
                                                      eval_device, deval_device_x,
                                                      deval_device_y, deval_device_z,
                                                      &syclQueue );

    std::vector<double> eval   ( nbf * npts ),
                        deval_x( nbf * npts ),
                        deval_y( nbf * npts ),
                        deval_z( nbf * npts );

    util::sycl_copy( nbf * npts, eval.data(),    eval_device,    syclQueue );
    util::sycl_copy( nbf * npts, deval_x.data(), deval_device_x, syclQueue );
    util::sycl_copy( nbf * npts, deval_y.data(), deval_device_y, syclQueue );
    util::sycl_copy( nbf * npts, deval_z.data(), deval_device_z, syclQueue );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_x[i] == Approx( d.deval_x[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_y[i] == Approx( d.deval_y[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_z[i] == Approx( d.deval_z[i] ) );
  }

  util::sycl_device_sync(syclQueue);

  util::sycl_free(shells_device, syclQueue);
  util::sycl_free(offs_device, syclQueue);
  util::sycl_free(pts_device, syclQueue);
  util::sycl_free(eval_device, syclQueue);
  util::sycl_free(deval_device_x, syclQueue);
  util::sycl_free(deval_device_y, syclQueue);
  util::sycl_free(deval_device_z, syclQueue);
}

void test_sycl_collocation_deriv1_masked( const BasisSet<double>& basis, std::ifstream& in_file,
                                          cl::sycl::queue& syclQueue) {

  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  auto shells_device  = util::sycl_malloc<Shell<double>>( basis.size(), syclQueue );
  auto offs_device    = util::sycl_malloc<size_t>( basis.size(), syclQueue );
  auto mask_device    = util::sycl_malloc<size_t>( basis.size(), syclQueue );
  auto pts_device     = util::sycl_malloc<double>( 3 * MAX_NPTS_CHECK, syclQueue );
  auto eval_device    = util::sycl_malloc<double>( basis.nbf() * MAX_NPTS_CHECK, syclQueue );
  auto deval_device_x = util::sycl_malloc<double>( basis.nbf() * MAX_NPTS_CHECK, syclQueue );
  auto deval_device_y = util::sycl_malloc<double>( basis.nbf() * MAX_NPTS_CHECK, syclQueue );
  auto deval_device_z = util::sycl_malloc<double>( basis.nbf() * MAX_NPTS_CHECK, syclQueue );

  std::vector<Shell<double>> shells( basis );
  util::sycl_copy( basis.size(), shells_device, shells.data(), syclQueue );

  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    util::sycl_copy( 3*npts, pts_device, pts.data()->data(), syclQueue );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::sycl_copy( mask.size(), mask_device, mask_ul.data(), syclQueue );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::sycl_copy( offs.size(), offs_device, offs.data(), syclQueue );

    integrator::sycl::eval_collocation_masked_deriv1( mask.size(), nbf, npts,
                                                      shells_device, mask_device,
                                                      offs_device, pts_device,
                                                      eval_device, deval_device_x,
                                                      deval_device_y, deval_device_z,
                                                      &syclQueue );

    std::vector<double> eval   ( nbf * npts ),
                        deval_x( nbf * npts ),
                        deval_y( nbf * npts ),
                        deval_z( nbf * npts );

    util::sycl_copy( nbf * npts, eval.data(),    eval_device,    syclQueue );
    util::sycl_copy( nbf * npts, deval_x.data(), deval_device_x, syclQueue );
    util::sycl_copy( nbf * npts, deval_y.data(), deval_device_y, syclQueue );
    util::sycl_copy( nbf * npts, deval_z.data(), deval_device_z, syclQueue );

    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( eval[i] == Approx( d.eval[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_x[i] == Approx( d.deval_x[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_y[i] == Approx( d.deval_y[i] ) );
    for( auto i = 0; i < npts * nbf; ++i )
      CHECK( deval_z[i] == Approx( d.deval_z[i] ) );

  }

  util::sycl_device_sync(syclQueue);

  util::sycl_free(shells_device, syclQueue);
  util::sycl_free(offs_device, syclQueue);
  util::sycl_free(pts_device, syclQueue);
  util::sycl_free(eval_device, syclQueue);
  util::sycl_free(deval_device_x, syclQueue);
  util::sycl_free(deval_device_y, syclQueue);
  util::sycl_free(deval_device_z, syclQueue);

}

void test_sycl_collocation_petite_combined_deriv1( const BasisSet<double>& basis, std::ifstream& in_file,
                                                   cl::sycl::queue& syclQueue) {

  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  std::vector< GauXC::sycl::XCTaskDevice<double> > tasks;

  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.nbe     = nbf;
    task.npts    = npts;
    task.nshells = mask.size();

    task.points     = util::sycl_malloc<double>( 3 * npts, syclQueue );
    task.shell_offs = util::sycl_malloc<size_t>( mask.size(), syclQueue );
    task.shells     = util::sycl_malloc<Shell<double>>(mask.size(), syclQueue);
    task.bf         = util::sycl_malloc<double>( nbf * npts, syclQueue );
    task.dbfx       = util::sycl_malloc<double>( nbf * npts, syclQueue );
    task.dbfy       = util::sycl_malloc<double>( nbf * npts, syclQueue );
    task.dbfz       = util::sycl_malloc<double>( nbf * npts, syclQueue );

    auto* pts_device = task.points;
    auto* offs_device = task.shell_offs;
    auto* shells_device = task.shells;

    util::sycl_copy( 3*npts, pts_device, pts.data()->data(), syclQueue );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::sycl_copy( offs.size(), offs_device, offs.data(), syclQueue );

    std::vector<Shell<double>> shells;
    for( auto idx : mask ) shells.emplace_back(basis[idx]);
    util::sycl_copy( shells.size(), shells_device, shells.data(), syclQueue );
  }

  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.nshells < b.nshells;
    })->nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::sycl_malloc<GauXC::sycl::XCTaskDevice<double>>( tasks.size(), syclQueue );
  util::sycl_copy( tasks.size(), tasks_device, tasks.data(), syclQueue );

  integrator::sycl::eval_collocation_petite_combined_deriv1( tasks.size(), npts_max,
    nshells_max, tasks_device, &syclQueue );

  util::sycl_device_sync(syclQueue);

  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    auto* ref_deval_x = ref_data[i].deval_x.data();
    auto* ref_deval_y = ref_data[i].deval_y.data();
    auto* ref_deval_z = ref_data[i].deval_z.data();

    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_x (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_y (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_z (tasks[i].nbe * tasks[i].npts);

    util::sycl_copy( eval.size(), eval.data(), tasks[i].bf, syclQueue );
    util::sycl_copy( eval.size(), deval_x.data(), tasks[i].dbfx, syclQueue );
    util::sycl_copy( eval.size(), deval_y.data(), tasks[i].dbfy, syclQueue );
    util::sycl_copy( eval.size(), deval_z.data(), tasks[i].dbfz, syclQueue );

    for( auto i = 0; i < eval.size(); ++i )
      CHECK( eval[i] == Approx( ref_eval[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_x[i] == Approx( ref_deval_x[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_y[i] == Approx( ref_deval_y[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_z[i] == Approx( ref_deval_z[i] ) );
  }

  for( auto& t : tasks ) {
      util::sycl_free( t.points, syclQueue );
      util::sycl_free( t.shell_offs, syclQueue );
      util::sycl_free( t.shells, syclQueue );
      util::sycl_free( t.bf, syclQueue );
      util::sycl_free( t.dbfx, syclQueue );
      util::sycl_free( t.dbfy, syclQueue );
      util::sycl_free( t.dbfz, syclQueue );
  }

  util::sycl_free( tasks_device, syclQueue );
}

void test_sycl_collocation_masked_combined_deriv1( const BasisSet<double>& basis, std::ifstream& in_file,
                                                   cl::sycl::queue& syclQueue) {

  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  std::vector< GauXC::sycl::XCTaskDevice<double> > tasks;

  auto shells_device  = util::sycl_malloc<Shell<double>>( basis.size(), syclQueue );
  std::vector<Shell<double>> shells( basis );
  util::sycl_copy( basis.size(), shells_device, shells.data(), syclQueue );

  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.nbe     = nbf;
    task.npts    = npts;
    task.nshells = mask.size();

    task.points     = util::sycl_malloc<double>( 3 * npts, syclQueue );
    task.shell_offs = util::sycl_malloc<size_t>( mask.size(), syclQueue );
    task.shell_list = util::sycl_malloc<size_t>( mask.size(), syclQueue );
    task.bf         = util::sycl_malloc<double>( nbf * npts, syclQueue );
    task.dbfx       = util::sycl_malloc<double>( nbf * npts, syclQueue );
    task.dbfy       = util::sycl_malloc<double>( nbf * npts, syclQueue );
    task.dbfz       = util::sycl_malloc<double>( nbf * npts, syclQueue );

    auto* pts_device = task.points;
    auto* offs_device = task.shell_offs;
    auto* mask_device = task.shell_list;

    util::sycl_copy( 3*npts, pts_device, pts.data()->data(), syclQueue );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::sycl_copy( mask.size(), mask_device, mask_ul.data(), syclQueue );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::sycl_copy( offs.size(), offs_device, offs.data(), syclQueue );
  }

  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.nshells < b.nshells;
    })->nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::sycl_malloc<GauXC::sycl::XCTaskDevice<double>>( tasks.size(), syclQueue );
  util::sycl_copy( tasks.size(), tasks_device, tasks.data(), syclQueue );

  integrator::sycl::eval_collocation_masked_combined_deriv1( tasks.size(), npts_max,
    nshells_max, shells_device, tasks_device, &syclQueue );

  util::sycl_device_sync(syclQueue);

  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    auto* ref_deval_x = ref_data[i].deval_x.data();
    auto* ref_deval_y = ref_data[i].deval_y.data();
    auto* ref_deval_z = ref_data[i].deval_z.data();

    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_x (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_y (tasks[i].nbe * tasks[i].npts);
    std::vector<double> deval_z (tasks[i].nbe * tasks[i].npts);

    util::sycl_copy( eval.size(), eval.data(), tasks[i].bf, syclQueue );
    util::sycl_copy( eval.size(), deval_x.data(), tasks[i].dbfx, syclQueue );
    util::sycl_copy( eval.size(), deval_y.data(), tasks[i].dbfy, syclQueue );
    util::sycl_copy( eval.size(), deval_z.data(), tasks[i].dbfz, syclQueue );

    for( auto i = 0; i < eval.size(); ++i )
      CHECK( eval[i] == Approx( ref_eval[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_x[i] == Approx( ref_deval_x[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_y[i] == Approx( ref_deval_y[i] ) );
    for( auto i = 0; i < eval.size(); ++i )
      CHECK( deval_z[i] == Approx( ref_deval_z[i] ) );
  }

  for( auto& t : tasks ) {
      util::sycl_free( t.points, syclQueue );
      util::sycl_free( t.shell_offs, syclQueue );
      util::sycl_free( t.shell_list, syclQueue );
      util::sycl_free( t.bf, syclQueue );
      util::sycl_free( t.dbfx, syclQueue );
      util::sycl_free( t.dbfy, syclQueue );
      util::sycl_free( t.dbfz, syclQueue );
  }
  util::sycl_free( tasks_device, syclQueue );
  util::sycl_free( shells_device, syclQueue );
}

#endif // GAUXC_ENABLE_SYCL





//#define GENERATE_TESTS

TEST_CASE( "Water / cc-pVDZ", "[collocation]" ) {

#ifdef GENERATE_TESTS
#ifdef GAUXC_ENABLE_MPI
  int world_size;
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
  if( world_size > 1 ) return;
#endif
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

  SECTION( "Host Eval Grad" ) {
    test_host_collocation_deriv1( basis, ref_data );
  }

#ifdef GAUXC_ENABLE_CUDA
  SECTION( "CUDA Eval: Petite Shell List" ) {
    test_cuda_collocation_petite( basis, ref_data );
  }
  SECTION( "CUDA Eval: Masked" ) {
    test_cuda_collocation_masked( basis, ref_data );
  }
  SECTION( "CUDA Eval: Petite Combined" ) {
    test_cuda_collocation_petite_combined( basis, ref_data );
  }
  SECTION( "CUDA Eval: Masked Combined" ) {
    test_cuda_collocation_masked_combined( basis, ref_data );
  }

  SECTION( "CUDA Eval Grad: Petite Shell List" ) {
    test_cuda_collocation_deriv1_petite( basis, ref_data );
  }
  SECTION( "CUDA Eval Grad: Masked" ) {
    test_cuda_collocation_deriv1_masked( basis, ref_data );
  }
  SECTION( "CUDA Eval Grad: Petite Combined" ) {
    test_cuda_collocation_petite_combined_deriv1( basis, ref_data );
  }
  SECTION( "CUDA Eval: Masked Combined" ) {
    test_cuda_collocation_masked_combined_deriv1( basis, ref_data );
  }
#endif // GAUXC_ENABLE_CUDA

#ifdef GAUXC_ENABLE_HIP
  SECTION( "HIP Eval: Petite Shell List" ) {
    test_hip_collocation_petite( basis, ref_data );
  }
  SECTION( "HIP Eval: Masked" ) {
    test_hip_collocation_masked( basis, ref_data );
  }
  SECTION( "HIP Eval: Petite Combined" ) {
    test_hip_collocation_petite_combined( basis, ref_data );
  }
  SECTION( "HIP Eval: Masked Combined" ) {
    test_hip_collocation_masked_combined( basis, ref_data );
  }

  SECTION( "HIP Eval Grad: Petite Shell List" ) {
    test_hip_collocation_deriv1_petite( basis, ref_data );
  }
  SECTION( "HIP Eval Grad: Masked" ) {
    test_hip_collocation_deriv1_masked( basis, ref_data );
  }
  SECTION( "HIP Eval Grad: Petite Combined" ) {
    test_hip_collocation_petite_combined_deriv1( basis, ref_data );
  }
  SECTION( "HIP Eval: Masked Combined" ) {
    test_hip_collocation_masked_combined_deriv1( basis, ref_data );
  }
#endif // GAUXC_ENABLE_HIP


#ifdef GAUXC_ENABLE_SYCL
  cl::sycl::gpu_selector device_selector;
  cl::sycl::queue syclQueue = cl::sycl::queue(device_selector,
                                              cl::sycl::property_list{cl::sycl::property::queue::in_order{}});

  SECTION( "SYCL Eval: Petite Shell List" ) {
    test_sycl_collocation_petite( basis, ref_data, syclQueue );
  }
  syclQueue.wait_and_throw();
  SECTION( "SYCL Eval: Masked" ) {
    test_sycl_collocation_masked( basis, ref_data, syclQueue );
  }
  syclQueue.wait_and_throw();
  SECTION( "SYCL Eval: Petite Combined" ) {
    test_sycl_collocation_petite_combined( basis, ref_data, syclQueue );
  }
  syclQueue.wait_and_throw();
  SECTION( "SYCL Eval: Masked Combined" ) {
    test_sycl_collocation_masked_combined( basis, ref_data, syclQueue );
  }
  syclQueue.wait_and_throw();

  SECTION( "SYCL Eval Grad: Petite Shell List" ) {
    test_sycl_collocation_deriv1_petite( basis, ref_data, syclQueue );
  }
  syclQueue.wait_and_throw();
  SECTION( "SYCL Eval Grad: Masked" ) {
    test_sycl_collocation_deriv1_masked( basis, ref_data, syclQueue );
  }
  syclQueue.wait_and_throw();
  SECTION( "SYCL Eval Grad: Petite Combined" ) {
    test_sycl_collocation_petite_combined_deriv1( basis, ref_data, syclQueue );
  }
  syclQueue.wait_and_throw();
  SECTION( "SYCL Eval: Masked Combined" ) {
    test_sycl_collocation_masked_combined_deriv1( basis, ref_data, syclQueue );
  }
  syclQueue.wait_and_throw();
#endif // GAUXC_ENABLE_SYCL

#endif

}
