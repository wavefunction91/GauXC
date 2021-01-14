#ifdef GAUXC_ENABLE_SYCL
#include "collocation_common.hpp"
#include <gauxc/exceptions/sycl_exception.hpp>
#include <gauxc/util/sycl_util.hpp>
#include "sycl/collocation_device.hpp"



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

    check_collocation( npts, nbf, d.eval.data(), eval.data() );

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

    check_collocation( npts, nbf, d.eval.data(), eval.data() );
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

    check_collocation( tasks[i].npts, tasks[i].nbe, ref_eval, eval.data() );
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

    check_collocation( tasks[i].npts, tasks[i].nbe, ref_eval, eval.data() );
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

    check_collocation( npts, nbf, d.eval.data(), eval.data() );
    check_collocation( npts, nbf, d.deval_x.data(), deval_x.data() );
    check_collocation( npts, nbf, d.deval_y.data(), deval_y.data() );
    check_collocation( npts, nbf, d.deval_z.data(), deval_z.data() );
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

    check_collocation( npts, nbf, d.eval.data(), eval.data() );
    check_collocation( npts, nbf, d.deval_x.data(), deval_x.data() );
    check_collocation( npts, nbf, d.deval_y.data(), deval_y.data() );
    check_collocation( npts, nbf, d.deval_z.data(), deval_z.data() );

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

    auto npts = tasks[i].npts;
    auto nbe  = tasks[i].nbe;
    check_collocation( npts, nbe, ref_eval, eval.data() );
    check_collocation( npts, nbe, ref_deval_x, deval_x.data() );
    check_collocation( npts, nbe, ref_deval_y, deval_y.data() );
    check_collocation( npts, nbe, ref_deval_z, deval_z.data() );
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

    auto npts = tasks[i].npts;
    auto nbe  = tasks[i].nbe;
    check_collocation( npts, nbe, ref_eval, eval.data() );
    check_collocation( npts, nbe, ref_deval_x, deval_x.data() );
    check_collocation( npts, nbe, ref_deval_y, deval_y.data() );
    check_collocation( npts, nbe, ref_deval_z, deval_z.data() );
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
