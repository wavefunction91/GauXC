#ifdef GAUXC_ENABLE_HIP
#include "collocation_common.hpp"
#include "device/common/collocation_device.hpp"
#include "device_specific/hip_util.hpp"


void test_hip_collocation( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  std::vector< XCDeviceTask > tasks;

  auto shells_device  = util::hip_malloc<Shell<double>>( basis.size() );
  std::vector<Shell<double>> shells( basis );
  util::hip_copy( basis.size(), shells_device, shells.data() );

  type_erased_queue stream( std::make_shared<util::hip_stream>() );
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

  auto* tasks_device = util::hip_malloc<XCDeviceTask>( tasks.size() );
  util::hip_copy( tasks.size(), tasks_device, tasks.data() );

  eval_collocation_masked_combined( tasks.size(), npts_max,
    nshells_max, shells_device, tasks_device, stream );

  util::hip_device_sync();


  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    util::hip_copy( eval.size(), eval.data(), tasks[i].bf );

    check_collocation_transpose( tasks[i].npts, tasks[i].nbe, ref_eval, eval.data() );
  }


  for( auto& t : tasks ) {
    util::hip_free( t.points, t.shell_offs, t.shell_list, t.bf );
  }
  util::hip_free( tasks_device, shells_device );
}

void test_hip_collocation_deriv1( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  std::vector< XCDeviceTask > tasks;

  auto shells_device  = util::hip_malloc<Shell<double>>( basis.size() );
  std::vector<Shell<double>> shells( basis );
  util::hip_copy( basis.size(), shells_device, shells.data() );

  type_erased_queue stream( std::make_shared<util::hip_stream>() );
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

  auto* tasks_device = util::hip_malloc<XCDeviceTask>( tasks.size() );
  util::hip_copy( tasks.size(), tasks_device, tasks.data() );

  eval_collocation_masked_combined_deriv1( tasks.size(), npts_max,
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


    auto npts = tasks[i].npts;
    auto nbe  = tasks[i].nbe;
    check_collocation_transpose( npts, nbe, ref_eval, eval.data() );
    check_collocation_transpose( npts, nbe, ref_deval_x, deval_x.data() );
    check_collocation_transpose( npts, nbe, ref_deval_y, deval_y.data() );
    check_collocation_transpose( npts, nbe, ref_deval_z, deval_z.data() );
  }


  for( auto& t : tasks ) {
    util::hip_free( t.points, t.shell_offs, t.shell_list, t.bf, t.dbfx, t.dbfy,
      t.dbfz );
  }
  util::hip_free( tasks_device, shells_device );
}
#endif // GAUXC_ENABLE_HIP

