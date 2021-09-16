#ifdef GAUXC_ENABLE_CUDA
#include "collocation_common.hpp"
#include "device/common/collocation_device.hpp"
#include <gauxc/util/cuda_util.hpp>


void test_cuda_collocation( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  std::vector< XCDeviceTask > tasks;

  auto shells_device  = util::cuda_malloc<Shell<double>>( basis.size() );
  std::vector<Shell<double>> shells( basis );
  util::cuda_copy( basis.size(), shells_device, shells.data() );

  type_erased_queue stream( std::make_shared<util::cuda_stream>() );
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

  auto* tasks_device = util::cuda_malloc<XCDeviceTask>( tasks.size() );
  util::cuda_copy( tasks.size(), tasks_device, tasks.data() );

  eval_collocation_masked_combined( tasks.size(), npts_max,
    nshells_max, shells_device, tasks_device, stream );

  util::cuda_device_sync();


  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    util::cuda_copy( eval.size(), eval.data(), tasks[i].bf );

    check_collocation_transpose( tasks[i].npts, tasks[i].nbe, ref_eval, eval.data() );
  }


  for( auto& t : tasks ) {
    util::cuda_free( t.points, t.shell_offs, t.shell_list, t.bf );
  }
  util::cuda_free( tasks_device, shells_device );
}

void test_cuda_collocation_deriv1( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  std::vector< XCDeviceTask > tasks;

  auto shells_device  = util::cuda_malloc<Shell<double>>( basis.size() );
  std::vector<Shell<double>> shells( basis );
  util::cuda_copy( basis.size(), shells_device, shells.data() );

  type_erased_queue stream( std::make_shared<util::cuda_stream>() );
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

  auto* tasks_device = util::cuda_malloc<XCDeviceTask>( tasks.size() );
  util::cuda_copy( tasks.size(), tasks_device, tasks.data() );

  eval_collocation_masked_combined_deriv1( tasks.size(), npts_max,
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


    auto npts = tasks[i].npts;
    auto nbe  = tasks[i].nbe;
    check_collocation_transpose( npts, nbe, ref_eval, eval.data() );
    check_collocation_transpose( npts, nbe, ref_deval_x, deval_x.data() );
    check_collocation_transpose( npts, nbe, ref_deval_y, deval_y.data() );
    check_collocation_transpose( npts, nbe, ref_deval_z, deval_z.data() );
  }


  for( auto& t : tasks ) {
    util::cuda_free( t.points, t.shell_offs, t.shell_list, t.bf, t.dbfx, t.dbfy,
      t.dbfz );
  }
  util::cuda_free( tasks_device, shells_device );
}
























void test_cuda_collocation_shell_to_task( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  //ref_data.erase( ref_data.begin()+1, ref_data.end() );

  std::vector< XCDeviceTask > tasks;

  auto shells_device  = util::cuda_malloc<Shell<double>>( basis.size() );
  std::vector<Shell<double>> shells( basis );
  util::cuda_copy( basis.size(), shells_device, shells.data() );

  std::vector< std::vector<int32_t> >
    shell_to_task_idx( basis.size() ),
    shell_to_task_off( basis.size() );
  std::vector< int32_t >
    shell_to_task_ntask( basis.size() ),
    shell_to_task_ls( basis.size() ),
    shell_to_task_pure( basis.size() );

  for( auto i = 0; i < basis.size(); ++i ) {
    shell_to_task_ls[i] = basis[i].l();
    shell_to_task_pure[i] = basis[i].pure();
  }


  type_erased_queue stream( std::make_shared<util::cuda_stream>() );
  int itask = 0;
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

    // Form shell -> task
    for( auto i = 0; i < mask.size(); ++i ) {
      auto ish = mask[i];
      shell_to_task_idx[ish].emplace_back(itask);
      shell_to_task_off[ish].emplace_back( offs[i] );
    }
    itask++;

  }

  
  auto* tasks_device = util::cuda_malloc<XCDeviceTask>( tasks.size() );
  util::cuda_copy( tasks.size(), tasks_device, tasks.data() );

  std::vector<int32_t*> shell_to_task_idx_device, shell_to_task_off_device;
  for( auto i = 0; i < basis.size(); ++i ) {
    const auto ntask = shell_to_task_idx[i].size();
    shell_to_task_ntask[i] = ntask;
    shell_to_task_idx_device.emplace_back( util::cuda_malloc<int32_t>( ntask ) );
    shell_to_task_off_device.emplace_back( util::cuda_malloc<int32_t>( ntask ) );

    util::cuda_copy( ntask, shell_to_task_idx_device.back(), shell_to_task_idx[i].data() );
    util::cuda_copy( ntask, shell_to_task_off_device.back(), shell_to_task_off[i].data() );
  }

  eval_collocation_shell_to_task( basis.size(), shells_device, 
    shell_to_task_idx_device.data(), shell_to_task_off_device.data(),
    shell_to_task_ntask.data(), shell_to_task_ls.data(), shell_to_task_pure.data(),
    tasks_device, stream );


  util::cuda_device_sync();

  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    util::cuda_copy( eval.size(), eval.data(), tasks[i].bf );

    check_collocation_transpose( tasks[i].npts, tasks[i].nbe, ref_eval, eval.data() );
  }
      
  for( auto& t : tasks ) {
    util::cuda_free( t.points, t.shell_offs, t.shell_list, t.bf );
  }
  util::cuda_free( tasks_device, shells_device );
  for( auto& p : shell_to_task_idx_device ) util::cuda_free( p );
  for( auto& p : shell_to_task_off_device ) util::cuda_free( p );
}



void test_cuda_collocation_shell_to_task_gradient( const BasisSet<double>& basis, std::ifstream& in_file) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  //ref_data.erase( ref_data.begin()+1, ref_data.end() );

  std::vector< XCDeviceTask > tasks;

  auto shells_device  = util::cuda_malloc<Shell<double>>( basis.size() );
  std::vector<Shell<double>> shells( basis );
  util::cuda_copy( basis.size(), shells_device, shells.data() );

  std::vector< std::vector<int32_t> >
    shell_to_task_idx( basis.size() ),
    shell_to_task_off( basis.size() );
  std::vector< int32_t >
    shell_to_task_ntask( basis.size() ),
    shell_to_task_ls( basis.size() ),
    shell_to_task_pure( basis.size() );

  for( auto i = 0; i < basis.size(); ++i ) {
    shell_to_task_ls[i] = basis[i].l();
    shell_to_task_pure[i] = basis[i].pure();
  }

  //std::cout << "NBF = " << basis.nbf() << std::endl;
  type_erased_queue stream( std::make_shared<util::cuda_stream>() );
  int itask = 0;
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

    // Form shell -> task
    for( auto i = 0; i < mask.size(); ++i ) {
      auto ish = mask[i];
      shell_to_task_idx[ish].emplace_back(itask);
      shell_to_task_off[ish].emplace_back( offs[i] );
    }
    itask++;

  }

  
  auto* tasks_device = util::cuda_malloc<XCDeviceTask>( tasks.size() );
  util::cuda_copy( tasks.size(), tasks_device, tasks.data() );

  std::vector<int32_t*> shell_to_task_idx_device, shell_to_task_off_device;
  for( auto i = 0; i < basis.size(); ++i ) {
    const auto ntask = shell_to_task_idx[i].size();
    shell_to_task_ntask[i] = ntask;
    shell_to_task_idx_device.emplace_back( util::cuda_malloc<int32_t>( ntask ) );
    shell_to_task_off_device.emplace_back( util::cuda_malloc<int32_t>( ntask ) );

    util::cuda_copy( ntask, shell_to_task_idx_device.back(), shell_to_task_idx[i].data() );
    util::cuda_copy( ntask, shell_to_task_off_device.back(), shell_to_task_off[i].data() );
  }

  eval_collocation_shell_to_task_gradient( basis.size(), shells_device, 
    shell_to_task_idx_device.data(), shell_to_task_off_device.data(),
    shell_to_task_ntask.data(), shell_to_task_ls.data(), shell_to_task_pure.data(),
    tasks_device, stream );


  util::cuda_device_sync();

#if 1
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


    auto npts = tasks[i].npts;
    auto nbe  = tasks[i].nbe;
    check_collocation_transpose( npts, nbe, ref_eval, eval.data() );
    check_collocation_transpose( npts, nbe, ref_deval_x, deval_x.data() );
    check_collocation_transpose( npts, nbe, ref_deval_y, deval_y.data() );
    check_collocation_transpose( npts, nbe, ref_deval_z, deval_z.data() );
  }
#endif
      
  for( auto& t : tasks ) {
    util::cuda_free( t.points, t.shell_offs, t.shell_list, t.bf, t.dbfx, t.dbfy,
      t.dbfz );
  }
  util::cuda_free( tasks_device, shells_device );
  for( auto& p : shell_to_task_idx_device ) util::cuda_free( p );
  for( auto& p : shell_to_task_off_device ) util::cuda_free( p );
}
#endif // GAUXC_ENABLE_SYCL

