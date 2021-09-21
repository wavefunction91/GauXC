#ifdef GAUXC_ENABLE_CUDA
#include "collocation_common.hpp"
#include "device/common/collocation_device.hpp"
#include "device_specific/cuda_util.hpp"


auto populate_device_cuda( const BasisSet<double>& basis,
                           const std::vector<ref_collocation_data>& ref_data,
                           bool pop_grad = false ) {

  std::vector< XCDeviceTask > tasks;

  auto shells_device  = util::cuda_malloc<Shell<double>>( basis.size() );
  std::vector<Shell<double>> shells( basis );
  util::cuda_copy( basis.size(), shells_device, shells.data() );

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
    if(pop_grad) {
      task.dbfx       = util::cuda_malloc<double>( nbf * npts );
      task.dbfy       = util::cuda_malloc<double>( nbf * npts );
      task.dbfz       = util::cuda_malloc<double>( nbf * npts );
    }

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

  return std::pair(shells_device,tasks);
}


void cuda_check_collocation( const std::vector<XCDeviceTask>& tasks,
                             const std::vector<ref_collocation_data>& ref_data,
                             bool check_grad = false ) {

  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    std::vector<double> eval (tasks[i].nbe * tasks[i].npts);
    util::cuda_copy( eval.size(), eval.data(), tasks[i].bf );

    check_collocation_transpose( tasks[i].npts, tasks[i].nbe, ref_eval, 
      eval.data(), "IT = " + std::to_string(i) + " BF EVAL" );

    if( check_grad ) {
      auto* ref_deval_x = ref_data[i].deval_x.data();
      auto* ref_deval_y = ref_data[i].deval_y.data();
      auto* ref_deval_z = ref_data[i].deval_z.data();

      std::vector<double> deval_x (tasks[i].nbe * tasks[i].npts);
      std::vector<double> deval_y (tasks[i].nbe * tasks[i].npts);
      std::vector<double> deval_z (tasks[i].nbe * tasks[i].npts);

      util::cuda_copy( eval.size(), deval_x.data(), tasks[i].dbfx );
      util::cuda_copy( eval.size(), deval_y.data(), tasks[i].dbfy );
      util::cuda_copy( eval.size(), deval_z.data(), tasks[i].dbfz );

      auto npts = tasks[i].npts;
      auto nbe  = tasks[i].nbe;
      check_collocation_transpose( npts, nbe, ref_deval_x, deval_x.data(), "IT = " + std::to_string(i) + " BFX EVAL" );
      check_collocation_transpose( npts, nbe, ref_deval_y, deval_y.data(), "IT = " + std::to_string(i) + " BFY EVAL" );
      check_collocation_transpose( npts, nbe, ref_deval_z, deval_z.data(), "IT = " + std::to_string(i) + " BFZ EVAL" );
    }

  }

}


    









void test_cuda_collocation_masked_combined( const BasisSet<double>& basis, std::ifstream& in_file, bool grad ) {



  std::vector<ref_collocation_data> ref_data;

  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }


  type_erased_queue stream( std::make_shared<util::cuda_stream>() );
  auto [shells_device,tasks] = populate_device_cuda( basis, ref_data, grad );


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

  if(grad)
    eval_collocation_masked_combined_deriv1( tasks.size(), npts_max,
      nshells_max, shells_device, tasks_device, stream );
  else
    eval_collocation_masked_combined( tasks.size(), npts_max,
      nshells_max, shells_device, tasks_device, stream );

  util::cuda_device_sync();

  cuda_check_collocation( tasks, ref_data, grad );


  for( auto& t : tasks ) {
    util::cuda_free( t.points, t.shell_offs, t.shell_list, t.bf );
    if(grad) util::cuda_free( t.dbfx, t.dbfy, t.dbfz );
  }
  util::cuda_free( tasks_device, shells_device );
}

void test_cuda_collocation( const BasisSet<double>& basis, 
  std::ifstream& in_file ) {

  test_cuda_collocation_masked_combined( basis, in_file, false );

}
void test_cuda_collocation_deriv1( const BasisSet<double>& basis, 
  std::ifstream& in_file ) {

  test_cuda_collocation_masked_combined( basis, in_file, true );

}
  














void test_cuda_collocation_shell_to_task( const BasisSet<double>& basis, 
  std::ifstream& in_file, bool grad) {

  // Load reference data
  std::vector<ref_collocation_data> ref_data;
  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  //for( auto i = 0; i < ref_data.size(); ++i ) {
  //  std::cout << "IT = " << i << std::endl;
  //  auto& pts = ref_data[i].pts;
  //  for( auto j = 0; j < pts.size(); ++j )
  //    std::cout << "  PT(" << j << ") = "
  //              << pts[j][0] << ", "
  //              << pts[j][1] << ", "
  //              << pts[j][2] << std::endl;
  //}

  // Populate base task information
  type_erased_queue stream( std::make_shared<util::cuda_stream>() );
  auto [shells_device,tasks] = populate_device_cuda( basis, ref_data, grad );

  // Send tasks to device
  auto* tasks_device = util::cuda_malloc<XCDeviceTask>( tasks.size() );
  util::cuda_copy( tasks.size(), tasks_device, tasks.data() );


  // Form Shell -> Task data structures 
  std::vector< std::vector<int32_t> >
    shell_to_task_idx( basis.size() ),
    shell_to_task_off( basis.size() );

  int itask = 0;
  for( auto& d : ref_data ) {
    const auto& mask = d.mask;

    // Reform offsets 
    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();

    // Form shell -> task
    for( auto i = 0; i < mask.size(); ++i ) {
      auto ish = mask[i];
      shell_to_task_idx[ish].emplace_back(itask);
      shell_to_task_off[ish].emplace_back( offs[i] );
    }
    itask++;

  }

  std::vector<ShellToTaskDevice> shell_to_task;
  for( auto ish = 0; ish < basis.size(); ++ish ) {
    shell_to_task.emplace_back();

    const auto ntask = shell_to_task_idx[ish].size();
    shell_to_task.back().ntask = ntask;
    shell_to_task.back().shell_device = shells_device + ish;
    shell_to_task.back().task_idx_device = util::cuda_malloc<int32_t>( ntask );
    shell_to_task.back().task_shell_offs_device =
      util::cuda_malloc<int32_t>( ntask );

    util::cuda_copy( ntask, shell_to_task.back().task_idx_device, 
      shell_to_task_idx[ish].data() );
    util::cuda_copy( ntask, shell_to_task.back().task_shell_offs_device, 
      shell_to_task_off[ish].data() );

  }


  std::vector<uint32_t> shell_idx( basis.size() );
  std::iota( shell_idx.begin(), shell_idx.end(), 0 );

  std::sort( shell_idx.begin(), shell_idx.end(),
    [&]( auto i, auto j ){ return basis.at(i).l() < basis.at(j).l(); } );

  {
  std::vector<ShellToTaskDevice> shell_to_task_sorted( basis.size() );
  for( auto i = 0; i < basis.size(); ++i ) 
    shell_to_task_sorted[i] = shell_to_task[shell_idx[i]];
  shell_to_task = std::move(shell_to_task_sorted);
  }


  auto max_l = std::max_element(basis.begin(),basis.end(),
    [](const auto&a, const auto& b){ return a.l() < b.l(); } )->l();
  std::vector<uint32_t> l_counts( max_l+1 );
  std::vector<bool>     l_pure( max_l+1);
  for( auto l = 0; l <= max_l; ++l ) {
    l_counts[l] = std::count_if( basis.begin(), basis.end(),
      [&](auto& s){ return s.l() == l; } );
    auto f = std::find_if( basis.begin(), basis.end(),
      [&](auto& s){ return s.l() == l; } );
    l_pure[l] = f->pure();
  }

  auto* shell_to_task_device = util::cuda_malloc<ShellToTaskDevice>(basis.size());
  util::cuda_copy( basis.size(), shell_to_task_device, shell_to_task.data() );
  util::cuda_device_sync();


#if 1
  auto* p = shell_to_task_device;
  for( auto l = 0; l <= max_l; ++l ) {
    if( grad ) eval_collocation_shell_to_task_gradient( l_counts[l], l, l_pure[l], p, tasks_device, stream );
    else       eval_collocation_shell_to_task         ( l_counts[l], l, l_pure[l], p ,tasks_device, stream );
    p += l_counts[l];
  }
#else
  for( auto i = 0; i < basis.size(); ++i ) {
    auto l = basis.at(shell_idx[i]).l();
    auto pure = basis.at(shell_idx[i]).pure();
    eval_collocation_shell_to_task( 1, l, pure, shell_to_task_device + i ,tasks_device, stream );
  }
#endif



  util::cuda_device_sync();
  cuda_check_collocation( tasks, ref_data, grad );

      
  for( auto& t : tasks ) {
    util::cuda_free( t.points, t.shell_offs, t.shell_list, t.bf );
    if(grad) util::cuda_free( t.dbfx, t.dbfy, t.dbfz );
  }
  util::cuda_free( tasks_device, shells_device, shell_to_task_device );
  for( auto& s : shell_to_task ) {
    util::cuda_free( s.task_idx_device, s.task_shell_offs_device );
  }
}



void test_cuda_collocation_shell_to_task( const BasisSet<double>& basis, 
  std::ifstream& in_file) {

  test_cuda_collocation_shell_to_task(basis,in_file,false);

}
void test_cuda_collocation_shell_to_task_gradient( const BasisSet<double>& basis, 
  std::ifstream& in_file) {

  test_cuda_collocation_shell_to_task(basis,in_file,true);

}





















#endif // GAUXC_ENABLE_CUDA

