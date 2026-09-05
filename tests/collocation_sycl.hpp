/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#ifdef GAUXC_HAS_SYCL
#include "collocation_common.hpp"
#include "hdf5_test_serialization.hpp"
#include "hdf5_test_serialization_impl.hpp"
#include "device/common/collocation_device.hpp"
#include "device_specific/sycl_util.hpp"
#include <gauxc/basisset_map.hpp>


auto populate_device_sycl( const BasisSet<double>& basis,
                           const std::vector<ref_collocation_data>& ref_data,
                           bool pop_grad, bool pop_hess, bool pop_lapl, bool pop_lapl_grad,
                           ::sycl::queue& q ) {

  std::vector< XCDeviceTask > tasks;

  auto shells_device  = util::sycl_malloc<Shell<double>>( basis.size(), q );
  std::vector<Shell<double>> shells( basis );
  util::sycl_copy( basis.size(), shells_device, shells.data(), q );

  for( auto& d : ref_data ) {
    const auto npts = d.pts.size();
    const auto nbf  = d.eval.size() / npts;

    const auto& mask = d.mask;
    const auto& pts  = d.pts;

    /// XXX: THIS DOES NOT POPULATE A VALID TASK, ONLY WHAT's REQUIRED FOR THIS
    //  TEST
    auto& task = tasks.emplace_back();
    task.npts    = npts;
    task.bfn_screening.nbe     = nbf;
    task.bfn_screening.nshells = mask.size();
    task.bfn_screening.shell_offs = util::sycl_malloc<size_t>( mask.size(), q );
    task.bfn_screening.shell_list = util::sycl_malloc<size_t>( mask.size(), q );

    task.points_x     = util::sycl_malloc<double>( npts, q );
    task.points_y     = util::sycl_malloc<double>( npts, q );
    task.points_z     = util::sycl_malloc<double>( npts, q );
    task.bf         = util::sycl_malloc<double>( nbf * npts, q );
    if(pop_grad) {
      task.dbfx = util::sycl_malloc<double>( nbf * npts, q );
      task.dbfy = util::sycl_malloc<double>( nbf * npts, q );
      task.dbfz = util::sycl_malloc<double>( nbf * npts, q );
    }

    if(pop_hess) {
      task.d2bfxx = util::sycl_malloc<double>( nbf * npts, q );
      task.d2bfxy = util::sycl_malloc<double>( nbf * npts, q );
      task.d2bfxz = util::sycl_malloc<double>( nbf * npts, q );
      task.d2bfyy = util::sycl_malloc<double>( nbf * npts, q );
      task.d2bfyz = util::sycl_malloc<double>( nbf * npts, q );
      task.d2bfzz = util::sycl_malloc<double>( nbf * npts, q );
    }

    if(pop_lapl) {
      task.d2bflapl = util::sycl_malloc<double>( nbf * npts, q );
    }

    if(pop_lapl_grad) {
      task.d3bflapl_x = util::sycl_malloc<double>( nbf * npts, q );
      task.d3bflapl_y = util::sycl_malloc<double>( nbf * npts, q );
      task.d3bflapl_z = util::sycl_malloc<double>( nbf * npts, q );
    }

    auto* pts_x_device = task.points_x;
    auto* pts_y_device = task.points_y;
    auto* pts_z_device = task.points_z;
    auto* offs_device = task.bfn_screening.shell_offs;
    auto* mask_device = task.bfn_screening.shell_list;


    std::vector<double> pts_x, pts_y, pts_z;
    for( auto pt : pts ) {
      pts_x.emplace_back(pt[0]);
      pts_y.emplace_back(pt[1]);
      pts_z.emplace_back(pt[2]);
    }
    util::sycl_copy( npts, pts_x_device, pts_x.data(), q );
    util::sycl_copy( npts, pts_y_device, pts_y.data(), q );
    util::sycl_copy( npts, pts_z_device, pts_z.data(), q );

    std::vector<size_t> mask_ul( mask.size() );
    std::copy( mask.begin(), mask.end(), mask_ul.begin() );
    util::sycl_copy( mask.size(), mask_device, mask_ul.data(), q );

    std::vector<size_t> offs( mask.size() );
    offs[0] = 0;
    for( int i = 1; i < mask.size(); ++i )
      offs[i] = offs[i-1] + basis[mask[i-1]].size();
    util::sycl_copy( offs.size(), offs_device, offs.data(), q );

  }

  return std::pair(shells_device,tasks);
}


void sycl_check_collocation( const std::vector<XCDeviceTask>& tasks,
                             const std::vector<ref_collocation_data>& ref_data,
                             bool check_grad, bool check_hess, bool check_lapl, bool check_lapl_grad,
                             ::sycl::queue& q) {

  for( int i = 0; i < tasks.size(); i++ ) {

    auto* ref_eval = ref_data[i].eval.data();
    std::vector<double> eval (tasks[i].bfn_screening.nbe * tasks[i].npts);
    util::sycl_copy( eval.size(), eval.data(), tasks[i].bf, q );

    check_collocation_transpose( tasks[i].npts, tasks[i].bfn_screening.nbe, ref_eval,
      eval.data(), "IT = " + std::to_string(i) + " BF EVAL" );

    if( check_grad ) {
      auto* ref_deval_x = ref_data[i].deval_x.data();
      auto* ref_deval_y = ref_data[i].deval_y.data();
      auto* ref_deval_z = ref_data[i].deval_z.data();

      std::vector<double> deval_x (tasks[i].bfn_screening.nbe * tasks[i].npts);
      std::vector<double> deval_y (tasks[i].bfn_screening.nbe * tasks[i].npts);
      std::vector<double> deval_z (tasks[i].bfn_screening.nbe * tasks[i].npts);

      util::sycl_copy( eval.size(), deval_x.data(), tasks[i].dbfx, q );
      util::sycl_copy( eval.size(), deval_y.data(), tasks[i].dbfy, q );
      util::sycl_copy( eval.size(), deval_z.data(), tasks[i].dbfz, q );

      auto npts = tasks[i].npts;
      auto nbe  = tasks[i].bfn_screening.nbe;
      check_collocation_transpose( npts, nbe, ref_deval_x, deval_x.data(), "IT = " + std::to_string(i) + " BFX EVAL" );
      check_collocation_transpose( npts, nbe, ref_deval_y, deval_y.data(), "IT = " + std::to_string(i) + " BFY EVAL" );
      check_collocation_transpose( npts, nbe, ref_deval_z, deval_z.data(), "IT = " + std::to_string(i) + " BFZ EVAL" );
    }

    if( check_hess ) {
      auto* ref_d2eval_xx = ref_data[i].d2eval_xx.data();
      auto* ref_d2eval_xy = ref_data[i].d2eval_xy.data();
      auto* ref_d2eval_xz = ref_data[i].d2eval_xz.data();
      auto* ref_d2eval_yy = ref_data[i].d2eval_yy.data();
      auto* ref_d2eval_yz = ref_data[i].d2eval_yz.data();
      auto* ref_d2eval_zz = ref_data[i].d2eval_zz.data();

      std::vector<double> d2eval_xx (tasks[i].bfn_screening.nbe * tasks[i].npts);
      std::vector<double> d2eval_xy (tasks[i].bfn_screening.nbe * tasks[i].npts);
      std::vector<double> d2eval_xz (tasks[i].bfn_screening.nbe * tasks[i].npts);
      std::vector<double> d2eval_yy (tasks[i].bfn_screening.nbe * tasks[i].npts);
      std::vector<double> d2eval_yz (tasks[i].bfn_screening.nbe * tasks[i].npts);
      std::vector<double> d2eval_zz (tasks[i].bfn_screening.nbe * tasks[i].npts);

      util::sycl_copy( eval.size(), d2eval_xx.data(), tasks[i].d2bfxx, q );
      util::sycl_copy( eval.size(), d2eval_xy.data(), tasks[i].d2bfxy, q );
      util::sycl_copy( eval.size(), d2eval_xz.data(), tasks[i].d2bfxz, q );
      util::sycl_copy( eval.size(), d2eval_yy.data(), tasks[i].d2bfyy, q );
      util::sycl_copy( eval.size(), d2eval_yz.data(), tasks[i].d2bfyz, q );
      util::sycl_copy( eval.size(), d2eval_zz.data(), tasks[i].d2bfzz, q );

      auto npts = tasks[i].npts;
      auto nbe  = tasks[i].bfn_screening.nbe;
      check_collocation_transpose( npts, nbe, ref_d2eval_xx, d2eval_xx.data(), "IT = " + std::to_string(i) + " BFXX EVAL" );
      check_collocation_transpose( npts, nbe, ref_d2eval_xy, d2eval_xy.data(), "IT = " + std::to_string(i) + " BFXY EVAL" );
      check_collocation_transpose( npts, nbe, ref_d2eval_xz, d2eval_xz.data(), "IT = " + std::to_string(i) + " BFXZ EVAL" );
      check_collocation_transpose( npts, nbe, ref_d2eval_yy, d2eval_yy.data(), "IT = " + std::to_string(i) + " BFYY EVAL" );
      check_collocation_transpose( npts, nbe, ref_d2eval_yz, d2eval_yz.data(), "IT = " + std::to_string(i) + " BFYZ EVAL" );
      check_collocation_transpose( npts, nbe, ref_d2eval_zz, d2eval_zz.data(), "IT = " + std::to_string(i) + " BFZZ EVAL" );
    }

    if( check_lapl ) {
      auto npts = tasks[i].npts;
      auto nbe  = tasks[i].bfn_screening.nbe;
      auto* ref_d2eval_lapl = ref_data[i].d2eval_lapl.data();
      std::vector<double> d2eval_lapl(npts * nbe);
      util::sycl_copy(eval.size(), d2eval_lapl.data(), tasks[i].d2bflapl, q);
      check_collocation_transpose(npts, nbe, ref_d2eval_lapl, d2eval_lapl.data(), "IT = " + std::to_string(i) + "BFLAPL EVAL" );
    }

    if( check_lapl_grad ) {
      auto npts = tasks[i].npts;
      auto nbe  = tasks[i].bfn_screening.nbe;
      auto* ref_d3eval_lapl_x = ref_data[i].d3eval_lapl_x.data();
      auto* ref_d3eval_lapl_y = ref_data[i].d3eval_lapl_y.data();
      auto* ref_d3eval_lapl_z = ref_data[i].d3eval_lapl_z.data();
      std::vector<double> d3eval_lapl_x(npts * nbe);
      std::vector<double> d3eval_lapl_y(npts * nbe);
      std::vector<double> d3eval_lapl_z(npts * nbe);
      util::sycl_copy(eval.size(), d3eval_lapl_x.data(), tasks[i].d3bflapl_x, q);
      util::sycl_copy(eval.size(), d3eval_lapl_y.data(), tasks[i].d3bflapl_y, q);
      util::sycl_copy(eval.size(), d3eval_lapl_z.data(), tasks[i].d3bflapl_z, q);
      check_collocation_transpose(npts, nbe, ref_d3eval_lapl_x, d3eval_lapl_x.data(), "IT = " + std::to_string(i) + "BFLAPL_X EVAL" );
      check_collocation_transpose(npts, nbe, ref_d3eval_lapl_y, d3eval_lapl_y.data(), "IT = " + std::to_string(i) + "BFLAPL_Y EVAL" );
      check_collocation_transpose(npts, nbe, ref_d3eval_lapl_z, d3eval_lapl_z.data(), "IT = " + std::to_string(i) + "BFLAPL_Z EVAL" );
    }

  }

}


void test_sycl_collocation_masked_combined( const BasisSet<double>& basis, const std::string& filename, bool grad ) {



  std::vector<ref_collocation_data> ref_data;
  read_collocation_data(ref_data, filename);


  device_queue stream( std::make_shared<util::sycl_queue>() );
  // USM alloc/copy are context-bound, so the raw sycl::queue& (as opposed
  // to the type-erased device_queue) is needed for the host-side setup below.
  ::sycl::queue& q = stream.queue_as<util::sycl_queue>().queue;
  auto [shells_device,tasks] = populate_device_sycl( basis, ref_data, grad, false, false, false, q );


  const auto nshells_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.bfn_screening.nshells < b.bfn_screening.nshells;
    })->bfn_screening.nshells;

  const auto npts_max = std::max_element( tasks.begin(), tasks.end(),
    []( const auto& a, const auto& b ) {
      return a.npts < b.npts;
    })->npts;

  auto* tasks_device = util::sycl_malloc<XCDeviceTask>( tasks.size(), q );
  util::sycl_copy( tasks.size(), tasks_device, tasks.data(), q );

  if(grad)
    eval_collocation_masked_combined_deriv1( tasks.size(), npts_max,
      nshells_max, shells_device, tasks_device, stream );
  else
    eval_collocation_masked_combined( tasks.size(), npts_max,
      nshells_max, shells_device, tasks_device, stream );

  util::sycl_device_sync(q);

  sycl_check_collocation( tasks, ref_data, grad, false, false, false, q );


  for( auto& t : tasks ) {
    util::sycl_free( q, t.points_x, t.points_y, t.points_z, t.bfn_screening.shell_offs, t.bfn_screening.shell_list, t.bf );
    if(grad) util::sycl_free( q, t.dbfx, t.dbfy, t.dbfz );
  }
  util::sycl_free( q, tasks_device, shells_device );
}

void test_sycl_collocation( const BasisSet<double>& basis,
  const std::string& filename ) {

  test_sycl_collocation_masked_combined( basis, filename, false );

}

void test_sycl_collocation_deriv1( const BasisSet<double>& basis,
  const std::string& filename ) {

  test_sycl_collocation_masked_combined( basis, filename, true );

}


void test_sycl_collocation_shell_to_task( const BasisSet<double>& basis,  const BasisSetMap& basis_map,
  const std::string& filename, bool grad, bool hess, bool lapl, bool lapl_grad) {

  // Load reference data
  std::vector<ref_collocation_data> ref_data;
  read_collocation_data(ref_data, filename);

  // Populate base task information
  device_queue stream( std::make_shared<util::sycl_queue>() );
  ::sycl::queue& q = stream.queue_as<util::sycl_queue>().queue;
  auto [shells_device,tasks] = populate_device_sycl( basis, ref_data, grad, hess, lapl, lapl_grad, q );

  // Send tasks to device
  auto* tasks_device = util::sycl_malloc<XCDeviceTask>( tasks.size(), q );
  util::sycl_copy( tasks.size(), tasks_device, tasks.data(), q );


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
    shell_to_task.back().task_idx_device = util::sycl_malloc<int32_t>( ntask, q );
    shell_to_task.back().task_shell_offs_device =
      util::sycl_malloc<int32_t>( ntask, q );

    util::sycl_copy( ntask, shell_to_task.back().task_idx_device,
      shell_to_task_idx[ish].data(), q );
    util::sycl_copy( ntask, shell_to_task.back().task_shell_offs_device,
      shell_to_task_off[ish].data(), q );

  }


  // Sort shells by L
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


  // Send Shell -> Task to device
  auto* shell_to_task_device = util::sycl_malloc<ShellToTaskDevice>(basis.size(), q);
  util::sycl_copy( basis.size(), shell_to_task_device, shell_to_task.data(), q );
  util::sycl_device_sync(q);

  // Form angular momentum batches for collocation eval
  auto max_l = std::max_element(basis.begin(),basis.end(),
    [](const auto&a, const auto& b){ return a.l() < b.l(); } )->l();
  std::vector<AngularMomentumShellToTaskBatch> l_batched_shell_to_task(max_l+1);
  {
  auto* p = shell_to_task_device;
  auto* h = shell_to_task.data();
  for( auto l = 0; l <= max_l; ++l ) {
    auto nsh = basis_map.nshells_with_l(l);
    auto pure = basis_map.l_purity(l);
    l_batched_shell_to_task[l].nshells_in_batch     = nsh;
    l_batched_shell_to_task[l].pure                 = pure;
    l_batched_shell_to_task[l].shell_to_task_device = p;

    size_t total_ntask = std::accumulate( h, h + nsh, 0ul,
      [](auto& a, auto& b){ return a + b.ntask; } );
    l_batched_shell_to_task[l].ntask_average = total_ntask / nsh;

    p += nsh;
    h += nsh;
  }
  }


  if( lapl_grad )
    eval_collocation_shell_to_task_lapgrad( max_l, l_batched_shell_to_task.data(),
      tasks_device, stream );
  else if( hess )
    eval_collocation_shell_to_task_hessian( max_l, l_batched_shell_to_task.data(),
      tasks_device, stream );
  else if( lapl )
    eval_collocation_shell_to_task_laplacian( max_l, l_batched_shell_to_task.data(),
      tasks_device, stream );
  else if( grad )
    eval_collocation_shell_to_task_gradient( max_l, l_batched_shell_to_task.data(),
      tasks_device, stream );
  else
    eval_collocation_shell_to_task( max_l, l_batched_shell_to_task.data(),
      tasks_device, stream );



  util::sycl_device_sync(q);
  sycl_check_collocation( tasks, ref_data, grad, hess, lapl, lapl_grad, q );


  for( auto& t : tasks ) {
    util::sycl_free( q, t.points_x, t.points_y, t.points_z, t.bfn_screening.shell_offs, t.bfn_screening.shell_list, t.bf );
    if(grad) util::sycl_free( q, t.dbfx, t.dbfy, t.dbfz );
    if(hess) util::sycl_free( q, t.d2bfxx, t.d2bfxy, t.d2bfxz, t.d2bfyy, t.d2bfyz, t.d2bfzz );
    if(lapl) util::sycl_free( q, t.d2bflapl );
    if(lapl_grad) util::sycl_free( q, t.d3bflapl_x, t.d3bflapl_y, t.d3bflapl_z );
  }
  util::sycl_free( q, tasks_device, shells_device, shell_to_task_device );
  for( auto& s : shell_to_task ) {
    util::sycl_free( q, s.task_idx_device, s.task_shell_offs_device );
  }
}



void test_sycl_collocation_shell_to_task( const BasisSet<double>& basis,
  const BasisSetMap& basis_map, const std::string& filename) {

  test_sycl_collocation_shell_to_task(basis,basis_map,filename,false, false, false, false);

}
void test_sycl_collocation_shell_to_task_gradient( const BasisSet<double>& basis,
  const BasisSetMap& basis_map, const std::string& filename) {

  test_sycl_collocation_shell_to_task(basis,basis_map,filename,true, false, false, false);

}
void test_sycl_collocation_shell_to_task_hessian( const BasisSet<double>& basis,
  const BasisSetMap& basis_map, const std::string& filename) {

  test_sycl_collocation_shell_to_task(basis,basis_map,filename,true, true, false, false);

}

void test_sycl_collocation_shell_to_task_laplacian( const BasisSet<double>& basis,
  const BasisSetMap& basis_map, const std::string& filename) {

  test_sycl_collocation_shell_to_task(basis,basis_map,filename,true, false, true, false);

}

void test_sycl_collocation_shell_to_task_lapgrad( const BasisSet<double>& basis,
  const BasisSetMap& basis_map, const std::string& filename) {

  test_sycl_collocation_shell_to_task(basis,basis_map,filename,true, true, true, true);

}

#endif // GAUXC_HAS_SYCL
