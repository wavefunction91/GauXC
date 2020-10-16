#ifdef GAUXC_ENABLE_HIP
#include "collocation_common.hpp"
#include <gauxc/exceptions/hip_exception.hpp>
#include <gauxc/util/hip_util.hpp>
#include "hip/collocation_device.hpp"




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

    check_collocation_transpose( npts, nbf, d.eval.data(), eval.data() );

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

    check_collocation_transpose( npts, nbf, d.eval.data(), eval.data() );

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

    check_collocation_transpose( tasks[i].npts, tasks[i].nbe, ref_eval, eval.data() );
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

    check_collocation_transpose( tasks[i].npts, tasks[i].nbe, ref_eval, eval.data() );
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

    check_collocation_transpose( npts, nbf, d.eval.data(), eval.data() );
    check_collocation_transpose( npts, nbf, d.deval_x.data(), deval_x.data() );
    check_collocation_transpose( npts, nbf, d.deval_y.data(), deval_y.data() );
    check_collocation_transpose( npts, nbf, d.deval_z.data(), deval_z.data() );

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

    check_collocation_transpose( npts, nbf, d.eval.data(), eval.data() );
    check_collocation_transpose( npts, nbf, d.deval_x.data(), deval_x.data() );
    check_collocation_transpose( npts, nbf, d.deval_y.data(), deval_y.data() );
    check_collocation_transpose( npts, nbf, d.deval_z.data(), deval_z.data() );

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

    auto npts = tasks[i].npts;
    auto nbe  = tasks[i].nbe;
    check_collocation_transpose( npts, nbe, ref_eval, eval.data() );
    check_collocation_transpose( npts, nbe, ref_deval_x, deval_x.data() );
    check_collocation_transpose( npts, nbe, ref_deval_y, deval_y.data() );
    check_collocation_transpose( npts, nbe, ref_deval_z, deval_z.data() );
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
#endif

