#include <gauxc/xc_integrator/xc_sycl_data.hpp>
#include "buffer_adaptor.hpp"
#include "integrator_common.hpp"

namespace GauXC {

auto gauxc_asynchandler = [] (cl::sycl::exception_list exceptions) {
    for (std::exception_ptr const& e : exceptions) {
        try {
            std::rethrow_exception(e);
        } catch (cl::sycl::exception const& ex) {
            std::cout << "Caught asynchronous SYCL exception:" << std::endl
            << ex.what() << std::endl;
        }
    }
};

template <typename F>
XCSyclData<F>::XCSyclData( size_t _natoms,
                           size_t _n_deriv,
                           size_t _nbf,
                           size_t _nshells,
                           bool _denpack_host,
                           bool _vxcinc_host,
                           bool _batch_l3_blas ):
  nshells(_nshells),
  nbf(_nbf),
  n_deriv(_n_deriv),
  natoms(_natoms),
  denpack_host(_denpack_host),
  vxcinc_host(_vxcinc_host),
  batch_l3_blas(_batch_l3_blas)
{

  // Create SYCL queue
  master_queue.reset( new cl::sycl::queue(cl::sycl::gpu_selector{},
                                          gauxc_asynchandler,
                                          cl::sycl::property_list{cl::sycl::property::queue::in_order{}}) );

  // Allocate up to fill_fraction
  size_t fill_sz = (*master_queue).get_device().get_info<cl::sycl::info::device::max_mem_alloc_size>();
  device_ptr = (void *)cl::sycl::malloc_device(fill_sz, *master_queue);

  //std::cout << "NS = " << nshells << ", NA = " << natoms << ", NBF = " << nbf << std::endl;

  //std::cout << "XCDeviceData has allocated " << fill_sz << " bytes of data"
  //          << std::endl;
  // Allocate static memory with proper alignment
  buffer_adaptor mem( device_ptr, fill_sz );

  shells_device     = mem.aligned_alloc<Shell<F>>( nshells );
  exc_device        = mem.aligned_alloc<F>( 1 );
  nel_device        = mem.aligned_alloc<F>( 1 );
  acc_scr_device    = mem.aligned_alloc<F>( 1 );
  rab_device        = mem.aligned_alloc<F>( natoms * natoms );
  coords_device     = mem.aligned_alloc<F>( 3 * natoms );

  if( not vxcinc_host )
    vxc_device = mem.aligned_alloc<F>( nbf * nbf );
  if( not denpack_host )
    dmat_device = mem.aligned_alloc<F>( nbf * nbf );

  // Get current stack location
  dynmem_ptr = mem.stack();
  dynmem_sz  = mem.nleft();

}



template <typename F>
XCSyclData<F>::~XCSyclData() noexcept {
  if( device_ptr ) util::sycl_free( device_ptr, *master_queue );
}


using task_iterator = std::vector< XCTask >::iterator;
template <typename F>
using device_task_container = std::vector< sycl::XCTaskDevice<F> >;

template <typename F>
std::tuple< task_iterator, device_task_container<F> >
  XCSyclData<F>::generate_buffers( const BasisSet<F>& basis,
                                   task_iterator      task_begin,
                                   task_iterator      task_end    ) {

  // Host data packing arrays
  std::vector< std::array<F,3> > points_pack;
  std::vector< F > weights_pack;
  std::vector< size_t > shell_list_pack;
  std::vector< size_t > shell_offs_pack;
  std::vector< std::array<int64_t, 3> > submat_cut_pack;
  std::vector< int32_t > iparent_pack;
  std::vector< F >  dist_nearest_pack;

  // Host copies for batched GEMM/SYRK arrays
  std::vector< F* > dmat_array, bf_array, zmat_array;
  std::vector< int64_t > m_array, n_array, k_array, ld_array;

  // abb: can get rid of these variables
  std::vector< F > alpha_array;
  std::vector< oneapi::mkl::transpose > trans_array, nontrans_array;
  std::vector< int64_t > groupsize_array;

  device_task_container tasks_device;


  auto concat_iterable = []( auto& a, const auto& b ) {
    a.insert( a.end(), b.begin(), b.end() );
  };


  size_t ntask    = 0;
  size_t total_npts     = 0;
  size_t total_nbe_nbe  = 0;
  size_t total_nbe_npts = 0;
  size_t total_nshells  = 0;
  size_t total_ncut     = 0;
  size_t memleft = dynmem_sz;

  // Offset memory by the static requirement of an extra pointer element
  // for each of the size batch arrays in MAGMA
  memleft -= 4 * sizeof(int64_t); //M,N,K,LDA[LDB,LDC]

  auto task_it = task_begin;
  while( task_it != task_end ) {

    auto iAtom      = task_it->iParent;
    auto points     = task_it->points    ;
    auto weights    = task_it->weights   ;
    auto shell_list = task_it->shell_list;
    auto nbe        = task_it->nbe;
    auto dist_nearest = task_it->dist_nearest;

    // Generate map from compressed to non-compressed matrices
    auto [submat_cut_32bit, foo] = 
      integrator::gen_compressed_submat_map( basis, shell_list, basis.nbf(), basis.nbf() );
    std::vector< std::array<int64_t,3> > submat_cut;
    submat_cut.reserve( submat_cut_32bit.size() );
    for( const auto& x : submat_cut_32bit ) 
      submat_cut.push_back( {(int64_t)x[0], (int64_t)x[1], (int64_t)x[2]} );

    size_t ncut     = submat_cut.size();
    size_t nshells  = shell_list.size();
    size_t npts     = points.size();


    size_t mem_points  = 3 * npts;
    size_t mem_weights = npts;

    size_t mem_shells     = nshells;
    size_t mem_shell_list = nshells;
    size_t mem_shell_offs = nshells;
    size_t mem_submat_cut = 3 * ncut;

    size_t mem_nbe_scr    = nbe * nbe;
    size_t mem_zmat       = nbe * npts;

    size_t mem_bf         = nbe * npts;
    size_t mem_dbfx       = mem_bf;
    size_t mem_dbfy       = mem_bf;
    size_t mem_dbfz       = mem_bf;

    size_t mem_den        = npts;
    size_t mem_denx       = npts;
    size_t mem_deny       = npts;
    size_t mem_denz       = npts;

    size_t mem_eps        = npts;
    size_t mem_gamma      = npts;
    size_t mem_vrho       = npts;
    size_t mem_vgamma     = npts;

    //size_t mem_partition_scr = natoms * npts;
    size_t mem_dist_scr      = natoms * npts;
    size_t mem_iparent       = npts;
    size_t mem_dist_nearest  = npts;

    size_t mem_batch_mat_arr = 3; // dmat/zmat/bf
    size_t mem_batch_sz_arr  = 4; // M/N/K/LDA/[LDB/LDC]
    size_t mem_task      = 1;


    size_t mem_req_batch =
      mem_points            * sizeof(F) +
      mem_weights           * sizeof(F) +
      mem_shells            * sizeof(Shell<F>) +
      mem_shell_list        * sizeof(size_t) +
      mem_shell_offs        * sizeof(size_t) +
      mem_submat_cut        * sizeof(int64_t) +
      mem_nbe_scr           * sizeof(F) +
      mem_zmat              * sizeof(F) +
      mem_bf                * sizeof(F) +
      mem_dbfx              * sizeof(F) +
      mem_dbfy              * sizeof(F) +
      mem_dbfz              * sizeof(F) +
      mem_den               * sizeof(F) +
      mem_denx              * sizeof(F) +
      mem_deny              * sizeof(F) +
      mem_denz              * sizeof(F) +
      mem_eps               * sizeof(F) +
      mem_gamma             * sizeof(F) +
      mem_vrho              * sizeof(F) +
      mem_vgamma            * sizeof(F) +
      //mem_partition_scr     * sizeof(F) +
      mem_dist_scr          * sizeof(F) +
      mem_iparent           * sizeof(int32_t) +
      mem_dist_nearest      * sizeof(F) +
      mem_batch_mat_arr     * sizeof(F*) +
      mem_batch_sz_arr      * sizeof(int64_t) +
      mem_task              * sizeof(sycl::XCTaskDevice<F>);

    //std::cout << "Memory requirement for task " << ntask+1 << " " << mem_req_batch << " memleft " << memleft << std::endl;

    if( mem_req_batch > memleft ) break;

    // Update memory and increment task iterator
    memleft -= mem_req_batch;
    ntask++;
    task_it++;

    // Update counters
    total_npts     += npts;
    total_nbe_nbe  += nbe*nbe;
    total_nbe_npts += nbe*npts;
    total_nshells  += nshells;
    total_ncut     += ncut;

    // Compute offsets
    std::vector< size_t > shell_offs( nshells );
    shell_offs.at(0) = 0;
    for( auto i = 1ul; i < nshells; ++i )
      shell_offs.at(i) = shell_offs.at(i-1) +
                           basis.at( shell_list.at(i-1) ).size();


    // Pack the data on host
    concat_iterable( points_pack,  points  );
    concat_iterable( weights_pack, weights );
    concat_iterable( shell_list_pack, shell_list );
    concat_iterable( shell_offs_pack, shell_offs );
    concat_iterable( submat_cut_pack, submat_cut );

    m_array.emplace_back( nbe  );
    n_array.emplace_back( npts );
    k_array.emplace_back( nbe  );

    ld_array.emplace_back( nbe  );

    alpha_array.emplace_back( 1. );
    nontrans_array.emplace_back( oneapi::mkl::transpose::nontrans );
    trans_array.emplace_back( oneapi::mkl::transpose::trans );
    groupsize_array.emplace_back( 1. );

    iparent_pack.insert( iparent_pack.end(), npts, iAtom );
    dist_nearest_pack.insert( dist_nearest_pack.end(), npts, dist_nearest );

    // Add task
    tasks_device.emplace_back();

    tasks_device.back().nbe          = nbe;
    tasks_device.back().npts         = npts;
    tasks_device.back().ncut         = ncut;
    tasks_device.back().nshells      = nshells;
    tasks_device.back().iParent      = iAtom;
    tasks_device.back().dist_nearest = dist_nearest;
  }


  std::cout << "XCDeviceData will stack allocate for " << tasks_device.size() << " tasks" << std::endl;

  // Allocate out of dynamic memory
  buffer_adaptor mem( dynmem_ptr, dynmem_sz );

  // (possibly) Large types
  important_shells_device = mem.aligned_alloc<Shell<F>>( total_nshells );
  device_tasks            = mem.aligned_alloc<sycl::XCTaskDevice<F>>( ntask );

  // 64-bit types
  nbe_scr_device     = mem.aligned_alloc<F>( total_nbe_nbe  );
  zmat_device        = mem.aligned_alloc<F>( total_nbe_npts );
  bf_eval_device     = mem.aligned_alloc<F>( total_nbe_npts );
  dbf_x_eval_device  = mem.aligned_alloc<F>( total_nbe_npts );
  dbf_y_eval_device  = mem.aligned_alloc<F>( total_nbe_npts );
  dbf_z_eval_device  = mem.aligned_alloc<F>( total_nbe_npts );

  den_eval_device   = mem.aligned_alloc<F>( total_npts );
  eps_eval_device   = mem.aligned_alloc<F>( total_npts );
  vrho_eval_device  = mem.aligned_alloc<F>( total_npts );

  den_x_eval_device  = mem.aligned_alloc<F>( total_npts );
  den_y_eval_device  = mem.aligned_alloc<F>( total_npts );
  den_z_eval_device  = mem.aligned_alloc<F>( total_npts );
  gamma_eval_device  = mem.aligned_alloc<F>( total_npts );
  vgamma_eval_device = mem.aligned_alloc<F>( total_npts );

  points_device_buffer     = mem.aligned_alloc<F>( 3 * total_npts );
  weights_device_buffer    = mem.aligned_alloc<F>( total_npts );
  shell_list_device_buffer = mem.aligned_alloc<size_t>( total_nshells );
  shell_offs_device_buffer = mem.aligned_alloc<size_t>( total_nshells );
  submat_cut_device_buffer = mem.aligned_alloc<int64_t>( 3 * total_ncut );

  dist_scratch_device = mem.aligned_alloc<F>( natoms * total_npts );
  dist_nearest_buffer = mem.aligned_alloc<F>( total_npts );

  dmat_array_device = mem.aligned_alloc<F*>( ntask );
  zmat_array_device = mem.aligned_alloc<F*>( ntask );
  bf_array_device   = mem.aligned_alloc<F*>( ntask );

  m_array_device   = mem.aligned_alloc<int64_t>( ntask );
  n_array_device   = mem.aligned_alloc<int64_t>( ntask );
  k_array_device   = mem.aligned_alloc<int64_t>( ntask );
  ld_array_device = mem.aligned_alloc<int64_t>( ntask );

  // following 3 arrays are required for batched oneMKL gemm
  alpha_array_device     = mem.aligned_alloc<F>( ntask );
  beta_array_device      = mem.aligned_alloc<F>( ntask );
  trans_array_device     = mem.aligned_alloc<oneapi::mkl::transpose>( ntask );
  nontrans_array_device  = mem.aligned_alloc<oneapi::mkl::transpose>( ntask );
  groupsize_array_device = mem.aligned_alloc<int64_t>( ntask );

  iparent_device_buffer = mem.aligned_alloc<int32_t>( total_npts );

  // Update tasks with allocated pointers
  {
  double* points_ptr  = points_device_buffer;
  double* weights_ptr = weights_device_buffer;

  size_t* shell_list_ptr  = shell_list_device_buffer;
  size_t* shell_offs_ptr  = shell_offs_device_buffer;
  int64_t* submat_cut_ptr = submat_cut_device_buffer;
  Shell<F>   * shells_ptr = important_shells_device;
  double*      nbe_ptr    = nbe_scr_device;
  double*      zmat_ptr   = zmat_device;

  double*      bf_ptr     = bf_eval_device;
  double*      dbfx_ptr   = dbf_x_eval_device;
  double*      dbfy_ptr   = dbf_y_eval_device;
  double*      dbfz_ptr   = dbf_z_eval_device;

  double*      den_ptr    = den_eval_device;
  double*      ddenx_ptr  = den_x_eval_device;
  double*      ddeny_ptr  = den_y_eval_device;
  double*      ddenz_ptr  = den_z_eval_device;

  double*      eps_ptr     = eps_eval_device;
  double*      gamma_ptr   = gamma_eval_device;
  double*      vrho_ptr    = vrho_eval_device;
  double*      vgamma_ptr  = vgamma_eval_device;


  double* dist_scratch_ptr      = dist_scratch_device;

  for( auto& task : tasks_device ) {

    task.points     = points_ptr;
    task.weights    = weights_ptr;
    task.shell_list = shell_list_ptr;
    task.shell_offs = shell_offs_ptr;
    task.submat_cut = submat_cut_ptr;

    task.shells  = shells_ptr;
    task.nbe_scr = nbe_ptr;
    task.zmat    = zmat_ptr;
    task.bf      = bf_ptr;
    task.dbfx    = dbfx_ptr;
    task.dbfy    = dbfy_ptr;
    task.dbfz    = dbfz_ptr;
    task.den     = den_ptr;
    task.ddenx   = ddenx_ptr;
    task.ddeny   = ddeny_ptr;
    task.ddenz   = ddenz_ptr;

    task.eps    = eps_ptr;
    task.gamma  = gamma_ptr;
    task.vrho   = vrho_ptr;
    task.vgamma = vgamma_ptr;

    task.dist_scratch      = dist_scratch_ptr;

    auto npts    = task.npts;
    auto nbe     = task.nbe;
    auto nshells = task.nshells;
    auto ncut    = task.ncut;

    points_ptr     += 3 * npts;
    weights_ptr    += npts;
    shell_list_ptr += nshells;
    shell_offs_ptr += nshells;
    submat_cut_ptr += 3 * ncut;

    shells_ptr += nshells;
    nbe_ptr    += nbe * nbe;
    zmat_ptr   += nbe * npts;

    bf_ptr     += nbe * npts;
    dbfx_ptr   += nbe * npts;
    dbfy_ptr   += nbe * npts;
    dbfz_ptr   += nbe * npts;

    den_ptr    += npts;
    ddenx_ptr  += npts;
    ddeny_ptr  += npts;
    ddenz_ptr  += npts;

    eps_ptr    += npts;
    gamma_ptr  += npts;
    vrho_ptr   += npts;
    vgamma_ptr += npts;

    dist_scratch_ptr += natoms * npts;



    // Batched LA
    dmat_array.emplace_back( task.nbe_scr );
    bf_array.emplace_back(   task.bf      );
    zmat_array.emplace_back( task.zmat    );
  }

  } // End task setup




  auto copy_rev = [&]( size_t n, const auto* src, auto* dest, cl::sycl::queue& queue,
                       std::string m ) {
    util::sycl_copy_async( n, dest, src, queue, m );
  };



  try {

  // Send the data to the device
  copy_rev( 3*points_pack.size(), points_pack.data()->data(),
                         points_device_buffer, *master_queue,
                         "send points buffer" );
  copy_rev( weights_pack.size(), weights_pack.data(),
                         weights_device_buffer, *master_queue,
                         "send weights buffer" );

  copy_rev( shell_list_pack.size(), shell_list_pack.data(),
                          shell_list_device_buffer, *master_queue,
                          "send_shell_list_buffer" );
  copy_rev( shell_offs_pack.size(), shell_offs_pack.data(),
                         shell_offs_device_buffer, *master_queue,
                         "send_shell_offs_buffer" );
  copy_rev( 3*submat_cut_pack.size(), submat_cut_pack.data()->data(),
                         submat_cut_device_buffer, *master_queue,
                         "send_submat_cut_buffer"  );

  copy_rev( tasks_device.size(), tasks_device.data(), device_tasks,
                          *master_queue, "send_tasks_device" );


  copy_rev( dmat_array.size(), dmat_array.data(), dmat_array_device,
                         *master_queue, "send dmat_array" );
  copy_rev( zmat_array.size(), zmat_array.data(), zmat_array_device,
                         *master_queue, "send zmat_array" );
  copy_rev( bf_array.size(), bf_array.data(), bf_array_device,
                         *master_queue, "send bf_array" );

  copy_rev( m_array.size(), m_array.data(), m_array_device,
                         *master_queue, "send m_array" );
  copy_rev( n_array.size(), n_array.data(), n_array_device,
                         *master_queue, "send n_array" );
  copy_rev( k_array.size(), k_array.data(), k_array_device,
                         *master_queue, "send k_array" );

  copy_rev( ld_array.size(), ld_array.data(), ld_array_device,
                         *master_queue, "send ld_array" );

  util::sycl_set_zero_async( ntask, beta_array_device, *master_queue, "betaZero" );

  // abb: uncomment this when `-sycl-std=2020` is supported
  // constexpr F alpha_pattern = 1.0;
  // constexpr oneapi::mkl::transpose nontrans_pattern = oneapi::mkl::transpose::nontrans;
  // constexpr oneapi::mkl::transpose trans_pattern    = oneapi::mkl::transpose::trans;
  // constexpr int64_t gs_pattern = 1.0;
  // master_queue->fill(alpha_array_device, alpha_pattern, ntask);
  // master_queue->fill(trans_array_device, trans_pattern, ntask);
  // master_queue->fill(nontrans_array_device, nontrans_pattern, ntask);
  // master_queue->fill(groupsize_array_device, gs_pattern, ntask);

  // abb: remove the next 4 copy statements when `-sycl-std=2020` is supported
  copy_rev( alpha_array.size(), alpha_array.data(), alpha_array_device,
            *master_queue, "send alpha_array" );
  copy_rev( trans_array.size(), trans_array.data(), trans_array_device,
            *master_queue, "send trans_array" );
  copy_rev( nontrans_array.size(), nontrans_array.data(), nontrans_array_device,
            *master_queue, "send nontrans_array" );
  copy_rev( groupsize_array.size(), groupsize_array.data(), groupsize_array_device,
            *master_queue, "send groupsize_array" );

  copy_rev( iparent_pack.size(), iparent_pack.data(),
                         iparent_device_buffer, *master_queue, "send iparent"  );
  copy_rev( dist_nearest_pack.size(), dist_nearest_pack.data(),
                         dist_nearest_buffer, *master_queue, "send dist_nearest" );

  } catch(...) {
    //teardown_();  throw;
    throw;
  }


  // To avoid packed vectors going out of scope
  util::sycl_device_sync( *master_queue );

  return std::make_tuple(task_it, tasks_device);
}


// Explicit Instantiations
template class XCSyclData<double>;

}
