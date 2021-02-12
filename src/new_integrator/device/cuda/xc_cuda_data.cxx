#include "device/cuda/xc_cuda_data.hpp"
#include <gauxc/util/div_ceil.hpp>

#include "device/buffer_adaptor.hpp"
#include "common/integrator_common.hpp"
#include "device/cuda/cuda_device_properties.hpp"

namespace GauXC {


namespace integrator::device {

template <typename T>
std::shared_ptr< XCDeviceData<T> > make_device_data() {
  return std::make_shared< XCCudaData<T> >();
}

template std::shared_ptr<XCDeviceData<double>> make_device_data();

}








template <typename F>
XCCudaData<F>::XCCudaData( bool _batch_l3_blas ):
#ifdef GAUXC_ENABLE_MAGMA
  batch_l3_blas(_batch_l3_blas)  
#else
  batch_l3_blas(false)  
#endif
{

  // TODO: Expose this
  double fill_fraction = 0.9;

  cudaError_t stat;

  // Get Total Available Memory
  size_t cuda_avail, cuda_total;
  stat = cudaMemGetInfo( &cuda_avail, &cuda_total );
  GAUXC_CUDA_ERROR( "MemInfo Failed", stat );

  // Allocate up to fill_fraction
  devmem_sz = fill_fraction * cuda_avail;
  stat = cudaMalloc( &device_ptr, devmem_sz );
  GAUXC_CUDA_ERROR( "CUDA Malloc Failed", stat );

  // Create CUDA Stream and CUBLAS Handles and make them talk to eachother
  master_stream = std::make_unique< util::cuda_stream >();
  master_handle = std::make_unique< util::cublas_handle >();

  cublasSetStream( *master_handle, *master_stream );

#ifdef GAUXC_ENABLE_MAGMA
  // Create MAGMA Queue from CUDA Stream and CUBLAS Handle
  master_magma_queue = 
    std::make_unique< util::magma_queue >( 0, *master_stream, *master_handle );
#endif

  if( not batch_l3_blas ) {

    // Create BLAS streams
    blas_streams.resize(4);
    blas_handles.resize(4);
    for( auto i = 0; i < 4; ++i )
      cublasSetStream( blas_handles[i], blas_streams[i] );

  }

}



template <typename F>
XCCudaData<F>::~XCCudaData() noexcept {
  if( device_ptr ) util::cuda_free( device_ptr );
} 







template <typename F>
void XCCudaData<F>::allocate_static_data( size_t _natoms,
                                          size_t _n_deriv, 
                                          size_t _nbf,
                                          size_t _nshells ) {


  // Save state
  nshells = _nshells;
  nbf     = _nbf; 
  n_deriv = _n_deriv; 
  natoms  = _natoms;

  LDatoms = util::div_ceil( natoms, cuda::weight_unroll ) * cuda::weight_unroll;

  // Allocate static memory with proper alignment
  buffer_adaptor mem( device_ptr, devmem_sz );

  shells_device     = mem.aligned_alloc<Shell<F>>( nshells );
  exc_device        = mem.aligned_alloc<F>( 1 );
  nel_device        = mem.aligned_alloc<F>( 1 );
  acc_scr_device    = mem.aligned_alloc<F>( 1 );
  rab_device        = mem.aligned_alloc<F>( LDatoms * natoms, sizeof(double2));
  coords_device     = mem.aligned_alloc<F>( 3 * natoms );

  vxc_device  = mem.aligned_alloc<F>( nbf * nbf );
  dmat_device = mem.aligned_alloc<F>( nbf * nbf );

  // Get current stack location
  dynmem_ptr = mem.stack();
  dynmem_sz  = mem.nleft(); 

}




using task_iterator = std::vector< XCTask >::iterator;
template <typename F>
using device_task_container = std::vector< cuda::XCTaskDevice<F> >;

template <typename F>
std::tuple< task_iterator, device_task_container<F> >
  XCCudaData<F>::generate_buffers( const BasisSet<F>& basis,
                                   task_iterator      task_begin,
                                   task_iterator      task_end    ) {

  // Host data packing arrays
  std::vector< std::array<double,3> > points_pack;
  std::vector< double > weights_pack;
  std::vector< size_t > shell_list_pack;
  std::vector< size_t > shell_offs_pack;
  std::vector< std::array<int32_t, 3> > submat_cut_pack;
  std::vector< int32_t > submat_block_pack;
  std::vector< int32_t > iparent_pack;
  std::vector< double >  dist_nearest_pack;

  // Host copies for batched GEMM/SYRK arrays
  std::vector< double* > dmat_array, bf_array, zmat_array;
  std::vector< int > m_array, n_array, k_array, lda_array, ldb_array, ldc_array;

  device_task_container tasks_device;


  auto concat_iterable = []( auto& a, const auto& b ) {
    a.insert( a.end(), b.begin(), b.end() );
  };


  size_t ntask          = 0;
  size_t total_npts     = 0;
  size_t total_nbe_nbe  = 0;
  size_t total_nbe_npts = 0;
  size_t total_nshells  = 0;
  size_t total_ncut     = 0;
  size_t total_nblock   = 0;
  size_t memleft = dynmem_sz;

  uint32_t submat_chunk_size = cuda::get_submat_cut_block(nbf, 0);

  // Offset memory by the static requirement of an extra pointer element 
  // for each of the size batch arrays in MAGMA
  memleft -= 6 * sizeof(int); //M,N,K,LDA,LDB,LDC

  auto task_it = task_begin;
  while( task_it != task_end ) {

    auto iAtom      = task_it->iParent;
    auto points     = task_it->points    ;
    auto weights    = task_it->weights   ;
    auto shell_list = task_it->shell_list;
    auto nbe        = task_it->nbe;
    auto dist_nearest = task_it->dist_nearest;

    // Generate map from compressed to non-compressed matrices
    auto [submat_cut, submat_block] = integrator::gen_compressed_submat_map( basis, shell_list, nbf, submat_chunk_size );
    size_t ncut     = submat_cut.size();
    size_t nblock   = submat_block.size();
    size_t nshells  = shell_list.size();
    size_t npts     = points.size();


    size_t mem_points  = 3 * npts; 
    size_t mem_weights = npts;     

    size_t mem_shells     = nshells;
    size_t mem_shell_list = nshells;
    size_t mem_shell_offs = nshells;
    size_t mem_submat_cut = 3 * ncut;
    size_t mem_submat_block = nblock;

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
    size_t mem_dist_scr      = LDatoms * npts;
    size_t mem_iparent       = npts;
    size_t mem_dist_nearest  = npts;

    size_t mem_batch_mat_arr = 3; // dmat/zmat/bf
    size_t mem_batch_sz_arr  = 6; // M/N/K/LDA/LDB/LDC
    size_t mem_task      = 1;


    size_t mem_req_batch = 
      mem_points            * sizeof(double) + 
      mem_weights           * sizeof(double) +    
      mem_shells            * sizeof(Shell<F>) +             
      mem_shell_list        * sizeof(size_t) +
      mem_shell_offs        * sizeof(size_t) + 
      mem_submat_cut        * sizeof(int32_t) +
      mem_submat_block      * sizeof(int32_t) +
      mem_nbe_scr           * sizeof(double) +
      mem_zmat              * sizeof(double) +
      mem_bf                * sizeof(double) +
      mem_dbfx              * sizeof(double) +
      mem_dbfy              * sizeof(double) +
      mem_dbfz              * sizeof(double) +
      mem_den               * sizeof(double) +
      mem_denx              * sizeof(double) +
      mem_deny              * sizeof(double) +
      mem_denz              * sizeof(double) +
      mem_eps               * sizeof(double) +
      mem_gamma             * sizeof(double) +
      mem_vrho              * sizeof(double) +
      mem_vgamma            * sizeof(double) +
      //mem_partition_scr     * sizeof(double) +
      mem_dist_scr          * sizeof(double) +
      mem_iparent           * sizeof(int32_t) +
      mem_dist_nearest      * sizeof(double) +
      mem_batch_mat_arr     * sizeof(double*) +
      mem_batch_sz_arr      * sizeof(int32_t) +
      mem_task              * sizeof(cuda::XCTaskDevice<F>);

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
    total_nblock   += nblock;

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
    concat_iterable( submat_block_pack, submat_block );

    m_array.emplace_back( npts  );
    n_array.emplace_back( nbe );
    k_array.emplace_back( nbe  );

    lda_array.emplace_back( nbe  );
    ldb_array.emplace_back( npts );
    ldc_array.emplace_back( npts );

    iparent_pack.insert( iparent_pack.end(), npts, iAtom );
    dist_nearest_pack.insert( dist_nearest_pack.end(), npts, dist_nearest );

    // Add task
    tasks_device.emplace_back();

    tasks_device.back().nbe          = nbe;
    tasks_device.back().npts         = npts;
    tasks_device.back().ncut         = ncut;
    tasks_device.back().nblock       = nblock;
    tasks_device.back().nshells      = nshells;
    tasks_device.back().iParent      = iAtom;
    tasks_device.back().dist_nearest = dist_nearest;
  }


  std::cout << "XCDeviceData will stack allocate for " << tasks_device.size() << " tasks"; 
  std::cout << " Using chunk size of " << submat_chunk_size << std::endl;

  // Allocate out of dynamic memory
  buffer_adaptor mem( dynmem_ptr, dynmem_sz );

  // (possibly) Large types
  important_shells_device = mem.aligned_alloc<Shell<F>>( total_nshells );
  device_tasks            = mem.aligned_alloc<cuda::XCTaskDevice<F>>( ntask );

  // 64-bit types
  nbe_scr_device     = mem.aligned_alloc<double>( total_nbe_nbe  );
  zmat_device        = mem.aligned_alloc<double>( total_nbe_npts );
  bf_eval_device     = mem.aligned_alloc<double>( total_nbe_npts );
  dbf_x_eval_device  = mem.aligned_alloc<double>( total_nbe_npts );
  dbf_y_eval_device  = mem.aligned_alloc<double>( total_nbe_npts );
  dbf_z_eval_device  = mem.aligned_alloc<double>( total_nbe_npts );

  den_eval_device   = mem.aligned_alloc<double>( total_npts );
  eps_eval_device   = mem.aligned_alloc<double>( total_npts );
  vrho_eval_device  = mem.aligned_alloc<double>( total_npts );

  den_x_eval_device  = mem.aligned_alloc<double>( total_npts );
  den_y_eval_device  = mem.aligned_alloc<double>( total_npts );
  den_z_eval_device  = mem.aligned_alloc<double>( total_npts );
  gamma_eval_device  = mem.aligned_alloc<double>( total_npts );
  vgamma_eval_device = mem.aligned_alloc<double>( total_npts );

  points_device_buffer     = mem.aligned_alloc<double>( 3 * total_npts );
  weights_device_buffer    = mem.aligned_alloc<double>( total_npts );
  shell_list_device_buffer = mem.aligned_alloc<size_t>( total_nshells );
  shell_offs_device_buffer = mem.aligned_alloc<size_t>( total_nshells );
  submat_cut_device_buffer = mem.aligned_alloc<int32_t>( 3 * total_ncut );
  submat_block_device_buffer = mem.aligned_alloc<int32_t>( total_nblock );

  dist_scratch_device = mem.aligned_alloc<double>( LDatoms * total_npts, 2 * sizeof(double) );
  dist_nearest_buffer = mem.aligned_alloc<double>( total_npts );

  dmat_array_device = mem.aligned_alloc<double*>( ntask );
  zmat_array_device = mem.aligned_alloc<double*>( ntask );
  bf_array_device   = mem.aligned_alloc<double*>( ntask );

  // 32-bit types
  m_array_device   = mem.aligned_alloc<int32_t>( ntask + 1 );
  n_array_device   = mem.aligned_alloc<int32_t>( ntask + 1 );
  k_array_device   = mem.aligned_alloc<int32_t>( ntask + 1 );
  lda_array_device = mem.aligned_alloc<int32_t>( ntask + 1 );
  ldb_array_device = mem.aligned_alloc<int32_t>( ntask + 1 );
  ldc_array_device = mem.aligned_alloc<int32_t>( ntask + 1 );

  iparent_device_buffer = mem.aligned_alloc<int32_t>( total_npts );


  // Update tasks with allocated pointers
  {
  double* points_ptr  = points_device_buffer;
  double* weights_ptr = weights_device_buffer;

  size_t* shell_list_ptr  = shell_list_device_buffer;
  size_t* shell_offs_ptr  = shell_offs_device_buffer;
  int32_t* submat_cut_ptr = submat_cut_device_buffer;
  int32_t* submat_block_ptr = submat_block_device_buffer;
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
    task.submat_block = submat_block_ptr;
    
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
    auto nblock  = task.nblock;

    points_ptr     += 3 * npts;
    weights_ptr    += npts;
    shell_list_ptr += nshells;
    shell_offs_ptr += nshells;
    submat_cut_ptr += 3 * ncut;
    submat_block_ptr += nblock;
    
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

    dist_scratch_ptr += LDatoms * npts;



    // Batched LA
    dmat_array.emplace_back( task.nbe_scr );
    bf_array.emplace_back(   task.bf      );
    zmat_array.emplace_back( task.zmat    );
  }

  } // End task setup




  auto copy_rev = [&]( size_t n, const auto* src, auto* dest, cudaStream_t stream,
                       std::string m ) {
    util::cuda_copy_async( n, dest, src, stream, m );
  };



  try {

  // Send the data to the device
  copy_rev( 3*points_pack.size(), points_pack.data()->data(), 
                         points_device_buffer, *master_stream, 
                         "send points buffer" ); 
  copy_rev( weights_pack.size(), weights_pack.data(), 
                         weights_device_buffer, *master_stream, 
                         "send weights buffer" ); 

  copy_rev( shell_list_pack.size(), shell_list_pack.data(), 
                          shell_list_device_buffer, *master_stream, 
                          "send_shell_list_buffer" );
  copy_rev( shell_offs_pack.size(), shell_offs_pack.data(), 
                         shell_offs_device_buffer, *master_stream, 
                         "send_shell_offs_buffer" );
//  std::cout << "Element size " << sizeof(std::get<0>(submat_cut_pack[0]) << std::endl;
  copy_rev( 3 * submat_cut_pack.size(), submat_cut_pack.data()->data(), 
                         submat_cut_device_buffer, *master_stream, 
                         "send_submat_cut_buffer"  ); 
  copy_rev( submat_block_pack.size(), submat_block_pack.data(), 
                         submat_block_device_buffer, *master_stream, 
                         "send_submat_block_buffer"  ); 

  copy_rev( tasks_device.size(), tasks_device.data(), device_tasks, 
                          *master_stream, "send_tasks_device" );


  copy_rev( dmat_array.size(), dmat_array.data(), dmat_array_device, 
                         *master_stream, "send dmat_array" );
  copy_rev( zmat_array.size(), zmat_array.data(), zmat_array_device, 
                         *master_stream, "send zmat_array" );
  copy_rev( bf_array.size(), bf_array.data(), bf_array_device, 
                         *master_stream, "send bf_array" );

  copy_rev( m_array.size(), m_array.data(), m_array_device, 
                         *master_stream, "send m_array" );
  copy_rev( n_array.size(), n_array.data(), n_array_device, 
                         *master_stream, "send n_array" );
  copy_rev( k_array.size(), k_array.data(), k_array_device, 
                         *master_stream, "send k_array" );

  copy_rev( lda_array.size(), lda_array.data(), lda_array_device, 
                         *master_stream, "send lda_array" );
  copy_rev( ldb_array.size(), ldb_array.data(), ldb_array_device, 
                         *master_stream, "send ldb_array" );
  copy_rev( ldc_array.size(), ldc_array.data(), ldc_array_device, 
                         *master_stream, "send ldc_array" );

  copy_rev( iparent_pack.size(), iparent_pack.data(), 
                         iparent_device_buffer, *master_stream, "send iparent"  );
  copy_rev( dist_nearest_pack.size(), dist_nearest_pack.data(), 
                         dist_nearest_buffer, *master_stream, "send dist_nearest" );

  } catch(...) {
    //teardown_();  throw;
    throw;
  }


  // To avoid packed vectors going out of scope
  cudaStreamSynchronize( *master_stream );

  return std::make_tuple(task_it, tasks_device);
}


// Explicit Instantiations
template class XCCudaData<double>;

}
