#include "cuda_aos_scheme1.hpp"
#include "buffer_adaptor.hpp"
#include <gauxc/util/div_ceil.hpp>

namespace GauXC {
namespace detail {

CudaAoSScheme1::Data::~Data() noexcept = default;

CudaAoSScheme1::Data::Data(bool batch_l3_blas) {
  if( batch_l3_blas ) {
    // Create BLAS streams
    blas_streams.resize(4);
    blas_handles.resize(4);
    for( auto i = 0; i < 4; ++i )
      cublasSetStream( blas_handles[i], blas_streams[i] );
  }
}

size_t CudaAoSScheme1::Data::get_ldatoms() {
  constexpr auto weight_unroll = CudaAoSScheme1::weight_unroll;
  return util::div_ceil( natoms, weight_unroll ) * weight_unroll;
}

size_t CudaAoSScheme1::Data::get_static_mem_requirement() {
  return 0;
}

size_t CudaAoSScheme1::Data::get_mem_req( const host_task_type& task,
  const BasisSetMap& basis_map ) {

  const auto ldatoms = get_ldatoms();
  const auto mem_dist_scr = ldatoms * task.npts;
  const auto mem_dist_ner = task.npts;
  const auto mem_iparent  = task.npts;

  size_t base_size = base_type::get_mem_req(task, basis_map);
  return base_size + 
    (mem_dist_scr + mem_dist_ner) * sizeof(double) + 
    mem_iparent * sizeof(int32_t);
}

CudaAoSScheme1::Data::device_buffer_t CudaAoSScheme1::Data::add_extra_to_indirection( 
  std::vector<XCDeviceTask>& tasks,
  device_buffer_t buf ) {
  
  // Host Packing Arrays
  std::vector< int32_t > iparent_pack;
  std::vector< double >  dist_nearest_pack;

  // Allocate additional device memory 
  auto [ ptr, sz ] = buf;
  buffer_adaptor mem( ptr, sz );


  const auto ldatoms = get_ldatoms();
  dist_scratch_device = 
    mem.aligned_alloc<double>( ldatoms * total_npts_task_batch, sizeof(double2) );
  dist_nearest_device = mem.aligned_alloc<double>( total_npts_task_batch );
  iparent_device = mem.aligned_alloc<int32_t>( total_npts_task_batch );

  double* dist_scratch_ptr = dist_scratch_device;
  // Pack additional host data and send
  for( auto& task : tasks ) {
    iparent_pack.insert( iparent_pack.end(), task.npts, task.iParent );
    dist_nearest_pack.insert( dist_nearest_pack.end(), task.npts, 
      task.dist_nearest );

    // Extra indirection for dist scratch
    task.dist_scratch  = dist_scratch_ptr;
    dist_scratch_ptr   += ldatoms * task.npts;
  }

  copy_async( iparent_pack.size(), iparent_pack.data(), 
              iparent_device, "send iparent"  );
  copy_async( dist_nearest_pack.size(), dist_nearest_pack.data(), 
              dist_nearest_device, "send dist_nearest" );

  master_queue_synchronize(); 

  return device_buffer_t{ mem.stack(), mem.nleft() };
}

void CudaAoSScheme1::Data::allocate_rab() {
  
  const auto ldatoms = get_ldatoms();

  buffer_adaptor mem( dynmem_ptr, dynmem_sz );
  rab_device = mem.aligned_alloc<double>( natoms * ldatoms, sizeof(double2) );

  // Update
  dynmem_ptr = mem.stack();
  dynmem_sz  = mem.nleft(); 

}

void CudaAoSScheme1::Data::send_rab(const MolMeta& meta) {

  const auto ldatoms = get_ldatoms();
  std::vector<double> rab_inv(natoms*natoms);
  for( int i = 0 ; i < (int)natoms; ++i ) rab_inv[i] = 1./meta.rab().data()[i];

  util::cuda_copy_2d( rab_device,     ldatoms * sizeof(double),
                      rab_inv.data(), natoms * sizeof(double),
                      natoms * sizeof(double), natoms, "RAB H2D" );

}

size_t CudaAoSScheme1::Data::get_submat_chunk_size(int32_t LDA, int32_t dev_id) {

  constexpr auto max_submat_blocks = CudaAoSScheme1::max_submat_blocks;

  int l2_cache_size;
  cudaDeviceGetAttribute(&l2_cache_size, cudaDevAttrL2CacheSize, dev_id);

  int l2_block_size = (int) sqrt(0.75 * ((double) l2_cache_size / 8));
  int min_block_size = LDA / max_submat_blocks;

  int block_size = std::max(l2_block_size, min_block_size);
  block_size = std::min(block_size, LDA);

  return block_size;

}


}
}
