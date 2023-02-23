/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#include "device/common/exx_ek_screening.hpp"
#include <gauxc/util/div_ceil.hpp>
#include <gauxc/shell.hpp>
#include "device_specific/cuda_util.hpp"
#include "cuda_extensions.hpp"
#include "device_specific/cuda_device_constants.hpp"
#include <cub/device/device_scan.cuh>
#include "buffer_adaptor.hpp"

namespace GauXC {

__global__ void exx_ek_screening_bfn_stats_kernel( size_t ntasks, 
                                                   double      * max_bfn_sum_device,
                                                   double      * bfn_max_device,
                                                   size_t        LDBFM,
                                                   XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.x;
  if( batch_idx >= ntasks ) return;
  
  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;
  auto* basis_eval_device = task.bf;
  const auto* weights_device    = task.weights;
  //double* bfn_max_device  = task.bfn_max;

  const int warp_lane = threadIdx.x % cuda::warp_size;
  const int warp_id   = threadIdx.x / cuda::warp_size;
  const int nwarp  = blockDim.x  / cuda::warp_size;

#if 0
  for(int ibf = warp_id; ibf < nbf; ibf += nwarp) {
    double tmp = 0.0;
    for(int ipt = warp_lane; ipt < npts; ipt += cuda::warp_size) {
      tmp = fmax( tmp,
        std::sqrt(weights_device[ipt]) *
        std::abs( basis_eval_device[ ipt + ibf*npts ] )
      );
    }
    bfn_max_device[ibf] = tmp;
  }
#endif

  //printf("[GPU] ITASK %d TID %d WL %d WID %d\n", batch_idx, threadIdx.x, warp_lane, warp_id);

  // Compute Max Bfn Sum
#if 0
  double max_bfn_sum = 0.0;
  for(int ipt = warp_lane; ipt < npts; ipt += cuda::warp_size) {
    double tmp = 0.0;
    for(int ibf = warp_id; ibf < nbf; ibf += nwarp) {
      tmp += std::abs( basis_eval_device[ ipt + ibf*npts ] );
    }
    max_bfn_sum = fmax( max_bfn_sum, std::sqrt(weights_device[ipt]) * tmp );
  }
  max_bfn_sum = cuda::warp_reduce_max<cuda::warp_size>(max_bfn_sum);
#else


  // First scale the basis functions by the weights
  for(int ipt = warp_lane; ipt < npts; ipt += cuda::warp_size) {
    const auto w = std::sqrt(weights_device[ipt]);
    for(int ibf = warp_id; ibf < nbf; ibf += nwarp) {
      const auto val = basis_eval_device[ ipt + ibf*npts ];
      basis_eval_device[ ipt + ibf*npts ] = w * std::abs(val);
    }
  }
  __syncthreads();




  __shared__ double bf_shared[32][32 + 1];
  __shared__ double bfn_sum_shared[32];
  bfn_sum_shared[warp_lane] = 0.0;
  __syncthreads();

  const int npts_chunks = GauXC::util::div_ceil(npts, cuda::warp_size);
  const int  nbf_chunks = GauXC::util::div_ceil(nbf,  nwarp);
  for(int ipts_chunk = 0; ipts_chunk < npts_chunks; ++ipts_chunk) {
    double tmp_bfn_sum = 0.0;
    const int ipt = ipts_chunk * cuda::warp_size + warp_lane;
  for(int  ibf_chunk = 0;  ibf_chunk <  nbf_chunks; ++ ibf_chunk) {
    const int ibf =  ibf_chunk * nwarp + warp_id;

    bf_shared[warp_id][warp_lane] = 0.0;

    // Load in a block of basis functions
    // Warp lane is the point index and warp ID is the bfn idx
    if(ipt < npts and ibf < nbf) 
      bf_shared[warp_id][warp_lane]  = basis_eval_device[ipt + ibf*npts];
    __syncthreads();

    // Do transpose
    // Warp lane is the bfn index and warp ID is the point idx
    auto tmp = bf_shared[warp_lane][warp_id];
    __syncthreads();
 
    // Do a sum reduce over basis functions for the chunk
    // Now every warp has the local bfn chunk sum in lane 0
    // corresponding to the point represented by the warp Id
    tmp_bfn_sum += cuda::warp_reduce_sum<cuda::warp_size>(tmp);
    
  }
    // At this point, every warp contains the total bfn sum
    // for the point corresponding to its warp id. Update the 
    // local value accordingly
    if(warp_lane == 0) {
      double val = bfn_sum_shared[warp_id];
      bfn_sum_shared[warp_id] = fmax( val, tmp_bfn_sum );
    }
    __syncthreads();

  }

  // Get global maximum
  double max_bfn_sum;
  if(warp_id == 0) {
    auto tmp = bfn_sum_shared[warp_lane];
    max_bfn_sum = cuda::warp_reduce_max<cuda::warp_size>(tmp);
  }

#endif
  if(threadIdx.x == 0) {
    //task.max_bfn_sum = max_bfn_sum;
    max_bfn_sum_device[batch_idx] =  max_bfn_sum;
    //printf("[GPU] ITASK = %d MAX_SUM = %.6e PTR = %x\n", batch_idx, max_bfn_sum, task.bfn_shell_indirection);
    //printf("[GPU] ITASK = %d NBE = %lu NPTS = %lu \n", batch_idx, nbf, npts);
  }


  __syncthreads();
  for(int ibf = warp_id; ibf < nbf; ibf += nwarp) {
    double max_bf = 0;
    for(int ipt = warp_lane; ipt < npts; ipt += cuda::warp_size) {
      const auto val = basis_eval_device[ipt + ibf*npts];
      max_bf = fmax( max_bf, val );
    }

    // Warp reduce bf max
    max_bf = cuda::warp_reduce_max<cuda::warp_size>(max_bf);
    if(warp_lane == 0) {
      //printf("[GPU] ITASK = %d MAX_BFN(0) = %.6e\n", batch_idx, max_bf);
      bfn_max_device[batch_idx + task.bfn_shell_indirection[ibf]*LDBFM] =
        max_bf; 
    }
  }

}


__global__ void exx_ek_collapse_fmax_to_shells_kernel(
  int                  ntask,
  int                  nshells,
  const Shell<double>* shells_device,
  const int32_t*       shell_to_bf,
  const double*        fmax_bfn_device,
  size_t               LDF_bfn,
  double*              fmax_shell_device,
  size_t               LDF_shell
) {

  const int total_nwarp_x   = (blockDim.x * gridDim.x) / cuda::warp_size;
  const int tid_x           = threadIdx.x + blockIdx.x*blockDim.x;
  const int warp_lane       = tid_x % cuda::warp_size;
  const int warp_id_x       = tid_x / cuda::warp_size;


  double sh_buffer[10];

  // Each warp gets a shell
  for(int ish = warp_id_x; ish < nshells; ish += total_nwarp_x) {

    const int sh_sz = shells_device[ish].size();
    const int sh_st = shell_to_bf[ish];

    // Read in tasks in warp-sized chunks
    for(int i_task = warp_lane; i_task < ntask; i_task += cuda::warp_size) {

      // Get shell max
      double sh_max = 0.0;
      int sh_rem = sh_sz;      
      while(sh_rem > 0) {
        int ndo = min(sh_sz, 10);
        // Load in batches to break up dependency tree
        for(int ii = 0; ii < ndo; ++ii) {
          sh_buffer[ii] = fmax_bfn_device[i_task + (ii + sh_st)*LDF_bfn];
        }
        
        for(int ii = 0; ii < ndo; ++ii) {
          sh_max = fmax(sh_max, sh_buffer[ii]);
        }
        sh_rem -= ndo;
      }
      
      // Write to main memory
      fmax_shell_device[i_task + ish*LDF_shell] = sh_max;
      
    }

  }
}





__global__ void exx_ek_shellpair_collision_kernel(
  int32_t       ntasks,
  int32_t       nshells,
  const double* V_max_device,
  size_t        LDV,
  const double* F_max_shl_device,
  size_t        LDF,
  const double* max_bf_sum_device,
  double        eps_E,
  double        eps_K,
  uint32_t*     collisions,
  int           LD_coll,
  uint32_t*     rc_collisions,
  int           LD_rc,
  uint32_t*     counts,
  uint32_t*     rc_counts
) {



  const int tid_x = threadIdx.x + blockIdx.x * blockDim.x;
  const int nt_x  = blockDim.x * gridDim.x;

  for(int i_task = tid_x; i_task < ntasks; i_task += nt_x) {

    const auto max_bf_sum = max_bf_sum_device[i_task];

    for(int ij = 0; ij < LD_coll; ++ij) 
      collisions[i_task * LD_coll + ij] = 0;
    for(int ij = 0; ij < LD_rc; ++ij) 
      rc_collisions[i_task * LD_rc + ij] = 0;

    for(int i_shell = 0; i_shell < nshells;  ++i_shell)
    for(int j_shell = 0; j_shell <= i_shell; ++j_shell) {

      const auto V_ij = V_max_device[i_shell + j_shell*LDV];
      const auto F_i  = F_max_shl_device[i_task + i_shell * LDF];
      const auto F_j  = F_max_shl_device[i_task + j_shell * LDF];

      const double eps_E_compare = F_i * F_j * V_ij;
      const double eps_K_compare = fmax(F_i, F_j) * V_ij * max_bf_sum;
      const bool comp = (eps_K_compare > eps_K or eps_E_compare > eps_E);

      const int ij = i_shell+ ((2*nshells - j_shell - 1)*j_shell)/2;
      const int ij_block = ij / 32;
      const int ij_local = ij % 32;
      collisions[i_task * LD_coll + ij_block] |= comp ? (1u << ij_local) : 0;

      const int i_block = i_shell / 32;
      const int i_local = i_shell % 32;
      rc_collisions[i_task * LD_rc + i_block] |= comp ? (1u << i_local) : 0;
      const int j_block = j_shell / 32;
      const int j_local = j_shell % 32;
      rc_collisions[i_task * LD_rc + j_block] |= comp ? (1u << j_local) : 0;
    }

    uint32_t count = 0;
    for(int ij = 0; ij < LD_coll; ++ij)  count += __popc(collisions[i_task * LD_coll + ij]);
    counts[i_task] = count;

    count = 0;
    for(int ij = 0; ij < LD_rc; ++ij)  count += __popc(rc_collisions[i_task * LD_rc + ij]);
    rc_counts[i_task] = count;
  }

}

__global__ void print_coll(size_t ntasks, size_t nshells, uint32_t* collisions,
  size_t LD_coll) {


  for(auto i_task = 0 ; i_task < ntasks; ++i_task) {

    printf("[GPU] ITASK %d: ", i_task);
    int count = 0;
    for(int i_shell = 0, ij = 0; i_shell < nshells;  ++i_shell      )
    for(int j_shell = 0;         j_shell <= i_shell; ++j_shell, ij++) {

      const int ij_block = ij / 32;
      const int ij_local = ij % 32;
      if( collisions[i_task * LD_coll + ij_block] & (1u << ij_local) ) {
        //printf("(%d, %d) ", i_shell, j_shell);
        count++;
      }
    }
    printf("%d\n", count);

  }
}

__global__ void print_counts(size_t ntasks, uint32_t* counts) {


  for(auto i_task = 0 ; i_task < ntasks; ++i_task) {

    printf("[GPU] ITASK %d: %d\n", i_task,counts[i_task]);

  }
}




template <int32_t buffer_size, typename buffer_type = uint32_t>
__global__ void bitvector_to_position_list_shellpair(
  size_t ntasks,
  size_t nsp,
  size_t LD_bit,
  const uint32_t* collisions,
  const uint32_t* counts,
  uint32_t*       position_list
) {

  constexpr auto warp_size = cuda::warp_size;
  constexpr auto element_size = CHAR_BIT * sizeof(buffer_type);
  constexpr auto buffer_size_bits = element_size * buffer_size;
  __shared__ buffer_type collisions_buffer[warp_size][warp_size][buffer_size];

  // We are converting a large number of small bitvectors into position lists. For this reason, I am assigning a single thread to each bitvector
  // This avoids having to do popcounts and warp wide reductions, but hurts the memory access pattern

  // All threads in a warp must be active to do shared memory loads, so we seperate out the threadId.x
  for (int i_base = threadIdx.y * blockDim.x + blockIdx.x * blockDim.x * blockDim.y; i_base < ntasks; i_base += blockDim.x * blockDim.y * gridDim.x) {
    const int i = i_base + threadIdx.x;
    auto* out = position_list;
    if (i != 0 && i < ntasks) {
      out += counts[i-1];
    } 

    int current = 0;
    size_t nsp_blocks = (nsp + buffer_size_bits - 1) / buffer_size_bits;
    for (int j_block = 0; j_block < nsp_blocks; j_block++) {
      // Each thread has a buffer of length BUFFER_SIZE. All the threads in the warp work to 
      // load this data in a coalesced way (at least as much as possible)
      for (int buffer_loop = 0; buffer_loop < warp_size; buffer_loop += warp_size/buffer_size) {
        const int t_id_x        = threadIdx.x % buffer_size;
        const int buffer_thread = threadIdx.x / buffer_size;
        const int buffer_idx    = buffer_thread + buffer_loop;
        if (j_block * buffer_size_bits + t_id_x * element_size < nsp && i_base + buffer_idx < ntasks) {
          collisions_buffer[threadIdx.y][buffer_idx][t_id_x] = collisions[(i_base + buffer_idx) * LD_bit + j_block * buffer_size + t_id_x];
        }
      }

      __syncwarp();
      if (i < ntasks) {  // Once the data has been loaded, we exclude the threads not corresponding to a bitvector
        // We have loaded in BUFFER_SIZE_BITS elements to be processed by each warp
        for (int j_inner = 0; j_inner < buffer_size_bits && j_block * buffer_size_bits + j_inner < nsp; j_inner++) {
          const int j = buffer_size_bits * j_block + j_inner;
          const int j_int = j_inner / element_size;
          const int j_bit = j_inner % element_size;
          if( collisions_buffer[threadIdx.y][threadIdx.x][j_int] & (1 << (j_bit)) ) {
            out[current++] = j;
          }
        }
      }
      __syncwarp();
    }
  }

}





template <int32_t buffer_size, typename buffer_type = uint32_t>
__global__ void bitvector_to_position_list_shells( 
           size_t  ntasks, 
           size_t  nshells, 
           size_t  LD_bit,
    const uint32_t* collisions, 
    const uint32_t* counts, 
    const int32_t* shell_size,
          uint32_t* position_list, 
           size_t* nbe_list
) {
  constexpr auto warp_size = cuda::warp_size;
  constexpr auto element_size = CHAR_BIT * sizeof(buffer_type);
  constexpr auto buffer_size_bits = element_size * buffer_size;
  __shared__ buffer_type collisions_buffer[warp_size][warp_size][buffer_size];

  // We are converting a large number of small bitvectors into position lists. For this reason, I am assigning a single thread to each bitvector
  // This avoids having to do popcounts and warp wide reductions, but hurts the memory access pattern

  // All threads in a warp must be active to do shared memory loads, so we seperate out the threadId.x
  for (int i_base = threadIdx.y * blockDim.x + blockIdx.x * blockDim.x * blockDim.y; i_base < ntasks; i_base += blockDim.x * blockDim.y * gridDim.x) {
    const int i = i_base + threadIdx.x;
    auto* out = position_list;
    if (i != 0 && i < ntasks) {
      out += counts[i-1];
    } 

    int current = 0;
    size_t nbe = 0;
    size_t nsphere_blocks = (nshells + buffer_size_bits - 1) / buffer_size_bits;
    for (int j_block = 0; j_block < nsphere_blocks; j_block++) {
      // Each thread has a buffer of length BUFFER_SIZE. All the threads in the warp work to 
      // load this data in a coalesced way (at least as much as possible)
      for (int buffer_loop = 0; buffer_loop < warp_size; buffer_loop += warp_size/buffer_size) {
        const int t_id_x        = threadIdx.x % buffer_size;
        const int buffer_thread = threadIdx.x / buffer_size;
        const int buffer_idx    = buffer_thread + buffer_loop;
        if (j_block * buffer_size_bits + t_id_x * element_size < nshells && i_base + buffer_idx < ntasks) {
          collisions_buffer[threadIdx.y][buffer_idx][t_id_x] = collisions[(i_base + buffer_idx) * LD_bit + j_block * buffer_size + t_id_x];
        }
      }

      __syncwarp();
      if (i < ntasks) {  // Once the data has been loaded, we exclude the threads not corresponding to a bitvector
        // We have loaded in BUFFER_SIZE_BITS elements to be processed by each warp
        for (int j_inner = 0; j_inner < buffer_size_bits && j_block * buffer_size_bits + j_inner < nshells; j_inner++) {
          const int j = buffer_size_bits * j_block + j_inner;
          const int j_int = j_inner / element_size;
          const int j_bit = j_inner % element_size;
          if( collisions_buffer[threadIdx.y][threadIdx.x][j_int] & (1 << (j_bit)) ) {
            out[current++] = j;
            nbe += shell_size[j];
          }
        }
      }
      __syncwarp();
    }
    if (i < ntasks) {
      nbe_list[i] = nbe;
    }
  }
}






void exx_ek_screening_bfn_stats( size_t        ntasks,
                                 XCDeviceTask* tasks_device,
                                 double      * max_bfn_sum_device,
                                 double      * bfn_max_device,
                                 size_t        LDBFM,
                                 device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>() ;
  dim3 threads = 1024;//cuda::max_threads_per_thread_block;
  dim3 blocks  = ntasks;
  exx_ek_screening_bfn_stats_kernel<<<blocks, threads, 0, stream >>>(
    ntasks, max_bfn_sum_device, bfn_max_device, LDBFM, tasks_device );

}


void exx_ek_collapse_fmax_to_shells(
  int                  ntask,
  int                  nshells,
  const Shell<double>* shells_device,
  const int32_t*       shell_to_bf,
  const double*        fmax_bfn_device,
  size_t               LDF_bfn,
  double*              fmax_shell_device,
  size_t               LDF_shell,
  device_queue         queue
) {


  cudaStream_t stream = queue.queue_as<util::cuda_stream>() ;
  dim3 threads = 1024;//cuda::max_threads_per_thread_block;
  dim3 blocks  = ntask / cuda::warp_size;
  exx_ek_collapse_fmax_to_shells_kernel<<<blocks, threads, 0, stream >>>(
    ntask, nshells, shells_device, shell_to_bf, fmax_bfn_device, LDF_bfn,
    fmax_shell_device, LDF_shell );

}

void exx_ek_shellpair_collision(
  int32_t       ntasks,
  int32_t       nshells,
  const double* V_max_device,
  size_t        LDV,
  const double* F_max_shl_device,
  size_t        LDF,
  const double* max_bf_sum_device,
  const int32_t* shell_sizes_device,
  double        eps_E,
  double        eps_K,
  uint32_t*     collisions,
  int           LD_coll,
  uint32_t*     counts,
  void*         dyn_stack,
  size_t        dyn_size,
  host_task_iterator tb,
  host_task_iterator te,
  device_queue  queue
) {

  buffer_adaptor stack(dyn_stack, dyn_size);

  size_t LD_rc = util::div_ceil(nshells, 32);
  auto rc_collisions = stack.aligned_alloc<uint32_t>(ntasks * LD_rc);
  auto rc_counts = stack.aligned_alloc<uint32_t>(ntasks);

  

  cudaStream_t stream = queue.queue_as<util::cuda_stream>() ;
  dim3 threads = 1024;//cuda::max_threads_per_thread_block;
  dim3 blocks  = GauXC::util::div_ceil(ntasks,threads.x);
  exx_ek_shellpair_collision_kernel<<<blocks, threads, 0 , stream>>>(
    ntasks, nshells, V_max_device, LDV, F_max_shl_device, LDF, 
    max_bf_sum_device, eps_E, eps_K, collisions, LD_coll, 
    rc_collisions, LD_rc, counts, rc_counts);



  cudaError_t stat;

  size_t prefix_sum_bytes = 0;
  stat = cub::DeviceScan::InclusiveSum( NULL, prefix_sum_bytes,
    counts, counts, ntasks, stream );


  void* prefix_sum_storage = stack.aligned_alloc<char>(prefix_sum_bytes, 16);
  

  // Get inclusive sums
  stat = cub::DeviceScan::InclusiveSum( prefix_sum_storage, prefix_sum_bytes,
    counts, counts, ntasks, stream );
  stat = cub::DeviceScan::InclusiveSum( prefix_sum_storage, prefix_sum_bytes,
    rc_counts.ptr, rc_counts.ptr, ntasks, stream );

  // Get counts after prefix sum
  std::vector<uint32_t> counts_host(ntasks), rc_counts_host(ntasks);
  util::cuda_copy(ntasks, counts_host.data(), counts);
  util::cuda_copy(ntasks, rc_counts_host.data(), rc_counts.ptr);

  uint32_t total_sp_count = counts_host.back();
  std::cout << "TOTAL SPCOUNT = " << total_sp_count << std::endl;
  uint32_t total_s_count = rc_counts_host.back();
  std::cout << "TOTAL SCOUNT = " << total_s_count << std::endl;



  size_t nsp_dense = (nshells * (nshells+1))/2;
  auto position_sp_list_device = stack.aligned_alloc<uint32_t>(total_sp_count);
  auto position_s_list_device  = stack.aligned_alloc<uint32_t>(total_s_count);
  auto nbe_list                = stack.aligned_alloc<size_t>(ntasks);
  {
  dim3 threads(32,32);
  dim3 blocks( util::div_ceil(ntasks, 1024) );
  bitvector_to_position_list_shellpair<8><<<blocks, threads, 0, stream>>>(
    ntasks, nsp_dense, LD_coll, collisions, counts, position_sp_list_device
  );
  bitvector_to_position_list_shells<8><<<blocks, threads, 0, stream>>>(
    ntasks, nshells, LD_rc, rc_collisions.ptr, rc_counts.ptr, shell_sizes_device,
    position_s_list_device.ptr, nbe_list.ptr
  );
  }

  std::vector<uint32_t> position_sp_list(total_sp_count);
  util::cuda_copy(total_sp_count, position_sp_list.data(), position_sp_list_device.ptr, "Position List ShellPair");



  std::vector<std::pair<int,int>> lt_indices(nsp_dense);
  {
  auto it = lt_indices.begin();
  for(int j = 0; j < nshells; ++j)
  for(int i = j; i < nshells; ++i) {
    *(it++) = std::make_pair(i,j);
  }
  }

  std::vector<uint32_t> position_s_list(total_s_count);
  std::vector<size_t> nbe_list_host(ntasks);
  util::cuda_copy(total_s_count, position_s_list.data(), position_s_list_device.ptr, "Position List Shell");
  util::cuda_copy(ntasks, nbe_list_host.data(), nbe_list.ptr, "NBE List");

  for( auto it = tb; it != te; ++it ) {
    {
    size_t begin = (it == tb) ? 0 : counts_host[std::distance(tb,it)-1];
    size_t end   = counts_host[std::distance(tb,it)];

    it->cou_screening.shell_pair_list.resize(end - begin);
    for( auto ij = begin, idx = 0ul; ij < end; ++ij, ++idx) {
      it->cou_screening.shell_pair_list[idx] = 
        lt_indices[position_sp_list[ij]];
    }
    }

    {
    size_t begin = (it == tb) ? 0 : rc_counts_host[std::distance(tb,it)-1];
    size_t end   = rc_counts_host[std::distance(tb,it)];

    it->cou_screening.shell_list.resize(end - begin);
    it->cou_screening.nbe = nbe_list_host[std::distance(tb,it)];
    for( auto ij = begin, idx = 0ul; ij < end; ++ij, ++idx) {
      it->cou_screening.shell_list[idx] = position_s_list[ij]; 
    }
    }
    
  }

  

  
}

}
