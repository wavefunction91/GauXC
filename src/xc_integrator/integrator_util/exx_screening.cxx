/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#include "exx_screening.hpp"
#include "host/blas.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "exceptions/cuda_exception.hpp"
#include <chrono>

namespace std {
template <typename T>
ostream& operator<<( ostream& out, const vector<T>& v ) {
  for( auto _v : v ) out << _v << " ";
  return out;
}
}

namespace GauXC {

void exx_ek_screening( 
  const BasisSet<double>& basis, const BasisSetMap& basis_map,
  const double* P_abs, size_t ldp, const double* V_shell_max, size_t ldv,
  double eps_E, double eps_K, LocalHostWorkDriver* lwd, 
  exx_detail::host_task_iterator task_begin,
  exx_detail::host_task_iterator task_end ) {

  const size_t nbf     = basis.nbf();
  const size_t nshells = basis.nshells();
  const size_t ntasks  = std::distance(task_begin, task_end);

  std::vector<double> task_max_bf_sum(ntasks);
  std::vector<double> task_max_bfn(nbf * ntasks);

  using hrt_t = std::chrono::high_resolution_clock;
  using dur_t = std::chrono::duration<double>;

  auto coll_st = hrt_t::now();
  #pragma omp parallel
  { // Scope temp mem
  std::vector<double> basis_eval;
  std::vector<double> bfn_max_grid(nbf);

  #pragma omp for schedule(dynamic)
  for(size_t i_task = 0; i_task < ntasks; ++i_task) {
    //std::cout << "ITASK = " << i_task << std::endl;

    const auto& task = *(task_begin + i_task);
    const auto npts = task.points.size();

    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();

    // Basis function shell list
    auto shell_list_bfn_ = task.bfn_screening.shell_list;
    int32_t* shell_list_bfn = shell_list_bfn_.data();
    size_t nshells_bfn = shell_list_bfn_.size();
    size_t nbe_bfn     = 
      basis.nbf_subset( shell_list_bfn_.begin(), shell_list_bfn_.end() );

    // Resize scratch
    basis_eval.resize( nbe_bfn * npts );


    // Evaluate basis functions
    lwd->eval_collocation( npts, nshells_bfn, nbe_bfn, points, basis,
      shell_list_bfn, basis_eval.data() );

    // Compute max bfn sum
    // MBFS = max_i sqrt(W[i]) * \sum_mu B(mu,i)
    double max_bfn_sum = 0.;
    for( auto ipt = 0ul; ipt < npts; ++ipt ) {
      double tmp = 0.;
      for( auto ibf = 0ul; ibf < nbe_bfn; ++ibf ) {
      //if(ipt == 0 and ibf == 0) printf("[HOST] ITASK = %lu BF(0,%lu) = %.6e\n", i_task,ibf,basis_eval[ ibf + ipt*nbe_bfn ]);
        tmp += std::abs( basis_eval[ ibf + ipt*nbe_bfn ] );
      }
      //printf("[HOST] ITASK = %lu SUM(%lu) = %.6e\n", i_task,ipt, tmp;
      max_bfn_sum = std::max( max_bfn_sum, std::sqrt(weights[ipt])*tmp );
    }
    task_max_bf_sum[i_task] = max_bfn_sum;
    //printf("[HOST] ITASK = %lu MAX_SUM = %.6e\n", i_task, max_bfn_sum);
    //printf("[HOST] ITASK = %lu NBE = %lu NPTS = %lu \n", i_task, nbe_bfn, npts);

    // Compute max value for each bfn over grid
    bfn_max_grid.resize(nbe_bfn);
    for( auto ibf = 0ul; ibf < nbe_bfn; ++ibf ) {
      double tmp = 0.;
      for( auto ipt = 0ul; ipt < npts; ++ipt ) {
        tmp = std::max(tmp,
          std::sqrt(weights[ipt]) *
          std::abs(basis_eval[ibf + ipt*nbe_bfn])
        );
      }
      bfn_max_grid[ibf] = tmp;
    }
    //printf("[HOST] ITASK = %d MAX_BFN(0) = %.6e\n", i_task, bfn_max_grid[0]);

    // Place max bfn into larger array
    auto task_max_bfn_it = task_max_bfn.data() + i_task*nbf;
    size_t ibf = 0ul;
    for( auto i = 0ul; i < nshells_bfn; ++i ) {
      const auto ish = shell_list_bfn[i];
      const auto sh_sz = basis_map.shell_size(ish);
      const auto sh_off = basis_map.shell_to_first_ao(ish);

      for( auto j = 0; j < sh_sz; ++j ) {
        task_max_bfn_it[j + sh_off] = bfn_max_grid[j + ibf];
      }

      ibf += sh_sz;
    }

  } // Loop over tasks
  } // Memory Scope
  auto coll_en = hrt_t::now();
  std::cout << "... done " << dur_t(coll_en-coll_st).count() << std::endl;
  

  // Compute approx F_i^(k) = |P_ij| * B_j^(k) 
  auto gemm_st = hrt_t::now();
  std::vector<double> task_approx_f( nbf * ntasks );
  blas::gemm( 'N', 'N', nbf, ntasks, nbf, 1., P_abs, ldp,
    task_max_bfn.data(), nbf, 0., task_approx_f.data(), nbf );

  auto gemm_en = hrt_t::now();
  std::cout << "... done " << dur_t(gemm_en-gemm_st).count() << std::endl;

  //std::cout << "CPU FMAX BFN = ";
  //for( auto x : task_approx_f ) std::cout << x << " ";
  //std::cout << std::endl;
#if 0
  // Compute EK shells list for each task
  std::vector<std::vector<int32_t>> ek_shell_task( ntasks );
  std::vector<double> max_F_shells(nshells);
  size_t i_task = 0;
  for( auto task_it = task_begin; task_it != task_end; task_it++, i_task++ ) {
    //std::cout << "ITASK = " << i_task << std::endl;
    std::set<int32_t> ek_shells;

    // Collapse max_F over shells
    const double* max_F_approx_bfn = task_approx_f.data() + i_task*nbf;
    for( auto ish = 0ul, ibf = 0ul; ish < nshells; ++ish) {
      const auto sh_sz = basis[ish].size();
      double tmp = 0.;
      for( auto i = 0; i < sh_sz; ++i ) {
        tmp = std::max( tmp, std::abs(max_F_approx_bfn[ibf + i]) );
      }
      max_F_shells[ish] = tmp;
      ibf += sh_sz;
    }

    // Compute important shell set
    const double max_bf_sum = task_max_bf_sum[i_task];
    for( auto j = 0ul; j < nshells; ++j )
    for( auto i = j;   i < nshells; ++i ) {
      const auto V_ij = V_shell_max[i + j*ldv];
      const auto F_i  = max_F_shells[i];
      const auto F_j  = max_F_shells[j];

      const double eps_E_compare = F_i * F_j * V_ij;
      const double eps_K_compare = std::max(F_i, F_j) * V_ij * max_bf_sum;
      if( eps_K_compare > eps_K or eps_E_compare > eps_E)  {
        ek_shells.insert(i); 
        ek_shells.insert(j); 
      }
    }

    // Append to list
    task_it->cou_screening.shell_list =
      std::vector<int32_t>(ek_shells.begin(), ek_shells.end());
    task_it->cou_screening.nbe = 
      basis.nbf_subset( ek_shells.begin(), ek_shells.end() );

    std::cout << "I_TASK " << i_task << " " << task_it->cou_screening.shell_list.size();
    std::cout << std::endl;
  } // Loop over tasks
#else

  //std::cout << "CPU FMAX SHELLS = ";
  auto list_st = hrt_t::now();
  #pragma omp parallel for schedule(dynamic)
  for(size_t i_task = 0; i_task < ntasks; ++i_task) {
    //std::cout << "ITASK = " << i_task << std::endl;
    std::vector<uint32_t> task_ek_shells(util::div_ceil(nshells,32),0);
    std::vector<double> max_F_shells(nshells);

    // Collapse max_F over shells
    const double* max_F_approx_bfn = task_approx_f.data() + i_task*nbf;
    for( auto ish = 0ul, ibf = 0ul; ish < nshells; ++ish) {
      const auto sh_sz = basis[ish].size();
      double tmp = 0.;
      for( auto i = 0; i < sh_sz; ++i ) {
        tmp = std::max( tmp, std::abs(max_F_approx_bfn[ibf + i]) );
      }
      max_F_shells[ish] = tmp;
      ibf += sh_sz;
    }
    //for(auto x : max_F_shells) std::cout << x << " ";

    auto task_it = task_begin + i_task;
    // Compute important shell set
    const double max_bf_sum = task_max_bf_sum[i_task];
    for( auto i = 0ul; i < nshells; ++i )
    for( auto j = 0ul; j <= i;      ++j ) {
      const auto V_ij = V_shell_max[i + j*ldv];
      const auto F_i  = max_F_shells[i];
      const auto F_j  = max_F_shells[j];

      const double eps_E_compare = F_i * F_j * V_ij;
      const double eps_K_compare = std::max(F_i, F_j) * V_ij * max_bf_sum;
      if( eps_K_compare > eps_K or eps_E_compare > eps_E)  {
        size_t i_block = i / 32;
        size_t j_block = j / 32;
        size_t i_local = i % 32;
        size_t j_local = j % 32;

        task_ek_shells[i_block] |= (1u << i_local); 
        task_ek_shells[j_block] |= (1u << j_local); 
        task_it->cou_screening.shell_pair_list.emplace_back(i,j);
      }
    }

    uint32_t total_shells = 0;
    for( auto x : task_ek_shells ) total_shells += __builtin_popcount(x);

    std::vector<uint32_t> ek_shells; ek_shells.reserve(total_shells);
    for( auto i_block = 0u; i_block < util::div_ceil(nshells,32); ++i_block ) {
    for( unsigned i_local = 0; i_local < 32; ++i_local ) 
    if( task_ek_shells[i_block] & (1u << i_local) ) {
      ek_shells.emplace_back(i_local + i_block*32);
    }
    }

    // Append to list
    task_it->cou_screening.shell_list =
      std::vector<int32_t>(ek_shells.begin(), ek_shells.end());
    task_it->cou_screening.nbe = 
      basis.nbf_subset( ek_shells.begin(), ek_shells.end() );
    //size_t nspt = (nshells*(nshells+1))/2;
    //size_t nsp   = task_it->cou_screening.shell_pair_list.size();
    //std::cout << "I_TASK " << i_task << " " << task_it->cou_screening.shell_list.size() << " " << nsp << " " << nspt << " " << nsp / double(nspt) << std::endl;

  } // Loop over tasks
  //std::cout << std::endl;
  auto list_en = hrt_t::now();
  //for(auto task_it = task_begin; task_it != task_end; ++task_it) {
  //  std::cout << "I_TASK " << std::distance(task_begin,task_it) << " " << task_it->cou_screening.shell_list.size();
  //  std::cout << std::endl;
  //}
  std::cout << "... done " << dur_t(list_en-list_st).count() << std::endl;

#endif
}


void exx_ek_screening( 
  const BasisSet<double>& basis, const BasisSetMap& basis_map,
  const double* P_abs, size_t ldp, const double* V_shell_max, size_t ldv,
  double eps_E, double eps_K, XCDeviceData& device_data, 
  LocalDeviceWorkDriver* lwd, 
  exx_detail::host_task_iterator task_begin,
  exx_detail::host_task_iterator task_end ) {

  const size_t nbf = basis.nbf();
  const auto nshells = basis.nshells();
  const size_t ntasks  = std::distance(task_begin, task_end);

  // Setup EXX EK Screening memory on the device
  device_data.reset_allocations();
  device_data.allocate_static_data_exx_ek_screening( ntasks, nbf, nshells, 
    basis_map.max_l() );
  device_data.send_static_data_density_basis( P_abs, ldp, basis );
  device_data.send_static_data_exx_ek_screening( V_shell_max, ldv, basis_map );

  device_data.zero_exx_ek_screening_intermediates();

  integrator_term_tracker enabled_terms;
  enabled_terms.exx_ek_screening = true;

  // Loop over tasks and form basis-related buffers
  auto task_it = task_begin;
  while( task_it != task_end ) {

    // Determine next task patch, send relevant data (EXX_EK only)
    task_it = 
      device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

    // Evaluate collocation
    lwd->eval_collocation( &device_data );

    // Evaluate EXX EK Screening Basis Statistics
    lwd->eval_exx_ek_screening_bfn_stats( &device_data );

  }

  // Compute approximate F Max
  lwd->eval_exx_ek_screening_approx_fmax( &device_data );
  lwd->exx_ek_collapse_fmat_to_shells( &device_data );

  // Retreive to host
  std::vector<double> task_f_bfn_max(ntasks * nbf);
  std::vector<double> task_f_shl_max(ntasks * nshells);
  device_data.retrieve_exx_ek_approx_fmax_bfn( task_f_bfn_max.data(), nbf );
  device_data.retrieve_exx_ek_approx_fmax_shell( task_f_shl_max.data(), nshells );

  GAUXC_CUDA_ERROR("End Sync", cudaDeviceSynchronize());

  //std::cout << "GPU FMAX BFN = ";
  //for(int i = 0; i < ntasks; ++i)
  //for(int j = 0; j < nbf;    ++j) std::cout << task_f_bfn_max[i + j*ntasks] << " ";
  //std::cout << std::endl;

  //std::cout << "GPU FMAX SHELLS = ";
  //for(int i = 0; i < ntasks; ++i)
  //for(int j = 0; j < nshells;    ++j) std::cout << task_f_shl_max[i + j*ntasks] << " ";
  //std::cout << std::endl;
  
  
}


}
