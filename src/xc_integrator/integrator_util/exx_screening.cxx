#include "exx_screening.hpp"
#include "host/blas.hpp"
#include <set>

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

  { // Scope temp mem
  std::vector<double> basis_eval;
  std::vector<double> bfn_max_grid(nbf);

  size_t i_task = 0;
  for( auto task_it = task_begin; task_it != task_end; task_it++, i_task++ ) {

    const auto& task = *task_it;
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
        tmp += std::abs( basis_eval[ ibf + ipt*nbe_bfn ] );
      }
      max_bfn_sum = std::max( max_bfn_sum, std::sqrt(weights[ipt])*tmp );
    }
    task_max_bf_sum[i_task] = max_bfn_sum;

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

  // Compute approx F_i^(k) = |P_ij| * B_j^(k) 
  std::vector<double> task_approx_f( nbf * ntasks );
  blas::gemm( 'N', 'N', nbf, ntasks, nbf, 1., P_abs, ldp,
    task_max_bfn.data(), nbf, 0., task_approx_f.data(), nbf );


  // Compute EK shells list for each task
  std::vector<std::vector<int32_t>> ek_shell_task( ntasks );
  std::vector<double> max_F_shells(nshells);
  size_t i_task = 0;
  for( auto task_it = task_begin; task_it != task_end; task_it++, i_task++ ) {
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

  } // Loop over tasks
}


}
