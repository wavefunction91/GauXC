/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "host/reference/weights.hpp"
#include "common/integrator_constants.hpp"

#include <gauxc/molgrid/defaults.hpp>

namespace GauXC {

// Reference Becke weights impl
void reference_becke_weights_host(
  const Molecule&        mol,
  const MolMeta&         meta,
  task_iterator          task_begin,
  task_iterator          task_end
) {

  // Becke partition functions
  auto hBecke = [](double x) {return 1.5 * x - 0.5 * x * x * x;}; // Eq. 19
  auto gBecke = [&](double x) {return hBecke(hBecke(hBecke(x)));}; // Eq. 20 f_3

  const size_t ntasks = std::distance(task_begin,task_end);
  const size_t natoms = mol.natoms();

  const auto&  RAB    = meta.rab();

  std::vector<double> slater_radii;
  for( auto& atom : mol ) {
    slater_radii.emplace_back( default_atomic_radius(atom.Z) );
  }


#if 0
  // TODO: Add a pathway for this
  std::vector<double> size_adj(natoms * natoms);
  for( auto i = 0; i < natoms; ++i ) 
  for( auto j = 0; j < natoms; ++j ) {
    const auto si  = slater_radii[i];
    const auto sj  = slater_radii[j];
    const auto chi = std::sqrt(si/sj);
    const auto u   = (chi-1.)/(chi+1.);
    const auto a   = u / (u*u-1.);

    size_adj[i + j*natoms] = a;
  }
#endif

  #pragma omp parallel 
  {

  std::vector<double> partitionScratch( natoms );
  std::vector<double> atomDist( natoms );

  #pragma omp for
  for( size_t iT = 0; iT < ntasks;                  ++iT )
  for( size_t i  = 0; i  < (task_begin+iT)->points.size(); ++i  ) {

    auto&       task   = *(task_begin+iT);
    auto&       weight = task.weights[i];
    const auto& point  = task.points[i];

    // Compute distances of each center to point
    for(size_t iA = 0; iA < natoms; iA++) {

      const double da_x = point[0] - mol[iA].x;
      const double da_y = point[1] - mol[iA].y;
      const double da_z = point[2] - mol[iA].z;

      atomDist[iA] = std::sqrt(da_x*da_x + da_y*da_y + da_z*da_z);

    }

    // Evaluate unnormalized partition functions 
    std::fill(partitionScratch.begin(),partitionScratch.end(),1.);
    for( size_t iA = 0; iA < natoms; iA++ ) 
    for( size_t jA = 0; jA < iA;     jA++ ){

      double mu  = (atomDist[iA] - atomDist[jA]) / RAB[jA + iA*natoms];

#if 0
      // Size Adjustment
      const double a = size_adj[iA + jA*natoms];
      mu = mu + a * ( 1. - mu*mu );
#endif

      const double g = gBecke(mu);

      partitionScratch[iA] *= 0.5 * (1. - g);
      partitionScratch[jA] *= 0.5 * (1. + g);
    }

    // Normalization
    double sum = 0.;
    for( size_t iA = 0; iA < natoms; iA++ )  sum += partitionScratch[iA];

    // Update Weights
    weight *= partitionScratch[task.iParent] / sum;

  } // Collapsed loop over tasks and points

  } // OMP context


}

void reference_ssf_weights_host(
  const Molecule&        mol,
  const MolMeta&         meta,
  task_iterator          task_begin,
  task_iterator          task_end
) {

  auto gFrisch = [&](double x) {
    const double s_x  = x / integrator::magic_ssf_factor<>;
    const double s_x2 = s_x  * s_x;
    const double s_x3 = s_x  * s_x2;
    const double s_x5 = s_x3 * s_x2;
    const double s_x7 = s_x5 * s_x2;

    return (35.*(s_x - s_x3) + 21.*s_x5 - 5.*s_x7) / 16.;
  };

  const size_t ntasks = std::distance(task_begin,task_end);
  const size_t natoms = mol.natoms();

  const auto&  RAB    = meta.rab();

  #pragma omp parallel 
  {

  std::vector<double> partitionScratch( natoms );
  std::vector<double> atomDist( natoms );

  #pragma omp for
  for( size_t iT = 0; iT < ntasks;                  ++iT )
  for( size_t i  = 0; i  < (task_begin+iT)->points.size(); ++i  ) {

    auto&       task   = *(task_begin+iT);
    auto&       weight = task.weights[i];
    const auto& point  = task.points[i];

    const auto dist_cutoff = 0.5 * (1-integrator::magic_ssf_factor<>) * task.dist_nearest;

    // Compute dist to parent atom
    {
      const double da_x = point[0] - mol[task.iParent].x;
      const double da_y = point[1] - mol[task.iParent].y;
      const double da_z = point[2] - mol[task.iParent].z;

      atomDist[task.iParent] = std::sqrt(da_x*da_x + da_y*da_y + da_z*da_z);
    }

    if( atomDist[task.iParent] < dist_cutoff ) continue; // Partition weight = 1

    // Compute distances of each center to point
    for(size_t iA = 0; iA < natoms; iA++) {

      if( iA == (size_t)task.iParent ) continue;

      const double da_x = point[0] - mol[iA].x;
      const double da_y = point[1] - mol[iA].y;
      const double da_z = point[2] - mol[iA].z;

      atomDist[iA] = std::sqrt(da_x*da_x + da_y*da_y + da_z*da_z);

    }

    // Evaluate unnormalized partition functions 
    std::fill(partitionScratch.begin(),partitionScratch.end(),1.);
#if 1
    for( size_t iA = 0; iA < natoms; iA++ ) 
    for( size_t jA = 0; jA < iA;     jA++ )
    if( partitionScratch[iA] > integrator::ssf_weight_tol or 
        partitionScratch[jA] > integrator::ssf_weight_tol ) {

      const double mu = (atomDist[iA] - atomDist[jA]) / RAB[jA + iA*natoms];

      if( mu <= -integrator::magic_ssf_factor<> ) {

        partitionScratch[jA] = 0.;

      } else if (mu >= integrator::magic_ssf_factor<>) {

        partitionScratch[iA] = 0.;

      } else {

        double g = 0.5 * ( 1. - gFrisch(mu) );
        partitionScratch[iA] *= g;
        partitionScratch[jA] *= 1. - g;

      }

    }
#else
    for(size_t iA = 0; iA < natoms; ++iA)
    for(size_t jA = 0; jA < natoms; ++jA) 
    if(iA != jA and partitionScratch[iA] > integrator::ssf_weight_tol) {
      const double mu = (atomDist[iA] - atomDist[jA]) / RAB[jA + iA*natoms];
      if( fabs(mu) < integrator::magic_ssf_factor<> ) {
        double g = 0.5 * (1. - gFrisch(mu));
        partitionScratch[iA] *= g;
      } else if(mu >= integrator::magic_ssf_factor<>) {
        partitionScratch[iA] = 0.0;
      }
    }

    if(partitionScratch[task.iParent] < std::numeric_limits<double>::epsilon()) {
      weight = 0;
      continue;
    }
#endif

    // Normalization
    double sum = 0.;
    for( size_t iA = 0; iA < natoms; iA++ )  sum += partitionScratch[iA];

    // Update Weights
    weight *= partitionScratch[task.iParent] / sum;

  } // Collapsed loop over tasks and points

  } // OMP context


}

void reference_lko_weights_host(
  const Molecule&        mol,
  const MolMeta&         meta,
  task_iterator          task_begin,
  task_iterator          task_end
) {


  // Sort on atom index
  std::stable_sort( task_begin, task_end, 
    [](const auto& a, const auto&b ) { return a.iParent < b.iParent; } );

  // Becke partition functions
  auto hBecke = [](double x) {return 1.5 * x - 0.5 * x * x * x;}; // Eq. 19
  auto gBecke = [&](double x) {return hBecke(hBecke(hBecke(x)));}; // Eq. 20 f_3
  constexpr double R_cutoff = 5;

  const size_t natoms = mol.natoms();

  const auto&  RAB    = meta.rab();

  #pragma omp parallel 
  {

  std::vector<double> partitionScratch( natoms );
  std::vector<double> atomDist( natoms );
  std::vector<size_t> inter_atom_dist_idx( natoms );
  std::vector<size_t> point_dist_idx( natoms );

  #pragma omp for schedule(dynamic)
  for( auto iAtom = 0ul; iAtom < natoms; ++iAtom ) {

    auto atom_begin = std::find_if( task_begin, task_end,
      [&](const auto& t){ return t.iParent == (int)iAtom; } );
    auto atom_end = std::find_if( task_begin, task_end,
      [&](const auto& t){ return t.iParent == (int)(iAtom+1); } );

    auto* RAB_parent = RAB.data() + iAtom*natoms;

    std::iota( inter_atom_dist_idx.begin(), inter_atom_dist_idx.end(), 0 );
    std::sort( inter_atom_dist_idx.begin(), inter_atom_dist_idx.end(),
      [&](auto i, auto j){ return RAB_parent[i] < RAB_parent[j]; } );

  for( auto task_it = atom_begin; task_it != atom_end; ++task_it ) {

    auto& points  = task_it->points;
    auto& weights = task_it->weights;
    const auto npts = points.size();

  for( auto ipt = 0ul; ipt < npts; ++ipt ) {

    auto& weight = weights[ipt];
    const auto point = points[ipt];

    std::fill( atomDist.begin(), atomDist.end(), std::numeric_limits<double>::infinity() );
    // Parent distance
    {
      const double da_x = point[0] - mol[iAtom].x;
      const double da_y = point[1] - mol[iAtom].y;
      const double da_z = point[2] - mol[iAtom].z;

      atomDist[iAtom] = std::sqrt(da_x*da_x + da_y*da_y + da_z*da_z);
    }

    double r_parent  = atomDist[iAtom];
    double r_nearest = r_parent;
    size_t natoms_keep = 1;
    // Compute distances of each center to point
    for(size_t iA = 1; iA < natoms; iA++) {
      auto idx = inter_atom_dist_idx[iA];
      if( RAB_parent[idx] > (r_parent + r_nearest + 2*R_cutoff) ) break;

      const double da_x = point[0] - mol[idx].x;
      const double da_y = point[1] - mol[idx].y;
      const double da_z = point[2] - mol[idx].z;

      const auto r = std::sqrt(da_x*da_x + da_y*da_y + da_z*da_z);
      r_nearest = std::min( r_nearest, r );
      atomDist[idx] = r;
      ++natoms_keep;
    }

     // Partition weight is 0
    if( r_parent > r_nearest + R_cutoff ) {
      weight = 0.;
      continue;
    }

    // Partiton atom indices into a petite list of non-negligible centers
    std::iota( point_dist_idx.begin(), point_dist_idx.end(), 0 );
    auto atom_keep_end = std::partition( point_dist_idx.begin(), point_dist_idx.end(), 
      [&](auto i){ return atomDist[i] < std::numeric_limits<double>::infinity(); } );

    // Only sort over non-negligible cetners
    std::sort( point_dist_idx.begin(), atom_keep_end,
      [&](auto i, auto j){ return atomDist[i] < atomDist[j]; } );

    // Get parent index
    auto parent_it  = std::find( point_dist_idx.begin(), atom_keep_end, iAtom );
    auto parent_idx = std::distance( point_dist_idx.begin(), parent_it );

    // Sort atom distances for contiguous reads in weight loop
    auto atom_dist_end = std::partition( atomDist.begin(), atomDist.end(),
      [](auto x){ return x < std::numeric_limits<double>::infinity(); } );
    std::sort( atomDist.begin(), atom_dist_end );



    // Evaluate unnormalized partition functions 
    std::fill_n(partitionScratch.begin(),natoms_keep,0.);
    for( auto i = 0ul; i < natoms_keep; ++i ) {
      auto idx_i = point_dist_idx[i];
      auto r_i = atomDist[i];
      if( r_i > (r_nearest + R_cutoff) ) { break; }
      partitionScratch[i] = 1.;

      const auto* RAB_i_idx = RAB.data() + idx_i*natoms;

    for( auto j = 0ul; j < i; ++j ) {
      auto idx_j = point_dist_idx[j];
      auto r_j = atomDist[j];
      if( r_j > (r_i + R_cutoff) ) { break; }

      const double mu = 
        (r_i - r_j) / std::min(RAB_i_idx[idx_j], R_cutoff);

      const double g = gBecke(mu);
      const auto   s_ij = 0.5 * (1. - g);
      partitionScratch[i] *= s_ij;
      partitionScratch[j] *= 1. - s_ij;
      
    }
    }

    // Normalization
    double sum = 0.;
    for( size_t iA = 0; iA < natoms_keep; iA++ )  sum += partitionScratch[iA];

    // Update Weights
    weight *= partitionScratch[parent_idx] / sum;

  } // Loop over points 
  } // Loop over tasks
  } // Loop over atoms

  } // OMP context

}

}
