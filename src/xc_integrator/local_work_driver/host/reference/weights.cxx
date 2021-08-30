#include "host/reference/weights.hpp"
#include "common/integrator_constants.hpp"

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
      const double mu = (atomDist[iA] - atomDist[jA]) / RAB[jA + iA*natoms];
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

    // Normalization
    double sum = 0.;
    for( size_t iA = 0; iA < natoms; iA++ )  sum += partitionScratch[iA];

    // Update Weights
    weight *= partitionScratch[task.iParent] / sum;

  } // Collapsed loop over tasks and points

  } // OMP context


}

}
