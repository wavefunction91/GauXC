#pragma once
#include <gauxc/molgrid.hpp>
#include <gauxc/load_balancer.hpp>


using namespace GauXC;

struct ref_weights_data {
  Molecule                  mol;
  std::shared_ptr<MolMeta>  meta;
  std::vector< XCTask > tasks_unm;
  std::vector< XCTask > tasks_mod; // This is only the weights

  template <typename Archive>
  void load( Archive& ar ) {
    ar( mol, tasks_unm, tasks_mod );
    meta = std::make_shared<MolMeta>(mol);
  }
  template <typename Archive>
  void save( Archive& ar ) const {
    ar( mol, tasks_unm, tasks_mod );
  }
};


#ifdef GAUXC_ENABLE_HOST
#include "host/reference/weights.hpp"

void generate_weights_data( const Molecule& mol, const BasisSet<double>& basis,
                                std::ofstream& out_file, size_t ntask_save = 15 ) {


  MolGrid mg(AtomicGridSizeDefault::FineGrid, mol);
  LoadBalancerFactory lb_factory(ExecutionSpace::Host, "Default");
#ifdef GAUXC_ENABLE_MPI
  auto lb = lb_factory.get_instance(MPI_COMM_WORLD, mol, mg, basis);
#else
  auto lb = lb_factory.get_instance(mol, mg, basis);
#endif
  auto& tasks = lb.get_tasks();

  ref_weights_data   ref_data;
  ref_data.mol       = mol;

  auto abs_comparator = []( const auto& a, const auto& b ) {
    return std::abs(a) < std::abs(b);
  };

  std::sort( tasks.begin(), tasks.end(),
    [&]( const auto& a, const auto& b ) {
      auto a_max =
        *std::max_element( a.weights.begin(), a.weights.end(),
                           abs_comparator );
      auto b_max =
        *std::max_element( b.weights.begin(), b.weights.end(),
                           abs_comparator );

      return a_max < b_max;
    });

  if( tasks.size() > ntask_save )
    tasks.erase( tasks.begin() + ntask_save, tasks.end() );

  ref_data.tasks_unm = tasks; // Make a copy of un modified tasks

  reference_ssf_weights_host( 
    mol, lb.molmeta(), tasks.begin(), tasks.end() );

  // Clear out unneeded data
  for( auto& task : tasks ) {
    task.points.clear();
    task.shell_list.clear();
  }
  ref_data.tasks_mod = tasks;

  {
    cereal::BinaryOutputArchive ar( out_file );
    ar( ref_data );
  }

}
#endif
