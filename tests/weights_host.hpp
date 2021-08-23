#pragma once
#include "weights_generate.hpp"
#include <fstream>
#include <string>

#ifdef GAUXC_ENABLE_HOST
#include "host/reference/weights.hpp"
using namespace GauXC;

void test_host_weights( std::ifstream& in_file ) {

  ref_weights_data ref_data;
  {
    cereal::BinaryInputArchive ar( in_file );
    ar( ref_data );
  }

  reference_ssf_weights_host( 
    ref_data.mol, *ref_data.meta, ref_data.tasks_unm.begin(), 
    ref_data.tasks_unm.end() );


  size_t ntasks = ref_data.tasks_unm.size();
  for( size_t itask = 0; itask < ntasks; ++itask ) {
    auto& task     = ref_data.tasks_unm.at(itask);
    auto& ref_task = ref_data.tasks_mod.at(itask);

    size_t npts = task.weights.size();
    for( size_t i = 0; i < npts; ++i ) {
      CHECK( task.weights.at(i) ==
             Approx(ref_task.weights.at(i)) );
    }
  }

}
#endif
