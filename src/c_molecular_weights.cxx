/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include <gauxc/molecular_weights.h>
#include <gauxc/molecular_weights.hpp>
#include <gauxc/load_balancer.h>
#include <gauxc/load_balancer.hpp>
#include <gauxc/util/c_molecular_weights.hpp>
#include <gauxc/util/c_load_balancer.hpp>
#include <gauxc/util/c_status.hpp>

namespace GauXC::C {
extern "C" {

GauXCMolecularWeightsFactory gauxc_molecular_weights_factory_new(
  GauXCStatus* status,
  enum GauXC_ExecutionSpace ex,
  const char* kernel_name,
  const GauXCMolecularWeightsSettings settings
) {
  GauXCMolecularWeightsFactory mwf{};
  mwf.hdr = GauXCHeader{GauXC_Type_MolecularWeightsFactory};
  mwf.ptr = nullptr;

  try {
    mwf.ptr = new MolecularWeightsFactory(
      static_cast<ExecutionSpace>(ex),
      std::string(kernel_name),
      detail::convert_molecular_weights_settings( settings )
    );
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return mwf;
}

GauXCMolecularWeights gauxc_molecular_weights_factory_get_instance(
  GauXCStatus* status,
  const GauXCMolecularWeightsFactory mwf
) {
  GauXCMolecularWeights mw{};
  mw.hdr = GauXCHeader{GauXC_Type_MolecularWeights};
  mw.ptr = nullptr;
  mw.owned = true;

  try {
    MolecularWeights mw_instance = detail::get_molecular_weights_factory_ptr(mwf)->get_instance();
    mw.ptr = new MolecularWeights( std::move(mw_instance) );
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return mw;
}

GauXCMolecularWeights gauxc_molecular_weights_factory_get_shared_instance(
  GauXCStatus* status,
  const GauXCMolecularWeightsFactory mwf
) {
  GauXCMolecularWeights mw{};
  mw.hdr = GauXCHeader{GauXC_Type_MolecularWeights};
  mw.ptr = nullptr;
  mw.owned = false;

  try {
    auto mw_instance_ptr = detail::get_molecular_weights_factory_ptr(mwf)->get_shared_instance();
    mw.ptr = new std::shared_ptr<MolecularWeights>( std::move(mw_instance_ptr) );
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return mw;
}

void gauxc_molecular_weights_modify_weights(
  GauXCStatus* status,
  const GauXCMolecularWeights mw,
  const GauXCLoadBalancer lb
) {
  try {
    if (mw.owned) {
      if (lb.owned)
        detail::get_molecular_weights_ptr(mw)->modify_weights(
          *detail::get_load_balancer_ptr(lb)
        );
      else
        detail::get_molecular_weights_ptr(mw)->modify_weights(
          **detail::get_load_balancer_shared(lb)
        );
    } else {
      if (lb.owned)
        detail::get_molecular_weights_shared(mw)->get()->modify_weights(
          *detail::get_load_balancer_ptr(lb)
        );
      else
        detail::get_molecular_weights_shared(mw)->get()->modify_weights(
          **detail::get_load_balancer_shared(lb)
        );
    }
    status->code = 0;
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
}


void gauxc_molecular_weights_delete(
  GauXCStatus* status,
  GauXCMolecularWeights* mw
) {
  status->code = 0;
  if(mw == nullptr) return;
  if(mw->ptr != nullptr) {
    if (mw->owned)
      delete detail::get_molecular_weights_ptr(*mw);
    else
      delete detail::get_molecular_weights_shared(*mw);
  }
  mw->ptr = nullptr;
}


void gauxc_molecular_weights_factory_delete(
  GauXCStatus* status,
  GauXCMolecularWeightsFactory* mwf
) {
  status->code = 0;
  if(mwf == nullptr) return;
  if(mwf->ptr != nullptr)
    delete detail::get_molecular_weights_factory_ptr(*mwf);
  mwf->ptr = nullptr;
}

} // extern "C"
} // namespace GauXC::C