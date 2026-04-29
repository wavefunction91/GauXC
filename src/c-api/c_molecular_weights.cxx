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
#include <gauxc/c/molecular_weights.h>
#include <gauxc/c/load_balancer.h>

#include <gauxc/molecular_weights.hpp>
#include <gauxc/load_balancer.hpp>

#include "c_molecular_weights.hpp"
#include "c_load_balancer.hpp"
#include "c_status.hpp"

namespace GauXC::C {
extern "C" {

GauXCMolecularWeightsFactory gauxc_molecular_weights_factory_new(
  GauXCStatus* status,
  enum GauXC_ExecutionSpace ex,
  const char* kernel_name,
  const GauXCMolecularWeightsSettings settings
) {
  detail::gauxc_status_init(status);
  GauXCMolecularWeightsFactory mwf{};
  mwf.hdr = GauXCHeader{GauXC_Type_MolecularWeightsFactory};
  mwf.ptr = nullptr;
  if (kernel_name == nullptr) {
    detail::gauxc_status_handle(status, 1, "Kernel name string cannot be null");
    return mwf;
  }

  try {
    mwf.ptr = new MolecularWeightsFactory(
      static_cast<ExecutionSpace>(ex),
      std::string(kernel_name),
      detail::convert_molecular_weights_settings( settings )
    );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return mwf;
}

GauXCMolecularWeights gauxc_molecular_weights_factory_get_instance(
  GauXCStatus* status,
  const GauXCMolecularWeightsFactory mwf
) {
  detail::gauxc_status_init(status);
  GauXCMolecularWeights mw{};
  mw.hdr = GauXCHeader{GauXC_Type_MolecularWeights};
  mw.ptr = nullptr;
  if (mwf.ptr == nullptr || mwf.hdr.type != GauXC_Type_MolecularWeightsFactory) {
    detail::gauxc_status_handle(status, 1, "Invalid MolecularWeightsFactory handle");
    return mw;
  }

  try {
    MolecularWeights mw_instance = detail::get_molecular_weights_factory_ptr(mwf)->get_instance();
    mw.ptr = new MolecularWeights( std::move(mw_instance) );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return mw;
}

void gauxc_molecular_weights_modify_weights(
  GauXCStatus* status,
  const GauXCMolecularWeights mw,
  const GauXCLoadBalancer lb
) {
  detail::gauxc_status_init(status);
  if (mw.ptr == nullptr || mw.hdr.type != GauXC_Type_MolecularWeights) {
    detail::gauxc_status_handle(status, 1, "Invalid MolecularWeights handle");
    return;
  }
  if (lb.ptr == nullptr || lb.hdr.type != GauXC_Type_LoadBalancer) {
    detail::gauxc_status_handle(status, 1, "Invalid LoadBalancer handle");
    return;
  }
  try {
    detail::get_molecular_weights_ptr(mw)->modify_weights(
      **detail::get_load_balancer_ptr(lb)
    );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}


void gauxc_molecular_weights_delete(
  GauXCStatus* status,
  GauXCMolecularWeights* mw
) {
  detail::gauxc_status_init(status);
  if(mw == nullptr) return;
  if(mw->ptr != nullptr) {
    delete detail::get_molecular_weights_ptr(*mw);
  }
  mw->ptr = nullptr;
}


void gauxc_molecular_weights_factory_delete(
  GauXCStatus* status,
  GauXCMolecularWeightsFactory* mwf
) {
  detail::gauxc_status_init(status);
  if(mwf == nullptr) return;
  if(mwf->ptr != nullptr)
    delete detail::get_molecular_weights_factory_ptr(*mwf);
  mwf->ptr = nullptr;
}

} // extern "C"
} // namespace GauXC::C