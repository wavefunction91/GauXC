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
#include <gauxc/load_balancer.h>
#include <gauxc/load_balancer.hpp>
#include <gauxc/util/c_molecule.hpp>
#include <gauxc/util/c_basisset.hpp>
#include <gauxc/util/c_molgrid.hpp>
#include <gauxc/util/c_runtime_environment.hpp>
#include <gauxc/util/c_load_balancer.hpp>
#include <gauxc/util/c_status.hpp>

namespace GauXC::C {
extern "C" {

void gauxc_load_balancer_delete(
  GauXCStatus* status,
  GauXCLoadBalancer* lb
) {
  status->code = 0;
  if(lb == nullptr) return;
  if(lb->ptr != nullptr) {
    if (lb->owned)
      delete detail::get_load_balancer_ptr(*lb);
    else
      delete detail::get_load_balancer_shared(*lb);
  }
  lb->ptr = nullptr;
}

GauXCLoadBalancerFactory gauxc_load_balancer_factory_new(
  GauXCStatus* status,
  enum GauXC_ExecutionSpace ex,
  const char* kernel_name
) {
  status->code = 0;
  GauXCLoadBalancerFactory lbf{};
  lbf.hdr = GauXCHeader{GauXC_Type_LoadBalancerFactory};
  lbf.ptr = nullptr;

  try {
    lbf.ptr = new LoadBalancerFactory(
      static_cast<ExecutionSpace>(ex),
      std::string(kernel_name)
    );
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return lbf;
}

void gauxc_load_balancer_factory_delete(
  GauXCStatus* status,
  GauXCLoadBalancerFactory* lbf
) {
  status->code = 0;
  if(lbf == nullptr) return;
  if(lbf->ptr != nullptr)
    delete detail::get_load_balancer_factory_ptr(*lbf);
  lbf->ptr = nullptr;
}

GauXCLoadBalancer gauxc_load_balancer_factory_get_instance(
  GauXCStatus* status,
  const GauXCLoadBalancerFactory lbf,
  const GauXCRuntimeEnvironment rt,
  const GauXCMolecule mol,
  const GauXCMolGrid mg,
  const GauXCBasisSet bs
) {
  status->code = 0;
  GauXCLoadBalancer lb{};
  lb.hdr = GauXCHeader{GauXC_Type_LoadBalancer};
  lb.ptr = nullptr;
  lb.owned = true;

  try {
    LoadBalancer lb_instance = detail::get_load_balancer_factory_ptr(lbf)->get_instance(
      *detail::get_runtime_environment_ptr(rt),
      *detail::get_molecule_ptr(mol),
      *detail::get_molgrid_ptr(mg),
      *detail::get_basisset_ptr(bs)
    );
    lb.ptr = new LoadBalancer( std::move(lb_instance) );
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return lb;
}

GauXCLoadBalancer gauxc_load_balancer_factory_get_shared_instance(
  GauXCStatus* status,
  const GauXCLoadBalancerFactory lbf,
  const GauXCRuntimeEnvironment rt,
  const GauXCMolecule mol,
  const GauXCMolGrid mg,
  const GauXCBasisSet bs
) {
  status->code = 0;
  GauXCLoadBalancer lb{};
  lb.hdr = GauXCHeader{GauXC_Type_LoadBalancer};
  lb.ptr = nullptr;
  lb.owned = false;

  try {
    auto lb_instance_ptr = detail::get_load_balancer_factory_ptr(lbf)->get_shared_instance(
      *detail::get_runtime_environment_ptr(rt),
      *detail::get_molecule_ptr(mol),
      *detail::get_molgrid_ptr(mg),
      *detail::get_basisset_ptr(bs)
    );
    lb.ptr = new std::shared_ptr<LoadBalancer>( std::move(lb_instance_ptr) );
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
  return lb;
}

} // extern "C"
} // namespace GauXC::C