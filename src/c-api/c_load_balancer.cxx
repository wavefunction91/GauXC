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
#include <gauxc/c/load_balancer.h>

#include <gauxc/load_balancer.hpp>

#include "c_molecule.hpp"
#include "c_basisset.hpp"
#include "c_molgrid.hpp"
#include "c_runtime_environment.hpp"
#include "c_load_balancer.hpp"
#include "c_status.hpp"

namespace GauXC::C {
extern "C" {

void gauxc_load_balancer_delete(
  GauXCStatus* status,
  GauXCLoadBalancer* lb
) {
  detail::gauxc_status_init(status);
  if(lb == nullptr) return;
  if(lb->ptr != nullptr) {
    delete detail::get_load_balancer_ptr(*lb);
  }
  lb->ptr = nullptr;
}

GauXCLoadBalancerFactory gauxc_load_balancer_factory_new(
  GauXCStatus* status,
  enum GauXC_ExecutionSpace ex,
  const char* kernel_name
) {
  detail::gauxc_status_init(status);
  GauXCLoadBalancerFactory lbf{};
  lbf.hdr = GauXCHeader{GauXC_Type_LoadBalancerFactory};
  lbf.ptr = nullptr;
  if (kernel_name == nullptr) {
    detail::gauxc_status_handle(status, 1, "Kernel name string cannot be null");
    return lbf;
  }

  try {
    lbf.ptr = new LoadBalancerFactory(
      static_cast<ExecutionSpace>(ex),
      std::string(kernel_name)
    );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return lbf;
}

void gauxc_load_balancer_factory_delete(
  GauXCStatus* status,
  GauXCLoadBalancerFactory* lbf
) {
  detail::gauxc_status_init(status);
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
  detail::gauxc_status_init(status);
  GauXCLoadBalancer lb{};
  lb.hdr = GauXCHeader{GauXC_Type_LoadBalancer};
  lb.ptr = nullptr;
  if (lbf.ptr == nullptr || lbf.hdr.type != GauXC_Type_LoadBalancerFactory) {
    detail::gauxc_status_handle(status, 1, "Invalid LoadBalancerFactory handle");
    return lb;
  }
  if (rt.hdr.type != GauXC_Type_RuntimeEnvironment) {
    detail::gauxc_status_handle(status, 1, "Invalid RuntimeEnvironment handle");
    return lb;
  }
  if (mol.ptr == nullptr || mol.hdr.type != GauXC_Type_Molecule) {
    detail::gauxc_status_handle(status, 1, "Invalid Molecule handle");
    return lb;
  }
  if (mg.ptr == nullptr || mg.hdr.type != GauXC_Type_MolGrid) {
    detail::gauxc_status_handle(status, 1, "Invalid MolGrid handle");
    return lb;
  }
  if (bs.ptr == nullptr || bs.hdr.type != GauXC_Type_BasisSet) {
    detail::gauxc_status_handle(status, 1, "Invalid BasisSet handle");
    return lb;
  }

  try {
    auto lb_instance_ptr = detail::get_load_balancer_factory_ptr(lbf)->get_shared_instance(
      *detail::get_runtime_environment_ptr(rt),
      *detail::get_molecule_ptr(mol),
      *detail::get_molgrid_ptr(mg),
      *detail::get_basisset_ptr(bs)
    );
    lb.ptr = new std::shared_ptr<LoadBalancer>( std::move(lb_instance_ptr) );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return lb;
}

} // extern "C"
} // namespace GauXC::C