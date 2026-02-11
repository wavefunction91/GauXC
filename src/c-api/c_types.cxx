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

#include <gauxc/types.h>
#include <gauxc/status.h>
#include <gauxc/molecule.h>
#include <gauxc/basisset.h>
#include <gauxc/molgrid.h>
#include <gauxc/runtime_environment.h>
#include <gauxc/load_balancer.h>
#include <gauxc/molecular_weights.h>
#include <gauxc/functional.h>
#include <gauxc/xc_integrator.h>
#include <gauxc/matrix.h>

#include "c_status.hpp"

namespace GauXC::C {
extern "C" {

void gauxc_object_delete(GauXCStatus* status, void** obj) {
   detail::gauxc_status_init(status);
   if(obj == nullptr) return;

   struct GauXCObject {
      GauXCHeader hdr;
   };

   GauXCHeader* header = &reinterpret_cast<struct GauXCObject*>(*obj)->hdr;

   switch(header->type) {
     case GauXC_Type_Molecule: {
       GauXCMolecule* mol = reinterpret_cast<GauXCMolecule*>(*obj);
       gauxc_molecule_delete(status, mol);
       break;
     }
     case GauXC_Type_BasisSet: {
       GauXCBasisSet* basis = reinterpret_cast<GauXCBasisSet*>(*obj);
       gauxc_basisset_delete(status, basis);
       break;
     }
     case GauXC_Type_MolGrid: {
       GauXCMolGrid* mg = reinterpret_cast<GauXCMolGrid*>(*obj);
       gauxc_molgrid_delete(status, mg);
       break;
     }
     case GauXC_Type_RuntimeEnvironment: {
       GauXCRuntimeEnvironment* re = reinterpret_cast<GauXCRuntimeEnvironment*>(*obj);
       gauxc_runtime_environment_delete(status, re);
       break;
     }
     case GauXC_Type_LoadBalancer: {
       GauXCLoadBalancer* lb = reinterpret_cast<GauXCLoadBalancer*>(*obj);
       gauxc_load_balancer_delete(status, lb);
       break;
     }
     case GauXC_Type_LoadBalancerFactory: {
       GauXCLoadBalancerFactory* lbf = reinterpret_cast<GauXCLoadBalancerFactory*>(*obj);
       gauxc_load_balancer_factory_delete(status, lbf);
       break;
     }
     case GauXC_Type_MolecularWeights: {
       GauXCMolecularWeights* mw = reinterpret_cast<GauXCMolecularWeights*>(*obj);
       gauxc_molecular_weights_delete(status, mw);
       break;
     }
     case GauXC_Type_MolecularWeightsFactory: {
       GauXCMolecularWeightsFactory* mwf = reinterpret_cast<GauXCMolecularWeightsFactory*>(*obj);
       gauxc_molecular_weights_factory_delete(status, mwf);
       break;
     }
     case GauXC_Type_Functional: {
       GauXCFunctional* func = reinterpret_cast<GauXCFunctional*>(*obj);
       gauxc_functional_delete(status, func);
       break;
     }
     case GauXC_Type_Integrator: {
       GauXCIntegrator* integrator = reinterpret_cast<GauXCIntegrator*>(*obj);
       gauxc_integrator_delete(status, integrator);
       break;
     }
     case GauXC_Type_IntegratorFactory: { 
       GauXCIntegratorFactory* integrator_factory = reinterpret_cast<GauXCIntegratorFactory*>(*obj);
       gauxc_integrator_factory_delete(status, integrator_factory);
       break;
     }
     case GauXC_Type_Matrix: {
       GauXCMatrix* matrix = reinterpret_cast<GauXCMatrix*>(*obj);
       gauxc_matrix_delete(status, matrix);
       break;
     }
     default: {
       detail::gauxc_status_handle(status, 1, "Unknown object type in gauxc_object_delete");
       break;
     }
   }
}

void gauxc_objects_delete(
  GauXCStatus* status,
  void** ptrs,
  size_t nptrs
) {
   detail::gauxc_status_init(status);
   for(void** ptr = ptrs; ptr < ptrs + nptrs; ++ptr) {
      if(*ptr != nullptr) {
         gauxc_object_delete(status, ptr);
         if(status->code != 0) return;
      }
   }
}

} // extern "C"
} // namespace GauXC::C