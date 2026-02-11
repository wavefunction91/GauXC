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
#pragma once

#include <gauxc/molecular_weights.h>
#include <gauxc/molecular_weights.hpp>

namespace GauXC::detail {
static inline MolecularWeights* get_molecular_weights_ptr(C::GauXCMolecularWeights mw) noexcept {
  return static_cast<MolecularWeights*>(mw.ptr);
}
static inline std::shared_ptr<MolecularWeights>* get_molecular_weights_shared(C::GauXCMolecularWeights mw) noexcept {
  return static_cast<std::shared_ptr<MolecularWeights>*>(mw.ptr);
}
static inline MolecularWeightsFactory* get_molecular_weights_factory_ptr(C::GauXCMolecularWeightsFactory mwf) noexcept {
  return static_cast<MolecularWeightsFactory*>(mwf.ptr);
}
static inline MolecularWeightsSettings
convert_molecular_weights_settings( const GauXC::C::GauXCMolecularWeightsSettings& c_settings ) {
  MolecularWeightsSettings settings;
  settings.weight_alg = static_cast<XCWeightAlg>( c_settings.weight_alg );
  settings.becke_size_adjustment = c_settings.becke_size_adjustment;
  return settings;
}
}