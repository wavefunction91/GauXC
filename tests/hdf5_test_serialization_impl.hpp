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

#include "collocation_common.hpp"
#include "hdf5_test_serialization.hpp"

namespace GauXC {

// Write ref_collocation_data to HDF5
inline void write_collocation_data(const std::vector<ref_collocation_data>& ref_data,
                                     const std::string& filename) {
  HighFive::File file(filename, HighFive::File::Truncate);
  
  // Write number of data entries
  {
    HighFive::DataSpace space(1);
    auto dset = file.createDataSet<size_t>("/nentries", space);
    size_t nentries = ref_data.size();
    dset.write_raw(&nentries);
  }
  
  // Write each collocation data entry
  for (size_t i = 0; i < ref_data.size(); ++i) {
    std::string group = "/entry_" + std::to_string(i);
    const auto& d = ref_data[i];
    
    // Write all the data fields
    write_int32_vector(d.mask, file, group + "/mask");
    write_points_array(d.pts, file, group + "/pts");
    write_double_vector(d.eval, file, group + "/eval");
    write_double_vector(d.deval_x, file, group + "/deval_x");
    write_double_vector(d.deval_y, file, group + "/deval_y");
    write_double_vector(d.deval_z, file, group + "/deval_z");
    write_double_vector(d.d2eval_xx, file, group + "/d2eval_xx");
    write_double_vector(d.d2eval_xy, file, group + "/d2eval_xy");
    write_double_vector(d.d2eval_xz, file, group + "/d2eval_xz");
    write_double_vector(d.d2eval_yy, file, group + "/d2eval_yy");
    write_double_vector(d.d2eval_yz, file, group + "/d2eval_yz");
    write_double_vector(d.d2eval_zz, file, group + "/d2eval_zz");
    write_double_vector(d.d2eval_lapl, file, group + "/d2eval_lapl");
    write_double_vector(d.d3eval_lapl_x, file, group + "/d3eval_lapl_x");
    write_double_vector(d.d3eval_lapl_y, file, group + "/d3eval_lapl_y");
    write_double_vector(d.d3eval_lapl_z, file, group + "/d3eval_lapl_z");
  }
}

// Read ref_collocation_data from HDF5
inline void read_collocation_data(std::vector<ref_collocation_data>& ref_data,
                                    const std::string& filename) {
  HighFive::File file(filename, HighFive::File::ReadOnly);
  
  // Read number of data entries
  size_t nentries;
  {
    auto dset = file.getDataSet("/nentries");
    dset.read(&nentries);
  }
  
  ref_data.resize(nentries);
  
  // Read each collocation data entry
  for (size_t i = 0; i < nentries; ++i) {
    std::string group = "/entry_" + std::to_string(i);
    auto& d = ref_data[i];
    
    // Read all the data fields
    read_int32_vector(d.mask, file, group + "/mask");
    read_points_array(d.pts, file, group + "/pts");
    read_double_vector(d.eval, file, group + "/eval");
    read_double_vector(d.deval_x, file, group + "/deval_x");
    read_double_vector(d.deval_y, file, group + "/deval_y");
    read_double_vector(d.deval_z, file, group + "/deval_z");
    read_double_vector(d.d2eval_xx, file, group + "/d2eval_xx");
    read_double_vector(d.d2eval_xy, file, group + "/d2eval_xy");
    read_double_vector(d.d2eval_xz, file, group + "/d2eval_xz");
    read_double_vector(d.d2eval_yy, file, group + "/d2eval_yy");
    read_double_vector(d.d2eval_yz, file, group + "/d2eval_yz");
    read_double_vector(d.d2eval_zz, file, group + "/d2eval_zz");
    read_double_vector(d.d2eval_lapl, file, group + "/d2eval_lapl");
    read_double_vector(d.d3eval_lapl_x, file, group + "/d3eval_lapl_x");
    read_double_vector(d.d3eval_lapl_y, file, group + "/d3eval_lapl_y");
    read_double_vector(d.d3eval_lapl_z, file, group + "/d3eval_lapl_z");
  }
}

// Write ref_weights_data to HDF5
inline void write_weights_data(const ref_weights_data& ref_data,
                                 const std::string& filename) {
  HighFive::File file(filename, HighFive::File::Truncate);
  
  // Write molecule
  write_hdf5_record(ref_data.mol, filename, "/MOLECULE");
  
  // Write tasks_unm
  write_xctask_vector(ref_data.tasks_unm, filename, "/tasks_unm");
  
  // Write tasks_mod
  write_xctask_vector(ref_data.tasks_mod, filename, "/tasks_mod");
}

// Read ref_weights_data from HDF5
inline void read_weights_data(ref_weights_data& ref_data,
                                const std::string& filename) {
  HighFive::File file(filename, HighFive::File::ReadOnly);
  
  // Read molecule
  read_hdf5_record(ref_data.mol, filename, "/MOLECULE");
  
  // Reconstruct meta from molecule
  ref_data.meta = std::make_shared<MolMeta>(ref_data.mol);
  
  // Read tasks_unm
  read_xctask_vector(ref_data.tasks_unm, filename, "/tasks_unm");
  
  // Read tasks_mod
  read_xctask_vector(ref_data.tasks_mod, filename, "/tasks_mod");
}

} // namespace GauXC
