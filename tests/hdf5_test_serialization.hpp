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

#include <gauxc/external/hdf5.hpp>
#include <highfive/H5File.hpp>
#include <gauxc/xc_task.hpp>
#include <gauxc/molecule.hpp>
#include <gauxc/molmeta.hpp>
#include <vector>
#include <array>
#include <string>

namespace GauXC {

// Forward declarations removed - actual definitions are in collocation_common.hpp
// which must be included before this file or hdf5_test_serialization_impl.hpp

// Helper function to write std::array<double,3> vectors as 2D arrays
inline void write_points_array(const std::vector<std::array<double,3>>& pts,
                                 HighFive::File& file, const std::string& dset_name) {
  if (pts.empty()) return;
  
  // Convert to 2D vector for HighFive
  std::vector<std::vector<double>> pts_2d(pts.size(), std::vector<double>(3));
  for (size_t i = 0; i < pts.size(); ++i) {
    pts_2d[i][0] = pts[i][0];
    pts_2d[i][1] = pts[i][1];
    pts_2d[i][2] = pts[i][2];
  }
  
  auto dset = file.createDataSet<double>(dset_name, HighFive::DataSpace::From(pts_2d));
  dset.write(pts_2d);
}

// Helper function to read std::array<double,3> vectors from 2D arrays
inline void read_points_array(std::vector<std::array<double,3>>& pts,
                                HighFive::File& file, const std::string& dset_name) {
  // Check if dataset exists
  if (!file.exist(dset_name)) {
    pts.clear();
    return;
  }
  
  auto dset = file.getDataSet(dset_name);
  std::vector<std::vector<double>> pts_2d;
  dset.read(pts_2d);
  
  pts.resize(pts_2d.size());
  for (size_t i = 0; i < pts_2d.size(); ++i) {
    pts[i][0] = pts_2d[i][0];
    pts[i][1] = pts_2d[i][1];
    pts[i][2] = pts_2d[i][2];
  }
}

// Helper to write vector<double> with proper dimensionality
inline void write_double_vector(const std::vector<double>& vec,
                                  HighFive::File& file, const std::string& dset_name) {
  if (vec.empty()) {
    // Write empty dataset
    HighFive::DataSpace space({0});
    auto dset = file.createDataSet<double>(dset_name, space);
    return;
  }
  
  HighFive::DataSpace space(vec.size());
  auto dset = file.createDataSet<double>(dset_name, space);
  dset.write(vec);
}

// Helper to read vector<double>
inline void read_double_vector(std::vector<double>& vec,
                                 HighFive::File& file, const std::string& dset_name) {
  auto dset = file.getDataSet(dset_name);
  auto dims = dset.getDimensions();
  
  if (dims.empty() || dims[0] == 0) {
    vec.clear();
    return;
  }
  
  vec.resize(dims[0]);
  dset.read(vec);
}

// Helper to write vector<int32_t>
inline void write_int32_vector(const std::vector<int32_t>& vec,
                                 HighFive::File& file, const std::string& dset_name) {
  if (vec.empty()) {
    HighFive::DataSpace space({0});
    auto dset = file.createDataSet<int32_t>(dset_name, space);
    return;
  }
  
  HighFive::DataSpace space(vec.size());
  auto dset = file.createDataSet<int32_t>(dset_name, space);
  dset.write(vec);
}

// Helper to read vector<int32_t>
inline void read_int32_vector(std::vector<int32_t>& vec,
                                HighFive::File& file, const std::string& dset_name) {
  auto dset = file.getDataSet(dset_name);
  auto dims = dset.getDimensions();
  
  if (dims.empty() || dims[0] == 0) {
    vec.clear();
    return;
  }
  
  vec.resize(dims[0]);
  dset.read(vec);
}

// Forward declare - implementations need the full struct definitions
// These will be defined after including collocation_common.hpp and weights_generate.hpp

// Write XCTask vector to HDF5
inline void write_xctask_vector(const std::vector<XCTask>& tasks,
                                  const std::string& filename,
                                  const std::string& group_name = "/tasks") {
  HighFive::File file(filename, HighFive::File::OpenOrCreate);
  
  // Write number of tasks
  {
    HighFive::DataSpace space(1);
    auto dset = file.createDataSet<size_t>(group_name + "/ntasks", space);
    size_t ntasks = tasks.size();
    dset.write(&ntasks);
  }
  
  // Write each task
  for (size_t i = 0; i < tasks.size(); ++i) {
    std::string task_group = group_name + "/task_" + std::to_string(i);
    
    const auto& task = tasks[i];
    
    // Write scalar fields
    {
      HighFive::DataSpace space(1);
      auto dset_iparent = file.createDataSet<int32_t>(task_group + "/iParent", space);
      dset_iparent.write(&task.iParent);
      
      auto dset_npts = file.createDataSet<int32_t>(task_group + "/npts", space);
      dset_npts.write(&task.npts);
      
      auto dset_dist = file.createDataSet<double>(task_group + "/dist_nearest", space);
      dset_dist.write(&task.dist_nearest);
      
      auto dset_maxw = file.createDataSet<double>(task_group + "/max_weight", space);
      dset_maxw.write(&task.max_weight);
      
      auto dset_nbe = file.createDataSet<int32_t>(task_group + "/bfn_screening_nbe", space);
      dset_nbe.write(&task.bfn_screening.nbe);
    }
    
    // Write points and weights
    write_points_array(task.points, file, task_group + "/points");
    write_double_vector(task.weights, file, task_group + "/weights");
    
    // Write screening data
    write_int32_vector(task.bfn_screening.shell_list, file, task_group + "/shell_list");
  }
}

// Read XCTask vector from HDF5
inline void read_xctask_vector(std::vector<XCTask>& tasks,
                                 const std::string& filename,
                                 const std::string& group_name = "/tasks") {
  HighFive::File file(filename, HighFive::File::ReadOnly);
  
  // Read number of tasks
  size_t ntasks;
  {
    auto dset = file.getDataSet(group_name + "/ntasks");
    dset.read(&ntasks);
  }
  
  tasks.resize(ntasks);
  
  // Read each task
  for (size_t i = 0; i < ntasks; ++i) {
    std::string task_group = group_name + "/task_" + std::to_string(i);
    
    auto& task = tasks[i];
    
    // Read scalar fields
    {
      auto dset_iparent = file.getDataSet(task_group + "/iParent");
      dset_iparent.read(&task.iParent);
      
      auto dset_npts = file.getDataSet(task_group + "/npts");
      dset_npts.read(&task.npts);
      
      auto dset_dist = file.getDataSet(task_group + "/dist_nearest");
      dset_dist.read(&task.dist_nearest);
      
      auto dset_maxw = file.getDataSet(task_group + "/max_weight");
      dset_maxw.read(&task.max_weight);
      
      auto dset_nbe = file.getDataSet(task_group + "/bfn_screening_nbe");
      dset_nbe.read(&task.bfn_screening.nbe);
    }
    
    // Read points and weights
    read_points_array(task.points, file, task_group + "/points");
    read_double_vector(task.weights, file, task_group + "/weights");
    
    // Read screening data
    read_int32_vector(task.bfn_screening.shell_list, file, task_group + "/shell_list");
  }
}

} // namespace GauXC
