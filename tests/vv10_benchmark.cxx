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
#include <gauxc/external/hdf5.hpp>
#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/molecular_weights.hpp>
#include <gauxc/runtime_environment.hpp>
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>

#include <highfive/H5File.hpp>

#define EIGEN_DONT_VECTORIZE
#define EIGEN_NO_CUDA
#include <Eigen/Core>

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using matrix_type = Eigen::MatrixXd;

struct Options {
  std::vector<std::string> files;
  int repeats = 5;
  int warmup = 1;
  std::size_t batch_size = 512;
  double basis_tol = 1e-10;
  std::string grid = "GM3";
  std::string pruning = "UNPRUNED";
  std::string rad_quad = "MURAKNOWLES";
  std::string exec_space = "Host";
  std::string nlc_math_mode = "FP64";
  std::string quantity = "VNLC";
  std::string gradient_mode = "FULL";
  bool force_rks = false;
  std::string integrator_kernel = "Default";
  std::string lwd_kernel = "Default";
  std::string reduction_kernel = "Default";
};

void uppercase(std::string& value) {
  std::transform(value.begin(), value.end(), value.begin(), ::toupper);
}

void print_usage(const char* exe) {
  std::cout << "Usage: " << exe << " [options] file1.h5 [file2.h5 ...]\n"
            << "Options:\n"
            << "  --repeats N              Timed repetitions (default: 5)\n"
            << "  --warmup N               Untimed warmup repetitions (default: 1)\n"
            << "  --batch-size N           Grid batch size (default: 512)\n"
            << "  --basis-tol X            Shell screening tolerance (default: 1e-10)\n"
            << "  --grid NAME              FINE, ULTRAFINE, SUPERFINE, GM3, GM5 (default: GM3)\n"
            << "  --pruning NAME           UNPRUNED, ROBUST, TREUTLER (default: UNPRUNED)\n"
            << "  --rad-quad NAME          MK/MURAKNOWLES, BECKE, TA, MHL (default: MURAKNOWLES)\n"
            << "  --exec-space NAME        Host or Device (default: Host)\n"
            << "  --nlc-math-mode NAME     FP64 or FLOATPAIR for exact device VV10 (default: FP64)\n"
            << "  --quantity NAME          VNLC, GRAD, or FNLC (default: VNLC)\n"
            << "  --gradient-mode NAME     FULL or HF for --quantity GRAD (default: FULL)\n"
            << "  --force-rks              Use scalar density as RKS even if spin datasets are present\n"
            << "  --integrator-kernel NAME Integrator kernel (default: Default)\n"
            << "  --lwd-kernel NAME        Local work driver kernel (default: Default)\n"
            << "  --reduction-kernel NAME  Reduction driver kernel (default: Default)\n";
}

Options parse_options(int argc, char** argv) {
  Options options;
  for(int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto require_value = [&](const char* name) -> std::string {
      if(i + 1 >= argc) throw std::runtime_error(std::string("Missing value for ") + name);
      return argv[++i];
    };

    if(arg == "--help" or arg == "-h") {
      print_usage(argv[0]);
      std::exit(0);
    } else if(arg == "--repeats") {
      options.repeats = std::stoi(require_value("--repeats"));
    } else if(arg == "--warmup") {
      options.warmup = std::stoi(require_value("--warmup"));
    } else if(arg == "--batch-size") {
      options.batch_size = static_cast<std::size_t>(std::stoull(require_value("--batch-size")));
    } else if(arg == "--basis-tol") {
      options.basis_tol = std::stod(require_value("--basis-tol"));
    } else if(arg == "--grid") {
      options.grid = require_value("--grid");
    } else if(arg == "--pruning") {
      options.pruning = require_value("--pruning");
    } else if(arg == "--rad-quad") {
      options.rad_quad = require_value("--rad-quad");
    } else if(arg == "--exec-space") {
      options.exec_space = require_value("--exec-space");
    } else if(arg == "--nlc-math-mode") {
      options.nlc_math_mode = require_value("--nlc-math-mode");
    } else if(arg == "--quantity") {
      options.quantity = require_value("--quantity");
    } else if(arg == "--gradient-mode") {
      options.gradient_mode = require_value("--gradient-mode");
    } else if(arg == "--force-rks") {
      options.force_rks = true;
    } else if(arg == "--integrator-kernel") {
      options.integrator_kernel = require_value("--integrator-kernel");
    } else if(arg == "--lwd-kernel") {
      options.lwd_kernel = require_value("--lwd-kernel");
    } else if(arg == "--reduction-kernel") {
      options.reduction_kernel = require_value("--reduction-kernel");
    } else {
      options.files.push_back(arg);
    }
  }

  uppercase(options.grid);
  uppercase(options.pruning);
  uppercase(options.rad_quad);
  uppercase(options.nlc_math_mode);
  uppercase(options.quantity);
  uppercase(options.gradient_mode);
  if(options.files.empty()) throw std::runtime_error("At least one HDF5 benchmark file is required");
  if(options.repeats <= 0) throw std::runtime_error("--repeats must be positive");
  if(options.warmup < 0) throw std::runtime_error("--warmup must be non-negative");
  if(options.quantity != "VNLC" and options.quantity != "GRAD" and options.quantity != "FNLC") {
    throw std::runtime_error("Unknown --quantity: " + options.quantity);
  }
  if(options.nlc_math_mode != "FP64" and options.nlc_math_mode != "FLOATPAIR") {
    throw std::runtime_error("Unknown --nlc-math-mode: " + options.nlc_math_mode);
  }
  if(options.gradient_mode != "FULL" and options.gradient_mode != "HF") {
    throw std::runtime_error("Unknown --gradient-mode: " + options.gradient_mode);
  }
  return options;
}

matrix_type read_matrix(HighFive::File& file, const std::string& name) {
  auto dset = file.getDataSet(name);
  auto dims = dset.getDimensions();
  if(dims.size() != 2 or dims[0] != dims[1]) {
    throw std::runtime_error("Expected square matrix dataset: " + name);
  }
  matrix_type matrix(dims[0], dims[1]);
  dset.read(matrix.data());
  return matrix;
}

double mpi_max_seconds(double local_seconds) {
#ifdef GAUXC_HAS_MPI
  double global_seconds = 0.0;
  MPI_Reduce(&local_seconds, &global_seconds, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  return global_seconds;
#else
  return local_seconds;
#endif
}

template <typename F>
double time_call(F&& func) {
#ifdef GAUXC_HAS_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  const auto start = std::chrono::high_resolution_clock::now();
  func();
#ifdef GAUXC_HAS_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  const auto stop = std::chrono::high_resolution_clock::now();
  return std::chrono::duration<double>(stop - start).count();
}

void print_result(int rank, const std::string& file, const std::string& mode,
                  const std::string& quantity,
                  int natoms, int nbf, const std::vector<double>& seconds,
                  double value, double norm) {
  if(rank != 0) return;
  auto sorted = seconds;
  std::sort(sorted.begin(), sorted.end());
  const auto sum = std::accumulate(sorted.begin(), sorted.end(), 0.0);
  const auto min = sorted.front();
  const auto median = sorted[sorted.size() / 2];
  const auto avg = sum / static_cast<double>(sorted.size());
  std::cout << std::left << std::setw(28) << file
            << std::right << std::setw(8) << mode
            << std::setw(10) << quantity
            << std::setw(8) << natoms
            << std::setw(8) << nbf
            << std::setw(22) << min
            << std::setw(22) << median
            << std::setw(22) << avg
            << std::setw(24) << value
            << std::setw(24) << norm
            << '\n';
}

template <typename RuntimeType>
void run_benchmarks(const Options& options, RuntimeType& runtime) {
  using namespace GauXC;
  using namespace ExchCXX;

  const int rank = runtime.comm_rank();

  std::map<std::string, AtomicGridSizeDefault> grid_map = {
    {"FINE", AtomicGridSizeDefault::FineGrid},
    {"ULTRAFINE", AtomicGridSizeDefault::UltraFineGrid},
    {"SUPERFINE", AtomicGridSizeDefault::SuperFineGrid},
    {"GM3", AtomicGridSizeDefault::GM3},
    {"GM5", AtomicGridSizeDefault::GM5}
  };
  std::map<std::string, PruningScheme> pruning_map = {
    {"UNPRUNED", PruningScheme::Unpruned},
    {"ROBUST", PruningScheme::Robust},
    {"TREUTLER", PruningScheme::Treutler}
  };
  std::map<std::string, RadialQuad> rad_quad_map = {
    {"BECKE", RadialQuad::Becke},
    {"MURAKNOWLES", RadialQuad::MuraKnowles},
    {"MK", RadialQuad::MuraKnowles},
    {"TREUTLERAHLRICHS", RadialQuad::TreutlerAhlrichs},
    {"TA", RadialQuad::TreutlerAhlrichs},
    {"MURRAYHANDYLAMING", RadialQuad::MurrayHandyLaming},
    {"MHL", RadialQuad::MurrayHandyLaming}
  };

  if(rank == 0) {
    std::cout << std::scientific << std::setprecision(12);
    std::cout << std::left << std::setw(28) << "file"
              << std::right << std::setw(8) << "mode"
              << std::setw(10) << "quantity"
              << std::setw(8) << "natoms"
              << std::setw(8) << "nbf"
              << std::setw(22) << "min_s"
              << std::setw(22) << "median_s"
              << std::setw(22) << "avg_s"
              << std::setw(24) << "value"
              << std::setw(24) << "out_norm"
              << '\n';
  }

  for(const auto& file_name : options.files) {
    Molecule mol;
    BasisSet<double> basis;
    read_hdf5_record(mol, file_name, "/MOLECULE");
    read_hdf5_record(basis, file_name, "/BASIS");
    for(auto& shell : basis) shell.set_shell_tolerance(options.basis_tol);

    HighFive::File file(file_name, HighFive::File::ReadOnly);
    const bool has_spin_density = file.exist("/DENSITY_Z");
    const bool uks = has_spin_density and not options.force_rks;
    const auto density_name = has_spin_density ? std::string("/DENSITY_SCALAR") : std::string("/DENSITY");
    auto P = read_matrix(file, density_name);
    matrix_type Pz;
    if(uks) Pz = read_matrix(file, "/DENSITY_Z");

    const auto molgrid = MolGridFactory::create_default_molgrid(
      mol, pruning_map.at(options.pruning), BatchSize(options.batch_size),
      rad_quad_map.at(options.rad_quad), grid_map.at(options.grid));

    LoadBalancerFactory lb_factory(ExecutionSpace::Host, "Default");
    auto load_balancer = lb_factory.get_shared_instance(runtime, mol, molgrid, basis);

    const auto exec_space = options.exec_space == "Device" ? ExecutionSpace::Device : ExecutionSpace::Host;
    MolecularWeightsFactory mw_factory(exec_space, "Default", MolecularWeightsSettings{});
    auto molecular_weights = mw_factory.get_instance();
    molecular_weights.modify_weights(*load_balancer);

    auto functional = functional_type(Backend::builtin, Functional::BLYP,
      uks ? Spin::Polarized : Spin::Unpolarized);
    XCIntegratorFactory<matrix_type> integrator_factory(exec_space, "Replicated",
      options.integrator_kernel, options.lwd_kernel, options.reduction_kernel);
    auto integrator = integrator_factory.get_instance(functional, load_balancer);

    IntegratorSettingsNLC settings;
    settings.math_mode = options.nlc_math_mode == "FLOATPAIR"
      ? NLCMathMode::FloatPair : NLCMathMode::NativeFP64;
    settings.include_weight_derivatives = options.gradient_mode == "FULL";
    double value = 0.0;
    double norm = 0.0;

    auto evaluate = [&]() {
      if(options.quantity == "VNLC") {
        if(uks) {
          auto result = integrator.eval_nlc_vnlc(P, Pz, settings);
          value = std::get<0>(result);
          norm = std::get<1>(result).norm();
        } else {
          auto result = integrator.eval_nlc_vnlc(P, settings);
          value = std::get<0>(result);
          norm = std::get<1>(result).norm();
        }
      } else if(options.quantity == "GRAD") {
        std::vector<double> gradient;
        if(uks) gradient = integrator.eval_nlc_grad(P, Pz, settings);
        else gradient = integrator.eval_nlc_grad(P, settings);
        value = 0.0;
        norm = std::sqrt(std::inner_product(gradient.begin(), gradient.end(), gradient.begin(), 0.0));
      } else {
        if(uks) {
          auto result = integrator.eval_nlc_fnlc_contraction(P, Pz, P, Pz, settings);
          value = 0.0;
          norm = std::get<0>(result).norm() + std::get<1>(result).norm();
        } else {
          auto result = integrator.eval_nlc_fnlc_contraction(P, P, settings);
          value = 0.0;
          norm = result.norm();
        }
      }
    };

    for(int i = 0; i < options.warmup; ++i) {
      evaluate();
    }

    std::vector<double> seconds;
    seconds.reserve(options.repeats);
    for(int i = 0; i < options.repeats; ++i) {
      const auto local_seconds = time_call(evaluate);
      seconds.push_back(mpi_max_seconds(local_seconds));
    }

    print_result(rank, file_name, uks ? "UKS" : "RKS", options.quantity,
      mol.natoms(), basis.nbf(), seconds, value, norm);
  }
}

} // namespace

int main(int argc, char** argv) {
#ifdef GAUXC_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  int status = 0;
  try {
    const auto options = parse_options(argc, argv);
#ifdef GAUXC_HAS_DEVICE
    if(options.exec_space == "Device") {
      auto runtime = GauXC::DeviceRuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD,) 0.9);
      run_benchmarks(options, runtime);
    } else
#endif
    {
      auto runtime = GauXC::RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
      run_benchmarks(options, runtime);
    }
  } catch(const std::exception& ex) {
    std::cerr << "vv10_benchmark error: " << ex.what() << std::endl;
    status = 1;
  }

#ifdef GAUXC_HAS_MPI
  MPI_Finalize();
#endif
  return status;
}