// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <mpi.h>
#include <omp.h>
#include <chrono>
#include <iostream>
#include <string>
#include "Eigen/Dense"
// #include "nlohmann/json.hpp"

#include "examples/implicit_surface_reconstruction/basic_mesh.h"
#include "examples/implicit_surface_reconstruction/binary.h"
#include "examples/implicit_surface_reconstruction/initialization.h"
#include "examples/implicit_surface_reconstruction/reconstruction.h"
#include "examples/implicit_surface_reconstruction/solver.h"
#include "examples/implicit_surface_reconstruction/surface_select.h"

// using json = nlohmann::json;

#define USE_MPI

static void print_usage_and_exit(int rank) {
  if (rank == 0) {
    std::cerr
        << "Usage:\n"
           "  mpirun -n P app --mode=init        "
           "--shape=<sphere|ellipsoid|...> --nx=<N>  --outdir=<dir>\n"
           "  mpirun -n P app --mode=metrics     --shape=<...>                 "
           "  --nx=<N>  --factor=<F> --outdir=<dir> --method=<Jibben|...>\n"
           "  mpirun -n P app --mode=convergence --shape=<...>                 "
           "  --nx=<N>               --outdir=<dir> --method=<Jibben|...> ";
  }
  MPI_Abort(MPI_COMM_WORLD, 1);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::string mode, shape, method, outdir;
  int Nx, factor;

  // parsing inputs
  for (int i = 1; i < argc; ++i) {
    const char* a = argv[i];
    if (std::strncmp(a, "--mode=", 7) == 0)
      mode = a + 7;
    else if (std::strncmp(a, "--shape=", 8) == 0)
      shape = a + 8;
    else if (std::strncmp(a, "--nx=", 5) == 0)
      Nx = std::stoi(a + 5);
    else if (std::strncmp(a, "--factor=", 9) == 0)
      factor = std::stoi(a + 9);
    else if (std::strncmp(a, "--outdir=", 9) == 0)
      outdir = a + 9;
    else if (std::strncmp(a, "--method=", 9) == 0)
      method = a + 9;
    else {
      if (rank == 0) std::cerr << "Unknown arg: " << a << "\n";
      print_usage_and_exit(rank);
    }
  }

  constexpr std::size_t VM_ORDER = 2;
  constexpr std::size_t SM_ORDER = 2;

  auto t0 = std::chrono::steady_clock::now();

  if (mode == "init") {
    // generating binary file with implicit surface moments
    run_initialization<VM_ORDER, SM_ORDER>(shape, Nx, outdir);
  } else if (mode == "metrics") {
    // only for paraboloids
    // reading binary file and computing metrics for one coarse mesh
    if (Nx % factor != 0) {
      if (rank == 0)
        std::fprintf(stderr,
                     "ERROR: --nx (%d) must be divisible by --factor (%d)\n",
                     Nx, factor);
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
    const MomentDiffNorms norms =
        computeReconstructionMetricsFromBin<VM_ORDER, SM_ORDER>(
            factor, method, shape, Nx, outdir);
  } else if (mode == "convergence") {
    // only for paraboloids
    run_convergence<VM_ORDER, SM_ORDER>(shape, Nx, method, outdir);
  } else if (mode == "curvednessConvergence") {
    runCurvednessConvergence<VM_ORDER, SM_ORDER>(shape, Nx, method, outdir);
  } else if (mode == "outputInterface") {
    output_interfaces<VM_ORDER, SM_ORDER>(shape, Nx, factor, method, outdir);
  } else {
    if (rank == 0)
      std::fprintf(stderr, "ERROR: Unknown --mode=%s\n", mode.c_str());
    print_usage_and_exit(rank);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  auto t1 = std::chrono::steady_clock::now();
  if (rank == 0) {
    const double runtime = std::chrono::duration<double>(t1 - t0).count();
    std::printf("Total wall time: %20.6f s\n", runtime);
  }

  MPI_Finalize();
  return 0;
}

// // ------- FOR TESTING PURPOSES -----------
// int main(int argc, char** argv) {
//   MPI_Init(&argc, &argv);
//   int rank = 0, size = 1;
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   MPI_Comm_size(MPI_COMM_WORLD, &size);
//   // std::string shape = "sphere";
//   // BasicMesh mesh(32, 32, 32, 1);
//   // SurfaceVariant surf = makeSurface(shape, mesh);

//   // double area = std::visit([](auto& s) { return s.surfaceArea(); }, surf);
//   // std::cout << area << std::endl;

//   // testing curvedness
//   constexpr std::size_t VM_ORDER = 2, SM_ORDER = 2;
//   const std::string shape = "ellipsoid";
//   const int Nx_fine = 256;
//   const int factor = 16;
//   const std::string reconstruction_method = "Jibben";
//   const std::string output_dir =
//       "/home/parinht2/Desktop/ppic paper/reconstruction_convergence/sweep_4";
//   // std::pair<double, double> curvedness_norms =
//   //     getCurvednessMetrics<VM_ORDER, SM_ORDER>(
//   //         shape, Nx_fine, factor, reconstruction_method, output_dir);

//   // std::cout << "Linf norm = " << curvedness_norms.first << std::endl;
//   // std::cout << "L2 norm = " << curvedness_norms.second << std::endl;

//   // outputting corresponding interface
//   output_interfaces<VM_ORDER, SM_ORDER>(shape, Nx_fine, factor,
//                                         reconstruction_method, output_dir);

//   // running curvedness convergence
//   // runCurvednessConvergence<VM_ORDER, SM_ORDER>(
//   //     shape, Nx_fine, reconstruction_method, output_dir);

//   MPI_Finalize();

//   return 0;
// }