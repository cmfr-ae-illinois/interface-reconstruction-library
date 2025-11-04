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
    if (rank == 0) {
      std::printf("[metrics] shape=%s Nx=%d factor=%d\n", shape.c_str(), Nx,
                  factor);
      std::printf("VOL  : M0  Linf=%.8e  L2=%.8e\n", norms.vol_M0_Linf,
                  norms.vol_M0_L2);
      std::printf("VOL  : M1  Linf=%.8e  L2=%.8e\n", norms.vol_M1_Linf,
                  norms.vol_M1_L2);
      std::printf("VOL  : M2  Linf=%.8e  L2=%.8e\n", norms.vol_M2_Linf,
                  norms.vol_M2_L2);
      std::printf("SURF : M0  Linf=%.8e  L2=%.8e\n", norms.surf_M0_Linf,
                  norms.surf_M0_L2);
      std::printf("SURF : M1  Linf=%.8e  L2=%.8e\n", norms.surf_M1_Linf,
                  norms.surf_M1_L2);
      std::printf("SURF : M2  Linf=%.8e  L2=%.8e\n", norms.surf_M2_Linf,
                  norms.surf_M2_L2);
    }
  } else if (mode == "convergence") {
    run_convergence<VM_ORDER, SM_ORDER>(shape, Nx, method, outdir);
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