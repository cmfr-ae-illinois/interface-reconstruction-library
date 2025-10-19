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
#include "nlohmann/json.hpp"

#include "examples/implicit_surface_reconstruction/basic_mesh.h"
#include "examples/implicit_surface_reconstruction/initialization.h"
#include "examples/implicit_surface_reconstruction/solver.h"

using json = nlohmann::json;

#define USE_MPI

// int main(int argc, char* argv[]) {
//   const int Nx = 32;
//   const std::string reconstruction_method = "Jibben";

//   // setting mesh
//   BasicMesh mesh(Nx, Nx, Nx, 1);
//   mesh.setCellBoundaries(IRL::Pt(-0.18, -0.18, -0.18),
//                          IRL::Pt(0.18, 0.18, 0.18));

// #ifdef USE_MPI
//   MPI_Init(&argc, &argv);
//   int rank = 0;
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// #endif

//   // cell status and mixed cells list
//   Data<int> cell_status(&mesh);
//   std::vector<std::tuple<int, int, int>> mixed_cells_list =
//       getCellStatus(&cell_status);

//   // moments
//   constexpr std::size_t SM_ORDER = 2, VM_ORDER = 2;
//   Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
//                  IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
//       moments(&mesh);

//   auto t0 = std::chrono::steady_clock::now();
//   getInitializedField(cell_status, mixed_cells_list, &moments);
//   auto t1 = std::chrono::steady_clock::now();
//   const double runtime = std::chrono::duration<double>(t1 - t0).count();

// #ifdef USE_MPI
//   if (rank == 0) {
//     std::printf("Total run time: %20.6f s\n", runtime);
//   }
//   MPI_Finalize();
// #else
//   std::printf("Total run time: %20.6f s\n", runtime);
// #endif

//   return 0;
// }

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  const int Nx = 32;
  const std::string reconstruction_method = "Jibben";

  // setting mesh
  BasicMesh mesh(Nx, Nx, Nx, 1);
  mesh.setCellBoundaries(IRL::Pt(-0.18, -0.18, -0.18),
                         IRL::Pt(0.18, 0.18, 0.18));

  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  Data<int> cell_status(&mesh);
  std::vector<std::tuple<int, int, int>> mixed_cells_list;
  if (rank == 0) {
    mixed_cells_list = getCellStatus(&cell_status);
  }

  // AMR params
  constexpr std::size_t MAX_REFINE = 4;
  constexpr std::size_t VM_ORDER = 2;
  constexpr std::size_t SM_ORDER = 2;

  Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                 IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
      moments(&mesh);

  // Parallel compute. Choose ONE of the two output patterns below:
  // A) per-rank binary outputs (simple, scalable):
  auto t0 = std::chrono::steady_clock::now();
  getInitializedField<MAX_REFINE, VM_ORDER, SM_ORDER>(
      cell_status, mixed_cells_list, &moments);
  auto t1 = std::chrono::steady_clock::now();
  const double runtime = std::chrono::duration<double>(t1 - t0).count();

  if (rank == 0) {
    std::printf("Total run time: %20.6f s\n", runtime);
  }

  MPI_Finalize();
  return 0;
}
