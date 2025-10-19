// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_INITIALIZATION_TPP_
#define EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_INITIALIZATION_TPP_

#include "examples/implicit_surface_reconstruction/initialization.h"

struct Range {
  int begin;
  int end;
};

inline Range block_partition(int N, int rank, int size) {
  const int q = N / size;
  const int r = N % size;
  const int start = rank * q + std::min(rank, r);
  const int count = q + (rank < r ? 1 : 0);
  return {start, start + count};
}

template <std::size_t MAX_REFINE, std::size_t VM_ORDER, std::size_t SM_ORDER>
void getInitializedField(
    const Data<int>& cell_status,
    std::vector<std::tuple<int, int, int>> mixed_cells_list_root,
    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>* moments) {
  int rank = 0, size = 1;
  (void)MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  (void)MPI_Comm_size(MPI_COMM_WORLD, &size);

  // broadcasting mixed cells list to all ranks
  std::vector<int> flat;
  if (rank == 0) {
    flat.reserve(3 * mixed_cells_list_root.size());
    for (auto& t : mixed_cells_list_root) {
      flat.push_back(std::get<0>(t));
      flat.push_back(std::get<1>(t));
      flat.push_back(std::get<2>(t));
    }
  }
  int Ntriples =
      (rank == 0) ? static_cast<int>(mixed_cells_list_root.size()) : 0;
  MPI_Bcast(&Ntriples, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank != 0) flat.resize(3 * Ntriples);
  if (Ntriples > 0) {
    MPI_Bcast(flat.data(), 3 * Ntriples, MPI_INT, 0, MPI_COMM_WORLD);
  }

  // rebuilding mixed cells list on all ranks
  std::vector<std::tuple<int, int, int>> mixed_cells_list;
  mixed_cells_list.reserve(Ntriples);
  for (int n = 0; n < Ntriples; ++n) {
    mixed_cells_list.emplace_back(flat[3 * n + 0], flat[3 * n + 1],
                                  flat[3 * n + 2]);
  }

  const BasicMesh& mesh = cell_status.getMesh();

  IRL::Sphere<double, MAX_REFINE> sphere(0., 0., 0., 0.15);

  // partitioning the workload
  const int N = static_cast<int>(mixed_cells_list.size());
  Range rnge = block_partition(N, rank, size);

  // moment calculation
  double local_area = 0.0;
  for (int n = rnge.begin; n < rnge.end; ++n) {
    const auto [i, j, k] = mixed_cells_list[n];

    IRL::Pt x0(mesh.x(i), mesh.y(j), mesh.z(k));
    IRL::Pt x1(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
    IRL::RectangularCuboid cell =
        IRL::RectangularCuboid::fromBoundingPts(x0, x1);

    IRL::ImplicitSurfaceCutter<IRL::Sphere<double, MAX_REFINE>,
                               IRL::GeneralMoments3D<VM_ORDER>>
        cutter(sphere, cell);

    auto vol = cutter.computeVolumeMoments();
    auto surf = cutter.template computeSurfaceMoments<SM_ORDER>(
        false, Eigen::Integrator<double, 2>::GaussKronrod15, 5);

    (*moments)(i, j, k).first = vol;
    (*moments)(i, j, k).second = surf;

    local_area += surf[0];
  }

  MPI_Barrier(MPI_COMM_WORLD);

  double global_area = 0.0;

  MPI_Reduce(&local_area, &global_area, 1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);

  if (rank == 0) {
    std::cout << "Total reconstructed surface area = " << global_area
              << std::endl;
  }
}

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_INITIALIZATION_TPP_
