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

template <std::size_t MAX_REFINE, std::size_t VM_ORDER, std::size_t SM_ORDER>
void getInitializedField(
    const Data<int>& cell_status,
    std::vector<std::tuple<int, int, int>> mixed_cells_list_root,
    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>* moments) {
  int rank = 0, size = 1;
  (void)MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  (void)MPI_Comm_size(MPI_COMM_WORLD, &size);

  // broadcast mixed (i,j,k) to all ranks
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

  std::vector<std::tuple<int, int, int>> mixed_cells_list;
  mixed_cells_list.reserve(Ntriples);
  for (int n = 0; n < Ntriples; ++n) {
    mixed_cells_list.emplace_back(flat[3 * n + 0], flat[3 * n + 1],
                                  flat[3 * n + 2]);
  }

  const BasicMesh& mesh = cell_status.getMesh();
  IRL::Sphere<double, MAX_REFINE> sphere(0., 0., 0., 0.15);

  // partition mixed cells across ranks
  const int N = static_cast<int>(mixed_cells_list.size());
  Range rnge = block_partition(N, rank, size);

  // computing local moment results
  constexpr int NV = (VM_ORDER + 1) * (VM_ORDER + 2) * (VM_ORDER + 3) / 6;
  constexpr int NS = (SM_ORDER + 1) * (SM_ORDER + 2) * (SM_ORDER + 3) / 6;

  struct CellMomentsRec {
    int i, j, k;
    double vol[NV];
    double surf[NS];
  };

  std::vector<CellMomentsRec> local;
  local.reserve(std::max(0, rnge.end - rnge.begin));

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

    CellMomentsRec rec{};
    rec.i = i;
    rec.j = j;
    rec.k = k;
    for (int t = 0; t < NV; ++t) rec.vol[t] = vol[t];
    for (int t = 0; t < NS; ++t) rec.surf[t] = surf[t];
    local.push_back(rec);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  int local_n = static_cast<int>(local.size());

  std::vector<int> counts(size, 0);
  MPI_Gather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT,
             /*root=*/0, MPI_COMM_WORLD);

  std::vector<int> displs, counts_bytes, displs_bytes;
  int total_n = 0;
  const int rec_bytes = static_cast<int>(sizeof(CellMomentsRec));

  if (rank == 0) {
    displs.resize(size);
    counts_bytes.resize(size);
    displs_bytes.resize(size);
    displs[0] = 0;
    for (int r = 1; r < size; ++r) displs[r] = displs[r - 1] + counts[r - 1];
    total_n = displs[size - 1] + counts[size - 1];
    for (int r = 0; r < size; ++r) {
      counts_bytes[r] = counts[r] * rec_bytes;
      displs_bytes[r] = displs[r] * rec_bytes;
    }
  }

  std::vector<CellMomentsRec> all;
  if (rank == 0) all.resize(total_n);

  MPI_Gatherv(reinterpret_cast<const void*>(local.data()), local_n * rec_bytes,
              MPI_BYTE,
              rank == 0 ? reinterpret_cast<void*>(all.data()) : nullptr,
              rank == 0 ? counts_bytes.data() : nullptr,
              rank == 0 ? displs_bytes.data() : nullptr, MPI_BYTE, /*root=*/0,
              MPI_COMM_WORLD);

  // writing all moments to rank 0
  if (rank == 0) {
    for (const auto& rec : all) {
      IRL::GeneralMoments3D<VM_ORDER> vol{};
      IRL::GeneralSurfaceMoments3D<SM_ORDER> surf{};
      for (int t = 0; t < NV; ++t) vol[t] = rec.vol[t];
      for (int t = 0; t < NS; ++t) surf[t] = rec.surf[t];
      (*moments)(rec.i, rec.j, rec.k).first = vol;
      (*moments)(rec.i, rec.j, rec.k).second = surf;
    }
    // volume moments for cells below
    for (int i = mesh.imin(); i <= mesh.imax(); i++) {
      for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
        for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
          if (cell_status(i, j, k) == -1) {
            IRL::Pt x0(mesh.x(i), mesh.y(j), mesh.z(k));
            IRL::Pt x1(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
            IRL::RectangularCuboid cell =
                IRL::RectangularCuboid::fromBoundingPts(x0, x1);
            IRL::GeneralMoments3D<VM_ORDER> cell_moment =
                IRL::getVolumeMoments<IRL::GeneralMoments3D<VM_ORDER>>(cell);
            (*moments)(i, j, k).first = cell_moment;
          }
        }
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_INITIALIZATION_TPP_
