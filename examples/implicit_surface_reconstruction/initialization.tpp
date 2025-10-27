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

// wrapper to set refine level to be 0 for a surface
template <class S>
struct ZeroRefineSurface : S {
  using S::S;
  ZeroRefineSurface(const S& s) : S(s) {}
  static constexpr std::size_t getMaxRefineLevel() { return 0; }
};

// finding mixed cells ----------------------------------------------
template <class SurfaceType>
std::vector<std::tuple<int, int, int>> getCellStatus(
    Data<int>* cell_status, const SurfaceType& surface) {
  using S0 = ZeroRefineSurface<std::decay_t<SurfaceType>>;
  S0 surface0(surface);

  std::vector<std::tuple<int, int, int>> mixed_cells_list;
  const BasicMesh& mesh = cell_status->getMesh();

  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        IRL::Pt x0(mesh.x(i), mesh.y(j), mesh.z(k));
        IRL::Pt x1(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
        IRL::RectangularCuboid cell =
            IRL::RectangularCuboid::fromBoundingPts(x0, x1);

        IRL::ImplicitSurfaceCutter<S0, IRL::Volume> cutter(surface0, cell);

        const int status = cutter.getBaseCellStatus();
        (*cell_status)(i, j, k) = status;
        if (status == 0) mixed_cells_list.emplace_back(i, j, k);
      }
    }
  }

  return mixed_cells_list;
}

// finding implicit surface moments -----------------------------------
template <class SurfaceType, std::size_t VM_ORDER, std::size_t SM_ORDER>
void getInitializedField(
    const Data<int>& cell_status,
    std::vector<std::tuple<int, int, int>> mixed_cells_list_root,
    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>* moments,
    const SurfaceType& surface) {
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

    IRL::ImplicitSurfaceCutter<SurfaceType, IRL::GeneralMoments3D<VM_ORDER>>
        cutter(surface, cell);

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
  MPI_Gather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0,
             MPI_COMM_WORLD);

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

// initializing moments and writing to binary ------------------------------
template <class SurfaceType, std::size_t VM_ORDER, std::size_t SM_ORDER>
Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
               IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
initializeMomentsAndWriteBin(const BasicMesh& mesh, const SurfaceType& surface,
                             const std::string& bin_path) {
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  Data<int> cell_status(&mesh);
  std::vector<std::tuple<int, int, int>> mixed_cells_list;
  if (rank == 0) {
    mixed_cells_list = getCellStatus(&cell_status, surface);
  }

  Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                 IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
      moments(&mesh);
  getInitializedField<SurfaceType, VM_ORDER, SM_ORDER>(
      cell_status, mixed_cells_list, &moments, surface);

  // writing moments to binary
  if (rank == 0) {
    writeMomentsToBinary<VM_ORDER, SM_ORDER>(moments, bin_path);
  }

  return moments;
}

// running initialization for generic shape ---------------------------
template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void run_initialization(const std::string& shape, int Nx,
                        const std::string& output_dir) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  BasicMesh mesh(Nx, Nx, Nx, 1);
  SurfaceVariant surf = makeSurface(shape, mesh);
  std::string bin_path = binary_filename(output_dir, shape, Nx);

  std::visit(
      [&](const auto& surface) {
        using S = std::decay_t<decltype(surface)>;

        auto t0 = std::chrono::steady_clock::now();
        auto data = initializeMomentsAndWriteBin<S, VM_ORDER, SM_ORDER>(
            mesh, surface, bin_path);
        auto t1 = std::chrono::steady_clock::now();
        if (rank == 0) {
          const double dt = std::chrono::duration<double>(t1 - t0).count();
          std::cout << "✅ Initialization complete for " << shape << " in "
                    << dt << " s\n";
        }
      },
      surf);
}

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_INITIALIZATION_TPP_
