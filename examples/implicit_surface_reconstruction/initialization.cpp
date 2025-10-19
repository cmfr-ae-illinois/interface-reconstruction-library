// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/implicit_surface_reconstruction/initialization.h"

// recentering M1
Eigen::Vector3d recenter_M1(const Eigen::Vector3d& M1, const double& M0,
                            const Eigen::Vector3d& x_ref) {
  return M1 - (M0 * x_ref);
}

// recentering M2
Eigen::Matrix3d recenter_M2(const Eigen::Matrix3d& M2,
                            const Eigen::Vector3d& M1, const double& M0,
                            const Eigen::Vector3d& x_ref) {
  return M2 - (x_ref * M1.transpose()) - (M1 * x_ref.transpose()) +
         (M0 * x_ref * x_ref.transpose());
}

// finding mixed cells
std::vector<std::tuple<int, int, int>> getCellStatus(Data<int>* cell_status) {
  std::vector<std::tuple<int, int, int>> mixed_cells_list;

  const BasicMesh& mesh = cell_status->getMesh();
  IRL::Sphere<double, 0> sphere(0., 0., 0., 0.15);

  for (int i = mesh.imin(); i <= mesh.imax(); i++) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); j++) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); k++) {
        IRL::Pt x0(mesh.x(i), mesh.y(j), mesh.z(k));
        IRL::Pt x1(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
        IRL::RectangularCuboid cell =
            IRL::RectangularCuboid::fromBoundingPts(x0, x1);
        IRL::ImplicitSurfaceCutter<IRL::Sphere<double, 0>, IRL::Volume> cutter(
            sphere, cell);
        (*cell_status)(i, j, k) = cutter.getBaseCellStatus();
        if ((*cell_status)(i, j, k) == 0)
          mixed_cells_list.emplace_back(i, j, k);
      }
    }
  }

  return mixed_cells_list;
}

// void getInitializedField(
//     const Data<int>& cell_status,
//     const std::vector<std::tuple<int, int, int>>& mixed_cells_list,
//     Data<std::pair<IRL::GeneralMoments3D<2>,
//     IRL::GeneralSurfaceMoments3D<2>>>*
//         moments) {
//   constexpr std::size_t MAX_REFINE = 4;
//   constexpr std::size_t SM_ORDER = 2;
//   constexpr std::size_t VM_ORDER = 2;

//   const BasicMesh& mesh = cell_status.getMesh();
//   IRL::Sphere<double, MAX_REFINE> sphere(0., 0., 0., 0.15);

//   int nmixed_global = mixed_cells_list.size();
//   double area = 0.;
//   for (int n = 0; n < mixed_cells_list.size(); n++) {
//     const auto [i, j, k] = mixed_cells_list[n];
//     IRL::Pt x0(mesh.x(i), mesh.y(j), mesh.z(k));
//     IRL::Pt x1(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
//     IRL::RectangularCuboid cell =
//         IRL::RectangularCuboid::fromBoundingPts(x0, x1);

//     IRL::ImplicitSurfaceCutter<IRL::Sphere<double, MAX_REFINE>,
//                                IRL::GeneralMoments3D<VM_ORDER>>
//         cutter(sphere, cell);

//     auto vol = cutter.computeVolumeMoments();
//     auto surf = cutter.template computeSurfaceMoments<SM_ORDER>(
//         false, Eigen::Integrator<double, 2>::GaussKronrod15, 5);

//     (*moments)(i, j, k).first = vol;
//     (*moments)(i, j, k).second = surf;
//     area += surf[0];
//   }
//   std::cout << area << std::endl;
//   std::cout << sphere.surfaceArea() << std::endl;
// }

// ---------------------------------------------------------------------------

// struct Range {
//   int begin;
//   int end;
// };  // [begin, end)
// inline Range block_partition(int N, int rank, int size) {
//   const int q = N / size;
//   const int r = N % size;
//   const int start = rank * q + std::min(rank, r);
//   const int count = q + (rank < r ? 1 : 0);
//   return {start, start + count};
// }

// template <std::size_t MAX_REFINE, std::size_t VM_ORDER, std::size_t SM_ORDER>
// void getInitializedField_MPI(
//     const Data<int>& cell_status_global_or_local,
//     std::vector<std::tuple<int, int, int>>
//         mixed_cells_list_root,  // valid only on rank 0
//     Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
//                    IRL::GeneralSurfaceMoments3D<SM_ORDER>>>* moments) {
//   int rank = 0, size = 1;
//   (void)MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   (void)MPI_Comm_size(MPI_COMM_WORLD, &size);

//   // ---- 1) Broadcast mixed_cells_list to all ranks ----
//   // Pack rank0 vector<int,int,int> as flat ints
//   std::vector<int> flat;
//   if (rank == 0) {
//     flat.reserve(3 * mixed_cells_list_root.size());
//     for (auto& t : mixed_cells_list_root) {
//       flat.push_back(std::get<0>(t));
//       flat.push_back(std::get<1>(t));
//       flat.push_back(std::get<2>(t));
//     }
//   }
//   int Ntriples =
//       (rank == 0) ? static_cast<int>(mixed_cells_list_root.size()) : 0;
//   MPI_Bcast(&Ntriples, 1, MPI_INT, 0, MPI_COMM_WORLD);
//   if (rank != 0) flat.resize(3 * Ntriples);
//   if (Ntriples > 0) {
//     MPI_Bcast(flat.data(), 3 * Ntriples, MPI_INT, 0, MPI_COMM_WORLD);
//   }

//   // Rebuild mixed_cells_list on all ranks
//   std::vector<std::tuple<int, int, int>> mixed_cells_list;
//   mixed_cells_list.reserve(Ntriples);
//   for (int n = 0; n < Ntriples; ++n) {
//     mixed_cells_list.emplace_back(flat[3 * n + 0], flat[3 * n + 1],
//                                   flat[3 * n + 2]);
//   }

//   // ---- 2) Everyone has mesh (we read from the local Data<> we were passed)
//   // ----
//   const BasicMesh& mesh = cell_status_global_or_local.getMesh();

//   // The surface definition (identical everywhere)
//   IRL::Sphere<double, MAX_REFINE> sphere(0., 0., 0., 0.15);

//   // ---- 3) Partition the workload ----
//   const int N = static_cast<int>(mixed_cells_list.size());
//   Range rnge = block_partition(N, rank, size);

//   // ---- 4) Compute locally on the assigned slice ----
//   // We'll hold local results in a simple struct for Option A or B.
//   struct Item {
//     int i, j, k;
//     IRL::GeneralMoments3D<VM_ORDER> vol;
//     IRL::GeneralSurfaceMoments3D<SM_ORDER> surf;
//   };
//   std::vector<Item> local;
//   local.reserve(std::max(0, rnge.end - rnge.begin));

//   for (int n = rnge.begin; n < rnge.end; ++n) {
//     const auto [i, j, k] = mixed_cells_list[n];

//     IRL::Pt x0(mesh.x(i), mesh.y(j), mesh.z(k));
//     IRL::Pt x1(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
//     IRL::RectangularCuboid cell =
//         IRL::RectangularCuboid::fromBoundingPts(x0, x1);

//     IRL::ImplicitSurfaceCutter<IRL::Sphere<double, MAX_REFINE>,
//                                IRL::GeneralMoments3D<VM_ORDER>>
//         cutter(sphere, cell);

//     auto vol = cutter.computeVolumeMoments();
//     auto surf = cutter.template computeSurfaceMoments<SM_ORDER>(
//         /*useAdaptive=*/false, Eigen::Integrator<double, 2>::GaussKronrod15,
//         /*npts=*/5);

//     // Fill local Data<> too (so every rank has its slice placed locally)
//     // NOTE: Data<> is process-local; there’s no shared memory across ranks.
//     (*moments)(i, j, k).first = vol;
//     (*moments)(i, j, k).second = surf;

//     local.push_back({i, j, k, std::move(vol), std::move(surf)});
//   }

//   MPI_Barrier(MPI_COMM_WORLD);
// }

// -----------------------------------------------------------------------------

// template <std::size_t MAX_REFINE, std::size_t VM_ORDER, std::size_t SM_ORDER>
// void getInitializedField(
//     const Data<int>& cell_status_global_or_local,
//     std::vector<std::tuple<int, int, int>> mixed_cells_list_root,
//     Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
//                    IRL::GeneralSurfaceMoments3D<SM_ORDER>>>* moments) {
// #ifdef USE_MPI
//   int rank = 0, size = 1;
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   MPI_Comm_size(MPI_COMM_WORLD, &size);
// #else
//   int rank = 0, size = 1;
// #endif

//   // ---- 1) Broadcast mixed_cells_list to all ranks ----
//   std::vector<int> flat;
// #ifdef USE_MPI
//   if (rank == 0) {
//     flat.reserve(3 * mixed_cells_list_root.size());
//     for (auto& t : mixed_cells_list_root) {
//       flat.push_back(std::get<0>(t));
//       flat.push_back(std::get<1>(t));
//       flat.push_back(std::get<2>(t));
//     }
//   }

//   int Ntriples =
//       (rank == 0) ? static_cast<int>(mixed_cells_list_root.size()) : 0;
//   MPI_Bcast(&Ntriples, 1, MPI_INT, 0, MPI_COMM_WORLD);

//   if (rank != 0) flat.resize(3 * Ntriples);
//   if (Ntriples > 0)
//     MPI_Bcast(flat.data(), 3 * Ntriples, MPI_INT, 0, MPI_COMM_WORLD);
// #else
//   // In serial mode, just use the root’s data directly
//   flat.reserve(3 * mixed_cells_list_root.size());
//   for (auto& t : mixed_cells_list_root) {
//     flat.push_back(std::get<0>(t));
//     flat.push_back(std::get<1>(t));
//     flat.push_back(std::get<2>(t));
//   }
//   int Ntriples = static_cast<int>(mixed_cells_list_root.size());
// #endif

//   // ---- 2) Rebuild mixed_cells_list ----
//   std::vector<std::tuple<int, int, int>> mixed_cells_list;
//   mixed_cells_list.reserve(Ntriples);
//   for (int n = 0; n < Ntriples; ++n) {
//     mixed_cells_list.emplace_back(flat[3 * n + 0], flat[3 * n + 1],
//                                   flat[3 * n + 2]);
//   }

//   // ---- 3) Mesh setup ----
//   const BasicMesh& mesh = cell_status_global_or_local.getMesh();
//   IRL::Sphere<double, MAX_REFINE> sphere(0., 0., 0., 0.15);

//   // ---- 4) Work partition ----
//   const int N = static_cast<int>(mixed_cells_list.size());
//   Range rnge = block_partition(N, rank, size);

//   struct Item {
//     int i, j, k;
//     IRL::GeneralMoments3D<VM_ORDER> vol;
//     IRL::GeneralSurfaceMoments3D<SM_ORDER> surf;
//   };
//   std::vector<Item> local;
//   local.reserve(std::max(0, rnge.end - rnge.begin));

//   for (int n = rnge.begin; n < rnge.end; ++n) {
//     const auto [i, j, k] = mixed_cells_list[n];
//     IRL::Pt x0(mesh.x(i), mesh.y(j), mesh.z(k));
//     IRL::Pt x1(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
//     IRL::RectangularCuboid cell =
//         IRL::RectangularCuboid::fromBoundingPts(x0, x1);

//     IRL::ImplicitSurfaceCutter<IRL::Sphere<double, MAX_REFINE>,
//                                IRL::GeneralMoments3D<VM_ORDER>>
//         cutter(sphere, cell);

//     auto vol = cutter.computeVolumeMoments();
//     auto surf = cutter.template computeSurfaceMoments<SM_ORDER>(
//         /*useAdaptive=*/false, Eigen::Integrator<double, 2>::GaussKronrod15,
//         /*npts=*/5);

//     (*moments)(i, j, k).first = vol;
//     (*moments)(i, j, k).second = surf;
//     local.push_back({i, j, k, std::move(vol), std::move(surf)});
//   }

// #ifdef USE_MPI
//   MPI_Barrier(MPI_COMM_WORLD);
// #endif
// }