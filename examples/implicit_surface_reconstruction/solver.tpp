// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_SOLVER_TPP_
#define EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_SOLVER_TPP_

#include "examples/implicit_surface_reconstruction/solver.h"

template <class ReturnType>
ReturnType recenterMoments(const ReturnType& moments, const IRL::Pt& xc) {
  ReturnType centered_moments;

  Eigen::Vector3d x_ref(xc[0], xc[1], xc[2]);

  double M0 = moments[0];
  Eigen::Vector3d M1(moments[1], moments[2], moments[3]);
  Eigen::Matrix3d M2 = Eigen::Matrix3d::Zero();
  M2(0, 0) = moments[4];
  M2(0, 1) = moments[5];
  M2(0, 2) = moments[6];
  M2(1, 0) = M2(0, 1);
  M2(1, 1) = moments[7];
  M2(1, 2) = moments[8];
  M2(2, 0) = M2(0, 2);
  M2(2, 1) = M2(1, 2);
  M2(2, 2) = moments[9];

  Eigen::Vector3d M1_recentered = M1 - (M0 * x_ref);
  Eigen::Matrix3d M2_recentered = M2 - (x_ref * M1.transpose()) -
                                  (M1 * x_ref.transpose()) +
                                  (M0 * x_ref * x_ref.transpose());

  centered_moments[0] = M0;
  centered_moments[1] = M1_recentered[0];
  centered_moments[2] = M1_recentered[1];
  centered_moments[3] = M1_recentered[2];
  centered_moments[4] = M2_recentered(0, 0);
  centered_moments[5] = M2_recentered(0, 1);
  centered_moments[6] = M2_recentered(0, 2);
  centered_moments[7] = M2_recentered(1, 1);
  centered_moments[8] = M2_recentered(1, 2);
  centered_moments[9] = M2_recentered(2, 2);

  return centered_moments;
}

template <std::size_t MAX_REFINE, std::size_t VM_ORDER, std::size_t SM_ORDER>
Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
               IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
initializeMoments(const BasicMesh& mesh) {
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  Data<int> cell_status(&mesh);
  std::vector<std::tuple<int, int, int>> mixed_cells_list;
  if (rank == 0) {
    mixed_cells_list = getCellStatus(&cell_status);
  }

  Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                 IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
      moments(&mesh);
  getInitializedField<MAX_REFINE, VM_ORDER, SM_ORDER>(
      cell_status, mixed_cells_list, &moments);

  // writing moments to binary
  if (rank == 0) {
    const int Nx = mesh.getNx();
    std::string outdir =
        "/home/parinht2/Desktop/ppic paper/reconstruction_convergence";
    std::string filename = outdir + "/moments_Nx" + std::to_string(Nx) + ".bin";
    writeMomentsToBinary<VM_ORDER, SM_ORDER>(moments, filename);
  }

  return moments;
}

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
struct ReconRec {
  int i, j, k;
  static constexpr int NV =
      (int)((VM_ORDER + 1) * (VM_ORDER + 2) * (VM_ORDER + 3) / 6);
  static constexpr int NS =
      (int)((SM_ORDER + 1) * (SM_ORDER + 2) * (SM_ORDER + 3) / 6);
  double vol[NV];
  double surf[NS];
};

template <std::size_t MAX_REFINE, std::size_t VM_ORDER, std::size_t SM_ORDER>
void computeReconstructionMetrics(const int& factor,
                                  const std::string& reconstruction_method,
                                  const BasicMesh& mesh) {
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // ------------------ reeading binary file -------------------------------
  Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                 IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
      moments_fine(&mesh);
  std::string filename =
      "/home/parinht2/Desktop/ppic "
      "paper/reconstruction_convergence/sphere_moments_Nx128.bin";

  if (rank == 0) {
    auto t0 = std::chrono::steady_clock::now();
    readMomentsFromBinary<VM_ORDER, SM_ORDER>(filename, &moments_fine);
    auto t1 = std::chrono::steady_clock::now();
    std::printf("time taken to read binary file: %20.6f s\n",
                std::chrono::duration<double>(t1 - t0).count());
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // ------------------ coarsening the data -------------------------------
  const int Nx_fine = mesh.getNx();
  const int Nx_coarse = Nx_fine / factor;

  BasicMesh mesh_coarse(Nx_coarse, Nx_coarse, Nx_coarse, 1);
  mesh_coarse.setCellBoundaries(IRL::Pt(-0.18, -0.18, -0.18),
                                IRL::Pt(0.18, 0.18, 0.18));

  Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                 IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
      moments_coarse(&mesh_coarse);

  if (rank == 0) {
    auto t0 = std::chrono::steady_clock::now();
    coarsenMomentsFromBinary<VM_ORDER, SM_ORDER>(filename, factor,
                                                 &moments_coarse);
    auto t1 = std::chrono::steady_clock::now();
    std::printf("time taken to coarsen moments: %20.6f s\n",
                std::chrono::duration<double>(t1 - t0).count());
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // -------- interface reconstruction from coarse moments  ------------
  Data<double> velU(&mesh_coarse), velV(&mesh_coarse), velW(&mesh_coarse);
  Data<IRL::VolumeMoments> liq_moments(&mesh_coarse), gas_moments(&mesh_coarse);
  Data<IRL::SeparatorVariant> interface(&mesh_coarse);

  if (rank == 0) {
    auto t0 = std::chrono::steady_clock::now();

    for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
      for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
        for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
          IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j), mesh_coarse.z(k));
          IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
                     mesh_coarse.z(k + 1));
          IRL::RectangularCuboid cell =
              IRL::RectangularCuboid::fromBoundingPts(x0, x1);

          double m0 = moments_coarse(i, j, k).first[0];
          IRL::Pt m1(moments_coarse(i, j, k).first[1],
                     moments_coarse(i, j, k).first[2],
                     moments_coarse(i, j, k).first[3]);

          liq_moments(i, j, k) = IRL::VolumeMoments(m0, m1);
          gas_moments(i, j, k) = cell.calculateMoments() - liq_moments(i, j, k);
        }
      }
    }

    getReconstruction(reconstruction_method, liq_moments, gas_moments, 0.0,
                      velU, velV, velW, &interface);

    auto t1 = std::chrono::steady_clock::now();
    std::printf("time taken to reconstruct interface: %20.6f s\n",
                std::chrono::duration<double>(t1 - t0).count());
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // ------- moments of inside cells and build mixed list ---------
  Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                 IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
      moments_reconstruction(&mesh_coarse);

  std::vector<std::tuple<int, int, int>> mixed_cells_list;

  if (rank == 0) {
    auto t0 = std::chrono::steady_clock::now();

    for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
      for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
        for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
          double vf = liq_moments(i, j, k).volume() / mesh_coarse.cell_volume();

          if (vf > IRL::global_constants::VF_HIGH) {
            const IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j),
                             mesh_coarse.z(k));
            const IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
                             mesh_coarse.z(k + 1));
            const IRL::RectangularCuboid cell =
                IRL::RectangularCuboid::fromBoundingPts(x0, x1);
            moments_reconstruction(i, j, k).first =
                IRL::getVolumeMoments<IRL::GeneralMoments3D<VM_ORDER>>(cell);
          }

          if (vf > IRL::global_constants::VF_LOW &&
              vf < IRL::global_constants::VF_HIGH) {
            mixed_cells_list.emplace_back(i, j, k);
          }
        }
      }
    }

    auto t1 = std::chrono::steady_clock::now();
    std::printf(
        "time taken to compute moments for cells inside and finding mixed "
        "list: %20.6f s\n",
        std::chrono::duration<double>(t1 - t0).count());
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // --------------- mixed cell moments using mpi ------------------------
  using VM = IRL::GeneralMoments3D<VM_ORDER>;
  using SM = IRL::GeneralSurfaceMoments3D<SM_ORDER>;
  using VMS = IRL::AddSurfaceOutput<IRL::VolumeMoments,
                                    IRL::ParaboloidParametrizedSurfaceOutput>;
  using RecOut = ReconRec<VM_ORDER, SM_ORDER>;

  auto t0 = std::chrono::steady_clock::now();

  // 0) Count mixed cells (rank 0), build I/J/K arrays, and serialize all
  // mixed-cell paraboloids contiguously.
  int Nmixed = 0;
  std::vector<int> I, J, K;  // (i,j,k) for mixed cells, global order
  IRL::ByteBuffer buf_all;   // all mixed-cell paraboloids, back-to-back
  int size_paraboloid = 0;   // fixed byte-size per paraboloid

  if (rank == 0) {
    // determine fixed serialized size of a paraboloid (one-shot)
    {
      IRL::Paraboloid dummy;
      IRL::ByteBuffer probe;
      probe.resize(0);
      probe.resetBufferPointer();
      IRL::serializeAndPack(dummy, &probe);
      size_paraboloid = probe.size();
    }

    // scan grid once: collect mixed (i,j,k) + pack paraboloid for each
    // (we already have mixed_cells_list from your prior step; use it to avoid
    // rescan if you prefer)
    for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); ++i) {
      for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); ++j) {
        for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); ++k) {
          const double vf =
              liq_moments(i, j, k).volume() / mesh_coarse.cell_volume();
          if (vf > IRL::global_constants::VF_LOW &&
              vf < IRL::global_constants::VF_HIGH) {
            // must be a paraboloid in the interface for mixed cells
            const auto* parab =
                std::get_if<IRL::Paraboloid>(&interface(i, j, k));
            if (!parab) continue;  // (or throw)
            I.push_back(i);
            J.push_back(j);
            K.push_back(k);
            ++Nmixed;
          }
        }
      }
    }

    // Pre-size a single contiguous buffer to hold all serialized paraboloids
    buf_all.resize(Nmixed * size_paraboloid);
    buf_all.resetBufferPointer();
    for (int n = 0; n < Nmixed; ++n) {
      const int i = I[n], j = J[n], k = K[n];
      const auto& parab = std::get<IRL::Paraboloid>(interface(i, j, k));
      IRL::serializeAndPack(parab, &buf_all);
    }
  }

  // 1) Broadcast counts / sizes so all ranks know how many bytes to receive.
  MPI_Bcast(&Nmixed, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&size_paraboloid, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> counts(size), displs(size);
  if (rank == 0) {
    displs[0] = 0;
    for (int r = 0; r < size; ++r) {
      auto rr = block_partition(Nmixed, r, size);
      counts[r] = rr.end - rr.begin;
      if (r > 0) displs[r] = displs[r - 1] + counts[r - 1];
    }
  }

  const auto bounds = block_partition(Nmixed, rank, size);
  const int myN = bounds.end - bounds.begin;

  std::vector<int> Ii(myN), Jj(myN), Kk(myN);
  MPI_Scatterv(rank == 0 ? I.data() : nullptr, counts.data(), displs.data(),
               MPI_INT, Ii.data(), myN, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(rank == 0 ? J.data() : nullptr, counts.data(), displs.data(),
               MPI_INT, Jj.data(), myN, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(rank == 0 ? K.data() : nullptr, counts.data(), displs.data(),
               MPI_INT, Kk.data(), myN, MPI_INT, 0, MPI_COMM_WORLD);

  // 4) Scatter the serialized paraboloids as raw bytes (one big contiguous
  // buffer).
  std::vector<int> countsB(size), displsB(size);
  if (rank == 0) {
    for (int r = 0; r < size; ++r) {
      countsB[r] = counts[r] * size_paraboloid;
      displsB[r] = displs[r] * size_paraboloid;
    }
  }
  IRL::ByteBuffer buf_local;
  buf_local.resize(myN * size_paraboloid);
  buf_local.resetBufferPointer();

  MPI_Scatterv(rank == 0 ? buf_all.data() : nullptr, countsB.data(),
               displsB.data(), MPI_BYTE, buf_local.data(),
               myN * size_paraboloid, MPI_BYTE, 0, MPI_COMM_WORLD);

  std::vector<RecOut> local;
  local.reserve(std::max(1, myN));

  // reset read pointer for sequential unpack (matches pack order during
  // scatter)
  buf_local.resetBufferPointer();
  for (int t = 0; t < myN; ++t) {
    IRL::Paraboloid P;
    IRL::unpackAndStore(&P, &buf_local);  // <-- same API as your code

    const int i = Ii[t], j = Jj[t], k = Kk[t];

    const IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j), mesh_coarse.z(k));
    const IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
                     mesh_coarse.z(k + 1));
    const IRL::RectangularCuboid cell =
        IRL::RectangularCuboid::fromBoundingPts(x0, x1);

    auto surface = IRL::getVolumeMoments<VMS>(cell, P).getSurface();
    auto vm = IRL::getVolumeMoments<VM>(cell, P);
    auto sm = surface.template getSurfaceMoments<SM_ORDER>();

    RecOut r{};
    r.i = i;
    r.j = j;
    r.k = k;
    for (int n = 0; n < RecOut::NV; ++n) r.vol[n] = vm[n];
    for (int n = 0; n < RecOut::NS; ++n) r.surf[n] = sm[n];
    local.push_back(r);
  }

  // 6) Gather results to rank 0 and write into moments_reconstruction
  const int rec_bytes = (int)sizeof(RecOut);
  int local_n = (int)local.size();

  std::vector<int> rcnts(size), rdisp(size), rcntsB(size), rdispB(size);
  MPI_Gather(&local_n, 1, MPI_INT, rcnts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<RecOut> all;
  int total = 0;
  if (rank == 0) {
    rdisp[0] = 0;
    for (int r = 1; r < size; ++r) rdisp[r] = rdisp[r - 1] + rcnts[r - 1];
    total = rdisp[size - 1] + rcnts[size - 1];
    all.resize(total);
    for (int r = 0; r < size; ++r) {
      rcntsB[r] = rcnts[r] * rec_bytes;
      rdispB[r] = rdisp[r] * rec_bytes;
    }
  }

  MPI_Gatherv(local.data(), local_n * rec_bytes, MPI_BYTE,
              rank == 0 ? all.data() : nullptr,
              rank == 0 ? rcntsB.data() : nullptr,
              rank == 0 ? rdispB.data() : nullptr, MPI_BYTE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    for (const auto& rec : all) {
      auto& dst = moments_reconstruction(rec.i, rec.j, rec.k);
      for (int n = 0; n < RecOut::NV; ++n) dst.first[n] = rec.vol[n];
      for (int n = 0; n < RecOut::NS; ++n) dst.second[n] = rec.surf[n];
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  auto t1 = std::chrono::steady_clock::now();
  if (rank == 0) {
    std::printf(
        "time taken to compute moments using interface [MPI/pack]: %20.6f s\n",
        std::chrono::duration<double>(t1 - t0).count());
  }

  // ------------ recentering moments and finding metrics  --------------------
  if (rank == 0) {
    auto t0 = std::chrono::steady_clock::now();

    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
        moments_reconstruction_recentered(&mesh_coarse);
    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
        moments_coarse_recentered(&mesh_coarse);

    for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
      for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
        for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
          const IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j),
                           mesh_coarse.z(k));
          const IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
                           mesh_coarse.z(k + 1));
          const IRL::Pt xc = 0.5 * (x0 + x1);

          IRL::GeneralMoments3D<VM_ORDER> vm = moments_coarse(i, j, k).first;
          IRL::GeneralSurfaceMoments3D<SM_ORDER> sm =
              moments_coarse(i, j, k).second;
          IRL::GeneralMoments3D<VM_ORDER> vmr =
              moments_reconstruction(i, j, k).first;
          IRL::GeneralSurfaceMoments3D<SM_ORDER> smr =
              moments_reconstruction(i, j, k).second;

          moments_reconstruction_recentered(i, j, k).first =
              recenterMoments<IRL::GeneralMoments3D<VM_ORDER>>(vmr, xc);
          moments_reconstruction_recentered(i, j, k).second =
              recenterMoments<IRL::GeneralSurfaceMoments3D<SM_ORDER>>(smr, xc);

          moments_coarse_recentered(i, j, k).first =
              recenterMoments<IRL::GeneralMoments3D<VM_ORDER>>(vm, xc);
          moments_coarse_recentered(i, j, k).second =
              recenterMoments<IRL::GeneralSurfaceMoments3D<SM_ORDER>>(sm, xc);
        }
      }
    }

    auto t1 = std::chrono::steady_clock::now();
    std::printf("time taken to recenter moments: %20.6f s\n",
                std::chrono::duration<double>(t1 - t0).count());

    // metrics
    MomentDiffNorms norms = compute_moment_diff_norms<VM_ORDER, SM_ORDER>(
        moments_coarse_recentered, moments_reconstruction_recentered);

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
  MPI_Barrier(MPI_COMM_WORLD);
}

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_SOLVER_TPP_