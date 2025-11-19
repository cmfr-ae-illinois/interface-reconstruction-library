
#ifndef EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_RECONSTRUCTION_METRICS_TPP_
#define EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_RECONSTRUCTION_METRICS_TPP_

#include "examples/implicit_surface_reconstruction/reconstruction.h"

// computing norms of moment metrics -----------------------------------------
template <std::size_t VM_ORDER, std::size_t SM_ORDER>
MomentDiffNorms compute_moment_diff_norms(
    const Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                         IRL::GeneralSurfaceMoments3D<SM_ORDER>>>& A,
    const Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                         IRL::GeneralSurfaceMoments3D<SM_ORDER>>>& B) {
  const BasicMesh& mA = A.getMesh();
  const BasicMesh& mB = B.getMesh();
  if (mA.imin() != mB.imin() || mA.imax() != mB.imax() ||
      mA.jmin() != mB.jmin() || mA.jmax() != mB.jmax() ||
      mA.kmin() != mB.kmin() || mA.kmax() != mB.kmax()) {
    throw std::runtime_error("Meshes do not match for A and B.");
  }

  // Accumulators
  MomentDiffNorms out{};
  out.vol_M0_Linf = out.vol_M1_Linf = out.vol_M2_Linf = 0.0;
  out.surf_M0_Linf = out.surf_M1_Linf = out.surf_M2_Linf = 0.0;

  long long Ncells = 0;
  long double vol_M0_L2_acc = 0.0L, vol_M1_L2_acc = 0.0L, vol_M2_L2_acc = 0.0L;
  long double surf_M0_L2_acc = 0.0L, surf_M1_L2_acc = 0.0L,
              surf_M2_L2_acc = 0.0L;

  for (int i = mA.imin(); i <= mA.imax(); ++i) {
    for (int j = mA.jmin(); j <= mA.jmax(); ++j) {
      for (int k = mA.kmin(); k <= mA.kmax(); ++k) {
        const auto& a = A(i, j, k);
        const auto& b = B(i, j, k);

        // M0
        const double dM0v = std::abs(a.first[M0] - b.first[M0]);
        out.vol_M0_Linf = std::max(out.vol_M0_Linf, dM0v);
        vol_M0_L2_acc += static_cast<long double>(dM0v) * dM0v;

        // M1
        const double dM1xv = a.first[M1x] - b.first[M1x];
        const double dM1yv = a.first[M1y] - b.first[M1y];
        const double dM1zv = a.first[M1z] - b.first[M1z];
        const double nM1v =
            std::sqrt(dM1xv * dM1xv + dM1yv * dM1yv + dM1zv * dM1zv);
        out.vol_M1_Linf = std::max(out.vol_M1_Linf, nM1v);
        vol_M1_L2_acc += static_cast<long double>(nM1v) * nM1v;

        // M2
        const double dMxxv = a.first[Mxx] - b.first[Mxx];
        const double dMxyv = a.first[Mxy] - b.first[Mxy];
        const double dMxzv = a.first[Mxz] - b.first[Mxz];
        const double dMyyv = a.first[Myy] - b.first[Myy];
        const double dMyzv = a.first[Myz] - b.first[Myz];
        const double dMzzv = a.first[Mzz] - b.first[Mzz];
        const double nM2v_sq =
            dMxxv * dMxxv + dMyyv * dMyyv + dMzzv * dMzzv +
            2.0 * (dMxyv * dMxyv + dMxzv * dMxzv + dMyzv * dMyzv);
        const double nM2v = std::sqrt(nM2v_sq);
        out.vol_M2_Linf = std::max(out.vol_M2_Linf, nM2v);
        vol_M2_L2_acc += static_cast<long double>(nM2v_sq);

        // surface
        const double dM0s = std::abs(a.second[M0] - b.second[M0]);
        out.surf_M0_Linf = std::max(out.surf_M0_Linf, dM0s);
        surf_M0_L2_acc += static_cast<long double>(dM0s) * dM0s;

        const double dM1xs = a.second[M1x] - b.second[M1x];
        const double dM1ys = a.second[M1y] - b.second[M1y];
        const double dM1zs = a.second[M1z] - b.second[M1z];
        const double nM1s =
            std::sqrt(dM1xs * dM1xs + dM1ys * dM1ys + dM1zs * dM1zs);
        out.surf_M1_Linf = std::max(out.surf_M1_Linf, nM1s);
        surf_M1_L2_acc += static_cast<long double>(nM1s) * nM1s;

        const double dMxxs = a.second[Mxx] - b.second[Mxx];
        const double dMxys = a.second[Mxy] - b.second[Mxy];
        const double dMxzs = a.second[Mxz] - b.second[Mxz];
        const double dMyys = a.second[Myy] - b.second[Myy];
        const double dMyzs = a.second[Myz] - b.second[Myz];
        const double dMzzs = a.second[Mzz] - b.second[Mzz];
        const double nM2s_sq =
            dMxxs * dMxxs + dMyys * dMyys + dMzzs * dMzzs +
            2.0 * (dMxys * dMxys + dMxzs * dMxzs + dMyzs * dMyzs);
        const double nM2s = std::sqrt(nM2s_sq);
        out.surf_M2_Linf = std::max(out.surf_M2_Linf, nM2s);
        surf_M2_L2_acc += static_cast<long double>(nM2s_sq);

        ++Ncells;
      }
    }
  }

  // L2
  const long double invN =
      (Ncells > 0) ? 1.0L / static_cast<long double>(Ncells) : 0.0L;

  out.vol_M0_L2 = std::sqrt(static_cast<double>(vol_M0_L2_acc * invN));
  out.vol_M1_L2 = std::sqrt(static_cast<double>(vol_M1_L2_acc * invN));
  out.vol_M2_L2 = std::sqrt(static_cast<double>(vol_M2_L2_acc * invN));
  out.surf_M0_L2 = std::sqrt(static_cast<double>(surf_M0_L2_acc * invN));
  out.surf_M1_L2 = std::sqrt(static_cast<double>(surf_M1_L2_acc * invN));
  out.surf_M2_L2 = std::sqrt(static_cast<double>(surf_M2_L2_acc * invN));

  return out;
}

// recentering moments ---------------------------------------------------
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

// computing metrics by reading binary file -----------------------------
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

// template <std::size_t VM_ORDER, std::size_t SM_ORDER>
// MomentDiffNorms computeReconstructionMetricsFromBin(
//     const int& factor, const std::string& reconstruction_method,
//     const std::string& shape, int Nx_fine, const std::string& output_dir) {
//   int rank = 0, size = 1;
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   MPI_Comm_size(MPI_COMM_WORLD, &size);

//   BasicMesh mesh_fine(Nx_fine, Nx_fine, Nx_fine, 1);
//   SurfaceVariant surf = makeSurface(shape, mesh_fine);
//   (void)surf;

//   const std::string bin_path = binary_filename(output_dir, shape, Nx_fine);

//   // ------------------ coarsening the data -------------------------------
//   const int Nx_coarse = Nx_fine / factor;

//   BasicMesh mesh_coarse(Nx_coarse, Nx_coarse, Nx_coarse, 1);
//   mesh_coarse.setCellBoundaries(
//       IRL::Pt(mesh_fine.x(mesh_fine.imin()), mesh_fine.y(mesh_fine.jmin()),
//               mesh_fine.z(mesh_fine.kmin())),
//       IRL::Pt(mesh_fine.x(mesh_fine.imax() + 1),
//               mesh_fine.y(mesh_fine.jmax() + 1),
//               mesh_fine.z(mesh_fine.kmax() + 1)));

//   Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
//                  IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
//       moments_coarse(&mesh_coarse);

//   if (rank == 0) {
//     auto t0 = std::chrono::steady_clock::now();
//     coarsenMomentsFromBinary<VM_ORDER, SM_ORDER>(bin_path, factor,
//                                                  &moments_coarse);
//     auto t1 = std::chrono::steady_clock::now();
//     std::printf("time taken to coarsen moments: %20.6f s\n",
//                 std::chrono::duration<double>(t1 - t0).count());
//   }
//   MPI_Barrier(MPI_COMM_WORLD);

//   // -------- interface reconstruction from coarse moments  ------------
//   Data<double> velU(&mesh_coarse), velV(&mesh_coarse), velW(&mesh_coarse);
//   Data<IRL::VolumeMoments> liq_moments(&mesh_coarse),
//   gas_moments(&mesh_coarse); Data<IRL::SeparatorVariant> interface(
//       &mesh_coarse);  // PUT THIS WITHIN RANK LOOP

//   if (rank == 0) {
//     auto t0 = std::chrono::steady_clock::now();

//     for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
//       for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
//         for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
//           IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j), mesh_coarse.z(k));
//           IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
//                      mesh_coarse.z(k + 1));
//           IRL::RectangularCuboid cell =
//               IRL::RectangularCuboid::fromBoundingPts(x0, x1);

//           double m0 = moments_coarse(i, j, k).first[0];
//           IRL::Pt m1(moments_coarse(i, j, k).first[1],
//                      moments_coarse(i, j, k).first[2],
//                      moments_coarse(i, j, k).first[3]);

//           liq_moments(i, j, k) = IRL::VolumeMoments(m0, m1);
//           gas_moments(i, j, k) = cell.calculateMoments() - liq_moments(i, j,
//           k);
//         }
//       }
//     }

//     getReconstruction(reconstruction_method, liq_moments, gas_moments, 0.0,
//                       velU, velV, velW, &interface);

//     auto t1 = std::chrono::steady_clock::now();
//     std::printf("time taken to reconstruct interface: %20.6f s\n",
//                 std::chrono::duration<double>(t1 - t0).count());
//   }
//   MPI_Barrier(MPI_COMM_WORLD);

//   // ------- moments of inside cells and build mixed list ---------
//   Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
//                  IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
//       moments_reconstruction(&mesh_coarse);

//   if (rank == 0) {
//     auto t0 = std::chrono::steady_clock::now();

//     for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
//       for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
//         for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
//           double vf = liq_moments(i, j, k).volume() /
//           mesh_coarse.cell_volume();

//           if (vf > IRL::global_constants::VF_HIGH) {
//             const IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j),
//                              mesh_coarse.z(k));
//             const IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
//                              mesh_coarse.z(k + 1));
//             const IRL::RectangularCuboid cell =
//                 IRL::RectangularCuboid::fromBoundingPts(x0, x1);
//             moments_reconstruction(i, j, k).first =
//                 IRL::getVolumeMoments<IRL::GeneralMoments3D<VM_ORDER>>(cell);
//           }
//         }
//       }
//     }

//     auto t1 = std::chrono::steady_clock::now();
//     std::printf(
//         "time taken to compute moments for cells inside and finding mixed "
//         "list: %20.6f s\n",
//         std::chrono::duration<double>(t1 - t0).count());
//   }
//   MPI_Barrier(MPI_COMM_WORLD);

//   // --------------- mixed cell moments using mpi ------------------------
//   using VM = IRL::GeneralMoments3D<VM_ORDER>;
//   using SM = IRL::GeneralSurfaceMoments3D<SM_ORDER>;
//   using VMS = IRL::AddSurfaceOutput<IRL::VolumeMoments,
//                                     IRL::ParaboloidParametrizedSurfaceOutput>;
//   using RecOut = ReconRec<VM_ORDER, SM_ORDER>;

//   auto t0 = std::chrono::steady_clock::now();

//   // 0) Count mixed cells (rank 0), build I/J/K arrays, and serialize all
//   // mixed-cell paraboloids contiguously.
//   int Nmixed = 0;
//   std::vector<int> I, J, K;  // (i,j,k) for mixed cells, global order
//   IRL::ByteBuffer buf_all;   // all mixed-cell paraboloids, back-to-back
//   int size_paraboloid = 0;   // fixed byte-size per paraboloid

//   if (rank == 0) {
//     // determine fixed serialized size of a paraboloid (one-shot)
//     {
//       IRL::Paraboloid dummy;
//       IRL::ByteBuffer probe;
//       probe.resize(0);
//       probe.resetBufferPointer();
//       IRL::serializeAndPack(dummy, &probe);
//       size_paraboloid = probe.size();
//     }

//     for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); ++i) {
//       for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); ++j) {
//         for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); ++k) {
//           const double vf =
//               liq_moments(i, j, k).volume() / mesh_coarse.cell_volume();
//           if (vf > IRL::global_constants::VF_LOW &&
//               vf < IRL::global_constants::VF_HIGH) {
//             // must be a paraboloid in the interface for mixed cells
//             const auto* parab =
//                 std::get_if<IRL::Paraboloid>(&interface(i, j, k));
//             if (!parab) continue;  // (or throw)
//             I.push_back(i);
//             J.push_back(j);
//             K.push_back(k);
//             ++Nmixed;
//           }
//         }
//       }
//     }

//     // Pre-size a single contiguous buffer to hold all serialized paraboloids
//     buf_all.resize(Nmixed * size_paraboloid);
//     buf_all.resetBufferPointer();
//     for (int n = 0; n < Nmixed; ++n) {
//       const int i = I[n], j = J[n], k = K[n];
//       const auto& parab = std::get<IRL::Paraboloid>(interface(i, j, k));
//       IRL::serializeAndPack(parab, &buf_all);
//     }
//   }

//   // 1) Broadcast counts / sizes so all ranks know how many bytes to receive.
//   MPI_Bcast(&Nmixed, 1, MPI_INT, 0, MPI_COMM_WORLD);
//   MPI_Bcast(&size_paraboloid, 1, MPI_INT, 0, MPI_COMM_WORLD);

//   std::vector<int> counts(size), displs(size);
//   if (rank == 0) {
//     displs[0] = 0;
//     for (int r = 0; r < size; ++r) {
//       auto rr = block_partition(Nmixed, r, size);
//       counts[r] = rr.end - rr.begin;
//       if (r > 0) displs[r] = displs[r - 1] + counts[r - 1];
//     }
//   }

//   const auto bounds = block_partition(Nmixed, rank, size);
//   const int myN = bounds.end - bounds.begin;

//   std::vector<int> Ii(myN), Jj(myN), Kk(myN);
//   MPI_Scatterv(rank == 0 ? I.data() : nullptr, counts.data(), displs.data(),
//                MPI_INT, Ii.data(), myN, MPI_INT, 0, MPI_COMM_WORLD);
//   MPI_Scatterv(rank == 0 ? J.data() : nullptr, counts.data(), displs.data(),
//                MPI_INT, Jj.data(), myN, MPI_INT, 0, MPI_COMM_WORLD);
//   MPI_Scatterv(rank == 0 ? K.data() : nullptr, counts.data(), displs.data(),
//                MPI_INT, Kk.data(), myN, MPI_INT, 0, MPI_COMM_WORLD);

//   // 4) Scatter the serialized paraboloids as raw bytes (one big contiguous
//   // buffer).
//   std::vector<int> countsB(size), displsB(size);
//   if (rank == 0) {
//     for (int r = 0; r < size; ++r) {
//       countsB[r] = counts[r] * size_paraboloid;
//       displsB[r] = displs[r] * size_paraboloid;
//     }
//   }
//   IRL::ByteBuffer buf_local;
//   buf_local.resize(myN * size_paraboloid);
//   buf_local.resetBufferPointer();

//   MPI_Scatterv(rank == 0 ? buf_all.data() : nullptr, countsB.data(),
//                displsB.data(), MPI_BYTE, buf_local.data(),
//                myN * size_paraboloid, MPI_BYTE, 0, MPI_COMM_WORLD);

//   std::vector<RecOut> local;
//   local.reserve(std::max(1, myN));

//   // reset read pointer for sequential unpack (matches pack order during
//   // scatter)
//   buf_local.resetBufferPointer();
//   for (int t = 0; t < myN; ++t) {
//     IRL::Paraboloid P;
//     IRL::unpackAndStore(&P, &buf_local);  // <-- same API as your code

//     const int i = Ii[t], j = Jj[t], k = Kk[t];

//     const IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j), mesh_coarse.z(k));
//     const IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
//                      mesh_coarse.z(k + 1));
//     const IRL::RectangularCuboid cell =
//         IRL::RectangularCuboid::fromBoundingPts(x0, x1);

//     auto surface = IRL::getVolumeMoments<VMS>(cell, P).getSurface();
//     auto vm = IRL::getVolumeMoments<VM>(cell, P);
//     auto sm = surface.template getSurfaceMoments<SM_ORDER>();

//     RecOut r{};
//     r.i = i;
//     r.j = j;
//     r.k = k;
//     for (int n = 0; n < RecOut::NV; ++n) r.vol[n] = vm[n];
//     for (int n = 0; n < RecOut::NS; ++n) r.surf[n] = sm[n];
//     local.push_back(r);
//   }

//   // 6) Gather results to rank 0 and write into moments_reconstruction
//   const int rec_bytes = (int)sizeof(RecOut);
//   int local_n = (int)local.size();

//   std::vector<int> rcnts(size), rdisp(size), rcntsB(size), rdispB(size);
//   MPI_Gather(&local_n, 1, MPI_INT, rcnts.data(), 1, MPI_INT, 0,
//   MPI_COMM_WORLD);

//   std::vector<RecOut> all;
//   int total = 0;
//   if (rank == 0) {
//     rdisp[0] = 0;
//     for (int r = 1; r < size; ++r) rdisp[r] = rdisp[r - 1] + rcnts[r - 1];
//     total = rdisp[size - 1] + rcnts[size - 1];
//     all.resize(total);
//     for (int r = 0; r < size; ++r) {
//       rcntsB[r] = rcnts[r] * rec_bytes;
//       rdispB[r] = rdisp[r] * rec_bytes;
//     }
//   }

//   MPI_Gatherv(local.data(), local_n * rec_bytes, MPI_BYTE,
//               rank == 0 ? all.data() : nullptr,
//               rank == 0 ? rcntsB.data() : nullptr,
//               rank == 0 ? rdispB.data() : nullptr, MPI_BYTE, 0,
//               MPI_COMM_WORLD);

//   if (rank == 0) {
//     for (const auto& rec : all) {
//       auto& dst = moments_reconstruction(rec.i, rec.j, rec.k);
//       for (int n = 0; n < RecOut::NV; ++n) dst.first[n] = rec.vol[n];
//       for (int n = 0; n < RecOut::NS; ++n) dst.second[n] = rec.surf[n];
//     }
//   }

//   MPI_Barrier(MPI_COMM_WORLD);
//   auto t1 = std::chrono::steady_clock::now();
//   if (rank == 0) {
//     std::printf(
//         "time taken to compute moments using interface [MPI/pack]: %20.6f
//         s\n", std::chrono::duration<double>(t1 - t0).count());
//   }

//   // ------------ recentering moments and finding metrics
//   -------------------- MomentDiffNorms norms{}; if (rank == 0) {
//     auto t0 = std::chrono::steady_clock::now();

//     Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
//                    IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
//         moments_reconstruction_recentered(&mesh_coarse);
//     Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
//                    IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
//         moments_coarse_recentered(&mesh_coarse);

//     for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
//       for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
//         for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
//           const IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j),
//                            mesh_coarse.z(k));
//           const IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
//                            mesh_coarse.z(k + 1));
//           const IRL::Pt xc = 0.5 * (x0 + x1);

//           IRL::GeneralMoments3D<VM_ORDER> vm = moments_coarse(i, j, k).first;
//           IRL::GeneralSurfaceMoments3D<SM_ORDER> sm =
//               moments_coarse(i, j, k).second;
//           IRL::GeneralMoments3D<VM_ORDER> vmr =
//               moments_reconstruction(i, j, k).first;
//           IRL::GeneralSurfaceMoments3D<SM_ORDER> smr =
//               moments_reconstruction(i, j, k).second;

//           moments_reconstruction_recentered(i, j, k).first =
//               recenterMoments<IRL::GeneralMoments3D<VM_ORDER>>(vmr, xc);
//           moments_reconstruction_recentered(i, j, k).second =
//               recenterMoments<IRL::GeneralSurfaceMoments3D<SM_ORDER>>(smr,
//               xc);

//           moments_coarse_recentered(i, j, k).first =
//               recenterMoments<IRL::GeneralMoments3D<VM_ORDER>>(vm, xc);
//           moments_coarse_recentered(i, j, k).second =
//               recenterMoments<IRL::GeneralSurfaceMoments3D<SM_ORDER>>(sm,
//               xc);
//         }
//       }
//     }

//     auto t1 = std::chrono::steady_clock::now();
//     std::printf("time taken to recenter moments: %20.6f s\n",
//                 std::chrono::duration<double>(t1 - t0).count());

//     // metrics
//     norms = compute_moment_diff_norms<VM_ORDER, SM_ORDER>(
//         moments_coarse_recentered, moments_reconstruction_recentered);

//     std::printf("VOL  : M0  Linf=%.8e  L2=%.8e\n", norms.vol_M0_Linf,
//                 norms.vol_M0_L2);
//     std::printf("VOL  : M1  Linf=%.8e  L2=%.8e\n", norms.vol_M1_Linf,
//                 norms.vol_M1_L2);
//     std::printf("VOL  : M2  Linf=%.8e  L2=%.8e\n", norms.vol_M2_Linf,
//                 norms.vol_M2_L2);
//     std::printf("SURF : M0  Linf=%.8e  L2=%.8e\n", norms.surf_M0_Linf,
//                 norms.surf_M0_L2);
//     std::printf("SURF : M1  Linf=%.8e  L2=%.8e\n", norms.surf_M1_Linf,
//                 norms.surf_M1_L2);
//     std::printf("SURF : M2  Linf=%.8e  L2=%.8e\n", norms.surf_M2_Linf,
//                 norms.surf_M2_L2);
//   }
//   double buf[12];
//   if (rank == 0) {
//     buf[0] = norms.vol_M0_Linf;
//     buf[1] = norms.vol_M0_L2;
//     buf[2] = norms.vol_M1_Linf;
//     buf[3] = norms.vol_M1_L2;
//     buf[4] = norms.vol_M2_Linf;
//     buf[5] = norms.vol_M2_L2;
//     buf[6] = norms.surf_M0_Linf;
//     buf[7] = norms.surf_M0_L2;
//     buf[8] = norms.surf_M1_Linf;
//     buf[9] = norms.surf_M1_L2;
//     buf[10] = norms.surf_M2_Linf;
//     buf[11] = norms.surf_M2_L2;
//   }
//   MPI_Bcast(buf, 12, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//   if (rank != 0) {
//     norms.vol_M0_Linf = buf[0];
//     norms.vol_M0_L2 = buf[1];
//     norms.vol_M1_Linf = buf[2];
//     norms.vol_M1_L2 = buf[3];
//     norms.vol_M2_Linf = buf[4];
//     norms.vol_M2_L2 = buf[5];
//     norms.surf_M0_Linf = buf[6];
//     norms.surf_M0_L2 = buf[7];
//     norms.surf_M1_Linf = buf[8];
//     norms.surf_M1_L2 = buf[9];
//     norms.surf_M2_Linf = buf[10];
//     norms.surf_M2_L2 = buf[11];
//   }

//   return norms;
// }

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
MomentDiffNorms computeReconstructionMetricsFromBin(
    const int& factor, const std::string& reconstruction_method,
    const std::string& shape, int Nx_fine, const std::string& output_dir) {
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  BasicMesh mesh_fine(Nx_fine, Nx_fine, Nx_fine, 1);
  SurfaceVariant surf = makeSurface(shape, mesh_fine);
  (void)surf;

  const std::string bin_path = binary_filename(output_dir, shape, Nx_fine);

  // ------------------ coarsening the data -------------------------------
  const int Nx_coarse = Nx_fine / factor;

  BasicMesh mesh_coarse(Nx_coarse, Nx_coarse, Nx_coarse, 1);
  mesh_coarse.setCellBoundaries(
      IRL::Pt(mesh_fine.x(mesh_fine.imin()), mesh_fine.y(mesh_fine.jmin()),
              mesh_fine.z(mesh_fine.kmin())),
      IRL::Pt(mesh_fine.x(mesh_fine.imax() + 1),
              mesh_fine.y(mesh_fine.jmax() + 1),
              mesh_fine.z(mesh_fine.kmax() + 1)));

  using VM = IRL::GeneralMoments3D<VM_ORDER>;
  using SM = IRL::GeneralSurfaceMoments3D<SM_ORDER>;
  std::unique_ptr<Data<std::pair<VM, SM>>> moments_coarse;

  if (rank == 0) {
    moments_coarse = std::make_unique<Data<std::pair<VM, SM>>>(&mesh_coarse);
    auto t0 = std::chrono::steady_clock::now();
    coarsenMomentsFromBinary<VM_ORDER, SM_ORDER>(bin_path, factor,
                                                 moments_coarse.get());
    auto t1 = std::chrono::steady_clock::now();
    std::printf("time taken to coarsen moments: %20.6f s\n",
                std::chrono::duration<double>(t1 - t0).count());
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // -------- interface reconstruction from coarse moments  ------------
  std::unique_ptr<Data<double>> velU, velV, velW;
  std::unique_ptr<Data<IRL::VolumeMoments>> liq_moments, gas_moments;
  std::unique_ptr<Data<IRL::SeparatorVariant>> interface;

  if (rank == 0) {
    velU = std::make_unique<Data<double>>(&mesh_coarse);
    velV = std::make_unique<Data<double>>(&mesh_coarse);
    velW = std::make_unique<Data<double>>(&mesh_coarse);
    liq_moments = std::make_unique<Data<IRL::VolumeMoments>>(&mesh_coarse);
    gas_moments = std::make_unique<Data<IRL::VolumeMoments>>(&mesh_coarse);
    interface = std::make_unique<Data<IRL::SeparatorVariant>>(&mesh_coarse);

    auto t0 = std::chrono::steady_clock::now();
    for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
      for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
        for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
          IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j), mesh_coarse.z(k));
          IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
                     mesh_coarse.z(k + 1));
          IRL::RectangularCuboid cell =
              IRL::RectangularCuboid::fromBoundingPts(x0, x1);

          double m0 = (*moments_coarse)(i, j, k).first[0];
          IRL::Pt m1((*moments_coarse)(i, j, k).first[1],
                     (*moments_coarse)(i, j, k).first[2],
                     (*moments_coarse)(i, j, k).first[3]);

          (*liq_moments)(i, j, k) = IRL::VolumeMoments(m0, m1);
          (*gas_moments)(i, j, k) =
              cell.calculateMoments() - (*liq_moments)(i, j, k);
        }
      }
    }

    getReconstruction(reconstruction_method, *liq_moments, *gas_moments, 0.0,
                      *velU, *velV, *velW, interface.get());

    auto t1 = std::chrono::steady_clock::now();
    std::printf("time taken to reconstruct interface: %20.6f s\n",
                std::chrono::duration<double>(t1 - t0).count());
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // ------- moments of inside cells and prepare mixed list on root ---------
  std::unique_ptr<Data<std::pair<VM, SM>>> moments_reconstruction;
  if (rank == 0) {
    moments_reconstruction =
        std::make_unique<Data<std::pair<VM, SM>>>(&mesh_coarse);
    auto t0 = std::chrono::steady_clock::now();

    for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
      for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
        for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
          double vf =
              (*liq_moments)(i, j, k).volume() / mesh_coarse.cell_volume();
          if (vf > IRL::global_constants::VF_HIGH) {
            const IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j),
                             mesh_coarse.z(k));
            const IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
                             mesh_coarse.z(k + 1));
            const IRL::RectangularCuboid cell =
                IRL::RectangularCuboid::fromBoundingPts(x0, x1);
            (*moments_reconstruction)(i, j, k).first =
                IRL::getVolumeMoments<VM>(cell);
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
  using VMS = IRL::AddSurfaceOutput<IRL::VolumeMoments,
                                    IRL::ParaboloidParametrizedSurfaceOutput>;
  using RecOut = ReconRec<VM_ORDER, SM_ORDER>;

  auto t0 = std::chrono::steady_clock::now();

  int Nmixed = 0;
  std::vector<int> I, J, K;  // root-only lists of mixed cells
  IRL::ByteBuffer buf_all;   // root-only contiguous paraboloids
  int size_paraboloid = 0;   // fixed size

  if (rank == 0) {
    // Determine serialized size once.
    {
      IRL::Paraboloid dummy;
      IRL::ByteBuffer probe;
      probe.resize(0);
      probe.resetBufferPointer();
      IRL::serializeAndPack(dummy, &probe);
      size_paraboloid = probe.size();
    }

    // Build mixed list and pack
    I.reserve((mesh_coarse.imax() - mesh_coarse.imin() + 1) *
              (mesh_coarse.jmax() - mesh_coarse.jmin() + 1) *
              (mesh_coarse.kmax() - mesh_coarse.kmin() + 1) / 10);  // heuristic
    J.reserve(I.capacity());
    K.reserve(I.capacity());

    for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); ++i) {
      for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); ++j) {
        for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); ++k) {
          const double vf =
              (*liq_moments)(i, j, k).volume() / mesh_coarse.cell_volume();
          if (vf > IRL::global_constants::VF_LOW &&
              vf < IRL::global_constants::VF_HIGH) {
            const auto* parab =
                std::get_if<IRL::Paraboloid>(&(*interface)(i, j, k));
            if (!parab) continue;
            I.push_back(i);
            J.push_back(j);
            K.push_back(k);
          }
        }
      }
    }
    Nmixed = (int)I.size();

    buf_all.resize(Nmixed * size_paraboloid);
    buf_all.resetBufferPointer();
    for (int n = 0; n < Nmixed; ++n) {
      const int i = I[n], j = J[n], k = K[n];
      const auto& parab = std::get<IRL::Paraboloid>((*interface)(i, j, k));
      IRL::serializeAndPack(parab, &buf_all);
    }

    // Free root-only fields we no longer need to save RAM before the MPI stage
    // (interface content has been serialized).
    interface.reset();
    velU.reset();
    velV.reset();
    velW.reset();
    // liq_moments is still used for VF checks above; safe to free now:
    liq_moments.reset();
    gas_moments.reset();
  }

  // Broadcast the counts / sizes
  MPI_Bcast(&Nmixed, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&size_paraboloid, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Partition work
  std::vector<int> counts, displs;
  if (rank == 0) {
    counts.resize(size);
    displs.resize(size);
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
  MPI_Scatterv(rank == 0 ? I.data() : nullptr,
               rank == 0 ? counts.data() : nullptr,
               rank == 0 ? displs.data() : nullptr, MPI_INT, Ii.data(), myN,
               MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(rank == 0 ? J.data() : nullptr,
               rank == 0 ? counts.data() : nullptr,
               rank == 0 ? displs.data() : nullptr, MPI_INT, Jj.data(), myN,
               MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(rank == 0 ? K.data() : nullptr,
               rank == 0 ? counts.data() : nullptr,
               rank == 0 ? displs.data() : nullptr, MPI_INT, Kk.data(), myN,
               MPI_INT, 0, MPI_COMM_WORLD);

  // Scatter serialized paraboloids as bytes
  std::vector<int> countsB, displsB;
  if (rank == 0) {
    countsB.resize(size);
    displsB.resize(size);
    for (int r = 0; r < size; ++r) {
      countsB[r] = counts[r] * size_paraboloid;
      displsB[r] = displs[r] * size_paraboloid;
    }
  }

  IRL::ByteBuffer buf_local;
  buf_local.resize(myN * size_paraboloid);
  buf_local.resetBufferPointer();

  MPI_Scatterv(rank == 0 ? buf_all.data() : nullptr,
               rank == 0 ? countsB.data() : nullptr,
               rank == 0 ? displsB.data() : nullptr, MPI_BYTE, buf_local.data(),
               myN * size_paraboloid, MPI_BYTE, 0, MPI_COMM_WORLD);

  // Free large root-only vectors ASAP after scattering
  if (rank == 0) {
    std::vector<int>().swap(I);
    std::vector<int>().swap(J);
    std::vector<int>().swap(K);
    IRL::ByteBuffer tmp;
    std::swap(buf_all, tmp);
    // IRL::ByteBuffer().swap(buf_all);
    std::vector<int>().swap(counts);
    std::vector<int>().swap(displs);
    std::vector<int>().swap(countsB);
    std::vector<int>().swap(displsB);
  }

  // Local computation
  std::vector<ReconRec<VM_ORDER, SM_ORDER>> local;
  local.reserve(std::max(1, myN));

  buf_local.resetBufferPointer();
  for (int t = 0; t < myN; ++t) {
    IRL::Paraboloid P;
    IRL::unpackAndStore(&P, &buf_local);

    const int i = Ii[t], j = Jj[t], k = Kk[t];

    const IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j), mesh_coarse.z(k));
    const IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
                     mesh_coarse.z(k + 1));
    const IRL::RectangularCuboid cell =
        IRL::RectangularCuboid::fromBoundingPts(x0, x1);

    auto surface = IRL::getVolumeMoments<VMS>(cell, P).getSurface();
    auto vm = IRL::getVolumeMoments<VM>(cell, P);
    auto sm = surface.template getSurfaceMoments<SM_ORDER>();

    ReconRec<VM_ORDER, SM_ORDER> r{};
    r.i = i;
    r.j = j;
    r.k = k;
    for (int n = 0; n < ReconRec<VM_ORDER, SM_ORDER>::NV; ++n) r.vol[n] = vm[n];
    for (int n = 0; n < ReconRec<VM_ORDER, SM_ORDER>::NS; ++n)
      r.surf[n] = sm[n];
    local.push_back(r);
  }

  // Gather structured results to root (others allocate no big recv buffers).
  const int rec_bytes = (int)sizeof(RecOut);
  int local_n = (int)local.size();

  std::vector<int> rcnts, rdisp, rcntsB, rdispB;
  if (rank == 0) {
    rcnts.resize(size);
    rdisp.resize(size);
    rcntsB.resize(size);
    rdispB.resize(size);
  }

  MPI_Gather(&local_n, 1, MPI_INT, rank == 0 ? rcnts.data() : nullptr, 1,
             MPI_INT, 0, MPI_COMM_WORLD);

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

  // Root writes gathered results into moments_reconstruction
  if (rank == 0) {
    for (const auto& rec : all) {
      auto& dst = (*moments_reconstruction)(rec.i, rec.j, rec.k);
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

  // --------------------- recentering & metrics  ------------------------
  MomentDiffNorms norms{};
  if (rank == 0) {
    auto t0m = std::chrono::steady_clock::now();

    Data<std::pair<VM, SM>> moments_reconstruction_recentered(&mesh_coarse);
    Data<std::pair<VM, SM>> moments_coarse_recentered(&mesh_coarse);

    for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
      for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
        for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
          const IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j),
                           mesh_coarse.z(k));
          const IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
                           mesh_coarse.z(k + 1));
          const IRL::Pt xc = 0.5 * (x0 + x1);

          const VM& vm = (*moments_coarse)(i, j, k).first;
          const SM& sm = (*moments_coarse)(i, j, k).second;
          const VM& vmr = (*moments_reconstruction)(i, j, k).first;
          const SM& smr = (*moments_reconstruction)(i, j, k).second;

          moments_reconstruction_recentered(i, j, k).first =
              recenterMoments<VM>(vmr, xc);
          moments_reconstruction_recentered(i, j, k).second =
              recenterMoments<SM>(smr, xc);

          moments_coarse_recentered(i, j, k).first =
              recenterMoments<VM>(vm, xc);
          moments_coarse_recentered(i, j, k).second =
              recenterMoments<SM>(sm, xc);
        }
      }
    }

    auto t1m = std::chrono::steady_clock::now();
    std::printf("time taken to recenter moments: %20.6f s\n",
                std::chrono::duration<double>(t1m - t0m).count());

    norms = compute_moment_diff_norms<VM_ORDER, SM_ORDER>(
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

    moments_reconstruction.reset();
    moments_coarse.reset();
  }

  // Broadcast metrics to all ranks for the return value
  double buf[12];
  if (rank == 0) {
    buf[0] = norms.vol_M0_Linf;
    buf[1] = norms.vol_M0_L2;
    buf[2] = norms.vol_M1_Linf;
    buf[3] = norms.vol_M1_L2;
    buf[4] = norms.vol_M2_Linf;
    buf[5] = norms.vol_M2_L2;
    buf[6] = norms.surf_M0_Linf;
    buf[7] = norms.surf_M0_L2;
    buf[8] = norms.surf_M1_Linf;
    buf[9] = norms.surf_M1_L2;
    buf[10] = norms.surf_M2_Linf;
    buf[11] = norms.surf_M2_L2;
  }
  MPI_Bcast(buf, 12, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    norms.vol_M0_Linf = buf[0];
    norms.vol_M0_L2 = buf[1];
    norms.vol_M1_Linf = buf[2];
    norms.vol_M1_L2 = buf[3];
    norms.vol_M2_Linf = buf[4];
    norms.vol_M2_L2 = buf[5];
    norms.surf_M0_Linf = buf[6];
    norms.surf_M0_L2 = buf[7];
    norms.surf_M1_Linf = buf[8];
    norms.surf_M1_L2 = buf[9];
    norms.surf_M2_Linf = buf[10];
    norms.surf_M2_L2 = buf[11];
  }

  return norms;
}

// convergence for a given shape and reconstruction method --------------------
template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void run_convergence(const std::string& shape, int Nx_fine,
                     const std::string& reconstruction_method,
                     const std::string& output_dir) {
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::string csv_path = output_dir + "/convergence_" + shape + "_" +
                         reconstruction_method + "_Nx" +
                         std::to_string(Nx_fine) + ".csv";

  // factors for convergence
  const std::vector<int> factors = {1, 2, 4, 8, 16};  // 256, 128, 64, 32, 16

  std::ofstream ofs;
  if (rank == 0) {
    ofs.open(csv_path, std::ios::out | std::ios::trunc);
    if (!ofs) {
      std::fprintf(stderr, "ERROR: Cannot open output file: %s\n",
                   csv_path.c_str());
      MPI_Abort(MPI_COMM_WORLD, 3);
    }
    ofs << "shape,Nx_fine,factor,"
           "vol_M0_Linf,vol_M0_L2,vol_M1_Linf,vol_M1_L2,vol_M2_Linf,vol_M2_L2,"
           "surf_M0_Linf,surf_M0_L2,surf_M1_Linf,surf_M1_L2,surf_M2_Linf,surf_"
           "M2_L2\n";
    ofs.flush();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  for (int f : factors) {
    // Ensure Nx_fine divisible by factor on all ranks
    int ok = (Nx_fine % f == 0) ? 1 : 0;
    if (rank == 0 && !ok) {
      std::fprintf(stderr, "ERROR: Nx_fine (%d) not divisible by factor (%d)\n",
                   Nx_fine, f);
    }
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!ok) MPI_Abort(MPI_COMM_WORLD, 4);

    // find norms
    const MomentDiffNorms norms =
        computeReconstructionMetricsFromBin<VM_ORDER, SM_ORDER>(
            f, reconstruction_method, shape, Nx_fine, output_dir);

    if (rank == 0) {
      ofs << shape << "," << Nx_fine << "," << f << ",";
      ofs << std::scientific << std::setprecision(8) << norms.vol_M0_Linf << ","
          << norms.vol_M0_L2 << "," << norms.vol_M1_Linf << ","
          << norms.vol_M1_L2 << "," << norms.vol_M2_Linf << ","
          << norms.vol_M2_L2 << "," << norms.surf_M0_Linf << ","
          << norms.surf_M0_L2 << "," << norms.surf_M1_Linf << ","
          << norms.surf_M1_L2 << "," << norms.surf_M2_Linf << ","
          << norms.surf_M2_L2 << "\n";
      ofs.flush();

      std::printf(
          "[convergence] shape=%s Nx_fine=%d factor=%d  "
          "VOL(M0 Linf=%.3e L2=%.3e, M1 Linf=%.3e L2=%.3e, M2 Linf=%.3e "
          "L2=%.3e)  "
          "SURF(M0 Linf=%.3e L2=%.3e, M1 Linf=%.3e L2=%.3e, M2 Linf=%.3e "
          "L2=%.3e)\n",
          shape.c_str(), Nx_fine, f, norms.vol_M0_Linf, norms.vol_M0_L2,
          norms.vol_M1_Linf, norms.vol_M1_L2, norms.vol_M2_Linf,
          norms.vol_M2_L2, norms.surf_M0_Linf, norms.surf_M0_L2,
          norms.surf_M1_Linf, norms.surf_M1_L2, norms.surf_M2_Linf,
          norms.surf_M2_L2);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (rank == 0) {
    ofs.close();
    std::printf("✅ Convergence results written to %s\n", csv_path.c_str());
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

// outputting interfaces for viz ----------------------------------------------
template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void output_interfaces(const std::string& shape, int Nx_fine, const int& factor,
                       const std::string& reconstruction_method,
                       const std::string& output_dir) {
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0) {
    // setting fine mesh
    BasicMesh mesh_fine(Nx_fine, Nx_fine, Nx_fine, 1);
    SurfaceVariant surf = makeSurface(shape, mesh_fine);
    (void)surf;

    // reading and coarsening moment data from binary
    const std::string bin_path = binary_filename(output_dir, shape, Nx_fine);
    const int Nx_coarse = Nx_fine / factor;
    BasicMesh mesh_coarse(Nx_coarse, Nx_coarse, Nx_coarse, 1);
    mesh_coarse.setCellBoundaries(
        IRL::Pt(mesh_fine.x(mesh_fine.imin()), mesh_fine.y(mesh_fine.jmin()),
                mesh_fine.z(mesh_fine.kmin())),
        IRL::Pt(mesh_fine.x(mesh_fine.imax() + 1),
                mesh_fine.y(mesh_fine.jmax() + 1),
                mesh_fine.z(mesh_fine.kmax() + 1)));
    Data<std::pair<IRL::GeneralMoments3D<VM_ORDER>,
                   IRL::GeneralSurfaceMoments3D<SM_ORDER>>>
        moments_coarse(&mesh_coarse);
    coarsenMomentsFromBinary<VM_ORDER, SM_ORDER>(bin_path, factor,
                                                 &moments_coarse);

    // getting interface
    Data<double> velU(&mesh_coarse), velV(&mesh_coarse), velW(&mesh_coarse);
    Data<IRL::VolumeMoments> liq_moments(&mesh_coarse),
        gas_moments(&mesh_coarse);
    Data<IRL::SeparatorVariant> interface(&mesh_coarse);

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

    // outputting interface
    std::string vtk_output_dir = output_dir + "/" + shape + "_" +
                                 reconstruction_method + "_Nx" +
                                 std::to_string(Nx_coarse);
    VTKOutput vtk_io(vtk_output_dir, "viz", mesh_coarse);
    double sim_time = 0.0;
    writeInterfaceToFile(liq_moments, interface, sim_time, &vtk_io, true);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

// finding curvedness
template <std::size_t VM_ORDER, std::size_t SM_ORDER>
std::pair<double, double> getCurvednessMetrics(
    const std::string& shape, int Nx_fine, const int& factor,
    const std::string& reconstruction_method, const std::string& output_dir) {
  // int rank = 0, size = 1;
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::pair<double, double> norms = {0.0, 0.0};  // Linf, L2

  // if (rank == 0) {
  // ---------------- CSV OUTPUTTING -----------------------
  // const std::string csv_path = output_dir + "/" + shape + "_Nx" +
  //                              std::to_string(Nx_fine / factor) +
  //                              "_curvedness_data.csv";
  // std::ofstream csv(csv_path);
  // if (!csv) {
  //   std::fprintf(stderr, "ERROR: cannot open CSV for writing: %s\n",
  //                csv_path.c_str());
  // } else {
  //   csv << "i,j,k,C_parab,C_surface,"
  //          "centroid_x,centroid_y,centroid_z,"
  //          "proj_x,proj_y,proj_z\n";
  //   csv << std::scientific << std::setprecision(10);
  // }
  // -------------------------------------------------------

  // -- coarsen --
  BasicMesh mesh_fine(Nx_fine, Nx_fine, Nx_fine, 1);
  SurfaceVariant surf = makeSurface(shape, mesh_fine);
  const std::string bin_path = binary_filename(output_dir, shape, Nx_fine);

  const int Nx_coarse = Nx_fine / factor;
  BasicMesh mesh_coarse(Nx_coarse, Nx_coarse, Nx_coarse, 1);
  mesh_coarse.setCellBoundaries(
      IRL::Pt(mesh_fine.x(mesh_fine.imin()), mesh_fine.y(mesh_fine.jmin()),
              mesh_fine.z(mesh_fine.kmin())),
      IRL::Pt(mesh_fine.x(mesh_fine.imax() + 1),
              mesh_fine.y(mesh_fine.jmax() + 1),
              mesh_fine.z(mesh_fine.kmax() + 1)));

  using VM = IRL::GeneralMoments3D<VM_ORDER>;
  using SM = IRL::GeneralSurfaceMoments3D<SM_ORDER>;
  Data<std::pair<VM, SM>> moments_coarse(&mesh_coarse);

  coarsenMomentsFromBinary<VM_ORDER, SM_ORDER>(bin_path, factor,
                                               &moments_coarse);

  // -- reconstruction --
  Data<double> velU(&mesh_coarse), velV(&mesh_coarse), velW(&mesh_coarse);
  Data<IRL::VolumeMoments> liq_moments(&mesh_coarse), gas_moments(&mesh_coarse);
  Data<IRL::SeparatorVariant> interface(&mesh_coarse);

  for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
    for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
      for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
        IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j), mesh_coarse.z(k));
        IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
                   mesh_coarse.z(k + 1));
        IRL::RectangularCuboid cell =
            IRL::RectangularCuboid::fromBoundingPts(x0, x1);

        const double m0 = moments_coarse(i, j, k).first[0];

        const IRL::Pt m1(moments_coarse(i, j, k).first[1],
                         moments_coarse(i, j, k).first[2],
                         moments_coarse(i, j, k).first[3]);

        liq_moments(i, j, k) = IRL::VolumeMoments(m0, m1);
        gas_moments(i, j, k) = cell.calculateMoments() - liq_moments(i, j, k);
      }
    }
  }

  getReconstruction(reconstruction_method, liq_moments, gas_moments, 0.0, velU,
                    velV, velW, &interface);

  // -- curvedness --
  double L2_sum = 0.0;
  int N = 0;

  for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
    for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
      for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
        double vf = liq_moments(i, j, k).volume() / mesh_coarse.cell_volume();
        if (vf < IRL::global_constants::VF_LOW ||
            vf > IRL::global_constants::VF_HIGH)
          continue;

        const auto* parab = std::get_if<IRL::Paraboloid>(&interface(i, j, k));
        if (!parab) continue;

        IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j), mesh_coarse.z(k));
        IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
                   mesh_coarse.z(k + 1));
        IRL::RectangularCuboid cell =
            IRL::RectangularCuboid::fromBoundingPts(x0, x1);

        using VMS =
            IRL::AddSurfaceOutput<IRL::VolumeMoments,
                                  IRL::ParaboloidParametrizedSurfaceOutput>;
        auto surface = IRL::getVolumeMoments<VMS>(cell, *parab).getSurface();

        auto surface_moments = surface.template getSurfaceMoments<1>();
        Eigen::Vector3d surface_centroid(
            surface_moments[1] / surface_moments[0],
            surface_moments[2] / surface_moments[0],
            surface_moments[3] / surface_moments[0]);

        Eigen::Vector3d x_proj = std::visit(
            [&](auto& s) { return s.projectPointOnSurface(surface_centroid); },
            surf);

        double surface_curvedness = std::visit(
            [&](auto& s) {
              return s.curvedness(x_proj[0], x_proj[1], x_proj[2]);
            },
            surf);

        double parab_curvedness = surface.getCurvednessNonAligned(
            IRL::Pt(x_proj[0], x_proj[1], x_proj[2]));

        norms.first = std::max(norms.first,
                               std::abs(surface_curvedness - parab_curvedness));
        L2_sum += (surface_curvedness - parab_curvedness) *
                  (surface_curvedness - parab_curvedness);
        N++;

        // ------------------------ csv output -----------------------
        // if (csv) {
        //   csv << i << ',' << j << ',' << k << ',' << parab_curvedness <<
        //   ','
        //       << surface_curvedness << ',' << surface_centroid[0] << ','
        //       << surface_centroid[1] << ',' << surface_centroid[2] << ','
        //       << x_proj[0] << ',' << x_proj[1] << ',' << x_proj[2] << '\n';
        // }
        // -----------------------------------------------------------
      }
    }
  }
  // csv.close();
  // if (csv) {
  //   std::printf("curvedness data outputted to csv: %s\n",
  //   csv_path.c_str());
  // }

  norms.second = std::sqrt(1.0 / static_cast<double>(N) * L2_sum);
  // }

  // double buf[2];
  // if (rank == 0) {
  //   buf[0] = norms.first;
  //   buf[1] = norms.second;
  // }
  // MPI_Bcast(buf, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // if (rank != 0) {
  //   norms.first = buf[0];
  //   norms.second = buf[1];
  // }

  return norms;
}

template <std::size_t VM_ORDER, std::size_t SM_ORDER>
void runCurvednessConvergence(const std::string& shape, int Nx_fine,
                              const std::string& reconstruction_method,
                              const std::string& output_dir) {
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::string csv_path = output_dir + "/curvedness_convergence_" + shape + "_" +
                         reconstruction_method + "_Nx" +
                         std::to_string(Nx_fine) + ".csv";

  // factors for convergence
  const std::vector<int> factors = {1, 2, 4, 8, 16};  // {1, 2, 4, 8, 16}

  if (rank == 0) {
    std::ofstream csv(csv_path);
    if (!csv) {
      std::fprintf(stderr, "ERROR: cannot open CSV for writing: %s\n",
                   csv_path.c_str());
    } else {
      csv << "factor,Linf,L2\n";
      csv << std::scientific << std::setprecision(10);
    }
    for (int f : factors) {
      std::cout << "Running factor = " << f << std::endl;
      int ok = (Nx_fine % f == 0) ? 1 : 0;
      if (!ok) {
        std::fprintf(stderr,
                     "ERROR: Nx_fine (%d) not divisible by factor (%d)\n",
                     Nx_fine, f);
        if (!ok) MPI_Abort(MPI_COMM_WORLD, 4);
      }
      std::pair<double, double> norms =
          getCurvednessMetrics<VM_ORDER, SM_ORDER>(
              shape, Nx_fine, f, reconstruction_method, output_dir);
      if (csv) {
        csv << f << "," << norms.first << "," << norms.second << "\n";
      }
    }
    csv.close();
    std::printf("✅ Curvedness convergence results written to %s\n",
                csv_path.c_str());
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

// --------------------------- LVIRA norms --------------------------------
template <std::size_t VM_ORDER, std::size_t SM_ORDER>
MomentDiffNorms computePLICMetricsFromBin(
    const int& factor, const std::string& reconstruction_method,
    const std::string& shape, int Nx_fine, const std::string& output_dir) {
  if (reconstruction_method != "LVIRA" && reconstruction_method != "ELVIRA") {
    std::cerr << "ERROR: reconstruction method must be LVIRA or ELVIRA (got "
              << reconstruction_method << ")\n";
    return {};
  }

  // -- coarsen --
  BasicMesh mesh_fine(Nx_fine, Nx_fine, Nx_fine, 1);
  SurfaceVariant surf = makeSurface(shape, mesh_fine);
  const std::string bin_path = binary_filename(output_dir, shape, Nx_fine);

  const int Nx_coarse = Nx_fine / factor;
  BasicMesh mesh_coarse(Nx_coarse, Nx_coarse, Nx_coarse, 1);
  mesh_coarse.setCellBoundaries(
      IRL::Pt(mesh_fine.x(mesh_fine.imin()), mesh_fine.y(mesh_fine.jmin()),
              mesh_fine.z(mesh_fine.kmin())),
      IRL::Pt(mesh_fine.x(mesh_fine.imax() + 1),
              mesh_fine.y(mesh_fine.jmax() + 1),
              mesh_fine.z(mesh_fine.kmax() + 1)));

  using VM = IRL::GeneralMoments3D<VM_ORDER>;
  using SM = IRL::GeneralSurfaceMoments3D<SM_ORDER>;
  Data<std::pair<VM, SM>> moments_coarse(&mesh_coarse);

  coarsenMomentsFromBinary<VM_ORDER, SM_ORDER>(bin_path, factor,
                                               &moments_coarse);

  // -- reconstruction --
  Data<double> velU(&mesh_coarse), velV(&mesh_coarse), velW(&mesh_coarse);
  Data<IRL::VolumeMoments> liq_moments(&mesh_coarse), gas_moments(&mesh_coarse);
  Data<IRL::SeparatorVariant> interface(&mesh_coarse);

  for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
    for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
      for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
        IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j), mesh_coarse.z(k));
        IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
                   mesh_coarse.z(k + 1));
        IRL::RectangularCuboid cell =
            IRL::RectangularCuboid::fromBoundingPts(x0, x1);

        const double m0 = moments_coarse(i, j, k).first[0];

        const IRL::Pt m1(moments_coarse(i, j, k).first[1],
                         moments_coarse(i, j, k).first[2],
                         moments_coarse(i, j, k).first[3]);

        liq_moments(i, j, k) = IRL::VolumeMoments(m0, m1);
        gas_moments(i, j, k) = cell.calculateMoments() - liq_moments(i, j, k);
      }
    }
  }

  getReconstruction(reconstruction_method, liq_moments, gas_moments, 0.0, velU,
                    velV, velW, &interface);

  // -- computing PLIC moments --
  Data<std::pair<VM, SM>> moments_reconstruction(&mesh_coarse);
  for (int i = mesh_coarse.imin(); i <= mesh_coarse.imax(); i++) {
    for (int j = mesh_coarse.jmin(); j <= mesh_coarse.jmax(); j++) {
      for (int k = mesh_coarse.kmin(); k <= mesh_coarse.kmax(); k++) {
        double vf = liq_moments(i, j, k).volume() / mesh_coarse.cell_volume();
        const IRL::Pt x0(mesh_coarse.x(i), mesh_coarse.y(j), mesh_coarse.z(k));
        const IRL::Pt x1(mesh_coarse.x(i + 1), mesh_coarse.y(j + 1),
                         mesh_coarse.z(k + 1));
        const IRL::RectangularCuboid cell =
            IRL::RectangularCuboid::fromBoundingPts(x0, x1);
        if (vf > IRL::global_constants::VF_HIGH) {
          moments_reconstruction(i, j, k).first =
              IRL::getVolumeMoments<VM>(cell);
        }
        if (vf > IRL::global_constants::VF_LOW &&
            vf < IRL::global_constants::VF_HIGH) {
          const auto planar_separator =
              std::get<IRL::PlanarSeparator>(interface(i, j, k));
          IRL::Polygon polygon =
              IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                  cell, planar_separator, planar_separator[0]);
          // polygon.calculateAndSetPlaneOfExistence();
          moments_reconstruction(i, j, k).second =
              polygon.calculateGeneralMoments<SM_ORDER>();
          moments_reconstruction(i, j, k).first =
              IRL::getVolumeMoments<VM>(cell, planar_separator);
          // RECENTER MOMENTS
          // CHANGE MOMENT CALC
        }
      }
    }
  }

  // -- moment norms --
  return compute_moment_diff_norms(moments_coarse, moments_reconstruction);
}

#endif  // EXAMPLES_IMPLICIT_SURFACE_RECONSTRUCTION_RECONSTRUCTION_METRICS_TPP_
