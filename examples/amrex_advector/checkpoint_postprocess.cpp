// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

#include "examples/amrex_advector/amrcore_advector.h"

namespace {

constexpr int NUM_MOMENTS = 10;

using MomentArray = std::array<double, NUM_MOMENTS>;

// Moment ordering:
// 0 : M0
// 1 : Mx
// 2 : My
// 3 : Mz
// 4 : Mxx
// 5 : Mxy
// 6 : Mxz
// 7 : Myy
// 8 : Myz
// 9 : Mzz

int momentOrder(const int m) {
  if (m == 0) {
    return 0;
  }
  if (m <= 3) {
    return 1;
  }
  return 2;
}

double momentNormalization(const int m, const double dx) {
  const int order = momentOrder(m);
  if (order == 0) {
    return 1.0;
  }
  if (order == 1) {
    return 1.0 / dx;
  }
  return 1.0 / (dx * dx);
}

void validateMomentOrder(const int order, const std::string& name) {
  if (order < 0 || order > 2) {
    std::ostringstream oss;
    oss << name << " must be 0, 1, or 2; got " << order;
    amrex::Abort(oss.str());
  }
}

void copyVolumeMomentsToArray(const IRL::VolumeMoments& moments,
                              MomentArray& result) {
  result[0] = moments.volume();
  result[1] = moments.centroid()[0];
  result[2] = moments.centroid()[1];
  result[3] = moments.centroid()[2];
}

template <class MomentsType>
void copyGeneralMomentsToArray(const MomentsType& moments,
                               MomentArray& result) {
  for (IRL::UnsignedIndex_t m = 0; m < moments.size(); ++m) {
    result[m] = moments[m];
  }
}

// Recenter moments around cell center xc
MomentArray recenterMoments(const MomentArray& moments, const IRL::Pt& xc) {
  MomentArray centered{};

  const double x = xc[0];
  const double y = xc[1];
  const double z = xc[2];

  const double M0 = moments[0];

  const double Mx = moments[1];
  const double My = moments[2];
  const double Mz = moments[3];

  const double Mxx = moments[4];
  const double Mxy = moments[5];
  const double Mxz = moments[6];
  const double Myy = moments[7];
  const double Myz = moments[8];
  const double Mzz = moments[9];

  centered[0] = M0;
  centered[1] = Mx - M0 * x;
  centered[2] = My - M0 * y;
  centered[3] = Mz - M0 * z;
  centered[4] = Mxx - 2.0 * x * Mx + M0 * x * x;
  centered[5] = Mxy - x * My - y * Mx + M0 * x * y;
  centered[6] = Mxz - x * Mz - z * Mx + M0 * x * z;
  centered[7] = Myy - 2.0 * y * My + M0 * y * y;
  centered[8] = Myz - y * Mz - z * My + M0 * y * z;
  centered[9] = Mzz - 2.0 * z * Mz + M0 * z * z;

  return centered;
}

// Surface moments of a planar polygon through the requested order.
MomentArray calculatePolygonSurfaceMoments(const IRL::Polygon& polygon,
                                           const int max_order) {
  MomentArray moments{};

  const int n = polygon.getNumberOfVertices();

  if (n < 3) {
    return moments;
  }

  const IRL::Pt& v0 = polygon[0];

  for (int iv = 1; iv < n - 1; ++iv) {
    const IRL::Pt& v1 = polygon[iv];
    const IRL::Pt& v2 = polygon[iv + 1];

    IRL::Polygon triangle;

    triangle.addVertex(v0);
    triangle.addVertex(v1);
    triangle.addVertex(v2);

    triangle.calculateAndSetPlaneOfExistence();

    const double area = triangle.calculateVolume();

    // M0
    moments[0] += area;

    if (max_order == 0) {
      continue;
    }

    // M1
    const IRL::Pt centroid = (v0 + v1 + v2) * (1.0 / 3.0);

    moments[1] += centroid[0] * area;
    moments[2] += centroid[1] * area;
    moments[3] += centroid[2] * area;

    if (max_order == 1) {
      continue;
    }

    // M2
    const auto a = v1 - v0;
    const auto b = v2 - v0;

    const double factor = 2.0 * area;

    auto accumulate = [&](const IRL::Pt& u, const IRL::Pt& v,
                          const double scale) {
      moments[4] += scale * u[0] * v[0];
      moments[5] += scale * u[0] * v[1];
      moments[6] += scale * u[0] * v[2];
      moments[7] += scale * u[1] * v[1];
      moments[8] += scale * u[1] * v[2];
      moments[9] += scale * u[2] * v[2];
    };

    accumulate(v0, v0, factor * 0.5);
    accumulate(v0, a, factor / 6.0);
    accumulate(a, v0, factor / 6.0);
    accumulate(v0, b, factor / 6.0);
    accumulate(b, v0, factor / 6.0);
    accumulate(a, a, factor / 12.0);
    accumulate(b, b, factor / 12.0);
    accumulate(a, b, factor / 24.0);
    accumulate(b, a, factor / 24.0);
  }

  return moments;
}

// Volume moments through the requested order.
MomentArray getVolumeMoments(const IRL::RectangularCuboid& cell,
                             const IRL::SeparatorUnion& interface,
                             const double vf, const int max_order) {
  MomentArray result{};

  if (vf <= IRL::global_constants::VF_LOW) {
    return result;
  }

  if (vf >= IRL::global_constants::VF_HIGH) {
    if (max_order == 0) {
      result[0] = IRL::getVolumeMoments<IRL::Volume>(cell);
    } else if (max_order == 1) {
      copyVolumeMomentsToArray(IRL::getVolumeMoments<IRL::VolumeMoments>(cell),
                               result);
    } else {
      copyGeneralMomentsToArray(
          IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell), result);
    }
    return result;
  }

  if (interface.type() == IRL::SeparatorUnion::SeparatorType::OnePlane) {
    const auto separator =
        IRL::PlanarSeparator::fromOnePlane(interface.getPlane());
    if (max_order == 0) {
      result[0] = IRL::getVolumeMoments<IRL::Volume>(cell, separator);
    } else if (max_order == 1) {
      copyVolumeMomentsToArray(
          IRL::getVolumeMoments<IRL::VolumeMoments>(cell, separator), result);
    } else {
      copyGeneralMomentsToArray(
          IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, separator),
          result);
    }
    return result;
  }

  if (interface.type() == IRL::SeparatorUnion::SeparatorType::Paraboloid) {
    const auto paraboloid = interface.getParaboloid();
    if (max_order == 0) {
      result[0] = IRL::getVolumeMoments<IRL::Volume>(cell, paraboloid);
    } else if (max_order == 1) {
      copyVolumeMomentsToArray(
          IRL::getVolumeMoments<IRL::VolumeMoments>(cell, paraboloid), result);
    } else {
      copyGeneralMomentsToArray(
          IRL::getVolumeMoments<IRL::GeneralMoments3D<2>>(cell, paraboloid),
          result);
    }
    return result;
  }

  amrex::Abort("Mixed cell has unsupported interface type.");

  return result;
}

// Surface moments through the requested order.
MomentArray getSurfaceMoments(const IRL::RectangularCuboid& cell,
                              const IRL::SeparatorUnion& interface,
                              const double vf, const int max_order) {
  MomentArray result{};

  if (vf <= IRL::global_constants::VF_LOW ||
      vf >= IRL::global_constants::VF_HIGH) {
    return result;
  }

  if (interface.type() == IRL::SeparatorUnion::SeparatorType::OnePlane) {
    const auto separator =
        IRL::PlanarSeparator::fromOnePlane(interface.getPlane());
    IRL::Polygon polygon = IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
        cell, separator, separator[0]);
    if (polygon.getNumberOfVertices() > 2) {
      polygon.calculateAndSetPlaneOfExistence();
      result = calculatePolygonSurfaceMoments(polygon, max_order);
    }
    return result;
  }

  if (interface.type() == IRL::SeparatorUnion::SeparatorType::Paraboloid) {
    using VolumeMomentsAndSurface =
        IRL::AddSurfaceOutput<IRL::Volume,
                              IRL::ParaboloidParametrizedSurfaceOutput>;
    const auto paraboloid = interface.getParaboloid();
    auto volume_and_surface =
        IRL::getVolumeMoments<VolumeMomentsAndSurface>(cell, paraboloid);
    auto surface = volume_and_surface.getSurface();
    if (max_order == 0) {
      copyGeneralMomentsToArray(surface.getSurfaceMoments<0>(), result);
    } else if (max_order == 1) {
      copyGeneralMomentsToArray(surface.getSurfaceMoments<1>(), result);
    } else {
      copyGeneralMomentsToArray(surface.getSurfaceMoments<2>(), result);
    }
    return result;
  }

  amrex::Abort("Mixed cell has unsupported surface interface type.");

  return result;
}

// Error norms
void computeMomentErrors(const amrex::MultiFab& initial_moments,
                         const amrex::MultiFab& final_moments,
                         const amrex::SepUnionMultiFab& initial_interface,
                         const amrex::SepUnionMultiFab& final_interface,
                         const amrex::Geometry& geom,
                         const int volume_moment_order,
                         const int surface_moment_order) {
  using namespace amrex;

  if (initial_moments.boxArray() != final_moments.boxArray()) {
    amrex::Abort("Initial and final uniform BoxArrays differ.");
  }

  if (initial_moments.DistributionMap() != final_moments.DistributionMap()) {
    amrex::Abort("Initial and final uniform DistributionMappings differ.");
  }

  const auto dx = geom.CellSizeArray();
  const auto problo = geom.ProbLoArray();

  const double cell_volume = dx[0] * dx[1] * dx[2];
  const double h = dx[0];

  const bool compute_volume_M0 = volume_moment_order >= 0;
  const bool compute_volume_M1 = volume_moment_order >= 1;
  const bool compute_volume_M2 = volume_moment_order >= 2;

  const bool compute_surface_M0 = surface_moment_order >= 0;
  const bool compute_surface_M1 = surface_moment_order >= 1;
  const bool compute_surface_M2 = surface_moment_order >= 2;

  // Volume errors
  double volume_L1_M0 = 0.0;
  double volume_L1_M1 = 0.0;
  double volume_L1_M2 = 0.0;
  double volume_Linf_M0 = 0.0;
  double volume_Linf_M1 = 0.0;
  double volume_Linf_M2 = 0.0;

  // Surface errors
  double surface_L1_M0 = 0.0;
  double surface_L1_M1 = 0.0;
  double surface_L1_M2 = 0.0;
  double surface_Linf_M0 = 0.0;
  double surface_Linf_M1 = 0.0;
  double surface_Linf_M2 = 0.0;

  for (MFIter mfi(initial_moments); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();

    const auto initial_mom_arr = initial_moments.const_array(mfi);
    const auto final_mom_arr = final_moments.const_array(mfi);

    const auto initial_interface_arr = initial_interface.const_array(mfi);
    const auto final_interface_arr = final_interface.const_array(mfi);

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          const double x0 = problo[0] + i * dx[0];
          const double y0 = problo[1] + j * dx[1];
          const double z0 = problo[2] + k * dx[2];

          const IRL::Pt lower(x0, y0, z0);
          const IRL::Pt upper(x0 + dx[0], y0 + dx[1], z0 + dx[2]);

          const IRL::Pt center(x0 + 0.5 * dx[0], y0 + 0.5 * dx[1],
                               z0 + 0.5 * dx[2]);

          const auto cell =
              IRL::RectangularCuboid::fromBoundingPts(lower, upper);

          const double initial_vf = initial_mom_arr(i, j, k, comp_vf);
          const double final_vf = final_mom_arr(i, j, k, comp_vf);

          MomentArray initial_volume =
              getVolumeMoments(cell, initial_interface_arr(i, j, k), initial_vf,
                               volume_moment_order);
          MomentArray final_volume =
              getVolumeMoments(cell, final_interface_arr(i, j, k), final_vf,
                               volume_moment_order);

          // Recenter about cell center
          initial_volume = recenterMoments(initial_volume, center);
          final_volume = recenterMoments(final_volume, center);

          if (compute_volume_M0) {
            const double M0_error =
                std::abs(final_volume[0] - initial_volume[0]);
            volume_L1_M0 += M0_error;
            volume_Linf_M0 = std::max(volume_Linf_M0, M0_error);
          }

          if (compute_volume_M1) {
            const double dMx = final_volume[1] - initial_volume[1];
            const double dMy = final_volume[2] - initial_volume[2];
            const double dMz = final_volume[3] - initial_volume[3];
            const double M1_error =
                std::sqrt(dMx * dMx + dMy * dMy + dMz * dMz);
            volume_L1_M1 += M1_error;
            volume_Linf_M1 = std::max(volume_Linf_M1, M1_error);
          }

          if (compute_volume_M2) {
            const double dMxx = final_volume[4] - initial_volume[4];
            const double dMxy = final_volume[5] - initial_volume[5];
            const double dMxz = final_volume[6] - initial_volume[6];
            const double dMyy = final_volume[7] - initial_volume[7];
            const double dMyz = final_volume[8] - initial_volume[8];
            const double dMzz = final_volume[9] - initial_volume[9];
            const double M2_error =
                std::sqrt(dMxx * dMxx + dMyy * dMyy + dMzz * dMzz +
                          2.0 * (dMxy * dMxy + dMxz * dMxz + dMyz * dMyz));
            volume_L1_M2 += M2_error;
            volume_Linf_M2 = std::max(volume_Linf_M2, M2_error);
          }

          MomentArray initial_surface =
              getSurfaceMoments(cell, initial_interface_arr(i, j, k),
                                initial_vf, surface_moment_order);
          MomentArray final_surface =
              getSurfaceMoments(cell, final_interface_arr(i, j, k), final_vf,
                                surface_moment_order);

          // Recenter about cell center
          initial_surface = recenterMoments(initial_surface, center);
          final_surface = recenterMoments(final_surface, center);

          if (compute_surface_M0) {
            const double SM0_error =
                std::abs(final_surface[0] - initial_surface[0]);
            surface_L1_M0 += SM0_error;
            surface_Linf_M0 = std::max(surface_Linf_M0, SM0_error);
          }

          if (compute_surface_M1) {
            const double dSMx = final_surface[1] - initial_surface[1];
            const double dSMy = final_surface[2] - initial_surface[2];
            const double dSMz = final_surface[3] - initial_surface[3];
            const double surface_M1_error =
                std::sqrt(dSMx * dSMx + dSMy * dSMy + dSMz * dSMz);
            surface_L1_M1 += surface_M1_error;
            surface_Linf_M1 = std::max(surface_Linf_M1, surface_M1_error);
          }

          if (compute_surface_M2) {
            const double dSMxx = final_surface[4] - initial_surface[4];
            const double dSMxy = final_surface[5] - initial_surface[5];
            const double dSMxz = final_surface[6] - initial_surface[6];
            const double dSMyy = final_surface[7] - initial_surface[7];
            const double dSMyz = final_surface[8] - initial_surface[8];
            const double dSMzz = final_surface[9] - initial_surface[9];
            const double surface_M2_error = std::sqrt(
                dSMxx * dSMxx + dSMyy * dSMyy + dSMzz * dSMzz +
                2.0 * (dSMxy * dSMxy + dSMxz * dSMxz + dSMyz * dSMyz));
            surface_L1_M2 += surface_M2_error;
            surface_Linf_M2 = std::max(surface_Linf_M2, surface_M2_error);
          }
        }
      }
    }
  }

  // MPI reduction
  if (compute_volume_M0) ParallelDescriptor::ReduceRealSum(volume_L1_M0);
  if (compute_volume_M1) ParallelDescriptor::ReduceRealSum(volume_L1_M1);
  if (compute_volume_M2) ParallelDescriptor::ReduceRealSum(volume_L1_M2);
  if (compute_volume_M0) ParallelDescriptor::ReduceRealMax(volume_Linf_M0);
  if (compute_volume_M1) ParallelDescriptor::ReduceRealMax(volume_Linf_M1);
  if (compute_volume_M2) ParallelDescriptor::ReduceRealMax(volume_Linf_M2);

  if (compute_surface_M0) ParallelDescriptor::ReduceRealSum(surface_L1_M0);
  if (compute_surface_M1) ParallelDescriptor::ReduceRealSum(surface_L1_M1);
  if (compute_surface_M2) ParallelDescriptor::ReduceRealSum(surface_L1_M2);
  if (compute_surface_M0) ParallelDescriptor::ReduceRealMax(surface_Linf_M0);
  if (compute_surface_M1) ParallelDescriptor::ReduceRealMax(surface_Linf_M1);
  if (compute_surface_M2) ParallelDescriptor::ReduceRealMax(surface_Linf_M2);

  // Normalize by moment order
  if (compute_volume_M1) {
    volume_L1_M1 /= h;
    volume_Linf_M1 /= h;
  }
  if (compute_volume_M2) {
    volume_L1_M2 /= (h * h);
    volume_Linf_M2 /= (h * h);
  }

  if (compute_surface_M1) {
    surface_L1_M1 /= h;
    surface_Linf_M1 /= h;
  }
  if (compute_surface_M2) {
    surface_L1_M2 /= (h * h);
    surface_Linf_M2 /= (h * h);
  }

  // Print
  amrex::Print() << "\n"
                 << "=============================================\n"
                 << "VOLUME MOMENT ERRORS\n"
                 << "=============================================\n"
                 << std::scientific << std::setprecision(16);
  if (compute_volume_M0)
    amrex::Print() << "M0  L1 = " << volume_L1_M0
                   << "  Linf = " << volume_Linf_M0 << "\n";
  if (compute_volume_M1)
    amrex::Print() << "M1  L1 = " << volume_L1_M1
                   << "  Linf = " << volume_Linf_M1 << "\n";
  if (compute_volume_M2)
    amrex::Print() << "M2  L1 = " << volume_L1_M2
                   << "  Linf = " << volume_Linf_M2 << "\n";

  amrex::Print() << "\n"
                 << "=============================================\n"
                 << "SURFACE MOMENT ERRORS\n"
                 << "=============================================\n"
                 << std::scientific << std::setprecision(16);
  if (compute_surface_M0)
    amrex::Print() << "M0  L1 = " << surface_L1_M0
                   << "  Linf = " << surface_Linf_M0 << "\n";
  if (compute_surface_M1)
    amrex::Print() << "M1  L1 = " << surface_L1_M1
                   << "  Linf = " << surface_Linf_M1 << "\n";
  if (compute_surface_M2)
    amrex::Print() << "M2  L1 = " << surface_L1_M2
                   << "  Linf = " << surface_Linf_M2 << "\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  amrex::Initialize(argc, argv);

  {
    BL_PROFILE("checkpoint_postprocess()");
    std::string initial_checkpoint;
    std::string final_checkpoint;
    int volume_moment_order = 2;
    int surface_moment_order = 2;
    {
      amrex::ParmParse pp("postprocess");
      pp.get("initial_checkpoint", initial_checkpoint);
      pp.get("final_checkpoint", final_checkpoint);
      pp.query("volumetric_moment_order", volume_moment_order);
      pp.query("surface_moment_order", surface_moment_order);
    }

    validateMomentOrder(volume_moment_order,
                        "postprocess.volumetric_moment_order");
    validateMomentOrder(surface_moment_order,
                        "postprocess.surface_moment_order");

    const auto start_total = amrex::second();

    AmrCoreAdv amr_core_adv;

    // Initial checkpoint.

    amrex::MultiFab initial_uniform_moments;
    amrex::SepUnionMultiFab initial_uniform_interface;
    amr_core_adv.BuildUniformCheckpointState(
        initial_checkpoint, initial_uniform_moments, initial_uniform_interface);

    // Keep the geometry corresponding to this uniform mesh.
    const amrex::Geometry initial_geom = amr_core_adv.GetFinestGeometry();

    // Final checkpoint.

    amrex::MultiFab final_uniform_moments;
    amrex::SepUnionMultiFab final_uniform_interface;
    amr_core_adv.BuildUniformCheckpointState(
        final_checkpoint, final_uniform_moments, final_uniform_interface);

    // Redistribute final state onto the initial state's DistributionMapping
    {
      const auto& common_ba = initial_uniform_moments.boxArray();
      const auto& common_dm = initial_uniform_moments.DistributionMap();

      // Moments
      amrex::MultiFab aligned_final_moments(common_ba, common_dm,
                                            final_uniform_moments.nComp(), 0);
      aligned_final_moments.setVal(0.0);
      aligned_final_moments.ParallelCopy(final_uniform_moments, 0, 0,
                                         final_uniform_moments.nComp(), 0, 0);
      std::swap(final_uniform_moments, aligned_final_moments);
      // Interface
      amrex::SepUnionMultiFab aligned_final_interface(common_ba, common_dm, 1,
                                                      0);
      InitializeSepUnionMultiFab(aligned_final_interface);
      aligned_final_interface.ParallelCopy(final_uniform_interface, 0, 0, 1, 0,
                                           0);
      std::swap(final_uniform_interface, aligned_final_interface);
    }
    const amrex::Geometry final_geom = amr_core_adv.GetFinestGeometry();
    if (initial_geom.Domain() != final_geom.Domain()) {
      amrex::Abort("Initial and final finest geometries differ.");
    }

    // Calculate error norms.

    computeMomentErrors(initial_uniform_moments, final_uniform_moments,
                        initial_uniform_interface, final_uniform_interface,
                        final_geom, volume_moment_order, surface_moment_order);

    auto total_time = amrex::second() - start_total;

    amrex::ParallelDescriptor::ReduceRealMax(
        total_time, amrex::ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "\nCheckpoint postprocess time: " << total_time << "\n";
  }

  amrex::Finalize();

  return 0;
}
