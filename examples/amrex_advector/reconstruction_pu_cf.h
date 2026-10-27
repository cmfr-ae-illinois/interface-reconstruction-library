// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_PU_CF_H_
#define EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_PU_CF_H_

#include <cmath>
#include <limits>

#include "irl/amrex/sepunion_multifab.h"

#include "examples/amrex_advector/reconstruction_lvira.h"
#include "examples/amrex_advector/reconstruction_plicnet.h"
#include "irl/interface_reconstruction_methods/cf.h"

using namespace amrex;

struct PU_CF {
  static IRL::RectangularCuboid makeCell(
      const int i, const int j, const int k,
      const GpuArray<Real, AMREX_SPACEDIM>& problo,
      const GpuArray<Real, AMREX_SPACEDIM>& dx) {
    const double x = problo[0] + i * dx[0];
    const double y = problo[1] + j * dx[1];
    const double z = problo[2] + k * dx[2];

    return IRL::RectangularCuboid::fromBoundingPts(
        IRL::Pt(x, y, z), IRL::Pt(x + dx[0], y + dx[1], z + dx[2]));
  }

  static IRL::Polygon polygonFromPLIC(
      const int i, const int j, const int k, const IRL::PlanarSeparator& plic,
      const GpuArray<Real, AMREX_SPACEDIM>& problo,
      const GpuArray<Real, AMREX_SPACEDIM>& dx) {
    const auto cell = makeCell(i, j, k, problo, dx);

    return IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(cell, plic,
                                                                plic[0]);
  }

  static void GetReconstruction(
      SepUnionMultiFab& interface, SepUnionMultiFab& interface_with_ghost,
      const MultiFab& moments, const Geometry& geom,
      std::vector<InterfaceScalarField>* scalar_fields = nullptr,
      Real* reconstruction_loop_time = nullptr) {
    // plic
    LVIRA::GetReconstruction(interface, interface_with_ghost, moments, geom,
                             nullptr, reconstruction_loop_time);

    // some params
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    const Real cell_vol = dx[0] * dx[1] * dx[2];
    const Real vol_inv = 1.0 / cell_vol;

    const Real dx_avg = std::cbrt(dx[0] * dx[1] * dx[2]);

    constexpr int nlayers = 1;

    constexpr int nstencil =
        (1 + 2 * nlayers) * (1 + 2 * nlayers) * (1 + 2 * nlayers);

    // multifab to store cf interfaces
    const BoxArray& ba = interface_with_ghost.boxArray();
    const DistributionMapping& dm = interface_with_ghost.DistributionMap();
    SepUnionMultiFab cf_interface(ba, dm, interface_with_ghost.nComp(),
                                  interface_with_ghost.nGrowVect());
    cf_interface.LocalCopy(interface_with_ghost, 0, 0,
                           interface_with_ghost.nComp(),
                           interface_with_ghost.nGrowVect());
    cf_interface.FillBoundaryWithPeriodicShift(geom);

    // reconstruction
    ReconstructionLoopTimer loop_timer(reconstruction_loop_time);
    loop_timer.start();
    for (MFIter mfi(cf_interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> cf_interface_array = cf_interface.array(mfi);

      Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
          interface_with_ghost.const_array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const double vf = moments_array(i, j, k, comp_vf);

        if (vf < IRL::global_constants::VF_LOW) {
          cf_interface_array(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
          return;
        } else if (vf > IRL::global_constants::VF_HIGH) {
          cf_interface_array(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
          return;
        }

        const auto center_plic = IRL::PlanarSeparator::fromOnePlane(
            interface_with_ghost_array(i, j, k).getPlane());

        const auto center_polygon =
            polygonFromPLIC(i, j, k, center_plic, problo, dx);
        const double center_area = center_polygon.calculateVolume();

        if (center_polygon.getNumberOfVertices() <= 2 ||
            !std::isfinite(center_area) ||
            center_area <= std::numeric_limits<double>::epsilon()) {
          return;
        }

        IRL::JibbenNeighborhood neighborhood;
        neighborhood.reserve(nstencil);
        neighborhood.setDelta(2.5 * dx_avg);
        neighborhood.emptyNeighborhood();

        int count = 0;

        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              const double neighbor_vf = moments_array(ii, jj, kk, comp_vf);

              if (!std::isfinite(neighbor_vf)) {
                continue;
              }

              if (neighbor_vf < IRL::global_constants::VF_LOW ||
                  neighbor_vf > IRL::global_constants::VF_HIGH) {
                continue;
              }

              const auto neighbor_plic = IRL::PlanarSeparator::fromOnePlane(
                  interface_with_ghost_array(ii, jj, kk).getPlane());

              const auto neighbor_polygon =
                  polygonFromPLIC(ii, jj, kk, neighbor_plic, problo, dx);
              const double neighbor_area = neighbor_polygon.calculateVolume();

              if (neighbor_polygon.getNumberOfVertices() <= 2 ||
                  !std::isfinite(neighbor_area) ||
                  neighbor_area <= std::numeric_limits<double>::epsilon()) {
                continue;
              }

              neighborhood.addMember(neighbor_polygon, neighbor_vf);

              if (i == ii && j == jj && k == kk) {
                neighborhood.setCenterOfStencil(count);
              }

              ++count;
            }
          }
        }

        if (neighborhood.size() < 2) {
          return;
        }

        IRL::CircleFit_3D circle_fit;

        auto paraboloid = circle_fit.solve(&neighborhood, dx_avg);

        if (!std::isfinite(paraboloid.getDatum()[0]) ||
            !std::isfinite(paraboloid.getDatum()[1]) ||
            !std::isfinite(paraboloid.getDatum()[2]) ||
            !std::isfinite(paraboloid.getAlignedParaboloid().a()) ||
            !std::isfinite(paraboloid.getAlignedParaboloid().b())) {
          return;
        }

        cf_interface_array(i, j, k) = paraboloid;

        const auto cell = makeCell(i, j, k, problo, dx);

        IRL::setDistanceToMatchVolumeFraction(
            cell, vf, &cf_interface_array(i, j, k), 1.0e-14);
      });
    }
    loop_timer.stop();
    cf_interface.FillBoundaryWithPeriodicShift(geom);

    // PU fit now
    loop_timer.start();
    for (MFIter mfi(interface_with_ghost, TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      Array4<IRL::SeparatorUnion const> cf_interface_array =
          cf_interface.const_array(mfi);

      Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);

      Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
          interface_with_ghost.const_array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const double liq_vf = moments_array(i, j, k, comp_vf);

        if (liq_vf < IRL::global_constants::VF_LOW ||
            liq_vf > IRL::global_constants::VF_HIGH) {
          return;
        }

        const auto center_planar_separator = IRL::PlanarSeparator::fromOnePlane(
            interface_with_ghost_array(i, j, k).getPlane());
        const auto center_polygon =
            polygonFromPLIC(i, j, k, center_planar_separator, problo, dx);
        const double center_area = center_polygon.calculateVolume();
        const IRL::Normal center_normal = center_planar_separator[0].normal();

        if (center_polygon.getNumberOfVertices() <= 2 ||
            !std::isfinite(center_area) ||
            center_area <= std::numeric_limits<double>::epsilon()) {
          return;
        }

        using CellType = IRL::RectangularCuboid;
        IRL::PUNeighborhood<CellType> pu_neighborhood;
        pu_neighborhood.reserve(nstencil);
        pu_neighborhood.emptyNeighborhood();

        const CellType center_cell = makeCell(i, j, k, problo, dx);
        pu_neighborhood.setCenterCell(&center_cell);

        int pu_count = 0;

        // Build PU neighborhood
        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              const double liq_vf_neighbor = moments_array(ii, jj, kk, comp_vf);

              if (liq_vf_neighbor < IRL::global_constants::VF_LOW ||
                  liq_vf_neighbor > IRL::global_constants::VF_HIGH) {
                continue;
              }

              const auto planar_separator = IRL::PlanarSeparator::fromOnePlane(
                  interface_with_ghost_array(ii, jj, kk).getPlane());

              const auto polygon =
                  polygonFromPLIC(ii, jj, kk, planar_separator, problo, dx);

              if (polygon.getNumberOfVertices() <= 2) {
                continue;
              }

              const double area_weight =
                  polygon.calculateVolume() / (dx_avg * dx_avg);

              if (!std::isfinite(area_weight) ||
                  area_weight <= std::numeric_limits<double>::epsilon()) {
                continue;
              }

              double vfrac_weight = 1.0;

              if (liq_vf_neighbor < 0.1) {
                vfrac_weight =
                    0.5 - 0.5 * std::cos(10.0 * M_PI * liq_vf_neighbor);
              } else if (liq_vf_neighbor > 0.9) {
                vfrac_weight =
                    0.5 - 0.5 * std::cos(10.0 * M_PI * (1.0 - liq_vf_neighbor));
              }

              const double weight = area_weight * vfrac_weight;
              const IRL::Normal neighbor_normal = planar_separator[0].normal();
              const double normal_weight =
                  std::max(0.0, center_normal[0] * neighbor_normal[0] +
                                    center_normal[1] * neighbor_normal[1] +
                                    center_normal[2] * neighbor_normal[2]);
              const IRL::Pt centroid = polygon.calculateCentroid();

              const IRL::SeparatorVariant neighbor_interface =
                  cf_interface_array(ii, jj, kk);

              pu_neighborhood.addMember(&centroid, &neighbor_interface,
                                        weight * normal_weight);

              if (i == ii && j == jj && k == kk) {
                pu_neighborhood.setCenterOfStencil(pu_count);
              }

              ++pu_count;
            }
          }
        }

        // Neighbor count check
        if (pu_count < 2) {
          return;
        }

        const double kernel_size = 2.5 * dx_avg;

        IRL::PUParaboloid<CellType> pu_solver(pu_neighborhood, kernel_size,
                                              dx_avg);

        IRL::Paraboloid pu_paraboloid = pu_solver.solve();

        // Check for NaN/Inf datum
        if (!std::isfinite(pu_paraboloid.getDatum()[0]) ||
            !std::isfinite(pu_paraboloid.getDatum()[1]) ||
            !std::isfinite(pu_paraboloid.getDatum()[2])) {
          return;
        }

        // Revert to PLIC if curvature is too large
        const auto& aligned_paraboloid = pu_paraboloid.getAlignedParaboloid();

        if (std::fabs(aligned_paraboloid.a()) * dx_avg > 4.0 ||
            std::fabs(aligned_paraboloid.b()) * dx_avg > 4.0) {
          return;
        }

        // Volume fraction matching
        const auto new_datum = pu_paraboloid.getDatum();
        const auto new_frame = pu_paraboloid.getReferenceFrame();

        const auto cell = makeCell(i, j, k, problo, dx);

        IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
            solver_distance(cell, liq_vf, 1.0e-14, pu_paraboloid);

        pu_paraboloid.setDatum(
            IRL::Pt(new_datum + solver_distance.getDistance() * new_frame[2]));

        interface_array(i, j, k) = pu_paraboloid;
      });
    }
    loop_timer.stop();

    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());

    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_PU_CF_H_
