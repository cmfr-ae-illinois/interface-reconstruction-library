// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_PU_H_
#define EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_PU_H_

#include <cmath>
#include <limits>

#include "irl/amrex/sepunion_multifab.h"
#include "irl/generic_cutting/cut_polygon.h"

#include "examples/amrex_advector/reconstruction_lvira.h"
#include "examples/amrex_advector/reconstruction_plicnet.h"

using namespace amrex;

struct PU {
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

    return IRL::cutPlaneByHexahedron<IRL::Polygon>(cell, plic[0]);
  }

  static void GetReconstruction(
      SepUnionMultiFab& interface, SepUnionMultiFab& interface_with_ghost,
      const MultiFab& moments, const Geometry& geom,
      std::vector<InterfaceScalarField>* scalar_fields = nullptr) {
    // initial plic reconstruction
    LVIRA::GetReconstruction(interface, interface_with_ghost, moments, geom,
                             nullptr);

    // ---------------------------------------------------------------------
    // params
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    const Real cell_vol = dx[0] * dx[1] * dx[2];
    const Real vol_inv = 1.0 / cell_vol;

    const Real dx_avg = std::cbrt(dx[0] * dx[1] * dx[2]);

    constexpr int nlayers = 1;
    constexpr int nstencil =
        (1 + 2 * nlayers) * (1 + 2 * nlayers) * (1 + 2 * nlayers);

    constexpr double normal_scatter_threshold = 0.10;
    constexpr double volume_error_threshold = 0.05;

    // ---------------------------------------------------------------------
    // some fields

    const BoxArray& ba = interface_with_ghost.boxArray();
    const DistributionMapping& dm = interface_with_ghost.DistributionMap();

    SepUnionMultiFab jibben_interface(ba, dm, interface_with_ghost.nComp(),
                                      interface_with_ghost.nGrowVect());

    SepUnionMultiFab pu_interface(ba, dm, interface_with_ghost.nComp(),
                                  interface_with_ghost.nGrowVect());

    SepUnionMultiFab pu_neighborhood_interface(
        ba, dm, interface_with_ghost.nComp(), interface_with_ghost.nGrowVect());

    // SepUnionMultiFab final_interface(ba, dm, interface.nComp(),
    //                                  interface.nGrowVect());

    iMultiFab is_underresolved(ba, dm, 1, 0);
    is_underresolved.setVal(0);

    // ---------------------------------------------------------------------
    // all interface fields are initialized with plic

    jibben_interface.LocalCopy(interface_with_ghost, 0, 0,
                               interface_with_ghost.nComp(),
                               interface_with_ghost.nGrowVect());

    pu_interface.LocalCopy(interface_with_ghost, 0, 0,
                           interface_with_ghost.nComp(),
                           interface_with_ghost.nGrowVect());

    pu_neighborhood_interface.LocalCopy(interface_with_ghost, 0, 0,
                                        interface_with_ghost.nComp(),
                                        interface_with_ghost.nGrowVect());

    // final_interface.LocalCopy(interface, 0, 0, interface.nComp(),
    //                           interface.nGrowVect());

    jibben_interface.FillBoundaryWithPeriodicShift(geom);
    pu_interface.FillBoundaryWithPeriodicShift(geom);
    pu_neighborhood_interface.FillBoundaryWithPeriodicShift(geom);
    // final_interface.FillBoundaryWithPeriodicShift(geom);

    // ---------------------------------------------------------------------
    // Jibben reconstruction

    for (MFIter mfi(interface_with_ghost, TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      Array4<IRL::SeparatorUnion> jibben_array = jibben_interface.array(mfi);

      Array4<IRL::SeparatorUnion> pu_neigh_array =
          pu_neighborhood_interface.array(mfi);

      Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
          interface_with_ghost.const_array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);

      Array4<int> underresolved_array = is_underresolved.array(mfi);

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const double liq_vf = moments_array(i, j, k) * vol_inv;

        if (liq_vf < IRL::global_constants::VF_LOW ||
            liq_vf > IRL::global_constants::VF_HIGH ||
            interface_with_ghost_array(i, j, k).type() !=
                IRL::SeparatorUnion::SeparatorType::OnePlane) {
          return;
        }

        const auto center_planar_separator = IRL::PlanarSeparator::fromOnePlane(
            interface_with_ghost_array(i, j, k).getPlane());
        const auto center_polygon =
            polygonFromPLIC(i, j, k, center_planar_separator, problo, dx);
        const double center_area = center_polygon.calculateVolume();

        if (center_polygon.getNumberOfVertices() <= 2 ||
            !std::isfinite(center_area) ||
            center_area <= std::numeric_limits<double>::epsilon()) {
          underresolved_array(i, j, k) = 1;
          return;
        }

        IRL::JibbenNeighborhood jibben_neighborhood;
        jibben_neighborhood.reserve(nstencil);
        jibben_neighborhood.setDelta(2.5 * dx_avg);
        jibben_neighborhood.emptyNeighborhood();

        int jibben_count = 0;
        double vf_supercell = 0.0;

        // Build Jibben stencil
        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              const double liq_vf_neighbor =
                  moments_array(ii, jj, kk) * vol_inv;

              if (liq_vf_neighbor < IRL::global_constants::VF_LOW ||
                  liq_vf_neighbor > IRL::global_constants::VF_HIGH) {
                continue;
              }

              const auto planar_separator = IRL::PlanarSeparator::fromOnePlane(
                  interface_with_ghost_array(ii, jj, kk).getPlane());

              const auto polygon =
                  polygonFromPLIC(ii, jj, kk, planar_separator, problo, dx);

              const double polygon_area = polygon.calculateVolume();

              if (polygon.getNumberOfVertices() > 2 &&
                  std::isfinite(polygon_area) &&
                  polygon_area > std::numeric_limits<double>::epsilon()) {
                jibben_neighborhood.addMember(polygon);

                if (i == ii && j == jj && k == kk) {
                  jibben_neighborhood.setCenterOfStencil(jibben_count);
                }

                ++jibben_count;
                vf_supercell += liq_vf_neighbor;
              }
            }
          }
        }

        // Neighbor count check
        if (jibben_count < 2) {
          underresolved_array(i, j, k) = 1;
          return;
        }

        jibben_neighborhood.localize();

        IRL::Jibben_3D jibben(&jibben_neighborhood);

        // Jibben reconstruction
        IRL::Paraboloid jibben_paraboloid = jibben.solve2(&jibben_neighborhood);

        // Normal scatter check
        const double normal_scatter = jibben.getNormalEigenMetric();

        if (normal_scatter > normal_scatter_threshold) {
          underresolved_array(i, j, k) = 1;
        }

        // Squared volume error check
        const double volume_error = jibben.getVolumeErrorSquared(dx_avg);

        if (volume_error > volume_error_threshold) {
          underresolved_array(i, j, k) = 1;
        }

        // Supercell check
        if (liq_vf < 1.0e-2 && vf_supercell < 1.0e-1) {
          underresolved_array(i, j, k) = 1;
        }

        jibben_array(i, j, k) = jibben_paraboloid;

        // Match volume fraction
        const auto cell = makeCell(i, j, k, problo, dx);

        IRL::setDistanceToMatchVolumeFraction(cell, liq_vf,
                                              &jibben_array(i, j, k), 1.0e-14);

        // Update PU-neighborhood interface with acceptable paraboloids
        const auto& aligned_paraboloid =
            jibben_paraboloid.getAlignedParaboloid();

        if (std::fabs(aligned_paraboloid.a()) * dx_avg <= 1.0 &&
            std::fabs(aligned_paraboloid.b()) * dx_avg <= 1.0 &&
            underresolved_array(i, j, k) == 0) {
          pu_neigh_array(i, j, k) = jibben_array(i, j, k);
        }
      });
    }

    jibben_interface.FillBoundaryWithPeriodicShift(geom);
    pu_neighborhood_interface.FillBoundaryWithPeriodicShift(geom);

    // ---------------------------------------------------------------------
    // PU reconstruction

    for (MFIter mfi(interface_with_ghost, TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      Array4<IRL::SeparatorUnion> pu_array = pu_interface.array(mfi);

      Array4<IRL::SeparatorUnion const> pu_neigh_array =
          pu_neighborhood_interface.const_array(mfi);

      Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
          interface_with_ghost.const_array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const double liq_vf = moments_array(i, j, k) * vol_inv;

        if (liq_vf < IRL::global_constants::VF_LOW ||
            liq_vf > IRL::global_constants::VF_HIGH) {
          return;
        }

        const auto center_planar_separator = IRL::PlanarSeparator::fromOnePlane(
            interface_with_ghost_array(i, j, k).getPlane());
        const auto center_polygon =
            polygonFromPLIC(i, j, k, center_planar_separator, problo, dx);
        const double center_area = center_polygon.calculateVolume();

        if (center_polygon.getNumberOfVertices() <= 2 ||
            !std::isfinite(center_area) ||
            center_area <= std::numeric_limits<double>::epsilon()) {
          return;
        }

        IRL::PUNeighborhood pu_neighborhood;
        pu_neighborhood.reserve(nstencil);
        pu_neighborhood.emptyNeighborhood();

        int pu_count = 0;

        // Build PU neighborhood
        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              const double liq_vf_neighbor =
                  moments_array(ii, jj, kk) * vol_inv;

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
              const IRL::Pt centroid = polygon.calculateCentroid();

              const IRL::SeparatorVariant neighbor_interface =
                  pu_neigh_array(ii, jj, kk);

              pu_neighborhood.addMember(neighbor_interface, centroid, weight);

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

        const double delta = 2.5 * dx_avg;

        IRL::Paraboloid pu_paraboloid =
            IRL::reconstructionWithPU3D(pu_neighborhood, delta, dx_avg);

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

        pu_array(i, j, k) = pu_paraboloid;
      });
    }

    pu_interface.FillBoundaryWithPeriodicShift(geom);

    // ---------------------------------------------------------------------
    // pu clean up

    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> pu_array = pu_interface.array(mfi);

      Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
          interface_with_ghost.const_array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const double liq_vf = moments_array(i, j, k) * vol_inv;

        if (liq_vf < IRL::global_constants::VF_LOW ||
            liq_vf > IRL::global_constants::VF_HIGH) {
          return;
        }

        double vf_supercell = 0.0;

        for (int kk = k - 1; kk <= k + 1; ++kk) {
          for (int jj = j - 1; jj <= j + 1; ++jj) {
            for (int ii = i - 1; ii <= i + 1; ++ii) {
              vf_supercell += moments_array(ii, jj, kk) * vol_inv;
            }
          }
        }

        if (liq_vf < 1.0e-2 && vf_supercell < 1.0e-1) {
          pu_array(i, j, k) = interface_with_ghost_array(i, j, k);
        }
      });
    }

    pu_interface.FillBoundaryWithPeriodicShift(geom);

    // ---------------------------------------------------------------------
    // storing final interface

    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);

      Array4<IRL::SeparatorUnion const> pu_array =
          pu_interface.const_array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const double liq_vf = moments_array(i, j, k) * vol_inv;

        if (liq_vf < IRL::global_constants::VF_LOW ||
            liq_vf > IRL::global_constants::VF_HIGH) {
          return;
        }

        interface_array(i, j, k) = pu_array(i, j, k);
      });
    }

    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());

    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_PU_H_
