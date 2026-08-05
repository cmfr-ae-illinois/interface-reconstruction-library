// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_HYBRID_H_
#define EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_HYBRID_H_

#include <algorithm>
#include <cmath>
#include <limits>

#include "irl/amrex/sepunion_multifab.h"
#include "irl/generic_cutting/cut_polygon.h"

#include "examples/amrex_advector/reconstruction_lvira.h"
#include "examples/amrex_advector/reconstruction_plicnet.h"
#include "examples/amrex_advector/reconstruction_pu.h"
#include "irl/interface_reconstruction_methods/iterative_jibben.h"

using namespace amrex;

struct HYBRID {
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

  static IRL::Normal localToGlobalNormal(
      const IRL::Normal& a_local_normal,
      const IRL::ReferenceFrame& a_local_frame) {
    const IRL::Normal& e_xi = a_local_frame[0];
    const IRL::Normal& e_eta = a_local_frame[1];
    const IRL::Normal& e_zeta = a_local_frame[2];

    IRL::Normal global_normal(
        a_local_normal[0] * e_xi[0] + a_local_normal[1] * e_eta[0] +
            a_local_normal[2] * e_zeta[0],
        a_local_normal[0] * e_xi[1] + a_local_normal[1] * e_eta[1] +
            a_local_normal[2] * e_zeta[1],
        a_local_normal[0] * e_xi[2] + a_local_normal[1] * e_eta[2] +
            a_local_normal[2] * e_zeta[2]);

    global_normal.normalize();
    return global_normal;
  }

  static void GetReconstruction(
      SepUnionMultiFab& interface, SepUnionMultiFab& interface_with_ghost,
      const MultiFab& moments, const Geometry& geom,
      std::vector<InterfaceScalarField>* scalar_fields = nullptr) {
    // ---------------------------------------------------------------------
    // initial plic reconstruction
    // ---------------------------------------------------------------------

    LVIRA::GetReconstruction(interface, interface_with_ghost, moments, geom,
                             nullptr);

    // ---------------------------------------------------------------------
    // Some params
    // ---------------------------------------------------------------------

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

    constexpr double jibben_pu_neighborhood_curvature_threshold = 1.0;
    constexpr double pu_curvature_threshold = 4.0;
    constexpr double pu_delta_factor = 2.5;

    constexpr double pi = 3.141592653589793238462643383279502884;

    // ---------------------------------------------------------------------
    // Some fields
    // ---------------------------------------------------------------------

    const BoxArray& ba = interface_with_ghost.boxArray();
    const DistributionMapping& dm = interface_with_ghost.DistributionMap();

    SepUnionMultiFab jibben_interface(ba, dm, interface_with_ghost.nComp(),
                                      interface_with_ghost.nGrowVect());

    SepUnionMultiFab pu_interface(ba, dm, interface_with_ghost.nComp(),
                                  interface_with_ghost.nGrowVect());

    SepUnionMultiFab pu_neighborhood_interface(
        ba, dm, interface_with_ghost.nComp(), interface_with_ghost.nGrowVect());

    iMultiFab underresolved(ba, dm, 1, 0);
    underresolved.setVal(0);

    // Initialize everything with the initial PLIC interface.
    jibben_interface.LocalCopy(interface_with_ghost, 0, 0,
                               interface_with_ghost.nComp(),
                               interface_with_ghost.nGrowVect());

    pu_interface.LocalCopy(interface_with_ghost, 0, 0,
                           interface_with_ghost.nComp(),
                           interface_with_ghost.nGrowVect());

    pu_neighborhood_interface.LocalCopy(interface_with_ghost, 0, 0,
                                        interface_with_ghost.nComp(),
                                        interface_with_ghost.nGrowVect());

    jibben_interface.FillBoundaryWithPeriodicShift(geom);
    pu_interface.FillBoundaryWithPeriodicShift(geom);
    pu_neighborhood_interface.FillBoundaryWithPeriodicShift(geom);

    // ---------------------------------------------------------------------
    // To output scalar interface fields
    // ---------------------------------------------------------------------

    if (scalar_fields) {
      scalar_fields->clear();

      scalar_fields->emplace_back("is_underresolved", moments.boxArray(),
                                  moments.DistributionMap(), 0);

      scalar_fields->emplace_back("interface_type", moments.boxArray(),
                                  moments.DistributionMap(), 0);
    }

    // ---------------------------------------------------------------------
    // Volume-fitting reconstruction
    // ---------------------------------------------------------------------
    //
    // Accepted VF paraboloids are stored in:
    //   jibben_interface
    //
    // Also, accepted low-curvature VF paraboloids are stored in:
    //   pu_neighborhood_interface
    //
    // Rejected cells leave pu_neighborhood_interface as the original PLIC
    // plane. Therefore the later PU neighborhood can contain both
    // paraboloids
    // and planes.
    // ---------------------------------------------------------------------

    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> jibben_array = jibben_interface.array(mfi);

      Array4<IRL::SeparatorUnion> pu_neigh_array =
          pu_neighborhood_interface.array(mfi);

      Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
          interface_with_ghost.const_array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);

      Array4<int> underresolved_array = underresolved.array(mfi);

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const double liq_vf = moments_array(i, j, k) * vol_inv;

        if (liq_vf < IRL::global_constants::VF_LOW ||
            liq_vf > IRL::global_constants::VF_HIGH) {
          return;
        }

        if (interface_with_ghost_array(i, j, k).type() !=
            IRL::SeparatorUnion::SeparatorType::OnePlane) {
          underresolved_array(i, j, k) = 1;
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

        IRL::JibbenNeighborhood neighborhood;
        neighborhood.reserve(nstencil);
        neighborhood.setDelta(2.5 * dx_avg);
        neighborhood.emptyNeighborhood();

        double vf_supercell = 0.0;
        int jibben_count = 0;

        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              const double neighbor_vf = moments_array(ii, jj, kk) * vol_inv;

              if (neighbor_vf < IRL::global_constants::VF_LOW ||
                  neighbor_vf > IRL::global_constants::VF_HIGH) {
                continue;
              }

              if (interface_with_ghost_array(ii, jj, kk).type() !=
                  IRL::SeparatorUnion::SeparatorType::OnePlane) {
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
                neighborhood.addMember(polygon);

                if (i == ii && j == jj && k == kk) {
                  neighborhood.setCenterOfStencil(jibben_count);
                }

                ++jibben_count;
                vf_supercell += neighbor_vf;
              }
            }
          }
        }

        if (jibben_count < 2) {
          underresolved_array(i, j, k) = 1;
          return;
        }

        // Small isolated VF check
        if (liq_vf < 1.0e-2 && vf_supercell < 1.0e-1) {
          underresolved_array(i, j, k) = 1;
          return;
        }

        neighborhood.localize();
        IRL::Jibben_3D jibben_solver(&neighborhood);

        // Normal scatter check
        const double normal_scatter = jibben_solver.getNormalScatterMetric();
        if (normal_scatter > normal_scatter_threshold) {
          underresolved_array(i, j, k) = 1;
          return;
        }

        // solving for paraboloid
        IRL::Paraboloid jibben_paraboloid = jibben_solver.solve2(&neighborhood);

        // Squared volume error check
        const double volume_error = jibben_solver.getVolumeErrorSquared(dx_avg);
        if (volume_error > volume_error_threshold) {
          underresolved_array(i, j, k) = 1;
          return;
        }

        // Curvature check for accepting Jibben as final.
        const auto& aligned_paraboloid =
            jibben_paraboloid.getAlignedParaboloid();

        if (std::fabs(aligned_paraboloid.a()) * dx_avg > 4.0 ||
            std::fabs(aligned_paraboloid.b()) * dx_avg > 4.0) {
          underresolved_array(i, j, k) = 1;
          return;
        }

        const auto new_datum = jibben_paraboloid.getDatum();
        const auto new_frame = jibben_paraboloid.getReferenceFrame();

        const auto cell = makeCell(i, j, k, problo, dx);

        IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
            solver_distance(cell, liq_vf, 1.0e-14, jibben_paraboloid);

        jibben_paraboloid.setDatum(
            IRL::Pt(new_datum + solver_distance.getDistance() * new_frame[2]));

        jibben_array(i, j, k) = jibben_paraboloid;

        if (std::fabs(aligned_paraboloid.a()) * dx_avg <=
                jibben_pu_neighborhood_curvature_threshold &&
            std::fabs(aligned_paraboloid.b()) * dx_avg <=
                jibben_pu_neighborhood_curvature_threshold) {
          pu_neigh_array(i, j, k) = jibben_array(i, j, k);
        }
      });
    }

    jibben_interface.FillBoundaryWithPeriodicShift(geom);
    pu_neighborhood_interface.FillBoundaryWithPeriodicShift(geom);

    // ---------------------------------------------------------------------
    // PU on underresolved cells only
    // ---------------------------------------------------------------------

    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> pu_array = pu_interface.array(mfi);

      Array4<IRL::SeparatorUnion const> pu_neigh_array =
          pu_neighborhood_interface.const_array(mfi);

      Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
          interface_with_ghost.const_array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);

      Array4<int const> underresolved_array = underresolved.const_array(mfi);

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const double liq_vf = moments_array(i, j, k) * vol_inv;

        if (liq_vf < IRL::global_constants::VF_LOW ||
            liq_vf > IRL::global_constants::VF_HIGH) {
          return;
        }

        // PU is only needed for underresolved cells.
        if (underresolved_array(i, j, k) == 0) {
          return;
        }

        if (interface_with_ghost_array(i, j, k).type() !=
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
          return;
        }

        IRL::PUNeighborhood pu_neighborhood;
        pu_neighborhood.reserve(nstencil);
        pu_neighborhood.emptyNeighborhood();

        int pu_count = 0;

        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              const double neighbor_vf = moments_array(ii, jj, kk) * vol_inv;

              if (neighbor_vf < IRL::global_constants::VF_LOW ||
                  neighbor_vf > IRL::global_constants::VF_HIGH) {
                continue;
              }

              if (interface_with_ghost_array(ii, jj, kk).type() !=
                  IRL::SeparatorUnion::SeparatorType::OnePlane) {
                continue;
              }

              // Polygon/centroid/area from original PLIC.
              const auto planar_separator = IRL::PlanarSeparator::fromOnePlane(
                  interface_with_ghost_array(ii, jj, kk).getPlane());

              const auto polygon =
                  polygonFromPLIC(ii, jj, kk, planar_separator, problo, dx);

              const double polygon_area = polygon.calculateVolume();

              if (polygon.getNumberOfVertices() <= 2 ||
                  !std::isfinite(polygon_area) ||
                  polygon_area <= std::numeric_limits<double>::epsilon()) {
                continue;
              }

              const double area_weight = polygon_area / (dx_avg * dx_avg);

              if (!std::isfinite(area_weight) ||
                  area_weight <= std::numeric_limits<double>::epsilon()) {
                continue;
              }

              double vfrac_weight = 1.0;

              if (neighbor_vf < 0.1) {
                vfrac_weight = 0.5 - 0.5 * std::cos(10.0 * pi * neighbor_vf);
              } else if (neighbor_vf > 0.9) {
                vfrac_weight =
                    0.5 - 0.5 * std::cos(10.0 * pi * (1.0 - neighbor_vf));
              }

              const double weight = area_weight * vfrac_weight;

              const IRL::Pt centroid = polygon.calculateCentroid();

              // Important:
              // Neighbor interface comes from the mixed Jibben/PLIC field.
              // This can be a paraboloid or a plane.
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

        if (pu_count < 2) {
          return;
        }

        const double delta = pu_delta_factor * dx_avg;

        IRL::Paraboloid pu_paraboloid =
            IRL::reconstructionWithPU3D(pu_neighborhood, delta, dx_avg);

        if (!std::isfinite(pu_paraboloid.getDatum()[0]) ||
            !std::isfinite(pu_paraboloid.getDatum()[1]) ||
            !std::isfinite(pu_paraboloid.getDatum()[2])) {
          return;
        }

        const auto& aligned_paraboloid = pu_paraboloid.getAlignedParaboloid();

        if (std::fabs(aligned_paraboloid.a()) * dx_avg >
                pu_curvature_threshold ||
            std::fabs(aligned_paraboloid.b()) * dx_avg >
                pu_curvature_threshold) {
          return;
        }

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
    // 7. PU cleanup for tiny isolated cells
    // ---------------------------------------------------------------------

    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> pu_array = pu_interface.array(mfi);

      Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
          interface_with_ghost.const_array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);

      Array4<int const> underresolved_array = underresolved.const_array(mfi);

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const double liq_vf = moments_array(i, j, k) * vol_inv;

        if (liq_vf < IRL::global_constants::VF_LOW ||
            liq_vf > IRL::global_constants::VF_HIGH) {
          return;
        }

        if (underresolved_array(i, j, k) == 0) {
          return;
        }

        double vf_supercell = 0.0;

        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
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
    // PU if underresolved, else Jibben
    // ---------------------------------------------------------------------

    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);

      Array4<IRL::SeparatorUnion const> jibben_array =
          jibben_interface.const_array(mfi);

      Array4<IRL::SeparatorUnion const> pu_array =
          pu_interface.const_array(mfi);

      Array4<Real const> moments_array = moments.const_array(mfi);

      Array4<int const> underresolved_array = underresolved.const_array(mfi);

      Array4<Real> is_underresolved_poly;
      Array4<Real> is_underresolved_parab;
      Array4<Real> interface_type_poly;
      Array4<Real> interface_type_parab;

      if (scalar_fields) {
        is_underresolved_poly =
            (*scalar_fields)[0].polygon_scalar_data.array(mfi);
        is_underresolved_parab =
            (*scalar_fields)[0].paraboloid_scalar_data.array(mfi);

        interface_type_poly =
            (*scalar_fields)[1].polygon_scalar_data.array(mfi);
        interface_type_parab =
            (*scalar_fields)[1].paraboloid_scalar_data.array(mfi);
      }

      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const double liq_vf = moments_array(i, j, k) * vol_inv;

        if (liq_vf < IRL::global_constants::VF_LOW) {
          interface_array(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
          return;
        } else if (liq_vf > IRL::global_constants::VF_HIGH) {
          interface_array(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
          return;
        }

        const bool is_under = underresolved_array(i, j, k) != 0;

        if (is_under) {
          interface_array(i, j, k) = pu_array(i, j, k);
        } else {
          interface_array(i, j, k) = jibben_array(i, j, k);
        }

        if (scalar_fields) {
          if (interface_array(i, j, k).type() ==
              IRL::SeparatorUnion::SeparatorType::Paraboloid) {
            is_underresolved_parab(i, j, k) = is_under ? 1.0 : 0.0;

            // 1 = accepted Jibben paraboloid
            // 2 = underresolved, using PU paraboloid
            interface_type_parab(i, j, k) = is_under ? 2.0 : 1.0;
          } else {
            is_underresolved_poly(i, j, k) = is_under ? 1.0 : 0.0;

            // 0 = PLIC/plane
            interface_type_poly(i, j, k) = 0.0;
          }
        }
      });
    }

    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());

    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
  }

  // static void GetReconstruction(
  //     SepUnionMultiFab& interface, SepUnionMultiFab& interface_with_ghost,
  //     const MultiFab& moments, const Geometry& geom,
  //     std::vector<InterfaceScalarField>* scalar_fields = nullptr) {
  //   // ---------------------------------------------------------------------
  //   // 1. Initial PLIC reconstruction
  //   // ---------------------------------------------------------------------

  //   LVIRA::GetReconstruction(interface, interface_with_ghost, moments, geom,
  //                            nullptr);

  //   // ---------------------------------------------------------------------
  //   // 2. Parameters
  //   // ---------------------------------------------------------------------

  //   const auto dx = geom.CellSizeArray();
  //   const auto problo = geom.ProbLoArray();

  //   const Real cell_vol = dx[0] * dx[1] * dx[2];
  //   const Real vol_inv = 1.0 / cell_vol;
  //   const Real dx_avg = std::cbrt(dx[0] * dx[1] * dx[2]);

  //   constexpr int nlayers = 1;
  //   constexpr int nstencil =
  //       (1 + 2 * nlayers) * (1 + 2 * nlayers) * (1 + 2 * nlayers);

  //   // Existing HYBRID underresolution thresholds.
  //   constexpr double normal_scatter_threshold = 0.10;
  //   constexpr double volume_error_threshold = 0.05;

  //   // iVF-style well-resolved thresholds.
  //   constexpr int ivf_min_neighbors = 10;
  //   constexpr double alpha_deg = 50.0;
  //   constexpr double pi = 3.141592653589793238462643383279502884;
  //   constexpr double alpha_rad = alpha_deg * pi / 180.0;
  //   const double cos_alpha = std::cos(alpha_rad);

  //   // PU thresholds.
  //   constexpr double jibben_pu_neighborhood_curvature_threshold = 1.0;
  //   constexpr double pu_curvature_threshold = 4.0;
  //   constexpr double pu_delta_factor = 2.5;

  //   // ---------------------------------------------------------------------
  //   // 3. Temporary fields
  //   // ---------------------------------------------------------------------

  //   const BoxArray& ba = interface.boxArray();
  //   const DistributionMapping& dm = interface.DistributionMap();

  //   SepUnionMultiFab jibben_interface(ba, dm, interface.nComp(),
  //                                     interface.nGrowVect());

  //   SepUnionMultiFab pu_interface(ba, dm, interface.nComp(),
  //                                 interface.nGrowVect());

  //   SepUnionMultiFab pu_neighborhood_interface(
  //       ba, dm, interface_with_ghost.nComp(),
  //       interface_with_ghost.nGrowVect());

  //   iMultiFab underresolved(ba, dm, 1, 0);
  //   iMultiFab well_resolved(ba, dm, 1, 0);

  //   underresolved.setVal(0);
  //   well_resolved.setVal(0);

  //   // Initialize everything with the initial PLIC interface.
  //   jibben_interface.LocalCopy(interface, 0, 0, interface.nComp(),
  //                              interface.nGrowVect());

  //   pu_interface.LocalCopy(interface, 0, 0, interface.nComp(),
  //                          interface.nGrowVect());

  //   // This is the mixed field used by the PU fallback neighborhood.
  //   // Initially it is all PLIC. Accepted Jibben paraboloids overwrite it.
  //   pu_neighborhood_interface.LocalCopy(interface_with_ghost, 0, 0,
  //                                       interface_with_ghost.nComp(),
  //                                       interface_with_ghost.nGrowVect());

  //   jibben_interface.FillBoundaryWithPeriodicShift(geom);
  //   pu_interface.FillBoundaryWithPeriodicShift(geom);
  //   pu_neighborhood_interface.FillBoundaryWithPeriodicShift(geom);

  //   // ---------------------------------------------------------------------
  //   // 4. Scalar fields
  //   // ---------------------------------------------------------------------

  //   if (scalar_fields) {
  //     scalar_fields->clear();

  //     scalar_fields->emplace_back("is_underresolved", moments.boxArray(),
  //                                 moments.DistributionMap(), 0);

  //     scalar_fields->emplace_back("is_well_resolved", moments.boxArray(),
  //                                 moments.DistributionMap(), 0);

  //     scalar_fields->emplace_back("interface_type", moments.boxArray(),
  //                                 moments.DistributionMap(), 0);
  //   }

  //   // ---------------------------------------------------------------------
  //   // 5. Jibben pass with iVF-style well-resolved override
  //   // ---------------------------------------------------------------------
  //   //
  //   // This does the expensive Jibben fit only once.
  //   //
  //   // For each mixed cell:
  //   //   1. Build PLIC polygon neighborhood.
  //   //   2. Fit Jibben paraboloid.
  //   //   3. Compute iVF-style averaged fitted normal using iJibben.
  //   //   4. Compare initial PLIC normal with this averaged fitted normal.
  //   //   5. Filter neighbors using the averaged fitted normal.
  //   //   6. If iVF-style conditions pass:
  //   //        mark well_resolved = 1
  //   //        skip normal-scatter / volume-error / curvature underresolved
  //   //        checks
  //   //      else:
  //   //        use existing HYBRID underresolved checks
  //   // ---------------------------------------------------------------------

  //   for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
  //     Array4<IRL::SeparatorUnion> jibben_array = jibben_interface.array(mfi);

  //     Array4<IRL::SeparatorUnion> pu_neigh_array =
  //         pu_neighborhood_interface.array(mfi);

  //     Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
  //         interface_with_ghost.const_array(mfi);

  //     Array4<Real const> moments_array = moments.const_array(mfi);

  //     Array4<int> underresolved_array = underresolved.array(mfi);
  //     Array4<int> well_resolved_array = well_resolved.array(mfi);

  //     const Box& bx = mfi.tilebox();

  //     amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
  //       const double liq_vf = moments_array(i, j, k) * vol_inv;

  //       if (liq_vf < IRL::global_constants::VF_LOW ||
  //           liq_vf > IRL::global_constants::VF_HIGH) {
  //         return;
  //       }

  //       if (interface_with_ghost_array(i, j, k).type() !=
  //           IRL::SeparatorUnion::SeparatorType::OnePlane) {
  //         underresolved_array(i, j, k) = 1;
  //         return;
  //       }

  //       const auto center_planar_separator =
  //       IRL::PlanarSeparator::fromOnePlane(
  //           interface_with_ghost_array(i, j, k).getPlane());

  //       const auto center_polygon =
  //           polygonFromPLIC(i, j, k, center_planar_separator, problo, dx);

  //       const double center_area = center_polygon.calculateVolume();

  //       if (center_polygon.getNumberOfVertices() <= 2 ||
  //           !std::isfinite(center_area) ||
  //           center_area <= std::numeric_limits<double>::epsilon()) {
  //         underresolved_array(i, j, k) = 1;
  //         return;
  //       }

  //       //
  //       -------------------------------------------------------------------
  //       // Build the Jibben neighborhood.
  //       //
  //       -------------------------------------------------------------------

  //       IRL::JibbenNeighborhood neighborhood;
  //       neighborhood.reserve(nstencil);
  //       neighborhood.setDelta(2.5 * dx_avg);
  //       neighborhood.emptyNeighborhood();

  //       double vf_supercell = 0.0;
  //       int jibben_count = 0;

  //       for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
  //         for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
  //           for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
  //             const double neighbor_vf = moments_array(ii, jj, kk) * vol_inv;

  //             if (neighbor_vf < IRL::global_constants::VF_LOW ||
  //                 neighbor_vf > IRL::global_constants::VF_HIGH) {
  //               continue;
  //             }

  //             if (interface_with_ghost_array(ii, jj, kk).type() !=
  //                 IRL::SeparatorUnion::SeparatorType::OnePlane) {
  //               continue;
  //             }

  //             const auto planar_separator =
  //             IRL::PlanarSeparator::fromOnePlane(
  //                 interface_with_ghost_array(ii, jj, kk).getPlane());

  //             const auto polygon =
  //                 polygonFromPLIC(ii, jj, kk, planar_separator, problo, dx);

  //             const double polygon_area = polygon.calculateVolume();

  //             if (polygon.getNumberOfVertices() > 2 &&
  //                 std::isfinite(polygon_area) &&
  //                 polygon_area > std::numeric_limits<double>::epsilon()) {
  //               neighborhood.addMember(polygon);

  //               if (i == ii && j == jj && k == kk) {
  //                 neighborhood.setCenterOfStencil(jibben_count);
  //               }

  //               ++jibben_count;
  //               vf_supercell += neighbor_vf;
  //             }
  //           }
  //         }
  //       }

  //       if (jibben_count < 2) {
  //         underresolved_array(i, j, k) = 1;
  //         return;
  //       }

  //       neighborhood.localize();

  //       //
  //       -------------------------------------------------------------------
  //       // Primary HYBRID / Jibben fit.
  //       //
  //       -------------------------------------------------------------------

  //       IRL::Jibben_3D jibben_solver(&neighborhood);

  //       IRL::Paraboloid jibben_paraboloid =
  //       jibben_solver.solve2(&neighborhood);

  //       //
  //       -------------------------------------------------------------------
  //       // iVF-style averaged-normal diagnostic.
  //       //
  //       // This does NOT change the PLIC normal.
  //       // It only computes the same averaged normal that iVF would use as
  //       its
  //       // updated normal, and then checks whether this neighborhood is
  //       // coherent.
  //       //
  //       -------------------------------------------------------------------

  //       bool is_well_resolved_ivf_style = false;

  //       IRL::iJibben_3D ijibben(&neighborhood);
  //       ijibben.getParaboloidCoefficients2();

  //       const std::pair<double, IRL::Normal> fit_result =
  //           ijibben.computeAveragedCurvatureAndNormal();

  //       const IRL::Normal updated_local_n = fit_result.second;

  //       IRL::Normal updated_global_n = localToGlobalNormal(
  //           updated_local_n, neighborhood.getReferenceFrame());

  //       updated_global_n.normalize();

  //       IRL::Normal init_n =
  //           interface_with_ghost_array(i, j, k).getPlane().normal();

  //       init_n.normalize();

  //       const double target_alignment = init_n * updated_global_n;

  //       const bool target_normal_ok = target_alignment >= cos_alpha;

  //       int filtered_count = 0;
  //       int filtered_center_count = -1;

  //       if (target_normal_ok) {
  //         for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
  //           for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
  //             for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
  //               const double neighbor_vf = moments_array(ii, jj, kk) *
  //               vol_inv;

  //               if (neighbor_vf < IRL::global_constants::VF_LOW ||
  //                   neighbor_vf > IRL::global_constants::VF_HIGH) {
  //                 continue;
  //               }

  //               if (interface_with_ghost_array(ii, jj, kk).type() !=
  //                   IRL::SeparatorUnion::SeparatorType::OnePlane) {
  //                 continue;
  //               }

  //               IRL::Normal neighbor_init_n =
  //                   interface_with_ghost_array(ii, jj,
  //                   kk).getPlane().normal();

  //               neighbor_init_n.normalize();

  //               const double neighbor_alignment =
  //                   neighbor_init_n * updated_global_n;

  //               if (neighbor_alignment < cos_alpha) {
  //                 continue;
  //               }

  //               const auto neighbor_plic =
  //               IRL::PlanarSeparator::fromOnePlane(
  //                   interface_with_ghost_array(ii, jj, kk).getPlane());

  //               const auto neighbor_polygon =
  //                   polygonFromPLIC(ii, jj, kk, neighbor_plic, problo, dx);

  //               const double neighbor_area =
  //               neighbor_polygon.calculateVolume();

  //               if (neighbor_polygon.getNumberOfVertices() <= 2 ||
  //                   !std::isfinite(neighbor_area) ||
  //                   neighbor_area <= std::numeric_limits<double>::epsilon())
  //                   {
  //                 continue;
  //               }

  //               if (ii == i && jj == j && kk == k) {
  //                 filtered_center_count = filtered_count;
  //               }

  //               ++filtered_count;
  //             }
  //           }
  //         }
  //       }

  //       if (target_normal_ok && filtered_count >= ivf_min_neighbors &&
  //           filtered_center_count >= 0) {
  //         is_well_resolved_ivf_style = true;
  //         well_resolved_array(i, j, k) = 1;
  //         underresolved_array(i, j, k) = 0;
  //       }

  //       //
  //       -------------------------------------------------------------------
  //       // Existing HYBRID underresolved checks.
  //       //
  //       // These are skipped if the iVF-style well-resolved test passed.
  //       //
  //       -------------------------------------------------------------------

  //       const auto& aligned_paraboloid =
  //           jibben_paraboloid.getAlignedParaboloid();

  //       if (!is_well_resolved_ivf_style) {
  //         // Normal scatter check.
  //         const double normal_scatter =
  //         jibben_solver.getNormalScatterMetric();

  //         if (normal_scatter > normal_scatter_threshold) {
  //           underresolved_array(i, j, k) = 1;
  //         }

  //         // Squared volume error check.
  //         const double volume_error =
  //             jibben_solver.getVolumeErrorSquared(dx_avg);

  //         if (volume_error > volume_error_threshold) {
  //           underresolved_array(i, j, k) = 1;
  //         }

  //         // Small isolated VF check.
  //         if (liq_vf < 1.0e-2 && vf_supercell < 1.0e-1) {
  //           underresolved_array(i, j, k) = 1;
  //         }

  //         // Curvature check.
  //         if (std::fabs(aligned_paraboloid.a()) * dx_avg > 4.0 ||
  //             std::fabs(aligned_paraboloid.b()) * dx_avg > 4.0) {
  //           underresolved_array(i, j, k) = 1;
  //         }
  //       }

  //       // If still underresolved, keep jibben_interface as PLIC.
  //       // The PU fallback pass will try to improve it.
  //       if (underresolved_array(i, j, k) != 0) {
  //         return;
  //       }

  //       //
  //       -------------------------------------------------------------------
  //       // Volume-match accepted Jibben paraboloid.
  //       //
  //       -------------------------------------------------------------------

  //       const auto new_datum = jibben_paraboloid.getDatum();
  //       const auto new_frame = jibben_paraboloid.getReferenceFrame();

  //       const auto cell = makeCell(i, j, k, problo, dx);

  //       IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
  //           solver_distance(cell, liq_vf, 1.0e-14, jibben_paraboloid);

  //       jibben_paraboloid.setDatum(
  //           IRL::Pt(new_datum + solver_distance.getDistance() *
  //           new_frame[2]));

  //       jibben_array(i, j, k) = jibben_paraboloid;

  //       //
  //       -------------------------------------------------------------------
  //       // Update mixed PU-neighborhood field.
  //       //
  //       // If iVF-style well-resolved, trust this paraboloid in the later PU
  //       // neighborhood. Otherwise only use it if its curvature is mild.
  //       //
  //       -------------------------------------------------------------------

  //       if (is_well_resolved_ivf_style ||
  //           (std::fabs(aligned_paraboloid.a()) * dx_avg <=
  //                jibben_pu_neighborhood_curvature_threshold &&
  //            std::fabs(aligned_paraboloid.b()) * dx_avg <=
  //                jibben_pu_neighborhood_curvature_threshold)) {
  //         pu_neigh_array(i, j, k) = jibben_array(i, j, k);
  //       }
  //     });
  //   }

  //   jibben_interface.FillBoundaryWithPeriodicShift(geom);
  //   pu_neighborhood_interface.FillBoundaryWithPeriodicShift(geom);

  //   // ---------------------------------------------------------------------
  //   // 6. PU fallback pass
  //   // ---------------------------------------------------------------------
  //   //
  //   // PU is run only on cells still marked underresolved.
  //   //
  //   // Geometry/centroid/area are built from the original PLIC field.
  //   // The actual neighbor interface used by PU comes from
  //   // pu_neighborhood_interface, which can contain accepted Jibben
  //   paraboloids
  //   // and original PLIC planes.
  //   // ---------------------------------------------------------------------

  //   for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
  //     Array4<IRL::SeparatorUnion> pu_array = pu_interface.array(mfi);

  //     Array4<IRL::SeparatorUnion const> pu_neigh_array =
  //         pu_neighborhood_interface.const_array(mfi);

  //     Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
  //         interface_with_ghost.const_array(mfi);

  //     Array4<Real const> moments_array = moments.const_array(mfi);

  //     Array4<int const> underresolved_array = underresolved.const_array(mfi);

  //     const Box& bx = mfi.tilebox();

  //     amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
  //       const double liq_vf = moments_array(i, j, k) * vol_inv;

  //       if (liq_vf < IRL::global_constants::VF_LOW ||
  //           liq_vf > IRL::global_constants::VF_HIGH) {
  //         return;
  //       }

  //       if (underresolved_array(i, j, k) == 0) {
  //         return;
  //       }

  //       if (interface_with_ghost_array(i, j, k).type() !=
  //           IRL::SeparatorUnion::SeparatorType::OnePlane) {
  //         return;
  //       }

  //       const auto center_planar_separator =
  //       IRL::PlanarSeparator::fromOnePlane(
  //           interface_with_ghost_array(i, j, k).getPlane());

  //       const auto center_polygon =
  //           polygonFromPLIC(i, j, k, center_planar_separator, problo, dx);

  //       const double center_area = center_polygon.calculateVolume();

  //       if (center_polygon.getNumberOfVertices() <= 2 ||
  //           !std::isfinite(center_area) ||
  //           center_area <= std::numeric_limits<double>::epsilon()) {
  //         return;
  //       }

  //       IRL::PUNeighborhood pu_neighborhood;
  //       pu_neighborhood.reserve(nstencil);
  //       pu_neighborhood.emptyNeighborhood();

  //       int pu_count = 0;

  //       for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
  //         for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
  //           for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
  //             const double neighbor_vf = moments_array(ii, jj, kk) * vol_inv;

  //             if (neighbor_vf < IRL::global_constants::VF_LOW ||
  //                 neighbor_vf > IRL::global_constants::VF_HIGH) {
  //               continue;
  //             }

  //             if (interface_with_ghost_array(ii, jj, kk).type() !=
  //                 IRL::SeparatorUnion::SeparatorType::OnePlane) {
  //               continue;
  //             }

  //             // Polygon geometry comes from original PLIC.
  //             const auto planar_separator =
  //             IRL::PlanarSeparator::fromOnePlane(
  //                 interface_with_ghost_array(ii, jj, kk).getPlane());

  //             const auto polygon =
  //                 polygonFromPLIC(ii, jj, kk, planar_separator, problo, dx);

  //             const double polygon_area = polygon.calculateVolume();

  //             if (polygon.getNumberOfVertices() <= 2 ||
  //                 !std::isfinite(polygon_area) ||
  //                 polygon_area <= std::numeric_limits<double>::epsilon()) {
  //               continue;
  //             }

  //             const double area_weight = polygon_area / (dx_avg * dx_avg);

  //             if (!std::isfinite(area_weight) ||
  //                 area_weight <= std::numeric_limits<double>::epsilon()) {
  //               continue;
  //             }

  //             double vfrac_weight = 1.0;

  //             if (neighbor_vf < 0.1) {
  //               vfrac_weight = 0.5 - 0.5 * std::cos(10.0 * pi * neighbor_vf);
  //             } else if (neighbor_vf > 0.9) {
  //               vfrac_weight =
  //                   0.5 - 0.5 * std::cos(10.0 * pi * (1.0 - neighbor_vf));
  //             }

  //             const double weight = area_weight * vfrac_weight;

  //             const IRL::Pt centroid = polygon.calculateCentroid();

  //             // Mixed neighbor reconstruction:
  //             // this may be a Jibben paraboloid or a PLIC plane.
  //             const IRL::SeparatorVariant neighbor_interface =
  //                 pu_neigh_array(ii, jj, kk);

  //             pu_neighborhood.addMember(neighbor_interface, centroid,
  //             weight);

  //             if (i == ii && j == jj && k == kk) {
  //               pu_neighborhood.setCenterOfStencil(pu_count);
  //             }

  //             ++pu_count;
  //           }
  //         }
  //       }

  //       if (pu_count < 2) {
  //         return;
  //       }

  //       const double delta = pu_delta_factor * dx_avg;

  //       IRL::Paraboloid pu_paraboloid =
  //           IRL::reconstructionWithPU3D(pu_neighborhood, delta, dx_avg);

  //       if (!std::isfinite(pu_paraboloid.getDatum()[0]) ||
  //           !std::isfinite(pu_paraboloid.getDatum()[1]) ||
  //           !std::isfinite(pu_paraboloid.getDatum()[2])) {
  //         return;
  //       }

  //       const auto& aligned_paraboloid =
  //       pu_paraboloid.getAlignedParaboloid();

  //       if (std::fabs(aligned_paraboloid.a()) * dx_avg >
  //               pu_curvature_threshold ||
  //           std::fabs(aligned_paraboloid.b()) * dx_avg >
  //               pu_curvature_threshold) {
  //         return;
  //       }

  //       const auto new_datum = pu_paraboloid.getDatum();
  //       const auto new_frame = pu_paraboloid.getReferenceFrame();

  //       const auto cell = makeCell(i, j, k, problo, dx);

  //       IRL::ProgressiveDistanceSolverParaboloid<IRL::RectangularCuboid>
  //           solver_distance(cell, liq_vf, 1.0e-14, pu_paraboloid);

  //       pu_paraboloid.setDatum(
  //           IRL::Pt(new_datum + solver_distance.getDistance() *
  //           new_frame[2]));

  //       pu_array(i, j, k) = pu_paraboloid;
  //     });
  //   }

  //   pu_interface.FillBoundaryWithPeriodicShift(geom);

  //   // ---------------------------------------------------------------------
  //   // 7. PU cleanup for tiny isolated cells
  //   // ---------------------------------------------------------------------

  //   for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
  //     Array4<IRL::SeparatorUnion> pu_array = pu_interface.array(mfi);

  //     Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
  //         interface_with_ghost.const_array(mfi);

  //     Array4<Real const> moments_array = moments.const_array(mfi);

  //     Array4<int const> underresolved_array = underresolved.const_array(mfi);

  //     const Box& bx = mfi.tilebox();

  //     amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
  //       const double liq_vf = moments_array(i, j, k) * vol_inv;

  //       if (liq_vf < IRL::global_constants::VF_LOW ||
  //           liq_vf > IRL::global_constants::VF_HIGH) {
  //         return;
  //       }

  //       if (underresolved_array(i, j, k) == 0) {
  //         return;
  //       }

  //       double vf_supercell = 0.0;

  //       for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
  //         for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
  //           for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
  //             vf_supercell += moments_array(ii, jj, kk) * vol_inv;
  //           }
  //         }
  //       }

  //       if (liq_vf < 1.0e-2 && vf_supercell < 1.0e-1) {
  //         pu_array(i, j, k) = interface_with_ghost_array(i, j, k);
  //       }
  //     });
  //   }

  //   pu_interface.FillBoundaryWithPeriodicShift(geom);

  //   // ---------------------------------------------------------------------
  //   // 8. Final hybrid selection
  //   // ---------------------------------------------------------------------
  //   //
  //   // If underresolved:
  //   //   use PU result.
  //   //
  //   // Else:
  //   //   use accepted Jibben result.
  //   //
  //   // If the iVF-style well-resolved check passed, the cell is protected
  //   from
  //   // the usual underresolved checks and therefore uses Jibben.
  //   // ---------------------------------------------------------------------

  //   for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
  //     Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);

  //     Array4<IRL::SeparatorUnion const> jibben_array =
  //         jibben_interface.const_array(mfi);

  //     Array4<IRL::SeparatorUnion const> pu_array =
  //         pu_interface.const_array(mfi);

  //     Array4<Real const> moments_array = moments.const_array(mfi);

  //     Array4<int const> underresolved_array = underresolved.const_array(mfi);

  //     Array4<int const> well_resolved_array = well_resolved.const_array(mfi);

  //     Array4<Real> is_underresolved_poly;
  //     Array4<Real> is_underresolved_parab;
  //     Array4<Real> is_well_resolved_poly;
  //     Array4<Real> is_well_resolved_parab;
  //     Array4<Real> interface_type_poly;
  //     Array4<Real> interface_type_parab;

  //     if (scalar_fields) {
  //       is_underresolved_poly =
  //           (*scalar_fields)[0].polygon_scalar_data.array(mfi);
  //       is_underresolved_parab =
  //           (*scalar_fields)[0].paraboloid_scalar_data.array(mfi);

  //       is_well_resolved_poly =
  //           (*scalar_fields)[1].polygon_scalar_data.array(mfi);
  //       is_well_resolved_parab =
  //           (*scalar_fields)[1].paraboloid_scalar_data.array(mfi);

  //       interface_type_poly =
  //           (*scalar_fields)[2].polygon_scalar_data.array(mfi);
  //       interface_type_parab =
  //           (*scalar_fields)[2].paraboloid_scalar_data.array(mfi);
  //     }

  //     const Box& bx = mfi.tilebox();

  //     amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
  //       const double liq_vf = moments_array(i, j, k) * vol_inv;

  //       if (liq_vf < IRL::global_constants::VF_LOW) {
  //         interface_array(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
  //         return;
  //       } else if (liq_vf > IRL::global_constants::VF_HIGH) {
  //         interface_array(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
  //         return;
  //       }

  //       const bool is_under = underresolved_array(i, j, k) != 0;

  //       const bool is_well = well_resolved_array(i, j, k) != 0;

  //       if (is_under) {
  //         interface_array(i, j, k) = pu_array(i, j, k);
  //       } else {
  //         interface_array(i, j, k) = jibben_array(i, j, k);
  //       }

  //       if (scalar_fields) {
  //         if (interface_array(i, j, k).type() ==
  //             IRL::SeparatorUnion::SeparatorType::Paraboloid) {
  //           is_underresolved_parab(i, j, k) = is_under ? 1.0 : 0.0;
  //           is_well_resolved_parab(i, j, k) = is_well ? 1.0 : 0.0;

  //           // interface_type:
  //           // 1 = accepted Jibben
  //           // 2 = PU fallback
  //           // 3 = accepted Jibben because iVF-style well-resolved
  //           // interface_type_parab(i, j, k) =
  //           //     is_under ? 2.0 : (is_well ? 3.0 : 1.0);
  //           interface_type_parab(i, j, k) = is_under ? 2.0 : 1.0;
  //         } else {
  //           is_underresolved_poly(i, j, k) = is_under ? 1.0 : 0.0;
  //           is_well_resolved_poly(i, j, k) = is_well ? 1.0 : 0.0;

  //           // 0 = plane / PLIC-like fallback
  //           interface_type_poly(i, j, k) = 0.0;
  //         }
  //       }
  //     });
  //   }

  //   // ---------------------------------------------------------------------
  //   // 9. Final ghost fill / periodic datum correction
  //   // ---------------------------------------------------------------------

  //   interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
  //                                  interface.nGrowVect());

  //   interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
  // }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_HYBRID_H_
