// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_VF_H_
#define EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_VF_H_

#include "irl/amrex/sepunion_multifab.h"

#include "examples/amrex_advector/reconstruction_lvira.h"
#include "examples/amrex_advector/reconstruction_plicnet.h"

using namespace amrex;

struct VF {
  static void GetReconstruction(
      SepUnionMultiFab& interface, SepUnionMultiFab& interface_with_ghost,
      const MultiFab& moments, const Geometry& geom,
      std::vector<InterfaceScalarField>* scalar_fields = nullptr) {
    // Produce initial guess with PLICNet
    // PLICNet::GetReconstruction(interface, interface_with_ghost, moments,
    // geom);
    LVIRA::GetReconstruction(interface, interface_with_ghost, moments, geom,
                             scalar_fields);

    // Now compute Jibben's reconstruction
    const auto dx = geom.CellSizeArray();
    const auto dx_avg = std::cbrt(dx[0] * dx[1] * dx[2]);
    const auto problo = geom.ProbLoArray();
    const Real vol = dx[0] * dx[1] * dx[2];
    const Real vol_inv = 1.0 / vol;

    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);
      Array4<IRL::SeparatorUnion const> interface_with_ghost_array =
          interface_with_ghost.const_array(mfi);
      Array4<Real const> moments_array = moments.const_array(mfi);
      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // Compute cell volume fraction
        const double vfrac = moments_array(i, j, k) * vol_inv;
        // Skip cell if vfrac ~0 or vfrac~1
        if (vfrac < IRL::global_constants::VF_LOW) {
          interface_array(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
          return;
        } else if (vfrac > IRL::global_constants::VF_HIGH) {
          interface_array(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
          return;
        }

        // Fill in Jibben neighborhood
        IRL::JibbenNeighborhood neighborhood;
        const int nlayers = 1;
        const int nstencil =
            (1 + 2 * nlayers) * (1 + 2 * nlayers) * (1 + 2 * nlayers);
        neighborhood.reserve(nstencil);
        neighborhood.setDelta(2.5 * dx_avg);
        neighborhood.emptyNeighborhood();
        int count = 0;
        for (int kk = k - nlayers; kk <= k + nlayers; ++kk) {
          for (int jj = j - nlayers; jj <= j + nlayers; ++jj) {
            for (int ii = i - nlayers; ii <= i + nlayers; ++ii) {
              const double vfrac_local = moments_array(ii, jj, kk) * vol_inv;
              if (vfrac_local < IRL::global_constants::VF_LOW ||
                  vfrac_local > IRL::global_constants::VF_HIGH) {
                continue;
              }
              const double xx = problo[0] + ii * dx[0];
              const double yy = problo[1] + jj * dx[1];
              const double zz = problo[2] + kk * dx[2];
              const auto cell = IRL::RectangularCuboid::fromBoundingPts(
                  IRL::Pt(xx, yy, zz),
                  IRL::Pt(xx + dx[0], yy + dx[1], zz + dx[2]));
              const auto planar_separator = IRL::PlanarSeparator::fromOnePlane(
                  interface_with_ghost_array(ii, jj, kk).getPlane());
              const auto polygon =
                  IRL::getPlanePolygonFromReconstruction<IRL::Polygon>(
                      cell, planar_separator, planar_separator[0]);
              if (polygon.getNumberOfVertices() > 0) {
                neighborhood.addMember(polygon);
                if (i == ii && j == jj && k == kk) {
                  neighborhood.setCenterOfStencil(count);
                }
                count++;
              }
            }
          }
        }
        neighborhood.localize();
        // Reconstruct interface with Jibben and store result in
        // interface_array
        auto paraboloid = IRL::reconstructionWithJibben3D(neighborhood);

        // Construct local cell
        const double x = problo[0] + i * dx[0];
        const double y = problo[1] + j * dx[1];
        const double z = problo[2] + k * dx[2];
        const IRL::Pt lower_cell_pt(x, y, z);
        const IRL::Pt upper_cell_pt(x + dx[0], y + dx[1], z + dx[2]);
        const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            lower_cell_pt, upper_cell_pt);
        const double distance_to_datum =
            IRL::magnitude(paraboloid.getDatum() - cell_center);
        const double max_curvature =
            std::max(std::fabs(paraboloid.getAlignedParaboloid().a()),
                     std::fabs(paraboloid.getAlignedParaboloid().b()));

        // if (distance_to_datum < 2.0 * dx_avg && max_curvature * dx_avg < 0.5)
        // { Assign paraboloid to interface array
        interface_array(i, j, k) = paraboloid;
        // Match to volume fraction
        IRL::setDistanceToMatchVolumeFraction(
            cell, vfrac, &interface_array(i, j, k), 1.0e-14);
        // }
      });
    }  // end mfi

    // This copies the interface, update ghosts, and shifts datums across
    // periodic boundaries
    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());
    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_VF_H_
