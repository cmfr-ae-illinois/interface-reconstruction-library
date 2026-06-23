// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_LVIRA_H_
#define EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_LVIRA_H_

#include "irl/amrex/sepunion_multifab.h"

#include "examples/amrex_advector/reconstruction_elvira.h"

using namespace amrex;

struct LVIRA {
  static void GetReconstruction(SepUnionMultiFab& interface,
                                SepUnionMultiFab& interface_with_ghost,
                                const MultiFab& moments, const Geometry& geom) {
    // Produce initial guess with ELVIRA
    ELVIRA::GetReconstruction(interface, interface_with_ghost, moments, geom);

    // Now compute LVIRA reconstruction
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    const Real vol = dx[0] * dx[1] * dx[2];
    const Real vol_inv = 1.0 / vol;

    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);
      Array4<Real const> moments_array = moments.const_array(mfi);
      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // Compute cell volume fraction
        const double vfrac = moments_array(i, j, k) * vol_inv;
        // Skip cell if vfrac ~0 or vfrac~1
        if (vfrac < IRL::global_constants::VF_LOW ||
            vfrac > IRL::global_constants::VF_HIGH) {
          const double distance = std::copysign(1.0, vfrac - 0.5);
          interface_array(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
              IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          return;
        }

        // Fill in ELVIRA neighborhood
        IRL::LVIRANeighborhood<IRL::RectangularCuboid> lvira_neighborhood;
        lvira_neighborhood.resize(27);
        IRL::RectangularCuboid cells[27];
        std::array<double, 27> cells_vfrac;
        for (int kk = k - 1; kk < k + 2; ++kk) {
          for (int jj = j - 1; jj < j + 2; ++jj) {
            for (int ii = i - 1; ii < i + 2; ++ii) {
              const int neigh_id =
                  (kk - k + 1) * 9 + (jj - j + 1) * 3 + (ii - i + 1);
              const double xx = problo[0] + ii * dx[0];
              const double yy = problo[1] + jj * dx[1];
              const double zz = problo[2] + kk * dx[2];
              cells[neigh_id] = IRL::RectangularCuboid::fromBoundingPts(
                  IRL::Pt(xx, yy, zz),
                  IRL::Pt(xx + dx[0], yy + dx[1], zz + dx[2]));
              cells_vfrac[neigh_id] = moments_array(ii, jj, kk, 0) * vol_inv;
              lvira_neighborhood.setMember(
                  static_cast<IRL::UnsignedIndex_t>(neigh_id), &cells[neigh_id],
                  &cells_vfrac[neigh_id]);
              // Trap center cell
              if (ii == i && jj == j && kk == k) {
                lvira_neighborhood.setCenterOfStencil(neigh_id);
              }
            }
          }
        }
        // Reconstruct interface with LVIRA and store result in
        // interface_array
        const auto planar_separator = IRL::PlanarSeparator::fromOnePlane(
            interface_array(i, j, k).getPlane());
        interface_array(i, j, k) = IRL::reconstructionWithLVIRA3D(
            lvira_neighborhood, planar_separator);
      });
    }  // end mfi

    // This copies the interface, update ghosts, and shifts datums across
    // periodic boundaries
    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());
    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_LVIRA_H_
