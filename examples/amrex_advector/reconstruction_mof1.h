// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_MOF1_H_
#define EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_MOF1_H_

#include "irl/amrex/sepunion_multifab.h"

using namespace amrex;

struct MOF1 {
  static void GetReconstruction(SepUnionMultiFab& interface,
                                SepUnionMultiFab& interface_with_ghost,
                                const MultiFab& moments, const Geometry& geom) {
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

        // Construct local cell as IRL object
        const double x = problo[0] + i * dx[0];
        const double y = problo[1] + j * dx[1];
        const double z = problo[2] + k * dx[2];
        const IRL::Pt lower_cell_pt(x, y, z);
        const IRL::Pt upper_cell_pt(x + dx[0], y + dx[1], z + dx[2]);
        const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
        const auto cell = IRL::RectangularCuboid::fromBoundingPts(
            lower_cell_pt, upper_cell_pt);
        // Compute gas and liquid moments (M0 and M1)
        const double liq_m0 = moments_array(i, j, k, 0);
        const double gas_m0 = vol - liq_m0;
        const IRL::Pt liq_m1 =
            (1.0 / liq_m0) * IRL::Pt(moments_array(i, j, k, 1),
                                     moments_array(i, j, k, 2),
                                     moments_array(i, j, k, 3));
        const IRL::Pt gas_m1 =
            (1.0 / gas_m0) *
            (vol * cell_center - IRL::Pt(moments_array(i, j, k, 1),
                                         moments_array(i, j, k, 2),
                                         moments_array(i, j, k, 3)));
        const IRL::SeparatedMoments<IRL::VolumeMoments> svm(
            IRL::VolumeMoments(liq_m0, liq_m1),
            IRL::VolumeMoments(gas_m0, gas_m1));
        // Reconstruct interface with ELVIRA and store result in
        // interface_array
        interface_array(i, j, k) =
            IRL::reconstructionWithMOF3D(cell, svm, 0.5, 0.5);
      });
    }  // end mfi

    // This copies the interface, update ghosts, and shifts datums across
    // periodic boundaries
    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());
    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_MOF1_H_
