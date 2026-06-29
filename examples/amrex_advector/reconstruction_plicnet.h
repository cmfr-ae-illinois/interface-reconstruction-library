// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_PLICNET_H_
#define EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_PLICNET_H_

#include "irl/amrex/sepunion_multifab.h"

using namespace amrex;

struct PLICNet {
  static void GetReconstruction(
      SepUnionMultiFab& interface, SepUnionMultiFab& interface_with_ghost,
      const MultiFab& moments, const Geometry& geom,
      std::vector<InterfaceScalarField>* scalar_fields = nullptr) {
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

        // Fill in PLICNet neighborhood
        IRL::PLICNet plicnet;
        for (int kk = k - 1; kk <= k + 1; ++kk) {
          for (int jj = j - 1; jj <= j + 1; ++jj) {
            for (int ii = i - 1; ii <= i + 1; ++ii) {
              const double xx = problo[0] + ii * dx[0];
              const double yy = problo[1] + jj * dx[1];
              const double zz = problo[2] + kk * dx[2];
              const IRL::Pt lower_cell_pt = IRL::Pt(xx, yy, zz);
              const IRL::Pt upper_cell_pt =
                  IRL::Pt(xx + dx[0], yy + dx[1], zz + dx[2]);
              const IRL::Pt cell_center = 0.5 * (lower_cell_pt + upper_cell_pt);
              // Compute gas and liquid moments (M0 and M1)
              const double liq_m0 = moments_array(ii, jj, kk, 0);
              const double gas_m0 = vol - liq_m0;
              const IRL::Pt liq_m1 = (1.0 / IRL::safelyTiny(liq_m0)) *
                                     IRL::Pt(moments_array(ii, jj, kk, 1),
                                             moments_array(ii, jj, kk, 2),
                                             moments_array(ii, jj, kk, 3));
              const IRL::Pt gas_m1 =
                  (1.0 / IRL::safelyTiny(gas_m0)) *
                  IRL::Pt(vol * cell_center -
                          IRL::Pt(moments_array(ii, jj, kk, 1),
                                  moments_array(ii, jj, kk, 2),
                                  moments_array(ii, jj, kk, 3)));
              plicnet.setMember(lower_cell_pt, upper_cell_pt, liq_m0 * vol_inv,
                                liq_m1, gas_m1, ii - i, jj - j, kk - k);
            }
          }
        }
        // Reconstruct interface with PLICNet and store result in
        // interface_array
        interface_array(i, j, k) = plicnet.getPlanarSeparator();
      });
    }  // end mfi

    // This copies the interface, update ghosts, and shifts datums across
    // periodic boundaries
    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());
    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_PLICNET_H_
