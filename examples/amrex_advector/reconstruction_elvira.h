// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_ELVIRA_H_
#define EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_ELVIRA_H_

#include "irl/amrex/sepunion_multifab.h"

using namespace amrex;

class ReconstructionLoopTimer {
 public:
  explicit ReconstructionLoopTimer(amrex::Real* accumulator)
      : accumulator_m(accumulator) {}

  void start() {
    if (accumulator_m != nullptr) start_m = amrex::second();
  }

  void stop() {
    if (accumulator_m != nullptr) {
      *accumulator_m += amrex::second() - start_m;
    }
  }

 private:
  amrex::Real* accumulator_m;
  amrex::Real start_m = 0.0;
};

struct ELVIRA {
  static void GetReconstruction(
      SepUnionMultiFab& interface, SepUnionMultiFab& interface_with_ghost,
      const MultiFab& moments, const Geometry& geom,
      std::vector<InterfaceScalarField>* scalar_fields = nullptr,
      Real* reconstruction_loop_time = nullptr) {
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    ReconstructionLoopTimer loop_timer(reconstruction_loop_time);
    loop_timer.start();
    for (MFIter mfi(interface, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Array4<IRL::SeparatorUnion> interface_array = interface.array(mfi);
      Array4<Real const> moments_array = moments.const_array(mfi);
      const Box& bx = mfi.tilebox();

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // Compute cell volume fraction
        const double vfrac = moments_array(i, j, k, comp_vf);
        // Skip cell if vfrac ~0 or vfrac~1
        if (vfrac < IRL::global_constants::VF_LOW ||
            vfrac > IRL::global_constants::VF_HIGH) {
          const double distance = std::copysign(1.0, vfrac - 0.5);
          interface_array(i, j, k) = IRL::PlanarSeparator::fromOnePlane(
              IRL::Plane(IRL::Normal(0.0, 0.0, 0.0), distance));
          return;
        }

        // Fill in ELVIRA neighborhood
        IRL::ELVIRANeighborhood elvira_neighborhood;
        elvira_neighborhood.resize(27);
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
              cells_vfrac[neigh_id] = moments_array(ii, jj, kk, comp_vf);
              elvira_neighborhood.setMember(&cells[neigh_id],
                                            &cells_vfrac[neigh_id], ii - i,
                                            jj - j, kk - k);
            }
          }
        }
        // Reconstruct interface with ELVIRA and store result in interface_array
        interface_array(i, j, k) =
            IRL::reconstructionWithELVIRA3D(elvira_neighborhood);
      });
    }  // end mfi
    loop_timer.stop();

    // This copies the interface, update ghosts, and shifts datums across
    // periodic boundaries
    interface_with_ghost.LocalCopy(interface, 0, 0, interface.nComp(),
                                   interface.nGrowVect());
    interface_with_ghost.FillBoundaryWithPeriodicShift(geom);
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_RECONSTRUCTION_ELVIRA_H_
