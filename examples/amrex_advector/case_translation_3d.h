// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_TRANSLATION_3D_H_
#define EXAMPLES_AMREX_ADVECTOR_TRANSLATION_3D_H_

#include <AMReX_AmrCore.H>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/variant_reconstruction/separator_union.h"

using namespace amrex;

struct Translation3D {
  static Real get_max_vel(void) { return 3.0; }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void initialize_case(
      Box const& bx, Array4<Real> const& moments,
      Array4<IRL::SeparatorUnion> const& interface,
      GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
      GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx) {
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    const IRL::Pt sphere_center(0.5, 0.5, 0.5);
    const double sphere_radius = 0.25;
    const double diag =
        std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

    for (int k = lo.z; k <= hi.z; ++k) {
      const double z = problo[2] + k * dx[2];
      for (int j = lo.y; j <= hi.y; ++j) {
        const double y = problo[1] + j * dx[1];
        for (int i = lo.x; i <= hi.x; ++i) {
          const double x = problo[0] + i * dx[0];
          const IRL::Pt lower_cell_pt(x, y, z);
          const IRL::Pt upper_cell_pt(x + dx[0], y + dx[1], z + dx[2]);
          const IRL::Pt mid_pt = 0.5 * (lower_cell_pt + upper_cell_pt);
          IRL::Pt disp = mid_pt - sphere_center;
          const auto mag = IRL::magnitude(disp);
          if (mag < sphere_radius - 0.5 * diag) {
            interface(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
          } else if (mag > sphere_radius + 0.5 * diag) {
            interface(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
          } else {
            auto sphere_normal = IRL::Normal::fromPt(disp);
            sphere_normal.normalize();
            const double curvature = 1.0 / sphere_radius;
            const auto frame = IRL::ReferenceFrame::fromNormal(sphere_normal);
            interface(i, j, k) =
                IRL::Paraboloid(sphere_center + sphere_radius * sphere_normal,
                                frame, 0.5 * curvature, 0.5 * curvature);
          }
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              lower_cell_pt, upper_cell_pt);
          const auto cell_moments = IRL::getVolumeMoments<IRL::VolumeMoments>(
              cell, interface(i, j, k));
          moments(i, j, k, 0) = cell_moments.volume();
          moments(i, j, k, 1) = cell_moments.centroid()[0];
          moments(i, j, k, 2) = cell_moments.centroid()[1];
          moments(i, j, k, 3) = cell_moments.centroid()[2];
        }
      }
    }
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real
  get_face_velocity_x(const Real x = 0.0, const Real y = 0.0,
                      const Real z = 0.0, const Real time = 0.0) {
    return 1.0;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real
  get_face_velocity_y(const Real x = 0.0, const Real y = 0.0,
                      const Real z = 0.0, const Real time = 0.0) {
    return 2.0;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real
  get_face_velocity_z(const Real x = 0.0, const Real y = 0.0,
                      const Real z = 0.0, const Real time = 0.0) {
    return 3.0;
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_TRANSLATION_3D_H_
