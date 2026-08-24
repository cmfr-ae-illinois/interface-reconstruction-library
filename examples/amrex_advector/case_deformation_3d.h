// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_DEFORMATION_3D_H_
#define EXAMPLES_AMREX_ADVECTOR_DEFORMATION_3D_H_

#include <AMReX_AmrCore.H>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/variant_reconstruction/separator_union.h"

using namespace amrex;

struct Deformation3D {
  static Real get_max_vel(void) { return 2.0; }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void initialize_case(
      Box const& bx, Array4<Real> const& moments,
      Array4<IRL::SeparatorUnion> const& interface,
      GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
      GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx, bool const& transport_m1,
      bool const& transport_m2) {
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    const IRL::Pt sphere_center(0.35, 0.35, 0.35);
    const double sphere_radius = 0.15;
    const double diag =
        std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    const double vol = dx[0] * dx[1] * dx[2];

    for (int k = lo.z; k <= hi.z; ++k) {
      const double z = problo[2] + k * dx[2];
      for (int j = lo.y; j <= hi.y; ++j) {
        const double y = problo[1] + j * dx[1];
        for (int i = lo.x; i <= hi.x; ++i) {
          const double x = problo[0] + i * dx[0];
          const IRL::Pt lower_cell_pt(x, y, z);
          const IRL::Pt upper_cell_pt(x + dx[0], y + dx[1], z + dx[2]);
          const IRL::Pt mid = 0.5 * (lower_cell_pt + upper_cell_pt);
          IRL::Pt disp = mid - sphere_center;
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
          const auto cell_moments =
              IRL::getVolumeMoments<IRL::SeparatedMoments<IRL::VolumeMoments>>(
                  cell, interface(i, j, k));
          moments(i, j, k, 0) = cell_moments[0].volume();
          if (transport_m1) {
            moments(i, j, k, 1) = cell_moments[0].centroid()[0];
            moments(i, j, k, 2) = cell_moments[0].centroid()[1];
            moments(i, j, k, 3) = cell_moments[0].centroid()[2];
            moments(i, j, k, 4) = cell_moments[1].centroid()[0];
            moments(i, j, k, 5) = cell_moments[1].centroid()[1];
            moments(i, j, k, 6) = cell_moments[1].centroid()[2];
          }
          if (transport_m2) {
            const auto cell_general_moments = IRL::getVolumeMoments<
                IRL::SeparatedMoments<IRL::GeneralMoments3D<2>>>(
                cell, interface(i, j, k));
            moments(i, j, k, 7) = cell_general_moments[0][4];
            moments(i, j, k, 8) = cell_general_moments[0][5];
            moments(i, j, k, 9) = cell_general_moments[0][6];
            moments(i, j, k, 10) = cell_general_moments[0][7];
            moments(i, j, k, 11) = cell_general_moments[0][8];
            moments(i, j, k, 12) = cell_general_moments[0][9];
            moments(i, j, k, 13) = cell_general_moments[1][4];
            moments(i, j, k, 14) = cell_general_moments[1][5];
            moments(i, j, k, 15) = cell_general_moments[1][6];
            moments(i, j, k, 16) = cell_general_moments[1][7];
            moments(i, j, k, 17) = cell_general_moments[1][8];
            moments(i, j, k, 18) = cell_general_moments[1][9];
          }
        }
      }
    }
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real get_face_velocity_x(
      const Real x, const Real y, const Real z, const Real time) {
    Real sinpix = std::sin(M_PI * x);
    Real sin2pix = std::sin(2.0 * M_PI * x);
    Real sinpiy = std::sin(M_PI * y);
    Real sin2piy = std::sin(2.0 * M_PI * y);
    Real sinpiz = std::sin(M_PI * z);
    Real sin2piz = std::sin(2.0 * M_PI * z);
    Real cospit = std::cos(M_PI * time / 3.0);
    return 2.0 * sinpix * sinpix * sin2piy * sin2piz * cospit;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real get_face_velocity_y(
      const Real x, const Real y, const Real z, const Real time) {
    Real sinpix = std::sin(M_PI * x);
    Real sin2pix = std::sin(2.0 * M_PI * x);
    Real sinpiy = std::sin(M_PI * y);
    Real sin2piy = std::sin(2.0 * M_PI * y);
    Real sinpiz = std::sin(M_PI * z);
    Real sin2piz = std::sin(2.0 * M_PI * z);
    Real cospit = std::cos(M_PI * time / 3.0);
    return -sinpiy * sinpiy * sin2pix * sin2piz * cospit;
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real get_face_velocity_z(
      const Real x, const Real y, const Real z, const Real time) {
    Real sinpix = std::sin(M_PI * x);
    Real sin2pix = std::sin(2.0 * M_PI * x);
    Real sinpiy = std::sin(M_PI * y);
    Real sin2piy = std::sin(2.0 * M_PI * y);
    Real sinpiz = std::sin(M_PI * z);
    Real sin2piz = std::sin(2.0 * M_PI * z);
    Real cospit = std::cos(M_PI * time / 3.0);
    return -sinpiz * sinpiz * sin2pix * sin2piy * cospit;
  }

  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE IRL::Vec3<double>
  get_velocity(const Real x, const Real y, const Real z, const Real time) {
    const Real u = get_face_velocity_x(x, y, z, time);
    const Real v = get_face_velocity_y(x, y, z, time);
    const Real w = get_face_velocity_z(x, y, z, time);

    return IRL::Vec3<double>(u, v, w);
  }

  static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Eigen::Matrix3d
  get_velocity_gradient(const Real x, const Real y, const Real z,
                        const Real time) {
    const double sinpix = std::sin(M_PI * x), sinpiy = std::sin(M_PI * y),
                 sinpiz = std::sin(M_PI * z);
    const double cospix = std::cos(M_PI * x), cospiy = std::cos(M_PI * y),
                 cospiz = std::cos(M_PI * z);
    const double sin2pix = std::sin(2.0 * M_PI * x),
                 sin2piy = std::sin(2.0 * M_PI * y),
                 sin2piz = std::sin(2.0 * M_PI * z);
    const double cos2pix = std::cos(2.0 * M_PI * x),
                 cos2piy = std::cos(2.0 * M_PI * y),
                 cos2piz = std::cos(2.0 * M_PI * z);
    const double cospit = std::cos(M_PI * time / 3.0);
    return Eigen::Matrix3d(
        {{4.0 * M_PI * cospix * sinpix * sin2piy * sin2piz * cospit,
          4.0 * M_PI * sinpix * sinpix * cos2piy * sin2piz * cospit,
          4.0 * M_PI * sinpix * sinpix * sin2piy * cos2piz * cospit},
         {-2.0 * M_PI * sinpiy * sinpiy * cos2pix * sin2piz * cospit,
          -2.0 * M_PI * cospiy * sinpiy * sin2pix * sin2piz * cospit,
          -2.0 * M_PI * sinpiy * sinpiy * sin2pix * cos2piz * cospit},
         {-2.0 * M_PI * sinpiz * sinpiz * cos2pix * sin2piy * cospit,
          -2.0 * M_PI * sinpiz * sinpiz * sin2pix * cos2piy * cospit,
          -2.0 * M_PI * cospiz * sinpiz * sin2pix * sin2piy * cospit}});
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_DEFORMATION_3D_H_
