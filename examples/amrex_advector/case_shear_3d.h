// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_SHEAR_3D_H_
#define EXAMPLES_AMREX_ADVECTOR_SHEAR_3D_H_

#include <AMReX_AmrCore.H>

#include <cmath>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/variant_reconstruction/separator_union.h"

using namespace amrex;

struct Shear3D {
  static constexpr Real Umax = 1.0;
  static constexpr Real R = 0.5;
  static constexpr Real T = 3.0;

  static Real get_max_vel(void) { return 1.03443; }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE void initialize_case(
      Box const& bx, Array4<Real> const& moments,
      Array4<IRL::SeparatorUnion> const& interface,
      GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
      GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx, bool const& transport_m1,
      bool const& transport_m2) {
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    const IRL::Pt sphere_center(0.5, 0.75, 0.25);
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

          moments(i, j, k, comp_vf) = cell_moments[0].volume() / vol;
          moments(i, j, k, comp_m0) = cell_moments[0].volume();

          if (transport_m1) {
            moments(i, j, k, comp_m1_l) = cell_moments[0].centroid()[0];
            moments(i, j, k, comp_m1_l + 1) = cell_moments[0].centroid()[1];
            moments(i, j, k, comp_m1_l + 2) = cell_moments[0].centroid()[2];
            moments(i, j, k, comp_m1_g) = cell_moments[1].centroid()[0];
            moments(i, j, k, comp_m1_g + 1) = cell_moments[1].centroid()[1];
            moments(i, j, k, comp_m1_g + 2) = cell_moments[1].centroid()[2];
          }

          if (transport_m2) {
            const auto cell_general_moments = IRL::getVolumeMoments<
                IRL::SeparatedMoments<IRL::GeneralMoments3D<2>>>(
                cell, interface(i, j, k));
            moments(i, j, k, comp_m2_l) = cell_general_moments[0][4];
            moments(i, j, k, comp_m2_l + 1) = cell_general_moments[0][5];
            moments(i, j, k, comp_m2_l + 2) = cell_general_moments[0][6];
            moments(i, j, k, comp_m2_l + 3) = cell_general_moments[0][7];
            moments(i, j, k, comp_m2_l + 4) = cell_general_moments[0][8];
            moments(i, j, k, comp_m2_l + 5) = cell_general_moments[0][9];
            moments(i, j, k, comp_m2_g) = cell_general_moments[1][4];
            moments(i, j, k, comp_m2_g + 1) = cell_general_moments[1][5];
            moments(i, j, k, comp_m2_g + 2) = cell_general_moments[1][6];
            moments(i, j, k, comp_m2_g + 3) = cell_general_moments[1][7];
            moments(i, j, k, comp_m2_g + 4) = cell_general_moments[1][8];
            moments(i, j, k, comp_m2_g + 5) = cell_general_moments[1][9];
          }
        }
      }
    }
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real get_face_velocity_x(
      const Real x, const Real y, const Real z, const Real time) {
    const Real sinpix = std::sin(M_PI * x);
    return sinpix * sinpix * std::sin(2.0 * M_PI * y) *
           std::cos(M_PI * time / T);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real get_face_velocity_y(
      const Real x, const Real y, const Real z, const Real time) {
    const Real sinpiy = std::sin(M_PI * y);
    return -std::sin(2.0 * M_PI * x) * sinpiy * sinpiy *
           std::cos(M_PI * time / T);
  }

  static AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real get_face_velocity_z(
      const Real x, const Real y, const Real z, const Real time) {
    const Real rx = x - 0.5;
    const Real ry = y - 0.5;
    const Real r = std::sqrt(rx * rx + ry * ry);
    const Real profile = 1.0 - r / R;
    return Umax * profile * profile * std::cos(M_PI * time / T);
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
    const double sinpix = std::sin(M_PI * x);
    const double sinpiy = std::sin(M_PI * y);
    const double sin2pix = std::sin(2.0 * M_PI * x);
    const double sin2piy = std::sin(2.0 * M_PI * y);
    const double cos2pix = std::cos(2.0 * M_PI * x);
    const double cos2piy = std::cos(2.0 * M_PI * y);
    const double cospit = std::cos(M_PI * time / T);

    const double dudx = M_PI * sin2pix * sin2piy * cospit;
    const double dudy = 2.0 * M_PI * sinpix * sinpix * cos2piy * cospit;
    const double dudz = 0.0;

    const double dvdx = -2.0 * M_PI * cos2pix * sinpiy * sinpiy * cospit;
    const double dvdy = -M_PI * sin2pix * sin2piy * cospit;
    const double dvdz = 0.0;

    const double rx = x - 0.5;
    const double ry = y - 0.75;
    const double r = std::sqrt(rx * rx + ry * ry);

    double dwdx = 0.0;
    double dwdy = 0.0;
    if (r > 0.0) {
      const double factor = -2.0 * Umax * (1.0 - r / R) * cospit / (R * r);
      dwdx = factor * rx;
      dwdy = factor * ry;
    }
    // At r == 0 the classical radial gradient is undefined.
    const double dwdz = 0.0;

    return Eigen::Matrix3d(
        {{dudx, dudy, dudz}, {dvdx, dvdy, dvdz}, {dwdx, dwdy, dwdz}});
  }
};

#endif  // EXAMPLES_AMREX_ADVECTOR_SHEAR_3D_H_