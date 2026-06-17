// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_AMREX_ADVECTOR_ADVECTION_HELPERS_TPP_
#define EXAMPLES_AMREX_ADVECTOR_ADVECTION_HELPERS_TPP_

#include "irl/amrex/sepunion_multifab.h"

using namespace amrex;

inline IRL::Vec3<double> GetVelocity(const IRL::Pt& pt,
                                     Array4<Real const> const& vx,
                                     Array4<Real const> const& vy,
                                     Array4<Real const> const& vz,
                                     const Box& bx, const Geometry& a_geom) {
  const auto& dx = a_geom.CellSizeArray();
  const auto& prob_lo = a_geom.ProbLoArray();
  const auto& lo = lbound(bx);
  const auto& hi = ubound(bx);

  IRL::Vec3<double> interpolated_velocity =
      IRL::Vec3<double>::fromScalarConstant(0.0);

  {  // X-component
    int i = static_cast<int>(Math::floor((pt[0] - prob_lo[0]) / dx[0]));
    int j = static_cast<int>(Math::floor((pt[1] - prob_lo[1]) / dx[1] - 0.5));
    int k = static_cast<int>(Math::floor((pt[2] - prob_lo[2]) / dx[2] - 0.5));
    if (!bx.contains(i, j, k)) {
      std::ostringstream oss;
      oss << "Position (" << pt[0] << "," << pt[1] << "," << pt[2]
          << ") leads to index (" << i << "," << j << "," << k
          << ") outside of box (" << lo.x << "," << lo.y << "," << lo.z << ")x("
          << hi.x << "," << hi.y << "," << hi.z << ")";
      throw std::runtime_error(oss.str());
    }
    // Cell lo face positions
    Real xlo = prob_lo[0] + i * dx[0];
    Real ylo = prob_lo[1] + j * dx[1] + 0.5 * dx[1];
    Real zlo = prob_lo[2] + k * dx[2] + 0.5 * dx[2];
    // Normalized position within cell [0,1]
    Real tx = (pt[0] - xlo) / dx[0];
    Real ty = (pt[1] - ylo) / dx[1];
    Real tz = (pt[2] - zlo) / dx[2];
    // Trilinear interpolation
    interpolated_velocity[0] =
        tz * (ty * (tx * vx(i + 1, j + 1, k + 1) +
                    (1.0 - tx) * vx(i, j + 1, k + 1)) +
              (1.0 - ty) *
                  (tx * vx(i + 1, j, k + 1) + (1.0 - tx) * vx(i, j, k + 1))) +
        (1.0 - tz) *
            (ty * (tx * vx(i + 1, j + 1, k) + (1.0 - tx) * vx(i, j + 1, k)) +
             (1.0 - ty) * (tx * vx(i + 1, j, k) + (1.0 - tx) * vx(i, j, k)));
  }
  {  // Y-component
    int i = static_cast<int>(Math::floor((pt[0] - prob_lo[0]) / dx[0] - 0.5));
    int j = static_cast<int>(Math::floor((pt[1] - prob_lo[1]) / dx[1]));
    int k = static_cast<int>(Math::floor((pt[2] - prob_lo[2]) / dx[2] - 0.5));
    if (!bx.contains(i, j, k)) {
      std::ostringstream oss;
      oss << "Position (" << pt[0] << "," << pt[1] << "," << pt[2]
          << ") leads to index (" << i << "," << j << "," << k
          << ") outside of box (" << lo.x << "," << lo.y << "," << lo.z << ")x("
          << hi.x << "," << hi.y << "," << hi.z << ")";
      throw std::runtime_error(oss.str());
    }
    // Cell lo face positions
    Real xlo = prob_lo[0] + i * dx[0] + 0.5 * dx[0];
    Real ylo = prob_lo[1] + j * dx[1];
    Real zlo = prob_lo[2] + k * dx[2] + 0.5 * dx[2];
    // Normalized position within cell [0,1]
    Real tx = (pt[0] - xlo) / dx[0];
    Real ty = (pt[1] - ylo) / dx[1];
    Real tz = (pt[2] - zlo) / dx[2];
    // Trilinear interpolation
    interpolated_velocity[1] =
        tz * (ty * (tx * vy(i + 1, j + 1, k + 1) +
                    (1.0 - tx) * vy(i, j + 1, k + 1)) +
              (1.0 - ty) *
                  (tx * vy(i + 1, j, k + 1) + (1.0 - tx) * vy(i, j, k + 1))) +
        (1.0 - tz) *
            (ty * (tx * vy(i + 1, j + 1, k) + (1.0 - tx) * vy(i, j + 1, k)) +
             (1.0 - ty) * (tx * vy(i + 1, j, k) + (1.0 - tx) * vy(i, j, k)));
  }
  {  // Z-component
    int i = static_cast<int>(Math::floor((pt[0] - prob_lo[0]) / dx[0] - 0.5));
    int j = static_cast<int>(Math::floor((pt[1] - prob_lo[1]) / dx[1] - 0.5));
    int k = static_cast<int>(Math::floor((pt[2] - prob_lo[2]) / dx[2]));
    if (!bx.contains(i, j, k)) {
      std::ostringstream oss;
      oss << "Position (" << pt[0] << "," << pt[1] << "," << pt[2]
          << ") leads to index (" << i << "," << j << "," << k
          << ") outside of box (" << lo.x << "," << lo.y << "," << lo.z << ")x("
          << hi.x << "," << hi.y << "," << hi.z << ")";
      throw std::runtime_error(oss.str());
    }
    // Cell lo face positions
    Real xlo = prob_lo[0] + i * dx[0] + 0.5 * dx[0];
    Real ylo = prob_lo[1] + j * dx[1] + 0.5 * dx[1];
    Real zlo = prob_lo[2] + k * dx[2];
    // Normalized position within cell [0,1]
    Real tx = (pt[0] - xlo) / dx[0];
    Real ty = (pt[1] - ylo) / dx[1];
    Real tz = (pt[2] - zlo) / dx[2];
    // Trilinear interpolation
    interpolated_velocity[2] =
        tz * (ty * (tx * vz(i + 1, j + 1, k + 1) +
                    (1.0 - tx) * vz(i, j + 1, k + 1)) +
              (1.0 - ty) *
                  (tx * vz(i + 1, j, k + 1) + (1.0 - tx) * vz(i, j, k + 1))) +
        (1.0 - tz) *
            (ty * (tx * vz(i + 1, j + 1, k) + (1.0 - tx) * vz(i, j + 1, k)) +
             (1.0 - ty) * (tx * vz(i + 1, j, k) + (1.0 - tx) * vz(i, j, k)));
  }
  return interpolated_velocity;
}

inline IRL::Pt ProjectVertex(const IRL::Pt& pt, const double dt,
                             Array4<Real const> const& vx,
                             Array4<Real const> const& vy,
                             Array4<Real const> const& vz, const Box& bx,
                             const Geometry& a_geom) {
  using Pt = IRL::Pt;
  auto v1 = GetVelocity(pt, vx, vy, vz, bx, a_geom);
  auto v2 =
      GetVelocity(pt + Pt::fromVec3(0.5 * dt * v1), vx, vy, vz, bx, a_geom);
  auto v3 =
      GetVelocity(pt + Pt::fromVec3(0.5 * dt * v2), vx, vy, vz, bx, a_geom);
  auto v4 = GetVelocity(pt + Pt::fromVec3(dt * v3), vx, vy, vz, bx, a_geom);
  return pt + Pt::fromVec3(dt * (v1 + 2.0 * v2 + 2.0 * v3 + v4) / 6.0);
}

#endif  // EXAMPLES_AMREX_ADVECTOR_ADVECTION_HELPERS_TPP_
