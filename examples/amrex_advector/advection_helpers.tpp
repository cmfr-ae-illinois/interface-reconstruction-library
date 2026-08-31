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

template <class PtType>
inline void RestrictPtToBBox(PtType& pt, const PtType& lo, const PtType& hi) {
  for (int i = 0; i < 3; i++) {
    pt[i] = std::max(pt[i], lo[i]);
    pt[i] = std::min(pt[i], hi[i]);
  }
}

// general velocity accessor
inline IRL::Vec3<double> GetVelocity(
    const IRL::Pt& pt, const double time,
    const VelocityFieldType velocity_field_type, Array4<Real const> const& vx,
    Array4<Real const> const& vy, Array4<Real const> const& vz, const Box& bx,
    const Geometry& a_geom) {
  if (velocity_field_type == VelocityFieldType::Deformation) {
    return Deformation3D::get_velocity(pt[0], pt[1], pt[2], time);
  } else if (velocity_field_type == VelocityFieldType::Rotation) {
    return Rotation3D::get_velocity(pt[0], pt[1], pt[2], time);
  } else if (velocity_field_type == VelocityFieldType::Translation) {
    return Translation3D::get_velocity(pt[0], pt[1], pt[2], time);
  } else {
    return GetInterpolatedVelocity(pt, vx, vy, vz, bx, a_geom);
  }
}

inline IRL::Vec3<double> GetInterpolatedVelocity(const IRL::Pt& pt,
                                                 Array4<Real const> const& vx,
                                                 Array4<Real const> const& vy,
                                                 Array4<Real const> const& vz,
                                                 const Box& bx,
                                                 const Geometry& a_geom) {
  const auto& dx = a_geom.CellSizeArray();
  const auto& prob_lo = a_geom.ProbLoArray();
#ifndef NDEBUG
  const auto& lo = lbound(bx);
  const auto& hi = ubound(bx);
#endif

  IRL::Vec3<double> interpolated_velocity =
      IRL::Vec3<double>::fromScalarConstant(0.0);

  {  // X-component
    int i = static_cast<int>(Math::floor((pt[0] - prob_lo[0]) / dx[0]));
    int j = static_cast<int>(Math::floor((pt[1] - prob_lo[1]) / dx[1] - 0.5));
    int k = static_cast<int>(Math::floor((pt[2] - prob_lo[2]) / dx[2] - 0.5));
#ifndef NDEBUG
    if (i < lo.x || j < lo.y || k < lo.z || i > hi.x + 1 || j > hi.y ||
        k > hi.z) {
      std::ostringstream oss;
      oss << "Position X  (" << pt[0] << "," << pt[1] << "," << pt[2]
          << ") leads to index (" << i << "," << j << "," << k
          << ") outside of box (" << lo.x << "," << lo.y << "," << lo.z << ")x("
          << hi.x << "," << hi.y << "," << hi.z << ")";
      oss << " Problo = " << prob_lo[0] << ", " << prob_lo[1] << ", "
          << prob_lo[2];
      throw std::runtime_error(oss.str());
    }
#endif
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
#ifndef NDEBUG
    if (i < lo.x || j < lo.y || k < lo.z || i > hi.x || j > hi.y + 1 ||
        k > hi.z) {
      std::ostringstream oss;
      oss << "Position Y (" << pt[0] << "," << pt[1] << "," << pt[2]
          << ") leads to index (" << i << "," << j << "," << k
          << ") outside of box (" << lo.x << "," << lo.y << "," << lo.z << ")x("
          << hi.x << "," << hi.y << "," << hi.z << ")";
      oss << " Problo = " << prob_lo[0] << ", " << prob_lo[1] << ", "
          << prob_lo[2];
      throw std::runtime_error(oss.str());
    }
#endif
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
#ifndef NDEBUG
    if (i < lo.x || j < lo.y || k < lo.z || i > hi.x || j > hi.y ||
        k > hi.z + 1) {
      std::ostringstream oss;
      oss << "Position Z (" << pt[0] << "," << pt[1] << "," << pt[2]
          << ") leads to index (" << i << "," << j << "," << k
          << ") outside of box (" << lo.x << "," << lo.y << "," << lo.z << ")x("
          << hi.x << "," << hi.y << "," << hi.z << ")";
      oss << " Problo = " << prob_lo[0] << ", " << prob_lo[1] << ", "
          << prob_lo[2];
      throw std::runtime_error(oss.str());
    }
#endif
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

inline Eigen::Vector3d GetVelocity(const Eigen::Vector3d& pt, const double time,
                                   const VelocityFieldType velocity_field_type,
                                   Array4<Real const> const& vx,
                                   Array4<Real const> const& vy,
                                   Array4<Real const> const& vz, const Box& bx,
                                   const Geometry& a_geom) {
  const IRL::Pt pt_irl(pt[0], pt[1], pt[2]);
  const IRL::Vec3<double> velocity_irl =
      GetVelocity(pt_irl, time, velocity_field_type, vx, vy, vz, bx, a_geom);
  return Eigen::Vector3d({velocity_irl[0], velocity_irl[1], velocity_irl[2]});
}

// general velocity gradient accessor
inline Eigen::Matrix3d GetVelocityGradient(
    const Eigen::Vector3d& pt, const double time,
    const VelocityFieldType velocity_field_type, Array4<Real const> const& vx,
    Array4<Real const> const& vy, Array4<Real const> const& vz, const Box& bx,
    const Geometry& a_geom) {
  if (velocity_field_type == VelocityFieldType::Deformation) {
    return Deformation3D::get_velocity_gradient(pt[0], pt[1], pt[2], time);
  } else if (velocity_field_type == VelocityFieldType::Rotation) {
    return Rotation3D::get_velocity_gradient(pt[0], pt[1], pt[2], time);
  } else if (velocity_field_type == VelocityFieldType::Translation) {
    return Translation3D::get_velocity_gradient(pt[0], pt[1], pt[2], time);
  } else {
    return GetInterpolatedVelocityGradient(pt, vx, vy, vz, bx, a_geom);
  }
}

inline Eigen::Matrix3d GetInterpolatedVelocityGradient(
    const Eigen::Vector3d& pt, Array4<Real const> const& vx,
    Array4<Real const> const& vy, Array4<Real const> const& vz, const Box& bx,
    const Geometry& a_geom) {
  using Array3_1D = std::array<double, 3>;
  using Array2_3D = std::array<std::array<std::array<double, 2>, 2>, 2>;

  const auto& dx = a_geom.CellSizeArray();
  const Array3_1D invdx = {1.0 / dx[0], 1.0 / dx[1], 1.0 / dx[2]};
  const auto& prob_lo = a_geom.ProbLoArray();
  const auto& lo = lbound(bx);
  const auto& hi = ubound(bx);

  Eigen::Matrix3d interpolated_grad = Eigen::Matrix3d::Zero();

  auto trilinear_interpolation = [](double tx, double ty, double tz,
                                    Array2_3D& data) -> double {
    return tz * (ty * (tx * data[1][1][1] + (1.0 - tx) * data[0][1][1]) +
                 (1.0 - ty) *
                     (tx * data[1][0][1] + (1.0 - tx) * data[0][0][1])) +
           (1.0 - tz) *
               (ty * (tx * data[1][1][0] + (1.0 - tx) * data[0][1][0]) +
                (1.0 - ty) * (tx * data[1][0][0] + (1.0 - tx) * data[0][0][0]));
  };

  // Start with diagonal gradient components
  {
    int i = static_cast<int>(Math::floor((pt[0] - prob_lo[0]) / dx[0] - 0.5));
    int j = static_cast<int>(Math::floor((pt[1] - prob_lo[1]) / dx[1] - 0.5));
    int k = static_cast<int>(Math::floor((pt[2] - prob_lo[2]) / dx[2] - 0.5));
    // Cell lo face positions
    Real xlo = prob_lo[0] + i * dx[0] + 0.5 * dx[0];
    Real ylo = prob_lo[1] + j * dx[1] + 0.5 * dx[1];
    Real zlo = prob_lo[2] + k * dx[2] + 0.5 * dx[2];
    // Normalized position within cell [0,1]
    Real tx = (pt[0] - xlo) / dx[0];
    Real ty = (pt[1] - ylo) / dx[1];
    Real tz = (pt[2] - zlo) / dx[2];
    // Compute gradient at interpolation domain corners
    Array2_3D GradUx, GradVy, GradWz;
    for (int ii = i; ii <= i + 1; ii++) {
      for (int jj = j; jj <= j + 1; jj++) {
        for (int kk = k; kk <= k + 1; kk++) {
          GradUx[ii - i][jj - j][kk - k] =
              (vx(ii + 1, jj, kk) - vx(ii, jj, kk)) * invdx[0];
          GradVy[ii - i][jj - j][kk - k] =
              (vy(ii, jj + 1, kk) - vy(ii, jj, kk)) * invdx[1];
          GradWz[ii - i][jj - j][kk - k] =
              (vz(ii, jj, kk + 1) - vz(ii, jj, kk)) * invdx[2];
        }
      }
    }
    // Trilinear interpolation
    interpolated_grad(0, 0) = trilinear_interpolation(tx, ty, tz, GradUx);
    interpolated_grad(1, 1) = trilinear_interpolation(tx, ty, tz, GradVy);
    interpolated_grad(2, 2) = trilinear_interpolation(tx, ty, tz, GradWz);
  }

  // GradUy, GradVx
  {
    int i = static_cast<int>(Math::floor((pt[0] - prob_lo[0]) / dx[0]));
    int j = static_cast<int>(Math::floor((pt[1] - prob_lo[1]) / dx[1]));
    int k = static_cast<int>(Math::floor((pt[2] - prob_lo[2]) / dx[2] - 0.5));
    // Cell lo face positions
    Real xlo = prob_lo[0] + i * dx[0];
    Real ylo = prob_lo[1] + j * dx[1];
    Real zlo = prob_lo[2] + k * dx[2] + 0.5 * dx[2];
    // Normalized position within cell [0,1]
    Real tx = (pt[0] - xlo) / dx[0];
    Real ty = (pt[1] - ylo) / dx[1];
    Real tz = (pt[2] - zlo) / dx[2];
    // Compute gradient at interpolation domain corners
    Array2_3D GradUy, GradVx;
    for (int ii = i; ii <= i + 1; ii++) {
      for (int jj = j; jj <= j + 1; jj++) {
        for (int kk = k; kk <= k + 1; kk++) {
          GradUy[ii - i][jj - j][kk - k] =
              (vx(ii, jj, kk) - vx(ii, jj - 1, kk)) * invdx[1];
          GradVx[ii - i][jj - j][kk - k] =
              (vy(ii, jj, kk) - vy(ii - 1, jj, kk)) * invdx[0];
        }
      }
    }
    // Trilinear interpolation
    interpolated_grad(0, 1) = trilinear_interpolation(tx, ty, tz, GradUy);
    interpolated_grad(1, 0) = trilinear_interpolation(tx, ty, tz, GradVx);
  }

  // GradUz, GradWx
  {
    int i = static_cast<int>(Math::floor((pt[0] - prob_lo[0]) / dx[0]));
    int j = static_cast<int>(Math::floor((pt[1] - prob_lo[1]) / dx[1] - 0.5));
    int k = static_cast<int>(Math::floor((pt[2] - prob_lo[2]) / dx[2]));
    // Cell lo face positions
    Real xlo = prob_lo[0] + i * dx[0];
    Real ylo = prob_lo[1] + j * dx[1] + 0.5 * dx[1];
    Real zlo = prob_lo[2] + k * dx[2];
    // Normalized position within cell [0,1]
    Real tx = (pt[0] - xlo) / dx[0];
    Real ty = (pt[1] - ylo) / dx[1];
    Real tz = (pt[2] - zlo) / dx[2];
    // Compute gradient at interpolation domain corners
    Array2_3D GradUz, GradWx;
    for (int ii = i; ii <= i + 1; ii++) {
      for (int jj = j; jj <= j + 1; jj++) {
        for (int kk = k; kk <= k + 1; kk++) {
          GradUz[ii - i][jj - j][kk - k] =
              (vx(ii, jj, kk) - vx(ii, jj, kk - 1)) * invdx[2];
          GradWx[ii - i][jj - j][kk - k] =
              (vz(ii, jj, kk) - vz(ii - 1, jj, kk)) * invdx[0];
        }
      }
    }
    // Trilinear interpolation
    interpolated_grad(0, 2) = trilinear_interpolation(tx, ty, tz, GradUz);
    interpolated_grad(2, 0) = trilinear_interpolation(tx, ty, tz, GradWx);
  }

  // GradVz, GradWy
  {
    int i = static_cast<int>(Math::floor((pt[0] - prob_lo[0]) / dx[0] - 0.5));
    int j = static_cast<int>(Math::floor((pt[1] - prob_lo[1]) / dx[1]));
    int k = static_cast<int>(Math::floor((pt[2] - prob_lo[2]) / dx[2]));
    // Cell lo face positions
    Real xlo = prob_lo[0] + i * dx[0] + 0.5 * dx[0];
    Real ylo = prob_lo[1] + j * dx[1];
    Real zlo = prob_lo[2] + k * dx[2];
    // Normalized position within cell [0,1]
    Real tx = (pt[0] - xlo) / dx[0];
    Real ty = (pt[1] - ylo) / dx[1];
    Real tz = (pt[2] - zlo) / dx[2];
    // Compute gradient at interpolation domain corners
    Array2_3D GradVz, GradWy;
    for (int ii = i; ii <= i + 1; ii++) {
      for (int jj = j; jj <= j + 1; jj++) {
        for (int kk = k; kk <= k + 1; kk++) {
          GradVz[ii - i][jj - j][kk - k] =
              (vy(ii, jj, kk) - vy(ii, jj, kk - 1)) * invdx[2];
          GradWy[ii - i][jj - j][kk - k] =
              (vz(ii, jj, kk) - vz(ii, jj - 1, kk)) * invdx[1];
        }
      }
    }
    // Trilinear interpolation
    interpolated_grad(1, 2) = trilinear_interpolation(tx, ty, tz, GradVz);
    interpolated_grad(2, 1) = trilinear_interpolation(tx, ty, tz, GradWy);
  }

  return interpolated_grad;
}

IRL::Pt ProjectVertex(const IRL::Pt& pt, const double dt, const double time,
                      const VelocityFieldType velocity_field_type,
                      Array4<Real const> const& vx,
                      Array4<Real const> const& vy,
                      Array4<Real const> const& vz, const Box& bx,
                      const Geometry& a_geom) {
  const auto v1 =
      GetVelocity(pt, time, velocity_field_type, vx, vy, vz, bx, a_geom);
  const auto v2 =
      GetVelocity(pt + IRL::Pt::fromVec3(0.5 * dt * v1), time + 0.5 * dt,
                  velocity_field_type, vx, vy, vz, bx, a_geom);
  const auto v3 =
      GetVelocity(pt + IRL::Pt::fromVec3(0.5 * dt * v2), time + 0.5 * dt,
                  velocity_field_type, vx, vy, vz, bx, a_geom);
  const auto v4 = GetVelocity(pt + IRL::Pt::fromVec3(dt * v3), time + dt,
                              velocity_field_type, vx, vy, vz, bx, a_geom);
  return pt + IRL::Pt::fromVec3(dt * (v1 + 2.0 * v2 + 2.0 * v3 + v4) / 6.0);
}

#endif  // EXAMPLES_AMREX_ADVECTOR_ADVECTION_HELPERS_TPP_
