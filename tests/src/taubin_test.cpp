// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "gtest/gtest.h"

#include <Eigen/Dense>
#include <algorithm>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

#include "irl/geometry/general/pt.h"

namespace {

using namespace IRL;

TEST(TaubinTest, MatrixTerms) {
  // line segment
  Eigen::Vector2d p0(0.0, 0.0);
  Eigen::Vector2d p1(1.0, 1.0);

  // ==== point sampling formulation =====
  const int n_points = 1000000;
  std::vector<Eigen::Vector2d> points(n_points);
  for (int i = 0; i < n_points; ++i) {
    double t = static_cast<double>(i) / (static_cast<double>(n_points) - 1.0);
    points[i] = (1.0 - t) * p0 + t * p1;
    // std::cout << "Point " << i << ": " << points[i].transpose() << std::endl;
  }
  // moment matrix
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  for (int i = 0; i < n_points; ++i) {
    double xi = points[i].x(), yi = points[i].y();
    double zi = xi * xi + yi * yi;
    Eigen::Vector4d u;
    u << zi, xi, yi, 1.0;
    M += u * u.transpose();
  }

  // constraint matrix
  Eigen::Matrix4d C = Eigen::Matrix4d::Zero();
  C(0, 0) = 4.0 * M(0, 3);
  C(0, 1) = 2.0 * M(1, 3);
  C(0, 2) = 2.0 * M(2, 3);
  C(1, 0) = C(0, 1);
  C(1, 1) = static_cast<double>(n_points);
  C(2, 0) = C(0, 2);
  C(2, 2) = static_cast<double>(n_points);

  M /= static_cast<double>(n_points);
  std::cout << "Moment matrix M (point sampling):\n" << M << std::endl;

  C /= static_cast<double>(n_points);
  std::cout << "Constraint matrix C (point sampling):\n" << C << std::endl;

  // ==== line segment formulation =====
  const double x0 = p0.x();
  const double y0 = p0.y();
  const double dx = p1.x() - p0.x();
  const double dy = p1.y() - p0.y();

  // moment matrix
  std::vector<double> terms(10, 0.0);

  // Mzz
  terms[0] = (dx * dx * dx * dx) / 5. + (2. * (dx * dx) * (dy * dy)) / 5. +
             (dy * dy * dy * dy) / 5. + dx * dx * dx * x0 +
             dx * (dy * dy) * x0 + 2. * (dx * dx) * (x0 * x0) +
             (2. * (dy * dy) * (x0 * x0)) / 3. + 2. * dx * (x0 * x0 * x0) +
             x0 * x0 * x0 * x0 + dx * dx * dy * y0 + dy * dy * dy * y0 +
             (8. * dx * dy * x0 * y0) / 3. + 2. * dy * (x0 * x0) * y0 +
             (2. * (dx * dx) * (y0 * y0)) / 3. + 2. * (dy * dy) * (y0 * y0) +
             2. * dx * x0 * (y0 * y0) + 2. * (x0 * x0) * (y0 * y0) +
             2. * dy * (y0 * y0 * y0) + y0 * y0 * y0 * y0;

  // Mxz
  terms[1] = (dx * dx * dx) / 4. + (dx * (dy * dy)) / 4. + dx * dx * x0 +
             (dy * dy * x0) / 3. + (3. * dx * (x0 * x0)) / 2. + x0 * x0 * x0 +
             (2. * dx * dy * y0) / 3. + dy * x0 * y0 + (dx * (y0 * y0)) / 2. +
             x0 * (y0 * y0);

  // Myz
  terms[2] = (dx * dx * dy) / 4. + (dy * dy * dy) / 4. +
             (2. * dx * dy * x0) / 3. + (dy * (x0 * x0)) / 2. +
             (dx * dx * y0) / 3. + dy * dy * y0 + dx * x0 * y0 + x0 * x0 * y0 +
             (3. * dy * (y0 * y0)) / 2. + y0 * y0 * y0;

  // Mz
  terms[3] =
      (dx * dx) / 3. + (dy * dy) / 3. + dx * x0 + x0 * x0 + dy * y0 + y0 * y0;

  // Mxx
  terms[4] = (dx * dx) / 3. + dx * x0 + x0 * x0;

  // Mxy
  terms[5] = (dx * dy) / 3. + (dy * x0) / 2. + (dx * y0) / 2. + x0 * y0;

  // Mx
  terms[6] = dx / 2. + x0;

  // Myy
  terms[7] = (dy * dy) / 3. + dy * y0 + y0 * y0;

  // My
  terms[8] = dy / 2. + y0;

  // Le
  terms[9] = 1.0;

  double Mzz = terms[0], Mxz = terms[1], Myz = terms[2], Mz = terms[3],
         Mxx = terms[4], Mxy = terms[5], Mx = terms[6], Myy = terms[7],
         My = terms[8];

  Eigen::Matrix4d M_line = Eigen::Matrix4d::Zero();
  M_line(0, 0) = Mzz;
  M_line(0, 1) = Mxz;
  M_line(0, 2) = Myz;
  M_line(0, 3) = Mz;
  M_line(1, 0) = Mxz;
  M_line(1, 1) = Mxx;
  M_line(1, 2) = Mxy;
  M_line(1, 3) = Mx;
  M_line(2, 0) = Myz;
  M_line(2, 1) = Mxy;
  M_line(2, 2) = Myy;
  M_line(2, 3) = My;
  M_line(3, 0) = Mz;
  M_line(3, 1) = Mx;
  M_line(3, 2) = My;
  M_line(3, 3) = terms[9];

  std::cout << "Moment matrix M (line segment):\n" << M_line << std::endl;

  // constraint matrix
  Eigen::Matrix4d C_line = Eigen::Matrix4d::Zero();
  C_line(0, 0) = 4.0 * M_line(0, 3);
  C_line(0, 1) = 2.0 * M_line(1, 3);
  C_line(0, 2) = 2.0 * M_line(2, 3);
  C_line(1, 0) = C_line(0, 1);
  C_line(1, 1) = M_line(3, 3);
  C_line(2, 0) = C_line(0, 2);
  C_line(2, 2) = M_line(3, 3);
  std::cout << "Constraint matrix C (line segment):\n" << C_line << std::endl;
}

TEST(TaubinTest, CircleFit) {
  // line segment end points
  std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> end_points = {
      {Eigen::Vector2d(0.0867, -0.0469), Eigen::Vector2d(0.0938, -0.0564)},
      {Eigen::Vector2d(0.0915, -0.0313), Eigen::Vector2d(0.0904, -0.0469)},
      {Eigen::Vector2d(0.1017, -0.0625), Eigen::Vector2d(0.1094, -0.0666)},
      {Eigen::Vector2d(0.0938, -0.0595), Eigen::Vector2d(0.1024, -0.0625)},
      {Eigen::Vector2d(0.0940, -0.0313), Eigen::Vector2d(0.0938, -0.0329)},
      {Eigen::Vector2d(0.1094, -0.0623), Eigen::Vector2d(0.1250, -0.0545)},
      {Eigen::Vector2d(0.1250, -0.0530), Eigen::Vector2d(0.1309, -0.0469)},
      {Eigen::Vector2d(0.1307, -0.0469), Eigen::Vector2d(0.1398, -0.0313)}};

  // ==== point sampling formulation =====
  std::vector<Eigen::Vector2d> points;
  const int n_points_per_segment = 10;
  for (const auto& ep : end_points) {
    const Eigen::Vector2d& p0 = ep.first;
    const Eigen::Vector2d& p1 = ep.second;
    for (int i = 0; i < n_points_per_segment; ++i) {
      double t = static_cast<double>(i) /
                 (static_cast<double>(n_points_per_segment) - 1.0);
      Eigen::Vector2d p = (1.0 - t) * p0 + t * p1;
      points.push_back(p);
    }
  }

  // moment matrix
  Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
  const int n_points = static_cast<int>(points.size());
  for (int i = 0; i < n_points; ++i) {
    double xi = points[i].x(), yi = points[i].y();
    double zi = xi * xi + yi * yi;
    Eigen::Vector4d u;
    u << zi, xi, yi, 1.0;
    M += u * u.transpose();
  }

  // constraint matrix
  Eigen::Matrix4d C = Eigen::Matrix4d::Zero();
  C(0, 0) = 4.0 * M(0, 3);
  C(0, 1) = 2.0 * M(1, 3);
  C(0, 2) = 2.0 * M(2, 3);
  C(1, 0) = C(0, 1);
  C(1, 1) = static_cast<double>(n_points);
  C(2, 0) = C(0, 2);
  C(2, 2) = static_cast<double>(n_points);

  // normalize matrices
  M /= static_cast<double>(n_points);
  C /= static_cast<double>(n_points);

  std::cout << "M = " << M << std::endl;
  std::cout << "C = " << C << std::endl;

  // solving generalized eigenvalue problem
  {
    Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
    ges.compute(M, C);
    const auto& evals = ges.eigenvalues();
    const auto& evecs = ges.eigenvectors();
    int bestIndex = -1;
    double bestLambda = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 4; i++) {
      const auto lam = evals(i);
      if (std::abs(lam.imag()) > 1e-9) continue;
      const double lambdaReal = lam.real();
      if (lambdaReal <= 0.0) continue;
      if (lambdaReal < bestLambda) {
        bestLambda = lambdaReal;
        bestIndex = i;
      }
    }

    // extracting eigenvector components
    Eigen::Vector4cd v_c = evecs.col(bestIndex);
    Eigen::Vector4d a = v_c.real();  // imaginary parts should be tiny
    double A = a(0);
    double B = a(1);
    double Cc = a(2);
    double D = a(3);
    const double num = B * B + Cc * Cc - 4.0 * A * D;
    const double R = 0.5 * std::sqrt(num) / std::abs(A);
    double xc_local = -B / (2.0 * A);
    double yc_local = -Cc / (2.0 * A);
    std::cout << "Fitted circle center: (" << xc_local << ", " << yc_local
              << "), radius: " << R << std::endl;
  }

  // ==== line segment formulation =====
  Eigen::Matrix4d M_line = Eigen::Matrix4d::Zero();
  Eigen::Matrix4d C_line = Eigen::Matrix4d::Zero();

  const int n_line_segments = end_points.size();
  for (int i = 0; i < n_line_segments; i++) {
    const Eigen::Vector2d& p0 = end_points[i].first;
    const Eigen::Vector2d& p1 = end_points[i].second;
    const double x0 = p0.x();
    const double y0 = p0.y();
    const double dx = p1.x() - p0.x();
    const double dy = p1.y() - p0.y();

    // moment matrix terms
    std::vector<double> terms(10, 0.0);

    // Mzz
    terms[0] = (dx * dx * dx * dx) / 5. + (2. * (dx * dx) * (dy * dy)) / 5. +
               (dy * dy * dy * dy) / 5. + dx * dx * dx * x0 +
               dx * (dy * dy) * x0 + 2. * (dx * dx) * (x0 * x0) +
               (2. * (dy * dy) * (x0 * x0)) / 3. + 2. * dx * (x0 * x0 * x0) +
               x0 * x0 * x0 * x0 + dx * dx * dy * y0 + dy * dy * dy * y0 +
               (8. * dx * dy * x0 * y0) / 3. + 2. * dy * (x0 * x0) * y0 +
               (2. * (dx * dx) * (y0 * y0)) / 3. + 2. * (dy * dy) * (y0 * y0) +
               2. * dx * x0 * (y0 * y0) + 2. * (x0 * x0) * (y0 * y0) +
               2. * dy * (y0 * y0 * y0) + y0 * y0 * y0 * y0;

    // Mxz
    terms[1] = (dx * dx * dx) / 4. + (dx * (dy * dy)) / 4. + dx * dx * x0 +
               (dy * dy * x0) / 3. + (3. * dx * (x0 * x0)) / 2. + x0 * x0 * x0 +
               (2. * dx * dy * y0) / 3. + dy * x0 * y0 + (dx * (y0 * y0)) / 2. +
               x0 * (y0 * y0);

    // Myz
    terms[2] = (dx * dx * dy) / 4. + (dy * dy * dy) / 4. +
               (2. * dx * dy * x0) / 3. + (dy * (x0 * x0)) / 2. +
               (dx * dx * y0) / 3. + dy * dy * y0 + dx * x0 * y0 +
               x0 * x0 * y0 + (3. * dy * (y0 * y0)) / 2. + y0 * y0 * y0;

    // Mz
    terms[3] =
        (dx * dx) / 3. + (dy * dy) / 3. + dx * x0 + x0 * x0 + dy * y0 + y0 * y0;

    // Mxx
    terms[4] = (dx * dx) / 3. + dx * x0 + x0 * x0;

    // Mxy
    terms[5] = (dx * dy) / 3. + (dy * x0) / 2. + (dx * y0) / 2. + x0 * y0;

    // Mx
    terms[6] = dx / 2. + x0;

    // Myy
    terms[7] = (dy * dy) / 3. + dy * y0 + y0 * y0;

    // My
    terms[8] = dy / 2. + y0;

    // Le
    terms[9] = 1.0;

    double Mzz = terms[0], Mxz = terms[1], Myz = terms[2], Mz = terms[3],
           Mxx = terms[4], Mxy = terms[5], Mx = terms[6], Myy = terms[7],
           My = terms[8];

    M_line(0, 0) += Mzz;
    M_line(0, 1) += Mxz;
    M_line(0, 2) += Myz;
    M_line(0, 3) += Mz;
    M_line(1, 0) += Mxz;
    M_line(1, 1) += Mxx;
    M_line(1, 2) += Mxy;
    M_line(1, 3) += Mx;
    M_line(2, 0) += Myz;
    M_line(2, 1) += Mxy;
    M_line(2, 2) += Myy;
    M_line(2, 3) += My;
    M_line(3, 0) += Mz;
    M_line(3, 1) += Mx;
    M_line(3, 2) += My;
    M_line(3, 3) += terms[9];
  }
  // constraint matrix
  C_line(0, 0) = 4.0 * M_line(0, 3);
  C_line(0, 1) = 2.0 * M_line(1, 3);
  C_line(0, 2) = 2.0 * M_line(2, 3);
  C_line(1, 0) = C_line(0, 1);
  C_line(1, 1) = M_line(3, 3);
  C_line(2, 0) = C_line(0, 2);
  C_line(2, 2) = M_line(3, 3);

  M_line /= static_cast<double>(n_line_segments);
  C_line /= static_cast<double>(n_line_segments);

  std::cout << "M_line = " << M_line << std::endl;
  std::cout << "C_line = " << C_line << std::endl;

  {
    Eigen::GeneralizedEigenSolver<Eigen::Matrix4d> ges;
    ges.compute(M_line, C_line);
    const auto& evals = ges.eigenvalues();
    const auto& evecs = ges.eigenvectors();
    int bestIndex = -1;
    double bestLambda = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 4; i++) {
      const auto lam = evals(i);
      if (std::abs(lam.imag()) > 1e-9) continue;
      const double lambdaReal = lam.real();
      if (lambdaReal <= 0.0) continue;
      if (lambdaReal < bestLambda) {
        bestLambda = lambdaReal;
        bestIndex = i;
      }
    }

    // extracting eigenvector components
    Eigen::Vector4cd v_c = evecs.col(bestIndex);
    Eigen::Vector4d a = v_c.real();  // imaginary parts should be tiny
    double A = a(0);
    double B = a(1);
    double Cc = a(2);
    double D = a(3);
    const double num = B * B + Cc * Cc - 4.0 * A * D;
    const double R = 0.5 * std::sqrt(num) / std::abs(A);
    double xc_local = -B / (2.0 * A);
    double yc_local = -Cc / (2.0 * A);
    std::cout << "Fitted circle center: (" << xc_local << ", " << yc_local
              << "), radius: " << R << std::endl;
  }
}

}  // namespace
