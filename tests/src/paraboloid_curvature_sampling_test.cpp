// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "gtest/gtest.h"

#include "Eigen/Dense"

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

#include "external/NumericalIntegration/NumericalIntegration.h"

namespace {

using namespace IRL;

TEST(ParaboloidCurvatureSamplingTest, CurvatureSamplingTest) {
  // defining paraboloid
  const double a = 1.;
  const double b = -2.;
  double x0 = 0.0, y0 = 0.0, z0 = 0.0;
  Pt datum(x0, y0, z0);
  Normal T1(1, 0, 0), T2(0, 1, 0), N(0, 0, 1);  // orthonormal basis
  ReferenceFrame frame(T1, T2, N);

  const int nslices = 101;
  const double theta_min = -M_PI;
  const double theta_max = M_PI;

  std::vector<double> coeffs(nslices);
  Eigen::Matrix3d Mtilde = Eigen::Matrix3d::Zero();

  for (int i = 0; i < nslices; i++) {
    double w = static_cast<double>(i) / (nslices - 1);
    double theta = theta_min + w * (theta_max - theta_min);

    // local slize tangential direction
    Normal t_hat(std::cos(theta), std::sin(theta), 0.);

    // coefficient for parabola slice
    double c_theta = a * std::cos(theta) * std::cos(theta) +
                     b * std::sin(theta) * std::sin(theta);

    coeffs[i] = c_theta;

    // std::cout << "Slice " << i + 1 << " Coefficient: " << c_theta
    //           << " angle: " << theta * 180. / M_PI << std::endl;

    const double weight_n = 1.0 / static_cast<double>(nslices);
    const double curvature_n = -2.0 * c_theta;
    const Eigen::Vector3d Ttheta_n(t_hat[0], t_hat[1], t_hat[2]);
    Mtilde += weight_n * curvature_n * Ttheta_n * Ttheta_n.transpose();
  }

  // recovering principal curvatures from sampled curvatures
  ReferenceFrame polygon_frame = frame;
  IRL::Normal W_plus = IRL::Normal(1, 0, 0) + polygon_frame[2];
  IRL::Normal W_minus = IRL::Normal(1, 0, 0) - polygon_frame[2];
  IRL::Normal W_max =
      (IRL::squaredMagnitude(W_plus) > IRL::squaredMagnitude(W_minus))
          ? W_plus
          : W_minus;
  W_max.normalize();
  const Eigen::Vector3d W(W_max[0], W_max[1], W_max[2]);
  const Eigen::Matrix3d Q =
      Eigen::Matrix3d::Identity() - 2.0 * W * W.transpose();
  // Compute minor matrix
  const Eigen::Matrix3d QMQ = Q.transpose() * Mtilde * Q;
  const double m11 = QMQ(1, 1);
  const double m12 = 0.5 * (QMQ(1, 2) + QMQ(2, 1));
  const double m22 = QMQ(2, 2);
  // Compute eigenvalues and Givens rotation angle
  const double tmp_sqrt =
      std::sqrt((m11 - m22) * (m11 - m22) + 4.0 * m12 * m12);
  const double lambda1 = 0.5 * (m11 + m22 + tmp_sqrt);
  const double lambda2 = 0.5 * (m11 + m22 - tmp_sqrt);
  const double theta = 0.5 * std::atan2(2.0 * m12, m11 - m22);
  // Compute principal directions (Darboux frame)
  const IRL::UnitQuaternion givens_rotation(theta, polygon_frame[2]);
  const auto darboux_frame = givens_rotation * polygon_frame;
  // Compute principal curvatures
  const double k1 = 3.0 * lambda1 - lambda2;
  const double k2 = 3.0 * lambda2 - lambda1;
  std::cout << "a = " << 0.5 * k1 << std::endl;
  std::cout << "b = " << 0.5 * k2 << std::endl;
}

TEST(ParaboloidCurvatureSamplingTest, CurvatureSamplingTest2) {
  const double a_true = 3.0;
  const double b_true = -5.0;

  // sample directional slice coefficients c(theta) = a cos^2 theta + b sin^2
  // theta
  const int nslices = 101;
  const double theta_min = -M_PI;
  const double theta_max = M_PI;

  std::vector<double> thetas(nslices), coeffs(nslices);
  for (int i = 0; i < nslices; ++i) {
    double w = static_cast<double>(i) / (nslices - 1);
    double theta = theta_min + w * (theta_max - theta_min);
    thetas[i] = theta;

    double ctheta = a_true * std::cos(theta) * std::cos(theta) +
                    b_true * std::sin(theta) * std::sin(theta);
    coeffs[i] = ctheta;
  }

  Eigen::Matrix3d XtX = Eigen::Matrix3d::Zero();
  Eigen::Vector3d Xty = Eigen::Vector3d::Zero();

  for (int i = 0; i < nslices; ++i) {
    const double theta = thetas[i];
    const double x0 = 1.0;
    const double x1 = std::cos(2.0 * theta);
    const double x2 = std::sin(2.0 * theta);
    const double y = coeffs[i];

    // rank-1 update XtX += x x^T; Xty += x y
    XtX(0, 0) += x0 * x0;
    XtX(0, 1) += x0 * x1;
    XtX(0, 2) += x0 * x2;
    XtX(1, 0) += x1 * x0;
    XtX(1, 1) += x1 * x1;
    XtX(1, 2) += x1 * x2;
    XtX(2, 0) += x2 * x0;
    XtX(2, 1) += x2 * x1;
    XtX(2, 2) += x2 * x2;

    Xty(0) += x0 * y;
    Xty(1) += x1 * y;
    Xty(2) += x2 * y;
  }

  Eigen::Vector3d p = XtX.ldlt().solve(Xty);
  const double alpha = p(0);
  const double beta = p(1);
  const double gamma = p(2);

  const double R = std::sqrt(beta * beta + gamma * gamma);
  const double a_est = alpha + R;
  const double b_est = alpha - R;
  const double phi = 0.5 * std::atan2(gamma, beta);  // principal-axis rotation

  const double tol = 1e-12;
  EXPECT_NEAR(a_true, a_est, 1e-10);
  EXPECT_NEAR(b_true, b_est, 1e-10);
  EXPECT_NEAR(0.0, gamma, tol);

  std::cout << "Recovered a = " << a_est << ", b = " << b_est
            << ", phi (rad) = " << phi << std::endl;
}

}  // namespace
