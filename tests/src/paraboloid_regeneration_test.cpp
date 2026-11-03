// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cmath>

#include "Eigen/Dense"
#include "gtest/gtest.h"

#include "irl/paraboloid_reconstruction/paraboloid.h"

namespace {

using namespace IRL;

TEST(ParaboloidRegeneration, RegenerateAtLocation) {
  constexpr double tolerance = 1.0e-13;
  constexpr double coef_a = 0.35;
  constexpr double coef_b = 1.10;

  const ReferenceFrame frame(
      Normal(1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0), 0.0),
      Normal(-1.0 / std::sqrt(6.0), 1.0 / std::sqrt(6.0),
             2.0 * 1.0 / std::sqrt(6.0)),
      Normal(1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0),
             1.0 / std::sqrt(3.0)));
  const Pt datum(0.4, -0.7, 1.2);

  Paraboloid regenerated(datum, frame, coef_a, coef_b);

  const double local_x = 2.0;
  const double local_y = 0.0;
  const double local_z = 0.80;
  Pt point = datum;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    point[d] +=
        local_x * frame[0][d] + local_y * frame[1][d] + local_z * frame[2][d];
  }

  const double projected_local_z =
      -coef_a * local_x * local_x - coef_b * local_y * local_y;
  Pt projected_point = datum;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    projected_point[d] += local_x * frame[0][d] + local_y * frame[1][d] +
                          projected_local_z * frame[2][d];
  }

  const Eigen::Vector3d local_gradient(2.0 * coef_a * local_x,
                                       2.0 * coef_b * local_y, 1.0);
  Eigen::Matrix3d local_hessian = Eigen::Matrix3d::Zero();
  local_hessian(0, 0) = 2.0 * coef_a;
  local_hessian(1, 1) = 2.0 * coef_b;

  Eigen::Matrix3d local_to_world;
  for (UnsignedIndex_t i = 0; i < 3; ++i) {
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      local_to_world(i, j) = frame[j][i];
    }
  }
  const Eigen::Vector3d world_gradient = local_to_world * local_gradient;
  const Eigen::Matrix3d world_hessian =
      local_to_world * local_hessian * local_to_world.transpose();

  const Paraboloid from_derivatives = Paraboloid::fromDerivatives(
      projected_point, world_gradient, world_hessian);
  regenerated.regenerateAtLocation(point);

  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    EXPECT_NEAR(regenerated.getDatum()[d], from_derivatives.getDatum()[d],
                tolerance);
  }

  const auto& regenerated_frame = regenerated.getReferenceFrame();
  const auto& derivative_frame = from_derivatives.getReferenceFrame();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    EXPECT_NEAR(regenerated_frame[2][d], derivative_frame[2][d], tolerance);
  }

  const double direct_alignment =
      std::fabs(regenerated_frame[0] * derivative_frame[0]) +
      std::fabs(regenerated_frame[1] * derivative_frame[1]);
  const double swapped_alignment =
      std::fabs(regenerated_frame[0] * derivative_frame[1]) +
      std::fabs(regenerated_frame[1] * derivative_frame[0]);
  const bool directions_are_swapped = swapped_alignment > direct_alignment;

  const UnsignedIndex_t derivative_direction_for_a =
      directions_are_swapped ? 1 : 0;
  const UnsignedIndex_t derivative_direction_for_b =
      directions_are_swapped ? 0 : 1;
  EXPECT_NEAR(std::fabs(regenerated_frame[0] *
                        derivative_frame[derivative_direction_for_a]),
              1.0, tolerance);
  EXPECT_NEAR(std::fabs(regenerated_frame[1] *
                        derivative_frame[derivative_direction_for_b]),
              1.0, tolerance);

  const auto& regenerated_coefficients = regenerated.getAlignedParaboloid();
  const auto& derivative_coefficients = from_derivatives.getAlignedParaboloid();
  const double derivative_coefficient_for_a = directions_are_swapped
                                                  ? derivative_coefficients.b()
                                                  : derivative_coefficients.a();
  const double derivative_coefficient_for_b = directions_are_swapped
                                                  ? derivative_coefficients.a()
                                                  : derivative_coefficients.b();
  EXPECT_NEAR(regenerated_coefficients.a(), derivative_coefficient_for_a,
              tolerance);
  EXPECT_NEAR(regenerated_coefficients.b(), derivative_coefficient_for_b,
              tolerance);
}

}  // namespace
