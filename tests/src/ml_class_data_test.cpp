// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Dense>
#include "gtest/gtest.h"

#include "irl/ml_classification/data_gen.h"

namespace {

using namespace IRL;

std::vector<std::vector<std::vector<double>>> makeScalarField(
    const int stencil_size, const double value = 0.0) {
  return std::vector<std::vector<std::vector<double>>>(
      stencil_size,
      std::vector<std::vector<double>>(
          stencil_size, std::vector<double>(stencil_size, value)));
}

std::vector<std::vector<std::vector<Eigen::Vector3d>>> makeVectorField(
    const int stencil_size) {
  return std::vector<std::vector<std::vector<Eigen::Vector3d>>>(
      stencil_size,
      std::vector<std::vector<Eigen::Vector3d>>(
          stencil_size,
          std::vector<Eigen::Vector3d>(stencil_size,
                                       Eigen::Vector3d::Zero())));
}

double sumVolumeFractions(
    const std::vector<std::vector<std::vector<double>>>& vfrac) {
  double volume = 0.0;
  for (const auto& plane : vfrac) {
    for (const auto& row : plane) {
      for (const double value : row) {
        volume += value;
      }
    }
  }
  return volume;
}

Eigen::Vector3d sumFirstMoments(
    const std::vector<std::vector<std::vector<Eigen::Vector3d>>>&
        first_moment) {
  Eigen::Vector3d moment = Eigen::Vector3d::Zero();
  for (const auto& plane : first_moment) {
    for (const auto& row : plane) {
      for (const auto& value : row) {
        moment += value;
      }
    }
  }
  return moment;
}

TEST(EllipsoidGeneration, SpecificEllipsoidMatchesAnalyticalMoments) {
  constexpr int stencil_size = 5;

  // The ellipsoid is fully contained inside the 5^3 stencil, whose physical
  // bounds are [-2.5, 2.5]^3.
  const Eigen::Vector3d origin(0.20, -0.15, 0.10);

  // Use a non-axis-aligned ellipsoid so this also exercises the supplied
  // orientation vectors.
  const double inv_sqrt_2 = 1.0 / std::sqrt(2.0);
  const Eigen::Vector3d axis0_direction(inv_sqrt_2, inv_sqrt_2, 0.0);
  const Eigen::Vector3d axis1_direction(-inv_sqrt_2, inv_sqrt_2, 0.0);
  const Eigen::Vector3d axis2_direction(0.0, 0.0, 1.0);

  constexpr double axis0 = 1.40;
  constexpr double axis1 = 0.90;
  constexpr double axis2 = 0.65;

  // Analytical values for the geometry above:
  // V = 4/3*pi*a*b*c
  constexpr double expected_volume = 3.430619177720054;
  const Eigen::Vector3d expected_centroid(0.20, -0.15, 0.10);
  const Eigen::Vector3d expected_first_moment(
      0.686123835544011, -0.514592876658008, 0.343061917772005);

  // generateSpecificEllipsoid uses numerical geometric integration, so these
  // tolerances should reflect discretization error rather than roundoff error.
  constexpr double volume_tolerance = 2.0e-3;
  constexpr double centroid_tolerance = 2.0e-3;
  constexpr double first_moment_tolerance = 3.0e-3;

  Data_gen data_gen;
  data_gen.set_stencil_size(stencil_size);
  data_gen.setVisualize(false);

  auto vfrac = makeScalarField(stencil_size);
  auto first_moment = makeVectorField(stencil_size);
  auto area = makeScalarField(stencil_size);
  auto centroid = makeVectorField(stencil_size);

  std::vector<ParaboloidParametrizedSurfaceOutput> surfaces;
  std::vector<double> coarse_coords(stencil_size + 1, 0.0);

  data_gen.generateSpecificEllipsoid(
      vfrac, first_moment, area, centroid, surfaces, origin, axis0_direction,
      axis1_direction, axis2_direction, axis0, axis1, axis2, coarse_coords);

  const double computed_volume = sumVolumeFractions(vfrac);
  const Eigen::Vector3d computed_first_moment =
      sumFirstMoments(first_moment);

  ASSERT_GT(computed_volume, 1.0e-12);

  const Eigen::Vector3d computed_centroid =
      computed_first_moment / computed_volume;

  EXPECT_NEAR(computed_volume, expected_volume, volume_tolerance);

  for (int d = 0; d < 3; ++d) {
    EXPECT_NEAR(computed_centroid[d], expected_centroid[d],
                centroid_tolerance);
    EXPECT_NEAR(computed_first_moment[d], expected_first_moment[d],
                first_moment_tolerance);
  }
}

TEST(EllipsoidGeneration, EqualAxesAgreeWithSphere) {
  constexpr int stencil_size = 5;
  constexpr double radius = 0.90;

  const Eigen::Vector3d origin(0.15, -0.10, 0.20);
  const Eigen::Vector3d axis0_direction(1.0, 0.0, 0.0);
  const Eigen::Vector3d axis1_direction(0.0, 1.0, 0.0);
  const Eigen::Vector3d axis2_direction(0.0, 0.0, 1.0);

  // Analytical sphere volume: 4/3*pi*r^3.
  constexpr double expected_volume = 3.053628059289279;

  // These two generators use different geometric paths.  The cellwise
  // comparison therefore allows a somewhat larger error than the global
  // volume/centroid checks.
  constexpr double analytical_volume_tolerance = 2.0e-3;
  constexpr double global_volume_difference_tolerance = 2.0e-3;
  constexpr double centroid_tolerance = 2.0e-3;
  constexpr double max_cell_vfrac_difference = 1.0e-2;
  constexpr double rms_cell_vfrac_difference = 2.0e-3;
  constexpr double max_cell_first_moment_difference = 1.0e-2;
  constexpr double rms_cell_first_moment_difference = 2.0e-3;

  Data_gen data_gen;
  data_gen.set_stencil_size(stencil_size);
  data_gen.setVisualize(false);

  auto sphere_vfrac = makeScalarField(stencil_size);
  auto sphere_first_moment = makeVectorField(stencil_size);
  auto sphere_area = makeScalarField(stencil_size);
  auto sphere_centroid = makeVectorField(stencil_size);
  std::vector<ParaboloidParametrizedSurfaceOutput> sphere_surfaces;
  std::vector<double> sphere_coarse_coords(stencil_size + 1, 0.0);

  auto ellipsoid_vfrac = makeScalarField(stencil_size);
  auto ellipsoid_first_moment = makeVectorField(stencil_size);
  auto ellipsoid_area = makeScalarField(stencil_size);
  auto ellipsoid_centroid = makeVectorField(stencil_size);
  std::vector<ParaboloidParametrizedSurfaceOutput> ellipsoid_surfaces;
  std::vector<double> ellipsoid_coarse_coords(stencil_size + 1, 0.0);

  data_gen.generateSpecificSphere(
      sphere_vfrac, sphere_first_moment, sphere_area, sphere_centroid,
      sphere_surfaces, origin, radius, sphere_coarse_coords);

  data_gen.generateSpecificEllipsoid(
      ellipsoid_vfrac, ellipsoid_first_moment, ellipsoid_area,
      ellipsoid_centroid, ellipsoid_surfaces, origin, axis0_direction,
      axis1_direction, axis2_direction, radius, radius, radius,
      ellipsoid_coarse_coords);

  double max_vfrac_difference = 0.0;
  double sum_squared_vfrac_difference = 0.0;
  double max_first_moment_difference = 0.0;
  double sum_squared_first_moment_difference = 0.0;
  int cell_count = 0;

  for (int i = 0; i < stencil_size; ++i) {
    for (int j = 0; j < stencil_size; ++j) {
      for (int k = 0; k < stencil_size; ++k) {
        const double dv =
            ellipsoid_vfrac[i][j][k] - sphere_vfrac[i][j][k];
        max_vfrac_difference =
            std::max(max_vfrac_difference, std::abs(dv));
        sum_squared_vfrac_difference += dv * dv;

        const Eigen::Vector3d dm =
            ellipsoid_first_moment[i][j][k] -
            sphere_first_moment[i][j][k];
        const double dm_norm = dm.norm();
        max_first_moment_difference =
            std::max(max_first_moment_difference, dm_norm);
        sum_squared_first_moment_difference += dm_norm * dm_norm;

        ++cell_count;
      }
    }
  }

  const double rms_vfrac_difference =
      std::sqrt(sum_squared_vfrac_difference /
                static_cast<double>(cell_count));
  const double rms_first_moment_difference =
      std::sqrt(sum_squared_first_moment_difference /
                static_cast<double>(cell_count));

  const double sphere_volume = sumVolumeFractions(sphere_vfrac);
  const double ellipsoid_volume = sumVolumeFractions(ellipsoid_vfrac);

  const Eigen::Vector3d sphere_total_first_moment =
      sumFirstMoments(sphere_first_moment);
  const Eigen::Vector3d ellipsoid_total_first_moment =
      sumFirstMoments(ellipsoid_first_moment);

  ASSERT_GT(sphere_volume, 1.0e-12);
  ASSERT_GT(ellipsoid_volume, 1.0e-12);

  const Eigen::Vector3d sphere_recovered_centroid =
      sphere_total_first_moment / sphere_volume;
  const Eigen::Vector3d ellipsoid_recovered_centroid =
      ellipsoid_total_first_moment / ellipsoid_volume;

  // Each implementation should recover the known sphere volume.
  EXPECT_NEAR(sphere_volume, expected_volume, analytical_volume_tolerance);
  EXPECT_NEAR(ellipsoid_volume, expected_volume,
              analytical_volume_tolerance);

  // Equal ellipsoid semi-axes should reproduce the sphere generator.
  EXPECT_NEAR(ellipsoid_volume, sphere_volume,
              global_volume_difference_tolerance);
  EXPECT_LT(max_vfrac_difference, max_cell_vfrac_difference);
  EXPECT_LT(rms_vfrac_difference, rms_cell_vfrac_difference);
  EXPECT_LT(max_first_moment_difference,
            max_cell_first_moment_difference);
  EXPECT_LT(rms_first_moment_difference,
            rms_cell_first_moment_difference);

  for (int d = 0; d < 3; ++d) {
    EXPECT_NEAR(sphere_recovered_centroid[d], origin[d],
                centroid_tolerance);
    EXPECT_NEAR(ellipsoid_recovered_centroid[d], origin[d],
                centroid_tolerance);
  }
}

}  // namespace
