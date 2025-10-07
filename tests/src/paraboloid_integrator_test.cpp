// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "gtest/gtest.h"

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

#include "external/NumericalIntegration/NumericalIntegration.h"

namespace {

using namespace IRL;

TEST(ParaboloidIntegrator, ParaboloidCapMoments) {
  const std::size_t max_nsubdivisions = 128;
  const double epsabs = 10.0 * DBL_EPSILON;
  const double epsrel = 10.0 * DBL_EPSILON;
  const double a = 0.0;
  const double b = 10.0;

  // parameterized paraboloid surface
  using VolumeMomentsAndSurface =
      AddSurfaceOutput<VolumeMoments, ParaboloidParametrizedSurfaceOutput>;
  auto cell = RectangularCuboid::fromBoundingPts(IRL::Pt(0., 0, 0.),
                                                 IRL::Pt(1., 1., 1.));
  double x0 = 0.5, y0 = 0.5, z0 = 0.1, alpha = 5.;
  Pt datum(x0, y0, z0);
  ReferenceFrame frame(IRL::Normal(1, 0, 0), IRL::Normal(0, 1, 0),
                       IRL::Normal(0, 0, 1));
  Paraboloid paraboloid(datum, frame, alpha, alpha);
  auto surface =
      getVolumeMoments<VolumeMomentsAndSurface>(cell, paraboloid).getSurface();

  // M0 exact
  const double M0_exact =
      M_PI / (6. * alpha * alpha) * (-1. + std::pow(1. + 4. * alpha * z0, 1.5));

  // M1 exact
  double Mx_exact = (M_PI * x0 * (-1 + std::pow(1 + 4 * z0 * alpha, 1.5))) /
                    (6. * std::pow(alpha, 2));
  double My_exact = (M_PI * y0 * (-1 + std::pow(1 + 4 * z0 * alpha, 1.5))) /
                    (6. * std::pow(alpha, 2));
  double Mz_exact =
      (M_PI * (-1 + std::sqrt(1 + 4 * z0 * alpha) +
               2 * z0 * alpha *
                   (-5 + 4 * std::sqrt(1 + 4 * z0 * alpha) +
                    8 * z0 * alpha * std::sqrt(1 + 4 * z0 * alpha)))) /
      (60. * std::pow(alpha, 3));

  // M2 exact
  double Mxx_exact =
      (M_PI * (1 - std::pow(1 + 4 * z0 * alpha, 1.5) +
               6 * z0 * alpha * std::pow(1 + 4 * z0 * alpha, 1.5) +
               20 * std::pow(x0, 2) * std::pow(alpha, 2) *
                   (-1 + std::pow(1 + 4 * z0 * alpha, 1.5)))) /
      (120. * std::pow(alpha, 4));
  double Mxy_exact =
      (M_PI * x0 * y0 * (-1 + std::pow(1 + 4 * z0 * alpha, 1.5))) /
      (6. * std::pow(alpha, 2));
  double Mxz_exact = (M_PI * x0 *
                      (-1 + std::sqrt(1 + 4 * z0 * alpha) +
                       2 * z0 * alpha *
                           (-5 + 4 * std::sqrt(1 + 4 * z0 * alpha) +
                            8 * z0 * alpha * std::sqrt(1 + 4 * z0 * alpha)))) /
                     (60. * std::pow(alpha, 3));
  double Myy_exact =
      (M_PI * (1 - std::pow(1 + 4 * z0 * alpha, 1.5) +
               6 * z0 * alpha * std::pow(1 + 4 * z0 * alpha, 1.5) +
               20 * std::pow(y0, 2) * std::pow(alpha, 2) *
                   (-1 + std::pow(1 + 4 * z0 * alpha, 1.5)))) /
      (120. * std::pow(alpha, 4));
  double Myz_exact = (M_PI * y0 *
                      (-1 + std::sqrt(1 + 4 * z0 * alpha) +
                       2 * z0 * alpha *
                           (-5 + 4 * std::sqrt(1 + 4 * z0 * alpha) +
                            8 * z0 * alpha * std::sqrt(1 + 4 * z0 * alpha)))) /
                     (60. * std::pow(alpha, 3));
  double Mzz_exact =
      (M_PI *
       (-1 - 14 * z0 * alpha - 70 * std::pow(z0, 2) * std::pow(alpha, 2) +
        std::pow(1 + 4 * z0 * alpha, 1.5) +
        8 * z0 * alpha * std::pow(1 + 4 * z0 * alpha, 1.5) +
        16 * std::pow(z0, 2) * std::pow(alpha, 2) *
            std::pow(1 + 4 * z0 * alpha, 1.5))) /
      (420. * std::pow(alpha, 4));

  // numerical solution
  GeneralSurfaceMoments3D<2> m = surface.getSurfaceMoments<2>();

  EXPECT_NEAR(M0_exact, m[0], 10.0 * DBL_EPSILON);
  EXPECT_NEAR(Mx_exact, m[1], 10.0 * DBL_EPSILON);
  EXPECT_NEAR(My_exact, m[2], 10.0 * DBL_EPSILON);
  EXPECT_NEAR(Mz_exact, m[3], 10.0 * DBL_EPSILON);
  EXPECT_NEAR(Mxx_exact, m[4], 10.0 * DBL_EPSILON);
  EXPECT_NEAR(Mxy_exact, m[5], 10.0 * DBL_EPSILON);
  EXPECT_NEAR(Mxz_exact, m[6], 10.0 * DBL_EPSILON);
  EXPECT_NEAR(Myy_exact, m[7], 10.0 * DBL_EPSILON);
  EXPECT_NEAR(Myz_exact, m[8], 10.0 * DBL_EPSILON);
  EXPECT_NEAR(Mzz_exact, m[9], 10.0 * DBL_EPSILON);
}

}  // namespace
