// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
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

TEST(ParaboloidIntegrator, ParaboloidCapSurfaceArea) {
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

  // analytical surface area of capped paraboloid
  auto surfaceArea_ana = [](const double& alpha, const double& z0) {
    return M_PI / (6. * alpha * alpha) *
           (-1. + std::pow(1. + 4. * alpha * z0, 1.5));
  };
  const double exact_surfaceArea = surfaceArea_ana(alpha, z0);

  // numerical surface area
  double surfaceArea = surface.getIntegrator([](const Pt&) { return 1.0; });

  std::cout << "Capped Paraboloid Surface Area (exact)     = "
            << std::setprecision(16) << exact_surfaceArea << std::endl;

  std::cout << "Capped Paraboloid Surface Area (numerical)     = "
            << std::setprecision(16) << surfaceArea << std::endl;

  EXPECT_NEAR(surfaceArea, exact_surfaceArea,
              std::max(epsabs, epsrel * exact_surfaceArea));
}

}  // namespace
