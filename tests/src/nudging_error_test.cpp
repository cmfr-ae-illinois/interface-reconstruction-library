// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "gtest/gtest.h"

#include <algorithm>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

#include "irl/geometry/general/pt.h"

namespace {

using namespace IRL;

TEST(NudgingErrorTest, Nudging) {
  // defining paraboloid
  double a = 2.872483563951871284942096371195e-06;
  double b = 1.934238742397893506819173126132e+00;
  Pt datum(-1.033447268666018770133518955845e+00,
           9.999999999931639127481730611180e-01,
           -5.573730449772535955332841695054e-01);
  ReferenceFrame frame(IRL::Normal(-8.506464644355833382149967292207e-01,
                                   -5.812027452668506836974682272512e-12,
                                   5.257381406586761896093662471685e-01),
                       IRL::Normal(5.257381406447728666719854118128e-01,
                                   7.272588000451114639735138850396e-06,
                                   8.506464644130877772454368823674e-01),
                       IRL::Normal(-3.823481837114372293797738555643e-06,
                                   9.999999999735547095980336962384e-01,
                                   -6.186398214275882284846671405054e-06));
  Paraboloid interface(datum, frame, a, b);

  // defining cell
  Pt x0(-1.03369140625, 0.99951171875, -0.5576171875);
  Pt x1(-1.033203125, 1.0, -0.55712890625);
  RectangularCuboid cell = RectangularCuboid::fromBoundingPts(x0, x1);

  // computing moments
  auto moments = getVolumeMoments<VolumeMoments>(cell, interface);
  std::cout << "Moments = " << moments << std::endl;

  // viz
  {
    using VolumeAndSuface =
        AddSurfaceOutput<Volume, ParaboloidParametrizedSurfaceOutput>;

    AlignedParaboloid aligned_paraboloid(
        {2.42933958821834e-09, 0.00163583973427797});
    Pt datum(0, 0, 0);
    ReferenceFrame frame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
    Paraboloid paraboloid(datum, frame, aligned_paraboloid.a(),
                          aligned_paraboloid.b());
    RectangularCuboid viz_cell = RectangularCuboid::fromBoundingPts(
        -Pt(0.5, 0.5, 0.5), Pt(0.5, 0.5, 0.5));
    auto temp_surface_and_moments =
        getVolumeMoments<VolumeAndSuface>(viz_cell, paraboloid);
    auto temp_param_surface = temp_surface_and_moments.getSurface();
    auto temp_tri_surface = temp_param_surface.triangulate(0.01);
    temp_tri_surface.write("nudging-debug-surface");
  }
}

}  // namespace
