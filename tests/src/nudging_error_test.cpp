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
  GeneralMoments3D<2> moments =
      getVolumeMoments<GeneralMoments3D<2>>(cell, interface);
  std::cout << moments << std::endl;
}

}  // namespace
