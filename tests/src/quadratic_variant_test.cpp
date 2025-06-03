// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// #define DEBUG_CYL_IRL
// #define NUDGE_REGION

#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/rotations.h"

#include <sys/time.h>
#include <cmath>
#include <random>
#include <variant>

#include "gtest/gtest.h"

#include "irl/generic_cutting/generic_cutting.h"

#include "irl/generic_cutting/cylinder_intersection/cylinder_intersection.h"
#include "irl/moments/general_moments.h"

#include "irl/cylinder_reconstruction/cylinder.h"
#include "irl/interface_reconstruction_methods/progressive_radius_solver_cylinder.h"
#include "irl/quadratic_reconstruction/parametrized_surface.h"

#include "irl/variant_reconstruction/separator_variant.h"

namespace {
using namespace IRL;

TEST(QuadraticVariant, Test1) {
  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto cell =
      RectangularCuboid::fromBoundingPts(Pt(-1, -1, -1), Pt(1, 1, 1));

  SeparatorVariant interface;
  interface = Paraboloid(datum, frame, 0.5, 0.5);

  Volume volume;
  if (const auto plane = std::get_if<PlanarSeparator>(&interface))
    volume = getVolumeMoments<Volume>(cell, *plane);
  else if (const auto paraboloid = std::get_if<Paraboloid>(&interface))
    volume = getVolumeMoments<Volume>(cell, *paraboloid);
  else if (const auto cylinder = std::get_if<Cylinder>(&interface))
    volume = getVolumeMoments<Volume>(cell, *cylinder);
  EXPECT_NEAR(volume, 4.0 - 4.0 / 3.0, 1.0e-14);

  interface = Cylinder(datum, frame, 1.0, 0.25);
  if (const auto plane = std::get_if<PlanarSeparator>(&interface))
    volume = getVolumeMoments<Volume>(cell, *plane);
  else if (const auto paraboloid = std::get_if<Paraboloid>(&interface))
    volume = getVolumeMoments<Volume>(cell, *paraboloid);
  else if (const auto cylinder = std::get_if<Cylinder>(&interface))
    volume = getVolumeMoments<Volume>(cell, *cylinder);
  EXPECT_NEAR(volume, M_PI * 0.25 * 2.0, 1.0e-14);
}

TEST(QuadraticVariant, Test2) {
  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto cell =
      RectangularCuboid::fromBoundingPts(Pt(-1, -1, -1), Pt(1, 1, 1));

  SeparatorVariant interface;
  interface = Paraboloid(datum, frame, 0.5, 0.5);

  auto volume = getVolumeMoments<Volume>(cell, interface);
  EXPECT_NEAR(volume, 4.0 - 4.0 / 3.0, 1.0e-14);

  interface = Cylinder(datum, frame, 1.0, 0.25);
  volume = getVolumeMoments<Volume>(cell, interface);
  EXPECT_NEAR(volume, M_PI * 0.25 * 2.0, 1.0e-14);
}
}  // namespace
