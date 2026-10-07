// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Fabien Evrard <fa.evrard@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "gtest/gtest.h"

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

namespace {

using namespace IRL;

TEST(ParaboloidPolygonIntersection, TestWithoutNudge) {
  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));

  const auto paraboloid = Paraboloid(datum, frame, 0.5, 0.5);

  Polygon polygon;
  polygon.setNumberOfVertices(4);
  polygon[0] = Pt(0, 0, -1);
  polygon[1] = Pt(1, 0, -1);
  polygon[2] = Pt(1, 0, 1);
  polygon[3] = Pt(0, 0, 1);

  auto intersection = intersectPolygonWithParaboloid(polygon, paraboloid);
  const auto arcs = intersection.getArcs();

  std::cout << "Intersection contains " << arcs.size() << " arcs." << std::endl;
  for (const auto& arc : arcs) {
    std::cout << "Arc: " << arc << std::endl;
  }

  EXPECT_EQ(arcs.size(), 1);

  const double expected_arc_length =
      std::sqrt(2.0) / 2.0 + std::log(1.0 + std::sqrt(2.0)) / 2.0;
  EXPECT_NEAR(arcs[0].arc_length(), expected_arc_length, 1e-3);
}

TEST(ParaboloidPolygonIntersection, TestWithNudge) {
  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto paraboloid = Paraboloid(datum, frame, 1.0, 1.0);

  Polygon polygon;
  polygon.setNumberOfVertices(4);
  polygon[0] = Pt(-1, 0, -1);
  polygon[1] = Pt(1, 0, -1);
  polygon[2] = Pt(1, 0, 1);
  polygon[3] = Pt(-1, 0, 1);

  auto intersection = intersectPolygonWithParaboloid(polygon, paraboloid);
  const auto arcs = intersection.getArcs();

  std::cout << "Intersection contains " << arcs.size() << " arcs." << std::endl;
  double arc_length = 0.0;
  for (const auto& arc : arcs) {
    std::cout << "Arc: " << arc << std::endl;
    arc_length += arc.arc_length();
  }
  const double expected_arc_length =
      std::sqrt(5.0) + std::log(2.0 + std::sqrt(5.0)) / 2.0;
  EXPECT_NEAR(arc_length, expected_arc_length, 1e-3);
}

}  // namespace
