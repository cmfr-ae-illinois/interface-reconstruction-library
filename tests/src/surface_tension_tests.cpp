// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Ilia Kheirkhah <iliak2@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "gtest/gtest.h"

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/moments/general_moments.h"

#include "irl/interface_reconstruction_methods/pu.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/interface_reconstruction_methods/pu_paraboloid.h"
#include "irl/interface_reconstruction_methods/pust.h"
#include "irl/variant_reconstruction/separator_variant.h"

namespace {
using namespace IRL;

TEST(SurfaceTension, solveFaceParaboloid) {
  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto paraboloid = Paraboloid(datum, frame, 0.5, 0.5);
  const SeparatorVariant para = paraboloid;
  // Face to intersect
  Pt P0 = Pt(1, 0, 0);
  Pt P1 = Pt(1, 0, -1);
  Pt P2 = Pt(1, 2, -1);
  Pt P3 = Pt(1, 2, 0);

  // Options
  double ST = 1.0;
  double Pressure = 0.0;
  Normal Marangoni = Normal(0.0, 0.0, 0.0);

  // Filler Solver
  PUNeighborhood<RectangularCuboid> neighborhood;
  PUST<RectangularCuboid> solver(neighborhood, 2.0);

  // Solve the surface tension face force
  Normal result;
  result = solver.solveFaceParaboloid(ST, P0, P1, P2, P3, &para, Pressure,
                                      Marangoni);
  //   std::cout << "Result: " << result << std::endl;
  // Exact Values
  EXPECT_NEAR(result[0], 0.866025403784, 1e-8);
  EXPECT_NEAR(result[1], -0.317837245196, 1e-8);
  EXPECT_NEAR(result[2], -0.658478948462, 1e-8);
}

TEST(SurfaceTension, solveFaceParaboloid_Nudge) {
  const auto datum = Pt(0, 0, 0);
  const auto frame =
      ReferenceFrame(Normal(1, 0, 0), Normal(0, 1, 0), Normal(0, 0, 1));
  const auto paraboloid = Paraboloid(datum, frame, 0.5, 0.5);
  const SeparatorVariant para = paraboloid;
  // Face to intersect
  Pt P0 = Pt(1, 0, 0);
  Pt P1 = Pt(1, 0, -1);
  Pt P2 = Pt(1, 1, -1);
  Pt P3 = Pt(1, 1, 0);

  // Options
  double ST = 1.0;
  double Pressure = 0.0;
  Normal Marangoni = Normal(0.0, 0.0, 0.0);

  // Filler Solver
  PUNeighborhood<RectangularCuboid> neighborhood;
  PUST<RectangularCuboid> solver(neighborhood, 2.0);

  // Solve the surface tension face force
  Normal result;
  result = solver.solveFaceParaboloid(ST, P0, P1, P2, P3, &para, Pressure,
                                      Marangoni);
  //   std::cout << "Result: " << result << std::endl;
  // Exact Values
  EXPECT_NEAR(result[0], 0.866025403784, 1e-8);
  EXPECT_NEAR(result[1], -0.317837245196, 1e-8);
  EXPECT_NEAR(result[2], -0.658478948462, 1e-8);
}

}  // namespace