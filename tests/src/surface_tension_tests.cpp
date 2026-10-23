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

TEST(SurfaceTension, areaComputation) {
  using AreaAndSurface =
      AddSurfaceOutput<Area, ParaboloidParametrizedSurfaceOutput>;
  // ============================================ TAKEN FROM SIMULATION
  const auto datum =
      Pt(-0.00764114344206871, 0.100706330883016, 0.476561469895421);
  const auto frame = ReferenceFrame(
      Normal(0.974019255888702, 0.226455957217311, -0.00204660665759728),
      Normal(-0.00233572054042677, 0.00100874691917091, -0.999996763414367),
      Normal(-0.226453159765053, 0.974020883693181, 0.00151147675436535));
  const auto paraboloid =
      Paraboloid(datum, frame, 15.3235133743669, -0.0935337715595459);
  const SeparatorVariant para = paraboloid;
  // Face
  Pt P0 =
      Pt(-1.5625000000000000E-002, 0.10156250000000000, 0.46875000000000000);
  Pt P1 =
      Pt(-1.5625000000000000E-002, 0.10156250000000000, 0.48437500000000000);
  Pt P2 = Pt(0.0000000000000000, 0.10156250000000000, 0.48437500000000000);
  Pt P3 = Pt(0.0000000000000000, 0.10156250000000000, 0.46875000000000000);
  // ST
  double ST = 7.1999999999999995E-002;
  double Pressure = 1.0;
  Normal Marangoni = Normal(0.0, 0.0, 0.0);

  // ======================================= Construct Problem
  Polygon polygon;
  polygon.setNumberOfVertices(4);
  polygon[0] = P0;
  polygon[1] = P1;
  polygon[2] = P2;
  polygon[3] = P3;
  polygon.calculateAndSetPlaneOfExistence();
  // Filler Solver
  PUNeighborhood<RectangularCuboid> neighborhood;
  PUST<RectangularCuboid> solver(neighborhood, 2.0);

  // ===================================== Run results
  // Solve Face Output
  Normal result;
  result = solver.solveFaceParaboloid(ST, P0, P1, P2, P3, &para, Pressure,
                                      Marangoni, true);
  // Area Isolation
  AreaAndSurface area_and_surf =
      intersectPolygonWithParaboloid<AreaAndSurface>(polygon, paraboloid, true);
  double area = area_and_surf.getMoments().area();

  // Compute Area
  std::cout << "===========================\n";
  std::cout << "solveFaceParaboloid RESULT = " << result << std::endl;
  std::cout << "areaResult " << area << std::endl;
  std::cout << "===========================\n";
}

}  // namespace