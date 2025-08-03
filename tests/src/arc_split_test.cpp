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
#include "irl/quadratic_reconstruction/rational_bezier_arc.h"

namespace {

using namespace IRL;

// Expected location of control point
TEST(ArcSplit, ControlPoint) {
  RationalBezierArc arc(Pt(1.25, 1.2, 0.0), Pt(1.25, 0.5, 0.0),
                        Pt(2.0, 0.2, 0.0), 2.0);
  auto [arc_1, arc_2] = arc.split();

  // expected control points for each arc
  Pt arc_1_cp = 1. / (1. + arc.weight()) *
                (arc.start_point() + arc.weight() * arc.control_point());
  Pt arc_2_cp = 1. / (1. + arc.weight()) *
                (arc.end_point() + arc.weight() * arc.control_point());

  std::cout << "Expected Arc 1 Control point: " << arc_1_cp << std::endl;
  std::cout << "Arc 1 Control point   : " << arc_1.control_point() << std::endl;
  std::cout << "Expected Arc 2 Control point: " << arc_2_cp << std::endl;
  std::cout << "Arc 2 Control point   : " << arc_2.control_point() << std::endl;

  for (int i = 0; i < 2; i++) {
    EXPECT_NEAR(arc_1_cp[i], arc_1.control_point()[i], 10.0 * DBL_EPSILON);
    EXPECT_NEAR(arc_2_cp[i], arc_2.control_point()[i], 10.0 * DBL_EPSILON);
  }
}

// expected weight of control point (for a quadratic bezier curve representing a
// circular arc since analytical expressions are known for comparison)
TEST(ArcSplit, Weights) {
  // arc list
  auto cell = IRL::RectangularCuboid::fromBoundingPts(IRL::Pt(0., 0, 0.),
                                                      IRL::Pt(1., 1., 1.));
  IRL::Pt datum(0.5, 0.5, 0.1);
  IRL::ReferenceFrame frame(IRL::Normal(1, 0, 0), IRL::Normal(0, 1, 0),
                            IRL::Normal(0, 0, 1));
  IRL::Paraboloid paraboloid(datum, frame, 5.0, 5.0);
  using VolumeMomentsAndSurface =
      IRL::AddSurfaceOutput<IRL::VolumeMoments,
                            IRL::ParaboloidParametrizedSurfaceOutput>;
  auto surface =
      IRL::getVolumeMoments<VolumeMomentsAndSurface>(cell, paraboloid)
          .getSurface();
  auto arcList = surface.getArcs();

  // expected weight of each arc when splt
  double num_split_arcs = arcList.size() * 2.;
  double theta = 360. / num_split_arcs * M_PI / 180.;
  double w_split_arcs_exact = std::cos(theta / 2.);
  std::cout << "Expected weight of each split arc   = " << std::setprecision(16)
            << w_split_arcs_exact << std::endl;

  // weight obtained after splitting
  for (const auto& arc : arcList) {
    auto [arc_1, arc_2] = arc.split();
    std::cout << "Weight of split arc      = " << std::setprecision(16)
              << arc_1.weight() << std::endl;
    std::cout << "Weight of split arc      = " << std::setprecision(16)
              << arc_2.weight() << std::endl;
    EXPECT_NEAR(w_split_arcs_exact, arc_1.weight(), 10.0 * DBL_EPSILON);
    EXPECT_NEAR(w_split_arcs_exact, arc_2.weight(), 10.0 * DBL_EPSILON);
  }
}

}  // namespace
