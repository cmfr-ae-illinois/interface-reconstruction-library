// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cmath>
#include <iomanip>

#include "gtest/gtest.h"

#include "irl/quadratic_reconstruction/coons_patch.h"

namespace {

using namespace IRL;

TEST(CoonsPatch, JacobianTest) {
  // bounding arcs of parallelogram (twice the size of unit parallelogram)
  RationalBezierArc arc_0(Pt(4., 0., 0.), Pt(5., 0., 0.), Pt(6., 0., 0.), 1.0);
  RationalBezierArc arc_1(Pt(4., 0., 0.), Pt(5.5, 1., 0.), Pt(7., 2., 0.), 1.0);
  RationalBezierArc arc_2(Pt(5., 2., 0.), Pt(6., 2., 0.), Pt(7., 2., 0.), 1.0);
  RationalBezierArc arc_3(Pt(4., 0., 0.), Pt(4.5, 1., 0.), Pt(5., 2., 0.), 1.0);
  CoonsPatch coons(arc_0, arc_1, arc_2, arc_3);

  double expected_detJ = 4.0;

  int Nu = 100, Nv = 100;
  double u, v;
  for (int i = 0; i < Nu; i++) {
    for (int j = 0; j < Nv; j++) {
      u = 1. / Nu * i;
      v = 1. / Nv * j;
      double detJ = coons.detJacobian(u, v);
      EXPECT_NEAR(expected_detJ, detJ, 100.0 * DBL_EPSILON);
    }
  }
}

}  // namespace
