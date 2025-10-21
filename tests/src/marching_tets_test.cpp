// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.


#include "gtest/gtest.h"

#include "irl/surface_mesher/marching_tets.h"

namespace {
using namespace IRL;

TEST(MarchingTets, SphereOctant) {
  const int nsubdivision = 500;
  auto cube = RectangularCuboid::fromBoundingPts(Pt(0, 0, 0), Pt(1, 1, 1));
  auto unit_sphere = [](Pt p) {
    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - 1.0;
  };
  MarchingTets mt(cube, unit_sphere);
  auto surface = mt.triangulate(nsubdivision);
  surface.write("test_marching_tets_sphere_octant_correct_500");
  SUCCEED();
}

}  // namespace