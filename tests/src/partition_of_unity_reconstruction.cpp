// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/rotations.h"

#include <sys/time.h>
#include <cmath>
#include <random>
#include <variant>

#include "gtest/gtest.h"

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/moments/general_moments.h"

#include "irl/variant_reconstruction/separator_variant.h"

namespace {
using namespace IRL;

void getIntersectionPts(const Polygon& a_polygon, const Plane& a_cutting_plane,
                        StackVector<IRL::Pt, 2>* a_intersection_pts) {
  a_intersection_pts->resize(0);
  if (a_polygon.getNumberOfVertices() > 0) {
    double distance = a_cutting_plane.signedDistanceToPoint(a_polygon[0]);
    for (int n = 0; n < a_polygon.getNumberOfVertices(); ++n) {
      const int next_id = (n + 1) % a_polygon.getNumberOfVertices();
      double next_distance =
          a_cutting_plane.signedDistanceToPoint(a_polygon[next_id]);
      if (distance * next_distance < 0.0) {
        a_intersection_pts->push_back(Pt::fromEdgeIntersection(
            a_polygon[n], distance, a_polygon[next_id], next_distance));
        if (a_intersection_pts->size() == 2) {
          break;
        }
      }
      distance = next_distance;
    }
  }
}

TEST(PUReconstruction, Test1) {
  const int nlayers = 2;
  const int ncells = (1 + 2 * nlayers) * (1 + 2 * nlayers);

  std::vector<PlanarSeparator> planar_separator(ncells);
  std::vector<RectangularCuboid> cells(ncells);

  // Generate planar separator corresponding to circle centred at (0,0)
  const auto sphere_center = Pt(0.0, 0.0, 0.0);
  const double sphere_radius = 2.5;
  int count = 0;
  for (int i = 0; i < 1 + 2 * nlayers; ++i) {
    for (int j = 0; j < 1 + 2 * nlayers; ++j) {
      const Pt lower_cell_pt(static_cast<double>(i), static_cast<double>(j),
                             -1.0);
      const Pt upper_cell_pt(static_cast<double>(i + 1),
                             static_cast<double>(j + 1), 1.0);
      const Pt mid_pt = 0.5 * (lower_cell_pt + upper_cell_pt);
      Normal normal = mid_pt - sphere_center;
      normal.normalize();
      planar_separator[count] =
          PlanarSeparator::fromOnePlane(Plane(normal, sphere_radius));
      cells[count] =
          RectangularCuboid::fromBoundingPts(lower_cell_pt, upper_cell_pt);
      count++;
    }
  }

  // Compute end points of 2D PLIC
  count = 0;
  StackVector<Pt, 2> intersections;
  const auto xy_plane = Plane(Normal(0.0, 0.0, 1.0), 0.0);
  for (int i = 0; i < 1 + 2 * nlayers; ++i) {
    for (int j = 0; j < 1 + 2 * nlayers; ++j) {
      const Polygon polygon = getPlanePolygonFromReconstruction<Polygon>(
          cells[count], planar_separator[count], planar_separator[count][0]);
      count++;
      getIntersectionPts(polygon, xy_plane, &intersections);
      if (intersections.size() == 2) {
        std::cout << "Start point for cell " << i << ", " << j << " = "
                  << intersections[0] << std::endl;
        std::cout << "  End point for cell " << i << ", " << j << " = "
                  << intersections[1] << std::endl;
      }
    }
  }

  SUCCEED();
}

}  // namespace
