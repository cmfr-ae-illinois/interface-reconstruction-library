// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Ilia Kheirkhah <iliak2@illinois.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/rotations.h"

#include <sys/time.h>
#include <cmath>
#include <random>
#include <tuple>
#include <variant>

#include "gtest/gtest.h"

#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/moments/general_moments.h"

#include "irl/helpers/wendland.h"
#include "irl/interface_reconstruction_methods/pu_neighborhood.h"
#include "irl/interface_reconstruction_methods/pu_solve.h"
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

TEST(Wendland, WendlandTests) {
  // ComputeR Tests
  double result;

  IRL::Pt xi1(0.0, 0.0, 0.0);
  IRL::Pt x1(2.0, 0.0, 0.0);
  result = Wendland::computeR(xi1, x1);
  EXPECT_EQ(result, 2.0) << "Single Direction ComputeR Fail";

  IRL::Pt xi2(0.0, 0.0, 0.0);
  IRL::Pt x2(1.0, 1.0, 1.0);
  result = Wendland::computeR(xi2, x2);
  EXPECT_NEAR(result, std::sqrt(3.0), std::numeric_limits<double>::epsilon())
      << "Three Direction ComputeR Fail";

  IRL::Pt xi3(-1.0, 0.0, 1.0);
  IRL::Pt x3(0.0, 0.0, 0.0);
  result = Wendland::computeR(xi3, x3);
  EXPECT_NEAR(result, std::sqrt(2.0), std::numeric_limits<double>::epsilon())
      << "Two Direction ComputeR Fail";

  // Value, First, and Second Ders at r=d
  double delta = 2.0;
  double r = 2.0;
  EXPECT_NEAR(Wendland::eval(r, delta), 0.0,
              std::numeric_limits<double>::epsilon())
      << "Wendland Right Endpoint Value Fail";  // Value
  EXPECT_NEAR(Wendland::firstDer(r, delta), 0.0,
              std::numeric_limits<double>::epsilon())
      << "Wendland Right Endpoint First Derivative Fail";  // First Der
  EXPECT_NEAR(Wendland::secondDer(r, delta), 0.0,
              std::numeric_limits<double>::epsilon())
      << "Wendland Right Endpoint Second Derivative Fail";  // Second Der

  // Value, First, and Second Der at r=0
  r = 0.0;
  EXPECT_NEAR(Wendland::eval(r, delta), 1,
              std::numeric_limits<double>::epsilon())
      << "Wendland Left Endpoint Value Fail";  // Value
  EXPECT_NEAR(Wendland::firstDer(r, delta), 0.0,
              std::numeric_limits<double>::epsilon())
      << "Wendland Left Endpoint First Derivative Fail";  // First Der
  EXPECT_NEAR(Wendland::secondDer(r, delta), -5,
              std::numeric_limits<double>::epsilon())
      << "Wendland Left Endpoint Second Derivative Fail";  // Second Der

  // Value, First, and Second Der at r=1
  r = 1.0;
  EXPECT_NEAR(Wendland::eval(r, delta), 0.1875,
              std::numeric_limits<double>::epsilon())
      << "Wendland Midpoint Value Fail";  // Value
  EXPECT_NEAR(Wendland::firstDer(r, delta), -0.625,
              std::numeric_limits<double>::epsilon())
      << "Wendland Midpoint First Derivative Fail";  // First Der
  EXPECT_NEAR(Wendland::secondDer(r, delta), 1.25,
              std::numeric_limits<double>::epsilon())
      << "Wendland Midpoint Second Derivative Fail";  // Second Der

  // Value, First, and Second Der for d=4, r=0
  r = 0.0;
  delta = 4.0;
  EXPECT_NEAR(Wendland::eval(r, delta), 1.0,
              std::numeric_limits<double>::epsilon())
      << "Wendland Left Endpoint, New Delta Value Fail";  // Value
  EXPECT_NEAR(Wendland::firstDer(r, delta), 0.0,
              std::numeric_limits<double>::epsilon())
      << "Wendland Left Endpoint, New Delta First Derivative Fail";  // First
  EXPECT_NEAR(Wendland::secondDer(r, delta), -1.25,
              std::numeric_limits<double>::epsilon())
      << "Wendland Left Endpoint, New Delta Second Derivative Fail";  // Second

  // Full Function Evaluation
  double res1;
  std::pair<double, Eigen::Vector3d> res2;
  std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> res3;
  // Expected Values
  double expectedValue = 0.1875;
  Eigen::Vector3d expectedGradient(0.0, 0.0, -0.625);
  Eigen::Matrix3d expectedHessian = Eigen::Matrix3d::Zero();
  expectedHessian(0, 0) = -5.0 / 8.0;
  expectedHessian(1, 1) = -5.0 / 8.0;
  expectedHessian(2, 2) = 1.25;

  // Points
  IRL::Pt xi(0.0, 0.0, 0.0);
  delta = 2;
  IRL::Pt x(0.0, 0.0, 1.0);
  Eigen::Vector3d gradTemp;
  Eigen::Matrix3d hessTemp;

  // Value Only
  Wendland::evaluate(xi, delta, x, &res1);
  EXPECT_NEAR(res1, expectedValue, std::numeric_limits<double>::epsilon())
      << "Evaluate Value Only Fail - Value";
  // Value and Gradient
  Wendland::evaluate(xi, delta, x, &res2);
  EXPECT_NEAR(std::get<0>(res2), expectedValue,
              std::numeric_limits<double>::epsilon())
      << "Evaluate Value And Gradient Fail - Value";

  gradTemp = std::get<1>(res2);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(gradTemp(i), expectedGradient(i),
                std::numeric_limits<double>::epsilon())
        << "Evaluate Value and Gradient Fail - Gradient Index " << i;
  }
  // Value and Gradient and Hessian
  Wendland::evaluate(xi, delta, x, &res3);
  EXPECT_NEAR(std::get<0>(res3), expectedValue,
              std::numeric_limits<double>::epsilon())
      << "Evaluate Value, Gradient, Hessian Fail - Value";

  gradTemp = std::get<1>(res3);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(gradTemp(i), expectedGradient(i),
                std::numeric_limits<double>::epsilon())
        << "Evaluate Value, Gradient, Hessian Fail - Gradient Index " << i;
  }

  hessTemp = std::get<2>(res3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_NEAR(hessTemp(i, j), expectedHessian(i, j),
                  std::numeric_limits<double>::epsilon())
          << "Evaluate Value, Gradient, Hessian Fail - Hessian Index " << i
          << "," << j;
    }
  }

  SUCCEED();
}

TEST(PartitionOfUnityImplicitSurface, Test) {
  std::vector<IRL::Pt> centroids;
  std::vector<IRL::SeparatorVariant> variantSeps;
  std::vector<double> weights;
  double delta = 5.0;

  // Construct a Plane
  IRL::Normal nor(1.0, 0.0, 0.0);  // Vertical Plane
  IRL::Plane plane = IRL::Plane(nor, 0.0);
  IRL::PlanarSeparator sep = IRL::PlanarSeparator::fromOnePlane(plane);
  variantSeps.push_back(sep);
  variantSeps.push_back(sep);
  variantSeps.push_back(sep);

  // Construct an Implicit Surface out of 3 Planes that are the same.
  IRL::Pt p0(0, 0, 0);
  IRL::Pt p1(0, 2, 0);
  IRL::Pt p2(0, -1, 0);
  centroids.push_back(p0);
  centroids.push_back(p1);
  centroids.push_back(p2);
  weights.push_back(1.0);
  weights.push_back(1.0);
  weights.push_back(1.0);

  IRL::PUImplicitSurface planarSurface(centroids, variantSeps, weights, delta);

  // Test Evaluate method, which in turn tests the implicitSeparator methods
  double res1;
  std::pair<double, Eigen::Vector3d> res2;
  std::tuple<double, Eigen::Vector3d, Eigen::Matrix3d> res3;

  IRL::Pt x1(0, 0, 0);
  IRL::Pt x2(1, 1, 0);
  // Value
  planarSurface.evaluate(x1, &res1);
  EXPECT_NEAR(res1, 0.0, std::numeric_limits<double>::epsilon())
      << "Function Value on Planar Separator is Nonzero";
  planarSurface.evaluate(x2, &res1);
  EXPECT_NEAR(res1, 1.0, std::numeric_limits<double>::epsilon())
      << "Function Value off Planar Separator is Wrong";

  // Gradient
  Eigen::Vector3d expectedGrad(1.0, 0.0, 0.0);
  planarSurface.evaluate(x1, &res2);
  Eigen::Vector3d grad = std::get<1>(res2);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(grad(i), expectedGrad(i),
                std::numeric_limits<double>::epsilon())
        << "Function Value on Planar Separator is Wrong";
  }

  planarSurface.evaluate(x2, &res2);
  grad = std::get<1>(res2);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(grad(i), expectedGrad(i),
                std::numeric_limits<double>::epsilon())
        << "Function Gradient off Planar Separator is Wrong";
  }

  // Hessian
  Eigen::Matrix3d expectedHessian = Eigen::Matrix3d::Zero();
  planarSurface.evaluate(x1, &res3);
  Eigen::Matrix3d hess = std::get<2>(res3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_NEAR(hess(i), expectedHessian(i),
                  std::numeric_limits<double>::epsilon())
          << "Function Hessian on Planar Separator is Wrong";
    }
  }

  planarSurface.evaluate(x2, &res3);
  hess = std::get<2>(res3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_NEAR(hess(i), expectedHessian(i),
                  std::numeric_limits<double>::epsilon())
          << "Function Hessian on Planar Separator is Wrong";
    }
  }
  // Intersect Edge Test
  IRL::Pt x0e(-1, 1, 0);
  IRL::Pt x1e(2, 0, 0);
  IRL::Pt expectedIntersection(0, 2.0 / 3.0, 0);
  std::vector<IRL::Pt> intersections =
      planarSurface.intersectEdge(x0e, x1e, 10, 0.2);
  EXPECT_EQ(intersections.size(), 1)
      << "Wrong Number of Planar Surface Intersections Found";
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(intersections[0][i], expectedIntersection[i], 1e-11)
        << "Planar Intersection Index " << i << " Wrong";
  }

  // Paraboloid Implicit Testing
  // Construct a Probaboloid
  variantSeps.clear();
  IRL::Normal e0(1, 0, 0);
  IRL::Normal e1(0, 1, 0);
  IRL::Normal e2(0, 0, 1);
  IRL::ReferenceFrame refFrame(e0, e1, e2);
  double a = 1;
  double b = 0;
  IRL::Pt datum(1, 0, 0);
  IRL::Paraboloid para = IRL::Paraboloid(datum, refFrame, a, b);
  variantSeps.push_back(para);
  variantSeps.push_back(para);
  variantSeps.push_back(para);

  SUCCEED();
}

TEST(PUReconstruction, Test1) {
  const int nlayers = 2;
  const int ncells = (1 + 2 * nlayers) * (1 + 2 * nlayers);

  std::vector<SeparatorVariant> planar_separator(ncells);
  std::vector<RectangularCuboid> cells(ncells);
  std::vector<Pt> centroids;
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

  // Compute end points of 2D PLIC, use to calculate centroids, and add to
  // neighborhood

  count = 0;
  StackVector<Pt, 2> intersections;
  PUSTNeighborhood<RectangularCuboid> neighborhood;
  const auto xy_plane = Plane(Normal(0.0, 0.0, 1.0), 0.0);
  for (int i = 0; i < 1 + 2 * nlayers; ++i) {
    for (int j = 0; j < 1 + 2 * nlayers; ++j) {
      auto separatorPtr =
          std::get_if<PlanarSeparator>(&planar_separator[count]);
      auto separator = *separatorPtr;
      const Polygon polygon = getPlanePolygonFromReconstruction<Polygon>(
          cells[count], separator, separator[0]);

      getIntersectionPts(polygon, xy_plane, &intersections);
      if (intersections.size() == 2) {
        // std::cout << "Start point for cell " << i << ", " << j << " = "
        //           << intersections[0] << std::endl;
        // std::cout << "  End point for cell " << i << ", " << j << " = "
        //           << intersections[1] << std::endl;

        Pt cen = (intersections[0] + intersections[1]) * 0.5;
        centroids.push_back(cen);
        neighborhood.addMember(&cen, &planar_separator[count]);
      }
      count++;
    }
  }

  std::cout << neighborhood.size() << "\n";
  // Now that we have neighborhood, print out normals and centroids to use for
  // some analytical stuffs
  for (int i = 0; i < neighborhood.size(); ++i) {
    auto cent = neighborhood.getCentroid(i);
    auto separatorPtr =
        std::get_if<PlanarSeparator>(&(neighborhood.getSeparator(i)));
    PlanarSeparator sep = *separatorPtr;
    auto plane0 = sep[0];
    auto normal = plane0.normal();
  }
  // Calculate Intersections
  IRL::Pt x1(0, 2, 0);
  IRL::Pt x0(0, 3, 0);
  PUST solver(neighborhood);
  PUImplicitSurface semi = solver.neighborhoodToImplicitSurface(5.0);
  std::vector<IRL::Pt> inters = semi.intersectEdge(x0, x1, 10, 0.2);
  EXPECT_EQ(inters.size(), 1) << "Wrong Number of Intersections Found";
  // Order goes x=0,x=1,y=2,x=2,y=1,y=0
  // Intersection Points, by hand
  IRL::Pt inter1Expected(0, 2.93192553235, 0);
  IRL::Pt inter2Expected(1, 2.53478876547624, 0);
  IRL::Pt inter3Expected(1.7590388618752, 2, 0);
  IRL::Pt inter4Expected(2, 1.7590388618752, 0);
  IRL::Pt inter5Expected(2.53478876547624, 1, 0);
  IRL::Pt inter6Expected(2.93192553235, 0, 0);
  std::vector<IRL::Pt> intersSet = {inter1Expected, inter2Expected,
                                    inter3Expected, inter4Expected,
                                    inter5Expected, inter6Expected};

  // Gradients
  Eigen::Vector3d grad1Expected(0.277899498135, 0.968982554617, 0);
  Eigen::Vector3d grad2Expected(0.458376679255, 0.853190619123, 0);
  Eigen::Vector3d grad3Expected(0.641841562767, 0.704583445946, 0);
  Eigen::Vector3d grad4Expected(0.704583445946, 0.641841562767, 0);
  Eigen::Vector3d grad5Expected(0.853190619123, 0.458376679255, 0);
  Eigen::Vector3d grad6Expected(0.968982554617, 0.277899498135, 0);
  std::vector<Eigen::Vector3d> gradSet = {grad1Expected, grad2Expected,
                                          grad3Expected, grad4Expected,
                                          grad5Expected, grad6Expected};

  // Forces
  IRL::Normal force1Expected(-0.968982554617, 0.277899498135, 0);
  force1Expected.normalize();
  IRL::Normal force2Expected(-0.853190619123, 0.458376679255, 0);
  force2Expected.normalize();
  IRL::Normal force3Expected(-0.704583445946, 0.641841562767, 0);
  force3Expected.normalize();
  IRL::Normal force4Expected(-0.641841562767, 0.704583445946, 0);
  force4Expected.normalize();
  IRL::Normal force5Expected(-0.458376679255, 0.853190619123, 0);
  force5Expected.normalize();
  IRL::Normal force6Expected(-0.277899498135, 0.968982554617, 0);
  force6Expected.normalize();
  std::vector<IRL::Normal> forceSet = {force1Expected, force2Expected,
                                       force3Expected, force4Expected,
                                       force5Expected, force6Expected};

  // Net Force
  IRL::Normal forceNetExpected(-0.685567462155, -0.685567462155, 0);
  IRL::Normal forceNetCalc(0, 0, 0);

  // Now that we can calculate everything using our methods.
  std::vector<IRL::Pt> x0Set = {IRL::Pt(0, 3, 0), IRL::Pt(1, 3, 0),
                                IRL::Pt(2, 2, 0), IRL::Pt(2, 2, 0),
                                IRL::Pt(3, 1, 0), IRL::Pt(3, 0, 0)};
  std::vector<IRL::Pt> x1Set = {IRL::Pt(0, 2, 0), IRL::Pt(1, 2, 0),
                                IRL::Pt(1, 2, 0), IRL::Pt(2, 1, 0),
                                IRL::Pt(2, 1, 0), IRL::Pt(2, 0, 0)};
  for (int i = 0; i < x0Set.size(); i++) {  // Loop Over Edges and Solve
    inters = semi.intersectEdge(x0Set[i], x1Set[i], 10, 0.2);
    // Check Intersection
    EXPECT_EQ(inters.size(), 1)
        << "Wrong Number of Intersections for Edge " << i;
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(inters[0][j], intersSet[i][j], 1e-9)
          << "Intersection " << i << " Index " << j << " Wrong";
    }
    // Check Gradient
    std::pair<double, Eigen::Vector3d> retVal;
    semi.evaluate(inters[0], &retVal);
    Eigen::Vector3d tempGrad = std::get<1>(retVal);
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(tempGrad[j], gradSet[i][j], 1e-9)
          << "Gradient " << i << " Index " << j << " Wrong";
    }
    // Calculate Forces
    IRL::Normal M = {0.0, 0.0, 0.0};
    double Pres = 0.0;
    IRL::Normal result = solver.solveEdge(1, x0Set[i], x1Set[i], 5.0, Pres, M);
    for (int j = 0; j < 3; j++) {
      EXPECT_NEAR(result[j], forceSet[i][j], 1e-9)
          << "Force " << i << " Index " << j << " Wrong";
    }
    // Net Force Adding
    if (i == 0) {
      forceNetCalc = forceNetCalc + result;
    }
    if (i == 5) {
      forceNetCalc = forceNetCalc + (-1 * result);
    }
  }
  for (int j = 0; j < 3; j++) {
    EXPECT_NEAR(forceNetCalc[j], forceNetExpected[j], 1e-9)
        << "Net Force Index " << j << " Wrong";
  }
  SUCCEED();
}

}  // namespace
