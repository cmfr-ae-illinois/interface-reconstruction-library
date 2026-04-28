// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2026 Parin Trivedi <parin.trivedi@hotmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cmath>
#include <iomanip>
#include <Eigen/Dense>

#include "gtest/gtest.h"

#include "irl/quadratic_reconstruction/coons_patch.h"
#include "irl/variant_reconstruction/separator_variant.h"
#include "irl/interface_reconstruction_methods/pu.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"


namespace {

using namespace IRL;

void getIntersectionPts(const IRL::Polygon& a_polygon,
                        const IRL::Plane& a_cutting_plane,
                        IRL::StackVector<IRL::Pt, 2>* a_intersection_pts) {
  a_intersection_pts->resize(0);
  double distance = a_cutting_plane.signedDistanceToPoint(a_polygon[0]);
  for (int n = 0; n < a_polygon.getNumberOfVertices(); ++n) {
    const int next_id = (n + 1) % a_polygon.getNumberOfVertices();
    double next_distance =
        a_cutting_plane.signedDistanceToPoint(a_polygon[next_id]);
    if (distance * next_distance < 0.0) {
      a_intersection_pts->push_back(IRL::Pt::fromEdgeIntersection(
          a_polygon[n], distance, a_polygon[next_id], next_distance));
      if (a_intersection_pts->size() == 2) {
        break;
      }
    }
    distance = next_distance;
  }
}

TEST(PU, GradientAndHessianTest) {
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

  count = 0;
  StackVector<Pt, 2> intersections;
  PUNeighborhood neighborhood;
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
        Pt cen = (intersections[0] + intersections[1]) * 0.5;
        centroids.push_back(cen);
        neighborhood.addMember(planar_separator[count], cen, 1.0);
      }
      count++;
    }
  }

  // PU class
  PU pu(&neighborhood, 5.0);

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

  // Hessians
  Eigen::Matrix3d hess1Expected = Eigen::Matrix3d::Zero();
  hess1Expected(0, 0) = 0.132679691006;
  hess1Expected(0, 1) = -0.0787101181139;
  hess1Expected(1, 0) = -0.0787101181139;
  hess1Expected(1, 1) = 0.042814285345;
  Eigen::Matrix3d hess2Expected = Eigen::Matrix3d::Zero();
  hess2Expected(0, 0) = 0.165582686659;
  hess2Expected(0, 1) = -0.11415777827;
  hess2Expected(1, 0) = -0.11415777827;
  hess2Expected(1, 1) = 0.0785849650454;
  Eigen::Matrix3d hess3Expected = Eigen::Matrix3d::Zero();
  hess3Expected(0, 0) = 0.137487487701;
  hess3Expected(0, 1) = -0.130895703411;
  hess3Expected(1, 0) = -0.130895703411;
  hess3Expected(1, 1) = 0.123584037533;
  Eigen::Matrix3d hess4Expected = Eigen::Matrix3d::Zero();
  hess4Expected(0, 0) = 0.123584037533;
  hess4Expected(0, 1) = -0.130895703411;
  hess4Expected(1, 0) = -0.130895703411;
  hess4Expected(1, 1) = 0.137487487701;
  Eigen::Matrix3d hess5Expected = Eigen::Matrix3d::Zero();
  hess5Expected(0, 0) = 0.0785849650454;
  hess5Expected(0, 1) = -0.11415777827;
  hess5Expected(1, 0) = -0.11415777827;
  hess5Expected(1, 1) = 0.165582686659;
  Eigen::Matrix3d hess6Expected = Eigen::Matrix3d::Zero();
  hess6Expected(0, 0) = 0.042814285345;
  hess6Expected(0, 1) = -0.0787101181139;
  hess6Expected(1, 0) = -0.0787101181139;
  hess6Expected(1, 1) = 0.132679691006;
  std::vector<Eigen::Matrix3d> hessSet = {hess1Expected, hess2Expected,
                                          hess3Expected, hess4Expected,
                                          hess5Expected, hess6Expected};

  // solving for gradient and hessian at intersection points
  for (size_t i = 0; i < intersSet.size(); ++i) {
    const auto [PU_value, PU_grad, PU_hess] =
        pu.getPUAndGradAndHessian(intersSet[i]);
    // outputting gradient value
    std::cout << "PU grad at inter " << i + 1 << " is " << PU_grad.transpose()
              << std::endl;
    std::cout << "PU hess at inter " << i + 1 << " is " << std::endl
              << PU_hess << std::endl;
  }
}

}  // namespace
