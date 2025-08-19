// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/paraboloid_advector/stagnation_3d.h"

#include <float.h>
#include <chrono>
#include <cmath>
#include <iostream>

#include "irl/distributions/k_means.h"
#include "irl/distributions/partition_by_normal_vector.h"
#include "irl/generic_cutting/cut_polygon.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/moments/volume_moments.h"
#include "irl/parameters/constants.h"
#include "irl/planar_reconstruction/localized_separator_link.h"

#include "examples/paraboloid_advector/data.h"
#include "examples/paraboloid_advector/reconstruction_types.h"
#include "examples/paraboloid_advector/solver.h"
#include "examples/paraboloid_advector/vof_advection.h"

constexpr double T = 0.5;
constexpr int GC = 3;
constexpr IRL::Pt lower_domain(-0.5, -0.5, -0.5);
constexpr IRL::Pt upper_domain(0.5, 0.5, 0.5);

BasicMesh Stagnation3D::setMesh(const int a_nx) {
  BasicMesh mesh(a_nx, a_nx, a_nx, GC);
  IRL::Pt my_lower_domain = lower_domain;
  IRL::Pt my_upper_domain = upper_domain;
  mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
  return mesh;
}

void Stagnation3D::initialize(Data<double>* a_U, Data<double>* a_V,
                              Data<double>* a_W,
                              Data<IRL::Paraboloid>* a_interface,
                              const double a_time, const double final_time) {
  Stagnation3D::setVelocity(a_time, a_U, a_V, a_W);
  const BasicMesh& mesh = a_U->getMesh();
  const IRL::Pt sphere_center(0.0, 0.0, 0.0);
  const double sphere_radius = 0.25;

  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
        const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                    mesh.z(k + 1));
        const IRL::Pt mid_pt = 0.5 * (lower_cell_pt + upper_cell_pt);
        IRL::Pt disp = mid_pt - sphere_center;
        const auto mag = magnitude(disp);
        if (mag < sphere_radius - 2.0 * mesh.dx()) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysAbove();
        } else if (mag > sphere_radius + 2.0 * mesh.dx()) {
          (*a_interface)(i, j, k) = IRL::Paraboloid::createAlwaysBelow();
        } else {
          auto sphere_normal = IRL::Normal::fromPt(disp);
          sphere_normal.normalize();
          (*a_interface)(i, j, k) =
              details::fromSphere(sphere_center, sphere_radius, sphere_normal);
        }
      }
    }
  }
  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void Stagnation3D::setVelocity(const double a_time, Data<double>* a_U,
                               Data<double>* a_V, Data<double>* a_W) {
  const BasicMesh& mesh = a_U->getMesh();
  const IRL::Pt sphere_center(0.0, 0.0, 0.0);
  const double cos2pit = std::cos(M_PI * (a_time) / T);
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        const double x = mesh.xm(i) - sphere_center[0],
                     y = mesh.ym(j) - sphere_center[1],
                     z = mesh.zm(k) - sphere_center[2];
        (*a_U)(i, j, k) = 2.0 * x * cos2pit;
        (*a_V)(i, j, k) = -y * cos2pit;
        (*a_W)(i, j, k) = -z * cos2pit;
      }
    }
  }
}

const Eigen::Vector3d Stagnation3D::getExactVelocity(
    const Eigen::Vector3d& a_location, const double a_time) {
  const IRL::Pt sphere_center(0.0, 0.0, 0.0);
  const double x = a_location[0] - sphere_center[0],
               y = a_location[1] - sphere_center[1],
               z = a_location[2] - sphere_center[2];
  const double cos2pit = std::cos(M_PI * (a_time) / T);
  return Eigen::Vector3d({2.0 * x * cos2pit, -y * cos2pit, -z * cos2pit});
}

const Eigen::Matrix3d Stagnation3D::getExactVelocityGradient(
    const Eigen::Vector3d& a_location, const double a_time) {
  const double cos2pit = std::cos(M_PI * (a_time) / T);
  return Eigen::Matrix3d(
      {{2.0 * cos2pit, 0.0, 0.0}, {0.0, -cos2pit, 0.0}, {0.0, 0.0, -cos2pit}});
}

const Eigen::Matrix3d Stagnation3D::getExactVelocityHessianX(
    const Eigen::Vector3d& a_location, const double a_time) {
  return Eigen::Matrix3d::Zero();
}

const Eigen::Matrix3d Stagnation3D::getExactVelocityHessianY(
    const Eigen::Vector3d& a_location, const double a_time) {
  return Eigen::Matrix3d::Zero();
}

const Eigen::Matrix3d Stagnation3D::getExactVelocityHessianZ(
    const Eigen::Vector3d& a_location, const double a_time) {
  return Eigen::Matrix3d::Zero();
}