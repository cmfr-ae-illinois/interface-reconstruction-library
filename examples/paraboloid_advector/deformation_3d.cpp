// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/paraboloid_advector/deformation_3d.h"

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

constexpr double T = 0.1;
constexpr int GC = 3;
constexpr IRL::Pt lower_domain(0.0, 0.0, 0.0);
constexpr IRL::Pt upper_domain(1.0, 1.0, 1.0);

BasicMesh Deformation3D::setMesh(const int a_nx) {
  BasicMesh mesh(a_nx, a_nx, a_nx, GC);
  IRL::Pt my_lower_domain = lower_domain;
  IRL::Pt my_upper_domain = upper_domain;
  mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
  return mesh;
}

void Deformation3D::initialize(Data<double>* a_U, Data<double>* a_V,
                               Data<double>* a_W,
                               Data<IRL::Paraboloid>* a_interface,
                               const double a_time, const double final_time) {
  Deformation3D::setVelocity(0.0, a_U, a_V, a_W);
  const BasicMesh& mesh = a_U->getMesh();
  const IRL::Pt sphere_center(0.35, 0.35, 0.35);
  const double sphere_radius = 0.15;

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

void Deformation3D::setVelocity(const double a_time, Data<double>* a_U,
                                Data<double>* a_V, Data<double>* a_W) {
  const BasicMesh& mesh = a_U->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        const double sinpix = std::sin(M_PI * mesh.xm(i)),
                     sinpiy = std::sin(M_PI * mesh.ym(j)),
                     sinpiz = std::sin(M_PI * mesh.zm(k));
        const double sin2pix = std::sin(2.0 * M_PI * mesh.xm(i)),
                     sin2piy = std::sin(2.0 * M_PI * mesh.ym(j)),
                     sin2piz = std::sin(2.0 * M_PI * mesh.zm(k));
        const double cospit = std::cos(M_PI * (a_time) / T);
        (*a_U)(i, j, k) = 2.0 * sinpix * sinpix * sin2piy * sin2piz * cospit;
        (*a_V)(i, j, k) = -sinpiy * sinpiy * sin2pix * sin2piz * cospit;
        (*a_W)(i, j, k) = -sinpiz * sinpiz * sin2pix * sin2piy * cospit;
      }
    }
  }
}

const Eigen::Vector3d Deformation3D::getExactVelocity(
    const Eigen::Vector3d& a_location, const double a_time) {
  const double sinpix = std::sin(M_PI * a_location[0]),
               sinpiy = std::sin(M_PI * a_location[1]),
               sinpiz = std::sin(M_PI * a_location[2]);
  const double sin2pix = std::sin(2.0 * M_PI * a_location[0]),
               sin2piy = std::sin(2.0 * M_PI * a_location[1]),
               sin2piz = std::sin(2.0 * M_PI * a_location[2]);
  const double cospit = std::cos(M_PI * (a_time) / T);
  return Eigen::Vector3d({2.0 * sinpix * sinpix * sin2piy * sin2piz * cospit,
                          -sinpiy * sinpiy * sin2pix * sin2piz * cospit,
                          -sinpiz * sinpiz * sin2pix * sin2piy * cospit});
}

const Eigen::Matrix3d Deformation3D::getExactVelocityGradient(
    const Eigen::Vector3d& a_location, const double a_time) {
  const double sinpix = std::sin(M_PI * a_location[0]),
               sinpiy = std::sin(M_PI * a_location[1]),
               sinpiz = std::sin(M_PI * a_location[2]);
  const double cospix = std::cos(M_PI * a_location[0]),
               cospiy = std::cos(M_PI * a_location[1]),
               cospiz = std::cos(M_PI * a_location[2]);
  const double sin2pix = std::sin(2.0 * M_PI * a_location[0]),
               sin2piy = std::sin(2.0 * M_PI * a_location[1]),
               sin2piz = std::sin(2.0 * M_PI * a_location[2]);
  const double cos2pix = std::cos(2.0 * M_PI * a_location[0]),
               cos2piy = std::cos(2.0 * M_PI * a_location[1]),
               cos2piz = std::cos(2.0 * M_PI * a_location[2]);
  const double cospit = std::cos(M_PI * (a_time) / T);
  return Eigen::Matrix3d(
      {{4.0 * M_PI * cospix * sinpix * sin2piy * sin2piz * cospit,
        4.0 * M_PI * sinpix * sinpix * cos2piy * sin2piz * cospit,
        4.0 * M_PI * sinpix * sinpix * sin2piy * cos2piz * cospit},
       {-2.0 * M_PI * sinpiy * sinpiy * cos2pix * sin2piz * cospit,
        -2.0 * M_PI * cospiy * sinpiy * sin2pix * sin2piz * cospit,
        -2.0 * M_PI * sinpiy * sinpiy * sin2pix * cos2piz * cospit},
       {-2.0 * M_PI * sinpiz * sinpiz * cos2pix * sin2piy * cospit,
        -2.0 * M_PI * cospiz * sinpiz * sin2pix * sin2piy * cospit,
        -2.0 * M_PI * sinpiz * sinpiz * sin2pix * cos2piy * cospit}});
}

const Eigen::Matrix3d Deformation3D::getExactVelocityHessianX(
    const Eigen::Vector3d& a_location, const double a_time) {
  const double sinpix = std::sin(M_PI * a_location[0]),
               sinpiy = std::sin(M_PI * a_location[1]),
               sinpiz = std::sin(M_PI * a_location[2]);
  const double cospix = std::cos(M_PI * a_location[0]),
               cospiy = std::cos(M_PI * a_location[1]),
               cospiz = std::cos(M_PI * a_location[2]);
  const double sin2pix = std::sin(2.0 * M_PI * a_location[0]),
               sin2piy = std::sin(2.0 * M_PI * a_location[1]),
               sin2piz = std::sin(2.0 * M_PI * a_location[2]);
  const double cos2pix = std::cos(2.0 * M_PI * a_location[0]),
               cos2piy = std::cos(2.0 * M_PI * a_location[1]),
               cos2piz = std::cos(2.0 * M_PI * a_location[2]);
  const double cospit = std::cos(M_PI * (a_time) / T);
  return Eigen::Matrix3d(
      {{4.0 * M_PI * M_PI * (-sinpix * sinpix + cospix * cospix) * sin2piy *
            sin2piz * cospit,
        8.0 * M_PI * M_PI * cospix * sinpix * cos2piy * sin2piz * cospit,
        8.0 * M_PI * M_PI * cospix * sinpix * sin2piy * cos2piz * cospit},
       {8.0 * M_PI * M_PI * cospix * sinpix * cos2piy * sin2piz * cospit,
        -8.0 * M_PI * M_PI * sinpix * sinpix * sin2piy * sin2piz * cospit,
        8.0 * M_PI * M_PI * sinpix * sinpix * cos2piy * cos2piz * cospit},
       {8.0 * M_PI * M_PI * cospix * sinpix * sin2piy * cos2piz * cospit,
        8.0 * M_PI * M_PI * sinpix * sinpix * cos2piy * cos2piz * cospit,
        -8.0 * M_PI * M_PI * sinpix * sinpix * sin2piy * sin2piz * cospit}});
}

const Eigen::Matrix3d Deformation3D::getExactVelocityHessianY(
    const Eigen::Vector3d& a_location, const double a_time) {
  const double sinpix = std::sin(M_PI * a_location[0]),
               sinpiy = std::sin(M_PI * a_location[1]),
               sinpiz = std::sin(M_PI * a_location[2]);
  const double cospix = std::cos(M_PI * a_location[0]),
               cospiy = std::cos(M_PI * a_location[1]),
               cospiz = std::cos(M_PI * a_location[2]);
  const double sin2pix = std::sin(2.0 * M_PI * a_location[0]),
               sin2piy = std::sin(2.0 * M_PI * a_location[1]),
               sin2piz = std::sin(2.0 * M_PI * a_location[2]);
  const double cos2pix = std::cos(2.0 * M_PI * a_location[0]),
               cos2piy = std::cos(2.0 * M_PI * a_location[1]),
               cos2piz = std::cos(2.0 * M_PI * a_location[2]);
  const double cospit = std::cos(M_PI * (a_time) / T);
  return Eigen::Matrix3d(
      {{4.0 * M_PI * M_PI * sinpiy * sinpiy * sin2pix * sin2piz * cospit,
        -4.0 * M_PI * M_PI * cospiy * sinpiy * cos2pix * sin2piz * cospit,
        -4.0 * M_PI * M_PI * sinpiy * sinpiy * cos2pix * cos2piz * cospit},
       {-4.0 * M_PI * M_PI * cospiy * sinpiy * cos2pix * sin2piz * cospit,
        -2.0 * M_PI * M_PI * (-sinpiy * sinpiy + cospiy * cospiy) * sin2pix *
            sin2piz * cospit,
        -4.0 * M_PI * M_PI * cospiy * sinpiy * sin2pix * cos2piz * cospit},
       {-4.0 * M_PI * M_PI * sinpiy * sinpiy * cos2pix * cos2piz * cospit,
        -4.0 * M_PI * M_PI * cospiy * sinpiy * sin2pix * cos2piz * cospit,
        4.0 * M_PI * M_PI * sinpiy * sinpiy * sin2pix * sin2piz * cospit}});
}

const Eigen::Matrix3d Deformation3D::getExactVelocityHessianZ(
    const Eigen::Vector3d& a_location, const double a_time) {
  const double sinpix = std::sin(M_PI * a_location[0]),
               sinpiy = std::sin(M_PI * a_location[1]),
               sinpiz = std::sin(M_PI * a_location[2]);
  const double cospix = std::cos(M_PI * a_location[0]),
               cospiy = std::cos(M_PI * a_location[1]),
               cospiz = std::cos(M_PI * a_location[2]);
  const double sin2pix = std::sin(2.0 * M_PI * a_location[0]),
               sin2piy = std::sin(2.0 * M_PI * a_location[1]),
               sin2piz = std::sin(2.0 * M_PI * a_location[2]);
  const double cos2pix = std::cos(2.0 * M_PI * a_location[0]),
               cos2piy = std::cos(2.0 * M_PI * a_location[1]),
               cos2piz = std::cos(2.0 * M_PI * a_location[2]);
  const double cospit = std::cos(M_PI * (a_time) / T);
  return Eigen::Matrix3d(
      {{4.0 * M_PI * M_PI * sinpiz * sinpiz * sin2pix * sin2piy * cospit,
        -4.0 * M_PI * M_PI * sinpiz * sinpiz * cos2pix * cos2piy * cospit,
        -4.0 * M_PI * M_PI * cospiz * sinpiz * cos2pix * sin2piy * cospit},
       {-4.0 * M_PI * M_PI * sinpiz * sinpiz * cos2pix * cos2piy * cospit,
        4.0 * M_PI * M_PI * sinpiz * sinpiz * sin2pix * sin2piy * cospit,
        -4.0 * M_PI * M_PI * cospiz * sinpiz * sin2pix * cos2piy * cospit},
       {-4.0 * M_PI * M_PI * cospiz * sinpiz * cos2pix * sin2piy * cospit,
        -4.0 * M_PI * M_PI * cospiz * sinpiz * sin2pix * cos2piy * cospit,
        -2.0 * M_PI * M_PI * (-sinpiz * sinpiz + cospiz * cospiz) * sin2pix *
            sin2piy * cospit}});
}