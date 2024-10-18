// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2024 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/2d_advector/rotation_2d.h"

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
#include "irl/interface_reconstruction_methods/progressive_distance_solver_paraboloid.h"
#include "irl/moments/volume_moments.h"
#include "irl/parameters/constants.h"
#include "irl/planar_reconstruction/localized_separator_link.h"

#include "examples/2d_advector/data.h"
#include "examples/2d_advector/reconstruction_types.h"
#include "examples/2d_advector/solver.h"
#include "examples/2d_advector/vof_advection.h"

constexpr int GC = 3;
constexpr IRL::Pt lower_domain(-0.5, -0.5, -0.5);
constexpr IRL::Pt upper_domain(0.5, 0.5, 0.5);

BasicMesh Rotation2D::setMesh(const int a_nx) {
  BasicMesh mesh(a_nx, a_nx, 1, GC);
  IRL::Pt my_lower_domain = lower_domain;
  IRL::Pt my_upper_domain = upper_domain;
  mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
  return mesh;
}

void Rotation2D::initialize(Data<double>* a_U, Data<double>* a_V,
                            Data<double>* a_W,
                            Data<IRL::Paraboloid>* a_interface,
                            const double a_time) {
  Rotation2D::setVelocity(a_time, a_U, a_V, a_W);
  const BasicMesh& mesh = a_U->getMesh();
  const auto sphere_center = IRL::Pt(0.0, 0.0, 0.0);
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
          for (int ii = 0; ii < 10; ii++) {
            auto cell = IRL::RectangularCuboid::fromBoundingPts(lower_cell_pt,
                                                                upper_cell_pt);
            auto moments = IRL::getVolumeMoments<IRL::VolumeMoments>(
                cell, (*a_interface)(i, j, k));
            const double liquid_volume_fraction =
                moments.volume() / mesh.cell_volume();
            if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
                liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
              auto centroid = IRL::Pt(moments.centroid() / moments.volume() -
                                      sphere_center);
              sphere_normal = IRL::Normal::fromPt(centroid);
              sphere_normal.normalize();
              (*a_interface)(i, j, k) = details::fromSphere(
                  sphere_center, sphere_radius, sphere_normal);
            } else {
              break;
            }
          }
        }
      }
    }
  }
  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfacePlaneBorders(a_interface);
}

void Rotation2D::setVelocity(const double a_time, Data<double>* a_U,
                             Data<double>* a_V, Data<double>* a_W) {
  const BasicMesh& mesh = a_U->getMesh();
  const double vel_scale = 2.0 * M_PI;
  auto plane_normal = IRL::Normal(0, 0, 1);
  plane_normal.normalize();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        auto loc = IRL::Pt(mesh.xm(i), mesh.ym(j), mesh.zm(k));
        (*a_U)(i, j, k) = -vel_scale * loc[1];
        (*a_V)(i, j, k) = vel_scale * loc[0];
        (*a_W)(i, j, k) = 0.0;
      }
    }
  }
}

const std::array<double, 3> Rotation2D::getExactVelocity(
    const IRL::Pt& a_location, const double a_time) {
  const double vel_scale = 2.0 * M_PI;
  return {-vel_scale * a_location[1], vel_scale * a_location[0], 0.0};
}

const std::array<std::array<double, 3>, 3> Rotation2D::getExactVelocityGradient(
    const IRL::Pt& a_location, const double a_time) {
  const double vel_scale = 2.0 * M_PI;
  return {-vel_scale * a_location[1], vel_scale * a_location[0], 0.0};
  const auto dudx = IRL::Pt(0.0, vel_scale, 0.0);
  const auto dudy = IRL::Pt(-vel_scale, 0.0, 0.0);
  const auto dudz = IRL::Pt(0.0, 0.0, 0.0);
  return std::array<std::array<double, 3>, 3>({{{dudx[0], dudy[0], dudz[0]},
                                                {dudx[1], dudy[1], dudz[1]},
                                                {dudx[2], dudy[2], dudz[2]}}});
}

const IRL2D::Vec Rotation2D::getExactVelocity2D(double t, const IRL2D::Vec& P) {
  const double vel_scale = 2.0 * M_PI;
  return IRL2D::Vec{0.25 + P.y() * P.y() * P.y(), 0.25 + P.x() * P.x() * P.x()};
  // return IRL2D::Vec{1.0 - P.y() * P.y(), 1.0 + P.x() * P.x()};
  // return IRL2D::Vec{1.0 - P.y(), 1.0 + P.x()};
  // return IRL2D::Vec{1.0, 1.0};
}

const IRL2D::Mat Rotation2D::getExactVelocityGradient2D(double t,
                                                        const IRL2D::Vec& P) {
  const double vel_scale = 2.0 * M_PI;
  return IRL2D::Mat(IRL2D::Vec{0.0, 3.0 * P.y() * P.y()},
                    IRL2D::Vec{3.0 * P.x() * P.x(), 0.0});
  // return IRL2D::Mat(IRL2D::Vec{0.0, -2.0 * P.y()},
  //                   IRL2D::Vec{2.0 * P.x(), 0.0});
  // return IRL2D::Mat(IRL2D::Vec{0.0, -1.0}, IRL2D::Vec{1.0, 0.0});
  // return IRL2D::Mat(IRL2D::Vec{0.0, 0.0}, IRL2D::Vec{0.0, 0.0});
}
