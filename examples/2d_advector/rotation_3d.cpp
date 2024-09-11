// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2024 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/paraboloid_advector/rotation_3d.h"

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

#include "examples/paraboloid_advector/data.h"
#include "examples/paraboloid_advector/reconstruction_types.h"
#include "examples/paraboloid_advector/solver.h"
#include "examples/paraboloid_advector/vof_advection.h"

constexpr int GC = 3;
constexpr IRL::Pt lower_domain(-0.5, -0.5, -0.5);
constexpr IRL::Pt upper_domain(0.5, 0.5, 0.5);

BasicMesh Rotation3D::setMesh(const int a_nx) {
  BasicMesh mesh(a_nx, a_nx, a_nx, GC);
  IRL::Pt my_lower_domain = lower_domain;
  IRL::Pt my_upper_domain = upper_domain;
  mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
  return mesh;
}

void Rotation3D::initialize(Data<double>* a_U, Data<double>* a_V,
                            Data<double>* a_W,
                            Data<IRL::Paraboloid>* a_interface,
                            const double a_time) {
  Rotation3D::setVelocity(a_time, a_U, a_V, a_W);
  const BasicMesh& mesh = a_U->getMesh();
  const double sphere_center_radius = std::sqrt(3.0 / (5.0 * 5.0));
  auto plane_normal = IRL::Normal(1.0, 2.0, 3.0);
  plane_normal.normalize();
  auto sphere_center_dir = IRL::Normal(1.0 / 5.0, 1.0 / 5.0, -1.0 / 5.0);
  sphere_center_dir.normalize();
  IRL::UnitQuaternion rotation(a_time * 2.0 * M_PI, plane_normal);
  sphere_center_dir = rotation * sphere_center_dir;
  const auto sphere_center = IRL::Pt(sphere_center_radius * sphere_center_dir);
  const double sphere_radius = 0.125;

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

void Rotation3D::setVelocity(const double a_time, Data<double>* a_U,
                             Data<double>* a_V, Data<double>* a_W) {
  const BasicMesh& mesh = a_U->getMesh();
  const double vel_scale = 2.0 * M_PI;
  auto plane_normal = IRL::Normal(1.0, 2.0, 3.0);
  plane_normal.normalize();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        auto loc = IRL::Pt(mesh.xm(i), mesh.ym(j), mesh.zm(k));
        const auto plane_loc =
            IRL::Pt(loc - plane_normal * IRL::dotProduct(plane_normal, loc));
        auto vel_dir = IRL::Normal(IRL::crossProduct(plane_normal, plane_loc));
        vel_dir.normalize();
        const double vel_mag = vel_scale * IRL::magnitude(plane_loc);
        (*a_U)(i, j, k) = vel_mag * vel_dir[0];
        (*a_V)(i, j, k) = vel_mag * vel_dir[1];
        (*a_W)(i, j, k) = vel_mag * vel_dir[2];
      }
    }
  }
}

const std::array<double, 3> Rotation3D::getExactVelocity(
    const IRL::Pt& a_location, const double a_time) {
  const double vel_scale = 2.0 * M_PI;
  auto plane_normal = IRL::Normal(1.0, 2.0, 3.0);
  plane_normal.normalize();
  const auto plane_loc = IRL::Pt(
      a_location - plane_normal * IRL::dotProduct(plane_normal, a_location));
  auto vel_dir = IRL::Normal(IRL::crossProduct(plane_normal, plane_loc));
  vel_dir.normalize();
  const double vel_mag = vel_scale * IRL::magnitude(plane_loc);
  return {vel_mag * vel_dir[0], vel_mag * vel_dir[1], vel_mag * vel_dir[2]};
}

const std::array<std::array<double, 3>, 3> Rotation3D::getExactVelocityGradient(
    const IRL::Pt& a_location, const double a_time) {
  const double epsilon = std::sqrt(1.0e-15);
  const IRL::Pt dudx =
      (IRL::Pt::fromRawDoublePointer(
           Rotation3D::getExactVelocity(
               IRL::Pt(a_location + epsilon * IRL::Pt(1.0, 0.0, 0.0)), a_time)
               .data()) -
       IRL::Pt::fromRawDoublePointer(
           Rotation3D::getExactVelocity(
               IRL::Pt(a_location - epsilon * IRL::Pt(1.0, 0.0, 0.0)), a_time)
               .data())) *
      (1.0 / (2.0 * epsilon));
  const IRL::Pt dudy =
      (IRL::Pt::fromRawDoublePointer(
           Rotation3D::getExactVelocity(
               IRL::Pt(a_location + epsilon * IRL::Pt(0.0, 1.0, 0.0)), a_time)
               .data()) -
       IRL::Pt::fromRawDoublePointer(
           Rotation3D::getExactVelocity(
               IRL::Pt(a_location - epsilon * IRL::Pt(0.0, 1.0, 0.0)), a_time)
               .data())) *
      (1.0 / (2.0 * epsilon));
  const IRL::Pt dudz =
      (IRL::Pt::fromRawDoublePointer(
           Rotation3D::getExactVelocity(
               IRL::Pt(a_location + epsilon * IRL::Pt(0.0, 0.0, 1.0)), a_time)
               .data()) -
       IRL::Pt::fromRawDoublePointer(
           Rotation3D::getExactVelocity(
               IRL::Pt(a_location - epsilon * IRL::Pt(0.0, 0.0, 1.0)), a_time)
               .data())) *
      (1.0 / (2.0 * epsilon));
  return std::array<std::array<double, 3>, 3>({{{dudx[0], dudy[0], dudz[0]},
                                                {dudx[1], dudy[1], dudz[1]},
                                                {dudx[2], dudy[2], dudz[2]}}});
}