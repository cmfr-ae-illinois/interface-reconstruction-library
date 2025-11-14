// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/variant_advector/rotation_3d.h"

constexpr double T = 1.0;
constexpr int GC = 5;
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
                            Data<IRL::SeparatorVariant>* a_interface,
                            const double a_time, const double final_time) {
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
  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

void Rotation3D::setVelocity(const double a_time, Data<double>* a_U,
                             Data<double>* a_V, Data<double>* a_W) {
  const BasicMesh& mesh = a_U->getMesh();
  const double vel_scale = 2.0 * M_PI * T;
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

double Rotation3D::getTimeStep(const BasicMesh& a_mesh,
                               const double a_max_cfl) {
  const double vel_scale = 2.0 * M_PI * T;
  IRL::Normal plane_normal(1.0, 2.0, 3.0);
  plane_normal.normalize();

  double Umax = 0.0, Vmax = 0.0, Wmax = 0.0;

  for (int i = a_mesh.imin(); i <= a_mesh.imax(); i++) {
    for (int j = a_mesh.jmin(); j <= a_mesh.jmax(); j++) {
      for (int k = a_mesh.kmin(); k <= a_mesh.kmax(); k++) {
        auto loc = IRL::Pt(a_mesh.xm(i), a_mesh.ym(j), a_mesh.zm(k));
        const auto plane_loc =
            IRL::Pt(loc - plane_normal * IRL::dotProduct(plane_normal, loc));
        auto vel_dir = IRL::Normal(IRL::crossProduct(plane_normal, plane_loc));
        vel_dir.normalize();
        const double vel_mag = vel_scale * IRL::magnitude(plane_loc);
        const double u = vel_mag * vel_dir[0];
        const double v = vel_mag * vel_dir[1];
        const double w = vel_mag * vel_dir[2];

        Umax = std::max(Umax, std::abs(u));
        Vmax = std::max(Vmax, std::abs(v));
        Wmax = std::max(Wmax, std::abs(w));
      }
    }
  }

  return a_max_cfl * std::min(a_mesh.dx() / Umax,
                              std::min(a_mesh.dy() / Vmax, a_mesh.dz() / Wmax));
}
