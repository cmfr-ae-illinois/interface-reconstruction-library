// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/surface_advector/translation_3d.h"

constexpr int GC = 3;
constexpr IRL::Pt lower_domain(0.0, 0.0, 0.0);
constexpr IRL::Pt upper_domain(1.0, 1.0, 1.0);

BasicMesh Translation3D::setMesh(const int a_nx) {
  BasicMesh mesh(a_nx, a_nx, a_nx, GC);
  IRL::Pt my_lower_domain = lower_domain;
  IRL::Pt my_upper_domain = upper_domain;
  mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
  return mesh;
}

void Translation3D::initialize(Data<double>* a_U, Data<double>* a_V,
                               Data<double>* a_W,
                               Data<IRL::SeparatorVariant>* a_interface,
                               const double a_time, const double final_time) {
  Translation3D::setVelocity(a_time, a_U, a_V, a_W);
  const BasicMesh& mesh = a_U->getMesh();
  const IRL::Pt sphere_center(std::fmod(0.5 + a_time * 1.0, 1.0),
                              std::fmod(0.5 + a_time * 1.0 / 1.5, 1.0),
                              std::fmod(0.5 + a_time * 1.0 / 3.0, 1.0));
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
        const auto mag = IRL::magnitude(disp);
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
  correctInterfaceBorders(a_interface);
}

void Translation3D::setVelocity(const double a_time, Data<double>* a_U,
                                Data<double>* a_V, Data<double>* a_W) {
  const BasicMesh& mesh = a_U->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_U)(i, j, k) = 1.0;
        (*a_V)(i, j, k) = 1.0 / 1.5;
        (*a_W)(i, j, k) = 1.0 / 3.0;
      }
    }
  }
}

double Translation3D::getTimeStep(const BasicMesh& a_mesh,
                                  const double a_max_cfl) {
  const double Umax = 1.0, Vmax = 1.0 / 1.5, Wmax = 1.0 / 3.0;
  return a_max_cfl * std::min(a_mesh.dx() / Umax,
                              std::min(a_mesh.dy() / Vmax, a_mesh.dz() / Wmax));
}
