// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/variant_advector/deformation_3d.h"

constexpr double T = 3.0;
constexpr int GC = 5;
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
                               Data<IRL::SeparatorVariant>* a_interface,
                               const double a_time, const double final_time) {
  Deformation3D::setVelocity(a_time, a_U, a_V, a_W);
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

double Deformation3D::getTimeStep(const BasicMesh& a_mesh,
                                  const double a_max_cfl) {
  const double Umax = 2.0, Vmax = 1.0, Wmax = 1.0;
  return a_max_cfl * std::min(a_mesh.dx() / Umax,
                              std::min(a_mesh.dy() / Vmax, a_mesh.dz() / Wmax));
}
