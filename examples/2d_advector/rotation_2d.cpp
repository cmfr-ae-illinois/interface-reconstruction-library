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
#include "examples/2d_advector/irl2d.h"
#include "examples/2d_advector/reconstruction_types.h"
#include "examples/2d_advector/solver.h"
#include "examples/2d_advector/vof_advection.h"

constexpr int GC = 3;
constexpr IRL2D::Vec lower_domain(-0.5, -0.5);
constexpr IRL2D::Vec upper_domain(0.5, 0.5);

BasicMesh Rotation2D::setMesh(const int a_nx) {
  BasicMesh mesh(a_nx, a_nx, GC);
  IRL2D::Vec my_lower_domain = lower_domain;
  IRL2D::Vec my_upper_domain = upper_domain;
  mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
  return mesh;
}

void Rotation2D::initialize(Data<double>* a_U, Data<double>* a_V,
                            Data<IRL2D::Parabola>* a_interface,
                            const double a_time) {
  Rotation2D::setVelocity(a_time, a_U, a_V);
  const BasicMesh& mesh = a_U->getMesh();
  const auto circle_center = IRL2D::Vec(0.0, 0.0);
  const double circle_radius = 0.25;

  // Loop over cells in domain. Skip if cell is not mixed phase.
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      const IRL2D::Vec lower_cell_pt(mesh.x(i), mesh.y(j));
      const IRL2D::Vec upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1));
      const IRL2D::Vec mid_pt = 0.5 * (lower_cell_pt + upper_cell_pt);
      IRL2D::Vec disp = mid_pt - circle_center;
      const auto mag = disp.magnitude();
      if (mag < circle_radius - 2.0 * mesh.dx()) {
        (*a_interface)(i, j).markAsAlwaysAbove();
      } else if (mag > circle_radius + 2.0 * mesh.dx()) {
        (*a_interface)(i, j).markAsAlwaysBelow();
      } else {
        auto circle_normal = IRL2D::Vec(disp);
        circle_normal.normalize();
        const IRL2D::Vec datum = circle_center + circle_radius * circle_normal;
        const IRL2D::ReferenceFrame frame(
            IRL2D::Vec(circle_normal[1], -circle_normal[0]), circle_normal);
        const double coeff = 1.0 / (2.0 * circle_radius);
        (*a_interface)(i, j) = IRL2D::Parabola(datum, frame, coeff);
        for (int ii = 0; ii < 10; ii++) {
          auto cell = IRL2D::RectangleFromBounds(lower_cell_pt, upper_cell_pt);
          auto moments = IRL2D::ComputeMoments(cell, (*a_interface)(i, j));
          const double liquid_volume_fraction =
              moments.m0() / IRL2D::ComputeVolume(cell);
          if (liquid_volume_fraction >= IRL::global_constants::VF_LOW &&
              liquid_volume_fraction <= IRL::global_constants::VF_HIGH) {
            const auto new_disp = moments.m1() / moments.m0() - circle_center;
            auto new_circle_normal = IRL2D::Vec(new_disp);
            new_circle_normal.normalize();
            const IRL2D::Vec new_datum =
                circle_center + circle_radius * new_circle_normal;
            const IRL2D::ReferenceFrame new_frame(
                IRL2D::Vec(new_circle_normal[1], -new_circle_normal[0]),
                new_circle_normal);
            (*a_interface)(i, j) = IRL2D::Parabola(new_datum, new_frame, coeff);
          } else {
            break;
          }
        }
      }
    }
  }
  // Update border with simple ghost-cell fill and correct datum for
  // assumed periodic boundary
  a_interface->updateBorder();
  correctInterfaceBorders(a_interface);
}

void Rotation2D::setVelocity(const double a_time, Data<double>* a_U,
                             Data<double>* a_V) {
  const double vel_scale = 2.0 * M_PI;
  const BasicMesh& mesh = a_U->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      auto loc = IRL2D::Vec(mesh.xm(i), mesh.ym(j));
      (*a_U)(i, j) = -vel_scale * loc[1];
      (*a_V)(i, j) = vel_scale * loc[0];
    }
  }
}

const IRL2D::Vec Rotation2D::getExactVelocity2D(double t, const IRL2D::Vec& P) {
  const double vel_scale = 2.0 * M_PI;
  return IRL2D::Vec{-vel_scale * P.y(), vel_scale * P.x()};
  // return IRL2D::Vec{P.y() * P.y() * P.y(), P.x() * P.x() * P.x()};
  // return IRL2D::Vec{1.0 - P.y() * P.y(), 1.0 + P.x() * P.x()};
  // return IRL2D::Vec{1.0 - P.y(), 1.0 + P.x()};
  // return IRL2D::Vec{1.0, 1.0};
}

const IRL2D::Mat Rotation2D::getExactVelocityGradient2D(double t,
                                                        const IRL2D::Vec& P) {
  const double vel_scale = 2.0 * M_PI;
  return IRL2D::Mat(IRL2D::Vec{0.0, -vel_scale}, IRL2D::Vec{vel_scale, 0.0});
  // return IRL2D::Mat(IRL2D::Vec{0.0, 3.0 * P.y() * P.y()},
  //                   IRL2D::Vec{3.0 * P.x() * P.x(), 0.0});
  // return IRL2D::Mat(IRL2D::Vec{0.0, -2.0 * P.y()},
  //                   IRL2D::Vec{2.0 * P.x(), 0.0});
  // return IRL2D::Mat(IRL2D::Vec{0.0, -1.0}, IRL2D::Vec{1.0, 0.0});
  // return IRL2D::Mat(IRL2D::Vec{0.0, 0.0}, IRL2D::Vec{0.0, 0.0});
}
