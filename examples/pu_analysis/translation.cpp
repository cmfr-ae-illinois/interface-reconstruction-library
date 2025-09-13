// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2024 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/pu_analysis/translation.h"

#include <float.h>
#include <chrono>
#include <cmath>
#include <iostream>

#include "examples/pu_analysis/data.h"
#include "examples/pu_analysis/irl2d.h"
#include "examples/pu_analysis/reconstruction_types.h"
#include "examples/pu_analysis/vof_advection.h"
#include "irl/parameters/constants.h"

constexpr int GC = 5;

BasicMesh Translation::setMesh(const int a_nx) {
  BasicMesh mesh(a_nx, a_nx, GC);
  IRL2D::Vec my_lower_domain(-8, -8);
  IRL2D::Vec my_upper_domain(8, 8);
  mesh.setCellBoundaries(my_lower_domain, my_upper_domain);
  return mesh;
}

void Translation::setVelocity(const double a_time, Data<double>* a_U,
                              Data<double>* a_V) {
  const BasicMesh& mesh = a_U->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      (*a_U)(i, j) = 1.0;
      (*a_V)(i, j) = 0.0;
    }
  }
}

const IRL2D::Vec Translation::getExactVelocity2D(double t,
                                                 const IRL2D::Vec& P) {
  return IRL2D::Vec{1.0, 0.0};
}

const IRL2D::Mat Translation::getExactVelocityGradient2D(double t,
                                                         const IRL2D::Vec& P) {
  return IRL2D::Mat(IRL2D::Vec{0.0, 0.0}, IRL2D::Vec{0.0, 0.0});
}
