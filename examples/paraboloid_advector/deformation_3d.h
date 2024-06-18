// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_PARABOLOID_ADVECTOR_DEFORMATION_3D_H_
#define EXAMPLES_PARABOLOID_ADVECTOR_DEFORMATION_3D_H_

#include "irl/paraboloid_reconstruction/paraboloid.h"

#include "examples/paraboloid_advector/basic_mesh.h"
#include "examples/paraboloid_advector/data.h"

struct Deformation3D {
  static BasicMesh setMesh(const int a_nx);

  static void initialize(Data<double>* a_U, Data<double>* a_V,
                         Data<double>* a_W, Data<IRL::Paraboloid>* a_interface,
                         const double a_time);

  static void setVelocity(const double a_time, Data<double>* a_U,
                          Data<double>* a_V, Data<double>* a_W);

  static const std::array<double, 3> getExactVelocity(const IRL::Pt& a_location,
                                                      const double a_time);

  static const std::array<std::array<double, 3>, 3> getExactVelocityGradient(
      const IRL::Pt& a_location, const double a_time);
};

#endif  // EXAMPLES_PARABOLOID_ADVECTOR_DEFORMATION_2D_H_
