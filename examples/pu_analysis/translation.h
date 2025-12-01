// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2024 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_PU_ANALYSIS_TRANSLATION_H_
#define EXAMPLES_PU_ANALYSIS_TRANSLATION_H_

#include "irl/paraboloid_reconstruction/paraboloid.h"

#include "examples/pu_analysis/basic_mesh.h"
#include "examples/pu_analysis/data.h"
#include "examples/pu_analysis/irl2d.h"

struct Translation {
  static BasicMesh setMesh(const int a_nx);

  static void setVelocity(const double a_time, Data<double>* a_U,
                          Data<double>* a_V);

  static const IRL2D::Vec getExactVelocity2D(double t, const IRL2D::Vec& P);
  static const IRL2D::Mat getExactVelocityGradient2D(double t,
                                                     const IRL2D::Vec& P);
};

#endif  // EXAMPLES_PU_ANALYSIS_TRANSLATION_H_
