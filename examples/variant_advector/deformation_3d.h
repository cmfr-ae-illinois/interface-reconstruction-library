// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2025 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_VARIANT_ADVECTOR_DEFORMATION_3D_H_
#define EXAMPLES_VARIANT_ADVECTOR_DEFORMATION_3D_H_

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/moments/volume_moments.h"
#include "irl/variant_reconstruction/separator_variant.h"

#include "examples/variant_advector/basic_mesh.h"
#include "examples/variant_advector/data.h"
#include "examples/variant_advector/reconstruction_types.h"

struct Deformation3D {
  static BasicMesh setMesh(const int a_nx);

  static void initialize(Data<double>* a_U, Data<double>* a_V,
                         Data<double>* a_W,
                         Data<IRL::SeparatorVariant>* a_interface,
                         const double a_time, const double final_time);

  static void setVelocity(const double a_time, Data<double>* a_U,
                          Data<double>* a_V, Data<double>* a_W);

  static double getTimeStep(const BasicMesh& a_mesh, const double a_max_cfl);
};

#endif  // EXAMPLES_VARIANT_ADVECTOR_DEFORMATION_3D_H_
